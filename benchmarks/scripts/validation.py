#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["polars>=1.0", "matplotlib>=3.8", "numpy>=1.26", "rich>=13.0", "typer>=0.9.0"]
# ///
"""SASA validation: compare accuracy across tools.

Runs tools on a PDB directory, collects per-file SASA values,
and compares across tools -- completely independent of timing benchmarks.

Usage:
    # Run main validation (4 n_points levels, all zsasa variants + rustsasa)
    ./benchmarks/scripts/validation.py run \
        -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
        -n ecoli --threads 1

    # Quick run with single n_points
    ./benchmarks/scripts/validation.py run \
        -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
        -n ecoli --n-points 100 --threads 1

    # Re-analyze existing results
    ./benchmarks/scripts/validation.py compare \
        -d benchmarks/results/validation/ecoli/sr

Output:
    benchmarks/results/validation/<name>/<algorithm>/
    ├── config.json
    ├── results_100.csv
    ├── results_200.csv
    ├── results_500.csv
    ├── results_1000.csv
    ├── results_lahuta_128.csv    # if lahuta available
    ├── validation_grid.png       # 4x4 grid (rows=variants, cols=n_points)
    ├── validation_zsasa_f64.png  # per-tool 1xN scatter
    ├── validation_zsasa_f32.png
    ├── validation_zsasa_bitmask_f64.png
    ├── validation_zsasa_bitmask_f32.png
    ├── validation_quicklook.png  # zsasa_f64 n_points=100 only
    └── validation_lahuta.png     # lahuta_bitmask n_points=128
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Annotated

import typer

if TYPE_CHECKING:
    import matplotlib.axes
    import polars as pl
from rich.console import Console

from bench_common import (
    LAHUTA_BITMASK_POINTS,
    get_binary_path,
    get_system_info,
)

app = typer.Typer(help="SASA validation: compare accuracy across tools")
console = Console()

# Default n_points levels for main validation
DEFAULT_N_POINTS = "100,200,500,1000"

# Lahuta bitmask uses a fixed n_points (must be in LAHUTA_BITMASK_POINTS)
LAHUTA_N_POINTS = 128

# zsasa variants to compare (column name, precision, use_bitmask)
ZSASA_VARIANTS = [
    ("zsasa_f64", "f64", False),
    ("zsasa_f32", "f32", False),
    ("zsasa_bitmask_f64", "f64", True),
    ("zsasa_bitmask_f32", "f32", True),
]


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def _root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def _freesasa_batch_path() -> Path:
    """Get freesasa_batch binary path (symlinked in external/bin/)."""
    return _root_dir().joinpath("benchmarks", "external", "bin", "freesasa_batch")


def _count_pdb_files(input_dir: Path) -> int:
    """Count PDB/ENT files in a directory."""
    return sum(
        1 for f in input_dir.iterdir() if f.is_file() and f.suffix in {".pdb", ".ent"}
    )


# ---------------------------------------------------------------------------
# Tool runners
# ---------------------------------------------------------------------------


def run_zsasa(
    input_dir: Path,
    algorithm: str,
    precision: str,
    threads: int,
    n_points: int,
    use_bitmask: bool = False,
) -> dict[str, tuple[float, int]]:
    """Run zsasa in batch mode. Returns {stem: (total_sasa, n_atoms)}."""
    zsasa = get_binary_path("zig")
    if not zsasa.exists():
        console.print(f"[yellow][SKIP] zsasa not found: {zsasa}[/]")
        return {}

    results: dict[str, tuple[float, int]] = {}
    bitmask_label = " bitmask" if use_bitmask else ""

    with tempfile.TemporaryDirectory(prefix="validation_zsasa_") as tmp:
        out_dir = Path(tmp)
        cmd = [
            str(zsasa),
            "batch",
            str(input_dir),
            str(out_dir),
            f"--algorithm={algorithm}",
            f"--precision={precision}",
            f"--threads={threads}",
            f"--n-points={n_points}",
        ]
        if use_bitmask:
            cmd.append("--use-bitmask")
        console.print(
            f"  [dim]$ zsasa batch ({precision}{bitmask_label}, n_points={n_points})[/]"
        )
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        except subprocess.TimeoutExpired:
            console.print(f"[red]zsasa timed out ({precision}{bitmask_label})[/]")
            return {}
        if proc.returncode != 0:
            console.print(f"[red]zsasa failed: {proc.stderr[:500]}[/]")
            return {}

        for json_file in sorted(out_dir.glob("*.json")):
            try:
                data = json.loads(json_file.read_text())
                total = float(data["total_area"])
                n_atoms = len(data.get("atom_areas", []))
                results[json_file.stem] = (total, n_atoms)
            except (json.JSONDecodeError, KeyError, ValueError):
                continue

    return results


def run_freesasa(
    input_dir: Path,
    n_points: int,
    threads: int = 1,
) -> dict[str, float]:
    """Run FreeSASA batch binary. Returns {stem: total_sasa}."""
    binary = _freesasa_batch_path()
    if not binary.exists():
        console.print(f"[yellow][SKIP] freesasa_batch not found: {binary}[/]")
        return {}

    results: dict[str, float] = {}

    with tempfile.TemporaryDirectory(prefix="validation_freesasa_") as tmp:
        out_dir = Path(tmp)
        cmd = [
            str(binary),
            str(input_dir),
            str(out_dir),
            f"--n-threads={threads}",
            f"--n-points={n_points}",
        ]
        console.print(f"  [dim]$ freesasa_batch (n_points={n_points})[/]")
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        except subprocess.TimeoutExpired:
            console.print("[red]freesasa_batch timed out[/]")
            return {}
        if proc.returncode != 0:
            console.print(f"[red]freesasa_batch failed: {proc.stderr[:500]}[/]")
            return {}

        for txt_file in sorted(out_dir.glob("*.txt")):
            try:
                total = float(txt_file.read_text().strip())
                results[txt_file.stem] = total
            except ValueError:
                continue

    return results


def run_rustsasa(
    input_dir: Path,
    n_points: int,
    threads: int = 1,
) -> dict[str, float]:
    """Run RustSASA on a PDB directory. Returns {stem: total_sasa}.

    Output format: per-file JSON with {"Protein": {"global_total": float}}.
    """
    rustsasa = get_binary_path("rustsasa")
    if not rustsasa.exists():
        console.print(f"[yellow][SKIP] rustsasa not found: {rustsasa}[/]")
        return {}

    results: dict[str, float] = {}

    with tempfile.TemporaryDirectory(prefix="validation_rustsasa_") as tmp:
        out_dir = Path(tmp)
        cmd = [
            str(rustsasa),
            str(input_dir),
            str(out_dir),
            "--format",
            "json",
            "-t",
            str(threads),
            "-n",
            str(n_points),
            "--allow-vdw-fallback",
            "-o",
            "protein",
        ]
        console.print(f"  [dim]$ rustsasa (n_points={n_points})[/]")
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        except subprocess.TimeoutExpired:
            console.print("[red]rustsasa timed out[/]")
            return {}
        if proc.returncode != 0:
            console.print(f"[red]rustsasa failed: {proc.stderr[:500]}[/]")
            return {}

        for json_file in sorted(out_dir.glob("*.json")):
            try:
                data = json.loads(json_file.read_text())
                total = float(data["Protein"]["global_total"])
                results[json_file.stem] = total
            except (json.JSONDecodeError, KeyError, ValueError, TypeError):
                continue

    return results


def run_lahuta(
    input_dir: Path,
    n_points: int,
    threads: int = 1,
) -> dict[str, float]:
    """Run Lahuta bitmask on a PDB directory. Returns {stem: total_sasa}.

    Output format: JSONL with {"model": "<path>", "sasa": [per_atom_floats]}.
    Total SASA = sum(sasa).
    """
    if n_points not in LAHUTA_BITMASK_POINTS:
        console.print(
            f"[yellow][SKIP] lahuta_bitmask only supports "
            f"n_points in {sorted(LAHUTA_BITMASK_POINTS)}, got {n_points}[/]"
        )
        return {}

    lahuta = get_binary_path("lahuta")
    if not lahuta.exists():
        console.print(f"[yellow][SKIP] lahuta not found: {lahuta}[/]")
        return {}

    # Resolve to absolute paths since cwd changes to temp directory
    abs_input = input_dir.resolve()
    abs_lahuta = lahuta.resolve()

    results: dict[str, float] = {}

    with tempfile.TemporaryDirectory(prefix="validation_lahuta_") as tmp:
        out_file = Path(tmp).joinpath("sasa.jsonl")
        cmd = [
            str(abs_lahuta),
            "sasa-sr",
            "-d",
            str(abs_input),
            "--is_af2_model",
            "--points",
            str(n_points),
            "--use-bitmask",
            "-t",
            str(threads),
            "--progress",
            "0",
            "-o",
            str(out_file),
        ]
        console.print(f"  [dim]$ lahuta sasa-sr (bitmask, n_points={n_points})[/]")
        try:
            proc = subprocess.run(
                cmd, cwd=tmp, capture_output=True, text=True, timeout=3600
            )
        except subprocess.TimeoutExpired:
            console.print("[red]lahuta timed out[/]")
            return {}
        if proc.returncode != 0:
            console.print(f"[red]lahuta failed: {proc.stderr[:500]}[/]")
            return {}

        if not out_file.exists():
            console.print("[red]lahuta produced no output[/]")
            return {}

        for line in out_file.read_text().strip().splitlines():
            try:
                data = json.loads(line)
                model_path = Path(data["model"])
                stem = model_path.stem
                total = sum(data["sasa"])
                results[stem] = round(total, 2)
            except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                continue

    return results


# ---------------------------------------------------------------------------
# Data collection and CSV creation
# ---------------------------------------------------------------------------


def collect_results_for_npoints(
    input_dir: Path,
    algorithm: str,
    threads: int,
    n_points: int,
    include_rustsasa: bool = True,
) -> tuple[dict[str, dict[str, tuple[float, int]]], dict[str, float], dict[str, float]]:
    """Run all tools for a single n_points value.

    Returns (zsasa_runs, freesasa_results, rustsasa_results).
    """
    console.print(f"\n[bold]=== n_points = {n_points} ===[/]")

    # FreeSASA baseline
    console.print("[bold cyan]Running FreeSASA...[/]")
    freesasa_results = run_freesasa(input_dir, n_points, threads)
    console.print(f"  Got {len(freesasa_results)} results")

    # zsasa variants
    zsasa_runs: dict[str, dict[str, tuple[float, int]]] = {}
    for col_name, precision, use_bitmask in ZSASA_VARIANTS:
        console.print(f"[bold cyan]Running {col_name}...[/]")
        zsasa_runs[col_name] = run_zsasa(
            input_dir, algorithm, precision, threads, n_points, use_bitmask
        )
        console.print(f"  Got {len(zsasa_runs[col_name])} results")

    # RustSASA (reference only)
    rustsasa_results: dict[str, float] = {}
    if include_rustsasa:
        console.print("[bold cyan]Running RustSASA...[/]")
        rustsasa_results = run_rustsasa(input_dir, n_points, threads)
        console.print(f"  Got {len(rustsasa_results)} results")

    return zsasa_runs, freesasa_results, rustsasa_results


def build_csv(
    zsasa_runs: dict[str, dict[str, tuple[float, int]]],
    freesasa_results: dict[str, float],
    rustsasa_results: dict[str, float],
    csv_path: Path,
) -> None:
    """Build and save CSV from collected results."""
    import polars as pl

    all_stems: set[str] = set()
    for results in zsasa_runs.values():
        all_stems.update(results.keys())
    all_stems.update(freesasa_results.keys())
    all_stems.update(rustsasa_results.keys())

    if not all_stems:
        console.print("[red]No results collected[/]")
        return

    rows: list[dict] = []
    for stem in sorted(all_stems):
        row: dict = {"structure": stem}

        # n_atoms from any zsasa run
        for results in zsasa_runs.values():
            if stem in results:
                _, n_atoms = results[stem]
                row["n_atoms"] = n_atoms
                break

        # FreeSASA
        if stem in freesasa_results:
            row["freesasa"] = round(freesasa_results[stem], 2)

        # zsasa variants
        for col_name, results in zsasa_runs.items():
            if stem in results:
                sasa, _ = results[stem]
                row[col_name] = round(sasa, 2)

        # RustSASA
        if stem in rustsasa_results:
            row["rustsasa"] = round(rustsasa_results[stem], 2)

        rows.append(row)

    # Column order: structure, n_atoms, freesasa, zsasa variants, rustsasa
    columns = ["structure", "n_atoms", "freesasa"]
    columns.extend(zsasa_runs.keys())
    if rustsasa_results:
        columns.append("rustsasa")

    df_raw = pl.DataFrame(rows)
    available_cols = [c for c in columns if c in df_raw.columns]
    df = df_raw.select(available_cols)
    df.write_csv(csv_path)
    console.print(f"[green]Saved:[/] {csv_path} ({df.height} structures)")


def _build_lahuta_csv(
    zsasa_runs: dict[str, dict[str, tuple[float, int]]],
    freesasa_results: dict[str, float],
    rustsasa_results: dict[str, float],
    lahuta_results: dict[str, float],
    csv_path: Path,
) -> None:
    """Build CSV with all tools + lahuta_bitmask for comparison."""
    import polars as pl

    all_stems: set[str] = set()
    for results in zsasa_runs.values():
        all_stems.update(results.keys())
    all_stems.update(freesasa_results.keys())
    all_stems.update(rustsasa_results.keys())
    all_stems.update(lahuta_results.keys())

    if not all_stems:
        console.print("[red]No results collected for lahuta CSV[/]")
        return

    rows: list[dict] = []
    for stem in sorted(all_stems):
        row: dict = {"structure": stem}

        # n_atoms from any zsasa run
        for results in zsasa_runs.values():
            if stem in results:
                _, n_atoms = results[stem]
                row["n_atoms"] = n_atoms
                break

        if stem in freesasa_results:
            row["freesasa"] = round(freesasa_results[stem], 2)

        for col_name, results in zsasa_runs.items():
            if stem in results:
                sasa, _ = results[stem]
                row[col_name] = round(sasa, 2)

        if stem in rustsasa_results:
            row["rustsasa"] = round(rustsasa_results[stem], 2)

        if stem in lahuta_results:
            row["lahuta_bitmask"] = round(lahuta_results[stem], 2)

        rows.append(row)

    columns = ["structure", "n_atoms", "freesasa"]
    columns.extend(zsasa_runs.keys())
    if rustsasa_results:
        columns.append("rustsasa")
    if lahuta_results:
        columns.append("lahuta_bitmask")

    df_raw = pl.DataFrame(rows)
    available_cols = [c for c in columns if c in df_raw.columns]
    df = df_raw.select(available_cols)
    df.write_csv(csv_path)
    console.print(f"[green]Saved:[/] {csv_path} ({df.height} structures)")


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def compute_stats(x: list[float], y: list[float]) -> dict[str, float]:
    """Compute R^2, mean/max relative error between two lists."""
    import numpy as np

    xa = np.array(x)
    ya = np.array(y)

    # Filter out zeros to avoid division by zero
    mask = ya > 0
    xa = xa[mask]
    ya = ya[mask]

    if len(xa) < 2:
        return {"r_squared": 0.0, "mean_rel_error": 0.0, "max_rel_error": 0.0}

    correlation = np.corrcoef(xa, ya)[0, 1]
    if np.isnan(correlation):
        return {"r_squared": 0.0, "mean_rel_error": 0.0, "max_rel_error": 0.0}
    r_squared = float(correlation**2)

    rel_errors = np.abs(xa - ya) / ya * 100
    mean_rel_error = float(np.mean(rel_errors))
    max_rel_error = float(np.max(rel_errors))

    return {
        "r_squared": r_squared,
        "mean_rel_error": mean_rel_error,
        "max_rel_error": max_rel_error,
    }


def _stats_for_pair(
    df: pl.DataFrame, reference: str, col: str
) -> dict[str, float] | None:
    """Compute stats for a pair of columns, or None if insufficient data."""
    pair = df.select([reference, col]).drop_nulls()
    if pair.height < 2:
        return None
    return compute_stats(pair[col].to_list(), pair[reference].to_list())


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def _setup_style() -> None:
    """Set up matplotlib style."""
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def _scatter_cell(
    ax: matplotlib.axes.Axes,
    df: pl.DataFrame,
    reference: str,
    col: str,
    title: str,
    ref_r2: dict[str, float] | None = None,
) -> None:
    """Draw a single scatter cell: col (x) vs reference (y).

    ref_r2: optional dict of {tool_name: r_squared} for reference tools,
    displayed in a separate box at bottom-right.
    """
    import numpy as np

    pair = df.select([reference, col]).drop_nulls()
    if pair.height < 2:
        ax.set_title(f"{title}: insufficient data")
        return

    ref_arr = pair[reference].to_numpy()
    comp_arr = pair[col].to_numpy()

    ax.scatter(comp_arr, ref_arr, alpha=0.3, s=8, color="#3498db", edgecolors="none")

    max_val = max(float(np.max(ref_arr)), float(np.max(comp_arr)))
    ax.plot([0, max_val], [0, max_val], "r--", linewidth=1, alpha=0.7)

    stats = compute_stats(comp_arr.tolist(), ref_arr.tolist())

    ax.set_xlabel(f"{col} SASA")
    ax.set_ylabel(f"{reference} SASA")
    ax.set_title(title, fontsize=10)

    # Main stats box (top-left)
    stats_lines = [
        f"R² = {stats['r_squared']:.6f}",
        f"Mean err = {stats['mean_rel_error']:.4f}%",
        f"Max err = {stats['max_rel_error']:.4f}%",
    ]
    ax.text(
        0.05,
        0.95,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    # Reference R² box (bottom-right, separate)
    if ref_r2:
        ref_lines = [f"{name} R²={r2:.6f}" for name, r2 in ref_r2.items()]
        ax.text(
            0.95,
            0.05,
            f"vs {reference}:\n" + "\n".join(ref_lines),
            transform=ax.transAxes,
            fontsize=7,
            verticalalignment="bottom",
            horizontalalignment="right",
            bbox=dict(boxstyle="round", facecolor="lightskyblue", alpha=0.4),
        )

    ax.set_aspect("equal")


def _collect_ref_r2(
    df: pl.DataFrame, reference: str, ref_tools: list[str]
) -> dict[str, float]:
    """Collect R² values for reference tools against the baseline."""
    ref_r2: dict[str, float] = {}
    for tool in ref_tools:
        if tool in df.columns and reference in df.columns:
            stats = _stats_for_pair(df, reference, tool)
            if stats:
                ref_r2[tool] = stats["r_squared"]
    return ref_r2


def _load_csvs(results_dir: Path, n_points_list: list[int]) -> dict[int, pl.DataFrame]:
    """Load results CSVs keyed by n_points. Returns only existing ones."""
    import polars as pl

    csvs: dict[int, pl.DataFrame] = {}
    for n_pts in n_points_list:
        csv_path = results_dir.joinpath(f"results_{n_pts}.csv")
        if csv_path.exists():
            csvs[n_pts] = pl.read_csv(csv_path)
    return csvs


def generate_grid_plot(
    results_dir: Path,
    n_points_list: list[int],
    algorithm: str,
) -> None:
    """Generate grid plot: rows=zsasa variants, cols=n_points.

    Each cell shows scatter vs FreeSASA with R^2/error stats.
    RustSASA R^2 shown as annotation in each cell.
    """
    import matplotlib.pyplot as plt

    _setup_style()

    csvs = _load_csvs(results_dir, n_points_list)
    if not csvs:
        console.print("[yellow]No results CSVs found, skipping grid plot[/]")
        return

    available_npts = sorted(csvs.keys())
    n_rows = len(ZSASA_VARIANTS)
    n_cols = len(available_npts)

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows), squeeze=False
    )

    for col_idx, n_pts in enumerate(available_npts):
        df = csvs[n_pts]
        ref_r2 = _collect_ref_r2(df, "freesasa", ["rustsasa"])

        for row_idx, (col_name, _precision, _bitmask) in enumerate(ZSASA_VARIANTS):
            ax = axes[row_idx][col_idx]
            if col_name not in df.columns or "freesasa" not in df.columns:
                ax.set_title(f"{n_pts} points: {col_name} missing")
                continue

            title = f"{col_name} ({n_pts} points)"
            _scatter_cell(ax, df, "freesasa", col_name, title, ref_r2=ref_r2)

    fig.suptitle(
        f"{algorithm.upper()} Validation: zsasa variants vs FreeSASA",
        fontsize=14,
        y=1.01,
    )
    fig.tight_layout()
    out_path = results_dir.joinpath("validation_grid.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {out_path}")


def generate_per_tool_plots(
    results_dir: Path,
    n_points_list: list[int],
    algorithm: str,
) -> None:
    """Generate 1xN scatter plot per zsasa variant across n_points levels."""
    import matplotlib.pyplot as plt

    _setup_style()

    csvs = _load_csvs(results_dir, n_points_list)
    if not csvs:
        return

    available_npts = sorted(csvs.keys())

    for col_name, _precision, _bitmask in ZSASA_VARIANTS:
        n_cols = len(available_npts)
        fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 5), squeeze=False)

        for col_idx, n_pts in enumerate(available_npts):
            ax = axes[0][col_idx]
            df = csvs[n_pts]

            if col_name not in df.columns or "freesasa" not in df.columns:
                ax.set_title(f"{n_pts} points: no data")
                continue

            ref_r2 = _collect_ref_r2(df, "freesasa", ["rustsasa"])
            title = f"{n_pts} points"
            _scatter_cell(ax, df, "freesasa", col_name, title, ref_r2=ref_r2)

        fig.suptitle(
            f"{algorithm.upper()}: {col_name} vs FreeSASA",
            fontsize=13,
            y=1.02,
        )
        fig.tight_layout()
        out_path = results_dir.joinpath(f"validation_{col_name}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        console.print(f"[green]Saved:[/] {out_path}")


def generate_quicklook_plot(
    results_dir: Path,
    algorithm: str,
) -> None:
    """Generate single scatter plot: zsasa_f64 at n_points=100 vs FreeSASA."""
    import matplotlib.pyplot as plt
    import polars as pl

    _setup_style()

    csv_path = results_dir.joinpath("results_100.csv")
    if not csv_path.exists():
        console.print("[yellow]results_100.csv not found, skipping quicklook plot[/]")
        return

    df = pl.read_csv(csv_path)
    if "zsasa_f64" not in df.columns or "freesasa" not in df.columns:
        console.print("[yellow]Missing zsasa_f64 or freesasa columns for quicklook[/]")
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    ref_r2 = _collect_ref_r2(df, "freesasa", ["rustsasa"])
    _scatter_cell(
        ax,
        df,
        "freesasa",
        "zsasa_f64",
        f"{algorithm.upper()}: zsasa (f64) vs FreeSASA (100 points)",
        ref_r2=ref_r2,
    )

    fig.tight_layout()
    out_path = results_dir.joinpath("validation_quicklook.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {out_path}")


def generate_lahuta_plot(
    results_dir: Path,
    algorithm: str,
) -> None:
    """Generate 1x4 scatter plot: zsasa variants vs FreeSASA at LAHUTA_N_POINTS.

    Each cell shows a zsasa variant with lahuta_bitmask R² as reference annotation.
    """
    import matplotlib.pyplot as plt
    import polars as pl

    _setup_style()

    csv_path = results_dir.joinpath(f"results_lahuta_{LAHUTA_N_POINTS}.csv")
    if not csv_path.exists():
        console.print(
            f"[yellow]results_lahuta_{LAHUTA_N_POINTS}.csv not found, "
            f"skipping lahuta plot[/]"
        )
        return

    df = pl.read_csv(csv_path)
    if "freesasa" not in df.columns:
        console.print("[yellow]Missing freesasa column for lahuta plot[/]")
        return

    n_cols = len(ZSASA_VARIANTS)
    fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 5), squeeze=False)

    ref_r2 = _collect_ref_r2(df, "freesasa", ["lahuta_bitmask", "rustsasa"])

    for col_idx, (col_name, _precision, _bitmask) in enumerate(ZSASA_VARIANTS):
        ax = axes[0][col_idx]
        if col_name not in df.columns:
            ax.set_title(f"{col_name}: no data")
            continue

        title = f"{col_name} ({LAHUTA_N_POINTS} points)"
        _scatter_cell(ax, df, "freesasa", col_name, title, ref_r2=ref_r2)

    fig.suptitle(
        f"{algorithm.upper()} Validation ({LAHUTA_N_POINTS} points): "
        f"zsasa variants vs FreeSASA",
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()
    out_path = results_dir.joinpath("validation_lahuta.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {out_path}")


def print_stats_table(
    results_dir: Path,
    n_points_list: list[int],
    reference: str = "freesasa",
) -> None:
    """Print statistics table for all n_points levels."""
    import polars as pl
    from rich.table import Table

    table = Table(title=f"Validation Statistics (reference: {reference})")
    table.add_column("n_points", style="cyan", justify="right")
    table.add_column("Tool", style="green")
    table.add_column("N", justify="right")
    table.add_column("R\u00b2", justify="right")
    table.add_column("Mean Error %", justify="right")
    table.add_column("Max Error %", justify="right")

    for n_pts in n_points_list:
        csv_path = results_dir.joinpath(f"results_{n_pts}.csv")
        if not csv_path.exists():
            continue

        df = pl.read_csv(csv_path)
        if reference not in df.columns:
            continue

        skip_cols = {"structure", "n_atoms", reference}
        compare_cols = [c for c in df.columns if c not in skip_cols]

        for col in compare_cols:
            stats = _stats_for_pair(df, reference, col)
            if stats is None:
                continue
            pair = df.select([reference, col]).drop_nulls()
            table.add_row(
                str(n_pts),
                col,
                str(pair.height),
                f"{stats['r_squared']:.6f}",
                f"{stats['mean_rel_error']:.4f}",
                f"{stats['max_rel_error']:.4f}",
            )

    # Lahuta CSV stats (includes all tools at LAHUTA_N_POINTS)
    lahuta_csv = results_dir.joinpath(f"results_lahuta_{LAHUTA_N_POINTS}.csv")
    if lahuta_csv.exists():
        df = pl.read_csv(lahuta_csv)
        if reference in df.columns:
            skip_cols = {"structure", "n_atoms", reference}
            compare_cols = [c for c in df.columns if c not in skip_cols]
            for col in compare_cols:
                stats = _stats_for_pair(df, reference, col)
                if stats is None:
                    continue
                pair = df.select([reference, col]).drop_nulls()
                table.add_row(
                    str(LAHUTA_N_POINTS),
                    col,
                    str(pair.height),
                    f"{stats['r_squared']:.6f}",
                    f"{stats['mean_rel_error']:.4f}",
                    f"{stats['max_rel_error']:.4f}",
                )

    console.print()
    console.print(table)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


@app.command()
def run(
    input_dir: Annotated[
        Path,
        typer.Option(
            "--input",
            "-i",
            help="Input directory containing PDB files",
            exists=True,
            file_okay=False,
            dir_okay=True,
        ),
    ],
    name: Annotated[
        str,
        typer.Option(
            "--name",
            "-n",
            help="Dataset name (used for output directory)",
        ),
    ],
    algorithm: Annotated[
        str,
        typer.Option(
            "--algorithm",
            "-a",
            help="Algorithm: sr (shrake-rupley) or lr (lee-richards)",
        ),
    ] = "sr",
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            "-T",
            help="Number of threads (default: all CPU cores)",
        ),
    ] = os.cpu_count() or 1,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output directory (default: benchmarks/results/validation/<name>/<algo>)",
        ),
    ] = None,
    n_points_str: Annotated[
        str,
        typer.Option(
            "--n-points",
            "-N",
            help="Comma-separated n_points values (default: 100,200,500,1000)",
        ),
    ] = DEFAULT_N_POINTS,
    skip_rustsasa: Annotated[
        bool,
        typer.Option(
            "--skip-rustsasa",
            help="Skip RustSASA (reference) runs",
        ),
    ] = False,
    skip_lahuta: Annotated[
        bool,
        typer.Option(
            "--skip-lahuta",
            help="Skip lahuta_bitmask run",
        ),
    ] = False,
) -> None:
    """Run all tools on PDB directory and compare SASA values."""
    try:
        n_points_list = sorted(int(x.strip()) for x in n_points_str.split(","))
    except ValueError:
        console.print(
            f"[red]Invalid --n-points: '{n_points_str}'. "
            f"Expected comma-separated integers (e.g. 100,200,500,1000)[/]"
        )
        raise typer.Exit(1) from None
    if any(n <= 0 for n in n_points_list):
        console.print("[red]All n_points values must be positive[/]")
        raise typer.Exit(1)

    # Set up output
    if output_dir is None:
        results_dir = _root_dir().joinpath(
            "benchmarks", "results", "validation", name, algorithm
        )
    else:
        results_dir = output_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    n_files = _count_pdb_files(input_dir)

    # Save config
    config = {
        "timestamp": datetime.now().strftime("%Y-%m-%d_%H%M%S"),
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "input_dir": str(input_dir),
            "n_files": n_files,
            "algorithm": algorithm,
            "threads": threads,
            "n_points": n_points_list,
            "skip_rustsasa": skip_rustsasa,
            "skip_lahuta": skip_lahuta,
            "lahuta_n_points": LAHUTA_N_POINTS,
        },
    }
    results_dir.joinpath("config.json").write_text(json.dumps(config, indent=2))

    # Print header
    console.print(f"\n[bold]=== SASA Validation: {name} ===[/]")
    console.print(f"Input: {input_dir} ({n_files} PDB files)")
    console.print(f"Output: {results_dir}")
    console.print(f"Algorithm: {algorithm}, Threads: {threads}")
    console.print(f"n_points: {n_points_list}")
    zsasa_labels = [v[0] for v in ZSASA_VARIANTS]
    tools_str = f"freesasa, {', '.join(zsasa_labels)}"
    if not skip_rustsasa:
        tools_str += ", rustsasa (ref)"
    if not skip_lahuta:
        tools_str += f", lahuta_bitmask (n={LAHUTA_N_POINTS})"
    console.print(f"Tools: {tools_str}")

    # --- Main runs: each n_points level ---
    for n_pts in n_points_list:
        zsasa_runs, freesasa_results, rustsasa_results = collect_results_for_npoints(
            input_dir, algorithm, threads, n_pts, include_rustsasa=not skip_rustsasa
        )
        csv_path = results_dir.joinpath(f"results_{n_pts}.csv")
        build_csv(zsasa_runs, freesasa_results, rustsasa_results, csv_path)

    # --- Lahuta bitmask run (separate n_points, includes all tools for comparison) ---
    if not skip_lahuta:
        console.print(
            f"\n[bold]=== Lahuta bitmask + all tools "
            f"(n_points={LAHUTA_N_POINTS}) ===[/]"
        )

        # Run all standard tools at LAHUTA_N_POINTS for apples-to-apples comparison
        zsasa_runs, freesasa_results, rustsasa_results = collect_results_for_npoints(
            input_dir,
            algorithm,
            threads,
            LAHUTA_N_POINTS,
            include_rustsasa=not skip_rustsasa,
        )

        # Run lahuta_bitmask
        console.print("[bold cyan]Running lahuta_bitmask...[/]")
        lahuta_results = run_lahuta(input_dir, LAHUTA_N_POINTS, threads)
        console.print(f"  Got {len(lahuta_results)} results")

        # Build CSV with all tools + lahuta_bitmask
        lahuta_csv = results_dir.joinpath(f"results_lahuta_{LAHUTA_N_POINTS}.csv")
        _build_lahuta_csv(
            zsasa_runs, freesasa_results, rustsasa_results, lahuta_results, lahuta_csv
        )

    # --- Statistics and plots ---
    print_stats_table(results_dir, n_points_list)
    generate_grid_plot(results_dir, n_points_list, algorithm)
    generate_per_tool_plots(results_dir, n_points_list, algorithm)
    generate_quicklook_plot(results_dir, algorithm)
    generate_lahuta_plot(results_dir, algorithm)

    console.print(f"\n[bold green]=== Done! Results: {results_dir} ===[/]")


@app.command()
def compare(
    results_dir: Annotated[
        Path,
        typer.Option(
            "--dir",
            "-d",
            help="Directory containing results CSV files",
            exists=True,
            file_okay=False,
            dir_okay=True,
        ),
    ],
) -> None:
    """Re-analyze existing validation results (regenerate plots and stats)."""
    # Detect algorithm from config.json
    algorithm = "sr"
    config_path = results_dir.joinpath("config.json")
    if config_path.exists():
        try:
            config = json.loads(config_path.read_text())
            algorithm = config.get("parameters", {}).get("algorithm", "sr")
        except (json.JSONDecodeError, KeyError):
            pass

    # Discover available n_points CSVs
    import re

    n_points_list: list[int] = []
    for csv_file in sorted(results_dir.glob("results_*.csv")):
        match = re.match(r"results_(\d+)\.csv", csv_file.name)
        if match:
            n_points_list.append(int(match.group(1)))

    if not n_points_list:
        console.print(f"[red]No results_*.csv files found in {results_dir}[/]")
        raise typer.Exit(1)

    console.print(f"[bold]=== Re-analyzing: {results_dir} ===[/]")
    console.print(f"Algorithm: {algorithm}")
    console.print(f"Found n_points: {n_points_list}")

    print_stats_table(results_dir, n_points_list)
    generate_grid_plot(results_dir, n_points_list, algorithm)
    generate_per_tool_plots(results_dir, n_points_list, algorithm)
    generate_quicklook_plot(results_dir, algorithm)
    generate_lahuta_plot(results_dir, algorithm)

    console.print("\n[bold green]=== Done! ===[/]")


if __name__ == "__main__":
    app()

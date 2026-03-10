#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "cffi",
#     "matplotlib>=3.8",
#     "MDAnalysis",
#     "mdtraj",
#     "numpy>=1.26",
#     "polars>=1.0",
#     "rich>=13.0",
#     "typer>=0.9.0",
# ]
# ///
"""MD trajectory SASA validation: compare per-frame accuracy across tools.

Loads an XTC trajectory and compares per-frame total SASA values from:
- mdtraj native (md.shrake_rupley, reference)
- zsasa_mdtraj (zsasa Python wrapper via mdtraj)
- zsasa_mdanalysis (zsasa Python wrapper via MDAnalysis)
- zsasa CLI traj mode (Zig native XTC reader, f32/f64/bitmask variants)

This validates that XTC reading and coordinate handling produce
consistent SASA values across implementations.

Usage:
    # Run all tools (default n_points: 100,200,500,1000)
    ./benchmarks/scripts/validation_md.py run \
        --xtc benchmarks/md_data/6sup_A_protein/6sup_A_prod_R1_fit.xtc \
        --pdb benchmarks/md_data/6sup_A_protein/6sup_A.pdb \
        -n 6sup_R1

    # Custom n_points and stride
    ./benchmarks/scripts/validation_md.py run \
        --xtc benchmarks/md_data/6sup_A_protein/6sup_A_prod_R1_fit.xtc \
        --pdb benchmarks/md_data/6sup_A_protein/6sup_A.pdb \
        -n 6sup_R1 --n-points 100,500 --stride 10

    # Re-analyze existing results
    ./benchmarks/scripts/validation_md.py compare \
        -d benchmarks/results/validation_md/6sup_R1

Output:
    benchmarks/results/validation_md/<name>/
    ├── config.json
    ├── results_100.csv
    ├── results_200.csv
    ├── results_500.csv
    ├── results_1000.csv
    ├── validation_grid.png
    ├── validation_zsasa_mdtraj.png
    ├── validation_zsasa_mdanalysis.png
    ├── validation_zsasa_cli_f64.png
    ├── validation_zsasa_cli_f32.png
    ├── validation_zsasa_cli_bitmask_f64.png
    ├── validation_zsasa_cli_bitmask_f32.png
    ├── validation_xtc_comparison.png
    └── validation_bitmask_comparison.png
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Annotated

import numpy as np
import typer

if TYPE_CHECKING:
    import matplotlib.axes
    import polars as pl

from rich.console import Console

from bench_common import get_binary_path, get_system_info

# Add zsasa Python package to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.joinpath("python")))

app = typer.Typer(help="MD trajectory SASA validation: compare per-frame accuracy")
console = Console()

# Default n_points levels (same as validation.py)
DEFAULT_N_POINTS = "100,200,500,1000"

# zsasa MD variants: (column_name, source, precision, use_bitmask)
# source: "mdtraj" = Python mdtraj wrapper, "mdanalysis" = Python MDAnalysis wrapper,
#         "cli" = CLI traj mode
ZSASA_MD_VARIANTS: list[tuple[str, str, str, bool]] = [
    ("zsasa_mdtraj", "mdtraj", "", False),
    ("zsasa_mdanalysis", "mdanalysis", "", False),
    ("zsasa_cli_f64", "cli", "f64", False),
    ("zsasa_cli_f32", "cli", "f32", False),
    ("zsasa_cli_bitmask_f64", "cli", "f64", True),
    ("zsasa_cli_bitmask_f32", "cli", "f32", True),
]


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------


def _root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


# ---------------------------------------------------------------------------
# Tool runners (all return per-frame total SASA in Å²)
# ---------------------------------------------------------------------------


def run_mdtraj_native(
    xtc: Path, pdb: Path, n_points: int, stride: int
) -> dict[int, float]:
    """Run mdtraj.shrake_rupley. Returns {frame_idx: total_sasa_angstrom2}."""
    import mdtraj as md

    console.print("  Loading trajectory with mdtraj...")
    traj = md.load(str(xtc), top=str(pdb), stride=stride)
    console.print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms")

    console.print("  Computing SASA with md.shrake_rupley...")
    sasa = md.shrake_rupley(traj, n_sphere_points=n_points)
    # sasa shape: (n_frames, n_atoms) in nm², convert to Å²
    totals = sasa.sum(axis=1) * 100.0

    return {i: round(float(v), 2) for i, v in enumerate(totals)}


def run_zsasa_mdtraj(
    xtc: Path, pdb: Path, n_points: int, stride: int, threads: int
) -> dict[int, float]:
    """Run zsasa via mdtraj wrapper. Returns {frame_idx: total_sasa_angstrom2}."""
    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    console.print("  Loading trajectory with mdtraj...")
    traj = md.load(str(xtc), top=str(pdb), stride=stride)
    console.print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms")

    console.print(f"  Computing SASA with zsasa.mdtraj (threads={threads})...")
    # Returns nm², convert to Å²
    totals = compute_sasa(traj, n_points=n_points, n_threads=threads, mode="total")
    totals_a2 = totals * 100.0

    return {i: round(float(v), 2) for i, v in enumerate(totals_a2)}


def run_zsasa_mdanalysis(
    xtc: Path, pdb: Path, n_points: int, stride: int, threads: int
) -> dict[int, float]:
    """Run zsasa via MDAnalysis wrapper. Returns {frame_idx: total_sasa_angstrom2}."""
    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    console.print("  Loading trajectory with MDAnalysis...")
    u = mda.Universe(str(pdb), str(xtc))
    console.print(f"  {len(u.trajectory)} frames, {u.atoms.n_atoms} atoms")

    console.print(f"  Computing SASA with zsasa.mdanalysis (threads={threads})...")
    analysis = SASAAnalysis(u)
    analysis.run(step=stride, n_points=n_points, n_threads=threads)
    # Returns Å² (no conversion needed)
    totals = analysis.results.total_area

    return {i: round(float(v), 2) for i, v in enumerate(totals)}


def run_zsasa_cli(
    xtc: Path,
    pdb: Path,
    n_points: int,
    stride: int,
    threads: int,
    precision: str = "f64",
    use_bitmask: bool = False,
) -> dict[int, float]:
    """Run zsasa CLI traj mode. Returns {frame_idx: total_sasa_angstrom2}."""
    zsasa = get_binary_path("zig")
    if not zsasa.exists():
        console.print(f"[yellow][SKIP] zsasa not found: {zsasa}[/]")
        return {}

    bitmask_label = " bitmask" if use_bitmask else ""

    with tempfile.NamedTemporaryFile(suffix=".csv", prefix="zsasa_val_") as tmp:
        out_path = tmp.name
        cmd = [
            str(zsasa),
            "traj",
            str(xtc),
            str(pdb),
            "--include-hydrogens",
            f"--threads={threads}",
            f"--n-points={n_points}",
            f"--stride={stride}",
            f"--precision={precision}",
        ]
        if use_bitmask:
            cmd.append("--use-bitmask")
        cmd.extend(["-o", out_path, "-q"])
        console.print(
            f"  [dim]$ zsasa traj ({precision}{bitmask_label}, n_points={n_points})[/]"
        )

        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        except subprocess.TimeoutExpired:
            console.print(f"[red]zsasa traj timed out ({precision}{bitmask_label})[/]")
            return {}
        if proc.returncode != 0:
            console.print(f"[red]zsasa traj failed: {proc.stderr[:500]}[/]")
            return {}

        # Parse CSV: frame,time,total_sasa (Å²)
        results: dict[int, float] = {}
        with open(out_path) as f:
            reader = csv.DictReader(f)
            for i, row in enumerate(reader):
                results[i] = float(row["total_sasa"])

    return results


# ---------------------------------------------------------------------------
# Data collection and CSV creation
# ---------------------------------------------------------------------------


def collect_results_for_npoints(
    xtc: Path,
    pdb: Path,
    stride: int,
    threads: int,
    n_points: int,
) -> tuple[dict[int, float], dict[str, dict[int, float]]]:
    """Run all tools for a single n_points value.

    Returns (mdtraj_results, zsasa_runs).
    """
    console.print(f"\n[bold]=== n_points = {n_points} ===[/]")

    # mdtraj native (reference)
    console.print("[bold cyan]Running mdtraj native...[/]")
    mdtraj_results = run_mdtraj_native(xtc, pdb, n_points, stride)
    console.print(f"  Got {len(mdtraj_results)} frames")

    # zsasa variants
    zsasa_runs: dict[str, dict[int, float]] = {}
    for col_name, source, precision, use_bitmask in ZSASA_MD_VARIANTS:
        console.print(f"[bold cyan]Running {col_name}...[/]")
        if source == "mdtraj":
            zsasa_runs[col_name] = run_zsasa_mdtraj(xtc, pdb, n_points, stride, threads)
        elif source == "mdanalysis":
            zsasa_runs[col_name] = run_zsasa_mdanalysis(
                xtc, pdb, n_points, stride, threads
            )
        elif source == "cli":
            zsasa_runs[col_name] = run_zsasa_cli(
                xtc, pdb, n_points, stride, threads, precision, use_bitmask
            )
        console.print(f"  Got {len(zsasa_runs[col_name])} frames")

    return mdtraj_results, zsasa_runs


def build_csv(
    mdtraj_results: dict[int, float],
    zsasa_runs: dict[str, dict[int, float]],
    csv_path: Path,
) -> None:
    """Build and save CSV from collected results."""
    import polars as pl

    # Use frame indices from mdtraj as the master set
    all_frames = sorted(mdtraj_results.keys())
    if not all_frames:
        console.print("[red]No results collected[/]")
        return

    rows: list[dict] = []
    for frame in all_frames:
        row: dict = {"frame": frame}
        row["mdtraj"] = mdtraj_results.get(frame)

        for col_name in zsasa_runs:
            row[col_name] = zsasa_runs[col_name].get(frame)

        rows.append(row)

    # Column order: frame, mdtraj, zsasa variants
    columns = ["frame", "mdtraj"]
    columns.extend(zsasa_runs.keys())

    df_raw = pl.DataFrame(rows)
    available_cols = [c for c in columns if c in df_raw.columns]
    df = df_raw.select(available_cols)
    df.write_csv(csv_path)
    console.print(f"[green]Saved:[/] {csv_path} ({df.height} frames)")


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def compute_stats(x: list[float], y: list[float]) -> dict[str, float]:
    """Compute R^2, mean/max relative error between two lists."""
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


def compute_delta_r2(x: list[float], y: list[float]) -> float:
    """Compute R² of frame-to-frame ΔSASA between two series.

    Measures how well the tool tracks relative changes in SASA,
    independent of absolute value offset.
    """
    xa = np.array(x)
    ya = np.array(y)

    if len(xa) < 3:
        return 0.0

    dx = np.diff(xa)
    dy = np.diff(ya)

    correlation = np.corrcoef(dx, dy)[0, 1]
    if np.isnan(correlation):
        return 0.0
    return float(correlation**2)


def _stats_for_pair(
    df: pl.DataFrame, reference: str, col: str
) -> dict[str, float] | None:
    """Compute stats for a pair of columns, or None if insufficient data."""
    pair = df.select([reference, col]).drop_nulls()
    if pair.height < 2:
        return None
    ref_list = pair[reference].to_list()
    col_list = pair[col].to_list()
    stats = compute_stats(col_list, ref_list)
    stats["delta_r2"] = compute_delta_r2(col_list, ref_list)
    return stats


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
    data_lim: tuple[float, float] | None = None,
) -> None:
    """Draw a single scatter cell: col (x) vs reference (y).

    ref_r2: optional dict of {tool_name: r_squared} for reference tools,
    displayed in a separate box at bottom-right.
    data_lim: optional (min, max) to set consistent axis limits across cells.
    If None, auto-zoom to data range with 5% margin.
    """
    pair = df.select([reference, col]).drop_nulls()
    if pair.height < 2:
        ax.set_title(f"{title}: insufficient data")
        return

    ref_arr = pair[reference].to_numpy()
    comp_arr = pair[col].to_numpy()

    ax.scatter(comp_arr, ref_arr, alpha=0.3, s=8, color="#3498db", edgecolors="none")

    # y=x line across the visible range
    if data_lim is not None:
        lo, hi = data_lim
    else:
        all_vals = np.concatenate([ref_arr, comp_arr])
        lo, hi = float(np.min(all_vals)), float(np.max(all_vals))
        margin = (hi - lo) * 0.05
        lo -= margin
        hi += margin
    ax.plot([lo, hi], [lo, hi], "r--", linewidth=1, alpha=0.7)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    comp_list = comp_arr.tolist()
    ref_list = ref_arr.tolist()
    stats = compute_stats(comp_list, ref_list)
    delta_r2 = compute_delta_r2(comp_list, ref_list)

    ax.set_xlabel(f"{col} SASA")
    ax.set_ylabel(f"{reference} SASA")
    ax.set_title(title, fontsize=10)

    # Main stats box (top-left)
    stats_lines = [
        f"R² = {stats['r_squared']:.6f}",
        f"ΔR² = {delta_r2:.6f}",
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


def _compute_global_lim(
    csvs: dict[int, pl.DataFrame],
    columns: list[str],
) -> tuple[float, float]:
    """Compute global (min, max) across all CSVs and columns with 5% margin."""
    all_vals: list[float] = []
    for df in csvs.values():
        for col in columns:
            if col in df.columns:
                vals = df[col].drop_nulls().to_list()
                all_vals.extend(vals)
    if not all_vals:
        return (0.0, 1.0)
    lo, hi = min(all_vals), max(all_vals)
    margin = (hi - lo) * 0.05
    return (lo - margin, hi + margin)


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
) -> None:
    """Generate grid plot: rows=zsasa variants, cols=n_points.

    Each cell shows scatter vs mdtraj with R²/error stats.
    """
    import matplotlib.pyplot as plt

    _setup_style()

    csvs = _load_csvs(results_dir, n_points_list)
    if not csvs:
        console.print("[yellow]No results CSVs found, skipping grid plot[/]")
        return

    available_npts = sorted(csvs.keys())
    n_rows = len(ZSASA_MD_VARIANTS)
    n_cols = len(available_npts)

    # Compute global axis limits across all cells
    all_columns = ["mdtraj"] + [v[0] for v in ZSASA_MD_VARIANTS]
    data_lim = _compute_global_lim(csvs, all_columns)

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows), squeeze=False
    )

    for col_idx, n_pts in enumerate(available_npts):
        df = csvs[n_pts]

        for row_idx, (col_name, _source, _prec, _bitmask) in enumerate(
            ZSASA_MD_VARIANTS
        ):
            ax = axes[row_idx][col_idx]
            if col_name not in df.columns or "mdtraj" not in df.columns:
                ax.set_title(f"{n_pts} points: {col_name} missing")
                continue

            title = f"{col_name} ({n_pts} points)"
            _scatter_cell(ax, df, "mdtraj", col_name, title, data_lim=data_lim)

    fig.suptitle(
        "MD Validation: zsasa variants vs mdtraj",
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
) -> None:
    """Generate 1xN scatter plot per zsasa variant across n_points levels."""
    import matplotlib.pyplot as plt

    _setup_style()

    csvs = _load_csvs(results_dir, n_points_list)
    if not csvs:
        return

    available_npts = sorted(csvs.keys())

    for col_name, _source, _prec, _bitmask in ZSASA_MD_VARIANTS:
        # Compute global limits for this tool across all n_points
        tool_lim = _compute_global_lim(csvs, ["mdtraj", col_name])

        n_cols = len(available_npts)
        fig, axes = plt.subplots(1, n_cols, figsize=(5 * n_cols, 5), squeeze=False)

        for col_idx, n_pts in enumerate(available_npts):
            ax = axes[0][col_idx]
            df = csvs[n_pts]

            if col_name not in df.columns or "mdtraj" not in df.columns:
                ax.set_title(f"{n_pts} points: no data")
                continue

            title = f"{n_pts} points"
            _scatter_cell(ax, df, "mdtraj", col_name, title, data_lim=tool_lim)

        fig.suptitle(
            f"MD: {col_name} vs mdtraj",
            fontsize=13,
            y=1.02,
        )
        fig.tight_layout()
        out_path = results_dir.joinpath(f"validation_{col_name}.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        console.print(f"[green]Saved:[/] {out_path}")


def generate_xtc_comparison_plot(
    results_dir: Path,
) -> None:
    """Generate scatter plot: zsasa_mdtraj vs zsasa_cli_f64 at n_points=100.

    This compares XTC decompression libraries (C via mdtraj vs Zig native).
    """
    import matplotlib.pyplot as plt
    import polars as pl

    _setup_style()

    csv_path = results_dir.joinpath("results_100.csv")
    if not csv_path.exists():
        console.print(
            "[yellow]results_100.csv not found, skipping XTC comparison plot[/]"
        )
        return

    df = pl.read_csv(csv_path)
    if "zsasa_mdtraj" not in df.columns or "zsasa_cli_f64" not in df.columns:
        console.print(
            "[yellow]Missing zsasa_mdtraj or zsasa_cli_f64 columns "
            "for XTC comparison[/]"
        )
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    _scatter_cell(
        ax,
        df,
        "zsasa_mdtraj",
        "zsasa_cli_f64",
        "XTC comparison: zsasa_cli_f64 vs zsasa_mdtraj (100 points)",
    )

    fig.tight_layout()
    out_path = results_dir.joinpath("validation_xtc_comparison.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {out_path}")


def generate_bitmask_comparison_plot(
    results_dir: Path,
) -> None:
    """Generate scatter plot: zsasa_cli_f64 vs zsasa_cli_bitmask_f64 at n_points=100.

    This isolates the bitmask LUT quantization error by comparing
    the same algorithm with and without bitmask optimization.
    """
    import matplotlib.pyplot as plt
    import polars as pl

    _setup_style()

    csv_path = results_dir.joinpath("results_100.csv")
    if not csv_path.exists():
        console.print(
            "[yellow]results_100.csv not found, skipping bitmask comparison plot[/]"
        )
        return

    df = pl.read_csv(csv_path)
    if "zsasa_cli_f64" not in df.columns or "zsasa_cli_bitmask_f64" not in df.columns:
        console.print(
            "[yellow]Missing zsasa_cli_f64 or zsasa_cli_bitmask_f64 columns "
            "for bitmask comparison[/]"
        )
        return

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    _scatter_cell(
        ax,
        df,
        "zsasa_cli_f64",
        "zsasa_cli_bitmask_f64",
        "Bitmask comparison: bitmask_f64 vs f64 (100 points)",
    )

    fig.tight_layout()
    out_path = results_dir.joinpath("validation_bitmask_comparison.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {out_path}")


def print_stats_table(
    results_dir: Path,
    n_points_list: list[int],
    reference: str = "mdtraj",
) -> None:
    """Print statistics table for all n_points levels."""
    import polars as pl
    from rich.table import Table

    table = Table(title=f"MD Validation Statistics (reference: {reference})")
    table.add_column("n_points", style="cyan", justify="right")
    table.add_column("Tool", style="green")
    table.add_column("N", justify="right")
    table.add_column("R\u00b2", justify="right")
    table.add_column("\u0394R\u00b2", justify="right")
    table.add_column("Mean Error %", justify="right")
    table.add_column("Max Error %", justify="right")

    for n_pts in n_points_list:
        csv_path = results_dir.joinpath(f"results_{n_pts}.csv")
        if not csv_path.exists():
            continue

        df = pl.read_csv(csv_path)
        if reference not in df.columns:
            continue

        skip_cols = {"frame", reference}
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
                f"{stats['delta_r2']:.6f}",
                f"{stats['mean_rel_error']:.4f}",
                f"{stats['max_rel_error']:.4f}",
            )

    # Cross-tool comparison stats at n_points=100
    csv_100 = results_dir.joinpath("results_100.csv")
    if csv_100.exists():
        df = pl.read_csv(csv_100)

        # XTC comparison: zsasa_mdtraj vs zsasa_cli_f64
        if "zsasa_mdtraj" in df.columns and "zsasa_cli_f64" in df.columns:
            stats = _stats_for_pair(df, "zsasa_mdtraj", "zsasa_cli_f64")
            if stats:
                pair = df.select(["zsasa_mdtraj", "zsasa_cli_f64"]).drop_nulls()
                table.add_row(
                    "100 (XTC)",
                    "zsasa_cli_f64 vs zsasa_mdtraj",
                    str(pair.height),
                    f"{stats['r_squared']:.6f}",
                    f"{stats['delta_r2']:.6f}",
                    f"{stats['mean_rel_error']:.4f}",
                    f"{stats['max_rel_error']:.4f}",
                )

        # Bitmask comparison: zsasa_cli_f64 vs zsasa_cli_bitmask_f64
        if "zsasa_cli_f64" in df.columns and "zsasa_cli_bitmask_f64" in df.columns:
            stats = _stats_for_pair(df, "zsasa_cli_f64", "zsasa_cli_bitmask_f64")
            if stats:
                pair = df.select(
                    ["zsasa_cli_f64", "zsasa_cli_bitmask_f64"]
                ).drop_nulls()
                table.add_row(
                    "100 (bitmask)",
                    "bitmask_f64 vs cli_f64",
                    str(pair.height),
                    f"{stats['r_squared']:.6f}",
                    f"{stats['delta_r2']:.6f}",
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
    xtc: Annotated[
        Path,
        typer.Option(
            "--xtc",
            help="XTC trajectory file",
            exists=True,
            file_okay=True,
            dir_okay=False,
        ),
    ],
    pdb: Annotated[
        Path,
        typer.Option(
            "--pdb",
            help="Topology PDB file",
            exists=True,
            file_okay=True,
            dir_okay=False,
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
    n_points_str: Annotated[
        str | None,
        typer.Option(
            "--n-points",
            "-N",
            help="Comma-separated n_points values (default: 100,200,500,1000)",
        ),
    ] = None,
    stride: Annotated[
        int,
        typer.Option(
            "--stride",
            "-s",
            help="Frame stride (process every Nth frame)",
        ),
    ] = 1,
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
            help="Output directory (default: benchmarks/results/validation_md/<name>)",
        ),
    ] = None,
) -> None:
    """Run all tools on MD trajectory and compare per-frame SASA values."""
    # Parse n_points
    if n_points_str is None:
        n_points_str = DEFAULT_N_POINTS
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
            "benchmarks", "results", "validation_md", name
        )
    else:
        results_dir = output_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    # Save config
    variant_names = [v[0] for v in ZSASA_MD_VARIANTS]
    config = {
        "timestamp": datetime.now().strftime("%Y-%m-%d_%H%M%S"),
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "xtc": str(xtc),
            "pdb": str(pdb),
            "tools": ["mdtraj"] + variant_names,
            "n_points": n_points_list,
            "stride": stride,
            "threads": threads,
        },
    }
    results_dir.joinpath("config.json").write_text(json.dumps(config, indent=2))

    # Print header
    console.print(f"\n[bold]=== MD SASA Validation: {name} ===[/]")
    console.print(f"XTC: {xtc}")
    console.print(f"PDB: {pdb}")
    console.print(f"Output: {results_dir}")
    console.print(f"Tools: mdtraj (ref), {', '.join(variant_names)}")
    console.print(f"n_points: {n_points_list}, Stride: {stride}, Threads: {threads}")

    # --- Main runs: each n_points level ---
    for n_pts in n_points_list:
        mdtraj_results, zsasa_runs = collect_results_for_npoints(
            xtc, pdb, stride, threads, n_pts
        )
        csv_path = results_dir.joinpath(f"results_{n_pts}.csv")
        build_csv(mdtraj_results, zsasa_runs, csv_path)

    # --- Statistics and plots ---
    print_stats_table(results_dir, n_points_list)
    generate_grid_plot(results_dir, n_points_list)
    generate_per_tool_plots(results_dir, n_points_list)
    generate_xtc_comparison_plot(results_dir)
    generate_bitmask_comparison_plot(results_dir)

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
    """Re-analyze existing MD validation results (regenerate plots and stats)."""
    import re

    # Discover available n_points CSVs
    n_points_list: list[int] = []
    for csv_file in sorted(results_dir.glob("results_*.csv")):
        match = re.match(r"results_(\d+)\.csv", csv_file.name)
        if match:
            n_points_list.append(int(match.group(1)))

    if not n_points_list:
        console.print(f"[red]No results_*.csv files found in {results_dir}[/]")
        raise typer.Exit(1)

    console.print(f"[bold]=== Re-analyzing: {results_dir} ===[/]")
    console.print(f"Found n_points: {n_points_list}")

    print_stats_table(results_dir, n_points_list)
    generate_grid_plot(results_dir, n_points_list)
    generate_per_tool_plots(results_dir, n_points_list)
    generate_xtc_comparison_plot(results_dir)
    generate_bitmask_comparison_plot(results_dir)

    console.print("\n[bold green]=== Done! ===[/]")


if __name__ == "__main__":
    app()

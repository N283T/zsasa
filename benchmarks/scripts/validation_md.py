#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "cffi",
#     "matplotlib>=3.8",
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
- zsasa CLI traj mode (Zig native XTC reader)

This validates that XTC reading and coordinate handling produce
consistent SASA values across implementations.

Usage:
    # Run all tools
    ./benchmarks/scripts/validation_md.py run \
        --xtc benchmarks/md_data/6sup_A_protein/6sup_A_prod_R1_fit.xtc \
        --pdb benchmarks/md_data/6sup_A_protein/6sup_A.pdb \
        -n 6sup_R1

    # Compare n_points convergence
    ./benchmarks/scripts/validation_md.py run \
        --xtc benchmarks/md_data/6sup_A_protein/6sup_A_prod_R1_fit.xtc \
        --pdb benchmarks/md_data/6sup_A_protein/6sup_A.pdb \
        -n 6sup_R1_npoints --n-points 100,500,960 --stride 100

    # Re-analyze existing CSV
    ./benchmarks/scripts/validation_md.py compare \
        -d benchmarks/results/validation_md/6sup_R1

Output:
    benchmarks/results/validation_md/<name>/
    ├── config.json              # System info, parameters
    ├── results.csv              # Per-frame SASA values
    ├── validation_md.png        # Scatter plot (single n_points)
    └── validation_md_<N>.png    # Per-n_points scatter plots (multi n_points)
"""

from __future__ import annotations

import csv
import json
import os
import platform
import subprocess
import sys
import tempfile
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Annotated

import numpy as np
import typer
from rich.console import Console
from rich.table import Table

# Add zsasa Python package to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.joinpath("python")))

app = typer.Typer(help="MD trajectory SASA validation: compare per-frame accuracy")
console = Console()


# ---------------------------------------------------------------------------
# Shared utilities
# ---------------------------------------------------------------------------


class MdTool(str, Enum):
    mdtraj = "mdtraj"
    zsasa_mdtraj = "zsasa_mdtraj"
    zsasa_cli = "zsasa_cli"


ALL_MD_TOOLS = [MdTool.mdtraj, MdTool.zsasa_mdtraj, MdTool.zsasa_cli]


def parse_n_points(n_points_str: str) -> list[int]:
    """Parse n_points specification like '100,500,960'."""
    return sorted(int(x.strip()) for x in n_points_str.split(","))


def _detect_reference(results_dir: Path) -> str:
    """Detect reference tool base name from config.json."""
    config_path = results_dir.joinpath("config.json")
    if config_path.exists():
        try:
            config = json.loads(config_path.read_text())
            return config.get("parameters", {}).get("reference", "mdtraj")
        except (json.JSONDecodeError, KeyError):
            pass
    return "mdtraj"


def get_root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def get_system_info() -> dict:
    """Get system information."""
    info = {
        "os": platform.system(),
        "os_version": platform.release(),
        "arch": platform.machine(),
        "cpu_cores": os.cpu_count() or 1,
    }

    if platform.system() == "Darwin":
        try:
            result = subprocess.run(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                info["cpu_model"] = result.stdout.strip()
            result = subprocess.run(
                ["sysctl", "-n", "hw.memsize"], capture_output=True, text=True
            )
            if result.returncode == 0:
                info["memory_gb"] = int(result.stdout.strip()) // (1024**3)
        except Exception:
            pass

    return info


# ---------------------------------------------------------------------------
# MD tool runners (return per-frame total SASA in nm²)
# ---------------------------------------------------------------------------


def run_mdtraj_native(xtc: Path, pdb: Path, n_points: int, stride: int) -> np.ndarray:
    """Run mdtraj.shrake_rupley. Returns per-frame total SASA in nm²."""
    import mdtraj as md

    console.print("  Loading trajectory with mdtraj...")
    traj = md.load(str(xtc), top=str(pdb), stride=stride)
    console.print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms")

    console.print("  Computing SASA with md.shrake_rupley...")
    sasa = md.shrake_rupley(traj, n_sphere_points=n_points)
    # sasa shape: (n_frames, n_atoms) in nm²
    return sasa.sum(axis=1)


def run_zsasa_mdtraj(
    xtc: Path, pdb: Path, n_points: int, stride: int, threads: int
) -> np.ndarray:
    """Run zsasa via mdtraj wrapper. Returns per-frame total SASA in nm²."""
    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    console.print("  Loading trajectory with mdtraj...")
    traj = md.load(str(xtc), top=str(pdb), stride=stride)
    console.print(f"  {traj.n_frames} frames, {traj.n_atoms} atoms")

    console.print(f"  Computing SASA with zsasa.mdtraj (threads={threads})...")
    return compute_sasa(traj, n_points=n_points, n_threads=threads, mode="total")


def run_zsasa_cli(
    xtc: Path, pdb: Path, n_points: int, stride: int, threads: int
) -> np.ndarray:
    """Run zsasa CLI traj mode. Returns per-frame total SASA in nm²."""
    root = get_root_dir()
    zsasa = root.joinpath("zig-out", "bin", "zsasa")
    if not zsasa.exists():
        console.print(f"[red]zsasa not found: {zsasa}[/]")
        return np.array([])

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
            "-o",
            out_path,
            "-q",
        ]
        console.print(f"  [dim]$ {' '.join(cmd)}[/]")
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if proc.returncode != 0:
            console.print(f"[red]zsasa traj failed: {proc.stderr[:500]}[/]")
            return np.array([])

        # Parse CSV: frame,time,total_sasa (Å²)
        totals = []
        with open(out_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                totals.append(float(row["total_sasa"]))

    # Convert Å² to nm² (1 nm² = 100 Å²)
    return np.array(totals) / 100.0


# ---------------------------------------------------------------------------
# Statistics and plotting
# ---------------------------------------------------------------------------


def compute_stats(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Compute R², mean/max relative error between two arrays."""
    mask = y > 0
    x = x[mask]
    y = y[mask]

    if len(x) < 2:
        return {"r_squared": 0.0, "mean_rel_error": 0.0, "max_rel_error": 0.0}

    correlation = np.corrcoef(x, y)[0, 1]
    r_squared = float(correlation**2)

    rel_errors = np.abs(x - y) / y * 100
    mean_rel_error = float(np.mean(rel_errors))
    max_rel_error = float(np.max(rel_errors))

    return {
        "r_squared": r_squared,
        "mean_rel_error": mean_rel_error,
        "max_rel_error": max_rel_error,
    }


def _plot_scatter_panel(
    ax: object,
    ref_arr: np.ndarray,
    comp_arr: np.ndarray,
    ref_label: str,
    comp_label: str,
    n_frames: int,
) -> None:
    """Draw a single scatter panel with y=x line and stats box."""
    ax.scatter(comp_arr, ref_arr, alpha=0.3, s=10, color="#3498db")

    max_val = max(float(np.max(ref_arr)), float(np.max(comp_arr)))
    min_val = min(float(np.min(ref_arr)), float(np.min(comp_arr)))
    ax.plot([min_val, max_val], [min_val, max_val], "r--", linewidth=1.5, label="y = x")

    stats = compute_stats(comp_arr, ref_arr)

    ax.set_xlabel(f"{comp_label} SASA (nm²)")
    ax.set_ylabel(f"{ref_label} SASA (nm²)")
    ax.set_title(f"{comp_label} vs {ref_label} (n={n_frames:,})")
    ax.legend(fontsize=8)

    stats_text = (
        f"R² = {stats['r_squared']:.6f}\n"
        f"Mean err = {stats['mean_rel_error']:.4f}%\n"
        f"Max err = {stats['max_rel_error']:.4f}%"
    )
    ax.text(
        0.05,
        0.95,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax.set_aspect("equal")


def _extract_ref_base(reference: str) -> str:
    """Extract base tool name from reference column.

    Examples: 'mdtraj_960' -> 'mdtraj', 'mdtraj' -> 'mdtraj'
    """
    import re

    m = re.match(r"^(.+?)_(\d+)$", reference)
    return m.group(1) if m else reference


def _group_columns_by_npoints(columns: list[str]) -> dict[str, list[str]]:
    """Group data columns by n_points suffix.

    Returns dict mapping n_points label to list of column names.
    For single n_points (no suffix), key is "".
    """
    import re

    groups: dict[str, list[str]] = {}
    for col in columns:
        if col == "frame":
            continue
        m = re.search(r"_(\d+)$", col)
        key = m.group(1) if m else ""
        groups.setdefault(key, []).append(col)
    return groups


def _build_comparison_pairs_by_npoints(
    columns: list[str],
    reference: str,
) -> dict[str, list[tuple[str, str]]]:
    """Build comparison pairs grouped by n_points.

    Within each group:
    - reference tool vs each other tool at same n_points
    - zsasa_mdtraj vs zsasa_cli (XTC reader comparison)
    """
    ref_base = _extract_ref_base(reference)
    groups = _group_columns_by_npoints(columns)
    result: dict[str, list[tuple[str, str]]] = {}

    for npoints_key, cols in sorted(
        groups.items(), key=lambda x: int(x[0]) if x[0] else 0
    ):
        suffix = f"_{npoints_key}" if npoints_key else ""
        pairs: list[tuple[str, str]] = []

        ref_col = f"{ref_base}{suffix}" if npoints_key else ref_base

        # Reference vs each other tool (same n_points)
        if ref_col in cols:
            for other in sorted(cols):
                if other != ref_col:
                    pairs.append((ref_col, other))

        # XTC comparison: zsasa_mdtraj vs zsasa_cli
        zsasa_mdtraj_col = f"zsasa_mdtraj{suffix}"
        zsasa_cli_col = f"zsasa_cli{suffix}"
        if zsasa_mdtraj_col in cols and zsasa_cli_col in cols:
            pair = (zsasa_mdtraj_col, zsasa_cli_col)
            if pair not in pairs:
                pairs.append(pair)

        if pairs:
            result[npoints_key] = pairs

    return result


def generate_md_scatter_plot(
    results_dir: Path,
    csv_path: Path,
    reference: str = "mdtraj",
) -> None:
    """Generate scatter plots comparing MD tools, one image per n_points."""
    import math

    import matplotlib.pyplot as plt
    import polars as pl

    df = pl.read_csv(csv_path)

    groups = _build_comparison_pairs_by_npoints(df.columns, reference)
    if not groups:
        console.print("[yellow]No comparison pairs found, skipping plot[/]")
        return

    for npoints_key, pairs in groups.items():
        n = len(pairs)
        ncols = min(n, 3)
        nrows = math.ceil(n / ncols)

        fig, axes = plt.subplots(
            nrows, ncols, figsize=(6 * ncols, 6 * nrows), squeeze=False
        )

        for idx, (ref_col, comp_col) in enumerate(pairs):
            row, col_idx = divmod(idx, ncols)
            ax = axes[row][col_idx]

            pair_df = df.select([ref_col, comp_col]).drop_nulls()
            if pair_df.height < 2:
                ax.set_title(f"{comp_col}: insufficient data")
                continue

            ref_arr = pair_df[ref_col].to_numpy()
            comp_arr = pair_df[comp_col].to_numpy()

            _plot_scatter_panel(
                ax, ref_arr, comp_arr, ref_col, comp_col, pair_df.height
            )

        # Hide unused axes
        for idx in range(n, nrows * ncols):
            row, col_idx = divmod(idx, ncols)
            axes[row][col_idx].set_visible(False)

        fig.tight_layout()

        if npoints_key:
            out_path = results_dir.joinpath(f"validation_md_{npoints_key}.png")
        else:
            out_path = results_dir.joinpath("validation_md.png")

        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        console.print(f"[green]Saved:[/] {out_path}")


def print_md_stats_table(
    csv_path: Path,
    reference: str = "mdtraj",
) -> None:
    """Print statistics table comparing MD tools, grouped by n_points."""
    import polars as pl

    df = pl.read_csv(csv_path)

    groups = _build_comparison_pairs_by_npoints(df.columns, reference)

    table = Table(title="MD Validation Statistics")
    table.add_column("n_points", style="bold")
    table.add_column("Reference", style="dim")
    table.add_column("Tool", style="cyan")
    table.add_column("Frames", justify="right")
    table.add_column("R²", justify="right")
    table.add_column("Mean Error %", justify="right")
    table.add_column("Max Error %", justify="right")

    for npoints_key, pairs in groups.items():
        for ref_col, comp_col in pairs:
            pair = df.select([ref_col, comp_col]).drop_nulls()
            if pair.height < 2:
                continue

            ref_arr = pair[ref_col].to_numpy()
            comp_arr = pair[comp_col].to_numpy()
            stats = compute_stats(comp_arr, ref_arr)

            table.add_row(
                npoints_key if npoints_key else "-",
                ref_col,
                comp_col,
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
    tools: Annotated[
        list[MdTool] | None,
        typer.Option(
            "--tool",
            "-t",
            help="Tools to compare. Default: all",
        ),
    ] = None,
    n_points: Annotated[
        str,
        typer.Option(
            "--n-points",
            help="Test points per atom (comma-separated for convergence: 100,500,960)",
        ),
    ] = "100",
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
            help="Number of threads (zsasa tools)",
        ),
    ] = 1,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output directory (default: benchmarks/results/validation_md/<name>)",
        ),
    ] = None,
    reference: Annotated[
        str,
        typer.Option(
            "--reference",
            "-r",
            help="Reference tool for comparison",
        ),
    ] = "mdtraj",
) -> None:
    """Compare per-frame SASA values across MD trajectory tools."""
    selected_tools = tools if tools else ALL_MD_TOOLS
    n_points_list = parse_n_points(n_points)
    multi_npoints = len(n_points_list) > 1

    # Set up output
    root = get_root_dir()
    if output_dir is None:
        results_dir = root.joinpath("benchmarks", "results", "validation_md", name)
    else:
        results_dir = output_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    # Save config
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "xtc": str(xtc),
            "pdb": str(pdb),
            "tools": [t.value for t in selected_tools],
            "n_points": n_points_list,
            "stride": stride,
            "threads": threads,
            "reference": reference,
        },
    }
    config_path = results_dir.joinpath("config.json")
    config_path.write_text(json.dumps(config, indent=2))

    # Print header
    console.print(f"[bold]=== MD SASA Validation: {name} ===[/]")
    console.print(f"XTC: {xtc}")
    console.print(f"PDB: {pdb}")
    console.print(f"Tools: {', '.join(t.value for t in selected_tools)}")
    console.print(f"N-points: {n_points_list}, Stride: {stride}, Threads: {threads}")
    console.print()

    # Run each tool x n_points combination
    tool_results: dict[str, np.ndarray] = {}

    for np_ in n_points_list:
        suffix = f"_{np_}" if multi_npoints else ""
        console.print(f"[bold]=== n_points={np_} ===[/]")

        if MdTool.mdtraj in selected_tools:
            col = f"mdtraj{suffix}"
            console.print(f"[bold cyan]Running mdtraj native (n_points={np_})...[/]")
            tool_results[col] = run_mdtraj_native(xtc, pdb, np_, stride)
            console.print(f"  Got {len(tool_results[col])} frames")

        if MdTool.zsasa_mdtraj in selected_tools:
            col = f"zsasa_mdtraj{suffix}"
            console.print(f"[bold cyan]Running zsasa_mdtraj (n_points={np_})...[/]")
            tool_results[col] = run_zsasa_mdtraj(xtc, pdb, np_, stride, threads)
            console.print(f"  Got {len(tool_results[col])} frames")

        if MdTool.zsasa_cli in selected_tools:
            col = f"zsasa_cli{suffix}"
            console.print(f"[bold cyan]Running zsasa CLI traj (n_points={np_})...[/]")
            tool_results[col] = run_zsasa_cli(xtc, pdb, np_, stride, threads)
            console.print(f"  Got {len(tool_results[col])} frames")

    if not tool_results:
        console.print("[red]No results collected.[/]")
        raise typer.Exit(1)

    # Filter out tools with empty results
    tool_results = {k: v for k, v in tool_results.items() if len(v) > 0}
    if not tool_results:
        console.print("[red]All tools returned 0 frames.[/]")
        raise typer.Exit(1)

    # Find the minimum frame count (tools may differ slightly)
    n_frames = min(len(v) for v in tool_results.values())

    # Build CSV
    import polars as pl

    data: dict[str, list] = {"frame": list(range(n_frames))}
    for tool_name, values in tool_results.items():
        data[tool_name] = [round(float(v), 6) for v in values[:n_frames]]

    df = pl.DataFrame(data)
    csv_path = results_dir.joinpath("results.csv")
    df.write_csv(csv_path)
    console.print(f"\n[green]Saved:[/] {csv_path} ({n_frames} frames)")

    # Statistics and plot (reference is base tool name, per-group handled internally)
    print_md_stats_table(csv_path, reference=reference)
    generate_md_scatter_plot(results_dir, csv_path, reference=reference)

    console.print(f"\n[bold green]=== Done! Results: {results_dir} ===[/]")


@app.command()
def compare(
    results_dir: Annotated[
        Path,
        typer.Option(
            "--dir",
            "-d",
            help="Directory containing results.csv",
            exists=True,
            file_okay=False,
            dir_okay=True,
        ),
    ],
    reference: Annotated[
        str | None,
        typer.Option(
            "--reference",
            "-r",
            help="Reference tool (auto-detected from config.json if omitted)",
        ),
    ] = None,
) -> None:
    """Re-analyze existing MD validation results."""
    csv_path = results_dir.joinpath("results.csv")
    if not csv_path.exists():
        console.print(f"[red]results.csv not found in {results_dir}[/]")
        raise typer.Exit(1)

    # Auto-detect reference from config.json
    if reference is None:
        reference = _detect_reference(results_dir)

    console.print(f"[bold]=== Re-analyzing: {results_dir} ===[/]")
    console.print(f"Reference: {reference}")

    print_md_stats_table(csv_path, reference=reference)
    generate_md_scatter_plot(results_dir, csv_path, reference=reference)

    console.print(f"\n[bold green]=== Done! ===[/]")


if __name__ == "__main__":
    app()

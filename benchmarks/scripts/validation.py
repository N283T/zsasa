#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["polars>=1.0", "matplotlib>=3.8", "numpy>=1.26", "rich>=13.0", "typer>=0.9.0"]
# ///
"""SASA validation: compare accuracy across tools.

Runs tools on a PDB directory, collects per-file SASA values,
and compares across tools -- completely independent of timing benchmarks.

Usage:
    # Run tools and compare (both f32 and f64)
    ./benchmarks/scripts/validation.py run \
        -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
        -n ecoli \
        --algorithm sr --threads 1

    # Run with specific precision
    ./benchmarks/scripts/validation.py run \
        -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
        -n ecoli \
        --algorithm sr --precision f64 --threads 1

    # Re-analyze existing CSV
    ./benchmarks/scripts/validation.py compare \
        -d benchmarks/results/validation/ecoli/sr

Output:
    benchmarks/results/validation/<name>/<algorithm>/
    ├── config.json          # System info, parameters
    ├── results.csv          # Per-file SASA values
    └── validation_sr.png    # Scatter plot
"""

from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import tempfile
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="SASA validation: compare accuracy across tools")
console = Console()


# ---------------------------------------------------------------------------
# Shared utilities (adapted from bench_batch.py)
# ---------------------------------------------------------------------------


class Tool(str, Enum):
    zig = "zig"
    freesasa = "freesasa"


ALL_TOOLS = [Tool.zig, Tool.freesasa]


def get_root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def get_binary_paths() -> dict[str, Path]:
    """Get paths to tool binaries."""
    root = get_root_dir()
    return {
        "zsasa": root.joinpath("zig-out", "bin", "zsasa"),
        "freesasa": root.joinpath(
            "benchmarks", "external", "freesasa-bench", "src", "freesasa"
        ),
    }


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
# Tool runners
# ---------------------------------------------------------------------------


def run_zsasa(
    input_dir: Path,
    algorithm: str,
    precision: str,
    threads: int,
    binaries: dict[str, Path],
    n_points: int = 100,
    use_bitmask: bool = False,
) -> dict[str, tuple[float, int]]:
    """Run zsasa in batch mode. Returns {stem: (total_sasa, n_atoms)}."""
    zsasa = binaries["zsasa"]
    if not zsasa.exists():
        console.print(f"[red]zsasa not found: {zsasa}[/]")
        return {}

    results: dict[str, tuple[float, int]] = {}

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
        console.print(f"  [dim]$ {' '.join(cmd)}[/]")
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
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
    algorithm: str,
    binaries: dict[str, Path],
) -> dict[str, float]:
    """Run FreeSASA per-file via CLI. Returns {stem: total_sasa}."""
    binary = binaries["freesasa"]
    if not binary.exists():
        console.print(f"[red]freesasa not found: {binary}[/]")
        return {}

    pdb_files = sorted(input_dir.glob("*.pdb"))
    if not pdb_files:
        pdb_files = sorted(input_dir.glob("*.ent"))
    if not pdb_files:
        console.print("[red]No PDB files found in input directory[/]")
        return {}

    results: dict[str, float] = {}

    from rich.progress import Progress

    with Progress(console=console) as progress:
        task = progress.add_task("  FreeSASA", total=len(pdb_files))
        for pdb_file in pdb_files:
            cmd = [str(binary), str(pdb_file)]
            if algorithm == "sr":
                cmd.extend(["--shrake-rupley", "--resolution=100"])
            else:
                cmd.extend(["--lee-richards", "--resolution=20"])

            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                if proc.returncode != 0:
                    progress.advance(task)
                    continue

                for line in proc.stdout.split("\n"):
                    if line.startswith("Total"):
                        match = re.search(r"[\d.]+", line)
                        if match:
                            results[pdb_file.stem] = float(match.group())
                            break
            except (subprocess.TimeoutExpired, ValueError):
                pass

            progress.advance(task)

    return results


# ---------------------------------------------------------------------------
# Statistics and plotting
# ---------------------------------------------------------------------------


def compute_stats(x: list[float], y: list[float]) -> dict[str, float]:
    """Compute R², mean/max relative error between two lists."""
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
    r_squared = float(correlation**2)

    rel_errors = np.abs(xa - ya) / ya * 100
    mean_rel_error = float(np.mean(rel_errors))
    max_rel_error = float(np.max(rel_errors))

    return {
        "r_squared": r_squared,
        "mean_rel_error": mean_rel_error,
        "max_rel_error": max_rel_error,
    }


def generate_scatter_plot(
    results_dir: Path,
    csv_path: Path,
    algorithm: str,
    reference: str = "freesasa",
) -> None:
    """Generate scatter plot comparing tools against reference."""
    import matplotlib.pyplot as plt
    import numpy as np
    import polars as pl

    df = pl.read_csv(csv_path)

    if reference not in df.columns:
        console.print(
            f"[yellow]Reference column '{reference}' not in CSV, skipping plot[/]"
        )
        return

    ref_values = df[reference]

    # Find comparison columns (everything except structure, n_atoms, and reference)
    skip_cols = {"structure", "n_atoms", reference}
    compare_cols = [c for c in df.columns if c not in skip_cols]

    if not compare_cols:
        console.print("[yellow]No comparison columns found, skipping plot[/]")
        return

    fig, axes = plt.subplots(
        1, len(compare_cols), figsize=(8 * len(compare_cols), 8), squeeze=False
    )

    for idx, col in enumerate(compare_cols):
        ax = axes[0][idx]

        # Drop nulls for this pair
        pair = df.select([reference, col]).drop_nulls()
        if pair.height < 2:
            ax.set_title(f"{col}: insufficient data")
            continue

        ref_arr = pair[reference].to_numpy()
        comp_arr = pair[col].to_numpy()

        ax.scatter(comp_arr, ref_arr, alpha=0.3, s=10, color="#3498db")

        max_val = max(float(np.max(ref_arr)), float(np.max(comp_arr)))
        ax.plot([0, max_val], [0, max_val], "r--", linewidth=1.5, label="y = x")

        stats = compute_stats(comp_arr.tolist(), ref_arr.tolist())

        ax.set_xlabel(f"{col} SASA")
        ax.set_ylabel(f"{reference} SASA")
        ax.set_title(f"{algorithm.upper()}: {col} vs {reference} (n={pair.height:,})")
        ax.legend()

        stats_text = (
            f"R² = {stats['r_squared']:.6f}\n"
            f"Mean error = {stats['mean_rel_error']:.4f}%\n"
            f"Max error = {stats['max_rel_error']:.4f}%"
        )
        ax.text(
            0.05,
            0.95,
            stats_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )

        ax.set_aspect("equal")

    fig.tight_layout()
    out_path = results_dir.joinpath(f"validation_{algorithm}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    console.print(f"[green]Saved:[/] {out_path}")


def print_stats_table(
    csv_path: Path,
    reference: str = "freesasa",
) -> None:
    """Print statistics table comparing tools against reference."""
    import polars as pl

    df = pl.read_csv(csv_path)

    if reference not in df.columns:
        console.print(f"[yellow]Reference column '{reference}' not in CSV[/]")
        return

    skip_cols = {"structure", "n_atoms", reference}
    compare_cols = [c for c in df.columns if c not in skip_cols]

    table = Table(title=f"Validation Statistics (reference: {reference})")
    table.add_column("Tool", style="cyan")
    table.add_column("N", justify="right")
    table.add_column("R²", justify="right")
    table.add_column("Mean Error %", justify="right")
    table.add_column("Max Error %", justify="right")

    for col in compare_cols:
        pair = df.select([reference, col]).drop_nulls()
        if pair.height < 2:
            continue

        ref_arr = pair[reference].to_list()
        comp_arr = pair[col].to_list()
        stats = compute_stats(comp_arr, ref_arr)

        table.add_row(
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
    tools: Annotated[
        list[Tool] | None,
        typer.Option(
            "--tool",
            "-t",
            help="Tools to compare (can specify multiple). Default: all",
        ),
    ] = None,
    algorithm: Annotated[
        str,
        typer.Option(
            "--algorithm",
            "-a",
            help="Algorithm: sr (shrake-rupley) or lr (lee-richards)",
        ),
    ] = "sr",
    precision: Annotated[
        str | None,
        typer.Option(
            "--precision",
            "-p",
            help="zsasa precision: f32, f64, or omit for both",
        ),
    ] = None,
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            "-T",
            help="Number of threads (zsasa)",
        ),
    ] = 1,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output directory (default: benchmarks/results/validation/<name>)",
        ),
    ] = None,
    reference: Annotated[
        str,
        typer.Option(
            "--reference",
            "-r",
            help="Reference tool for comparison",
        ),
    ] = "freesasa",
    n_points: Annotated[
        int,
        typer.Option(
            "--n-points",
            "-N",
            help="Number of sphere test points per atom (default: 100)",
        ),
    ] = 100,
    use_bitmask: Annotated[
        bool,
        typer.Option(
            "--use-bitmask",
            help="Use bitmask neighbor list for zsasa",
        ),
    ] = False,
) -> None:
    """Run tools on PDB directory and compare SASA values."""
    selected_tools = tools if tools else ALL_TOOLS
    precisions = [precision] if precision else ["f32", "f64"]

    # Set up output
    root = get_root_dir()
    if output_dir is None:
        results_dir = root.joinpath(
            "benchmarks", "results", "validation", name, algorithm
        )
    else:
        results_dir = output_dir
    results_dir.mkdir(parents=True, exist_ok=True)

    binaries = get_binary_paths()

    # Count input files
    n_files = sum(
        1 for f in input_dir.iterdir() if f.is_file() and f.suffix in {".pdb", ".ent"}
    )

    # Save config
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "input_dir": str(input_dir),
            "n_files": n_files,
            "tools": [t.value for t in selected_tools],
            "algorithm": algorithm,
            "precision": precisions,
            "threads": threads,
            "n_points": n_points,
            "use_bitmask": use_bitmask,
            "reference": reference,
        },
    }
    config_path = results_dir.joinpath("config.json")
    config_path.write_text(json.dumps(config, indent=2))

    # Print header
    console.print(f"[bold]=== SASA Validation: {name} ===[/]")
    console.print(f"Input: {input_dir} ({n_files} PDB files)")
    console.print(f"Output: {results_dir}")
    console.print(f"Tools: {', '.join(t.value for t in selected_tools)}")
    console.print(
        f"Algorithm: {algorithm}, Precision: {', '.join(precisions)}, Threads: {threads}"
    )
    console.print()

    # Collect results from each tool x precision
    zsasa_runs: dict[str, dict[str, tuple[float, int]]] = {}

    if Tool.zig in selected_tools:
        for prec in precisions:
            col = f"zsasa_{prec}"
            console.print(f"[bold cyan]Running zsasa ({prec})...[/]")
            zsasa_runs[col] = run_zsasa(
                input_dir, algorithm, prec, threads, binaries, n_points, use_bitmask
            )
            console.print(f"  Got {len(zsasa_runs[col])} results")

    freesasa_results: dict[str, float] = {}

    if Tool.freesasa in selected_tools:
        console.print("[bold cyan]Running FreeSASA...[/]")
        freesasa_results = run_freesasa(input_dir, algorithm, binaries)
        console.print(f"  Got {len(freesasa_results)} results")

    # Merge by stem
    import polars as pl

    all_stems: set[str] = set()
    for results in zsasa_runs.values():
        all_stems.update(results.keys())
    if freesasa_results:
        all_stems.update(freesasa_results.keys())

    if not all_stems:
        console.print("[red]No results collected. Check tool output above.[/]")
        raise typer.Exit(1)

    rows: list[dict] = []
    for stem in sorted(all_stems):
        row: dict = {"structure": stem}

        # n_atoms from any zsasa run
        for results in zsasa_runs.values():
            if stem in results:
                _, n_atoms = results[stem]
                row["n_atoms"] = n_atoms
                break

        for col, results in zsasa_runs.items():
            if stem in results:
                sasa, _ = results[stem]
                row[col] = round(sasa, 2)

        if freesasa_results and stem in freesasa_results:
            row["freesasa"] = round(freesasa_results[stem], 2)

        rows.append(row)

    # Build DataFrame with consistent column order
    columns = ["structure", "n_atoms"]
    columns.extend(zsasa_runs.keys())
    if freesasa_results:
        columns.append("freesasa")

    df = pl.DataFrame(rows).select(
        [c for c in columns if c in pl.DataFrame(rows).columns]
    )

    csv_path = results_dir.joinpath("results.csv")
    df.write_csv(csv_path)
    console.print(f"\n[green]Saved:[/] {csv_path} ({df.height} structures)")

    # Statistics and plot
    print_stats_table(csv_path, reference=reference)
    generate_scatter_plot(results_dir, csv_path, algorithm, reference=reference)

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
        str,
        typer.Option(
            "--reference",
            "-r",
            help="Reference tool for comparison",
        ),
    ] = "freesasa",
) -> None:
    """Re-analyze existing validation results."""
    csv_path = results_dir.joinpath("results.csv")
    if not csv_path.exists():
        console.print(f"[red]results.csv not found in {results_dir}[/]")
        raise typer.Exit(1)

    # Detect algorithm from config.json if available
    algorithm = "sr"
    config_path = results_dir.joinpath("config.json")
    if config_path.exists():
        try:
            config = json.loads(config_path.read_text())
            algorithm = config.get("parameters", {}).get("algorithm", "sr")
        except (json.JSONDecodeError, KeyError):
            pass

    console.print(f"[bold]=== Re-analyzing: {results_dir} ===[/]")
    console.print(f"Algorithm: {algorithm}, Reference: {reference}")

    print_stats_table(csv_path, reference=reference)
    generate_scatter_plot(results_dir, csv_path, algorithm, reference=reference)

    console.print(f"\n[bold green]=== Done! ===[/]")


if __name__ == "__main__":
    app()

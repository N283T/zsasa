#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "cffi",
#     "MDAnalysis",
#     "mdsasa-bolt",
#     "mdtraj",
#     "numpy",
#     "rich>=13.0",
#     "typer>=0.9.0",
# ]
# ///
"""Benchmark SASA calculation: zsasa vs MDTraj vs mdsasa-bolt.

Compares SASA calculation performance across implementations:
- zsasa (Zig, configurable threads)
- MDTraj shrake_rupley (single-threaded)
- mdsasa-bolt (RustSASA, all cores)

Usage:
    ./benchmarks/scripts/benchmark_md.py trajectory.xtc topology.pdb
    ./benchmarks/scripts/benchmark_md.py trajectory.xtc topology.pdb -t "1,4,8"
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Annotated

import mdtraj as md
import numpy as np
import typer
from rich.console import Console
from rich.table import Table

# Add zsasa to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))

from zsasa.xtc import compute_sasa_trajectory

console = Console()
app = typer.Typer()


def benchmark_zsasa(
    xtc_path: Path, top_path: Path, n_threads: int = 0, n_runs: int = 3
) -> dict:
    """Benchmark zsasa SASA calculation."""
    traj = md.load_frame(str(xtc_path), 0, top=str(top_path))
    radii = np.array([1.7] * traj.n_atoms, dtype=np.float32)

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = compute_sasa_trajectory(xtc_path, radii, n_threads=n_threads)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = result.n_frames
        n_atoms = result.n_atoms

    return {
        "name": "zsasa",
        "threads": str(n_threads) if n_threads > 0 else "auto",
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


def benchmark_mdtraj(xtc_path: Path, top_path: Path, n_runs: int = 3) -> dict:
    """Benchmark MDTraj shrake_rupley (single-threaded)."""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        traj = md.load(str(xtc_path), top=str(top_path))
        _ = md.shrake_rupley(traj, n_sphere_points=100)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = traj.n_frames
        n_atoms = traj.n_atoms

    return {
        "name": "MDTraj",
        "threads": "1",
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


def benchmark_mdsasa_bolt(
    xtc_path: Path, top_path: Path, n_runs: int = 3
) -> dict | None:
    """Benchmark mdsasa-bolt (RustSASA, all cores)."""
    try:
        import MDAnalysis as mda
        from mdsasa_bolt import SASAAnalysis
    except ImportError:
        return None

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        u = mda.Universe(str(top_path), str(xtc_path))
        sasa = SASAAnalysis(u.atoms, n_points=100)
        sasa.run()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = len(u.trajectory)
        n_atoms = len(u.atoms)

    return {
        "name": "mdsasa-bolt",
        "threads": "all",
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


@app.command()
def benchmark(
    xtc_path: Annotated[Path, typer.Argument(help="Path to XTC trajectory file")],
    top_path: Annotated[Path, typer.Argument(help="Path to topology file (PDB/GRO)")],
    n_runs: Annotated[int, typer.Option("--runs", "-n", help="Number of runs")] = 3,
    threads: Annotated[
        str,
        typer.Option("--threads", "-t", help="Thread counts for zsasa (e.g., '1,4,8')"),
    ] = "1,4,8",
) -> None:
    """Benchmark SASA calculation across implementations."""
    if not xtc_path.exists():
        console.print(f"[red]Error: XTC file not found: {xtc_path}[/red]")
        raise typer.Exit(1)
    if not top_path.exists():
        console.print(f"[red]Error: Topology file not found: {top_path}[/red]")
        raise typer.Exit(1)

    # Get trajectory info
    traj = md.load_frame(str(xtc_path), 0, top=str(top_path))
    with md.open(str(xtc_path)) as f:
        n_frames = len(f)

    console.print(f"\n[bold]SASA Benchmark[/bold]")
    console.print(f"Trajectory: {xtc_path}")
    console.print(f"Frames: {n_frames}, Atoms: {traj.n_atoms}")
    console.print(f"Runs: {n_runs}\n")

    results = []

    # zsasa with different thread counts
    thread_counts = [int(t) for t in threads.split(",")]
    for n_threads in thread_counts:
        console.print(f"[cyan]zsasa (threads={n_threads})...[/cyan]")
        results.append(benchmark_zsasa(xtc_path, top_path, n_threads, n_runs))

    # MDTraj (single-threaded)
    console.print("[cyan]MDTraj (single-threaded)...[/cyan]")
    mdtraj_result = benchmark_mdtraj(xtc_path, top_path, n_runs)
    results.append(mdtraj_result)

    # mdsasa-bolt (all cores)
    console.print("[cyan]mdsasa-bolt (all cores)...[/cyan]")
    bolt_result = benchmark_mdsasa_bolt(xtc_path, top_path, n_runs)
    if bolt_result:
        results.append(bolt_result)
    else:
        console.print("[yellow]  mdsasa-bolt not installed, skipping[/yellow]")

    # Results table
    console.print()
    table = Table(title="SASA Benchmark Results")
    table.add_column("Method", style="cyan")
    table.add_column("Threads", justify="right")
    table.add_column("Time (s)", justify="right")
    table.add_column("fps", justify="right")
    table.add_column("vs MDTraj", justify="right")
    if bolt_result:
        table.add_column("vs mdsasa-bolt", justify="right")

    mdtraj_time = mdtraj_result["mean"]
    bolt_time = bolt_result["mean"] if bolt_result else None

    for r in results:
        fps = r["n_frames"] / r["mean"]
        vs_mdtraj = mdtraj_time / r["mean"]
        vs_mdtraj_str = f"{vs_mdtraj:.1f}x" if r["name"] != "MDTraj" else "-"

        row = [r["name"], r["threads"], f"{r['mean']:.2f}", f"{fps:.0f}", vs_mdtraj_str]

        if bolt_result:
            vs_bolt = bolt_time / r["mean"]
            vs_bolt_str = f"{vs_bolt:.1f}x" if r["name"] != "mdsasa-bolt" else "-"
            row.append(vs_bolt_str)

        table.add_row(*row)

    console.print(table)


def main(
    xtc_path: Annotated[Path, typer.Argument(help="Path to XTC trajectory file")],
    top_path: Annotated[Path, typer.Argument(help="Path to topology file (PDB/GRO)")],
    n_runs: Annotated[int, typer.Option("--runs", "-n", help="Number of runs")] = 3,
    threads: Annotated[
        str,
        typer.Option("--threads", "-t", help="Thread counts for zsasa (e.g., '1,4,8')"),
    ] = "1,4,8",
) -> None:
    """Benchmark SASA calculation across implementations."""
    benchmark(xtc_path, top_path, n_runs, threads)


if __name__ == "__main__":
    typer.run(main)

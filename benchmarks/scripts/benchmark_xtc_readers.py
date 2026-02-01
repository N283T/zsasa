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
"""Benchmark XTC readers and SASA calculation.

Compares:
1. XTC Reading: zsasa (Zig) vs MDTraj (C xdrfile)
2. SASA Calculation: zsasa vs MDTraj shrake_rupley vs mdsasa-bolt (RustSASA)

Usage:
    # XTC reading only
    uv run benchmarks/scripts/benchmark_xtc_readers.py trajectory.xtc topology.pdb

    # Include SASA calculation
    uv run benchmarks/scripts/benchmark_xtc_readers.py trajectory.xtc topology.pdb --sasa

    # Specify thread counts
    uv run benchmarks/scripts/benchmark_xtc_readers.py trajectory.xtc topology.pdb --sasa -t "1,4,8"
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Annotated

import mdtraj as md
import numpy as np
import typer
from rich.console import Console
from rich.table import Table

# Add zsasa to path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))

from zsasa.xtc import XtcReader, compute_sasa_trajectory

console = Console()


def benchmark_zsasa_read(xtc_path: Path, n_runs: int = 3) -> dict:
    """Benchmark zsasa native XTC reader."""
    times = []

    for _ in range(n_runs):
        start = time.perf_counter()
        frame_count = 0
        with XtcReader(xtc_path) as reader:
            natoms = reader.natoms
            for frame in reader:
                frame_count += 1
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "name": "zsasa (Zig)",
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": frame_count,
        "n_atoms": natoms,
    }


def benchmark_mdtraj_read(xtc_path: Path, top_path: Path, n_runs: int = 3) -> dict:
    """Benchmark MDTraj XTC reader."""
    times = []

    for _ in range(n_runs):
        start = time.perf_counter()
        traj = md.load(str(xtc_path), top=str(top_path))
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = traj.n_frames
        n_atoms = traj.n_atoms

    return {
        "name": "MDTraj (C xdrfile)",
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


def benchmark_zsasa_sasa(
    xtc_path: Path, top_path: Path, n_threads: int = 0, n_runs: int = 3
) -> dict:
    """Benchmark zsasa native XTC + SASA calculation."""
    # Get radii from topology
    traj = md.load_frame(str(xtc_path), 0, top=str(top_path))
    radii = np.array([1.7] * traj.n_atoms, dtype=np.float32)  # Simplified

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = compute_sasa_trajectory(xtc_path, radii, n_threads=n_threads)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = result.n_frames
        n_atoms = result.n_atoms

    return {
        "name": f"zsasa SASA (threads={n_threads or 'auto'})",
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


def benchmark_mdtraj_sasa(xtc_path: Path, top_path: Path, n_runs: int = 3) -> dict:
    """Benchmark MDTraj shrake_rupley SASA calculation.

    Note: MDTraj shrake_rupley doesn't support n_jobs parameter.
    It's single-threaded.
    """
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        traj = md.load(str(xtc_path), top=str(top_path))
        sasa = md.shrake_rupley(traj, n_sphere_points=100)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = traj.n_frames
        n_atoms = traj.n_atoms

    return {
        "name": "MDTraj shrake_rupley (single)",
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


def benchmark_mdsasa_bolt(
    xtc_path: Path, top_path: Path, n_runs: int = 3
) -> dict | None:
    """Benchmark mdsasa-bolt (RustSASA) if available.

    Note: mdsasa-bolt uses rayon with no thread control - always uses all cores.
    """
    try:
        import MDAnalysis as mda
        from mdsasa_bolt import SASAAnalysis as MdsasaAnalysis
    except ImportError:
        return None

    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        u = mda.Universe(str(top_path), str(xtc_path))
        sasa = MdsasaAnalysis(u.atoms, n_points=100)
        sasa.run()
        elapsed = time.perf_counter() - start
        times.append(elapsed)
        n_frames = len(u.trajectory)
        n_atoms = len(u.atoms)

    return {
        "name": "mdsasa-bolt (all cores)",
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": n_frames,
        "n_atoms": n_atoms,
    }


def benchmark_mdtraj_iterload(
    xtc_path: Path, top_path: Path, chunk: int = 100, n_runs: int = 3
) -> dict:
    """Benchmark MDTraj iterload (streaming)."""
    times = []

    for _ in range(n_runs):
        start = time.perf_counter()
        frame_count = 0
        for chunk_traj in md.iterload(str(xtc_path), top=str(top_path), chunk=chunk):
            frame_count += chunk_traj.n_frames
            n_atoms = chunk_traj.n_atoms
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        "name": f"MDTraj iterload (chunk={chunk})",
        "times": times,
        "mean": np.mean(times),
        "std": np.std(times),
        "n_frames": frame_count,
        "n_atoms": n_atoms,
    }


def main(
    xtc_path: Annotated[Path, typer.Argument(help="Path to XTC trajectory file")],
    top_path: Annotated[
        Path, typer.Argument(help="Path to topology file (PDB/GRO)")
    ] = None,
    n_runs: Annotated[int, typer.Option("--runs", "-n", help="Number of runs")] = 3,
    chunk: Annotated[
        int, typer.Option("--chunk", "-c", help="Chunk size for iterload")
    ] = 100,
    sasa: Annotated[
        bool, typer.Option("--sasa", "-s", help="Include SASA calculation benchmark")
    ] = False,
    threads: Annotated[
        str,
        typer.Option("--threads", "-t", help="Thread counts for SASA (e.g., '1,4,8')"),
    ] = "1,4,8",
) -> None:
    """Benchmark XTC readers."""
    if not xtc_path.exists():
        console.print(f"[red]Error: XTC file not found: {xtc_path}[/red]")
        raise typer.Exit(1)

    console.print(f"\n[bold]Benchmarking XTC readers[/bold]")
    console.print(f"File: {xtc_path}")
    console.print(f"Size: {xtc_path.stat().st_size / 1024 / 1024:.1f} MB")
    console.print(f"Runs: {n_runs}\n")

    results = []

    # zsasa native reader (no topology needed)
    console.print("[cyan]Running zsasa native reader...[/cyan]")
    zsasa_result = benchmark_zsasa_read(xtc_path, n_runs)
    results.append(zsasa_result)
    console.print(
        f"  {zsasa_result['n_frames']} frames, {zsasa_result['n_atoms']} atoms"
    )

    # MDTraj requires topology
    if top_path and top_path.exists():
        console.print("[cyan]Running MDTraj load...[/cyan]")
        mdtraj_result = benchmark_mdtraj_read(xtc_path, top_path, n_runs)
        results.append(mdtraj_result)

        console.print("[cyan]Running MDTraj iterload...[/cyan]")
        iterload_result = benchmark_mdtraj_iterload(xtc_path, top_path, chunk, n_runs)
        results.append(iterload_result)
    else:
        console.print(
            "[yellow]Skipping MDTraj benchmarks (no topology provided)[/yellow]"
        )

    # SASA calculation benchmarks
    if sasa and top_path and top_path.exists():
        console.print("\n[bold]SASA Calculation Benchmarks[/bold]\n")
        sasa_results = []

        thread_counts = [int(t) for t in threads.split(",")]
        for n_threads in thread_counts:
            console.print(f"[cyan]Running zsasa SASA (threads={n_threads})...[/cyan]")
            sasa_result = benchmark_zsasa_sasa(xtc_path, top_path, n_threads, n_runs)
            sasa_results.append(sasa_result)

        # MDTraj (single-threaded only - doesn't support n_jobs)
        console.print("[cyan]Running MDTraj shrake_rupley (single-threaded)...[/cyan]")
        mdtraj_sasa_result = benchmark_mdtraj_sasa(xtc_path, top_path, n_runs)
        sasa_results.append(mdtraj_sasa_result)

        # mdsasa-bolt (RustSASA) if available
        console.print("[cyan]Running mdsasa-bolt (RustSASA)...[/cyan]")
        bolt_result = benchmark_mdsasa_bolt(xtc_path, top_path, n_runs)
        if bolt_result:
            sasa_results.append(bolt_result)
        else:
            console.print("[yellow]  mdsasa not installed, skipping[/yellow]")

        # SASA results table
        console.print()
        sasa_table = Table(title="SASA Calculation Benchmark Results")
        sasa_table.add_column("Method", style="cyan")
        sasa_table.add_column("Mean (s)", justify="right")
        sasa_table.add_column("Std (s)", justify="right")
        sasa_table.add_column("Throughput (frames/s)", justify="right")
        sasa_table.add_column("vs baseline", justify="right")

        # Use mdsasa-bolt as baseline if available, else first MDTraj result
        baseline_result = bolt_result if bolt_result else sasa_results[-1]
        baseline = baseline_result["mean"]
        baseline_name = baseline_result["name"]

        console.print(f"[dim]Baseline: {baseline_name}[/dim]\n")

        for r in sasa_results:
            speedup = baseline / r["mean"] if r["mean"] > 0 else 0
            speedup_str = f"{speedup:.2f}x" if r["name"] != baseline_name else "-"
            throughput = r["n_frames"] / r["mean"]
            sasa_table.add_row(
                r["name"],
                f"{r['mean']:.3f}",
                f"{r['std']:.3f}",
                f"{throughput:.0f}",
                speedup_str,
            )

        console.print(sasa_table)

    # Results table
    console.print()
    table = Table(title="XTC Reader Benchmark Results")
    table.add_column("Reader", style="cyan")
    table.add_column("Mean (s)", justify="right")
    table.add_column("Std (s)", justify="right")
    table.add_column("Throughput (frames/s)", justify="right")
    table.add_column("Speedup", justify="right")

    baseline = results[0]["mean"]  # zsasa as baseline
    for r in results:
        speedup = baseline / r["mean"] if r["mean"] > 0 else 0
        speedup_str = f"{speedup:.2f}x" if r["name"] != "zsasa (Zig)" else "-"
        throughput = r["n_frames"] / r["mean"]
        table.add_row(
            r["name"],
            f"{r['mean']:.3f}",
            f"{r['std']:.3f}",
            f"{throughput:.0f}",
            speedup_str,
        )

    console.print(table)

    # File read speed
    file_size_mb = xtc_path.stat().st_size / 1024 / 1024
    console.print(f"\n[bold]File read speed:[/bold]")
    for r in results:
        speed = file_size_mb / r["mean"]
        console.print(f"  {r['name']}: {speed:.1f} MB/s")


if __name__ == "__main__":
    typer.run(main)

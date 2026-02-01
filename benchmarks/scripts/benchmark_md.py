#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "numpy>=2.0",
#     "mdtraj>=1.10",
#     "typer>=0.15",
#     "rich>=13.0",
#     "cffi>=1.16",
# ]
# ///
"""Benchmark MD trajectory SASA calculation."""

import time
from pathlib import Path

import mdtraj as md
import numpy as np
import typer
from rich.console import Console
from rich.table import Table

# Add python directory to path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "python"))

from zsasa import calculate_sasa_batch
from zsasa.mdtraj import compute_sasa

console = Console()
app = typer.Typer()


def time_func(func, *args, n_runs: int = 3, **kwargs):
    """Time a function and return mean/std."""
    times = []
    for _ in range(n_runs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        times.append(elapsed)
    return result, np.mean(times), np.std(times)


@app.command()
def benchmark(
    trajectory: Path = typer.Argument(..., help="Path to trajectory file (.xtc)"),
    topology: Path = typer.Argument(..., help="Path to topology file (.pdb or .tpr)"),
    n_frames: int = typer.Option(100, "--frames", "-f", help="Number of frames to use"),
    n_threads: int = typer.Option(
        0, "--threads", "-t", help="Number of threads (0=auto)"
    ),
    algorithm: str = typer.Option(
        "sr", "--algorithm", "-a", help="Algorithm: sr or lr"
    ),
    n_runs: int = typer.Option(3, "--runs", "-r", help="Number of runs for timing"),
):
    """Benchmark MD trajectory SASA calculation."""
    console.print(f"\n[bold]Loading trajectory:[/bold] {trajectory}")
    console.print(f"[bold]Topology:[/bold] {topology}")

    # Load trajectory
    traj = md.load(str(trajectory), top=str(topology))
    console.print(f"[bold]Total frames:[/bold] {traj.n_frames}")
    console.print(f"[bold]Atoms:[/bold] {traj.n_atoms}")

    # Subsample if needed
    if n_frames < traj.n_frames:
        step = traj.n_frames // n_frames
        traj = traj[::step][:n_frames]
        console.print(f"[bold]Using frames:[/bold] {traj.n_frames}")

    # Warm up
    console.print("\n[yellow]Warming up...[/yellow]")
    _ = compute_sasa(traj[:10], n_threads=n_threads, algorithm=algorithm)

    # Benchmark zsasa
    console.print(f"\n[green]Benchmarking zsasa ({algorithm.upper()})...[/green]")
    _, zsasa_mean, zsasa_std = time_func(
        compute_sasa,
        traj,
        n_threads=n_threads,
        algorithm=algorithm,
        n_runs=n_runs,
    )

    # Benchmark mdtraj.shrake_rupley for comparison
    console.print("\n[green]Benchmarking MDTraj shrake_rupley...[/green]")
    _, mdtraj_mean, mdtraj_std = time_func(
        md.shrake_rupley,
        traj,
        n_runs=n_runs,
    )

    # Results table
    table = Table(
        title=f"Benchmark Results ({traj.n_frames} frames, {traj.n_atoms} atoms)"
    )
    table.add_column("Implementation", style="cyan")
    table.add_column("Time (s)", justify="right")
    table.add_column("Per-frame (ms)", justify="right")
    table.add_column("Speedup", justify="right", style="green")

    mdtraj_per_frame = (mdtraj_mean / traj.n_frames) * 1000
    zsasa_per_frame = (zsasa_mean / traj.n_frames) * 1000
    speedup = mdtraj_mean / zsasa_mean

    table.add_row(
        "MDTraj",
        f"{mdtraj_mean:.3f} ± {mdtraj_std:.3f}",
        f"{mdtraj_per_frame:.2f}",
        "1.00x",
    )
    table.add_row(
        f"zsasa ({algorithm.upper()})",
        f"{zsasa_mean:.3f} ± {zsasa_std:.3f}",
        f"{zsasa_per_frame:.2f}",
        f"{speedup:.2f}x",
    )

    console.print(table)

    # Throughput info
    console.print(f"\n[bold]Throughput:[/bold]")
    console.print(f"  MDTraj: {traj.n_frames / mdtraj_mean:.1f} frames/s")
    console.print(f"  zsasa:  {traj.n_frames / zsasa_mean:.1f} frames/s")

    # Atoms per second
    atoms_per_frame = traj.n_atoms
    console.print(f"\n[bold]Atoms per second:[/bold]")
    console.print(
        f"  MDTraj: {atoms_per_frame * traj.n_frames / mdtraj_mean / 1e6:.2f}M atoms/s"
    )
    console.print(
        f"  zsasa:  {atoms_per_frame * traj.n_frames / zsasa_mean / 1e6:.2f}M atoms/s"
    )


@app.command()
def profile(
    trajectory: Path = typer.Argument(..., help="Path to trajectory file (.xtc)"),
    topology: Path = typer.Argument(..., help="Path to topology file (.pdb or .tpr)"),
    n_frames: int = typer.Option(100, "--frames", "-f", help="Number of frames to use"),
    n_threads: int = typer.Option(1, "--threads", "-t", help="Number of threads"),
    algorithm: str = typer.Option(
        "sr", "--algorithm", "-a", help="Algorithm: sr or lr"
    ),
):
    """Profile time breakdown of SASA calculation.

    This profiles individual component times:
    - Coordinate extraction
    - Radii extraction
    - SASA calculation (which includes neighbor list building)
    """
    console.print(f"\n[bold]Loading trajectory:[/bold] {trajectory}")

    # Load trajectory
    traj = md.load(str(trajectory), top=str(topology))

    # Subsample
    if n_frames < traj.n_frames:
        step = traj.n_frames // n_frames
        traj = traj[::step][:n_frames]

    console.print(
        f"[bold]Frames:[/bold] {traj.n_frames}, [bold]Atoms:[/bold] {traj.n_atoms}"
    )

    # Time coordinate extraction
    start = time.perf_counter()
    coords = traj.xyz * 10.0  # nm -> Angstrom
    coords = np.ascontiguousarray(coords.astype(np.float32))
    coord_time = time.perf_counter() - start

    # Time radii extraction
    from zsasa.mdtraj import _get_radii_from_topology

    start = time.perf_counter()
    radii = _get_radii_from_topology(traj.topology)
    radii_time = time.perf_counter() - start

    # Time SASA calculation
    start = time.perf_counter()
    result = calculate_sasa_batch(
        coords, radii, n_threads=n_threads, algorithm=algorithm
    )
    sasa_time = time.perf_counter() - start

    total_time = coord_time + radii_time + sasa_time

    # Results
    table = Table(title="Time Breakdown")
    table.add_column("Component", style="cyan")
    table.add_column("Time (ms)", justify="right")
    table.add_column("Percentage", justify="right")

    table.add_row(
        "Coordinate extraction",
        f"{coord_time * 1000:.1f}",
        f"{coord_time / total_time * 100:.1f}%",
    )
    table.add_row(
        "Radii extraction",
        f"{radii_time * 1000:.1f}",
        f"{radii_time / total_time * 100:.1f}%",
    )
    table.add_row(
        "SASA calculation",
        f"{sasa_time * 1000:.1f}",
        f"{sasa_time / total_time * 100:.1f}%",
    )
    table.add_row(
        "[bold]Total[/bold]",
        f"[bold]{total_time * 1000:.1f}[/bold]",
        "[bold]100%[/bold]",
    )

    console.print(table)

    # Per-frame breakdown
    console.print(f"\n[bold]Per-frame timing:[/bold]")
    console.print(f"  Total:       {total_time / traj.n_frames * 1000:.3f} ms/frame")
    console.print(f"  SASA only:   {sasa_time / traj.n_frames * 1000:.3f} ms/frame")


if __name__ == "__main__":
    app()

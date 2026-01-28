#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich", "typer"]
# ///
"""Compare file-level vs atom-level parallelism.

Usage:
    # Run only file-level
    ./benchmarks/scripts/compare_parallelism.py file \
        --input-dir benchmarks/inputs \
        --sample-file benchmarks/samples/large_20k_10k.json \
        --threads 10 --limit 100

    # Run only atom-level (after cooldown)
    ./benchmarks/scripts/compare_parallelism.py atom \
        --input-dir benchmarks/inputs \
        --sample-file benchmarks/samples/large_20k_10k.json \
        --threads 10 --limit 100

    # Run both (back-to-back, less accurate)
    ./benchmarks/scripts/compare_parallelism.py both \
        --input-dir benchmarks/inputs \
        --sample-file benchmarks/samples/large_20k_10k.json \
        --threads 10 --limit 100
"""

from __future__ import annotations

import json
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Compare parallelism strategies")
console = Console()


def get_zig_binary() -> Path:
    root = Path(__file__).parent.parent.parent
    return root / "zig-out" / "bin" / "freesasa_zig"


def run_batch(
    input_dir: Path,
    files: list[str],
    threads: int,
    precision: str,
    parallelism: str,
) -> dict:
    """Run batch mode with specified parallelism strategy."""
    binary = get_zig_binary()

    with tempfile.TemporaryDirectory() as work_dir:
        work = Path(work_dir)
        batch_input = work / "input"
        batch_output = work / "output"
        batch_input.mkdir()
        batch_output.mkdir()

        # Create symlinks
        for f in files:
            for ext in [".json.gz", ".json"]:
                src = input_dir / f"{f}{ext}"
                if src.exists():
                    (batch_input / f"{f}{ext}").symlink_to(src.resolve())
                    break

        cmd = [
            str(binary),
            str(batch_input),
            str(batch_output),
            "--algorithm=sr",
            f"--threads={threads}",
            f"--precision={precision}",
            f"--parallelism={parallelism}",
            "--format=compact",
            "--timing",
            "--quiet",
        ]

        start = time.perf_counter()
        result = subprocess.run(cmd, capture_output=True, text=True)
        wall_time = (time.perf_counter() - start) * 1000

        output = result.stderr + result.stdout

        # Parse timing
        sasa_time = 0.0
        total_time = 0.0
        for line in output.split("\n"):
            if line.startswith("BATCH_SASA_TIME_MS:"):
                sasa_time = float(line.split(":")[1])
            elif line.startswith("BATCH_TOTAL_TIME_MS:"):
                total_time = float(line.split(":")[1])

        return {
            "sasa_time_ms": sasa_time,
            "total_time_ms": total_time,
            "wall_time_ms": wall_time,
        }


def load_files(sample_file: Path, limit: int) -> list[str]:
    """Load sample file list."""
    with open(sample_file) as f:
        data = json.load(f)
    return data.get("samples", [])[:limit]


def print_result(strategy: str, result: dict, n_files: int) -> None:
    """Print result for a single strategy."""
    throughput = n_files / (result["wall_time_ms"] / 1000)
    console.print(f"\n[bold]{strategy}[/bold]")
    console.print(f"  Wall time:  {result['wall_time_ms']:.0f}ms")
    console.print(f"  SASA time:  {result['sasa_time_ms']:.0f}ms")
    console.print(f"  Throughput: {throughput:.1f} files/s")


@app.command()
def file(
    input_dir: Annotated[Path, typer.Option("--input-dir", "-i")],
    sample_file: Annotated[Path, typer.Option("--sample-file", "-S")],
    threads: Annotated[int, typer.Option("--threads", "-t")] = 10,
    precision: Annotated[str, typer.Option("--precision", "-p")] = "f64",
    limit: Annotated[
        int, typer.Option("--limit", "-l", help="Limit number of files")
    ] = 100,
):
    """Run file-level parallelism only (N files × 1 thread)."""
    files = load_files(sample_file, limit)
    console.print(f"[bold]File-level Parallelism[/bold]")
    console.print(f"Files: {len(files)}, Threads: {threads}, Precision: {precision}")

    result = run_batch(input_dir, files, threads, precision, "file")
    print_result("File-level (N files × 1 thread)", result, len(files))


@app.command()
def atom(
    input_dir: Annotated[Path, typer.Option("--input-dir", "-i")],
    sample_file: Annotated[Path, typer.Option("--sample-file", "-S")],
    threads: Annotated[int, typer.Option("--threads", "-t")] = 10,
    precision: Annotated[str, typer.Option("--precision", "-p")] = "f64",
    limit: Annotated[
        int, typer.Option("--limit", "-l", help="Limit number of files")
    ] = 100,
):
    """Run atom-level parallelism only (1 file × N threads)."""
    files = load_files(sample_file, limit)
    console.print(f"[bold]Atom-level Parallelism[/bold]")
    console.print(f"Files: {len(files)}, Threads: {threads}, Precision: {precision}")

    result = run_batch(input_dir, files, threads, precision, "atom")
    print_result("Atom-level (1 file × N threads)", result, len(files))


@app.command()
def pipeline(
    input_dir: Annotated[Path, typer.Option("--input-dir", "-i")],
    sample_file: Annotated[Path, typer.Option("--sample-file", "-S")],
    threads: Annotated[int, typer.Option("--threads", "-t")] = 10,
    precision: Annotated[str, typer.Option("--precision", "-p")] = "f64",
    limit: Annotated[
        int, typer.Option("--limit", "-l", help="Limit number of files")
    ] = 100,
):
    """Run pipelined parallelism (I/O prefetch + atom-level SASA)."""
    files = load_files(sample_file, limit)
    console.print(f"[bold]Pipelined Parallelism[/bold]")
    console.print(f"Files: {len(files)}, Threads: {threads}, Precision: {precision}")

    result = run_batch(input_dir, files, threads, precision, "pipeline")
    print_result("Pipelined (I/O prefetch + N threads)", result, len(files))


@app.command()
def both(
    input_dir: Annotated[Path, typer.Option("--input-dir", "-i")],
    sample_file: Annotated[Path, typer.Option("--sample-file", "-S")],
    threads: Annotated[int, typer.Option("--threads", "-t")] = 10,
    precision: Annotated[str, typer.Option("--precision", "-p")] = "f64",
    limit: Annotated[
        int, typer.Option("--limit", "-l", help="Limit number of files")
    ] = 100,
    cooldown: Annotated[
        int, typer.Option("--cooldown", "-c", help="Cooldown seconds between runs")
    ] = 0,
):
    """Run both strategies (back-to-back comparison)."""
    files = load_files(sample_file, limit)
    console.print("[bold]Parallelism Comparison[/bold]")
    console.print(f"Files: {len(files)}, Threads: {threads}, Precision: {precision}")
    if cooldown > 0:
        console.print(f"Cooldown: {cooldown}s between runs")
    console.print()

    # Run file-level parallelism
    console.print("[cyan]Running file-level parallelism...[/cyan]")
    file_result = run_batch(input_dir, files, threads, precision, "file")
    console.print(f"  Wall time: {file_result['wall_time_ms']:.0f}ms")

    # Cooldown
    if cooldown > 0:
        console.print(f"[yellow]Cooling down for {cooldown}s...[/yellow]")
        time.sleep(cooldown)

    # Run atom-level parallelism
    console.print("[cyan]Running atom-level parallelism...[/cyan]")
    atom_result = run_batch(input_dir, files, threads, precision, "atom")
    console.print(f"  Wall time: {atom_result['wall_time_ms']:.0f}ms")

    # Summary
    console.print()
    table = Table(title="Results")
    table.add_column("Strategy")
    table.add_column("Wall Time", justify="right")
    table.add_column("SASA Time", justify="right")
    table.add_column("Throughput", justify="right")

    file_tp = len(files) / (file_result["wall_time_ms"] / 1000)
    atom_tp = len(files) / (atom_result["wall_time_ms"] / 1000)

    table.add_row(
        "File-level (N files × 1 thread)",
        f"{file_result['wall_time_ms']:.0f}ms",
        f"{file_result['sasa_time_ms']:.0f}ms",
        f"{file_tp:.1f}/s",
    )
    table.add_row(
        "Atom-level (1 file × N threads)",
        f"{atom_result['wall_time_ms']:.0f}ms",
        f"{atom_result['sasa_time_ms']:.0f}ms",
        f"{atom_tp:.1f}/s",
    )

    console.print(table)

    # Winner
    if file_result["wall_time_ms"] < atom_result["wall_time_ms"]:
        speedup = atom_result["wall_time_ms"] / file_result["wall_time_ms"]
        console.print(f"\n[green]File-level is {speedup:.2f}x faster[/green]")
    else:
        speedup = file_result["wall_time_ms"] / atom_result["wall_time_ms"]
        console.print(f"\n[green]Atom-level is {speedup:.2f}x faster[/green]")


if __name__ == "__main__":
    app()

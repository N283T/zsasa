#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Batch SASA benchmark runner.

Runs benchmark using batch processing mode for each tool.
Unlike run.py which processes files one-by-one, this uses native batch modes:
  - Zig: directory input auto-detected
  - Rust: --json-dir flag
  - FreeSASA: freesasa_batch.sh wrapper

Usage:
    # Default dataset (9 structures)
    ./benchmarks/scripts/run_batch.py --tool zig --algorithm sr --threads 1,4,8

    # With stratified sample
    ./benchmarks/scripts/run_batch.py --tool zig --algorithm sr \
        --input-dir benchmarks/inputs \
        --sample-file benchmarks/samples/stratified_1k.json \
        --threads 1,4,8

Output:
    benchmarks/results/batch_{tool}_{algorithm}/
    ├── config.json   # System info and parameters
    └── results.csv   # Benchmark results (per-run timing)
"""

from __future__ import annotations

import csv
import json
import os
import platform
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn
from rich.table import Table

app = typer.Typer(help="Batch SASA benchmark runner")
console = Console()

TOOLS = ["zig", "rust", "freesasa"]
ALGORITHMS = ["sr", "lr"]


def get_system_info() -> dict:
    """Get system information for benchmark metadata."""
    return {
        "platform": platform.platform(),
        "processor": platform.processor(),
        "python_version": platform.python_version(),
        "cpu_count": os.cpu_count(),
    }


def get_binary_path(tool: str) -> Path:
    """Get path to tool binary."""
    root = Path(__file__).parent.parent.parent
    if tool == "zig":
        return root / "zig-out" / "bin" / "freesasa_zig"
    elif tool == "freesasa":
        return root / "benchmarks" / "scripts" / "freesasa_batch.sh"
    elif tool == "rust":
        return (
            root
            / "benchmarks"
            / "external"
            / "rustsasa-bench"
            / "target"
            / "release"
            / "rust-sasa"
        )
    raise ValueError(f"Unknown tool: {tool}")


def parse_threads(threads_str: str) -> list[int]:
    """Parse thread specification like '1-10' or '1,4,8'."""
    if "-" in threads_str and "," not in threads_str:
        start, end = map(int, threads_str.split("-"))
        return list(range(start, end + 1))
    return [int(x.strip()) for x in threads_str.split(",")]


def load_sample_file(sample_path: Path) -> list[str]:
    """Load sample IDs from sample file."""
    with open(sample_path) as f:
        data = json.load(f)
    return data.get("samples", [])


def prepare_input_dir(
    input_dir: Path,
    sample_ids: set[str] | None,
    work_dir: Path,
) -> tuple[Path, int]:
    """Prepare input directory with sampled files.

    If sample_ids is None, returns original input_dir.
    Otherwise, creates symlinks in work_dir for sampled files.

    Returns (input_path, file_count).
    """
    if sample_ids is None:
        # Count files in input_dir
        count = sum(
            1
            for f in os.scandir(input_dir)
            if f.is_file() and (f.name.endswith(".json.gz") or f.name.endswith(".json"))
        )
        return input_dir, count

    # Create symlinks for sampled files
    batch_input = work_dir / "input"
    batch_input.mkdir(parents=True, exist_ok=True)

    count = 0
    for sample_id in sample_ids:
        # Try .json.gz first, then .json
        for ext in [".json.gz", ".json"]:
            src = input_dir / f"{sample_id}{ext}"
            if src.exists():
                dst = batch_input / f"{sample_id}{ext}"
                # Use absolute path for symlink to work from any directory
                dst.symlink_to(src.resolve())
                count += 1
                break

    return batch_input, count


def run_zig_batch(
    input_dir: Path,
    output_dir: Path,
    algorithm: str,
    n_threads: int,
) -> dict:
    """Run Zig batch mode."""
    binary = get_binary_path("zig")
    if not binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {binary}")

    cmd = [
        str(binary),
        str(input_dir),
        str(output_dir),
        f"--algorithm={algorithm}",
        f"--threads={n_threads}",
        "--format=compact",
        "--timing",
        "--quiet",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stderr + result.stdout

    # Parse BATCH_* output
    return parse_batch_output(output)


def run_rust_batch(
    input_dir: Path,
    output_dir: Path,
    n_threads: int,
) -> dict:
    """Run Rust batch mode."""
    binary = get_binary_path("rust")
    if not binary.exists():
        raise FileNotFoundError(f"Rust binary not found: {binary}")

    cmd = [
        str(binary),
        f"--json-dir={input_dir}",
        f"--output-dir={output_dir}",
        f"-t={n_threads}",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stderr + result.stdout

    return parse_batch_output(output)


def run_freesasa_batch(
    input_dir: Path,
    output_dir: Path,
    n_jobs: int,
) -> dict:
    """Run FreeSASA batch mode via shell script."""
    script = get_binary_path("freesasa")
    if not script.exists():
        raise FileNotFoundError(f"FreeSASA batch script not found: {script}")

    cmd = [str(script), str(input_dir), str(output_dir), str(n_jobs)]

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stderr + result.stdout

    return parse_batch_output(output)


def parse_batch_output(output: str) -> dict:
    """Parse BATCH_* output from any tool."""
    result = {
        "sasa_time_ms": 0.0,
        "total_time_ms": 0.0,
        "files": 0,
        "successful": 0,
        "failed": 0,
    }

    for line in output.split("\n"):
        if line.startswith("BATCH_SASA_TIME_MS:"):
            result["sasa_time_ms"] = float(line.split(":")[1])
        elif line.startswith("BATCH_TOTAL_TIME_MS:"):
            result["total_time_ms"] = float(line.split(":")[1])
        elif line.startswith("BATCH_FILES:"):
            result["files"] = int(line.split(":")[1])
        elif line.startswith("BATCH_SUCCESS:") or line.startswith("BATCH_SUCCESSFUL:"):
            result["successful"] = int(line.split(":")[1])
        elif line.startswith("BATCH_FAILED:"):
            result["failed"] = int(line.split(":")[1])

    return result


def run_batch_benchmark(
    tool: str,
    input_dir: Path,
    output_dir: Path,
    algorithm: str,
    n_threads: int,
) -> dict:
    """Run batch benchmark for a specific tool."""
    if tool == "zig":
        return run_zig_batch(input_dir, output_dir, algorithm, n_threads)
    elif tool == "rust":
        return run_rust_batch(input_dir, output_dir, n_threads)
    elif tool == "freesasa":
        return run_freesasa_batch(input_dir, output_dir, n_threads)
    raise ValueError(f"Unknown tool: {tool}")


@app.command()
def main(
    tool: Annotated[
        str, typer.Option("--tool", "-t", help="Tool: zig, freesasa, rust")
    ],
    algorithm: Annotated[
        str, typer.Option("--algorithm", "-a", help="Algorithm: sr, lr")
    ] = "sr",
    threads: Annotated[
        str, typer.Option("--threads", "-T", help="Thread counts: '1-8' or '1,4,8'")
    ] = "1,4,8",
    runs: Annotated[
        int, typer.Option("--runs", "-r", help="Number of runs per configuration")
    ] = 3,
    output_dir: Annotated[
        Path | None, typer.Option("--output-dir", "-o", help="Output directory")
    ] = None,
    input_dir: Annotated[
        Path | None,
        typer.Option("--input-dir", "-i", help="Input directory with .json.gz files"),
    ] = None,
    sample_file: Annotated[
        Path | None,
        typer.Option("--sample-file", "-S", help="Sample file from sample.py"),
    ] = None,
):
    """Run batch SASA benchmark."""
    thread_counts = parse_threads(threads)

    # Validate tool
    if tool not in TOOLS:
        console.print(f"[red]Error:[/red] Unknown tool: {tool}")
        console.print(f"Available: {', '.join(TOOLS)}")
        raise typer.Exit(1)

    # Validate algorithm
    if algorithm not in ALGORITHMS:
        console.print(f"[red]Error:[/red] Unknown algorithm: {algorithm}")
        console.print(f"Available: {', '.join(ALGORITHMS)}")
        raise typer.Exit(1)

    # RustSASA only supports SR
    if tool == "rust" and algorithm == "lr":
        console.print("[red]Error:[/red] RustSASA only supports SR algorithm")
        raise typer.Exit(1)

    # Validate sample file requires input-dir
    sample_ids: set[str] | None = None
    if sample_file is not None:
        if input_dir is None:
            console.print("[red]Error:[/red] --sample-file requires --input-dir")
            raise typer.Exit(1)
        if not sample_file.exists():
            console.print(f"[red]Error:[/red] Sample file not found: {sample_file}")
            raise typer.Exit(1)
        sample_ids = set(load_sample_file(sample_file))
        console.print(
            f"Loaded [cyan]{len(sample_ids):,}[/cyan] samples from {sample_file}"
        )

    # Setup input directory
    if input_dir is None:
        input_dir = Path(__file__).parent.parent / "dataset"

    if not input_dir.exists():
        console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
        raise typer.Exit(1)

    # Setup output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    if output_dir is None:
        output_dir = (
            Path(__file__).parent.parent / "results" / f"batch_{tool}_{algorithm}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create work directory for batch output
    work_dir = Path(tempfile.mkdtemp(prefix="sasa_batch_"))

    try:
        # Prepare input (create symlinks if sampling)
        batch_input, file_count = prepare_input_dir(input_dir, sample_ids, work_dir)

        console.print(f"[bold]{tool.upper()} {algorithm.upper()} (Batch Mode)[/bold]")
        console.print(f"Files: {file_count:,}, Threads: {thread_counts}, Runs: {runs}")
        console.print()

        # Save config
        config = {
            "timestamp": timestamp,
            "system": get_system_info(),
            "parameters": {
                "tool": tool,
                "algorithm": algorithm,
                "thread_counts": thread_counts,
                "runs": runs,
                "file_count": file_count,
                "input_dir": str(input_dir),
                "sample_file": str(sample_file) if sample_file else None,
                "mode": "batch",
            },
        }
        config_path = output_dir / "config.json"
        config_path.write_text(json.dumps(config, indent=2))

        # Run benchmarks
        csv_path = output_dir / "results.csv"
        results = []

        total_runs = len(thread_counts) * runs

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running", total=total_runs)

            for n_threads in thread_counts:
                for run_num in range(1, runs + 1):
                    progress.update(task, description=f"t={n_threads} run={run_num}")

                    # Create output dir for this run
                    batch_output = work_dir / f"output_t{n_threads}_r{run_num}"
                    batch_output.mkdir(parents=True, exist_ok=True)

                    try:
                        result = run_batch_benchmark(
                            tool, batch_input, batch_output, algorithm, n_threads
                        )
                        results.append(
                            {
                                "tool": tool,
                                "algorithm": algorithm,
                                "threads": n_threads,
                                "run": run_num,
                                "files": result["files"],
                                "successful": result["successful"],
                                "failed": result["failed"],
                                "sasa_time_ms": result["sasa_time_ms"],
                                "total_time_ms": result["total_time_ms"],
                            }
                        )
                    except Exception as e:
                        console.print(f"[red]Error:[/red] {e}")
                        results.append(
                            {
                                "tool": tool,
                                "algorithm": algorithm,
                                "threads": n_threads,
                                "run": run_num,
                                "files": 0,
                                "successful": 0,
                                "failed": 0,
                                "sasa_time_ms": 0,
                                "total_time_ms": 0,
                            }
                        )

                    progress.advance(task)

                    # Clean up output dir
                    shutil.rmtree(batch_output, ignore_errors=True)

        # Write CSV
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "tool",
                    "algorithm",
                    "threads",
                    "run",
                    "files",
                    "successful",
                    "failed",
                    "sasa_time_ms",
                    "total_time_ms",
                ],
            )
            writer.writeheader()
            writer.writerows(results)

        console.print(f"\n[green]Results saved to {csv_path}[/green]")

        # Print summary table
        table = Table(title="Batch Benchmark Summary")
        table.add_column("Threads", justify="right")
        table.add_column("SASA Time (ms)", justify="right")
        table.add_column("Total Time (ms)", justify="right")
        table.add_column("Files/sec", justify="right")

        # Group by threads and calculate averages
        from collections import defaultdict

        by_threads = defaultdict(list)
        for r in results:
            by_threads[r["threads"]].append(r)

        for n_threads in sorted(by_threads.keys()):
            thread_results = by_threads[n_threads]
            avg_sasa = sum(r["sasa_time_ms"] for r in thread_results) / len(
                thread_results
            )
            avg_total = sum(r["total_time_ms"] for r in thread_results) / len(
                thread_results
            )
            files = thread_results[0]["files"]
            throughput = files / (avg_total / 1000) if avg_total > 0 else 0

            table.add_row(
                str(n_threads),
                f"{avg_sasa:.2f}",
                f"{avg_total:.2f}",
                f"{throughput:.1f}",
            )

        console.print(table)

    finally:
        # Cleanup work directory
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    app()

#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""SASA benchmark runner.

Runs benchmark for a single tool and algorithm with configurable thread counts.
Results are saved to CSV with execution config in JSON.

Usage:
    ./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10
    ./benchmarks/scripts/run.py --tool freesasa --algorithm lr --threads 1-10
    ./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1 --runs 1

Output:
    benchmarks/results/{tool}_{algorithm}/
    ├── config.json   # System info and parameters
    └── results.csv   # Benchmark results
"""

from __future__ import annotations

import csv
import gzip
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

app = typer.Typer(help="SASA benchmark with separated tool execution")
console = Console()

# All benchmark structures
STRUCTURES = [
    ("1crn", 327, "Crambin"),
    ("1ubq", 602, "Ubiquitin"),
    ("1a0q", 3183, "Lipid transfer protein"),
    ("3hhb", 4384, "Hemoglobin"),
    ("1aon", 58674, "GroEL-GroES"),
    ("4v6x", 237685, "Ribosome"),
]

TOOLS = ["zig", "freesasa", "rust"]
ALGORITHMS = ["sr", "lr"]


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


def parse_threads(threads_str: str) -> list[int]:
    """Parse thread specification like '1-10' or '1,4,8'."""
    result = []
    for part in threads_str.split(","):
        if "-" in part:
            start, end = part.split("-", 1)
            result.extend(range(int(start), int(end) + 1))
        else:
            result.append(int(part))
    return sorted(set(result))


def get_binary_path(tool: str) -> Path:
    """Get binary path for a tool."""
    # Project root (freesasa-zig/)
    root = Path(__file__).parent.parent.parent

    if tool == "zig":
        return root / "zig-out" / "bin" / "freesasa_zig"
    elif tool == "freesasa":
        return root / "benchmarks" / "external" / "freesasa-bench" / "src" / "freesasa"
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
    else:
        raise ValueError(f"Unknown tool: {tool}")


def run_zig(
    json_path: Path,
    algorithm: str,
    n_threads: int,
) -> tuple[float, float]:
    """Run Zig benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("zig")
    if not binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {binary}")

    input_path: Path | None = None
    output_path: Path | None = None
    cleanup_input = False

    try:
        # Decompress if needed
        if str(json_path).endswith(".gz"):
            with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
                input_path = Path(f.name)
            cleanup_input = True
            with gzip.open(json_path, "rb") as f_in:
                with open(input_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            input_path = json_path

        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            output_path = Path(f.name)

        cmd = [
            str(binary),
            "--timing",
            f"--algorithm={algorithm}",
            f"--threads={n_threads}",
            str(input_path),
            str(output_path),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        # Parse timing
        sasa_time_ms = 0.0
        for line in result.stderr.split("\n"):
            match = re.search(r"SASA calculation:\s*([\d.]+)\s*ms", line)
            if match:
                sasa_time_ms = float(match.group(1))
                break

        # Parse output
        with open(output_path) as f:
            data = json.load(f)

        return sasa_time_ms, data["total_area"]

    finally:
        if output_path is not None:
            output_path.unlink(missing_ok=True)
        if cleanup_input and input_path is not None:
            input_path.unlink(missing_ok=True)


def run_freesasa(
    json_path: Path,
    algorithm: str,
    n_threads: int,
    n_points: int = 100,
    n_slices: int = 20,
) -> tuple[float, float]:
    """Run FreeSASA C benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("freesasa")
    if not binary.exists():
        raise FileNotFoundError(f"FreeSASA binary not found: {binary}")

    cmd = [
        str(binary),
        "--json-input",
        str(json_path),
        f"--n-threads={n_threads}",
    ]

    if algorithm == "sr":
        cmd.extend(["--shrake-rupley", f"--resolution={n_points}"])
    else:
        cmd.extend(["--lee-richards", f"--resolution={n_slices}"])

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"FreeSASA failed: {result.stderr}")

    # Parse total SASA
    total_sasa = 0.0
    for line in result.stdout.split("\n"):
        if line.startswith("Total"):
            match = re.search(r"[\d.]+", line)
            if match:
                total_sasa = float(match.group())
                break

    # Parse timing
    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        match = re.search(r"SASA calculation time:\s*([\d.]+)\s*ms", line)
        if match:
            sasa_time_ms = float(match.group(1))
            break

    return sasa_time_ms, total_sasa


def run_rust(
    json_path: Path,
    n_threads: int,
    n_points: int = 100,
) -> tuple[float, float]:
    """Run RustSASA benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("rust")
    if not binary.exists():
        raise FileNotFoundError(f"RustSASA binary not found: {binary}")

    cmd = [
        str(binary),
        "-J",
        str(json_path),
        "-n",
        str(n_points),
        "-t",
        str(n_threads),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        raise RuntimeError(f"RustSASA failed: {result.stderr}")

    # Parse timing
    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        if "SASA_TIME_US:" in line:
            try:
                us = int(line.split(":")[1].strip())
                sasa_time_ms = us / 1000.0
            except (ValueError, IndexError):
                pass
            break

    # Parse total SASA
    total_sasa = 0.0
    for line in result.stdout.split("\n"):
        match = re.search(r"Total SASA:\s*([\d.]+)", line)
        if match:
            total_sasa = float(match.group(1))
            break

    return sasa_time_ms, total_sasa


def run_benchmark(
    tool: str,
    json_path: Path,
    algorithm: str,
    n_threads: int,
) -> tuple[float, float]:
    """Run benchmark for a specific tool. Returns (sasa_time_ms, total_sasa)."""
    if tool == "zig":
        return run_zig(json_path, algorithm, n_threads)
    elif tool == "freesasa":
        return run_freesasa(json_path, algorithm, n_threads)
    elif tool == "rust":
        if algorithm != "sr":
            raise ValueError("RustSASA only supports SR algorithm")
        return run_rust(json_path, n_threads)
    else:
        raise ValueError(f"Unknown tool: {tool}")


@app.command()
def main(
    tool: Annotated[
        str,
        typer.Option("--tool", "-t", help="Tool: zig, freesasa, rust"),
    ],
    algorithm: Annotated[
        str,
        typer.Option("--algorithm", "-a", help="Algorithm: sr, lr"),
    ],
    threads: Annotated[
        str,
        typer.Option("--threads", "-T", help="Thread counts: '1-10' or '1,4,8'"),
    ] = "1-10",
    runs: Annotated[
        int,
        typer.Option("--runs", "-r", help="Number of runs per configuration"),
    ] = 3,
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Output directory"),
    ] = None,
) -> None:
    """Run SASA benchmark for a single tool and algorithm."""

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

    # Setup output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    if output_dir is None:
        output_dir = Path(__file__).parent.parent / "results" / f"{tool}_{algorithm}"
    output_dir.mkdir(parents=True, exist_ok=True)

    json_dir = Path(__file__).parent.parent / "dataset"

    # Save config
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "tool": tool,
            "algorithm": algorithm,
            "thread_counts": thread_counts,
            "runs": runs,
            "structures": [s[0] for s in STRUCTURES],
        },
    }
    config_path = output_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2))

    # Calculate total runs
    total_runs = len(STRUCTURES) * len(thread_counts) * runs

    console.print(f"[bold]{tool.upper()} {algorithm.upper()}[/bold]")
    console.print(f"Threads: {thread_counts}, Runs: {runs}, Total: {total_runs}\n")

    # Run benchmarks
    csv_path = output_dir / "results.csv"

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "tool",
                "structure",
                "n_atoms",
                "algorithm",
                "threads",
                "run",
                "sasa_time_ms",
                "total_sasa",
            ],
        )
        writer.writeheader()

        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            console=console,
        ) as progress:
            task = progress.add_task("Running", total=total_runs)

            for n_threads in thread_counts:
                for pdb_id, n_atoms, _ in STRUCTURES:
                    json_path = json_dir / f"{pdb_id}.json.gz"
                    if not json_path.exists():
                        console.print(f"[yellow]Skip {pdb_id}: not found[/yellow]")
                        continue

                    for run_num in range(1, runs + 1):
                        desc = f"{pdb_id} t={n_threads} #{run_num}"
                        progress.update(task, description=desc)

                        try:
                            sasa_time, total_sasa = run_benchmark(
                                tool, json_path, algorithm, n_threads
                            )

                            writer.writerow(
                                {
                                    "tool": tool,
                                    "structure": pdb_id,
                                    "n_atoms": n_atoms,
                                    "algorithm": algorithm,
                                    "threads": n_threads,
                                    "run": run_num,
                                    "sasa_time_ms": sasa_time,
                                    "total_sasa": total_sasa,
                                }
                            )
                            f.flush()

                        except Exception as e:
                            err = f"{pdb_id} t={n_threads}"
                            console.print(f"[red]Error: {err}: {e}[/red]")

                        progress.advance(task)

    console.print(f"\n[green]Done![/green] Results saved to: {output_dir}")
    console.print(f"  - {csv_path.name}")
    console.print(f"  - {config_path.name}")


if __name__ == "__main__":
    app()

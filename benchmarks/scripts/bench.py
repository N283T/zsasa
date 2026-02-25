#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""SASA benchmark runner.

Runs benchmark for a single tool and algorithm with configurable thread counts.
Uses internal timing (SASA-only) with warmup runs for statistical reliability.
Results are saved to CSV with execution config in JSON.

Tools: zig_f64 (default), zig_f32, freesasa, rust (zig = zig_f64 alias)

Usage:
    # Default dataset (benchmarks/dataset/json/, ~2k structures)
    ./benchmarks/scripts/bench.py --tool zig_f64 --algorithm sr
    ./benchmarks/scripts/bench.py --tool zig_f32 --algorithm sr
    ./benchmarks/scripts/bench.py --tool freesasa --algorithm sr

    # "zig" is shorthand for "zig_f64"
    ./benchmarks/scripts/bench.py --tool zig --algorithm sr

    # Quick test (no warmup, 1 run)
    ./benchmarks/scripts/bench.py --tool zig --algorithm sr --threads 1 --warmup 0 --runs 1

Output:
    benchmarks/results/{tool}_{algorithm}/
    ├── config.json   # System info and parameters
    └── results.csv   # Benchmark results (warmup runs excluded)

    Examples: zig_f64_sr/, zig_f32_sr/, freesasa_sr/, rust_sr/
"""

from __future__ import annotations

import csv
import gzip
import json
import math
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

app = typer.Typer(help="SASA benchmark runner")
console = Console()

TOOLS = ["zig_f64", "zig_f32", "freesasa", "rust"]
TOOL_ALIASES = {"zig": "zig_f64"}
ALGORITHMS = ["sr", "lr"]


def parse_tool(tool: str) -> tuple[str, str, str]:
    """Parse tool name into (canonical, base, precision).

    Examples:
        "zig_f64" -> ("zig_f64", "zig", "f64")
        "zig_f32" -> ("zig_f32", "zig", "f32")
        "zig"     -> ("zig_f64", "zig", "f64")  # alias
        "freesasa" -> ("freesasa", "freesasa", "f64")
        "rust"     -> ("rust", "rust", "f64")
    """
    tool = TOOL_ALIASES.get(tool, tool)
    if tool.startswith("zig_f"):
        return tool, "zig", tool.split("_")[1]
    return tool, tool, "f64"


def get_n_atoms_from_json(json_path: Path) -> int:
    """Get number of atoms from a JSON file."""
    try:
        if str(json_path).endswith(".gz"):
            with gzip.open(json_path, "rt") as f:
                data = json.load(f)
        else:
            with open(json_path) as f:
                data = json.load(f)
        return len(data.get("x", []))
    except Exception:
        return 0


def scan_input_directory(input_dir: Path) -> list[tuple[str, int]]:
    """Scan directory for .json.gz files and return (id, n_atoms=0) list.

    Uses os.scandir for fast scanning. n_atoms is resolved lazily during run.
    """
    entries = []
    with os.scandir(input_dir) as it:
        for entry in it:
            if entry.is_file() and (
                entry.name.endswith(".json.gz") or entry.name.endswith(".json")
            ):
                entries.append(entry.name)

    entries.sort()

    structures = []
    for filename in entries:
        if filename.endswith(".json.gz"):
            pdb_id = filename[:-8]
        else:
            pdb_id = filename[:-5]
        structures.append((pdb_id, 0))

    return structures


def load_sample_file(sample_path: Path) -> list[str]:
    """Load sample file and return list of IDs.

    Supports both formats:
    - v1: {"samples": ["id1", "id2", ...]}
    - v2: {"samples": {"bin": [{"id": "...", ...}, ...], ...}}
    """
    with open(sample_path) as f:
        data = json.load(f)

    if "samples" not in data:
        raise ValueError(f"Invalid sample file: missing 'samples' key in {sample_path}")

    samples = data["samples"]

    # v1: flat list
    if isinstance(samples, list):
        return samples

    # v2: dict of bins
    ids = []
    for entries in samples.values():
        for entry in entries:
            ids.append(entry["id"])
    return ids


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
    root = Path(__file__).parent.parent.parent

    if tool == "zig":
        return root.joinpath("zig-out", "bin", "zsasa")
    elif tool == "freesasa":
        return root.joinpath(
            "benchmarks", "external", "freesasa-bench", "src", "freesasa"
        )
    elif tool == "rust":
        return root.joinpath(
            "benchmarks",
            "external",
            "rustsasa-bench",
            "target",
            "release",
            "rust-sasa",
        )
    else:
        raise ValueError(f"Unknown tool: {tool}")


def run_zig(
    json_path: Path,
    algorithm: str,
    n_threads: int,
    precision: str = "f64",
    n_points: int = 100,
    use_bitmask: bool = False,
) -> tuple[float, float]:
    """Run Zig benchmark. Returns (sasa_time_ms, total_sasa)."""
    binary = get_binary_path("zig")
    if not binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {binary}")

    input_path: Path | None = None
    output_path: Path | None = None
    cleanup_input = False

    try:
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
            "calc",
            "--timing",
            f"--algorithm={algorithm}",
            f"--threads={n_threads}",
            f"--precision={precision}",
            f"--n-points={n_points}",
            str(input_path),
            str(output_path),
        ]
        if use_bitmask:
            cmd.insert(-2, "--use-bitmask")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        sasa_time_ms = 0.0
        for line in result.stderr.split("\n"):
            match = re.search(r"SASA calculation:\s*([\d.]+)\s*ms", line)
            if match:
                sasa_time_ms = float(match.group(1))
                break

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

    total_sasa = 0.0
    for line in result.stdout.split("\n"):
        if line.startswith("Total"):
            match = re.search(r"[\d.]+", line)
            if match:
                total_sasa = float(match.group())
                break

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

    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        if "SASA_TIME_US:" in line:
            try:
                us = int(line.split(":")[1].strip())
                sasa_time_ms = us / 1000.0
            except (ValueError, IndexError):
                pass
            break

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
    precision: str = "f64",
    n_points: int = 100,
    use_bitmask: bool = False,
) -> tuple[float, float]:
    """Run benchmark for a specific tool. Returns (sasa_time_ms, total_sasa)."""
    if tool == "zig":
        return run_zig(
            json_path, algorithm, n_threads, precision, n_points, use_bitmask
        )
    elif tool == "freesasa":
        return run_freesasa(json_path, algorithm, n_threads, n_points)
    elif tool == "rust":
        if algorithm != "sr":
            raise ValueError("RustSASA only supports SR algorithm")
        return run_rust(json_path, n_threads, n_points)
    else:
        raise ValueError(f"Unknown tool: {tool}")


def print_summary(csv_path: Path, warmup: int, runs: int) -> None:
    """Print summary statistics from results CSV."""
    # Aggregate: per (threads, structure) → list of sasa_time_ms
    from collections import defaultdict

    by_threads: dict[int, list[float]] = defaultdict(list)

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            t = int(row["threads"])
            by_threads[t].append(float(row["sasa_time_ms"]))

    if not by_threads:
        return

    table = Table(title="Summary (SASA time per structure)")
    table.add_column("Threads", style="cyan", justify="right")
    table.add_column("N", justify="right")
    table.add_column("Mean (ms)", justify="right")
    table.add_column("Std (ms)", justify="right")
    table.add_column("Min (ms)", justify="right")
    table.add_column("Max (ms)", justify="right")
    table.add_column("Total (s)", justify="right")

    for t in sorted(by_threads):
        times = by_threads[t]
        n = len(times)
        mean = sum(times) / n
        variance = sum((x - mean) ** 2 for x in times) / n if n > 1 else 0.0
        std = math.sqrt(variance)
        total_s = sum(times) / 1000.0
        table.add_row(
            str(t),
            f"{n:,}",
            f"{mean:.3f}",
            f"{std:.3f}",
            f"{min(times):.3f}",
            f"{max(times):.3f}",
            f"{total_s:.1f}",
        )

    console.print()
    console.print(table)


@app.command()
def main(
    tool: Annotated[
        str,
        typer.Option(
            "--tool",
            "-t",
            help="Tool: zig_f64, zig_f32, freesasa, rust (zig = zig_f64)",
        ),
    ],
    algorithm: Annotated[
        str,
        typer.Option("--algorithm", "-a", help="Algorithm: sr, lr"),
    ],
    threads: Annotated[
        str,
        typer.Option("--threads", "-T", help="Thread counts: '1,4,10' or '1-10'"),
    ] = "1,4,10",
    runs: Annotated[
        int,
        typer.Option("--runs", "-r", help="Number of measured runs per configuration"),
    ] = 5,
    warmup: Annotated[
        int,
        typer.Option("--warmup", "-w", help="Warmup runs (not recorded)"),
    ] = 1,
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Output directory"),
    ] = None,
    input_dir: Annotated[
        Path | None,
        typer.Option("--input-dir", "-i", help="Input directory with .json.gz files"),
    ] = None,
    sample_file: Annotated[
        Path | None,
        typer.Option(
            "--sample-file",
            "-S",
            help="Sample file to filter structures (v1 or v2 format)",
        ),
    ] = None,
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
    """Run SASA benchmark for a single tool and algorithm."""

    thread_counts = parse_threads(threads)

    # Parse and validate tool
    all_valid = list(TOOLS) + list(TOOL_ALIASES.keys())
    if tool not in all_valid:
        console.print(f"[red]Error:[/red] Unknown tool: {tool}")
        console.print(f"Available: {', '.join(TOOLS)} (zig = zig_f64)")
        raise typer.Exit(1)

    tool_canonical, tool_base, precision = parse_tool(tool)

    # Validate algorithm
    if algorithm not in ALGORITHMS:
        console.print(f"[red]Error:[/red] Unknown algorithm: {algorithm}")
        console.print(f"Available: {', '.join(ALGORITHMS)}")
        raise typer.Exit(1)

    if tool_base == "rust" and algorithm == "lr":
        console.print("[red]Error:[/red] RustSASA only supports SR algorithm")
        raise typer.Exit(1)

    # Load sample filter
    sample_ids: set[str] | None = None
    if sample_file is not None:
        if not sample_file.exists():
            console.print(f"[red]Error:[/red] Sample file not found: {sample_file}")
            raise typer.Exit(1)
        sample_ids = set(load_sample_file(sample_file))
        console.print(
            f"Loaded [cyan]{len(sample_ids):,}[/cyan] samples from {sample_file}"
        )

    # Setup input directory
    if input_dir is not None:
        if not input_dir.exists():
            console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
            raise typer.Exit(1)
        json_dir = input_dir
    else:
        json_dir = Path(__file__).parent.parent.joinpath("dataset", "json")
        if not json_dir.exists():
            console.print(f"[red]Error:[/red] Default dataset not found: {json_dir}")
            console.print("Run sampling first or specify --input-dir")
            raise typer.Exit(1)

    structures = scan_input_directory(json_dir)
    if not structures:
        console.print(f"[red]Error:[/red] No .json.gz files found in {json_dir}")
        raise typer.Exit(1)

    # Filter by sample file if provided
    if sample_ids is not None:
        structures = [(pdb_id, n) for pdb_id, n in structures if pdb_id in sample_ids]
        if not structures:
            console.print("[red]Error:[/red] No matching structures found")
            raise typer.Exit(1)

    console.print(f"Found [cyan]{len(structures):,}[/cyan] structures in {json_dir}")

    # Setup output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    if output_dir is None:
        output_dir = Path(__file__).parent.parent.joinpath(
            "results", f"{tool_canonical}_{algorithm}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save config
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "tool": tool_canonical,
            "tool_base": tool_base,
            "algorithm": algorithm,
            "precision": precision,
            "thread_counts": thread_counts,
            "warmup": warmup,
            "runs": runs,
            "n_points": n_points,
            "use_bitmask": use_bitmask,
            "n_structures": len(structures),
            "input_dir": str(json_dir),
            "sample_file": str(sample_file) if sample_file else None,
        },
    }
    config_path = output_dir.joinpath("config.json")
    config_path.write_text(json.dumps(config, indent=2))

    # Calculate totals
    total_measured = len(structures) * len(thread_counts) * runs
    total_warmup = len(structures) * len(thread_counts) * warmup
    total_all = total_measured + total_warmup

    console.print(f"\n[bold]{tool_canonical.upper()} {algorithm.upper()}[/bold]")
    console.print(
        f"Threads: {thread_counts}, "
        f"Warmup: {warmup}, Runs: {runs}, "
        f"Total: {total_all:,} ({total_warmup:,} warmup + {total_measured:,} measured)\n"
    )

    # Run benchmarks
    csv_path = output_dir.joinpath("results.csv")

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "tool",
                "structure",
                "n_atoms",
                "algorithm",
                "precision",
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
            task = progress.add_task("Running", total=total_all)

            n_atoms_cache: dict[str, int] = {}

            for n_threads in thread_counts:
                for pdb_id, n_atoms in structures:
                    json_path = json_dir.joinpath(f"{pdb_id}.json.gz")
                    if not json_path.exists():
                        json_path = json_dir.joinpath(f"{pdb_id}.json")
                    if not json_path.exists():
                        console.print(f"[yellow]Skip {pdb_id}: not found[/yellow]")
                        continue

                    # Resolve n_atoms lazily
                    if n_atoms == 0:
                        if pdb_id not in n_atoms_cache:
                            n_atoms_cache[pdb_id] = get_n_atoms_from_json(json_path)
                        n_atoms = n_atoms_cache[pdb_id]

                    # Warmup runs (not recorded)
                    for w in range(warmup):
                        desc = f"{pdb_id} t={n_threads} warmup {w + 1}/{warmup}"
                        progress.update(task, description=desc)
                        try:
                            run_benchmark(
                                tool_base,
                                json_path,
                                algorithm,
                                n_threads,
                                precision,
                                n_points,
                                use_bitmask,
                            )
                        except Exception:
                            pass
                        progress.advance(task)

                    # Measured runs
                    for run_num in range(1, runs + 1):
                        desc = f"{pdb_id} t={n_threads} run {run_num}/{runs}"
                        progress.update(task, description=desc)

                        try:
                            sasa_time, total_sasa = run_benchmark(
                                tool_base,
                                json_path,
                                algorithm,
                                n_threads,
                                precision,
                                n_points,
                                use_bitmask,
                            )

                            writer.writerow(
                                {
                                    "tool": tool_base,
                                    "structure": pdb_id,
                                    "n_atoms": n_atoms,
                                    "algorithm": algorithm,
                                    "precision": precision,
                                    "threads": n_threads,
                                    "run": run_num,
                                    "sasa_time_ms": sasa_time,
                                    "total_sasa": total_sasa,
                                }
                            )
                            f.flush()

                        except Exception as e:
                            console.print(
                                f"[red]Error: {pdb_id} t={n_threads}: {e}[/red]"
                            )

                        progress.advance(task)

    console.print(f"\n[green]Done![/green] Results saved to: {output_dir}")
    console.print(f"  - {csv_path.name}")
    console.print(f"  - {config_path.name}")

    # Print summary
    print_summary(csv_path, warmup, runs)


if __name__ == "__main__":
    app()

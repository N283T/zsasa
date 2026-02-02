#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Batch SASA benchmark using hyperfine (RustSASA paper methodology).

Flexible batch benchmark script that works with any PDB directory.
Uses hyperfine for timing measurements (warmup + multiple runs).

Usage:
    # E. coli proteome benchmark
    ./batch_bench.py -i benchmarks/UP000000625_83333_ECOLI_v6/pdb -n ecoli

    # SwissProt benchmark
    ./batch_bench.py -i benchmarks/swissprot_pdb_v6 -n swissprot --threads 10

    # Single tool test
    ./batch_bench.py -i /path/to/pdb -n test --tool zig --runs 1

Output:
    benchmarks/results/batch/<name>/
    ├── bench_zsasa_f64_8t.json
    ├── bench_zsasa_f32_8t.json
    ├── bench_freesasa_1t.json
    ├── bench_rustsasa_8t.json
    └── ...
"""

from __future__ import annotations

import json
import shutil
import subprocess
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Batch SASA benchmark (hyperfine-based)")
console = Console()


class Tool(str, Enum):
    zig = "zig"
    freesasa = "freesasa"
    rustsasa = "rustsasa"
    all = "all"


def get_root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def get_binary_paths() -> dict[str, Path]:
    """Get paths to tool binaries."""
    root = get_root_dir()
    return {
        "zsasa": root.joinpath("zig-out", "bin", "zsasa"),
        "freesasa_batch": root.joinpath(
            "benchmarks", "external", "freesasa-batch", "sasa_batch"
        ),
        "rustsasa": root.joinpath(
            "benchmarks", "external", "rustsasa-bench", "target", "release", "rust-sasa"
        ),
    }


def check_hyperfine() -> bool:
    """Check if hyperfine is available."""
    return shutil.which("hyperfine") is not None


def run_benchmark(
    name: str,
    cmd: str,
    results_dir: Path,
    warmup: int,
    runs: int,
    dry_run: bool,
) -> dict | None:
    """Run a single benchmark with hyperfine."""
    json_out = results_dir.joinpath(f"bench_{name}.json")

    console.print(f">>> [bold cyan]{name}[/]")

    if dry_run:
        console.print(f"    {cmd}")
        return None

    hyperfine_cmd = [
        "hyperfine",
        "--warmup",
        str(warmup),
        "--runs",
        str(runs),
        "--export-json",
        str(json_out),
        cmd,
    ]

    try:
        subprocess.run(hyperfine_cmd, check=True)
        console.print()

        # Parse results
        if json_out.exists():
            with open(json_out) as f:
                data = json.load(f)
                if data.get("results"):
                    return data["results"][0]
    except subprocess.CalledProcessError as e:
        console.print(f"[red]Error running benchmark: {e}[/]")

    return None


def run_zig(
    input_dir: Path,
    temp_out: Path,
    results_dir: Path,
    threads: int,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
) -> list[dict]:
    """Run zsasa benchmarks."""
    zsasa = binaries["zsasa"]
    if not zsasa.exists():
        console.print("[yellow][SKIP] zsasa not found[/]")
        return []

    out_dir = temp_out.joinpath("zig")
    if not dry_run:
        shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir(parents=True, exist_ok=True)

    results = []

    # f64 precision
    result = run_benchmark(
        f"zsasa_f64_{threads}t",
        f"{zsasa} {input_dir} {out_dir} --threads={threads} --parallelism=file --precision=f64",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": f"zsasa_f64_{threads}t", **result})

    result = run_benchmark(
        "zsasa_f64_1t",
        f"{zsasa} {input_dir} {out_dir} --threads=1 --precision=f64",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": "zsasa_f64_1t", **result})

    # f32 precision
    result = run_benchmark(
        f"zsasa_f32_{threads}t",
        f"{zsasa} {input_dir} {out_dir} --threads={threads} --parallelism=file --precision=f32",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": f"zsasa_f32_{threads}t", **result})

    result = run_benchmark(
        "zsasa_f32_1t",
        f"{zsasa} {input_dir} {out_dir} --threads=1 --precision=f32",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": "zsasa_f32_1t", **result})

    return results


def run_freesasa(
    input_dir: Path,
    temp_out: Path,
    results_dir: Path,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
) -> list[dict]:
    """Run FreeSASA benchmark (single-threaded only)."""
    sasa_batch = binaries["freesasa_batch"]
    if not sasa_batch.exists():
        console.print("[yellow][SKIP] sasa_batch not found[/]")
        return []

    out_dir = temp_out.joinpath("freesasa")
    if not dry_run:
        shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir(parents=True, exist_ok=True)

    results = []

    result = run_benchmark(
        "freesasa_1t",
        f"{sasa_batch} {input_dir} {out_dir}",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": "freesasa_1t", **result})

    return results


def run_rustsasa(
    input_dir: Path,
    temp_out: Path,
    results_dir: Path,
    threads: int,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
) -> list[dict]:
    """Run RustSASA benchmarks."""
    rustsasa = binaries["rustsasa"]
    if not rustsasa.exists():
        console.print("[yellow][SKIP] RustSASA not found[/]")
        return []

    out_dir = temp_out.joinpath("rustsasa")
    if not dry_run:
        shutil.rmtree(out_dir, ignore_errors=True)
        out_dir.mkdir(parents=True, exist_ok=True)

    results = []

    result = run_benchmark(
        f"rustsasa_{threads}t",
        f"{rustsasa} {input_dir} {out_dir} --format json -t {threads}",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": f"rustsasa_{threads}t", **result})

    result = run_benchmark(
        "rustsasa_1t",
        f"{rustsasa} {input_dir} {out_dir} --format json -t 1",
        results_dir,
        warmup,
        runs,
        dry_run,
    )
    if result:
        results.append({"name": "rustsasa_1t", **result})

    return results


def print_summary(results_dir: Path) -> None:
    """Print summary table from result JSON files."""
    table = Table(title="Benchmark Results")
    table.add_column("Tool", style="cyan")
    table.add_column("Mean (s)", justify="right")
    table.add_column("Std Dev", justify="right")
    table.add_column("Min (s)", justify="right")
    table.add_column("Max (s)", justify="right")

    for json_file in sorted(results_dir.glob("bench_*.json")):
        try:
            with open(json_file) as f:
                data = json.load(f)
                if data.get("results"):
                    r = data["results"][0]
                    name = json_file.stem.replace("bench_", "")
                    table.add_row(
                        name,
                        f"{r['mean']:.3f}",
                        f"±{r['stddev']:.3f}",
                        f"{r['min']:.3f}",
                        f"{r['max']:.3f}",
                    )
        except (json.JSONDecodeError, KeyError):
            continue

    console.print()
    console.print(table)


@app.command()
def main(
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
            help="Benchmark name (used for output directory)",
        ),
    ],
    threads: Annotated[
        int,
        typer.Option(
            "--threads",
            "-T",
            help="Number of threads for multi-threaded benchmarks",
        ),
    ] = 8,
    runs: Annotated[
        int,
        typer.Option(
            "--runs",
            "-r",
            help="Number of benchmark runs",
        ),
    ] = 3,
    warmup: Annotated[
        int,
        typer.Option(
            "--warmup",
            "-w",
            help="Number of warmup runs",
        ),
    ] = 3,
    tool: Annotated[
        Tool,
        typer.Option(
            "--tool",
            "-t",
            help="Tool to benchmark (zig, freesasa, rustsasa, or all)",
        ),
    ] = Tool.all,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output directory (default: benchmarks/results/batch/<name>)",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run",
            help="Show commands without running",
        ),
    ] = False,
) -> None:
    """Run batch SASA benchmarks using hyperfine."""
    # Check hyperfine
    if not check_hyperfine():
        console.print("[red]Error: hyperfine not found. Please install it first.[/]")
        raise typer.Exit(1)

    # Set up paths
    root = get_root_dir()
    if output_dir is None:
        results_dir = root.joinpath("benchmarks", "results", "batch", name)
    else:
        results_dir = output_dir

    temp_out = results_dir.joinpath("temp_out")
    binaries = get_binary_paths()

    # Create directories
    if not dry_run:
        results_dir.mkdir(parents=True, exist_ok=True)
        temp_out.mkdir(parents=True, exist_ok=True)

    # Print header
    console.print(f"[bold]=== Batch SASA Benchmark: {name} ===[/]")
    console.print(f"Input: {input_dir}")
    console.print(f"Output: {results_dir}")
    console.print(f"Warmup: {warmup}, Runs: {runs}, Threads: {threads}")
    console.print()

    all_results = []

    # Run benchmarks
    if tool in (Tool.zig, Tool.all):
        results = run_zig(
            input_dir, temp_out, results_dir, threads, warmup, runs, dry_run, binaries
        )
        all_results.extend(results)

    if tool in (Tool.freesasa, Tool.all):
        results = run_freesasa(
            input_dir, temp_out, results_dir, warmup, runs, dry_run, binaries
        )
        all_results.extend(results)

    if tool in (Tool.rustsasa, Tool.all):
        results = run_rustsasa(
            input_dir, temp_out, results_dir, threads, warmup, runs, dry_run, binaries
        )
        all_results.extend(results)

    # Print summary
    console.print(f"[bold green]=== Done! Results: {results_dir} ===[/]")

    if not dry_run and results_dir.exists():
        print_summary(results_dir)


if __name__ == "__main__":
    app()

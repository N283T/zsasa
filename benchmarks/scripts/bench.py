#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""SR single-file benchmark runner (hyperfine).

Runs wall-clock benchmarks for a single tool using the Shrake-Rupley algorithm
with configurable thread counts. Uses hyperfine for timing (includes I/O).

Tools: zig_f64, zig_f32, freesasa, rust, lahuta

Usage:
    ./benchmarks/scripts/bench.py --tool zig_f64 --threads 1,4,8
    ./benchmarks/scripts/bench.py --tool freesasa --threads 1 --warmup 3 --runs 10

    # Quick test (1 structure, 1 run, no warmup)
    ./benchmarks/scripts/bench.py --tool zig_f64 --threads 1 --warmup 0 --runs 1 \
        --input benchmarks/dataset/pdb/1gyt.pdb

Output:
    benchmarks/results/single/<n_points>/<tool>_sr/
    ├── config.json   # System info and parameters
    └── results.csv   # Benchmark results (mean, stddev, min, max in seconds)
"""

from __future__ import annotations

import csv
import json
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

from bench_common import (
    get_binary_path,
    get_n_atoms_from_pdb,
    get_system_info,
    load_sample_file,
    parse_threads,
    scan_input_directory,
)

app = typer.Typer(help="SR single-file benchmark runner (hyperfine)")
console = Console()

SR_TOOLS = ["zig_f64", "zig_f32", "freesasa", "rust", "lahuta"]
TOOL_ALIASES = {"zig": "zig_f64"}


def _parse_tool(tool: str) -> tuple[str, str, str]:
    """Parse tool name into (canonical, base, precision).

    Returns:
        (canonical_name, base_tool, precision)
    """
    tool = TOOL_ALIASES.get(tool, tool)
    if tool.startswith("zig_f"):
        return tool, "zig", tool.split("_")[1]
    return tool, tool, "f64"


def _build_command(
    base: str, precision: str, pdb_path: Path, n_threads: int, n_points: int
) -> str | None:
    """Build shell command for a tool. Returns None if tool should be skipped."""
    binary = get_binary_path(base if base != "lahuta" else "lahuta")

    if base == "zig":
        return (
            f"{binary} calc --algorithm=sr --threads={n_threads}"
            f" --precision={precision} --n-points={n_points}"
            f" {pdb_path} /dev/null"
        )
    elif base == "freesasa":
        return (
            f"{binary} --shrake-rupley --resolution={n_points}"
            f" --n-threads={n_threads} {pdb_path}"
        )
    elif base == "rust":
        return f"{binary} {pdb_path} -n {n_points} -t {n_threads} -o protein"
    elif base == "lahuta":
        return (
            f"{binary} sasa-sr -f {pdb_path} --is_af2_model"
            f" --points {n_points} -t {n_threads} -o /dev/null"
        )
    else:
        return None


def _run_hyperfine(
    cmd: str, warmup: int, runs: int, json_path: Path
) -> dict | None:
    """Run hyperfine and return results dict, or None on failure."""
    hyperfine_cmd = [
        "hyperfine",
        "--warmup",
        str(warmup),
        "--runs",
        str(runs),
        "--export-json",
        str(json_path),
        cmd,
    ]
    try:
        subprocess.run(
            hyperfine_cmd, check=True, capture_output=True, text=True, timeout=600
        )
        if json_path.exists():
            with open(json_path) as f:
                data = json.load(f)
            if data.get("results"):
                return data["results"][0]
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        pass
    return None


def _print_summary(csv_path: Path) -> None:
    """Print summary table from results CSV."""
    from collections import defaultdict

    by_threads: dict[int, list[float]] = defaultdict(list)

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            by_threads[int(row["threads"])].append(float(row["mean_s"]))

    if not by_threads:
        return

    table = Table(title="Summary (wall-clock per structure)")
    table.add_column("Threads", style="cyan", justify="right")
    table.add_column("N", justify="right")
    table.add_column("Mean (s)", justify="right")
    table.add_column("Min (s)", justify="right")
    table.add_column("Max (s)", justify="right")

    for t in sorted(by_threads):
        times = by_threads[t]
        n = len(times)
        mean = sum(times) / n
        table.add_row(
            str(t),
            f"{n:,}",
            f"{mean:.4f}",
            f"{min(times):.4f}",
            f"{max(times):.4f}",
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
            help="Tool: zig_f64, zig_f32, freesasa, rust, lahuta (zig = zig_f64)",
        ),
    ],
    threads: Annotated[
        str,
        typer.Option("--threads", "-T", help="Thread counts: '1,4,8' or '1-10'"),
    ] = "1,4,10",
    runs: Annotated[
        int,
        typer.Option("--runs", "-r", help="Number of hyperfine runs per configuration"),
    ] = 5,
    warmup: Annotated[
        int,
        typer.Option("--warmup", "-w", help="Hyperfine warmup runs"),
    ] = 1,
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Output directory"),
    ] = None,
    input_dir: Annotated[
        Path | None,
        typer.Option("--input-dir", "-i", help="Input directory with .pdb files"),
    ] = None,
    input_file: Annotated[
        Path | None,
        typer.Option("--input", help="Single .pdb file to benchmark"),
    ] = None,
    sample_file: Annotated[
        Path | None,
        typer.Option("--sample-file", "-S", help="Sample file to filter structures"),
    ] = None,
    n_points: Annotated[
        int,
        typer.Option("--n-points", "-N", help="Number of sphere test points per atom"),
    ] = 100,
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite existing results"),
    ] = False,
) -> None:
    """Run SR single-file benchmark using hyperfine."""
    # Check hyperfine
    if not shutil.which("hyperfine"):
        console.print("[red]Error: hyperfine not found. Install it first.[/red]")
        raise typer.Exit(1)

    # Parse tool
    all_valid = SR_TOOLS + list(TOOL_ALIASES.keys())
    if tool not in all_valid:
        console.print(f"[red]Error:[/red] Unknown tool: {tool}")
        console.print(f"Available: {', '.join(SR_TOOLS)} (zig = zig_f64)")
        raise typer.Exit(1)

    tool_canonical, tool_base, precision = _parse_tool(tool)
    thread_counts = parse_threads(threads)

    # Check binary exists
    binary = get_binary_path(tool_base)
    if not binary.exists():
        console.print(f"[red]Error:[/red] Binary not found: {binary}")
        raise typer.Exit(1)

    # Resolve input: single file or directory
    if input_file is not None:
        if not input_file.exists():
            console.print(f"[red]Error:[/red] Input file not found: {input_file}")
            raise typer.Exit(1)
        pdb_dir = input_file.parent
        pdb_id = input_file.stem
        structures = [(pdb_id, 0)]
    else:
        if input_dir is not None:
            if not input_dir.exists():
                console.print(
                    f"[red]Error:[/red] Input directory not found: {input_dir}"
                )
                raise typer.Exit(1)
            pdb_dir = input_dir
        else:
            pdb_dir = Path(__file__).parent.parent.joinpath("dataset", "pdb")
            if not pdb_dir.exists():
                console.print(f"[red]Error:[/red] Default dataset not found: {pdb_dir}")
                raise typer.Exit(1)

        structures = scan_input_directory(pdb_dir)
        if not structures:
            console.print(f"[red]Error:[/red] No .pdb files found in {pdb_dir}")
            raise typer.Exit(1)

        # Filter by sample file
        if sample_file is not None:
            if not sample_file.exists():
                console.print(
                    f"[red]Error:[/red] Sample file not found: {sample_file}"
                )
                raise typer.Exit(1)
            sample_ids = set(load_sample_file(sample_file))
            structures = [
                (pdb_id, n) for pdb_id, n in structures if pdb_id in sample_ids
            ]
            if not structures:
                console.print("[red]Error:[/red] No matching structures found")
                raise typer.Exit(1)
            console.print(
                f"Loaded [cyan]{len(sample_ids):,}[/cyan] samples from {sample_file}"
            )

    console.print(f"Found [cyan]{len(structures):,}[/cyan] structures in {pdb_dir}")

    # Setup output directory
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    if output_dir is None:
        output_dir = Path(__file__).parent.parent.joinpath(
            "results", "single", str(n_points), f"{tool_canonical}_sr"
        )

    existing_csv = output_dir.joinpath("results.csv")
    if existing_csv.exists() and not force:
        console.print(f"[yellow]Warning:[/yellow] Results already exist: {output_dir}")
        console.print("Use [bold]--force[/bold] to overwrite")
        raise typer.Exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Save config
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "tool": tool_canonical,
            "tool_base": tool_base,
            "algorithm": "sr",
            "precision": precision,
            "thread_counts": thread_counts,
            "warmup": warmup,
            "runs": runs,
            "n_points": n_points,
            "n_structures": len(structures),
            "input_dir": str(pdb_dir),
            "sample_file": str(sample_file) if sample_file else None,
        },
    }
    config_path = output_dir.joinpath("config.json")
    config_path.write_text(json.dumps(config, indent=2))

    total = len(structures) * len(thread_counts)

    console.print(f"\n[bold]{tool_canonical.upper()} SR (hyperfine)[/bold]")
    console.print(
        f"Threads: {thread_counts}, Warmup: {warmup}, Runs: {runs}, "
        f"Structures: {len(structures)}, Total benchmarks: {total}\n"
    )

    # Run benchmarks
    csv_path = output_dir.joinpath("results.csv")
    n_atoms_cache: dict[str, int] = {}

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
                "mean_s",
                "stddev_s",
                "min_s",
                "max_s",
                "median_s",
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
            task = progress.add_task("Running", total=total)

            with tempfile.TemporaryDirectory(prefix="bench_sr_") as tmpdir:
                for n_threads in thread_counts:
                    for pdb_id, n_atoms in structures:
                        pdb_path = pdb_dir.joinpath(f"{pdb_id}.pdb")
                        if not pdb_path.exists():
                            console.print(
                                f"[yellow]Skip {pdb_id}: not found[/yellow]"
                            )
                            progress.advance(task)
                            continue

                        # Resolve n_atoms lazily
                        if n_atoms == 0:
                            if pdb_id not in n_atoms_cache:
                                n_atoms_cache[pdb_id] = get_n_atoms_from_pdb(pdb_path)
                            n_atoms = n_atoms_cache[pdb_id]

                        desc = f"{pdb_id} t={n_threads}"
                        progress.update(task, description=desc)

                        cmd = _build_command(
                            tool_base, precision, pdb_path, n_threads, n_points
                        )
                        if cmd is None:
                            progress.advance(task)
                            continue

                        json_path = Path(tmpdir).joinpath(
                            f"{pdb_id}_{n_threads}t.json"
                        )
                        result = _run_hyperfine(cmd, warmup, runs, json_path)

                        if result:
                            writer.writerow(
                                {
                                    "tool": tool_canonical,
                                    "structure": pdb_id,
                                    "n_atoms": n_atoms,
                                    "algorithm": "sr",
                                    "precision": precision,
                                    "threads": n_threads,
                                    "mean_s": result["mean"],
                                    "stddev_s": result["stddev"],
                                    "min_s": result["min"],
                                    "max_s": result["max"],
                                    "median_s": result["median"],
                                }
                            )
                            f.flush()
                        else:
                            console.print(
                                f"[red]Error: {pdb_id} t={n_threads}: "
                                f"hyperfine failed[/red]"
                            )

                        progress.advance(task)

    console.print(f"\n[green]Done![/green] Results saved to: {output_dir}")
    console.print(f"  - {csv_path.name}")
    console.print(f"  - {config_path.name}")

    _print_summary(csv_path)


if __name__ == "__main__":
    app()

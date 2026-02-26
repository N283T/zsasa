#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""LR benchmark runner.

Runs benchmark for a single tool using the Lee-Richards algorithm
with configurable thread counts.
Uses internal timing (SASA-only) with warmup runs for statistical reliability.
Results are saved to CSV with execution config in JSON.

Tools: zig_f64 (default), zig_f32, freesasa (zig = zig_f64 alias)
Note: RustSASA does not support LR algorithm.

Usage:
    ./benchmarks/scripts/bench_lr.py --tool zig_f64
    ./benchmarks/scripts/bench_lr.py --tool zig_f32
    ./benchmarks/scripts/bench_lr.py --tool freesasa

    # "zig" is shorthand for "zig_f64"
    ./benchmarks/scripts/bench_lr.py --tool zig

    # Quick test (no warmup, 1 run)
    ./benchmarks/scripts/bench_lr.py --tool zig --threads 1 --warmup 0 --runs 1

    # Replay from existing config
    ./benchmarks/scripts/bench_lr.py --config results/single_lr/20/zig_f64_lr/config.json

Output:
    benchmarks/results/single_lr/<n_slices>/<tool>_lr/
    ├── config.json   # System info and parameters
    └── results.csv   # Benchmark results (warmup runs excluded)

    Examples: single_lr/20/zig_f64_lr/, single_lr/20/freesasa_lr/
"""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn

from bench_common import (
    TOOL_ALIASES,
    TOOLS,
    get_n_atoms_from_json,
    get_system_info,
    load_sample_file,
    parse_threads,
    parse_tool,
    print_summary,
    run_benchmark,
    scan_input_directory,
)

app = typer.Typer(help="LR benchmark runner")
console = Console()

LR_TOOLS = ["zig_f64", "zig_f32", "freesasa"]


def _resolve_params(
    cfg: dict,
    *,
    tool: str | None,
    threads: str | None,
    runs: int | None,
    warmup: int | None,
    n_slices: int | None,
    input_dir: Path | None,
    sample_file: Path | None,
) -> dict:
    """Resolve parameters: CLI > config > defaults."""
    return {
        "tool": tool or cfg.get("tool"),
        "threads": threads
        or (
            ",".join(str(t) for t in cfg["thread_counts"])
            if "thread_counts" in cfg
            else "1,4,10"
        ),
        "runs": runs if runs is not None else cfg.get("runs", 5),
        "warmup": warmup if warmup is not None else cfg.get("warmup", 1),
        "n_slices": n_slices if n_slices is not None else cfg.get("n_slices", 20),
        "input_dir": input_dir
        or (Path(cfg["input_dir"]) if cfg.get("input_dir") else None),
        "sample_file": (
            sample_file
            or (Path(cfg["sample_file"]) if cfg.get("sample_file") else None)
        ),
    }


@app.command()
def main(
    config_file: Annotated[
        Path | None,
        typer.Option(
            "--config",
            "-c",
            help="Load parameters from existing config.json (CLI args override)",
        ),
    ] = None,
    tool: Annotated[
        str | None,
        typer.Option(
            "--tool",
            "-t",
            help="Tool: zig_f64, zig_f32, freesasa (zig = zig_f64)",
        ),
    ] = None,
    threads: Annotated[
        str | None,
        typer.Option("--threads", "-T", help="Thread counts: '1,4,10' or '1-10'"),
    ] = None,
    runs: Annotated[
        int | None,
        typer.Option("--runs", "-r", help="Number of measured runs per configuration"),
    ] = None,
    warmup: Annotated[
        int | None,
        typer.Option("--warmup", "-w", help="Warmup runs (not recorded)"),
    ] = None,
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
    n_slices: Annotated[
        int | None,
        typer.Option(
            "--n-slices",
            help="Number of slices per atom diameter (default: 20)",
        ),
    ] = None,
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite existing results"),
    ] = False,
) -> None:
    """Run LR benchmark for a single tool.

    Use --config to replay from an existing config.json. CLI args override config values.
    """

    # Load config if provided
    cfg: dict = {}
    if config_file is not None:
        if not config_file.exists():
            console.print(f"[red]Error:[/red] Config file not found: {config_file}")
            raise typer.Exit(1)
        cfg = json.loads(config_file.read_text()).get("parameters", {})
        console.print(f"[cyan]Loaded config:[/cyan] {config_file}")

    # Resolve: CLI > config > defaults
    p = _resolve_params(
        cfg,
        tool=tool,
        threads=threads,
        runs=runs,
        warmup=warmup,
        n_slices=n_slices,
        input_dir=input_dir,
        sample_file=sample_file,
    )
    tool = p["tool"]
    algorithm = "lr"
    runs = p["runs"]
    warmup = p["warmup"]
    n_slices = p["n_slices"]
    input_dir = p["input_dir"]
    sample_file = p["sample_file"]

    # Validate required params
    if not tool:
        console.print("[red]Error:[/red] --tool is required (or use --config)")
        raise typer.Exit(1)

    thread_counts = parse_threads(p["threads"])

    # Parse and validate tool
    all_valid = list(LR_TOOLS) + [k for k, v in TOOL_ALIASES.items() if v in LR_TOOLS]
    if tool not in all_valid:
        console.print(f"[red]Error:[/red] Unknown or unsupported tool for LR: {tool}")
        console.print(f"Available: {', '.join(LR_TOOLS)} (zig = zig_f64)")
        raise typer.Exit(1)

    tool_canonical, tool_base, precision = parse_tool(tool)

    if tool_base == "rust":
        console.print("[red]Error:[/red] RustSASA does not support LR algorithm")
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
            "results",
            "single_lr",
            str(n_slices),
            f"{tool_canonical}_lr",
        )

    # Overwrite protection
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
            "algorithm": algorithm,
            "precision": precision,
            "thread_counts": thread_counts,
            "warmup": warmup,
            "runs": runs,
            "n_slices": n_slices,
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

    console.print(f"\n[bold]{tool_canonical.upper()} LR[/bold]")
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
                                n_slices=n_slices,
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
                                n_slices=n_slices,
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

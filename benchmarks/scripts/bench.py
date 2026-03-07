#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""SR single-file benchmark runner (hyperfine).

Runs wall-clock benchmarks using the Shrake-Rupley algorithm with configurable
thread counts. Uses hyperfine for timing (includes I/O).

Tools: zig_f64, zig_f32, zig_f64_bitmask, zig_f32_bitmask, freesasa, rust
       (zig = zig_f64, zig_bitmask = zig_f64_bitmask)

Usage:
    # All tools (default)
    ./benchmarks/scripts/bench.py --threads 1,4,8

    # Specific tools
    ./benchmarks/scripts/bench.py --tool zig_f64 --tool freesasa --threads 1,4,8

    # Single tool
    ./benchmarks/scripts/bench.py --tool zig_f64 --threads 1,4,8

    # Quick test (single file, 1 hyperfine run, no warmup)
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
import re
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.progress import BarColumn, Progress, SpinnerColumn, TextColumn

from bench_common import (
    LAHUTA_BITMASK_POINTS,
    TOOL_ALIASES,
    check_hyperfine,
    ensure_zsasa_built,
    get_binary_path,
    get_n_atoms_from_pdb,
    get_system_info,
    load_sample_file,
    parse_threads,
    parse_tool,
    print_hyperfine_summary,
    quote_path,
    run_hyperfine,
    scan_input_directory,
)

app = typer.Typer(help="SR single-file benchmark runner (hyperfine)")
console = Console()

SR_TOOLS = [
    "zig_f64",
    "zig_f32",
    "zig_f64_bitmask",
    "zig_f32_bitmask",
    "freesasa",
    "rustsasa",
]


_TIMING_TOOLS = {"zig", "freesasa", "rustsasa"}


def _build_command(
    base: str,
    precision: str,
    pdb_path: Path,
    n_threads: int,
    n_points: int,
    use_bitmask: bool = False,
    *,
    timing: bool = False,
) -> str:
    """Build shell command for a tool.

    When timing=True, appends --timing flag (supported by zig, freesasa, rustsasa).

    Raises:
        ValueError: If tool base is not recognized, or if use_bitmask is True
            for a tool that does not support bitmask mode.
    """
    if use_bitmask and base not in ("zig", "lahuta"):
        raise ValueError(
            f"Bitmask mode is not supported for tool base '{base}'. "
            f"Only 'zig' and 'lahuta' support --use-bitmask."
        )

    binary = quote_path(get_binary_path(base))
    quoted = quote_path(pdb_path)
    bitmask_flag = " --use-bitmask" if use_bitmask else ""
    timing_flag = " --timing" if timing else ""

    if base == "zig":
        return (
            f"{binary} calc --algorithm=sr --threads={n_threads}"
            f" --precision={precision} --n-points={n_points}"
            f"{bitmask_flag}{timing_flag}"
            f" {quoted} /dev/null"
        )
    elif base == "freesasa":
        return (
            f"{binary} --shrake-rupley --resolution={n_points}"
            f" --n-threads={n_threads}{timing_flag} {quoted}"
        )
    elif base == "rustsasa":
        return (
            f"{binary} {quoted} /dev/null -n {n_points} -t {n_threads}"
            f" -o protein --allow-vdw-fallback{timing_flag}"
        )
    elif base == "lahuta":
        return (
            f"{binary} sasa-sr -f {quoted} --is_af2_model"
            f" --points {n_points}{bitmask_flag} -t {n_threads}"
            f" --stdout --progress 0 > /dev/null"
        )
    else:
        raise ValueError(f"No command builder for tool base: {base}")


_TIMING_RE = re.compile(r"^(\w+_TIME_MS):([0-9.]+)\s*$", re.MULTILINE)


def _run_timing(cmd: str, timeout: int = 60) -> dict[str, float] | None:
    """Run a command and parse PARSE_TIME_MS/SASA_TIME_MS/TOTAL_TIME_MS from stderr.

    Keys are lowercased: e.g. SASA_TIME_MS -> sasa_time_ms.
    """
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode != 0:
            snippet = (result.stderr or "").strip()[-200:]
            console.print(
                f"[yellow]Warning:[/yellow] timing command failed "
                f"(exit {result.returncode}): {cmd[:80]}"
            )
            if snippet:
                console.print(f"[dim]  {snippet}[/dim]")
            return None
        timing = {}
        for match in _TIMING_RE.finditer(result.stderr):
            timing[match.group(1).lower()] = float(match.group(2))
        if not timing:
            console.print(
                f"[yellow]Warning:[/yellow] no timing markers in stderr: {cmd[:80]}"
            )
        return timing if timing else None
    except subprocess.TimeoutExpired:
        console.print(
            f"[yellow]Warning:[/yellow] timing command timed out ({timeout}s)"
        )
        return None
    except OSError as e:
        console.print(f"[yellow]Warning:[/yellow] timing command error: {e}")
        return None


def _resolve_tools(tools: list[str] | None) -> list[str]:
    """Resolve tool list: None or ["all"] -> SR_TOOLS, otherwise validate."""
    if tools is None or tools == ["all"]:
        return list(SR_TOOLS)

    all_valid = SR_TOOLS + list(TOOL_ALIASES.keys())
    for t in tools:
        if t == "all":
            return list(SR_TOOLS)
        if t not in all_valid:
            console.print(f"[red]Error:[/red] Unknown tool: {t}")
            console.print(
                f"Available: all, {', '.join(SR_TOOLS)} "
                f"(zig = zig_f64, zig_bitmask = zig_f64_bitmask)"
            )
            raise typer.Exit(1)
    return list(tools)


def _run_tool(
    tool_name: str,
    *,
    structures: list[tuple[str, int]],
    pdb_dir: Path,
    thread_counts: list[int],
    n_points: int,
    warmup: int,
    runs: int,
    timeout: int,
    prepare: str | None,
    output_dir: Path,
    force: bool,
    timestamp: str,
    sample_file: Path | None,
) -> bool:
    """Run benchmarks for a single tool. Returns True on success."""
    tool_canonical, tool_base, precision, use_bitmask = parse_tool(tool_name)

    # Validate lahuta bitmask n_points
    if tool_base == "lahuta" and use_bitmask and n_points not in LAHUTA_BITMASK_POINTS:
        console.print(
            f"[red]Error:[/red] lahuta_bitmask only supports n_points "
            f"{sorted(LAHUTA_BITMASK_POINTS)}, got {n_points}."
        )
        return False

    # Check binary exists
    binary = get_binary_path(tool_base)
    if not binary.exists():
        console.print(f"[red]Error:[/red] Binary not found: {binary}")
        return False

    # Setup output directory
    existing_csv = output_dir.joinpath("results.csv")
    if existing_csv.exists() and not force:
        console.print(f"[yellow]Warning:[/yellow] Results already exist: {output_dir}")
        console.print("Use [bold]--force[/bold] to overwrite")
        return False

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
            "prepare": prepare,
        },
    }
    config_path = output_dir.joinpath("config.json")
    config_path.write_text(json.dumps(config, indent=2))

    total = len(structures) * len(thread_counts)
    prepare_info = f", Prepare: '{prepare}'" if prepare else ""
    console.print(f"\n[bold]{tool_canonical.upper()} SR (hyperfine)[/bold]")
    console.print(
        f"Threads: {thread_counts}, Warmup: {warmup}, Runs: {runs}, "
        f"Structures: {len(structures)}, Total benchmarks: {total}"
        f"{prepare_info}\n"
    )

    # Run benchmarks
    csv_path = output_dir.joinpath("results.csv")
    n_atoms_cache: dict[str, int | None] = {}
    n_success = 0
    n_failed = 0

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
                "user_s",
                "system_s",
                "memory_bytes",
                "parse_time_ms",
                "sasa_time_ms",
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
                            console.print(f"[yellow]Skip {pdb_id}: not found[/yellow]")
                            n_failed += 1
                            progress.advance(task)
                            continue

                        # lahuta only supports AlphaFold models
                        if tool_base == "lahuta" and not pdb_id.startswith("af-"):
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
                            tool_base,
                            precision,
                            pdb_path,
                            n_threads,
                            n_points,
                            use_bitmask,
                        )

                        # lahuta writes report files to cwd; redirect to tmpdir
                        if tool_base == "lahuta":
                            cmd = f"cd {quote_path(tmpdir)} && {cmd}"

                        json_path = Path(tmpdir).joinpath(f"{pdb_id}_{n_threads}t.json")
                        result = run_hyperfine(
                            cmd,
                            warmup,
                            runs,
                            json_path,
                            timeout=timeout,
                            prepare=prepare,
                        )

                        if result:
                            # Run internal timing (single additional run)
                            timing = None
                            if tool_base in _TIMING_TOOLS:
                                timing_cmd = _build_command(
                                    tool_base,
                                    precision,
                                    pdb_path,
                                    n_threads,
                                    n_points,
                                    use_bitmask,
                                    timing=True,
                                )
                                timing = _run_timing(timing_cmd, timeout=timeout)

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
                                    "user_s": result.get("user"),
                                    "system_s": result.get("system"),
                                    "memory_bytes": (
                                        result["memory_usage_byte"][0]
                                        if result.get("memory_usage_byte")
                                        else None
                                    ),
                                    "parse_time_ms": (
                                        timing.get("parse_time_ms") if timing else None
                                    ),
                                    "sasa_time_ms": (
                                        timing.get("sasa_time_ms") if timing else None
                                    ),
                                }
                            )
                            f.flush()
                            n_success += 1
                        else:
                            n_failed += 1

                        progress.advance(task)

    if n_failed > 0:
        console.print(
            f"\n[yellow]Warning:[/yellow] {n_failed}/{total} benchmarks failed"
        )

    if n_success == 0:
        console.print("[red]Error: all benchmarks failed, no results recorded[/red]")
        return False

    console.print(f"\n[green]Done![/green] {n_success}/{total} benchmarks completed")
    console.print(f"  Results: {output_dir}")
    console.print(f"  - {csv_path.name}")
    console.print(f"  - {config_path.name}")

    print_hyperfine_summary(csv_path)
    return True


@app.command()
def main(
    tools: Annotated[
        list[str] | None,
        typer.Option(
            "--tool",
            "-t",
            help=(
                "Tools to benchmark (can specify multiple: --tool X --tool Y). "
                "Default: all. "
                "Available: all, zig_f64, zig_f32, zig_f64_bitmask, zig_f32_bitmask, "
                "freesasa, rust "
                "(zig = zig_f64, zig_bitmask = zig_f64_bitmask)"
            ),
        ),
    ] = None,
    threads: Annotated[
        str,
        typer.Option("--threads", "-T", help="Thread counts: '1,4,10' or '1-10'"),
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
        typer.Option(
            "--output-dir",
            "-o",
            help="Output directory (per-tool subdirs when multiple tools)",
        ),
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
    timeout: Annotated[
        int,
        typer.Option(
            "--timeout", help="Timeout per benchmark in seconds (default: 600)"
        ),
    ] = 600,
    prepare: Annotated[
        str | None,
        typer.Option(
            "--prepare",
            "-p",
            help="Shell command to run before each timing run (passed to hyperfine --prepare). "
            "E.g. 'sync' or 'sudo purge' (macOS) to clear filesystem caches.",
        ),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Show commands without running"),
    ] = False,
    force: Annotated[
        bool,
        typer.Option("--force", "-f", help="Overwrite existing results"),
    ] = False,
) -> None:
    """Run SR single-file benchmark using hyperfine."""
    # Mutual exclusivity check
    if input_file is not None and input_dir is not None:
        console.print(
            "[red]Error:[/red] --input and --input-dir are mutually exclusive"
        )
        raise typer.Exit(1)

    # Resolve and validate tools
    selected_tools = _resolve_tools(tools)

    # Check hyperfine
    if not dry_run and not check_hyperfine():
        console.print("[red]Error: hyperfine not found. Install it first.[/red]")
        raise typer.Exit(1)

    thread_counts = parse_threads(threads)

    # Auto-build zsasa if any zig tool selected
    has_zig = any(parse_tool(t)[1] == "zig" for t in selected_tools)
    if not dry_run and has_zig:
        ensure_zsasa_built()

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
                console.print(f"[red]Error:[/red] Sample file not found: {sample_file}")
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
    console.print(f"Tools: [cyan]{', '.join(selected_tools)}[/cyan]")

    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")

    # Dry run: show commands for all tools and exit
    if dry_run:
        for tool_name in selected_tools:
            tool_canonical, tool_base, precision, use_bitmask = parse_tool(tool_name)
            console.print(f"\n[bold]{tool_canonical.upper()} SR[/bold]")
            for n_threads in thread_counts:
                for pdb_id, _n_atoms in structures:
                    pdb_path = pdb_dir.joinpath(f"{pdb_id}.pdb")
                    cmd = _build_command(
                        tool_base, precision, pdb_path, n_threads, n_points, use_bitmask
                    )
                    console.print(f"  [dim]{cmd}[/dim]")
        console.print("\n[bold green]Dry run complete.[/bold green]")
        return

    # Run benchmarks for each tool
    results_base = Path(__file__).parent.parent.joinpath(
        "results", "single", str(n_points)
    )
    n_succeeded = 0
    n_tools = len(selected_tools)

    for tool_name in selected_tools:
        tool_canonical = parse_tool(tool_name)[0]
        tool_output_dir = (
            output_dir
            if output_dir is not None and n_tools == 1
            else results_base.joinpath(f"{tool_canonical}_sr")
        )

        success = _run_tool(
            tool_name,
            structures=structures,
            pdb_dir=pdb_dir,
            thread_counts=thread_counts,
            n_points=n_points,
            warmup=warmup,
            runs=runs,
            timeout=timeout,
            prepare=prepare,
            output_dir=tool_output_dir,
            force=force,
            timestamp=timestamp,
            sample_file=sample_file,
        )
        if success:
            n_succeeded += 1

    if n_tools > 1:
        console.print(
            f"\n[bold]{'=' * 40}[/bold]"
            f"\n[bold green]Completed {n_succeeded}/{n_tools} tools[/bold green]"
        )


if __name__ == "__main__":
    app()

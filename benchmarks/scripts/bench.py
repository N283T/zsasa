#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""SR single-file benchmark runner.

Two independent subcommands:
  - wall: Wall-clock timing via hyperfine. Saves per-run JSON + aggregated CSV.
  - sasa: Internal tool timing via --timing flag (parse_time_ms, sasa_time_ms).
          Updates existing CSV from a previous `wall` run.

Tools: zig_f64, zig_f32, zig_f64_bitmask, zig_f32_bitmask, freesasa, rust
       (zig = zig_f64, zig_bitmask = zig_f64_bitmask)

Usage:
    # Wall-clock benchmarks (hyperfine)
    ./benchmarks/scripts/bench.py wall --threads 1,4,8

    # Add internal timing to existing results
    ./benchmarks/scripts/bench.py sasa

    # Specific tools
    ./benchmarks/scripts/bench.py wall --tool zig_f64 --tool freesasa --threads 1,4,8

    # Quick test (single file, 1 hyperfine run, no warmup)
    ./benchmarks/scripts/bench.py wall --tool zig_f64 --threads 1 --warmup 0 --runs 1 \
        --input benchmarks/dataset/pdb/1gyt.pdb

Output:
    benchmarks/results/single/<n_points>/<tool>_sr/
    ├── config.json   # System info and parameters
    ├── runs/         # Per-run hyperfine JSON (individual times preserved)
    │   ├── <structure>_<threads>t.json
    │   └── ...
    └── results.csv   # Aggregated results (generated from runs/)
"""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import statistics
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

app = typer.Typer(help="SR single-file benchmark runner")
console = Console()


CSV_COLUMNS = [
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
]

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


def _aggregate_runs_to_csv(
    output_dir: Path,
    tool_canonical: str,
    precision: str,
) -> None:
    """Generate results.csv from individual hyperfine JSON files in runs/.

    Reads each JSON, extracts per-run times and metadata, computes aggregates.
    """
    runs_dir = output_dir.joinpath("runs")
    if not runs_dir.exists():
        return

    json_files = sorted(runs_dir.glob("*.json"))
    if not json_files:
        return

    csv_path = output_dir.joinpath("results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()

        for json_file in json_files:
            try:
                data = json.loads(json_file.read_text())
            except (json.JSONDecodeError, OSError):
                continue

            result = data.get("result")
            meta = data.get("meta", {})
            if not result:
                continue

            times = result.get("times", [])
            if not times:
                continue

            writer.writerow(
                {
                    "tool": tool_canonical,
                    "structure": meta.get("structure", ""),
                    "n_atoms": meta.get("n_atoms", 0),
                    "algorithm": "sr",
                    "precision": precision,
                    "threads": meta.get("threads", 1),
                    "mean_s": statistics.mean(times),
                    "stddev_s": statistics.stdev(times) if len(times) > 1 else 0,
                    "min_s": min(times),
                    "max_s": max(times),
                    "median_s": statistics.median(times),
                    "user_s": result.get("user"),
                    "system_s": result.get("system"),
                    "memory_bytes": (
                        result["memory_usage_byte"][0]
                        if result.get("memory_usage_byte")
                        else None
                    ),
                    "parse_time_ms": "",
                    "sasa_time_ms": "",
                }
            )

    console.print(f"[green]Aggregated:[/green] {csv_path}")


def _resolve_input(
    input_file: Path | None,
    input_dir: Path | None,
    sample_file: Path | None,
) -> tuple[Path, list[tuple[str, int]]]:
    """Resolve input files and return (pdb_dir, structures)."""
    if input_file is not None and input_dir is not None:
        console.print(
            "[red]Error:[/red] --input and --input-dir are mutually exclusive"
        )
        raise typer.Exit(1)

    if input_file is not None:
        if not input_file.exists():
            console.print(f"[red]Error:[/red] Input file not found: {input_file}")
            raise typer.Exit(1)
        return input_file.parent, [(input_file.stem, 0)]

    if input_dir is not None:
        if not input_dir.exists():
            console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
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

    if sample_file is not None:
        if not sample_file.exists():
            console.print(f"[red]Error:[/red] Sample file not found: {sample_file}")
            raise typer.Exit(1)
        sample_ids = set(load_sample_file(sample_file))
        structures = [(pdb_id, n) for pdb_id, n in structures if pdb_id in sample_ids]
        if not structures:
            console.print("[red]Error:[/red] No matching structures found")
            raise typer.Exit(1)
        console.print(
            f"Loaded [cyan]{len(sample_ids):,}[/cyan] samples from {sample_file}"
        )

    return pdb_dir, structures


# === Shared CLI Option Helpers ===

_TOOL_HELP = (
    "Tools to benchmark (can specify multiple: --tool X --tool Y). "
    "Default: all. "
    "Available: all, zig_f64, zig_f32, zig_f64_bitmask, zig_f32_bitmask, "
    "freesasa, rust "
    "(zig = zig_f64, zig_bitmask = zig_f64_bitmask)"
)


# === wall command ===


@app.command()
def wall(
    tools: Annotated[
        list[str] | None, typer.Option("--tool", "-t", help=_TOOL_HELP)
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
        int, typer.Option("--warmup", "-w", help="Hyperfine warmup runs")
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
        int, typer.Option("--timeout", help="Timeout per benchmark in seconds")
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
        bool, typer.Option("--dry-run", help="Show commands without running")
    ] = False,
    force: Annotated[
        bool, typer.Option("--force", "-f", help="Overwrite existing results")
    ] = False,
) -> None:
    """Run wall-clock benchmarks via hyperfine. Saves per-run JSON + aggregated CSV."""
    selected_tools = _resolve_tools(tools)

    if not dry_run and not check_hyperfine():
        console.print("[red]Error: hyperfine not found. Install it first.[/red]")
        raise typer.Exit(1)

    thread_counts = parse_threads(threads)
    pdb_dir, structures = _resolve_input(input_file, input_dir, sample_file)

    has_zig = any(parse_tool(t)[1] == "zig" for t in selected_tools)
    if not dry_run and has_zig:
        ensure_zsasa_built()

    console.print(f"Found [cyan]{len(structures):,}[/cyan] structures in {pdb_dir}")
    console.print(f"Tools: [cyan]{', '.join(selected_tools)}[/cyan]")

    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")

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

    results_base = Path(__file__).parent.parent.joinpath(
        "results", "single", str(n_points)
    )
    n_succeeded = 0
    n_tools = len(selected_tools)

    for tool_name in selected_tools:
        tool_canonical, tool_base, precision, use_bitmask = parse_tool(tool_name)
        tool_output_dir = (
            output_dir
            if output_dir is not None and n_tools == 1
            else results_base.joinpath(f"{tool_canonical}_sr")
        )

        success = _run_wall(
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


def _run_wall(
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
    """Run hyperfine benchmarks for a single tool. Returns True on success.

    Saves individual hyperfine JSON to runs/ directory, then aggregates to CSV.
    """
    tool_canonical, tool_base, precision, use_bitmask = parse_tool(tool_name)

    if tool_base == "lahuta" and use_bitmask and n_points not in LAHUTA_BITMASK_POINTS:
        console.print(
            f"[red]Error:[/red] lahuta_bitmask only supports n_points "
            f"{sorted(LAHUTA_BITMASK_POINTS)}, got {n_points}."
        )
        return False

    binary = get_binary_path(tool_base)
    if not binary.exists():
        console.print(f"[red]Error:[/red] Binary not found: {binary}")
        return False

    runs_dir = output_dir.joinpath("runs")
    if runs_dir.exists() and not force:
        console.print(f"[yellow]Warning:[/yellow] Results already exist: {output_dir}")
        console.print("Use [bold]--force[/bold] to overwrite")
        return False

    output_dir.mkdir(parents=True, exist_ok=True)
    if runs_dir.exists() and force:
        shutil.rmtree(runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)

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
    console.print(f"\n[bold]{tool_canonical.upper()} SR (wall)[/bold]")
    console.print(
        f"Threads: {thread_counts}, Warmup: {warmup}, Runs: {runs}, "
        f"Structures: {len(structures)}, Total benchmarks: {total}"
        f"{prepare_info}\n"
    )

    n_atoms_cache: dict[str, int | None] = {}
    n_success = 0
    n_failed = 0

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

                    if tool_base == "lahuta" and not pdb_id.startswith("af-"):
                        progress.advance(task)
                        continue

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

                    if tool_base == "lahuta":
                        cmd = f"cd {quote_path(tmpdir)} && {cmd}"

                    # Run hyperfine, save JSON to temp
                    hf_json_path = Path(tmpdir).joinpath(f"{pdb_id}_{n_threads}t.json")
                    result = run_hyperfine(
                        cmd,
                        warmup,
                        runs,
                        hf_json_path,
                        timeout=timeout,
                        prepare=prepare,
                    )

                    if result:
                        # Save enriched JSON to runs/ with metadata
                        run_data = {
                            "meta": {
                                "structure": pdb_id,
                                "n_atoms": n_atoms,
                                "threads": n_threads,
                                "tool": tool_canonical,
                                "precision": precision,
                                "n_points": n_points,
                            },
                            "result": result,
                        }
                        dest = runs_dir.joinpath(f"{pdb_id}_{n_threads}t.json")
                        dest.write_text(json.dumps(run_data, indent=2))
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

    # Aggregate runs/ to results.csv
    _aggregate_runs_to_csv(output_dir, tool_canonical, precision)

    console.print(f"\n[green]Done![/green] {n_success}/{total} benchmarks completed")
    console.print(f"  Results: {output_dir}")
    console.print(f"  - runs/ ({n_success} JSON files)")
    console.print(f"  - results.csv")
    console.print(f"  - config.json")

    print_hyperfine_summary(output_dir.joinpath("results.csv"))
    return True


# === sasa command ===


@app.command()
def sasa(
    tools: Annotated[
        list[str] | None, typer.Option("--tool", "-t", help=_TOOL_HELP)
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
        int, typer.Option("--timeout", help="Timeout per benchmark in seconds")
    ] = 600,
) -> None:
    """Run internal SASA timing. Updates parse_time_ms/sasa_time_ms in existing CSV."""
    selected_tools = _resolve_tools(tools)
    pdb_dir, _ = _resolve_input(input_file, input_dir, sample_file)

    has_zig = any(parse_tool(t)[1] == "zig" for t in selected_tools)
    if has_zig:
        ensure_zsasa_built()

    results_base = Path(__file__).parent.parent.joinpath(
        "results", "single", str(n_points)
    )
    n_succeeded = 0
    n_tools = len(selected_tools)

    for tool_name in selected_tools:
        tool_canonical = parse_tool(tool_name)[0]
        tool_output_dir = results_base.joinpath(f"{tool_canonical}_sr")

        success = _run_sasa(
            tool_name,
            pdb_dir=pdb_dir,
            n_points=n_points,
            timeout=timeout,
            output_dir=tool_output_dir,
        )
        if success:
            n_succeeded += 1

    if n_tools > 1:
        console.print(
            f"\n[bold]{'=' * 40}[/bold]"
            f"\n[bold green]Completed {n_succeeded}/{n_tools} tools[/bold green]"
        )


def _run_sasa(
    tool_name: str,
    *,
    pdb_dir: Path,
    n_points: int,
    timeout: int,
    output_dir: Path,
) -> bool:
    """Run internal timing for a single tool, updating existing CSV."""
    tool_canonical, tool_base, precision, use_bitmask = parse_tool(tool_name)

    if tool_base not in _TIMING_TOOLS:
        console.print(
            f"[dim]Skipping timing for {tool_canonical} (not supported)[/dim]"
        )
        return True

    csv_path = output_dir.joinpath("results.csv")
    if not csv_path.exists():
        console.print(
            f"[red]Error:[/red] No results.csv found in {output_dir}. "
            f"Run 'wall' command first."
        )
        return False

    binary = get_binary_path(tool_base)
    if not binary.exists():
        console.print(f"[red]Error:[/red] Binary not found: {binary}")
        return False

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        console.print(f"[yellow]Warning:[/yellow] Empty CSV: {csv_path}")
        return True

    tool_rows = [i for i, r in enumerate(rows) if r["tool"] == tool_canonical]
    if not tool_rows:
        console.print(
            f"[yellow]Warning:[/yellow] No rows for {tool_canonical} in {csv_path}"
        )
        return True

    console.print(f"\n[bold]{tool_canonical.upper()} SR (sasa)[/bold]")
    console.print(f"Structures × threads: {len(tool_rows)}, Timeout: {timeout}s\n")

    n_success = 0
    n_failed = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        console=console,
    ) as progress:
        task = progress.add_task("Timing", total=len(tool_rows))

        for idx in tool_rows:
            row = rows[idx]
            structure = row["structure"]
            n_threads = int(row["threads"])
            pdb_path = pdb_dir.joinpath(f"{structure}.pdb")

            desc = f"{structure} t={n_threads}"
            progress.update(task, description=desc)

            if not pdb_path.exists():
                console.print(f"[yellow]Skip {structure}: not found[/yellow]")
                n_failed += 1
                progress.advance(task)
                continue

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

            if timing:
                rows[idx]["parse_time_ms"] = timing.get("parse_time_ms", "")
                rows[idx]["sasa_time_ms"] = timing.get("sasa_time_ms", "")
                n_success += 1
            else:
                n_failed += 1

            progress.advance(task)

    # Atomic write
    tmp_path = csv_path.with_suffix(".csv.tmp")
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp_path, csv_path)

    if n_failed > 0:
        console.print(
            f"\n[yellow]Warning:[/yellow] {n_failed}/{len(tool_rows)} timing runs failed"
        )

    console.print(
        f"\n[green]Done![/green] {n_success}/{len(tool_rows)} timing runs completed"
    )
    console.print(f"  Updated: {csv_path}")
    return True


if __name__ == "__main__":
    app()

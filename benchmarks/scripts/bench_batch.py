#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Batch SASA benchmark using hyperfine (RustSASA paper methodology).

Flexible batch benchmark script that works with any PDB directory.
Uses hyperfine for timing measurements (warmup + multiple runs).

Usage:
    # E. coli proteome benchmark
    ./benchmarks/scripts/bench_batch.py -i benchmarks/UP000000625_83333_ECOLI_v6/pdb -n ecoli

    # SwissProt benchmark
    ./benchmarks/scripts/bench_batch.py -i benchmarks/swissprot_pdb_v6 -n swissprot --threads 10

    # Single tool test
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --tool zig --runs 1

    # Bitmask variant
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --tool zig_bitmask

    # Multiple tools (skip freesasa)
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --tool zig --tool rustsasa

    # f64 only (skip f32)
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --precision f64

    # Multiple thread counts
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --threads 1,8,10

    # Clear filesystem cache before each run (macOS)
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --prepare 'sudo purge'

    # Sync filesystem buffers before each run
    ./benchmarks/scripts/bench_batch.py -i /path/to/pdb -n test --prepare sync

Output:
    benchmarks/results/batch/<n_points>/<name>/
    ├── config.json             # System info and parameters
    ├── bench_zsasa_f64_8t.json
    ├── bench_zsasa_f32_8t.json
    ├── bench_zsasa_f64_bitmask_8t.json  # --tool zig_bitmask
    ├── bench_freesasa_8t.json
    ├── bench_rustsasa_8t.json
    ├── bench_lahuta_8t.json
    └── bench_lahuta_bitmask_8t.json     # --tool lahuta_bitmask
"""

from __future__ import annotations

import tempfile
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from bench_common import (
    LAHUTA_BITMASK_POINTS,
    Precision,
    check_hyperfine,
    ensure_zsasa_built,
    get_binary_path,
    get_system_info,
    parse_threads,
    print_benchmark_summary,
    quote_path,
    run_benchmark,
    save_config,
)

app = typer.Typer(help="Batch SASA benchmark (hyperfine-based)")
console = Console()


class Tool(str, Enum):
    zig = "zig"
    zig_bitmask = "zig_bitmask"
    freesasa = "freesasa"
    rustsasa = "rustsasa"
    lahuta = "lahuta"
    lahuta_bitmask = "lahuta_bitmask"


ALL_TOOLS = [
    Tool.zig,
    Tool.zig_bitmask,
    Tool.freesasa,
    Tool.rustsasa,
    Tool.lahuta,
    Tool.lahuta_bitmask,
]


def get_root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def get_binary_paths() -> dict[str, Path]:
    """Get paths to tool binaries (all in external/bin/ via setup.sh)."""
    bin_dir = get_root_dir().joinpath("benchmarks", "external", "bin")
    return {
        "zsasa": get_binary_path("zig"),
        "freesasa_batch": bin_dir.joinpath("freesasa_batch"),
        "rustsasa": get_binary_path("rust"),
        "lahuta": get_binary_path("lahuta"),
    }


def run_zig(
    input_dir: Path,
    results_dir: Path,
    thread_counts: list[int],
    n_points: int,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
    use_bitmask: bool = False,
    timeout: int = 600,
    prepare: str | None = None,
    precisions: list[str] | None = None,
) -> list[dict]:
    """Run zsasa benchmarks."""
    zsasa = binaries["zsasa"]
    if not zsasa.exists():
        console.print("[yellow][SKIP] zsasa not found[/]")
        return []

    results = []
    bitmask_suffix = "_bitmask" if use_bitmask else ""

    for precision in precisions or ["f64", "f32"]:
        with tempfile.TemporaryDirectory(prefix=f"zsasa_{precision}_") as tmp:
            out_file = Path(tmp).joinpath("sasa.jsonl")

            for n_threads in thread_counts:
                bitmask_flag = " --use-bitmask" if use_bitmask else ""
                bench_name = f"zsasa_{precision}{bitmask_suffix}_{n_threads}t"
                result = run_benchmark(
                    bench_name,
                    f"{quote_path(zsasa)} batch {quote_path(input_dir)} --format=jsonl -o {quote_path(out_file)} --threads={n_threads} --precision={precision} --n-points={n_points}{bitmask_flag}",
                    results_dir,
                    warmup,
                    runs,
                    dry_run,
                    timeout,
                    prepare,
                )
                if result:
                    results.append({"name": bench_name, **result})

    return results


def run_freesasa(
    input_dir: Path,
    results_dir: Path,
    thread_counts: list[int],
    n_points: int,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
    timeout: int = 600,
    prepare: str | None = None,
) -> list[dict]:
    """Run FreeSASA benchmarks with file-level parallelism."""
    freesasa_batch = binaries["freesasa_batch"]
    if not freesasa_batch.exists():
        console.print("[yellow][SKIP] freesasa_batch not found[/]")
        return []

    results = []

    with tempfile.TemporaryDirectory(prefix="freesasa_") as tmp:
        out_dir = Path(tmp)

        for n_threads in thread_counts:
            result = run_benchmark(
                f"freesasa_{n_threads}t",
                f"{quote_path(freesasa_batch)} {quote_path(input_dir)} {quote_path(out_dir)} --n-threads={n_threads} --n-points={n_points}",
                results_dir,
                warmup,
                runs,
                dry_run,
                timeout,
                prepare,
            )
            if result:
                results.append({"name": f"freesasa_{n_threads}t", **result})

    return results


def run_rustsasa(
    input_dir: Path,
    results_dir: Path,
    thread_counts: list[int],
    n_points: int,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
    timeout: int = 600,
    prepare: str | None = None,
) -> list[dict]:
    """Run RustSASA benchmarks."""
    rustsasa = binaries["rustsasa"]
    if not rustsasa.exists():
        console.print("[yellow][SKIP] RustSASA not found[/]")
        return []

    results = []

    with tempfile.TemporaryDirectory(prefix="rustsasa_") as tmp:
        out_dir = Path(tmp)

        for n_threads in thread_counts:
            result = run_benchmark(
                f"rustsasa_{n_threads}t",
                f"{quote_path(rustsasa)} {quote_path(input_dir)} {quote_path(out_dir)} --format json -t {n_threads} -n {n_points} --allow-vdw-fallback",
                results_dir,
                warmup,
                runs,
                dry_run,
                timeout,
                prepare,
            )
            if result:
                results.append({"name": f"rustsasa_{n_threads}t", **result})

    return results


def run_lahuta(
    input_dir: Path,
    results_dir: Path,
    thread_counts: list[int],
    n_points: int,
    warmup: int,
    runs: int,
    dry_run: bool,
    binaries: dict[str, Path],
    use_bitmask: bool = False,
    single_tool: bool = False,
    timeout: int = 600,
    prepare: str | None = None,
) -> list[dict]:
    """Run Lahuta benchmarks (AF2 PDB only, Shrake-Rupley)."""
    # Resolve to absolute path since the command cd's to a temp directory
    input_dir = input_dir.resolve()
    lahuta = binaries["lahuta"]
    if not lahuta.exists():
        console.print("[yellow][SKIP] lahuta not found[/]")
        return []

    if use_bitmask and n_points not in LAHUTA_BITMASK_POINTS:
        if single_tool:
            console.print(
                f"[red]Error: lahuta bitmask only supports n_points "
                f"{sorted(LAHUTA_BITMASK_POINTS)}, got {n_points}.[/red]"
            )
            raise typer.Exit(1)
        console.print(
            f"[yellow][SKIP] lahuta_bitmask: n_points={n_points} not in "
            f"{sorted(LAHUTA_BITMASK_POINTS)}[/yellow]"
        )
        return []

    results = []
    bitmask_suffix = "_bitmask" if use_bitmask else ""
    bitmask_flag = " --use-bitmask" if use_bitmask else ""

    with tempfile.TemporaryDirectory(prefix="lahuta_") as tmp:
        out_file = Path(tmp).joinpath("sasa.jsonl")

        for n_threads in thread_counts:
            bench_name = f"lahuta{bitmask_suffix}_{n_threads}t"
            result = run_benchmark(
                bench_name,
                (
                    f"cd {quote_path(tmp)} && "
                    f"{quote_path(lahuta)} sasa-sr"
                    f" -d {quote_path(input_dir)}"
                    f" --is_af2_model"
                    f" --points {n_points}"
                    f"{bitmask_flag}"
                    f" -t {n_threads}"
                    f" --progress 0"
                    f" -o {quote_path(out_file)}"
                ),
                results_dir,
                warmup,
                runs,
                dry_run,
                timeout,
                prepare,
            )
            if result:
                results.append({"name": bench_name, **result})

    return results


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
        str,
        typer.Option(
            "--threads",
            "-T",
            help="Thread counts: '1,8,10' or '1-10'",
        ),
    ] = "1,8",
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
    tools: Annotated[
        list[Tool] | None,
        typer.Option(
            "--tool",
            "-t",
            help="Tools to benchmark (can specify multiple: --tool zig --tool rustsasa). Default: all",
        ),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output directory (default: benchmarks/results/batch/<name>)",
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
    timeout: Annotated[
        int,
        typer.Option(
            "--timeout",
            help="Timeout per benchmark in seconds (default: 3600)",
        ),
    ] = 3600,
    prepare: Annotated[
        str | None,
        typer.Option(
            "--prepare",
            "-p",
            help="Shell command to run before each timing run (passed to hyperfine --prepare). "
            "E.g. 'sync' or 'sudo purge' (macOS) to clear filesystem caches.",
        ),
    ] = None,
    precisions: Annotated[
        list[Precision] | None,
        typer.Option(
            "--precision",
            "-P",
            help="Precision to benchmark (can specify multiple: --precision f64 --precision f32). Default: f64,f32",
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

    thread_counts = parse_threads(threads)
    selected_precisions = (
        [p.value for p in precisions] if precisions else ["f64", "f32"]
    )

    # Default to all tools if none specified
    selected_tools = tools if tools else ALL_TOOLS

    # Rebuild zsasa in ReleaseFast mode to ensure benchmarks use optimized code
    zig_tools = {Tool.zig, Tool.zig_bitmask}
    if not dry_run and zig_tools & set(selected_tools):
        ensure_zsasa_built()

    # Set up paths
    root = get_root_dir()
    if output_dir is None:
        results_dir = root.joinpath(
            "benchmarks", "results", "batch", str(n_points), name
        )
    else:
        results_dir = output_dir

    binaries = get_binary_paths()

    # Create directories
    if not dry_run:
        results_dir.mkdir(parents=True, exist_ok=True)

    # Count PDB files in input directory
    n_files = sum(1 for f in input_dir.iterdir() if f.is_file() and f.suffix == ".pdb")

    # Save config
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "input_dir": str(input_dir),
            "n_files": n_files,
            "tools": [t.value for t in selected_tools],
            "thread_counts": thread_counts,
            "warmup": warmup,
            "runs": runs,
            "n_points": n_points,
            "precisions": selected_precisions,
            "prepare": prepare,
        },
    }
    if not dry_run:
        save_config(config, results_dir)

    # Print header
    console.print(f"[bold]=== Batch SASA Benchmark: {name} ===[/]")
    console.print(f"Input: {input_dir}")
    console.print(f"Output: {results_dir}")
    console.print(f"Tools: {', '.join(t.value for t in selected_tools)}")
    prepare_info = f", Prepare: '{prepare}'" if prepare else ""
    precision_info = f", Precision: {','.join(selected_precisions)}"
    console.print(
        f"Warmup: {warmup}, Runs: {runs}, Threads: {thread_counts}, Points: {n_points}{precision_info}{prepare_info}"
    )
    console.print()

    all_results = []

    # Run benchmarks
    if Tool.zig in selected_tools:
        results = run_zig(
            input_dir,
            results_dir,
            thread_counts,
            n_points,
            warmup,
            runs,
            dry_run,
            binaries,
            use_bitmask=False,
            timeout=timeout,
            prepare=prepare,
            precisions=selected_precisions,
        )
        all_results.extend(results)

    if Tool.zig_bitmask in selected_tools:
        results = run_zig(
            input_dir,
            results_dir,
            thread_counts,
            n_points,
            warmup,
            runs,
            dry_run,
            binaries,
            use_bitmask=True,
            timeout=timeout,
            prepare=prepare,
            precisions=selected_precisions,
        )
        all_results.extend(results)

    if Tool.freesasa in selected_tools:
        results = run_freesasa(
            input_dir,
            results_dir,
            thread_counts,
            n_points,
            warmup,
            runs,
            dry_run,
            binaries,
            timeout=timeout,
            prepare=prepare,
        )
        all_results.extend(results)

    if Tool.rustsasa in selected_tools:
        results = run_rustsasa(
            input_dir,
            results_dir,
            thread_counts,
            n_points,
            warmup,
            runs,
            dry_run,
            binaries,
            timeout=timeout,
            prepare=prepare,
        )
        all_results.extend(results)

    if Tool.lahuta in selected_tools:
        results = run_lahuta(
            input_dir,
            results_dir,
            thread_counts,
            n_points,
            warmup,
            runs,
            dry_run,
            binaries,
            use_bitmask=False,
            single_tool=(len(selected_tools) == 1),
            timeout=timeout,
            prepare=prepare,
        )
        all_results.extend(results)

    if Tool.lahuta_bitmask in selected_tools:
        results = run_lahuta(
            input_dir,
            results_dir,
            thread_counts,
            n_points,
            warmup,
            runs,
            dry_run,
            binaries,
            use_bitmask=True,
            single_tool=(len(selected_tools) == 1),
            timeout=timeout,
            prepare=prepare,
        )
        all_results.extend(results)

    # Print summary
    console.print(f"[bold green]=== Done! Results: {results_dir} ===[/]")
    if not dry_run:
        console.print(f"  - config.json")

    if not dry_run and results_dir.exists():
        print_benchmark_summary(results_dir)


if __name__ == "__main__":
    app()

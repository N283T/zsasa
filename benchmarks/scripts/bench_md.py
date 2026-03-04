#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Batch MD trajectory SASA benchmark using hyperfine.

Compares SASA calculation performance across implementations:
- zsasa CLI (Zig, traj mode, f32/f64, with optional bitmask neighbor list)
- zsasa.mdtraj (Python wrapper)
- zsasa.mdanalysis (Python wrapper)
- MDTraj shrake_rupley (native, single-threaded)
- mdsasa-bolt (RustSASA, all cores)

Usage:
    # All tools, default threads
    ./benchmarks/scripts/bench_md.py \\
        --xtc benchmarks/md_data/5vz0_A_protein/5vz0_A_prod_R1_fit.xtc \\
        --pdb benchmarks/md_data/5vz0_A_protein/5vz0_A.pdb \\
        --name 5vz0_R1

    # Specific tools and threads
    ./benchmarks/scripts/bench_md.py \\
        --xtc benchmarks/md_data/5vz0_A_protein/5vz0_A_prod_R1_fit.xtc \\
        --pdb benchmarks/md_data/5vz0_A_protein/5vz0_A.pdb \\
        --name 5vz0_R1 \\
        --tool zig --tool mdtraj \\
        --threads 1,4,8

    # Bitmask variant
    ./benchmarks/scripts/bench_md.py \\
        --xtc traj.xtc --pdb top.pdb \\
        --name test --tool zig_bitmask

    # Quick test with stride
    ./benchmarks/scripts/bench_md.py \\
        --xtc traj.xtc --pdb top.pdb \\
        --name test --runs 1 --warmup 0 --stride 100

    # f64 only (skip f32)
    ./benchmarks/scripts/bench_md.py \\
        --xtc traj.xtc --pdb top.pdb \\
        --name test --precision f64

    # Clear filesystem cache before each run (macOS)
    ./benchmarks/scripts/bench_md.py \\
        --xtc traj.xtc --pdb top.pdb \\
        --name test --prepare 'sudo purge'

Output:
    benchmarks/results/md/<n_points>/<name>/
    ├── config.json
    ├── bench_zig_f32_4t.json
    ├── bench_zig_f64_4t.json
    ├── bench_zig_f32_bitmask_4t.json   # --tool zig_bitmask
    ├── bench_zig_f64_bitmask_4t.json   # --tool zig_bitmask
    ├── bench_zsasa_mdtraj_4t.json
    ├── bench_zsasa_mdanalysis_4t.json
    ├── bench_mdtraj_1t.json
    └── bench_mdsasa_bolt_all.json
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
    Precision,
    check_hyperfine,
    ensure_zsasa_built,
    get_binary_path,
    get_system_info,
    parse_threads,
    print_benchmark_summary,
    run_benchmark,
    save_config,
)

app = typer.Typer(help="MD trajectory SASA benchmark (hyperfine-based)")
console = Console()


class Tool(str, Enum):
    zig = "zig"
    zig_bitmask = "zig_bitmask"
    zsasa_mdtraj = "zsasa_mdtraj"
    zsasa_mdanalysis = "zsasa_mdanalysis"
    mdtraj = "mdtraj"
    mdsasa_bolt = "mdsasa_bolt"


ALL_TOOLS = [
    Tool.zig,
    Tool.zig_bitmask,
    Tool.zsasa_mdtraj,
    Tool.zsasa_mdanalysis,
    Tool.mdtraj,
    Tool.mdsasa_bolt,
]


def get_root_dir() -> Path:
    """Get project root directory."""
    return Path(__file__).parent.parent.parent


def get_runner_path() -> Path:
    """Get path to bench_md_runner.py."""
    return Path(__file__).parent.joinpath("bench_md_runner.py")


def get_trajectory_info(xtc: Path, pdb: Path) -> dict:
    """Get trajectory metadata (atom_count, total_frames).

    Uses mdtraj if available, otherwise returns empty dict.
    """
    try:
        import mdtraj as md

        traj = md.load(str(xtc), top=str(pdb))
        return {"total_frames": traj.n_frames, "atom_count": traj.n_atoms}
    except Exception:
        return {}


def run_zig(
    xtc: Path,
    pdb: Path,
    results_dir: Path,
    thread_counts: list[int],
    warmup: int,
    runs: int,
    dry_run: bool,
    stride: int,
    n_points: int,
    use_bitmask: bool = False,
    timeout: int = 600,
    prepare: str | None = None,
    precisions: list[str] | None = None,
) -> list[dict]:
    """Run zsasa CLI traj benchmarks (precision x threads)."""
    zsasa = get_binary_path("zig")
    if not zsasa.exists():
        console.print("[yellow][SKIP] zsasa binary not found[/]")
        return []

    results = []
    bitmask_suffix = "_bitmask" if use_bitmask else ""

    for precision in precisions or ["f32", "f64"]:
        with tempfile.NamedTemporaryFile(suffix=".csv", prefix="zsasa_") as tmp:
            out_path = tmp.name

            for n_threads in thread_counts:
                bitmask_flag = " --use-bitmask" if use_bitmask else ""
                bench_name = f"zig_{precision}{bitmask_suffix}_{n_threads}t"
                cmd = (
                    f"{zsasa} traj {xtc} {pdb}"
                    f" --include-hydrogens"
                    f" --threads={n_threads}"
                    f" --precision={precision}"
                    f" --stride={stride}"
                    f" --n-points={n_points}"
                    f"{bitmask_flag}"
                    f" -o {out_path} -q"
                )
                result = run_benchmark(
                    bench_name,
                    cmd,
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


def run_python_tool(
    tool_name: str,
    label: str,
    xtc: Path,
    pdb: Path,
    results_dir: Path,
    runner: Path,
    warmup: int,
    runs: int,
    dry_run: bool,
    n_threads: int | None = None,
    stride: int = 1,
    n_points: int = 100,
    timeout: int = 600,
    prepare: str | None = None,
) -> dict | None:
    """Run a Python-based benchmark tool via the runner script."""
    cmd_parts = [
        f"uv run --script {runner}",
        f"--tool {tool_name}",
        f"--xtc {xtc}",
        f"--pdb {pdb}",
        f"--n-points {n_points}",
        f"--stride {stride}",
    ]
    if n_threads is not None:
        cmd_parts.append(f"--threads {n_threads}")

    cmd = " ".join(cmd_parts)

    return run_benchmark(
        label, cmd, results_dir, warmup, runs, dry_run, timeout, prepare
    )


def run_threaded_python_tool(
    tool_name: str,
    xtc: Path,
    pdb: Path,
    results_dir: Path,
    runner: Path,
    thread_counts: list[int],
    warmup: int,
    runs: int,
    dry_run: bool,
    stride: int,
    n_points: int,
    timeout: int = 600,
    prepare: str | None = None,
) -> list[dict]:
    """Run a Python-based benchmark tool across multiple thread counts.

    Used for zsasa_mdtraj and zsasa_mdanalysis which share the same
    loop-over-threads pattern.
    """
    results = []
    for n_threads in thread_counts:
        label = f"{tool_name}_{n_threads}t"
        result = run_python_tool(
            tool_name,
            label,
            xtc,
            pdb,
            results_dir,
            runner,
            warmup,
            runs,
            dry_run,
            n_threads=n_threads,
            stride=stride,
            n_points=n_points,
            timeout=timeout,
            prepare=prepare,
        )
        if result:
            results.append({"name": label, **result})
    return results


def run_mdtraj(
    xtc: Path,
    pdb: Path,
    results_dir: Path,
    runner: Path,
    warmup: int,
    runs: int,
    dry_run: bool,
    stride: int,
    n_points: int,
    timeout: int = 600,
    prepare: str | None = None,
) -> list[dict]:
    """Run native MDTraj benchmark (single-threaded)."""
    label = "mdtraj_1t"
    result = run_python_tool(
        "mdtraj",
        label,
        xtc,
        pdb,
        results_dir,
        runner,
        warmup,
        runs,
        dry_run,
        stride=stride,
        n_points=n_points,
        timeout=timeout,
        prepare=prepare,
    )
    if result:
        return [{"name": label, **result}]
    return []


def run_mdsasa_bolt(
    xtc: Path,
    pdb: Path,
    results_dir: Path,
    runner: Path,
    warmup: int,
    runs: int,
    dry_run: bool,
    stride: int,
    n_points: int,
    timeout: int = 600,
    prepare: str | None = None,
) -> list[dict]:
    """Run mdsasa-bolt benchmark (all cores)."""
    label = "mdsasa_bolt_all"
    result = run_python_tool(
        "mdsasa_bolt",
        label,
        xtc,
        pdb,
        results_dir,
        runner,
        warmup,
        runs,
        dry_run,
        stride=stride,
        n_points=n_points,
        timeout=timeout,
        prepare=prepare,
    )
    if result:
        return [{"name": label, **result}]
    return []


@app.command()
def main(
    xtc: Annotated[
        Path,
        typer.Option(
            "--xtc",
            help="XTC trajectory file",
            exists=True,
            file_okay=True,
            dir_okay=False,
        ),
    ],
    pdb: Annotated[
        Path,
        typer.Option(
            "--pdb",
            help="Topology PDB file",
            exists=True,
            file_okay=True,
            dir_okay=False,
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
            help="Thread counts: '1,4,8' or '1-8' (mdtraj always 1t, mdsasa-bolt always all)",
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
    ] = 1,
    stride: Annotated[
        int,
        typer.Option(
            "--stride",
            "-s",
            help="Frame stride (process every Nth frame)",
        ),
    ] = 1,
    n_points: Annotated[
        int,
        typer.Option(
            "--n-points",
            help="Test points per atom (default 100 for fair comparison)",
        ),
    ] = 100,
    tools: Annotated[
        list[Tool] | None,
        typer.Option(
            "--tool",
            "-t",
            help="Tools to benchmark (can specify multiple). Default: all",
        ),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output directory (default: benchmarks/results/md/<n_points>/<name>)",
        ),
    ] = None,
    timeout: Annotated[
        int,
        typer.Option(
            "--timeout",
            help="Timeout per benchmark in seconds (default: 600)",
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
    precisions: Annotated[
        list[Precision] | None,
        typer.Option(
            "--precision",
            "-P",
            help="Precision for zig benchmarks (can specify multiple: --precision f64 --precision f32). Default: f32,f64",
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
    """Run MD trajectory SASA benchmarks using hyperfine."""
    if not check_hyperfine():
        console.print("[red]Error: hyperfine not found. Please install it first.[/]")
        raise typer.Exit(1)

    thread_counts = parse_threads(threads)
    selected_precisions = (
        [p.value for p in precisions] if precisions else ["f32", "f64"]
    )

    # Set up paths
    root = get_root_dir()
    if output_dir is None:
        results_dir = root.joinpath("benchmarks", "results", "md", str(n_points), name)
    else:
        results_dir = output_dir

    runner = get_runner_path()

    if not dry_run:
        results_dir.mkdir(parents=True, exist_ok=True)

    selected_tools = tools if tools else ALL_TOOLS

    # Rebuild zsasa in ReleaseFast mode to ensure benchmarks use optimized code
    zig_tools = {Tool.zig, Tool.zig_bitmask}
    if not dry_run and zig_tools & set(selected_tools):
        ensure_zsasa_built()

    # Get trajectory metadata
    traj_info = get_trajectory_info(xtc, pdb)

    # Save config (merge with existing if present)
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    config = {
        "timestamp": timestamp,
        "system": get_system_info(),
        "parameters": {
            "name": name,
            "xtc": str(xtc),
            "pdb": str(pdb),
            "tools": [t.value for t in selected_tools],
            "thread_counts": thread_counts,
            "warmup": warmup,
            "runs": runs,
            "stride": stride,
            "n_points": n_points,
            "precisions": selected_precisions,
            "prepare": prepare,
            **traj_info,
        },
    }
    if not dry_run:
        save_config(config, results_dir)

    # Print header
    console.print(f"[bold]=== MD Trajectory SASA Benchmark: {name} ===[/]")
    console.print(f"XTC: {xtc}")
    console.print(f"PDB: {pdb}")
    console.print(f"Tools: {', '.join(t.value for t in selected_tools)}")
    prepare_info = f", Prepare: '{prepare}'" if prepare else ""
    precision_info = f", Precision: {','.join(selected_precisions)}"
    console.print(
        f"Warmup: {warmup}, Runs: {runs}, Threads: {thread_counts}, "
        f"Stride: {stride}, N-points: {n_points}{precision_info}{prepare_info}"
    )
    console.print()

    all_results = []

    # Run benchmarks
    if Tool.zig in selected_tools:
        results = run_zig(
            xtc,
            pdb,
            results_dir,
            thread_counts,
            warmup,
            runs,
            dry_run,
            stride,
            n_points,
            use_bitmask=False,
            timeout=timeout,
            prepare=prepare,
            precisions=selected_precisions,
        )
        all_results.extend(results)

    if Tool.zig_bitmask in selected_tools:
        results = run_zig(
            xtc,
            pdb,
            results_dir,
            thread_counts,
            warmup,
            runs,
            dry_run,
            stride,
            n_points,
            use_bitmask=True,
            timeout=timeout,
            prepare=prepare,
            precisions=selected_precisions,
        )
        all_results.extend(results)

    for python_tool in [Tool.zsasa_mdtraj, Tool.zsasa_mdanalysis]:
        if python_tool in selected_tools:
            results = run_threaded_python_tool(
                python_tool.value,
                xtc,
                pdb,
                results_dir,
                runner,
                thread_counts,
                warmup,
                runs,
                dry_run,
                stride,
                n_points,
                timeout,
                prepare,
            )
            all_results.extend(results)

    if Tool.mdtraj in selected_tools:
        results = run_mdtraj(
            xtc,
            pdb,
            results_dir,
            runner,
            warmup,
            runs,
            dry_run,
            stride,
            n_points,
            timeout,
            prepare,
        )
        all_results.extend(results)

    if Tool.mdsasa_bolt in selected_tools:
        results = run_mdsasa_bolt(
            xtc,
            pdb,
            results_dir,
            runner,
            warmup,
            runs,
            dry_run,
            stride,
            n_points,
            timeout,
            prepare,
        )
        all_results.extend(results)

    # Print summary
    console.print(f"[bold green]=== Done! Results: {results_dir} ===[/]")
    if not dry_run:
        console.print("  - config.json")

    if not dry_run and results_dir.exists():
        print_benchmark_summary(results_dir)


if __name__ == "__main__":
    app()

#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Fair SASA benchmark using pre-processed JSON input.

Compares Zig, FreeSASA C, and RustSASA using identical JSON input files,
eliminating CIF parsing overhead for accurate SASA calculation timing.

JSON inputs are pre-generated with ProtOr radii (protein atoms only).

Requirements:
    - Zig binary: zig build -Doptimize=ReleaseFast
    - FreeSASA C binary with --json-input support (N283T/freesasa-bench)
    - RustSASA binary with -J support (N283T/rustsasa-bench)
"""

from __future__ import annotations

import json
import os
import platform
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="Fair SASA benchmark using JSON input")
console = Console()

# Benchmark structures (same as quick_compare.py)
STRUCTURES = [
    ("1crn", "tiny", "Crambin"),
    ("1ubq", "small", "Ubiquitin"),
    ("1a0q", "medium", "Lipid transfer protein"),
    ("3hhb", "medium", "Hemoglobin"),
    ("1aon", "large", "GroEL-GroES"),
    ("4v6x", "xlarge", "Ribosome"),
]


def get_cpu_count() -> int:
    """Get the number of CPU cores available."""
    return os.cpu_count() or 1


def get_system_info() -> dict[str, str]:
    """Get system information for benchmark context."""
    info = {
        "os": f"{platform.system()} {platform.release()}",
        "arch": platform.machine(),
        "cpu_cores": str(get_cpu_count()),
    }

    system = platform.system()
    if system == "Darwin":
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
                mem_bytes = int(result.stdout.strip())
                info["memory"] = f"{mem_bytes // (1024**3)} GB"
        except Exception:
            pass
    elif system == "Linux":
        try:
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        info["cpu_model"] = line.split(":")[1].strip()
                        break
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal"):
                        mem_kb = int(line.split()[1])
                        info["memory"] = f"{mem_kb // (1024**2)} GB"
                        break
        except Exception:
            pass

    return info


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    pdb_id: str
    n_atoms: int
    algorithm: str
    tool: str
    sasa_time_ms: float
    total_area: float


def run_zig_benchmark(
    json_path: Path,
    algorithm: str = "sr",
    n_threads: int = 0,
) -> tuple[float, float]:
    """Run Zig benchmark with JSON input. Returns (sasa_time_ms, total_area)."""
    import gzip
    import shutil

    zig_binary = Path(__file__).parent.parent / "zig-out" / "bin" / "freesasa_zig"

    if not zig_binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {zig_binary}")

    # Decompress gzipped JSON if needed (Zig doesn't support .gz)
    if str(json_path).endswith(".gz"):
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            input_path = Path(f.name)
        with gzip.open(json_path, "rb") as f_in:
            with open(input_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        cleanup_input = True
    else:
        input_path = json_path
        cleanup_input = False

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    try:
        # Note: No --classifier since JSON already has pre-computed radii
        cmd = [
            str(zig_binary),
            "--timing",
            f"--algorithm={algorithm}",
            str(input_path),
            str(output_path),
        ]
        if n_threads > 0:
            cmd.insert(3, f"--threads={n_threads}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        # Parse SASA calculation time from stderr
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
        output_path.unlink(missing_ok=True)
        if cleanup_input:
            input_path.unlink(missing_ok=True)


def run_freesasa_c_benchmark(
    json_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    n_threads: int = 1,
) -> tuple[float, float]:
    """Run FreeSASA C benchmark with JSON input. Returns (sasa_time_ms, total_area)."""
    fs_c_binary = (
        Path(__file__).parent.parent
        / "benchmarks"
        / "external"
        / "freesasa-bench"
        / "src"
        / "freesasa"
    )

    if not fs_c_binary.exists():
        raise FileNotFoundError(f"FreeSASA C binary not found: {fs_c_binary}")

    cmd = [
        str(fs_c_binary),
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
        raise RuntimeError(f"FreeSASA C failed: {result.stderr}")

    # Parse total area from stdout
    total_area = 0.0
    for line in result.stdout.split("\n"):
        if line.startswith("Total"):
            match = re.search(r"[\d.]+", line)
            if match:
                total_area = float(match.group())
                break

    # Parse SASA time from stderr
    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        match = re.search(r"SASA calculation time:\s*([\d.]+)\s*ms", line)
        if match:
            sasa_time_ms = float(match.group(1))
            break

    return sasa_time_ms, total_area


def run_rustsasa_benchmark(
    json_path: Path,
    n_points: int = 100,
    n_threads: int = 1,
) -> tuple[float, float]:
    """Run RustSASA benchmark with JSON input. Returns (sasa_time_ms, total_area)."""
    rust_binary = (
        Path(__file__).parent.parent
        / "benchmarks"
        / "external"
        / "rustsasa-bench"
        / "target"
        / "release"
        / "rust-sasa"
    )

    if not rust_binary.exists():
        raise FileNotFoundError(f"RustSASA binary not found: {rust_binary}")

    cmd = [
        str(rust_binary),
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

    # Parse SASA_TIME_US from stderr
    sasa_time_ms = 0.0
    for line in result.stderr.split("\n"):
        if "SASA_TIME_US:" in line:
            try:
                us = int(line.split(":")[1].strip())
                sasa_time_ms = us / 1000.0
            except (ValueError, IndexError):
                pass
            break

    # Parse total area from stdout
    total_area = 0.0
    for line in result.stdout.split("\n"):
        match = re.search(r"Total SASA:\s*([\d.]+)", line)
        if match:
            total_area = float(match.group(1))
            break

    return sasa_time_ms, total_area


def print_config(n_runs: int, n_threads: int) -> None:
    """Print benchmark configuration and system info."""
    sys_info = get_system_info()

    table = Table(title="Benchmark Configuration", show_header=False, box=None)
    table.add_column("Key", style="bold cyan")
    table.add_column("Value")

    if "cpu_model" in sys_info:
        table.add_row("CPU", sys_info["cpu_model"])
    table.add_row("CPU Cores", sys_info["cpu_cores"])
    if "memory" in sys_info:
        table.add_row("Memory", sys_info["memory"])
    table.add_row("OS", sys_info["os"])
    table.add_row("Arch", sys_info["arch"])
    table.add_row("", "")

    actual_threads = n_threads if n_threads > 0 else get_cpu_count()
    table.add_row("Runs per benchmark", str(n_runs))
    table.add_row(
        "Threads", f"{actual_threads}" + (" (auto)" if n_threads == 0 else "")
    )
    table.add_row("Input format", "JSON (ProtOr radii, protein only)")

    console.print()
    console.print(table)
    console.print()


def run_benchmarks(
    structures: list[tuple[str, str, str]],
    json_dir: Path,
    n_runs: int = 3,
    n_threads: int = 0,
) -> list[BenchmarkResult]:
    """Run all benchmarks and return results."""
    results: list[BenchmarkResult] = []
    algorithms = ["sr", "lr"]

    for pdb_id, category, description in structures:
        json_path = json_dir / f"{pdb_id}.json.gz"

        if not json_path.exists():
            console.print(f"  [yellow]{pdb_id}: Skipping (JSON not found)[/yellow]")
            continue

        # Get atom count from JSON
        import gzip

        with gzip.open(json_path, "rt") as f:
            data = json.load(f)
            n_atoms = len(data.get("x", []))

        console.print(f"\n[bold]{pdb_id.upper()}[/bold] ({n_atoms:,} atoms)")

        for algo in algorithms:
            # Zig benchmark
            zig_times = []
            zig_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area = run_zig_benchmark(json_path, algo, n_threads)
                    zig_times.append(t)
                    zig_area = area
                except Exception as e:
                    console.print(f"    [red]Zig {algo.upper()}: ERROR[/red] - {e}")
                    break

            if zig_times:
                avg_time = sum(zig_times) / len(zig_times)
                results.append(
                    BenchmarkResult(pdb_id, n_atoms, algo, "zig", avg_time, zig_area)
                )
                algo_str = algo.upper()
                console.print(f"  [green]Zig[/green]   {algo_str}: {avg_time:8.2f} ms")

            # FreeSASA C benchmark
            fsc_times = []
            fsc_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area = run_freesasa_c_benchmark(
                        json_path,
                        algo,
                        n_threads=n_threads if n_threads > 0 else 1,
                    )
                    fsc_times.append(t)
                    fsc_area = area
                except FileNotFoundError:
                    break
                except Exception as e:
                    console.print(f"    [red]FS-C {algo.upper()}: ERROR[/red] - {e}")
                    break

            if fsc_times:
                avg_time = sum(fsc_times) / len(fsc_times)
                results.append(
                    BenchmarkResult(
                        pdb_id, n_atoms, algo, "freesasa-c", avg_time, fsc_area
                    )
                )
                console.print(
                    f"  [yellow]FS-C[/yellow]  {algo.upper()}: {avg_time:8.2f} ms"
                )

            # RustSASA benchmark (SR only)
            if algo == "sr":
                rust_times = []
                rust_area = 0.0
                for _ in range(n_runs):
                    try:
                        t, area = run_rustsasa_benchmark(
                            json_path,
                            n_threads=n_threads if n_threads > 0 else 1,
                        )
                        rust_times.append(t)
                        rust_area = area
                    except FileNotFoundError:
                        break
                    except Exception as e:
                        console.print(f"    [red]Rust SR: ERROR[/red] - {e}")
                        break

                if rust_times:
                    avg_time = sum(rust_times) / len(rust_times)
                    results.append(
                        BenchmarkResult(
                            pdb_id, n_atoms, algo, "rustsasa", avg_time, rust_area
                        )
                    )
                    console.print(f"  [magenta]Rust[/magenta]  SR: {avg_time:8.2f} ms")

    return results


def print_summary(results: list[BenchmarkResult]) -> None:
    """Print benchmark summary table."""
    pdbs = sorted(
        set(r.pdb_id for r in results),
        key=lambda p: next(r.n_atoms for r in results if r.pdb_id == p),
    )

    has_rust = any(r.tool == "rustsasa" for r in results)

    for algo, algo_name in [("sr", "Shrake-Rupley"), ("lr", "Lee-Richards")]:
        show_rust = has_rust and algo == "sr"

        table = Table(
            title=f"[bold]{algo_name}[/bold] (SASA calculation time only)",
            show_header=True,
            header_style="bold cyan",
        )

        table.add_column("PDB", style="bold")
        table.add_column("Atoms", justify="right")
        table.add_column("Zig (ms)", justify="right", style="green")
        if show_rust:
            table.add_column("RustSASA (ms)", justify="right", style="magenta")
        table.add_column("FreeSASA C (ms)", justify="right", style="yellow")
        table.add_column("Zig vs FS-C", justify="right", style="bold cyan")
        if show_rust:
            table.add_column("Zig vs Rust", justify="right", style="bold green")

        for pdb in pdbs:
            pdb_results = [r for r in results if r.pdb_id == pdb]
            n_atoms = pdb_results[0].n_atoms if pdb_results else 0

            zig_r = next(
                (r for r in pdb_results if r.algorithm == algo and r.tool == "zig"),
                None,
            )
            fsc_r = next(
                (
                    r
                    for r in pdb_results
                    if r.algorithm == algo and r.tool == "freesasa-c"
                ),
                None,
            )
            rust_r = (
                next(
                    (
                        r
                        for r in pdb_results
                        if r.algorithm == algo and r.tool == "rustsasa"
                    ),
                    None,
                )
                if show_rust
                else None
            )

            if zig_r and fsc_r:
                zig_time = zig_r.sasa_time_ms
                fsc_time = fsc_r.sasa_time_ms
                rust_time = rust_r.sasa_time_ms if rust_r else None
                speedup_fsc = fsc_time / zig_time if zig_time > 0 else 0
                speedup_rust = (
                    rust_time / zig_time if rust_time and zig_time > 0 else None
                )

                row = [
                    pdb.upper(),
                    f"{n_atoms:,}",
                    f"{zig_time:.2f}",
                ]
                if show_rust:
                    row.append(f"{rust_time:.2f}" if rust_time else "-")
                row.append(f"{fsc_time:.2f}")
                row.append(f"{speedup_fsc:.2f}x")
                if show_rust:
                    row.append(f"{speedup_rust:.2f}x" if speedup_rust else "-")

                table.add_row(*row)

        console.print()
        console.print(table)


@app.command()
def main(
    runs: Annotated[
        int, typer.Option("--runs", "-r", help="Number of runs per benchmark")
    ] = 3,
    threads: Annotated[
        int, typer.Option("--threads", "-t", help="Number of threads (0=auto)")
    ] = 0,
    structure: Annotated[
        str | None,
        typer.Option("--structure", "-s", help="Single structure (e.g., 4v6x)"),
    ] = None,
) -> None:
    """Run fair SASA benchmark using JSON input (no CIF parsing overhead)."""
    structures = list(STRUCTURES)

    if structure:
        structure_lower = structure.lower()
        structures = [(p, c, d) for p, c, d in structures if p == structure_lower]
        if not structures:
            console.print(f"[red]Error:[/red] Unknown structure: {structure}")
            raise typer.Exit(1)

    json_dir = Path(__file__).parent.parent / "benchmarks" / "inputs_json"

    if not json_dir.exists():
        console.print(f"[red]Error:[/red] JSON input directory not found: {json_dir}")
        console.print("Generate JSON inputs first with cif_to_protor_json.py")
        raise typer.Exit(1)

    console.rule("[bold]SASA Benchmark (JSON Input)[/bold]")
    print_config(n_runs=runs, n_threads=threads)

    results = run_benchmarks(structures, json_dir, runs, threads)
    print_summary(results)


if __name__ == "__main__":
    app()

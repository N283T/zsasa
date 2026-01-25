#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0", "typer>=0.9.0"]
# ///
"""Benchmark comparing Zig SASA implementation with FreeSASA C and RustSASA.

Measures SASA calculation performance across multiple structure sizes.
Compares Zig CLI vs FreeSASA C vs RustSASA (all native, with multi-threading).

Requirements:
    - Zig binary built with: zig build -Doptimize=ReleaseFast
    - FreeSASA C binary built with thread support (see README)
    - RustSASA binary built with: cargo build --release --features cli
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

app = typer.Typer(help="Benchmark Zig SASA vs FreeSASA C vs RustSASA")
console = Console()


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

    # Try to get CPU model
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
            # Get memory
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


def print_config(
    n_runs: int,
    n_threads: int,
    n_points: int = 100,
    n_slices: int = 20,
    classifier: str = "protor",
    probe_radius: float = 1.4,
) -> None:
    """Print benchmark configuration and system info."""
    sys_info = get_system_info()

    table = Table(title="Benchmark Configuration", show_header=False, box=None)
    table.add_column("Key", style="bold cyan")
    table.add_column("Value")

    # System info
    if "cpu_model" in sys_info:
        table.add_row("CPU", sys_info["cpu_model"])
    table.add_row("CPU Cores", sys_info["cpu_cores"])
    if "memory" in sys_info:
        table.add_row("Memory", sys_info["memory"])
    table.add_row("OS", sys_info["os"])
    table.add_row("Arch", sys_info["arch"])
    table.add_row("", "")  # Separator

    # Benchmark settings
    actual_threads = n_threads if n_threads > 0 else get_cpu_count()
    table.add_row("Runs per benchmark", str(n_runs))
    table.add_row(
        "Threads", f"{actual_threads}" + (" (auto)" if n_threads == 0 else "")
    )
    table.add_row("SR n_points", str(n_points))
    table.add_row("LR n_slices", str(n_slices))
    table.add_row("Classifier", classifier)
    table.add_row("Probe radius", f"{probe_radius} Å")

    console.print()
    console.print(table)
    console.print()


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    pdb_id: str
    n_atoms: int
    algorithm: str
    tool: str
    time_ms: float
    total_area: float
    sasa_only_ms: float | None = None


def run_zig_benchmark(
    input_path: Path,
    algorithm: str = "sr",
    n_threads: int = 0,
    classifier: str = "protor",
) -> tuple[float, float, dict[str, float]]:
    """Run Zig benchmark with timing.

    Returns (time_ms, total_area, timing_breakdown).
    """
    zig_binary = Path(__file__).parent.parent / "zig-out" / "bin" / "freesasa_zig"

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    try:
        cmd = [
            str(zig_binary),
            "--timing",
            f"--algorithm={algorithm}",
            f"--classifier={classifier}",
            str(input_path),
            str(output_path),
        ]
        if n_threads > 0:
            cmd.insert(3, f"--threads={n_threads}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        # Parse timing from stderr
        timing = {}
        for line in result.stderr.split("\n"):
            if ":" in line and "ms" in line:
                match = re.match(r"\s*(.+?):\s*([\d.]+)\s*ms", line)
                if match:
                    key = match.group(1).strip().lower().replace(" ", "_")
                    timing[key] = float(match.group(2))

        # Get total time
        total_time = timing.get("total", 0.0)

        # Parse output
        with open(output_path) as f:
            data = json.load(f)

        return total_time, data["total_area"], timing

    finally:
        output_path.unlink(missing_ok=True)


def run_freesasa_c_benchmark(
    cif_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    n_threads: int = 1,
    fs_c_binary: Path | None = None,
) -> tuple[float, float, float | None]:
    """Run FreeSASA C benchmark. Returns (time_ms, total_area, sasa_only_ms).

    Requires FreeSASA C binary compiled with thread support.
    The sasa_only_ms is parsed from stderr if available (requires patched binary).
    """
    if fs_c_binary is None:
        # Check new location first, then legacy
        new_path = (
            Path(__file__).parent.parent
            / "benchmarks"
            / "external"
            / "freesasa-bench"
            / "src"
            / "freesasa"
        )
        legacy_path = Path(__file__).parent.parent / "freesasa-c" / "src" / "freesasa"
        if new_path.exists():
            fs_c_binary = new_path
        else:
            fs_c_binary = legacy_path

    if not fs_c_binary.exists():
        raise FileNotFoundError(f"FreeSASA C binary not found: {fs_c_binary}")

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    cmd = [
        str(fs_c_binary),
        "--cif",
        "--radii=protor",
        f"--n-threads={n_threads}",
    ]

    if algorithm == "sr":
        cmd.extend(["--shrake-rupley", f"--resolution={n_points}"])
    else:
        cmd.extend(["--lee-richards", f"--resolution={n_slices}"])

    cmd.append(str(cif_path))

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

    # Parse SASA-only time from stderr (if patched binary)
    sasa_only_ms = None
    for line in result.stderr.split("\n"):
        match = re.search(r"SASA calculation time:\s*([\d.]+)\s*ms", line)
        if match:
            sasa_only_ms = float(match.group(1))
            break

    # Use SASA-only time as the primary time if available
    elapsed_ms = sasa_only_ms if sasa_only_ms is not None else 0.0

    return elapsed_ms, total_area, sasa_only_ms


def run_rustsasa_benchmark(
    cif_path: Path,
    n_points: int = 100,
    n_threads: int = 1,
    rust_binary: Path | None = None,
) -> tuple[float, float, float | None]:
    """Run RustSASA benchmark. Returns (time_ms, total_area, sasa_only_ms).

    Requires RustSASA binary with timing patch (outputs SASA_TIME_US to stderr).
    RustSASA only supports Shrake-Rupley algorithm.
    """
    if rust_binary is None:
        # Check new location first, then legacy
        new_path = (
            Path(__file__).parent.parent
            / "benchmarks"
            / "external"
            / "rustsasa-bench"
            / "target"
            / "release"
            / "rust-sasa"
        )
        legacy_path = (
            Path(__file__).parent.parent
            / "rust-sasa-bench"
            / "target"
            / "release"
            / "rust-sasa"
        )
        if new_path.exists():
            rust_binary = new_path
        else:
            rust_binary = legacy_path

    if not rust_binary.exists():
        raise FileNotFoundError(f"RustSASA binary not found: {rust_binary}")

    if not cif_path.exists():
        raise FileNotFoundError(f"CIF file not found: {cif_path}")

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    try:
        cmd = [
            str(rust_binary),
            str(cif_path),
            str(output_path),
            "-o",
            "atom",
            "-n",
            str(n_points),
            "-t",
            str(n_threads),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            raise RuntimeError(f"RustSASA failed: {result.stderr}")

        # Parse SASA_TIME_US from stderr
        sasa_only_ms = None
        for line in result.stderr.split("\n"):
            if "SASA_TIME_US:" in line:
                try:
                    us = int(line.split(":")[1].strip())
                    sasa_only_ms = us / 1000.0
                except (ValueError, IndexError):
                    pass  # Skip malformed timing output
                break

        # Parse total area from output JSON
        with open(output_path) as f:
            data = json.load(f)
        total_area = sum(data.get("Atom", []))

        return sasa_only_ms or 0.0, total_area, sasa_only_ms

    finally:
        output_path.unlink(missing_ok=True)


def download_cif_if_needed(pdb_id: str, cif_dir: Path) -> Path:
    """Download CIF file from RCSB if not present."""
    cif_path = cif_dir / f"{pdb_id}.cif"
    if cif_path.exists():
        return cif_path

    import urllib.error
    import urllib.request

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    console.print(f"  Downloading {pdb_id}.cif...")
    cif_dir.mkdir(parents=True, exist_ok=True)

    try:
        urllib.request.urlretrieve(url, cif_path)
    except (urllib.error.URLError, urllib.error.HTTPError) as e:
        cif_path.unlink(missing_ok=True)
        raise FileNotFoundError(f"Failed to download {pdb_id}.cif: {e}") from e

    return cif_path


def run_benchmarks(
    structures: list[tuple[str, str, str]],
    base_dir: Path,
    n_runs: int = 3,
    n_threads: int = 0,
    fs_c_path: Path | None = None,
    rust_path: Path | None = None,
) -> list[BenchmarkResult]:
    """Run all benchmarks and return results."""
    results = []
    algorithms = ["sr", "lr"]
    cif_dir = base_dir / "cif"

    for pdb_id, category, description in structures:
        input_path = base_dir / "inputs_protor" / f"{pdb_id}.json"

        if not input_path.exists():
            console.print(f"  {pdb_id}: Skipping (input not found)")
            continue

        # Get atom count
        with open(input_path) as f:
            data = json.load(f)
            n_atoms = len(data.get("x", []))

        # Download CIF for FreeSASA C
        try:
            cif_path = download_cif_if_needed(pdb_id, cif_dir)
        except FileNotFoundError as e:
            console.print(f"  {pdb_id}: Skipping ({e})")
            continue

        console.print(f"\n[bold]{pdb_id.upper()}[/bold] ({n_atoms:,} atoms)")

        for algo in algorithms:
            # Zig benchmark
            zig_times = []
            zig_sasa_times = []
            zig_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area, timing = run_zig_benchmark(
                        input_path, algo, n_threads, "protor"
                    )
                    zig_times.append(t)
                    zig_area = area
                    if "sasa_calculation" in timing:
                        zig_sasa_times.append(timing["sasa_calculation"])
                except Exception as e:
                    console.print(f"    [red]Zig {algo.upper()}: ERROR[/red] - {e}")
                    break

            if zig_times:
                avg_time = sum(zig_times) / len(zig_times)
                avg_sasa_time = (
                    sum(zig_sasa_times) / len(zig_sasa_times)
                    if zig_sasa_times
                    else None
                )
                results.append(
                    BenchmarkResult(
                        pdb_id,
                        n_atoms,
                        algo,
                        "zig",
                        avg_time,
                        zig_area,
                        avg_sasa_time,
                    )
                )
                t = avg_sasa_time or avg_time
                console.print(f"  [green]Zig[/green]   {algo.upper():2s}: {t:8.2f} ms")

            # FreeSASA C benchmark
            fsc_times = []
            fsc_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area, sasa_ms = run_freesasa_c_benchmark(
                        cif_path,
                        algo,
                        n_threads=n_threads if n_threads > 0 else 1,
                        fs_c_binary=fs_c_path,
                    )
                    fsc_times.append(t)
                    fsc_area = area
                except Exception as e:
                    console.print(f"    [red]FS-C {algo.upper()}: ERROR[/red] - {e}")
                    break

            if fsc_times:
                avg_time = sum(fsc_times) / len(fsc_times)
                results.append(
                    BenchmarkResult(
                        pdb_id,
                        n_atoms,
                        algo,
                        "freesasa-c",
                        avg_time,
                        fsc_area,
                        avg_time,  # Already SASA-only from patched binary
                    )
                )
                console.print(
                    f"  [yellow]FS-C[/yellow]  {algo.upper():2s}: {avg_time:8.2f} ms"
                )

            # RustSASA benchmark (SR only - RustSASA doesn't support LR)
            if algo == "sr":
                rust_times = []
                rust_area = 0.0
                for _ in range(n_runs):
                    try:
                        t, area, sasa_ms = run_rustsasa_benchmark(
                            cif_path,
                            n_threads=n_threads if n_threads > 0 else 1,
                            rust_binary=rust_path,
                        )
                        rust_times.append(t)
                        rust_area = area
                    except FileNotFoundError:
                        # RustSASA binary not available, skip silently
                        break
                    except Exception as e:
                        console.print(f"    [red]Rust SR: ERROR[/red] - {e}")
                        break

                if rust_times:
                    avg_time = sum(rust_times) / len(rust_times)
                    results.append(
                        BenchmarkResult(
                            pdb_id,
                            n_atoms,
                            algo,
                            "rustsasa",
                            avg_time,
                            rust_area,
                            avg_time,  # Already SASA-only from timing patch
                        )
                    )
                    console.print(
                        f"  [magenta]Rust[/magenta]  {algo.upper():2s}: {avg_time:8.2f} ms"
                    )

    return results


def print_summary(results: list[BenchmarkResult]) -> None:
    """Print benchmark summary table using rich."""
    pdbs = sorted(
        set(r.pdb_id for r in results),
        key=lambda p: next(r.n_atoms for r in results if r.pdb_id == p),
    )

    # Check if RustSASA results are available
    has_rust = any(r.tool == "rustsasa" for r in results)

    for algo, algo_name in [("sr", "Shrake-Rupley"), ("lr", "Lee-Richards")]:
        # For LR, RustSASA is not available
        show_rust = has_rust and algo == "sr"

        table = Table(
            title=f"[bold]{algo_name}[/bold] (SASA calculation time)",
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
                zig_time = zig_r.sasa_only_ms or zig_r.time_ms
                fsc_time = fsc_r.sasa_only_ms or fsc_r.time_ms
                rust_time = (rust_r.sasa_only_ms or rust_r.time_ms) if rust_r else None
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
        typer.Option(
            "--structure", "-s", help="Single structure to benchmark (e.g., 4v6x)"
        ),
    ] = None,
    fs_c_path: Annotated[
        Path | None, typer.Option("--fs-c-path", help="Path to FreeSASA C binary")
    ] = None,
    rust_path: Annotated[
        Path | None, typer.Option("--rust-path", help="Path to RustSASA binary")
    ] = None,
) -> None:
    """Run SASA benchmarks comparing Zig vs FreeSASA C vs RustSASA."""
    structures = [
        ("1crn", "tiny", "Crambin"),
        ("1ubq", "small", "Ubiquitin"),
        ("1a0q", "medium", "Lipid transfer protein"),
        ("3hhb", "medium", "Hemoglobin"),
        ("1aon", "large", "GroEL-GroES"),
        ("4v6x", "xlarge", "Ribosome"),
    ]

    if structure:
        structure_lower = structure.lower()
        structures = [(p, c, d) for p, c, d in structures if p == structure_lower]
        if not structures:
            console.print(f"[red]Error:[/red] Unknown structure: {structure}")
            raise typer.Exit(1)

    base_dir = Path(__file__).parent.parent / "benchmarks"

    console.rule("[bold]SASA Benchmark: Zig vs FreeSASA C vs RustSASA[/bold]")
    print_config(n_runs=runs, n_threads=threads)

    results = run_benchmarks(structures, base_dir, runs, threads, fs_c_path, rust_path)
    print_summary(results)


if __name__ == "__main__":
    app()

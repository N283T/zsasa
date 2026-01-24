#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["rich>=13.0"]
# ///
"""Benchmark comparing Zig SASA implementation with FreeSASA C.

Measures SASA calculation performance across multiple structure sizes.
Compares Zig CLI vs FreeSASA C (both native, both with multi-threading).

Usage:
    ./scripts/benchmark.py [--runs=N] [--threads=N] [--structure=PDB]
    ./scripts/benchmark.py --fs-c-path=PATH [--threads=N]

Examples:
    ./scripts/benchmark.py                      # Run benchmark
    ./scripts/benchmark.py --runs=5             # 5 runs per benchmark
    ./scripts/benchmark.py --structure=4v6x     # Single structure only
    ./scripts/benchmark.py --threads=4          # Use 4 threads

Requirements:
    - Zig binary built with: zig build -Doptimize=ReleaseFast
    - FreeSASA C binary built with thread support (see README)
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


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
    """Run Zig benchmark with timing. Returns (time_ms, total_area, timing_breakdown)."""
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
        fs_c_binary = Path(__file__).parent.parent / "freesasa-c" / "src" / "freesasa"

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


def download_cif_if_needed(pdb_id: str, cif_dir: Path) -> Path:
    """Download CIF file from RCSB if not present."""
    cif_path = cif_dir / f"{pdb_id}.cif"
    if cif_path.exists():
        return cif_path

    import urllib.error
    import urllib.request

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    print(f"  Downloading {pdb_id}.cif...")
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
) -> list[BenchmarkResult]:
    """Run all benchmarks and return results."""
    results = []
    algorithms = ["sr", "lr"]
    cif_dir = base_dir / "cif"

    for pdb_id, category, description in structures:
        input_path = base_dir / "inputs_protor" / f"{pdb_id}.json"

        if not input_path.exists():
            print(f"  {pdb_id}: Skipping (input not found)")
            continue

        # Get atom count
        with open(input_path) as f:
            data = json.load(f)
            n_atoms = len(data.get("x", []))

        # Download CIF for FreeSASA C
        try:
            cif_path = download_cif_if_needed(pdb_id, cif_dir)
        except FileNotFoundError as e:
            print(f"  {pdb_id}: Skipping ({e})")
            continue

        print(f"\n{pdb_id.upper()} ({n_atoms} atoms):")

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
                    print(f"    Zig {algo.upper()}: ERROR - {e}")
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
                print(f"  Zig   {algo.upper():2s}: {avg_sasa_time or avg_time:8.2f} ms")

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
                    print(f"    FS-C {algo.upper()}: ERROR - {e}")
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
                print(f"  FS-C  {algo.upper():2s}: {avg_time:8.2f} ms")

    return results


def print_summary(results: list[BenchmarkResult]) -> None:
    """Print benchmark summary table using rich."""
    from rich.console import Console
    from rich.table import Table

    console = Console()

    pdbs = sorted(
        set(r.pdb_id for r in results),
        key=lambda p: next(r.n_atoms for r in results if r.pdb_id == p),
    )

    for algo, algo_name in [("sr", "Shrake-Rupley"), ("lr", "Lee-Richards")]:
        table = Table(
            title=f"[bold]{algo_name}[/bold] (SASA calculation time)",
            show_header=True,
            header_style="bold cyan",
        )

        table.add_column("PDB", style="bold")
        table.add_column("Atoms", justify="right")
        table.add_column("Zig (ms)", justify="right", style="green")
        table.add_column("FreeSASA C (ms)", justify="right", style="yellow")
        table.add_column("Speedup", justify="right", style="bold magenta")

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

            if zig_r and fsc_r:
                zig_time = zig_r.sasa_only_ms or zig_r.time_ms
                fsc_time = fsc_r.sasa_only_ms or fsc_r.time_ms
                speedup = fsc_time / zig_time if zig_time > 0 else 0
                table.add_row(
                    pdb.upper(),
                    f"{n_atoms:,}",
                    f"{zig_time:.2f}",
                    f"{fsc_time:.2f}",
                    f"{speedup:.2f}x",
                )

        console.print()
        console.print(table)


def main() -> int:
    n_runs = 3
    n_threads = 0
    structure_filter = None
    fs_c_path = None

    for arg in sys.argv[1:]:
        if arg.startswith("--runs="):
            n_runs = int(arg.split("=")[1])
        elif arg.startswith("--threads="):
            n_threads = int(arg.split("=")[1])
        elif arg.startswith("--structure="):
            structure_filter = arg.split("=")[1].lower()
        elif arg.startswith("--fs-c-path="):
            fs_c_path = Path(arg.split("=")[1])

    structures = [
        ("1crn", "tiny", "Crambin"),
        ("1ubq", "small", "Ubiquitin"),
        ("1a0q", "medium", "Lipid transfer protein"),
        ("3hhb", "medium", "Hemoglobin"),
        ("1aon", "large", "GroEL-GroES"),
        ("4v6x", "xlarge", "Ribosome"),
    ]

    if structure_filter:
        structures = [(p, c, d) for p, c, d in structures if p == structure_filter]
        if not structures:
            print(f"Error: Unknown structure: {structure_filter}")
            return 1

    base_dir = Path(__file__).parent.parent / "benchmarks"

    print("=" * 70)
    print(
        f"SASA Benchmark: Zig vs FreeSASA C (runs={n_runs}, threads={n_threads or 'auto'})"
    )
    print("=" * 70)

    results = run_benchmarks(structures, base_dir, n_runs, n_threads, fs_c_path)
    print_summary(results)

    return 0


if __name__ == "__main__":
    sys.exit(main())

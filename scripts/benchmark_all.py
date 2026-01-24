#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
# ]
# ///
"""Unified benchmark comparing Zig implementation with FreeSASA across all structures.

Measures performance across multiple structure sizes with both algorithms.

Usage:
    ./benchmark_all.py [--runs=N] [--threads=N] [--structure=PDB]

Examples:
    ./benchmark_all.py                      # Run all benchmarks
    ./benchmark_all.py --runs=5             # 5 runs per benchmark
    ./benchmark_all.py --structure=1a0q     # Single structure only
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import freesasa
import gemmi


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    pdb_id: str
    n_atoms: int
    algorithm: str
    tool: str
    time_ms: float
    total_area: float
    sasa_only_ms: float | None = (
        None  # SASA calculation time only (for fair comparison)
    )


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


def run_freesasa_benchmark(
    structure_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
) -> tuple[float, float]:
    """Run FreeSASA benchmark. Returns (time_ms, total_area)."""
    # Load structure
    st = gemmi.read_structure(str(structure_path))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    # Write temporary PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
        tmp_pdb = Path(f.name)

    try:
        st.write_pdb(str(tmp_pdb))

        # Configure FreeSASA
        params = freesasa.Parameters()
        if algorithm == "sr":
            params.setAlgorithm(freesasa.ShrakeRupley)
            params.setNPoints(n_points)
        else:
            params.setAlgorithm(freesasa.LeeRichards)
            params.setNSlices(n_slices)

        # Load structure
        structure = freesasa.Structure(str(tmp_pdb))

        # Time the calculation
        start = time.perf_counter()
        result = freesasa.calc(structure, params)
        elapsed = (time.perf_counter() - start) * 1000  # ms

        return elapsed, result.totalArea()

    finally:
        tmp_pdb.unlink(missing_ok=True)


def run_benchmarks(
    structures: list[tuple[str, str, str]],
    base_dir: Path,
    n_runs: int = 3,
    n_threads: int = 0,
) -> list[BenchmarkResult]:
    """Run all benchmarks and return results."""
    results = []
    algorithms = ["sr", "lr"]

    for pdb_id, category, description in structures:
        input_path = base_dir / "inputs" / f"{pdb_id}.json"
        structure_path = base_dir / "structures" / f"{pdb_id}.cif.gz"

        if not input_path.exists() or not structure_path.exists():
            print(f"  {pdb_id}: Skipping (files not found)")
            continue

        # Get atom count
        with open(input_path) as f:
            data = json.load(f)
            n_atoms = len(data.get("x", []))

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
                    # Extract SASA-only time from timing breakdown
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
                        pdb_id, n_atoms, algo, "zig", avg_time, zig_area, avg_sasa_time
                    )
                )
                sasa_str = f", SASA: {avg_sasa_time:.2f}" if avg_sasa_time else ""
                print(
                    f"  Zig {algo.upper():2s}: {avg_time:8.2f} ms{sasa_str} (area: {zig_area:.2f})"
                )

            # FreeSASA benchmark
            fs_times = []
            fs_area = 0.0
            for _ in range(n_runs):
                try:
                    t, area = run_freesasa_benchmark(structure_path, algo)
                    fs_times.append(t)
                    fs_area = area
                except Exception as e:
                    print(f"    FreeSASA {algo.upper()}: ERROR - {e}")
                    break

            if fs_times:
                avg_time = sum(fs_times) / len(fs_times)
                results.append(
                    BenchmarkResult(
                        pdb_id, n_atoms, algo, "freesasa", avg_time, fs_area
                    )
                )
                print(
                    f"  FS  {algo.upper():2s}: {avg_time:8.2f} ms (area: {fs_area:.2f})"
                )

    return results


def print_summary(results: list[BenchmarkResult]) -> None:
    """Print benchmark summary table."""
    print("\n" + "=" * 80)
    print("BENCHMARK SUMMARY (SASA-only for fair comparison)")
    print("=" * 80)

    # Group by PDB
    pdbs = sorted(
        set(r.pdb_id for r in results),
        key=lambda p: next(r.n_atoms for r in results if r.pdb_id == p),
    )

    print(
        f"\n{'PDB':<8} {'Atoms':>8} {'Algo':<4} {'Zig SASA':>10} {'FS (ms)':>10} {'Speedup':>10}"
    )
    print("-" * 60)

    for pdb in pdbs:
        pdb_results = [r for r in results if r.pdb_id == pdb]
        n_atoms = pdb_results[0].n_atoms if pdb_results else 0

        for algo in ["sr", "lr"]:
            zig_r = next(
                (r for r in pdb_results if r.algorithm == algo and r.tool == "zig"),
                None,
            )
            fs_r = next(
                (
                    r
                    for r in pdb_results
                    if r.algorithm == algo and r.tool == "freesasa"
                ),
                None,
            )

            if zig_r and fs_r:
                # Use SASA-only time for fair comparison (or total if not available)
                zig_time = zig_r.sasa_only_ms if zig_r.sasa_only_ms else zig_r.time_ms
                speedup = fs_r.time_ms / zig_time
                print(
                    f"{pdb:<8} {n_atoms:>8} {algo.upper():<4} "
                    f"{zig_time:>10.2f} {fs_r.time_ms:>10.2f} {speedup:>9.2f}x"
                )


def main() -> int:
    # Parse arguments
    n_runs = 3
    n_threads = 0
    structure_filter = None

    for arg in sys.argv[1:]:
        if arg.startswith("--runs="):
            n_runs = int(arg.split("=")[1])
        elif arg.startswith("--threads="):
            n_threads = int(arg.split("=")[1])
        elif arg.startswith("--structure="):
            structure_filter = arg.split("=")[1].lower()

    # Available structures
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

    print("=" * 80)
    print(f"SASA Benchmark (runs={n_runs}, threads={n_threads or 'auto'})")
    print("=" * 80)

    results = run_benchmarks(structures, base_dir, n_runs, n_threads)
    print_summary(results)

    return 0


if __name__ == "__main__":
    sys.exit(main())

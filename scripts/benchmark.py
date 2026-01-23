#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
# ]
# ///
"""Benchmark freesasa-zig against FreeSASA Python bindings.

Usage:
    ./benchmark.py <input_file> [--runs N] [--zig-binary PATH] [--threads N]

Examples:
    ./benchmark.py examples/1A0Q.cif.gz
    ./benchmark.py examples/input_1a0q.json --runs 5
    ./benchmark.py examples/1A0Q.cif.gz --threads 4
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, stdev

import freesasa
import gemmi


@dataclass
class BenchmarkResult:
    """Benchmark results for a single implementation."""

    name: str
    times: list[float]
    total_area: float
    n_atoms: int

    @property
    def mean_time(self) -> float:
        return mean(self.times)

    @property
    def stdev_time(self) -> float:
        return stdev(self.times) if len(self.times) > 1 else 0.0

    @property
    def min_time(self) -> float:
        return min(self.times)

    @property
    def max_time(self) -> float:
        return max(self.times)


def prepare_inputs(
    input_path: Path,
) -> tuple[Path, Path | None]:
    """Prepare input files for both implementations.

    Returns (json_path, pdb_path) for zig and freesasa respectively.
    """
    # Check if input is already JSON
    if input_path.suffix == ".json":
        # Need to create PDB for FreeSASA
        # This is a limitation - we'd need the original structure
        # For now, just return the JSON path and None for PDB
        return input_path, None

    # Input is a structure file (CIF/PDB)
    st = gemmi.read_structure(str(input_path))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    # Create JSON for Zig
    xs, ys, zs, rs = [], [], [], []
    for cra in st[0].all():
        atom = cra.atom
        xs.append(atom.pos.x)
        ys.append(atom.pos.y)
        zs.append(atom.pos.z)
        rs.append(atom.element.vdw_r)

    json_tmp = tempfile.NamedTemporaryFile(suffix=".json", delete=False)
    json.dump({"x": xs, "y": ys, "z": zs, "r": rs}, open(json_tmp.name, "w"))

    # Create PDB for FreeSASA
    pdb_tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    st.write_pdb(pdb_tmp.name)

    return Path(json_tmp.name), Path(pdb_tmp.name)


def benchmark_zig(
    json_path: Path,
    zig_binary: Path,
    runs: int,
    threads: int | None = None,
) -> BenchmarkResult:
    """Benchmark the Zig implementation."""
    times = []
    total_area = 0.0
    n_atoms = 0

    # Build command
    cmd = [str(zig_binary)]
    if threads is not None:
        cmd.append(f"--threads={threads}")
    cmd.extend([str(json_path)])

    for i in range(runs):
        output_tmp = tempfile.NamedTemporaryFile(suffix=".json", delete=False)

        start = time.perf_counter()
        result = subprocess.run(
            cmd + [output_tmp.name],
            capture_output=True,
            text=True,
        )
        elapsed = time.perf_counter() - start

        if result.returncode != 0:
            print(f"Zig error: {result.stderr}", file=sys.stderr)
            sys.exit(1)

        times.append(elapsed)

        # Parse output on first run
        if i == 0:
            with open(output_tmp.name) as f:
                data = json.load(f)
                total_area = data["total_area"]
                n_atoms = len(data["atom_areas"])

        Path(output_tmp.name).unlink()

    name = "Zig"
    if threads is not None:
        name = f"Zig ({threads} thread{'s' if threads != 1 else ''})"

    return BenchmarkResult(
        name=name,
        times=times,
        total_area=total_area,
        n_atoms=n_atoms,
    )


def benchmark_freesasa(
    pdb_path: Path | None,
    json_path: Path | None,
    runs: int,
) -> BenchmarkResult | None:
    """Benchmark the FreeSASA Python implementation."""
    if pdb_path is None:
        return None

    times = []
    total_area = 0.0
    n_atoms = 0

    for i in range(runs):
        start = time.perf_counter()
        structure = freesasa.Structure(str(pdb_path))
        result = freesasa.calc(structure)
        elapsed = time.perf_counter() - start

        times.append(elapsed)

        if i == 0:
            total_area = result.totalArea()
            n_atoms = result.nAtoms()

    return BenchmarkResult(
        name="FreeSASA (Python)",
        times=times,
        total_area=total_area,
        n_atoms=n_atoms,
    )


def print_results(results: list[BenchmarkResult], runs: int) -> None:
    """Print benchmark results in a formatted table."""
    print(f"\n{'=' * 60}")
    print(f"Benchmark Results ({runs} runs)")
    print(f"{'=' * 60}\n")

    # Table header
    print(
        f"{'Implementation':<20} {'Mean (s)':<12} {'Std (s)':<12} {'Min (s)':<12} {'Max (s)':<12}"
    )
    print("-" * 68)

    for r in results:
        print(
            f"{r.name:<20} {r.mean_time:<12.4f} {r.stdev_time:<12.4f} "
            f"{r.min_time:<12.4f} {r.max_time:<12.4f}"
        )

    print()

    # SASA comparison
    print(f"{'Implementation':<20} {'Atoms':<10} {'Total SASA (Å²)':<20}")
    print("-" * 50)

    for r in results:
        print(f"{r.name:<20} {r.n_atoms:<10} {r.total_area:<20.2f}")

    if len(results) == 2:
        diff = abs(results[0].total_area - results[1].total_area)
        diff_pct = diff / results[1].total_area * 100
        speedup = results[1].mean_time / results[0].mean_time

        print()
        print(f"SASA difference: {diff:.2f} Å² ({diff_pct:.2f}%)")
        print(f"Speedup: {speedup:.2f}x")


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_file", type=Path, help="Input file (CIF, PDB, or JSON)")
    parser.add_argument("--runs", type=int, default=3, help="Number of benchmark runs")
    parser.add_argument(
        "--zig-binary",
        type=Path,
        default=Path("zig-out/bin/freesasa_zig"),
        help="Path to Zig binary",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Number of threads for Zig (default: auto-detect)",
    )
    args = parser.parse_args()

    if not args.input_file.exists():
        print(f"Error: File not found: {args.input_file}", file=sys.stderr)
        return 1

    if not args.zig_binary.exists():
        print(f"Error: Zig binary not found: {args.zig_binary}", file=sys.stderr)
        print("Run 'zig build' first.", file=sys.stderr)
        return 1

    print(f"Input: {args.input_file}")
    print(f"Runs: {args.runs}")
    if args.threads is not None:
        print(f"Threads: {args.threads}")

    # Prepare inputs
    json_path, pdb_path = prepare_inputs(args.input_file)

    results = []

    # Benchmark Zig
    thread_info = f" ({args.threads} threads)" if args.threads else ""
    print(f"\nBenchmarking Zig implementation{thread_info}...")
    zig_result = benchmark_zig(json_path, args.zig_binary, args.runs, args.threads)
    results.append(zig_result)

    # Benchmark FreeSASA (if PDB available)
    if pdb_path:
        print("Benchmarking FreeSASA (Python)...")
        freesasa_result = benchmark_freesasa(pdb_path, json_path, args.runs)
        if freesasa_result:
            results.append(freesasa_result)
        pdb_path.unlink()

    # Cleanup temp JSON if we created it
    if json_path != args.input_file:
        json_path.unlink()

    print_results(results, args.runs)

    return 0


if __name__ == "__main__":
    sys.exit(main())

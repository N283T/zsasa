#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "numpy>=1.24",
# ]
# ///
"""Benchmark freesasa-zig Python bindings vs FreeSASA Python bindings.

This benchmark compares the performance of freesasa-zig Python bindings
against the official FreeSASA Python library, both using direct library calls
(no subprocess overhead).

Usage:
    ./benchmark_python.py [--iterations=5] [--warmup=2]

Examples:
    ./benchmark_python.py                     # Default settings
    ./benchmark_python.py --iterations=10     # More iterations for accuracy
"""

from __future__ import annotations

import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import freesasa
import numpy as np

# Add the python package to path
sys.path.insert(0, str(Path(__file__).parent.parent / "python"))
from freesasa_zig import calculate_sasa  # noqa: E402


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""

    pdb_id: str
    n_atoms: int
    zig_sr_ms: float
    zig_lr_ms: float
    freesasa_sr_ms: float
    freesasa_lr_ms: float
    zig_sr_sasa: float
    zig_lr_sasa: float
    freesasa_sr_sasa: float
    freesasa_lr_sasa: float


def load_input(input_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load benchmark input data."""
    with open(input_path) as f:
        data = json.load(f)

    coords = np.column_stack([data["x"], data["y"], data["z"]])
    radii = np.array(data["r"])
    return coords, radii


def benchmark_zig(
    coords: np.ndarray,
    radii: np.ndarray,
    algorithm: str,
    n_points: int = 100,
    n_slices: int = 20,
    iterations: int = 5,
) -> tuple[float, float]:
    """Benchmark zig SASA calculation. Returns (avg_time_ms, sasa)."""
    times = []
    sasa = 0.0

    for _ in range(iterations):
        start = time.perf_counter()
        if algorithm == "sr":
            result = calculate_sasa(coords, radii, algorithm="sr", n_points=n_points)
        else:
            result = calculate_sasa(coords, radii, algorithm="lr", n_slices=n_slices)
        end = time.perf_counter()
        times.append((end - start) * 1000)
        sasa = result.total_area

    return np.mean(times), sasa


def benchmark_freesasa(
    coords: np.ndarray,
    radii: np.ndarray,
    algorithm: str,
    n_points: int = 100,
    n_slices: int = 20,
    iterations: int = 5,
) -> tuple[float, float]:
    """Benchmark FreeSASA calculation. Returns (avg_time_ms, sasa)."""
    times = []
    sasa = 0.0

    # Flatten coords to (x1, y1, z1, x2, y2, z2, ...)
    coord_flat = coords.flatten().tolist()
    radii_list = radii.tolist()

    # Configure parameters
    params = freesasa.Parameters()
    params.setProbeRadius(1.4)
    if algorithm == "sr":
        params.setAlgorithm(freesasa.ShrakeRupley)
        params.setNPoints(n_points)
    else:
        params.setAlgorithm(freesasa.LeeRichards)
        params.setNSlices(n_slices)

    for _ in range(iterations):
        start = time.perf_counter()
        result = freesasa.calcCoord(coord_flat, radii_list, params)
        end = time.perf_counter()
        times.append((end - start) * 1000)
        sasa = result.totalArea()

    return np.mean(times), sasa


def main() -> int:
    # Parse arguments
    iterations = 5
    warmup = 2

    for arg in sys.argv[1:]:
        if arg.startswith("--iterations="):
            iterations = int(arg.split("=")[1])
        elif arg.startswith("--warmup="):
            warmup = int(arg.split("=")[1])

    # Setup paths
    base_dir = Path(__file__).parent.parent / "benchmarks"
    inputs_dir = base_dir / "inputs"

    structures = [
        ("1crn", "tiny"),
        ("1ubq", "small"),
        ("1a0q", "medium"),
        ("3hhb", "medium"),
        ("1aon", "large"),
    ]

    print("=" * 80)
    print("Python Bindings Benchmark: freesasa-zig vs FreeSASA")
    print(f"Iterations: {iterations}, Warmup: {warmup}")
    print("=" * 80)

    results = []

    for pdb_id, category in structures:
        input_path = inputs_dir / f"{pdb_id}.json"
        if not input_path.exists():
            print(f"{pdb_id}: Input not found, skipping")
            continue

        print(f"\n{pdb_id} ({category}):")
        coords, radii = load_input(input_path)
        n_atoms = len(coords)
        print(f"  Atoms: {n_atoms}")

        # Warmup
        for _ in range(warmup):
            calculate_sasa(coords, radii, algorithm="sr")
            calculate_sasa(coords, radii, algorithm="lr")

        # Benchmark
        zig_sr_ms, zig_sr_sasa = benchmark_zig(coords, radii, "sr", iterations=iterations)
        zig_lr_ms, zig_lr_sasa = benchmark_zig(coords, radii, "lr", iterations=iterations)
        freesasa_sr_ms, freesasa_sr_sasa = benchmark_freesasa(
            coords, radii, "sr", iterations=iterations
        )
        freesasa_lr_ms, freesasa_lr_sasa = benchmark_freesasa(
            coords, radii, "lr", iterations=iterations
        )

        # Calculate speedup
        sr_speedup = freesasa_sr_ms / zig_sr_ms if zig_sr_ms > 0 else 0
        lr_speedup = freesasa_lr_ms / zig_lr_ms if zig_lr_ms > 0 else 0

        print(
            f"  SR: Zig {zig_sr_ms:7.2f}ms vs FreeSASA {freesasa_sr_ms:7.2f}ms ({sr_speedup:.1f}x)"
        )
        print(
            f"  LR: Zig {zig_lr_ms:7.2f}ms vs FreeSASA {freesasa_lr_ms:7.2f}ms ({lr_speedup:.1f}x)"
        )
        print(
            f"  SASA diff: SR {abs(zig_sr_sasa - freesasa_sr_sasa):.2f} Å², LR {abs(zig_lr_sasa - freesasa_lr_sasa):.2f} Å²"
        )

        results.append(
            BenchmarkResult(
                pdb_id=pdb_id,
                n_atoms=n_atoms,
                zig_sr_ms=zig_sr_ms,
                zig_lr_ms=zig_lr_ms,
                freesasa_sr_ms=freesasa_sr_ms,
                freesasa_lr_ms=freesasa_lr_ms,
                zig_sr_sasa=zig_sr_sasa,
                zig_lr_sasa=zig_lr_sasa,
                freesasa_sr_sasa=freesasa_sr_sasa,
                freesasa_lr_sasa=freesasa_lr_sasa,
            )
        )

    # Summary
    if results:
        print("\n" + "=" * 80)
        print("Summary")
        print("=" * 80)
        print(f"{'PDB':<8} {'Atoms':>8} {'Zig SR':>10} {'FS SR':>10} {'SR Speedup':>12}")
        print("-" * 50)

        total_zig_sr = sum(r.zig_sr_ms for r in results)
        total_freesasa_sr = sum(r.freesasa_sr_ms for r in results)

        for r in results:
            speedup = r.freesasa_sr_ms / r.zig_sr_ms if r.zig_sr_ms > 0 else 0
            print(
                f"{r.pdb_id:<8} {r.n_atoms:>8} {r.zig_sr_ms:>9.2f}ms {r.freesasa_sr_ms:>9.2f}ms {speedup:>11.1f}x"
            )

        print("-" * 50)
        overall_speedup = total_freesasa_sr / total_zig_sr if total_zig_sr > 0 else 0
        print(
            f"{'Total':<8} {'':<8} {total_zig_sr:>9.2f}ms {total_freesasa_sr:>9.2f}ms {overall_speedup:>11.1f}x"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())

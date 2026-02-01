#!/usr/bin/env python3
"""Basic SASA calculation example.

This example demonstrates the core SASA calculation functionality using
numpy arrays for coordinates and radii.

Requirements:
    pip install zsasa
"""

import numpy as np

from zsasa import calculate_sasa, get_version


def main() -> None:
    print(f"zsasa version: {get_version()}")
    print("=" * 50)

    # Example 1: Single atom (isolated sphere)
    print("\n1. Single atom SASA")
    print("-" * 30)

    coords = np.array([[0.0, 0.0, 0.0]])
    radii = np.array([1.5])  # van der Waals radius

    result = calculate_sasa(coords, radii, probe_radius=1.4)

    # Expected: 4 * pi * (1.5 + 1.4)^2 = 105.68 A^2
    print(f"Coordinates: {coords[0]}")
    print(f"Radius: {radii[0]} A")
    print(f"Probe radius: 1.4 A")
    print(f"Total SASA: {result.total_area:.2f} A^2")
    print(f"(Expected: ~105.68 A^2 for isolated sphere)")

    # Example 2: Two atoms (partial overlap)
    print("\n2. Two atoms with partial overlap")
    print("-" * 30)

    coords = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],  # 3 A apart
    ])
    radii = np.array([1.5, 1.5])

    result = calculate_sasa(coords, radii, probe_radius=1.4)

    print(f"Atom 1: {coords[0]}, radius={radii[0]} A")
    print(f"Atom 2: {coords[1]}, radius={radii[1]} A")
    print(f"Distance: 3.0 A")
    print(f"Total SASA: {result.total_area:.2f} A^2")
    print(f"Atom 1 SASA: {result.atom_areas[0]:.2f} A^2")
    print(f"Atom 2 SASA: {result.atom_areas[1]:.2f} A^2")

    # Example 3: Algorithm comparison (SR vs LR)
    print("\n3. Algorithm comparison (Shrake-Rupley vs Lee-Richards)")
    print("-" * 30)

    # Random 100-atom system
    np.random.seed(42)
    n_atoms = 100
    coords = np.random.uniform(-10, 10, (n_atoms, 3))
    radii = np.random.uniform(1.2, 2.0, n_atoms)

    # Shrake-Rupley algorithm
    result_sr = calculate_sasa(
        coords, radii,
        algorithm="sr",
        n_points=100,
        probe_radius=1.4,
    )

    # Lee-Richards algorithm
    result_lr = calculate_sasa(
        coords, radii,
        algorithm="lr",
        n_slices=20,
        probe_radius=1.4,
    )

    print(f"Number of atoms: {n_atoms}")
    print(f"Shrake-Rupley (100 points): {result_sr.total_area:.2f} A^2")
    print(f"Lee-Richards (20 slices):   {result_lr.total_area:.2f} A^2")
    print(f"Difference: {abs(result_sr.total_area - result_lr.total_area):.2f} A^2")

    # Example 4: Multi-threaded calculation
    print("\n4. Multi-threaded calculation")
    print("-" * 30)

    import time

    # Larger system for timing
    n_atoms = 1000
    coords = np.random.uniform(-20, 20, (n_atoms, 3))
    radii = np.random.uniform(1.2, 2.0, n_atoms)

    # Single-threaded
    start = time.perf_counter()
    result_1t = calculate_sasa(coords, radii, n_threads=1)
    time_1t = time.perf_counter() - start

    # Multi-threaded (auto-detect)
    start = time.perf_counter()
    result_mt = calculate_sasa(coords, radii, n_threads=0)
    time_mt = time.perf_counter() - start

    print(f"Number of atoms: {n_atoms}")
    print(f"Single-threaded: {time_1t*1000:.2f} ms")
    print(f"Multi-threaded:  {time_mt*1000:.2f} ms")
    print(f"Speedup: {time_1t/time_mt:.2f}x")
    print(f"Results match: {np.allclose(result_1t.atom_areas, result_mt.atom_areas)}")


if __name__ == "__main__":
    main()

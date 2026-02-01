#!/usr/bin/env python3
"""Basic SASA calculation example.

This example demonstrates the core SASA calculation functionality using
numpy arrays for coordinates and radii.

Requirements:
    pip install zsasa

Examples:
    # Run all examples
    python basic_sasa.py

    # Use in your own code
    >>> from zsasa import calculate_sasa
    >>> import numpy as np
    >>> coords = np.array([[0.0, 0.0, 0.0]])
    >>> radii = np.array([1.5])
    >>> result = calculate_sasa(coords, radii)
    >>> print(f"SASA: {result.total_area:.2f} A^2")
"""

import time

import numpy as np

from zsasa import calculate_sasa, calculate_sasa_batch, get_version


def example_single_atom() -> None:
    """Calculate SASA for a single isolated atom.

    For an isolated atom, the SASA equals the surface area of a sphere
    with radius = van_der_waals_radius + probe_radius.

    Formula: SASA = 4 * π * (r_vdw + r_probe)²

    For r_vdw=1.5 Å and r_probe=1.4 Å:
    SASA = 4 * π * (2.9)² ≈ 105.68 Å²
    """
    print("\n1. Single atom SASA")
    print("-" * 30)

    # Define a single atom at the origin
    coords = np.array([[0.0, 0.0, 0.0]])
    radii = np.array([1.5])  # Typical carbon van der Waals radius

    # Calculate SASA with standard water probe radius (1.4 Å)
    result = calculate_sasa(coords, radii, probe_radius=1.4)

    # Display results
    print(f"Coordinates: {coords[0]}")
    print(f"van der Waals radius: {radii[0]} Å")
    print(f"Probe radius: 1.4 Å")
    print(f"Total SASA: {result.total_area:.2f} Å²")
    print(f"(Expected: ~105.68 Å² for isolated sphere)")

    # Verify the calculation
    expected = 4 * np.pi * (1.5 + 1.4) ** 2
    print(f"Theoretical: {expected:.2f} Å²")


def example_two_atoms() -> None:
    """Calculate SASA for two overlapping atoms.

    When atoms are close together, their solvent-accessible surfaces
    overlap, reducing the total SASA compared to isolated atoms.

    The burial effect depends on:
    - Distance between atoms
    - Atomic radii
    - Probe radius
    """
    print("\n2. Two atoms with partial overlap")
    print("-" * 30)

    # Two atoms 3 Å apart (typical C-C bond distance ~1.5 Å)
    coords = np.array(
        [
            [0.0, 0.0, 0.0],  # Atom 1 at origin
            [3.0, 0.0, 0.0],  # Atom 2 along x-axis
        ]
    )
    radii = np.array([1.5, 1.5])  # Same radius for both

    result = calculate_sasa(coords, radii, probe_radius=1.4)

    # Two isolated atoms would have SASA = 2 * 105.68 = 211.36 Å²
    # With overlap, total SASA is reduced
    print(f"Atom 1: {coords[0]}, radius={radii[0]} Å")
    print(f"Atom 2: {coords[1]}, radius={radii[1]} Å")
    print(f"Distance: 3.0 Å")
    print()
    print(f"Per-atom SASA:")
    print(f"  Atom 1: {result.atom_areas[0]:.2f} Å²")
    print(f"  Atom 2: {result.atom_areas[1]:.2f} Å²")
    print(f"Total SASA: {result.total_area:.2f} Å²")
    print(f"(Two isolated atoms would be ~211.4 Å²)")
    print(f"Buried area: {211.4 - result.total_area:.2f} Å²")


def example_algorithm_comparison() -> None:
    """Compare Shrake-Rupley and Lee-Richards algorithms.

    zsasa supports two SASA algorithms:

    Shrake-Rupley (SR):
    - Places test points on sphere surface
    - Counts exposed points
    - n_points controls accuracy (default: 100)
    - Good for small systems

    Lee-Richards (LR):
    - Uses circular slices through atoms
    - Analytically calculates arc lengths
    - n_slices controls accuracy (default: 20)
    - Faster for large systems

    Both algorithms should give similar results when properly configured.
    """
    print("\n3. Algorithm comparison (Shrake-Rupley vs Lee-Richards)")
    print("-" * 30)

    # Create a random 100-atom system for testing
    np.random.seed(42)
    n_atoms = 100
    coords = np.random.uniform(-10, 10, (n_atoms, 3))  # Random positions
    radii = np.random.uniform(1.2, 2.0, n_atoms)  # Random radii

    # Shrake-Rupley with 100 test points per atom
    result_sr = calculate_sasa(
        coords,
        radii,
        algorithm="sr",  # Shrake-Rupley
        n_points=100,  # Test points per atom
        probe_radius=1.4,
    )

    # Lee-Richards with 20 slices per atom
    result_lr = calculate_sasa(
        coords,
        radii,
        algorithm="lr",  # Lee-Richards
        n_slices=20,  # Slices per atom
        probe_radius=1.4,
    )

    print(f"Number of atoms: {n_atoms}")
    print()
    print(f"Shrake-Rupley (100 points): {result_sr.total_area:.2f} Å²")
    print(f"Lee-Richards (20 slices):   {result_lr.total_area:.2f} Å²")
    print()
    diff = abs(result_sr.total_area - result_lr.total_area)
    print(f"Absolute difference: {diff:.2f} Å²")
    print(f"Relative difference: {diff / result_sr.total_area:.2%}")


def example_multithreaded() -> None:
    """Demonstrate multi-threaded SASA calculation.

    zsasa automatically parallelizes calculations across CPU cores:

    n_threads parameter:
    - 0: Auto-detect (use all available cores) [default]
    - 1: Single-threaded (useful for debugging)
    - N: Use exactly N threads

    Speedup depends on:
    - Number of atoms (larger = better scaling)
    - Number of available CPU cores
    - System load
    """
    print("\n4. Multi-threaded calculation")
    print("-" * 30)

    # Larger system to demonstrate threading benefit
    np.random.seed(42)
    n_atoms = 1000
    coords = np.random.uniform(-20, 20, (n_atoms, 3))
    radii = np.random.uniform(1.2, 2.0, n_atoms)

    # Single-threaded for baseline
    start = time.perf_counter()
    result_1t = calculate_sasa(coords, radii, n_threads=1)
    time_1t = time.perf_counter() - start

    # Multi-threaded (auto-detect cores)
    start = time.perf_counter()
    result_mt = calculate_sasa(coords, radii, n_threads=0)
    time_mt = time.perf_counter() - start

    print(f"Number of atoms: {n_atoms}")
    print()
    print(f"Single-threaded: {time_1t * 1000:.2f} ms")
    print(f"Multi-threaded:  {time_mt * 1000:.2f} ms")
    print(f"Speedup: {time_1t / time_mt:.2f}x")
    print()

    # Verify results are identical
    match = np.allclose(result_1t.atom_areas, result_mt.atom_areas)
    print(f"Results match: {match}")


def example_precision_comparison() -> None:
    """Compare f64 (double) and f32 (float) precision modes.

    zsasa supports two precision modes via calculate_sasa_batch:

    f64 (default):
    - Double precision (64-bit)
    - Higher accuracy
    - Standard for most applications

    f32:
    - Single precision (32-bit)
    - Faster on some hardware
    - Sufficient for many use cases
    - Lower memory usage for large systems
    """
    print("\n5. Precision comparison (f64 vs f32)")
    print("-" * 30)

    # Create test system
    np.random.seed(42)
    n_atoms = 500
    coords = np.random.uniform(-15, 15, (n_atoms, 3))
    radii = np.random.uniform(1.2, 2.0, n_atoms)

    # Reshape for batch API: (n_frames, n_atoms, 3)
    coords_batch = coords.reshape(1, n_atoms, 3)

    # Calculate with f64 precision (default)
    result_f64 = calculate_sasa_batch(coords_batch, radii, precision="f64")

    # Calculate with f32 precision
    result_f32 = calculate_sasa_batch(coords_batch, radii, precision="f32")

    print(f"Number of atoms: {n_atoms}")
    print()
    print(f"f64 (double): {result_f64.total_areas[0]:.6f} Å²")
    print(f"f32 (float):  {result_f32.total_areas[0]:.6f} Å²")
    print()

    diff = abs(result_f64.total_areas[0] - result_f32.total_areas[0])
    print(f"Absolute difference: {diff:.6f} Å²")
    print(f"Relative difference: {diff / result_f64.total_areas[0]:.6%}")


def example_batch_calculation() -> None:
    """Process multiple frames in a single call.

    calculate_sasa_batch accepts 3D arrays:
    - Shape: (n_frames, n_atoms, 3)
    - All frames processed efficiently
    - Automatic parallelization across frames

    This is useful for:
    - MD trajectory analysis
    - Conformational ensembles
    - Multiple structure comparison
    """
    print("\n6. Batch calculation (multiple frames)")
    print("-" * 30)

    # Create a multi-frame dataset (simulating MD trajectory)
    np.random.seed(42)
    n_frames = 5
    n_atoms = 100

    # Base coordinates
    base_coords = np.random.uniform(-10, 10, (n_atoms, 3))
    radii = np.random.uniform(1.2, 2.0, n_atoms)

    # Create frames with small perturbations (simulating dynamics)
    coords_batch = np.zeros((n_frames, n_atoms, 3))
    for i in range(n_frames):
        perturbation = np.random.normal(0, 0.5, (n_atoms, 3))
        coords_batch[i] = base_coords + perturbation

    # Calculate SASA for all frames at once
    result = calculate_sasa_batch(coords_batch, radii)

    print(f"Number of frames: {n_frames}")
    print(f"Atoms per frame: {n_atoms}")
    print()
    print("Per-frame total SASA:")
    for i, area in enumerate(result.total_areas):
        print(f"  Frame {i}: {area:.2f} Å²")
    print()
    print(f"Mean SASA: {np.mean(result.total_areas):.2f} Å²")
    print(f"Std SASA:  {np.std(result.total_areas):.2f} Å²")
    print()
    print(f"Per-atom SASA shape: {result.atom_areas.shape}")  # (n_frames, n_atoms)


def main() -> None:
    """Run all SASA calculation examples."""
    print(f"zsasa version: {get_version()}")
    print("=" * 50)
    print("Basic SASA Calculation Examples")
    print("=" * 50)

    # Run each example
    example_single_atom()
    example_two_atoms()
    example_algorithm_comparison()
    example_multithreaded()
    example_precision_comparison()
    example_batch_calculation()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

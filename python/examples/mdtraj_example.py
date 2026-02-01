#!/usr/bin/env python3
"""MDTraj integration example.

This example demonstrates SASA calculation for MD trajectories
using the MDTraj library.

zsasa provides a drop-in replacement for MDTraj's shrake_rupley:
- Same API and output format
- 10-100x faster (Zig + SIMD + multi-threading)
- Automatic parallelization across frames

Requirements:
    pip install zsasa mdtraj

Note:
    This example uses a PDB file as a single-frame "trajectory".
    For real MD analysis, you would use:
        traj = md.load('trajectory.xtc', top='topology.pdb')

Examples:
    # Run all examples
    python mdtraj_example.py

    # Use in your own code
    >>> import mdtraj as md
    >>> from zsasa.mdtraj import compute_sasa
    >>> traj = md.load('trajectory.xtc', top='topology.pdb')
    >>> sasa = compute_sasa(traj)  # Shape: (n_frames, n_atoms)
"""

from pathlib import Path
from typing import Optional

import numpy as np

# Find example structure files
EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def get_structure_file() -> Optional[Path]:
    """Find an available example structure file."""
    pdb_file = EXAMPLES_DIR / "1ubq.pdb"
    if pdb_file.exists():
        return pdb_file
    return None


def basic_sasa() -> None:
    """Calculate SASA for a trajectory.

    The compute_sasa function is a drop-in replacement for
    MDTraj's shrake_rupley function.

    Output:
    - Shape: (n_frames, n_atoms)
    - Units: nm² (same as MDTraj convention)
    """
    print("\n1. Basic SASA calculation")
    print("-" * 30)

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    # Load structure (single frame)
    traj = md.load(str(structure_file))
    print(f"Loaded: {traj.n_frames} frame(s), {traj.n_atoms} atoms")

    # Calculate SASA (returns nm², same as MDTraj convention)
    sasa = compute_sasa(traj)

    print(f"SASA shape: {sasa.shape}")  # (n_frames, n_atoms)
    print(f"Total SASA: {sasa.sum():.4f} nm²")
    print(f"Total SASA: {sasa.sum() * 100:.2f} Å²")  # Convert to Å²


def mdtraj_comparison() -> None:
    """Compare with MDTraj's built-in shrake_rupley.

    zsasa should give very similar results to MDTraj's implementation,
    but significantly faster.
    """
    print("\n2. Comparison with MDTraj's shrake_rupley")
    print("-" * 30)

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    traj = md.load(str(structure_file))

    # MDTraj's built-in (slower Python implementation)
    sasa_mdtraj = md.shrake_rupley(traj)

    # zsasa (fast Zig implementation)
    sasa_zsasa = compute_sasa(traj)

    print(f"MDTraj total: {sasa_mdtraj.sum():.6f} nm²")
    print(f"zsasa total:  {sasa_zsasa.sum():.6f} nm²")
    print()

    diff = abs(sasa_mdtraj.sum() - sasa_zsasa.sum())
    print(f"Absolute difference: {diff:.6f} nm²")

    # Per-atom correlation
    corr = np.corrcoef(sasa_mdtraj.flatten(), sasa_zsasa.flatten())[0, 1]
    print(f"Per-atom correlation: {corr:.6f}")


def algorithm_options() -> None:
    """Use different SASA algorithms.

    Shrake-Rupley (SR):
    - Point-based method
    - n_points controls accuracy (default: 100)

    Lee-Richards (LR):
    - Slice-based method
    - n_slices controls accuracy (default: 20)
    """
    print("\n3. Algorithm options")
    print("-" * 30)

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    traj = md.load(str(structure_file))

    # Shrake-Rupley with different point counts
    sasa_sr_100 = compute_sasa(traj, algorithm="sr", n_points=100)
    sasa_sr_500 = compute_sasa(traj, algorithm="sr", n_points=500)

    # Lee-Richards
    sasa_lr = compute_sasa(traj, algorithm="lr", n_slices=20)

    print(f"SR (100 points): {sasa_sr_100.sum():.6f} nm²")
    print(f"SR (500 points): {sasa_sr_500.sum():.6f} nm²")
    print(f"LR (20 slices):  {sasa_lr.sum():.6f} nm²")


def per_residue_sasa() -> None:
    """Aggregate per-atom SASA to per-residue.

    MDTraj provides topology information to map atoms to residues.
    """
    print("\n4. Per-residue SASA")
    print("-" * 30)

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    traj = md.load(str(structure_file))
    sasa = compute_sasa(traj)

    # Get atom-to-residue mapping
    topology = traj.topology
    residue_indices = np.array([atom.residue.index for atom in topology.atoms])

    # Aggregate per-atom SASA to per-residue
    n_residues = topology.n_residues
    residue_sasa = np.zeros(n_residues)

    for i, res_idx in enumerate(residue_indices):
        residue_sasa[res_idx] += sasa[0, i]

    # Show first 10 residues
    print(f"{'Residue':<12} {'SASA (nm²)':>12} {'SASA (Å²)':>12}")
    print("-" * 38)

    for i in range(min(10, n_residues)):
        res = topology.residue(i)
        print(f"{res.name}{res.resSeq:<8} {residue_sasa[i]:>12.6f} {residue_sasa[i]*100:>12.2f}")


def selection() -> None:
    """Calculate SASA for atom selections.

    MDTraj supports selection syntax to extract subsets of atoms.
    """
    print("\n5. Selection-based SASA")
    print("-" * 30)

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    traj = md.load(str(structure_file))

    # Protein heavy atoms only (no hydrogens)
    protein_atoms = traj.topology.select("protein and not element H")
    print(f"Protein heavy atoms: {len(protein_atoms)}")

    # Extract subset and calculate SASA
    traj_protein = traj.atom_slice(protein_atoms)
    sasa_protein = compute_sasa(traj_protein)
    print(f"Protein SASA: {sasa_protein.sum():.6f} nm²")

    # Backbone atoms
    backbone_atoms = traj.topology.select("backbone")
    print(f"Backbone atoms: {len(backbone_atoms)}")

    traj_backbone = traj.atom_slice(backbone_atoms)
    sasa_backbone = compute_sasa(traj_backbone)
    print(f"Backbone SASA: {sasa_backbone.sum():.6f} nm²")


def multiframe_trajectory() -> None:
    """Process multi-frame trajectories.

    zsasa automatically parallelizes across frames for
    efficient trajectory processing.
    """
    print("\n6. Multi-frame trajectory (simulated)")
    print("-" * 30)

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    traj = md.load(str(structure_file))

    # Create a fake multi-frame trajectory by adding small perturbations
    n_frames = 10
    coords_original = traj.xyz[0].copy()
    xyz_perturbed = np.zeros((n_frames, traj.n_atoms, 3))

    np.random.seed(42)
    for i in range(n_frames):
        # Add 0.01 nm (0.1 Å) RMS perturbation
        perturbation = np.random.normal(0, 0.01, coords_original.shape)
        xyz_perturbed[i] = coords_original + perturbation

    # Create new trajectory
    traj_multi = md.Trajectory(xyz_perturbed, traj.topology)
    print(f"Created trajectory: {traj_multi.n_frames} frames")

    # Calculate SASA for all frames
    sasa_multi = compute_sasa(traj_multi)
    print(f"SASA shape: {sasa_multi.shape}")

    # Per-frame total SASA
    total_per_frame = sasa_multi.sum(axis=1)
    print()
    print("Per-frame total SASA (nm²):")
    print(f"  Mean: {total_per_frame.mean():.6f}")
    print(f"  Std:  {total_per_frame.std():.6f}")
    print(f"  Min:  {total_per_frame.min():.6f}")
    print(f"  Max:  {total_per_frame.max():.6f}")


def performance_note() -> None:
    """Performance comparison notes.

    zsasa provides significant speedups over MDTraj's shrake_rupley.
    """
    print("\n7. Performance note")
    print("-" * 30)

    print("zsasa uses optimized Zig code with SIMD and multi-threading.")
    print()
    print("Typical speedups over MDTraj's shrake_rupley:")
    print("  - Small structures (< 1000 atoms): 5-10x")
    print("  - Medium structures (1000-10000 atoms): 10-50x")
    print("  - Large structures (> 10000 atoms): 50-100x")
    print()
    print("For batch processing of many frames, zsasa automatically")
    print("uses parallel processing across frames.")


def main() -> None:
    """Run all MDTraj integration examples."""
    print("MDTraj Integration Examples")
    print("=" * 50)

    structure_file = get_structure_file()
    if not structure_file:
        print(f"Example file not found: {EXAMPLES_DIR / '1ubq.pdb'}")
        return

    print(f"Using structure file: {structure_file.name}")

    # Run each example
    basic_sasa()
    mdtraj_comparison()
    algorithm_options()
    per_residue_sasa()
    selection()
    multiframe_trajectory()
    performance_note()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

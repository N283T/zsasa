#!/usr/bin/env python3
"""MDTraj integration example.

This example demonstrates SASA calculation for MD trajectories
using the MDTraj library.

Requirements:
    pip install zsasa mdtraj

Note:
    This example uses a PDB file as a single-frame "trajectory".
    For real MD analysis, you would use:
        traj = md.load('trajectory.xtc', top='topology.pdb')
"""

from pathlib import Path

import numpy as np


def main() -> None:
    print("MDTraj Integration Examples")
    print("=" * 50)

    # Find example structure file
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    pdb_file = examples_dir / "1ubq.pdb"

    if not pdb_file.exists():
        print(f"Example file not found: {pdb_file}")
        return

    print(f"Using structure file: {pdb_file.name}")

    import mdtraj as md

    from zsasa.mdtraj import compute_sasa

    # Example 1: Basic SASA calculation
    print("\n1. Basic SASA calculation")
    print("-" * 30)

    # Load structure (single frame)
    traj = md.load(str(pdb_file))
    print(f"Loaded: {traj.n_frames} frame(s), {traj.n_atoms} atoms")

    # Calculate SASA (returns nm^2, same as MDTraj convention)
    sasa = compute_sasa(traj)

    print(f"SASA shape: {sasa.shape}")  # (n_frames, n_atoms)
    print(f"Total SASA: {sasa.sum():.2f} nm^2")
    print(f"Total SASA: {sasa.sum() * 100:.2f} A^2")  # Convert to A^2

    # Example 2: Compare with MDTraj's built-in SASA
    print("\n2. Comparison with MDTraj's shrake_rupley")
    print("-" * 30)

    # MDTraj's built-in (slower Python implementation)
    sasa_mdtraj = md.shrake_rupley(traj)

    # zsasa (fast Zig implementation)
    sasa_zsasa = compute_sasa(traj)

    print(f"MDTraj total: {sasa_mdtraj.sum():.4f} nm^2")
    print(f"zsasa total:  {sasa_zsasa.sum():.4f} nm^2")
    print(f"Difference:   {abs(sasa_mdtraj.sum() - sasa_zsasa.sum()):.4f} nm^2")

    # Per-atom correlation
    corr = np.corrcoef(sasa_mdtraj.flatten(), sasa_zsasa.flatten())[0, 1]
    print(f"Per-atom correlation: {corr:.6f}")

    # Example 3: Algorithm options
    print("\n3. Algorithm options")
    print("-" * 30)

    # Shrake-Rupley with different point counts
    sasa_sr_100 = compute_sasa(traj, algorithm="sr", n_points=100)
    sasa_sr_500 = compute_sasa(traj, algorithm="sr", n_points=500)

    # Lee-Richards
    sasa_lr = compute_sasa(traj, algorithm="lr", n_slices=20)

    print(f"SR (100 points): {sasa_sr_100.sum():.4f} nm^2")
    print(f"SR (500 points): {sasa_sr_500.sum():.4f} nm^2")
    print(f"LR (20 slices):  {sasa_lr.sum():.4f} nm^2")

    # Example 4: Per-residue SASA
    print("\n4. Per-residue SASA")
    print("-" * 30)

    # Get atom-to-residue mapping
    topology = traj.topology
    residue_indices = np.array([atom.residue.index for atom in topology.atoms])

    # Aggregate per-atom SASA to per-residue
    n_residues = topology.n_residues
    residue_sasa = np.zeros(n_residues)

    for i, res_idx in enumerate(residue_indices):
        residue_sasa[res_idx] += sasa[0, i]

    # Show first 10 residues
    print(f"{'Residue':<12} {'SASA (nm^2)':>12} {'SASA (A^2)':>12}")
    print("-" * 38)
    for i in range(min(10, n_residues)):
        res = topology.residue(i)
        print(f"{res.name}{res.resSeq:<8} {residue_sasa[i]:>12.4f} {residue_sasa[i]*100:>12.2f}")

    # Example 5: Selection-based SASA
    print("\n5. Selection-based SASA")
    print("-" * 30)

    # Protein heavy atoms only
    protein_atoms = traj.topology.select("protein and not element H")
    print(f"Protein heavy atoms: {len(protein_atoms)}")

    # Extract subset
    traj_protein = traj.atom_slice(protein_atoms)
    sasa_protein = compute_sasa(traj_protein)
    print(f"Protein SASA: {sasa_protein.sum():.4f} nm^2")

    # Backbone atoms
    backbone_atoms = traj.topology.select("backbone")
    print(f"Backbone atoms: {len(backbone_atoms)}")

    traj_backbone = traj.atom_slice(backbone_atoms)
    sasa_backbone = compute_sasa(traj_backbone)
    print(f"Backbone SASA: {sasa_backbone.sum():.4f} nm^2")

    # Example 6: Multi-frame trajectory simulation
    print("\n6. Multi-frame trajectory (simulated)")
    print("-" * 30)

    # Create a fake multi-frame trajectory by adding small perturbations
    n_frames = 10
    coords_original = traj.xyz[0].copy()
    xyz_perturbed = np.zeros((n_frames, traj.n_atoms, 3))

    np.random.seed(42)
    for i in range(n_frames):
        perturbation = np.random.normal(0, 0.01, coords_original.shape)  # 0.01 nm = 0.1 A
        xyz_perturbed[i] = coords_original + perturbation

    # Create new trajectory
    traj_multi = md.Trajectory(xyz_perturbed, traj.topology)
    print(f"Created trajectory: {traj_multi.n_frames} frames")

    # Calculate SASA for all frames
    sasa_multi = compute_sasa(traj_multi)
    print(f"SASA shape: {sasa_multi.shape}")

    # Per-frame total SASA
    total_per_frame = sasa_multi.sum(axis=1)
    print(f"\nPer-frame total SASA (nm^2):")
    print(f"  Mean: {total_per_frame.mean():.4f}")
    print(f"  Std:  {total_per_frame.std():.4f}")
    print(f"  Min:  {total_per_frame.min():.4f}")
    print(f"  Max:  {total_per_frame.max():.4f}")

    # Example 7: Performance comparison
    print("\n7. Performance note")
    print("-" * 30)
    print("zsasa uses optimized Zig code with SIMD and multi-threading.")
    print("For large trajectories, it can be 10-100x faster than MDTraj's")
    print("pure Python implementation.")
    print()
    print("For batch processing of many frames, zsasa automatically")
    print("uses parallel processing across frames.")


if __name__ == "__main__":
    main()

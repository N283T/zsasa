#!/usr/bin/env python3
"""MDAnalysis integration example.

This example demonstrates SASA calculation for MD trajectories
using the MDAnalysis library with zsasa's SASAAnalysis class.

Requirements:
    pip install zsasa MDAnalysis

Note:
    This example uses a PDB file as a single-frame "trajectory".
    For real MD analysis, you would use:
        u = mda.Universe('topology.pdb', 'trajectory.xtc')
"""

from pathlib import Path

import numpy as np


def main() -> None:
    print("MDAnalysis Integration Examples")
    print("=" * 50)

    # Find example structure file
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    pdb_file = examples_dir / "1ubq.pdb"

    if not pdb_file.exists():
        print(f"Example file not found: {pdb_file}")
        return

    print(f"Using structure file: {pdb_file.name}")

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    # Example 1: Basic SASA calculation
    print("\n1. Basic SASA calculation")
    print("-" * 30)

    # Load structure
    u = mda.Universe(str(pdb_file))
    print(f"Loaded: {u.trajectory.n_frames} frame(s), {u.atoms.n_atoms} atoms")

    # Create and run SASA analysis
    sasa = SASAAnalysis(u, select="protein")
    sasa.run()

    # Results are in A^2 (MDAnalysis convention)
    print(f"Total SASA: {sasa.results.total_area[0]:.2f} A^2")
    print(f"Per-atom SASA shape: {sasa.results.atom_area.shape}")

    # Example 2: Selection-based analysis
    print("\n2. Selection-based analysis")
    print("-" * 30)

    # All atoms
    sasa_all = SASAAnalysis(u, select="all")
    sasa_all.run()
    print(f"All atoms:    {sasa_all.results.total_area[0]:.2f} A^2")

    # Protein only (recommended for most analyses)
    sasa_protein = SASAAnalysis(u, select="protein")
    sasa_protein.run()
    print(f"Protein only: {sasa_protein.results.total_area[0]:.2f} A^2")

    # Example 3: Per-residue analysis
    print("\n3. Per-residue SASA")
    print("-" * 30)

    # Run analysis on protein
    sasa = SASAAnalysis(u, select="protein")
    sasa.run()

    # Get atom group for mapping
    protein = u.select_atoms("protein")

    # Aggregate to residues
    residue_sasa = {}
    for atom, area in zip(protein.atoms, sasa.results.atom_area[0]):
        res_key = (atom.resname, atom.resid, atom.segid)
        if res_key not in residue_sasa:
            residue_sasa[res_key] = 0.0
        residue_sasa[res_key] += area

    # Show first 10 residues
    print(f"{'Residue':<12} {'SASA (A^2)':>12}")
    print("-" * 26)
    for i, ((resname, resid, segid), sasa_val) in enumerate(residue_sasa.items()):
        if i >= 10:
            break
        print(f"{resname}{resid:<8} {sasa_val:>12.2f}")

    # Example 4: Algorithm options (passed to run())
    print("\n4. Algorithm options")
    print("-" * 30)

    # Shrake-Rupley with different point counts
    sasa_sr = SASAAnalysis(u, select="protein")
    sasa_sr.run(algorithm="sr", n_points=100)
    print(f"SR (100 points): {sasa_sr.results.total_area[0]:.2f} A^2")

    sasa_sr_high = SASAAnalysis(u, select="protein")
    sasa_sr_high.run(algorithm="sr", n_points=500)
    print(f"SR (500 points): {sasa_sr_high.results.total_area[0]:.2f} A^2")

    # Lee-Richards
    sasa_lr = SASAAnalysis(u, select="protein")
    sasa_lr.run(algorithm="lr", n_slices=20)
    print(f"LR (20 slices):  {sasa_lr.results.total_area[0]:.2f} A^2")

    # Example 5: Probe radius effect (passed to run())
    print("\n5. Probe radius effect")
    print("-" * 30)

    for probe in [1.2, 1.4, 1.6, 2.0]:
        sasa = SASAAnalysis(u, select="protein")
        sasa.run(probe_radius=probe)
        print(f"Probe {probe:.1f} A: {sasa.results.total_area[0]:.2f} A^2")

    # Example 6: Multi-frame trajectory (simulated)
    print("\n6. Multi-frame trajectory (simulated)")
    print("-" * 30)

    # MDAnalysis doesn't easily support in-memory trajectory modification
    # So we'll demonstrate the API usage pattern

    print("For real trajectories, use:")
    print()
    print("  u = mda.Universe('topology.pdb', 'trajectory.xtc')")
    print("  sasa = SASAAnalysis(u, select='protein')")
    print("  sasa.run()")
    print()
    print("  # Per-frame total SASA")
    print("  print(sasa.results.total_area)  # Shape: (n_frames,)")
    print()
    print("  # Per-frame, per-atom SASA")
    print("  print(sasa.results.atom_area)  # Shape: (n_frames, n_atoms)")

    # Example 7: Chain analysis
    print("\n7. Chain/Segment analysis")
    print("-" * 30)

    # Get unique segments
    segments = u.segments
    print(f"Number of segments: {len(segments)}")

    for seg in segments:
        sasa = SASAAnalysis(u, select=f"segid {seg.segid}")
        sasa.run()
        n_atoms = len(u.select_atoms(f"segid {seg.segid}"))
        print(f"Segment {seg.segid}: {n_atoms} atoms, {sasa.results.total_area[0]:.2f} A^2")

    # Example 8: Using compute_sasa convenience function
    print("\n8. Using compute_sasa convenience function")
    print("-" * 30)

    from zsasa.mdanalysis import compute_sasa

    # Simple usage - returns per-atom SASA
    atom_sasa = compute_sasa(u, select="protein")
    print(f"Per-atom SASA shape: {atom_sasa.shape}")
    print(f"Total SASA: {atom_sasa.sum():.2f} A^2")

    # Get per-residue SASA
    residue_sasa = compute_sasa(u, select="protein", mode="residue")
    print(f"Per-residue SASA shape: {residue_sasa.shape}")

    # Get just total SASA
    total_sasa = compute_sasa(u, select="protein", mode="total")
    print(f"Total SASA (scalar): {total_sasa[0]:.2f} A^2")


if __name__ == "__main__":
    main()

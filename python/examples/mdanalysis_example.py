#!/usr/bin/env python3
"""MDAnalysis integration example.

This example demonstrates SASA calculation for MD trajectories
using the MDAnalysis library with zsasa's SASAAnalysis class.

zsasa integrates with MDAnalysis via:
- SASAAnalysis: AnalysisBase subclass for trajectory analysis
- compute_sasa: Convenience function for quick calculations

Requirements:
    pip install zsasa MDAnalysis

Note:
    This example uses a PDB file as a single-frame "trajectory".
    For real MD analysis, you would use:
        u = mda.Universe('topology.pdb', 'trajectory.xtc')

Examples:
    # Run all examples
    python mdanalysis_example.py

    # Use in your own code
    >>> import MDAnalysis as mda
    >>> from zsasa.mdanalysis import SASAAnalysis
    >>> u = mda.Universe('topology.pdb', 'trajectory.xtc')
    >>> sasa = SASAAnalysis(u, select='protein')
    >>> sasa.run()
    >>> print(sasa.results.total_area)  # Per-frame total SASA
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
    """Calculate SASA using SASAAnalysis class.

    SASAAnalysis follows MDAnalysis conventions:
    - Create analysis object with Universe and selection
    - Call run() to process trajectory
    - Access results via results attribute

    Results (in Å², MDAnalysis convention):
    - results.total_area: Per-frame total SASA, shape (n_frames,)
    - results.atom_area: Per-frame per-atom SASA, shape (n_frames, n_atoms)
    """
    print("\n1. Basic SASA calculation")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    # Load structure
    u = mda.Universe(str(structure_file))
    print(f"Loaded: {u.trajectory.n_frames} frame(s), {u.atoms.n_atoms} atoms")

    # Create and run SASA analysis
    sasa = SASAAnalysis(u, select="protein")
    sasa.run()

    # Results are in Å² (MDAnalysis convention)
    print(f"Total SASA: {sasa.results.total_area[0]:.2f} Å²")
    print(f"Per-atom SASA shape: {sasa.results.atom_area.shape}")


def selection_based() -> None:
    """Use different atom selections.

    MDAnalysis selection syntax allows flexible atom selection:
    - "protein": Standard amino acids
    - "all": All atoms
    - "backbone": Backbone atoms only
    - "resname ALA": Specific residues
    """
    print("\n2. Selection-based analysis")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    u = mda.Universe(str(structure_file))

    # All atoms
    sasa_all = SASAAnalysis(u, select="all")
    sasa_all.run()
    print(f"All atoms:    {sasa_all.results.total_area[0]:.2f} Å²")

    # Protein only (recommended for most analyses)
    sasa_protein = SASAAnalysis(u, select="protein")
    sasa_protein.run()
    print(f"Protein only: {sasa_protein.results.total_area[0]:.2f} Å²")


def per_residue_sasa() -> None:
    """Aggregate per-atom SASA to per-residue.

    Map atom-level SASA to residue-level for analysis.
    """
    print("\n3. Per-residue SASA")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    u = mda.Universe(str(structure_file))

    # Run analysis on protein
    sasa = SASAAnalysis(u, select="protein")
    sasa.run()

    # Get atom group for mapping
    protein = u.select_atoms("protein")

    # Aggregate to residues
    residue_sasa: dict[tuple[str, int, str], float] = {}
    for atom, area in zip(protein.atoms, sasa.results.atom_area[0]):
        res_key = (atom.resname, atom.resid, atom.segid)
        if res_key not in residue_sasa:
            residue_sasa[res_key] = 0.0
        residue_sasa[res_key] += area

    # Show first 10 residues
    print(f"{'Residue':<12} {'SASA (Å²)':>12}")
    print("-" * 26)

    for i, ((resname, resid, segid), sasa_val) in enumerate(residue_sasa.items()):
        if i >= 10:
            break
        print(f"{resname}{resid:<8} {sasa_val:>12.2f}")


def algorithm_options() -> None:
    """Use different SASA algorithms.

    Algorithm options are passed to run():
    - algorithm: "sr" (Shrake-Rupley) or "lr" (Lee-Richards)
    - n_points: Points per atom for SR (default: 100)
    - n_slices: Slices per atom for LR (default: 20)
    """
    print("\n4. Algorithm options")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    u = mda.Universe(str(structure_file))

    # Shrake-Rupley with different point counts
    sasa_sr = SASAAnalysis(u, select="protein")
    sasa_sr.run(algorithm="sr", n_points=100)
    print(f"SR (100 points): {sasa_sr.results.total_area[0]:.2f} Å²")

    sasa_sr_high = SASAAnalysis(u, select="protein")
    sasa_sr_high.run(algorithm="sr", n_points=500)
    print(f"SR (500 points): {sasa_sr_high.results.total_area[0]:.2f} Å²")

    # Lee-Richards
    sasa_lr = SASAAnalysis(u, select="protein")
    sasa_lr.run(algorithm="lr", n_slices=20)
    print(f"LR (20 slices):  {sasa_lr.results.total_area[0]:.2f} Å²")


def probe_radius() -> None:
    """Effect of probe radius on SASA.

    Probe radius represents the water molecule size.
    Standard value is 1.4 Å, but can be adjusted for:
    - Larger probes: More realistic solvent exclusion
    - Smaller probes: Access smaller cavities
    """
    print("\n5. Probe radius effect")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    u = mda.Universe(str(structure_file))

    for probe in [1.2, 1.4, 1.6, 2.0]:
        sasa = SASAAnalysis(u, select="protein")
        sasa.run(probe_radius=probe)
        print(f"Probe {probe:.1f} Å: {sasa.results.total_area[0]:.2f} Å²")


def trajectory_usage() -> None:
    """Usage pattern for real trajectories.

    This shows the typical usage for MD trajectory analysis.
    """
    print("\n6. Multi-frame trajectory (usage pattern)")
    print("-" * 30)

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


def chain_analysis() -> None:
    """Analyze SASA by chain/segment.

    MDAnalysis segments can be used to analyze multi-chain structures.
    """
    print("\n7. Chain/Segment analysis")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import SASAAnalysis

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    u = mda.Universe(str(structure_file))

    # Get unique segments
    segments = u.segments
    print(f"Number of segments: {len(segments)}")

    for seg in segments:
        sasa = SASAAnalysis(u, select=f"segid {seg.segid}")
        sasa.run()
        n_atoms = len(u.select_atoms(f"segid {seg.segid}"))
        print(f"Segment {seg.segid}: {n_atoms} atoms, {sasa.results.total_area[0]:.2f} Å²")


def convenience_function() -> None:
    """Use compute_sasa convenience function.

    For quick one-off calculations without creating an analysis object.
    """
    print("\n8. Using compute_sasa convenience function")
    print("-" * 30)

    import MDAnalysis as mda

    from zsasa.mdanalysis import compute_sasa

    structure_file = get_structure_file()
    if not structure_file:
        print("Example file not found")
        return

    u = mda.Universe(str(structure_file))

    # Simple usage - returns per-atom SASA
    atom_sasa = compute_sasa(u, select="protein")
    print(f"Per-atom SASA shape: {atom_sasa.shape}")
    print(f"Total SASA: {atom_sasa.sum():.2f} Å²")

    # Get per-residue SASA
    residue_sasa = compute_sasa(u, select="protein", mode="residue")
    print(f"Per-residue SASA shape: {residue_sasa.shape}")

    # Get just total SASA
    total_sasa = compute_sasa(u, select="protein", mode="total")
    print(f"Total SASA (scalar): {total_sasa[0]:.2f} Å²")


def main() -> None:
    """Run all MDAnalysis integration examples."""
    print("MDAnalysis Integration Examples")
    print("=" * 50)

    structure_file = get_structure_file()
    if not structure_file:
        print(f"Example file not found: {EXAMPLES_DIR / '1ubq.pdb'}")
        return

    print(f"Using structure file: {structure_file.name}")

    # Run each example
    basic_sasa()
    selection_based()
    per_residue_sasa()
    algorithm_options()
    probe_radius()
    trajectory_usage()
    chain_analysis()
    convenience_function()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

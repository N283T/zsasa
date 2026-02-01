"""zsasa: Fast SASA calculation using Zig.

This package provides Python bindings for the zsasa library,
a high-performance implementation of Solvent Accessible Surface Area (SASA)
calculation algorithms.

Example:
    >>> import numpy as np
    >>> from zsasa import calculate_sasa
    >>>
    >>> # Single atom
    >>> coords = np.array([[0.0, 0.0, 0.0]])
    >>> radii = np.array([1.5])
    >>> result = calculate_sasa(coords, radii)
    >>> print(f"Total SASA: {result.total_area:.2f} Å²")

    >>> # Classify atoms
    >>> from zsasa import classify_atoms, get_radius
    >>> result = classify_atoms(["ALA", "ALA"], ["CA", "O"])
    >>> print(result.radii)  # [1.87, 1.4]

Integrations:
    For structure file support, use the gemmi integration:

    >>> # pip install zsasa[gemmi]
    >>> from zsasa.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> print(f"Total: {result.total_area:.1f} Å²")

Analysis:
    For per-residue aggregation and RSA calculation:

    >>> from zsasa import aggregate_from_result
    >>> from zsasa.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> residues = aggregate_from_result(result)
    >>> for res in residues:
    ...     if res.rsa is not None:
    ...         print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: RSA={res.rsa:.1%}")

MDTraj Integration:
    For MD trajectory analysis (requires mdtraj):

    >>> # pip install mdtraj
    >>> from zsasa.mdtraj import compute_sasa
    >>> import mdtraj as md
    >>> traj = md.load('trajectory.xtc', top='topology.pdb')
    >>> sasa = compute_sasa(traj)  # Returns (n_frames, n_atoms) in nm²

MDAnalysis Integration:
    For MD trajectory analysis with MDAnalysis (requires MDAnalysis):

    >>> # pip install MDAnalysis
    >>> import MDAnalysis as mda
    >>> from zsasa.mdanalysis import SASAAnalysis
    >>> u = mda.Universe('topology.pdb', 'trajectory.xtc')
    >>> sasa = SASAAnalysis(u, select='protein')
    >>> sasa.run()
    >>> print(sasa.results.total_area)  # Returns per-frame SASA in Å²

Native XTC Reader:
    For standalone XTC reading without MDTraj/MDAnalysis dependencies:

    >>> from zsasa.xtc import XtcReader, compute_sasa_trajectory
    >>> # Low-level reader API
    >>> with XtcReader("trajectory.xtc") as reader:
    ...     for frame in reader:
    ...         print(f"Step {frame.step}")
    >>>
    >>> # High-level SASA calculation (radii from topology)
    >>> result = compute_sasa_trajectory("trajectory.xtc", radii)
    >>> print(result.total_areas)  # Per-frame SASA in Å²
"""

from zsasa.analysis import (
    ResidueResult,
    aggregate_by_residue,
    aggregate_from_result,
)
from zsasa.core import (
    MAX_SASA,
    AtomClass,
    BatchSasaResult,
    ClassificationResult,
    ClassifierType,
    SasaResult,
    calculate_rsa,
    calculate_rsa_batch,
    calculate_sasa,
    calculate_sasa_batch,
    classify_atoms,
    get_atom_class,
    get_max_sasa,
    get_radius,
    get_version,
    guess_radius,
    guess_radius_from_atom_name,
)

__all__ = [
    # SASA calculation
    "calculate_sasa",
    "calculate_sasa_batch",
    "SasaResult",
    "BatchSasaResult",
    # Classifier
    "ClassifierType",
    "AtomClass",
    "ClassificationResult",
    "get_radius",
    "get_atom_class",
    "guess_radius",
    "guess_radius_from_atom_name",
    "classify_atoms",
    # RSA
    "MAX_SASA",
    "get_max_sasa",
    "calculate_rsa",
    "calculate_rsa_batch",
    # Analysis
    "ResidueResult",
    "aggregate_by_residue",
    "aggregate_from_result",
    # Utility
    "get_version",
]

__version__ = get_version()

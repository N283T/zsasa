"""freesasa-zig: Fast SASA calculation using Zig.

This package provides Python bindings for the freesasa-zig library,
a high-performance implementation of Solvent Accessible Surface Area (SASA)
calculation algorithms.

Example:
    >>> import numpy as np
    >>> from freesasa_zig import calculate_sasa
    >>>
    >>> # Single atom
    >>> coords = np.array([[0.0, 0.0, 0.0]])
    >>> radii = np.array([1.5])
    >>> result = calculate_sasa(coords, radii)
    >>> print(f"Total SASA: {result.total_area:.2f} Å²")

    >>> # Classify atoms
    >>> from freesasa_zig import classify_atoms, get_radius
    >>> result = classify_atoms(["ALA", "ALA"], ["CA", "O"])
    >>> print(result.radii)  # [1.87, 1.4]

Integrations:
    For structure file support, use the gemmi integration:

    >>> # pip install freesasa-zig[gemmi]
    >>> from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> print(f"Total: {result.total_area:.1f} Å²")

Analysis:
    For per-residue aggregation and RSA calculation:

    >>> from freesasa_zig import aggregate_from_result
    >>> from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> residues = aggregate_from_result(result)
    >>> for res in residues:
    ...     if res.rsa is not None:
    ...         print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: RSA={res.rsa:.1%}")

MDTraj Integration:
    For MD trajectory analysis (requires mdtraj):

    >>> # pip install mdtraj
    >>> from freesasa_zig.mdtraj import compute_sasa
    >>> import mdtraj as md
    >>> traj = md.load('trajectory.xtc', top='topology.pdb')
    >>> sasa = compute_sasa(traj)  # Returns (n_frames, n_atoms) in nm²
"""

from freesasa_zig.analysis import (
    ResidueResult,
    aggregate_by_residue,
    aggregate_from_result,
)
from freesasa_zig.core import (
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

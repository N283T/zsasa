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
"""

from freesasa_zig.core import (
    MAX_SASA,
    AtomClass,
    ClassificationResult,
    ClassifierType,
    SasaResult,
    calculate_rsa,
    calculate_rsa_batch,
    calculate_sasa,
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
    "SasaResult",
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
    # Utility
    "get_version",
]

__version__ = get_version()

"""Biotite integration for freesasa-zig.

This module provides convenience functions for calculating SASA
from Biotite AtomArray/AtomArrayStack objects.

Requires: pip install freesasa-zig[biotite]

Also works with AtomWorks, which is built on Biotite:
    >>> from atomworks.io.utils.io_utils import load_any
    >>> from freesasa_zig.integrations.biotite import calculate_sasa_from_atom_array
    >>> atom_array = load_any("protein.cif")
    >>> result = calculate_sasa_from_atom_array(atom_array)

Example:
    >>> from freesasa_zig.integrations.biotite import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.pdb")
    >>> print(f"Total SASA: {result.total_area:.2f} Å²")
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np

from freesasa_zig.core import (
    AtomClass,
    ClassifierType,
    calculate_sasa,
    classify_atoms,
)
from freesasa_zig.integrations._types import AtomData, SasaResultWithAtoms

if TYPE_CHECKING:
    from biotite.structure import AtomArray, AtomArrayStack

__all__ = [
    "AtomData",
    "SasaResultWithAtoms",
    "extract_atoms_from_atom_array",
    "calculate_sasa_from_atom_array",
    "calculate_sasa_from_structure",
]

# Hydrogen-like elements (includes deuterium and tritium)
# Note: More comprehensive than gemmi's Element("H") check
_HYDROGEN_ELEMENTS = frozenset({"H", "D", "T"})


def _import_biotite():
    """Import biotite with helpful error message."""
    try:
        import biotite.structure as struc
        import biotite.structure.io as strucio

        return struc, strucio
    except ImportError as e:
        msg = (
            "biotite is required for this functionality. "
            "Install with: pip install freesasa-zig[biotite]"
        )
        raise ImportError(msg) from e


def extract_atoms_from_atom_array(
    atom_array: AtomArray,
    *,
    include_hetatm: bool = True,
    include_hydrogens: bool = False,
) -> AtomData:
    """Extract atom data from a Biotite AtomArray.

    Args:
        atom_array: A Biotite AtomArray object.
        include_hetatm: Whether to include HETATM records (ligands, waters, etc.).
        include_hydrogens: Whether to include hydrogen atoms.

    Returns:
        AtomData containing coordinates and atom metadata.

    Example:
        >>> import biotite.structure.io as strucio
        >>> atom_array = strucio.load_structure("protein.pdb")
        >>> atoms = extract_atoms_from_atom_array(atom_array)
        >>> print(f"Extracted {len(atoms)} atoms")
    """
    _import_biotite()  # Ensure biotite is available

    # Build filter mask
    mask = np.ones(len(atom_array), dtype=bool)

    # Filter HETATM if not requested
    if not include_hetatm:
        mask &= ~atom_array.hetero

    # Filter hydrogens if not requested
    if not include_hydrogens:
        mask &= ~np.isin(atom_array.element, list(_HYDROGEN_ELEMENTS))

    # Apply mask
    filtered = atom_array[mask]

    if len(filtered) == 0:
        return AtomData(
            coords=np.empty((0, 3), dtype=np.float64),
            residue_names=[],
            atom_names=[],
            chain_ids=[],
            residue_ids=[],
            elements=[],
        )

    return AtomData(
        coords=np.asarray(filtered.coord, dtype=np.float64),
        residue_names=filtered.res_name.tolist(),
        atom_names=filtered.atom_name.tolist(),
        chain_ids=filtered.chain_id.tolist(),
        residue_ids=filtered.res_id.tolist(),
        elements=filtered.element.tolist(),
    )


def calculate_sasa_from_atom_array(
    atom_array: AtomArray,
    *,
    classifier: ClassifierType = ClassifierType.NACCESS,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    include_hetatm: bool = True,
    include_hydrogens: bool = False,
) -> SasaResultWithAtoms:
    """Calculate SASA from a Biotite AtomArray.

    This is a convenience function that extracts atoms, classifies them,
    and calculates SASA in one step.

    Args:
        atom_array: A Biotite AtomArray object.
        classifier: Classifier for atom radii. Default: NACCESS.
        algorithm: SASA algorithm ("sr" or "lr"). Default: "sr".
        n_points: Test points per atom (SR algorithm). Default: 100.
        n_slices: Slices per atom (LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads (0 = auto). Default: 0.
        include_hetatm: Include HETATM records. Default: True.
        include_hydrogens: Include hydrogen atoms. Default: False.

    Returns:
        SasaResultWithAtoms with SASA values and atom metadata.

    Example:
        >>> import biotite.structure.io as strucio
        >>> atom_array = strucio.load_structure("protein.pdb")
        >>> result = calculate_sasa_from_atom_array(atom_array)
        >>> print(f"Total: {result.total_area:.1f} Å²")
        >>> print(f"Polar: {result.polar_area:.1f} Å²")
        >>> print(f"Apolar: {result.apolar_area:.1f} Å²")
    """
    # Extract atoms
    atom_data = extract_atoms_from_atom_array(
        atom_array,
        include_hetatm=include_hetatm,
        include_hydrogens=include_hydrogens,
    )

    if len(atom_data) == 0:
        return SasaResultWithAtoms(
            total_area=0.0,
            atom_areas=np.array([], dtype=np.float64),
            atom_classes=np.array([], dtype=np.int32),
            atom_data=atom_data,
            polar_area=0.0,
            apolar_area=0.0,
        )

    # Classify atoms
    classification = classify_atoms(
        atom_data.residue_names,
        atom_data.atom_names,
        classifier,
    )

    # Calculate SASA
    sasa_result = calculate_sasa(
        atom_data.coords,
        classification.radii,
        algorithm=algorithm,
        n_points=n_points,
        n_slices=n_slices,
        probe_radius=probe_radius,
        n_threads=n_threads,
    )

    # Calculate polar/apolar areas
    polar_mask = classification.classes == AtomClass.POLAR
    apolar_mask = classification.classes == AtomClass.APOLAR

    polar_area = float(np.sum(sasa_result.atom_areas[polar_mask]))
    apolar_area = float(np.sum(sasa_result.atom_areas[apolar_mask]))

    return SasaResultWithAtoms(
        total_area=sasa_result.total_area,
        atom_areas=sasa_result.atom_areas,
        atom_classes=classification.classes,
        atom_data=atom_data,
        polar_area=polar_area,
        apolar_area=apolar_area,
    )


def calculate_sasa_from_structure(
    source: str | Path | AtomArray | AtomArrayStack,
    *,
    model_index: int = 0,
    classifier: ClassifierType = ClassifierType.NACCESS,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    include_hetatm: bool = True,
    include_hydrogens: bool = False,
) -> SasaResultWithAtoms:
    """Calculate SASA from a structure file or Biotite AtomArray/AtomArrayStack.

    This is the highest-level convenience function. It accepts either
    a file path (PDB, mmCIF, BinaryCIF) or a Biotite structure object.

    Args:
        source: Path to structure file or Biotite AtomArray/AtomArrayStack.
        model_index: Model index to use (for AtomArrayStack). Default: 0.
        classifier: Classifier for atom radii. Default: NACCESS.
        algorithm: SASA algorithm ("sr" or "lr"). Default: "sr".
        n_points: Test points per atom (SR algorithm). Default: 100.
        n_slices: Slices per atom (LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads (0 = auto). Default: 0.
        include_hetatm: Include HETATM records. Default: True.
        include_hydrogens: Include hydrogen atoms. Default: False.

    Returns:
        SasaResultWithAtoms with SASA values and atom metadata.

    Example:
        >>> from freesasa_zig.integrations.biotite import calculate_sasa_from_structure
        >>>
        >>> # From file path
        >>> result = calculate_sasa_from_structure("protein.pdb")
        >>> print(f"Total: {result.total_area:.1f} Å²")
        >>>
        >>> # From Biotite AtomArray
        >>> import biotite.structure.io as strucio
        >>> atom_array = strucio.load_structure("protein.cif")
        >>> result = calculate_sasa_from_structure(atom_array)
    """
    struc, strucio = _import_biotite()

    # Load structure if path is given
    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"Structure file not found: {path}")
        # load_structure returns AtomArray for single model, AtomArrayStack for multiple
        structure = strucio.load_structure(str(path))
    else:
        structure = source

    # Handle AtomArrayStack (multiple models)
    if isinstance(structure, struc.AtomArrayStack):
        n_models = structure.stack_depth()
        if model_index < 0 or model_index >= n_models:
            msg = f"Model index {model_index} out of range (structure has {n_models} models)"
            raise IndexError(msg)
        atom_array = structure[model_index]
    else:
        # Single AtomArray
        if model_index != 0:
            msg = f"Model index {model_index} out of range (structure has 1 model)"
            raise IndexError(msg)
        atom_array = structure

    return calculate_sasa_from_atom_array(
        atom_array,
        classifier=classifier,
        algorithm=algorithm,
        n_points=n_points,
        n_slices=n_slices,
        probe_radius=probe_radius,
        n_threads=n_threads,
        include_hetatm=include_hetatm,
        include_hydrogens=include_hydrogens,
    )

"""Gemmi integration for zsasa.

This module provides convenience functions for calculating SASA
from gemmi Structure/Model objects.

Requires: pip install zsasa[gemmi]

Example:
    >>> from zsasa.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> print(f"Total SASA: {result.total_area:.2f} Å²")
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np

from zsasa.classifier import AtomClass, ClassifierType, classify_atoms
from zsasa.integrations._types import AtomData, SasaResultWithAtoms
from zsasa.sasa import calculate_sasa

if TYPE_CHECKING:
    import gemmi

__all__ = [
    "AtomData",
    "SasaResultWithAtoms",
    "extract_atoms_from_model",
    "calculate_sasa_from_model",
    "calculate_sasa_from_structure",
]


def _import_gemmi() -> "gemmi":  # noqa: UP037
    """Import gemmi with helpful error message."""
    try:
        import gemmi

        return gemmi
    except ImportError as e:
        msg = "gemmi is required for this functionality. Install with: pip install zsasa[gemmi]"
        raise ImportError(msg) from e


def extract_atoms_from_model(
    model: "gemmi.Model",  # noqa: UP037
    *,
    include_hetatm: bool = True,
    include_hydrogens: bool = False,
) -> AtomData:
    """Extract atom data from a gemmi Model.

    Args:
        model: A gemmi Model object.
        include_hetatm: Whether to include HETATM records (ligands, waters, etc.).
        include_hydrogens: Whether to include hydrogen atoms.

    Returns:
        AtomData containing coordinates and atom metadata.

    Example:
        >>> import gemmi
        >>> structure = gemmi.read_structure("protein.cif")
        >>> atoms = extract_atoms_from_model(structure[0])
        >>> print(f"Extracted {len(atoms)} atoms")
    """
    gemmi = _import_gemmi()

    coords = []
    residue_names = []
    atom_names = []
    chain_ids = []
    residue_ids = []
    elements = []

    for chain in model:
        for residue in chain:
            # Skip HETATM if not requested
            if not include_hetatm and residue.het_flag == "H":
                continue

            for atom in residue:
                # Skip hydrogens if not requested
                if not include_hydrogens and atom.element == gemmi.Element("H"):
                    continue

                coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                residue_names.append(residue.name)
                atom_names.append(atom.name)
                chain_ids.append(chain.name)
                residue_ids.append(residue.seqid.num)
                elements.append(atom.element.name)

    return AtomData(
        coords=np.array(coords, dtype=np.float64),
        residue_names=residue_names,
        atom_names=atom_names,
        chain_ids=chain_ids,
        residue_ids=residue_ids,
        elements=elements,
    )


def calculate_sasa_from_model(
    model: "gemmi.Model",  # noqa: UP037
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
    """Calculate SASA from a gemmi Model.

    This is a convenience function that extracts atoms, classifies them,
    and calculates SASA in one step.

    Args:
        model: A gemmi Model object.
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
        >>> import gemmi
        >>> structure = gemmi.read_structure("protein.cif")
        >>> result = calculate_sasa_from_model(structure[0])
        >>> print(f"Total: {result.total_area:.1f} Å²")
        >>> print(f"Polar: {result.polar_area:.1f} Å²")
        >>> print(f"Apolar: {result.apolar_area:.1f} Å²")
    """
    # Extract atoms
    atom_data = extract_atoms_from_model(
        model,
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
    source: str | Path | "gemmi.Structure",  # noqa: UP037
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
    """Calculate SASA from a structure file or gemmi Structure.

    This is the highest-level convenience function. It accepts either
    a file path (mmCIF or PDB) or a gemmi Structure object.

    Args:
        source: Path to structure file (mmCIF/PDB) or gemmi Structure object.
        model_index: Model index to use. Default: 0 (first model).
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
        >>> from zsasa.integrations.gemmi import calculate_sasa_from_structure
        >>>
        >>> # From file path
        >>> result = calculate_sasa_from_structure("protein.cif")
        >>> print(f"Total: {result.total_area:.1f} Å²")
        >>>
        >>> # From gemmi Structure
        >>> import gemmi
        >>> structure = gemmi.read_structure("protein.pdb")
        >>> result = calculate_sasa_from_structure(structure)
    """
    gemmi = _import_gemmi()

    # Load structure if path is given
    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"Structure file not found: {path}")
        structure = gemmi.read_structure(str(path))
    else:
        structure = source

    if model_index < 0 or model_index >= len(structure):
        msg = f"Model index {model_index} out of range (structure has {len(structure)} models)"
        raise IndexError(msg)

    model = structure[model_index]

    return calculate_sasa_from_model(
        model,
        classifier=classifier,
        algorithm=algorithm,
        n_points=n_points,
        n_slices=n_slices,
        probe_radius=probe_radius,
        n_threads=n_threads,
        include_hetatm=include_hetatm,
        include_hydrogens=include_hydrogens,
    )

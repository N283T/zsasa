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

from zsasa.classifier import AtomClass, ClassifierType
from zsasa.integrations._altloc import SiteAtom, keep_auto_altloc
from zsasa.integrations._types import AtomData, SasaResultWithAtoms, classify_atom_data
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
    include_hetatm: bool = False,
    include_hydrogens: bool = False,
) -> AtomData:
    """Extract atom data from a gemmi Model.

    Where the model has alternate locations, one conformer is kept by the
    rules of ``--altloc=auto`` in the zsasa command line: an atom without an
    altloc ID first, then altloc ``A``, then the highest occupancy, and one
    whole residue where the alternates of a position are different residues
    (microheterogeneity). To choose conformers yourself, edit the model before
    passing it.

    Args:
        model: A gemmi Model object.
        include_hetatm: Whether to include HETATM records (ligands, waters, etc.).
        include_hydrogens: Whether to include hydrogen atoms. Deuterium counts
            as hydrogen.

    Returns:
        AtomData containing coordinates and atom metadata.

    Example:
        >>> import gemmi
        >>> structure = gemmi.read_structure("protein.cif")
        >>> atoms = extract_atoms_from_model(structure[0])
        >>> print(f"Extracted {len(atoms)} atoms")
    """
    _import_gemmi()  # Ensure gemmi is available

    coords = []
    residue_names = []
    atom_names = []
    chain_ids = []
    residue_ids = []
    insertion_codes = []
    elements = []
    site_atoms = []

    for chain in model:
        for residue in chain:
            # Skip HETATM if not requested
            if not include_hetatm and residue.het_flag == "H":
                continue

            position = (chain.name, residue.seqid.num, residue.seqid.icode)
            for atom in residue:
                # Skip hydrogens if not requested (deuterium is an isotope of H)
                if not include_hydrogens and atom.is_hydrogen():
                    continue

                coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                residue_names.append(residue.name)
                atom_names.append(atom.name)
                chain_ids.append(chain.name)
                residue_ids.append(residue.seqid.num)
                # gemmi writes a blank for a residue without an insertion code
                insertion_codes.append(residue.seqid.icode.strip())
                elements.append(atom.element.name)
                site_atoms.append(
                    SiteAtom(
                        position=position,
                        residue_name=residue.name,
                        atom_name=atom.name,
                        altloc=atom.altloc if atom.has_altloc() else "",
                        occupancy=atom.occ,
                    )
                )

    # Iterating a gemmi residue yields the atoms of every alternate conformer
    keep = keep_auto_altloc(site_atoms)
    if not all(keep):
        coords, residue_names, atom_names, chain_ids, residue_ids, insertion_codes, elements = (
            [value for value, kept in zip(values, keep, strict=True) if kept]
            for values in (
                coords,
                residue_names,
                atom_names,
                chain_ids,
                residue_ids,
                insertion_codes,
                elements,
            )
        )

    return AtomData(
        coords=np.array(coords, dtype=np.float64),
        residue_names=residue_names,
        atom_names=atom_names,
        chain_ids=chain_ids,
        residue_ids=residue_ids,
        elements=elements,
        insertion_codes=insertion_codes,
    )


def calculate_sasa_from_model(
    model: "gemmi.Model",  # noqa: UP037
    *,
    classifier: ClassifierType = ClassifierType.CCD,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    include_hetatm: bool = False,
    include_hydrogens: bool = False,
) -> SasaResultWithAtoms:
    """Calculate SASA from a gemmi Model.

    This is a convenience function that extracts atoms, classifies them,
    and calculates SASA in one step.

    Args:
        model: A gemmi Model object.
        classifier: Classifier for atom radii. Default: CCD.
        algorithm: SASA algorithm ("sr" or "lr"). Default: "sr".
        n_points: Test points per atom (SR algorithm). Default: 100.
        n_slices: Slices per atom (LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads (0 = auto). Default: 0.
        include_hetatm: Include HETATM records. Default: False.
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

    # Classify atoms, falling back to element-derived radii when available.
    classification = classify_atom_data(atom_data, classifier)

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
    classifier: ClassifierType = ClassifierType.CCD,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    include_hetatm: bool = False,
    include_hydrogens: bool = False,
) -> SasaResultWithAtoms:
    """Calculate SASA from a structure file or gemmi Structure.

    This is the highest-level convenience function. It accepts either
    a file path (mmCIF or PDB) or a gemmi Structure object.

    Args:
        source: Path to structure file (mmCIF/PDB) or gemmi Structure object.
        model_index: Model index to use. Default: 0 (first model).
        classifier: Classifier for atom radii. Default: CCD.
        algorithm: SASA algorithm ("sr" or "lr"). Default: "sr".
        n_points: Test points per atom (SR algorithm). Default: 100.
        n_slices: Slices per atom (LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads (0 = auto). Default: 0.
        include_hetatm: Include HETATM records. Default: False.
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

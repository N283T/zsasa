"""BioPython integration for freesasa-zig.

This module provides convenience functions for calculating SASA
from BioPython Structure/Model objects.

Requires: pip install freesasa-zig[biopython]

Example:
    >>> from freesasa_zig.integrations.biopython import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.pdb")
    >>> print(f"Total SASA: {result.total_area:.2f} A^2")
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
from freesasa_zig.integrations.gemmi import AtomData, SasaResultWithAtoms

if TYPE_CHECKING:
    from Bio.PDB.Atom import Atom
    from Bio.PDB.Model import Model
    from Bio.PDB.Structure import Structure

__all__ = [
    "extract_atoms_from_model",
    "calculate_sasa_from_model",
    "calculate_sasa_from_structure",
]


def _import_biopython():
    """Import BioPython with helpful error message."""
    try:
        from Bio.PDB import PDBParser

        return PDBParser
    except ImportError as e:
        msg = (
            "BioPython is required for this functionality. "
            "Install with: pip install freesasa-zig[biopython]"
        )
        raise ImportError(msg) from e


def _is_hydrogen(atom: Atom) -> bool:
    """Check if an atom is hydrogen."""
    element = getattr(atom, "element", None)
    if element:
        return element.strip().upper() == "H"
    # Fallback: check atom name
    name = atom.get_name().strip()
    return name.startswith("H") or (len(name) > 1 and name[0].isdigit() and name[1] == "H")


def _is_hetatm(residue_id: tuple) -> bool:
    """Check if a residue is a HETATM (hetero residue).

    BioPython residue ID format: (hetfield, resseq, icode)
    - hetfield is ' ' for standard residues
    - hetfield is 'H_xxx' for hetero residues
    - hetfield is 'W' for water
    """
    hetfield = residue_id[0]
    return hetfield != " "


def extract_atoms_from_model(
    model: Model,
    *,
    include_hetatm: bool = True,
    include_hydrogens: bool = False,
) -> AtomData:
    """Extract atom data from a BioPython Model.

    Args:
        model: A BioPython Model object.
        include_hetatm: Whether to include HETATM records (ligands, waters, etc.).
        include_hydrogens: Whether to include hydrogen atoms.

    Returns:
        AtomData containing coordinates and atom metadata.

    Example:
        >>> from Bio.PDB import PDBParser
        >>> parser = PDBParser(QUIET=True)
        >>> structure = parser.get_structure("protein", "protein.pdb")
        >>> atoms = extract_atoms_from_model(structure[0])
        >>> print(f"Extracted {len(atoms)} atoms")
    """
    _import_biopython()  # Ensure BioPython is available

    coords = []
    residue_names = []
    atom_names = []
    chain_ids = []
    residue_ids = []
    elements = []

    for chain in model:
        chain_id = chain.get_id()

        for residue in chain:
            res_id = residue.get_id()

            # Skip HETATM if not requested
            if not include_hetatm and _is_hetatm(res_id):
                continue

            res_name = residue.get_resname()
            res_seq = res_id[1]  # (hetfield, resseq, icode)

            for atom in residue:
                # Skip hydrogens if not requested
                if not include_hydrogens and _is_hydrogen(atom):
                    continue

                coord = atom.get_coord()
                coords.append([coord[0], coord[1], coord[2]])
                residue_names.append(res_name)
                atom_names.append(atom.get_name())
                chain_ids.append(chain_id)
                residue_ids.append(res_seq)
                elements.append(getattr(atom, "element", "") or "")

    return AtomData(
        coords=np.array(coords, dtype=np.float64) if coords else np.empty((0, 3), dtype=np.float64),
        residue_names=residue_names,
        atom_names=atom_names,
        chain_ids=chain_ids,
        residue_ids=residue_ids,
        elements=elements,
    )


def calculate_sasa_from_model(
    model: Model,
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
    """Calculate SASA from a BioPython Model.

    This is a convenience function that extracts atoms, classifies them,
    and calculates SASA in one step.

    Args:
        model: A BioPython Model object.
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
        >>> from Bio.PDB import PDBParser
        >>> parser = PDBParser(QUIET=True)
        >>> structure = parser.get_structure("protein", "protein.pdb")
        >>> result = calculate_sasa_from_model(structure[0])
        >>> print(f"Total: {result.total_area:.1f} A^2")
        >>> print(f"Polar: {result.polar_area:.1f} A^2")
        >>> print(f"Apolar: {result.apolar_area:.1f} A^2")
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
    source: str | Path | Structure,
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
    """Calculate SASA from a structure file or BioPython Structure.

    This is the highest-level convenience function. It accepts either
    a file path (PDB or mmCIF) or a BioPython Structure object.

    Args:
        source: Path to structure file (PDB/mmCIF) or BioPython Structure object.
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
        >>> from freesasa_zig.integrations.biopython import calculate_sasa_from_structure
        >>>
        >>> # From file path
        >>> result = calculate_sasa_from_structure("protein.pdb")
        >>> print(f"Total: {result.total_area:.1f} A^2")
        >>>
        >>> # From BioPython Structure
        >>> from Bio.PDB import PDBParser
        >>> parser = PDBParser(QUIET=True)
        >>> structure = parser.get_structure("protein", "protein.pdb")
        >>> result = calculate_sasa_from_structure(structure)
    """
    _import_biopython()
    from Bio.PDB import MMCIFParser, PDBParser

    # Load structure if path is given
    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"Structure file not found: {path}")

        # Choose parser based on file extension
        suffix = path.suffix.lower()
        parser = MMCIFParser(QUIET=True) if suffix in (".cif", ".mmcif") else PDBParser(QUIET=True)

        structure = parser.get_structure(path.stem, str(path))
    else:
        structure = source

    # Get list of models
    models = list(structure.get_models())
    if model_index < 0 or model_index >= len(models):
        msg = f"Model index {model_index} out of range (structure has {len(models)} models)"
        raise IndexError(msg)

    model = models[model_index]

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

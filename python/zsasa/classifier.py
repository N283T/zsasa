"""Atom classifier functions for radius and polarity assignment."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np
from numpy.typing import NDArray

from zsasa._ffi import (
    ZSASA_ATOM_CLASS_APOLAR,
    ZSASA_ATOM_CLASS_POLAR,
    ZSASA_ATOM_CLASS_UNKNOWN,
    ZSASA_CLASSIFIER_CCD,
    ZSASA_CLASSIFIER_NACCESS,
    ZSASA_CLASSIFIER_OONS,
    ZSASA_CLASSIFIER_PROTOR,
    ZSASA_ERROR_INVALID_INPUT,
    ZSASA_OK,
    _get_lib,
)


class ClassifierType(IntEnum):
    """Available classifier types for atom radius assignment.

    Attributes:
        NACCESS: NACCESS-compatible radii (default, most commonly used).
        PROTOR: ProtOr radii based on hybridization state.
        OONS: Ooi, Oobatake, Nemethy, Scheraga radii (older FreeSASA default).
        CCD: CCD-based radii derived from bond topology (ProtOr-compatible).
    """

    NACCESS = ZSASA_CLASSIFIER_NACCESS
    PROTOR = ZSASA_CLASSIFIER_PROTOR
    OONS = ZSASA_CLASSIFIER_OONS
    CCD = ZSASA_CLASSIFIER_CCD


class AtomClass(IntEnum):
    """Atom polarity classes.

    Attributes:
        POLAR: Polar atoms (e.g., N, O).
        APOLAR: Apolar/hydrophobic atoms (e.g., C).
        UNKNOWN: Unknown classification.
    """

    POLAR = ZSASA_ATOM_CLASS_POLAR
    APOLAR = ZSASA_ATOM_CLASS_APOLAR
    UNKNOWN = ZSASA_ATOM_CLASS_UNKNOWN


def get_radius(
    residue: str,
    atom: str,
    classifier_type: ClassifierType = ClassifierType.NACCESS,
) -> float | None:
    """Get van der Waals radius for an atom using the specified classifier.

    Args:
        residue: Residue name (e.g., "ALA", "GLY").
        atom: Atom name (e.g., "CA", "CB").
        classifier_type: Classifier to use. Default: ClassifierType.NACCESS.

    Returns:
        Radius in Angstroms, or None if atom is not found in classifier.

    Example:
        >>> from zsasa import get_radius, ClassifierType
        >>> get_radius("ALA", "CA")
        1.87
        >>> get_radius("ALA", "XX")  # Unknown atom
        None
    """
    _, lib = _get_lib()
    radius = lib.zsasa_classifier_get_radius(
        classifier_type,
        residue.encode("utf-8"),
        atom.encode("utf-8"),
    )
    if np.isnan(radius):
        return None
    return radius


def get_atom_class(
    residue: str,
    atom: str,
    classifier_type: ClassifierType = ClassifierType.NACCESS,
) -> int:
    """Get atom polarity class using the specified classifier.

    Args:
        residue: Residue name (e.g., "ALA", "GLY").
        atom: Atom name (e.g., "CA", "CB").
        classifier_type: Classifier to use. Default: ClassifierType.NACCESS.

    Returns:
        AtomClass constant (POLAR, APOLAR, or UNKNOWN).

    Example:
        >>> from zsasa import get_atom_class, AtomClass
        >>> get_atom_class("ALA", "CA") == AtomClass.APOLAR
        True
        >>> get_atom_class("ALA", "O") == AtomClass.POLAR
        True
    """
    _, lib = _get_lib()
    return lib.zsasa_classifier_get_class(
        classifier_type,
        residue.encode("utf-8"),
        atom.encode("utf-8"),
    )


def guess_radius(element: str) -> float | None:
    """Guess van der Waals radius from element symbol.

    Args:
        element: Element symbol (e.g., "C", "N", "FE").
                 Case-insensitive, whitespace is trimmed.

    Returns:
        Radius in Angstroms, or None if element is not recognized.

    Example:
        >>> from zsasa import guess_radius
        >>> guess_radius("C")
        1.7
        >>> guess_radius("FE")
        1.26
        >>> guess_radius("XX")  # Unknown
        None
    """
    _, lib = _get_lib()
    radius = lib.zsasa_guess_radius(element.encode("utf-8"))
    if np.isnan(radius):
        return None
    return radius


def guess_radius_from_atom_name(atom_name: str) -> float | None:
    """Guess van der Waals radius from PDB-style atom name.

    Extracts element symbol from atom name following PDB conventions:
    - Leading space indicates single-char element (e.g., " CA " = Carbon alpha)
    - No leading space may indicate 2-char element (e.g., "FE  " = Iron)

    Args:
        atom_name: PDB-style atom name (e.g., " CA ", "FE  ").

    Returns:
        Radius in Angstroms, or None if element cannot be determined.

    Example:
        >>> from zsasa import guess_radius_from_atom_name
        >>> guess_radius_from_atom_name(" CA ")  # Carbon alpha
        1.7
        >>> guess_radius_from_atom_name("FE  ")  # Iron
        1.26
    """
    _, lib = _get_lib()
    radius = lib.zsasa_guess_radius_from_atom_name(atom_name.encode("utf-8"))
    if np.isnan(radius):
        return None
    return radius


@dataclass
class ClassificationResult:
    """Result of batch atom classification.

    Attributes:
        radii: Per-atom van der Waals radii in Angstroms.
               NaN values indicate atoms not found in the classifier.
               Use np.isnan(result.radii) to find unknown atoms.
        classes: Per-atom polarity classes (AtomClass constants).

    Example:
        >>> result = classify_atoms(residues, atoms)
        >>> unknown_mask = np.isnan(result.radii)
        >>> if unknown_mask.any():
        ...     print(f"Found {unknown_mask.sum()} unknown atoms")
    """

    radii: NDArray[np.float64]
    classes: NDArray[np.int32]

    def __repr__(self) -> str:
        return f"ClassificationResult(n_atoms={len(self.radii)})"


def classify_atoms(
    residues: list[str],
    atoms: list[str],
    classifier_type: ClassifierType = ClassifierType.NACCESS,
    *,
    include_classes: bool = True,
) -> ClassificationResult:
    """Classify multiple atoms at once (batch operation).

    This is more efficient than calling get_radius for each atom individually.

    Args:
        residues: List of residue names.
        atoms: List of atom names (must be same length as residues).
        classifier_type: Classifier to use. Default: ClassifierType.NACCESS.
        include_classes: Whether to compute atom classes. Default: True.

    Returns:
        ClassificationResult with radii and classes arrays.
        Unknown atoms have NaN radius and UNKNOWN class.

    Raises:
        ValueError: If residues and atoms have different lengths.

    Example:
        >>> from zsasa import classify_atoms
        >>> result = classify_atoms(
        ...     ["ALA", "ALA", "GLY"],
        ...     ["CA", "O", "N"],
        ... )
        >>> result.radii
        array([1.87, 1.4 , 1.65])
    """
    ffi, lib = _get_lib()

    if len(residues) != len(atoms):
        msg = f"residues and atoms must have same length: {len(residues)} != {len(atoms)}"
        raise ValueError(msg)

    n_atoms = len(residues)
    if n_atoms == 0:
        return ClassificationResult(
            radii=np.array([], dtype=np.float64),
            classes=np.array([], dtype=np.int32),
        )

    # Encode strings and create cffi arrays
    # Keep references to prevent garbage collection
    residues_bytes = [ffi.new("char[]", r.encode("utf-8")) for r in residues]
    atoms_bytes = [ffi.new("char[]", a.encode("utf-8")) for a in atoms]

    residues_arr = ffi.new("char*[]", residues_bytes)
    atoms_arr = ffi.new("char*[]", atoms_bytes)

    # Allocate output arrays
    radii = np.zeros(n_atoms, dtype=np.float64)
    classes = np.zeros(n_atoms, dtype=np.int32) if include_classes else None

    radii_ptr = ffi.cast("double*", radii.ctypes.data)
    classes_ptr = ffi.cast("int*", classes.ctypes.data) if include_classes else ffi.NULL

    result = lib.zsasa_classify_atoms(
        classifier_type,
        residues_arr,
        atoms_arr,
        n_atoms,
        radii_ptr,
        classes_ptr,
    )

    if result == ZSASA_ERROR_INVALID_INPUT:
        msg = f"Invalid classifier type: {classifier_type}"
        raise ValueError(msg)
    elif result != ZSASA_OK:
        msg = f"Classification error: {result}"
        raise RuntimeError(msg)

    if not include_classes:
        classes = np.full(n_atoms, AtomClass.UNKNOWN, dtype=np.int32)

    return ClassificationResult(radii=radii, classes=classes)

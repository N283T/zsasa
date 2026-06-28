"""Shared types for integration modules.

This module contains data classes shared across different structure parsing
library integrations (BioPython, gemmi, etc.).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from zsasa.classifier import ClassificationResult, ClassifierType, classify_atoms, guess_radius
from zsasa.sasa import SasaResult

__all__ = [
    "AtomData",
    "SasaResultWithAtoms",
    "classify_atom_data",
]


@dataclass
class AtomData:
    """Extracted atom data from a structure.

    Attributes:
        coords: Atom coordinates as (N, 3) array in Angstroms.
        residue_names: Residue names for each atom.
        atom_names: Atom names for each atom.
        chain_ids: Chain IDs for each atom.
        residue_ids: Residue sequence numbers for each atom.
        elements: Element symbols for each atom.
    """

    coords: NDArray[np.float64]
    residue_names: list[str]
    atom_names: list[str]
    chain_ids: list[str]
    residue_ids: list[int]
    elements: list[str]

    def __len__(self) -> int:
        return len(self.residue_names)

    def __repr__(self) -> str:
        return f"AtomData(n_atoms={len(self)})"


@dataclass
class SasaResultWithAtoms(SasaResult):
    """SASA result with atom metadata.

    Extends SasaResult with additional atom information for analysis.

    Attributes:
        total_area: Total solvent accessible surface area in Angstroms squared.
        atom_areas: Per-atom SASA values in Angstroms squared.
        atom_classes: Per-atom polarity classes.
        atom_data: Original atom data from the structure.
        polar_area: Total polar SASA in Angstroms squared.
        apolar_area: Total apolar SASA in Angstroms squared.

    Note:
        polar_area + apolar_area may be less than total_area if some atoms
        have UNKNOWN classification.
    """

    atom_classes: NDArray[np.int32]
    atom_data: AtomData
    polar_area: float
    apolar_area: float

    def __repr__(self) -> str:
        return (
            f"SasaResultWithAtoms(total={self.total_area:.1f}, "
            f"polar={self.polar_area:.1f}, apolar={self.apolar_area:.1f}, "
            f"n_atoms={len(self.atom_areas)})"
        )


def _format_atom_identifier(atom_data: AtomData, index: int) -> str:
    """Format a clear atom identifier for integration error messages."""
    return (
        f"chain {atom_data.chain_ids[index]} residue {atom_data.residue_ids[index]} "
        f"{atom_data.residue_names[index]} atom {atom_data.atom_names[index]}"
    )


def classify_atom_data(
    atom_data: AtomData,
    classifier: ClassifierType = ClassifierType.CCD,
) -> ClassificationResult:
    """Classify integration atoms and resolve unknown radii safely.

    Classifier misses are filled from element-derived radii when an element is
    available. If neither the classifier nor element fallback can assign a
    radius, a ValueError identifies the problematic atom instead of allowing a
    NaN radius to reach the native calculator.
    """
    classification = classify_atoms(
        atom_data.residue_names,
        atom_data.atom_names,
        classifier,
    )

    missing = np.nonzero(~np.isfinite(classification.radii))[0]
    for index in missing:
        element = atom_data.elements[index].strip()
        radius = guess_radius(element) if element else None
        if radius is None:
            identifier = _format_atom_identifier(atom_data, int(index))
            msg = f"Unknown radius for {identifier} (element={element!r})"
            raise ValueError(msg)
        classification.radii[index] = radius

    return classification

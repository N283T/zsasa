"""Shared types for integration modules.

This module contains data classes shared across different structure parsing
library integrations (BioPython, gemmi, etc.).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from zsasa.sasa import SasaResult

__all__ = [
    "AtomData",
    "SasaResultWithAtoms",
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

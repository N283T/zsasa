"""Analysis functions for per-residue SASA aggregation.

This module provides functions to aggregate per-atom SASA values to
per-residue values and calculate RSA (Relative Solvent Accessibility).

Example:
    >>> from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
    >>> from freesasa_zig.analysis import aggregate_from_result
    >>>
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> residues = aggregate_from_result(result)
    >>>
    >>> for res in residues:
    ...     rsa_str = f"{res.rsa:.1%}" if res.rsa else "N/A"
    ...     print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: "
    ...           f"{res.total_area:.1f} A^2 (RSA: {rsa_str})")
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from freesasa_zig.core import MAX_SASA, AtomClass

if TYPE_CHECKING:
    from freesasa_zig.integrations.gemmi import SasaResultWithAtoms

__all__ = [
    "ResidueResult",
    "aggregate_by_residue",
    "aggregate_from_result",
]


@dataclass
class ResidueResult:
    """Per-residue SASA result.

    Attributes:
        chain_id: Chain identifier (e.g., "A", "B").
        residue_id: Residue sequence number.
        residue_name: 3-letter residue name (e.g., "ALA", "GLY").
        total_area: Total SASA for this residue in A^2.
        polar_area: Polar SASA (N, O atoms, etc.) in A^2.
        apolar_area: Apolar SASA (C atoms, etc.) in A^2.
        rsa: Relative Solvent Accessibility (0.0-1.0+), or None for
             non-standard amino acids.
        n_atoms: Number of atoms in this residue.
    """

    chain_id: str
    residue_id: int
    residue_name: str
    total_area: float
    polar_area: float
    apolar_area: float
    rsa: float | None
    n_atoms: int

    def __repr__(self) -> str:
        rsa_str = f"{self.rsa:.3f}" if self.rsa is not None else "None"
        return (
            f"ResidueResult({self.chain_id}:{self.residue_name}{self.residue_id}, "
            f"total={self.total_area:.1f}, rsa={rsa_str}, n_atoms={self.n_atoms})"
        )


def aggregate_by_residue(
    atom_areas: NDArray[np.float64],
    chain_ids: list[str],
    residue_ids: list[int],
    residue_names: list[str],
    atom_classes: NDArray[np.int32] | None = None,
) -> list[ResidueResult]:
    """Aggregate per-atom SASA values to per-residue.

    Groups atoms by (chain_id, residue_id) and sums their SASA values.
    Also calculates polar/apolar breakdown and RSA if atom classes are provided.

    Args:
        atom_areas: Per-atom SASA values in A^2.
        chain_ids: Chain ID for each atom.
        residue_ids: Residue sequence number for each atom.
        residue_names: Residue name for each atom.
        atom_classes: Optional per-atom polarity classes (AtomClass values).
                      If provided, polar_area and apolar_area will be calculated.

    Returns:
        List of ResidueResult objects, one per unique residue.
        Results are ordered by appearance in the input (preserves chain/residue order).

    Example:
        >>> import numpy as np
        >>> from freesasa_zig.analysis import aggregate_by_residue
        >>>
        >>> atom_areas = np.array([10.0, 20.0, 15.0, 25.0])
        >>> chain_ids = ["A", "A", "A", "A"]
        >>> residue_ids = [1, 1, 2, 2]
        >>> residue_names = ["ALA", "ALA", "GLY", "GLY"]
        >>>
        >>> residues = aggregate_by_residue(
        ...     atom_areas, chain_ids, residue_ids, residue_names
        ... )
        >>> len(residues)
        2
        >>> residues[0].total_area
        30.0
    """
    n_atoms = len(atom_areas)
    if n_atoms == 0:
        return []

    # Validate input lengths
    if len(chain_ids) != n_atoms:
        msg = f"chain_ids length ({len(chain_ids)}) != atom_areas length ({n_atoms})"
        raise ValueError(msg)
    if len(residue_ids) != n_atoms:
        msg = f"residue_ids length ({len(residue_ids)}) != atom_areas length ({n_atoms})"
        raise ValueError(msg)
    if len(residue_names) != n_atoms:
        msg = f"residue_names length ({len(residue_names)}) != atom_areas length ({n_atoms})"
        raise ValueError(msg)
    if atom_classes is not None and len(atom_classes) != n_atoms:
        msg = f"atom_classes length ({len(atom_classes)}) != atom_areas length ({n_atoms})"
        raise ValueError(msg)

    # Group atoms by (chain_id, residue_id)
    # Use dict to preserve insertion order (Python 3.7+)
    residue_data: dict[tuple[str, int], dict] = {}

    for i in range(n_atoms):
        key = (chain_ids[i], residue_ids[i])

        if key not in residue_data:
            residue_data[key] = {
                "residue_name": residue_names[i],
                "total_area": 0.0,
                "polar_area": 0.0,
                "apolar_area": 0.0,
                "n_atoms": 0,
            }

        data = residue_data[key]
        area = float(atom_areas[i])
        data["total_area"] += area
        data["n_atoms"] += 1

        if atom_classes is not None:
            atom_class = atom_classes[i]
            if atom_class == AtomClass.POLAR:
                data["polar_area"] += area
            elif atom_class == AtomClass.APOLAR:
                data["apolar_area"] += area

    # Build result list
    results = []
    for (chain_id, residue_id), data in residue_data.items():
        residue_name = data["residue_name"]
        total_area = data["total_area"]

        # Calculate RSA if this is a standard amino acid
        max_sasa = MAX_SASA.get(residue_name)
        rsa = total_area / max_sasa if max_sasa is not None else None

        results.append(
            ResidueResult(
                chain_id=chain_id,
                residue_id=residue_id,
                residue_name=residue_name,
                total_area=total_area,
                polar_area=data["polar_area"],
                apolar_area=data["apolar_area"],
                rsa=rsa,
                n_atoms=data["n_atoms"],
            )
        )

    return results


def aggregate_from_result(result: SasaResultWithAtoms) -> list[ResidueResult]:
    """Aggregate per-atom SASA from a SasaResultWithAtoms to per-residue.

    This is a convenience wrapper that extracts the necessary data from
    a SasaResultWithAtoms object returned by gemmi integration functions.

    Args:
        result: A SasaResultWithAtoms object from calculate_sasa_from_structure
                or calculate_sasa_from_model.

    Returns:
        List of ResidueResult objects, one per unique residue.

    Example:
        >>> from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
        >>> from freesasa_zig.analysis import aggregate_from_result
        >>>
        >>> result = calculate_sasa_from_structure("protein.cif")
        >>> residues = aggregate_from_result(result)
        >>>
        >>> # Print buried residues (RSA < 0.25)
        >>> for res in residues:
        ...     if res.rsa is not None and res.rsa < 0.25:
        ...         print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: {res.rsa:.1%}")
    """
    return aggregate_by_residue(
        atom_areas=result.atom_areas,
        chain_ids=result.atom_data.chain_ids,
        residue_ids=result.atom_data.residue_ids,
        residue_names=result.atom_data.residue_names,
        atom_classes=result.atom_classes,
    )

"""MDTraj integration for freesasa-zig.

This module provides functions to calculate SASA for MDTraj trajectories,
offering a drop-in replacement for mdtraj.shrake_rupley() with better performance.

Example:
    >>> import mdtraj as md
    >>> from freesasa_zig.mdtraj import compute_sasa
    >>>
    >>> traj = md.load('trajectory.xtc', top='topology.pdb')
    >>> sasa = compute_sasa(traj)  # Returns (n_frames, n_atoms) array
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from freesasa_zig.core import calculate_sasa_batch

if TYPE_CHECKING:
    import mdtraj as md

# Van der Waals radii from MDTraj (in nm, converted to Angstrom)
# These match MDTraj's _ATOMIC_RADII for consistency
_ATOMIC_RADII_NM = {
    "H": 0.120,
    "He": 0.140,
    "Li": 0.076,
    "Be": 0.059,
    "B": 0.192,
    "C": 0.170,
    "N": 0.155,
    "O": 0.152,
    "F": 0.147,
    "Ne": 0.154,
    "Na": 0.102,
    "Mg": 0.086,
    "Al": 0.184,
    "Si": 0.210,
    "P": 0.180,
    "S": 0.180,
    "Cl": 0.181,
    "Ar": 0.188,
    "K": 0.138,
    "Ca": 0.114,
    "Sc": 0.211,
    "Ti": 0.200,
    "V": 0.200,
    "Cr": 0.200,
    "Mn": 0.200,
    "Fe": 0.200,
    "Co": 0.200,
    "Ni": 0.163,
    "Cu": 0.140,
    "Zn": 0.139,
    "Ga": 0.187,
    "Ge": 0.211,
    "As": 0.185,
    "Se": 0.190,
    "Br": 0.185,
    "Kr": 0.202,
    "Rb": 0.303,
    "Sr": 0.249,
    "Y": 0.200,
    "Zr": 0.200,
    "Nb": 0.200,
    "Mo": 0.200,
    "Tc": 0.200,
    "Ru": 0.200,
    "Rh": 0.200,
    "Pd": 0.163,
    "Ag": 0.172,
    "Cd": 0.158,
    "In": 0.193,
    "Sn": 0.217,
    "Sb": 0.206,
    "Te": 0.206,
    "I": 0.198,
    "Xe": 0.216,
    "Cs": 0.167,
    "Ba": 0.149,
}


def _get_radii_from_topology(topology: md.Topology) -> NDArray[np.float32]:
    """Extract atomic radii from MDTraj topology.

    Args:
        topology: MDTraj Topology object.

    Returns:
        Array of atomic radii in Angstroms.
    """
    radii = []
    for atom in topology.atoms:
        symbol = atom.element.symbol if atom.element is not None else "C"
        radius_nm = _ATOMIC_RADII_NM.get(symbol, 0.200)  # Default 2.0 Å
        radii.append(radius_nm * 10.0)  # nm -> Angstrom
    return np.array(radii, dtype=np.float32)


def compute_sasa(
    traj: md.Trajectory,
    *,
    probe_radius: float = 1.4,
    n_points: int = 960,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    mode: Literal["atom", "residue"] = "atom",
) -> NDArray[np.float32]:
    """Compute SASA for MDTraj trajectory using freesasa-zig.

    This function is a drop-in replacement for mdtraj.shrake_rupley() with
    better performance through SIMD optimization and frame-level parallelization.

    Args:
        traj: MDTraj Trajectory object.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_points: Number of test points per atom (for SR algorithm). Default: 960.
                  Higher values give more accuracy but slower computation.
        algorithm: Algorithm to use: "sr" (Shrake-Rupley) or "lr" (Lee-Richards).
                   Default: "sr".
        n_slices: Number of slices per atom (for LR algorithm). Default: 20.
        n_threads: Number of threads. 0 = auto-detect. Default: 0.
        mode: Output mode:
              - "atom": Return per-atom SASA, shape (n_frames, n_atoms).
              - "residue": Return per-residue SASA, shape (n_frames, n_residues).
              Default: "atom".

    Returns:
        SASA values in nm² (matching MDTraj's output units).
        Shape: (n_frames, n_atoms) if mode='atom', (n_frames, n_residues) if mode='residue'.

    Note:
        - Input coordinates are automatically converted from nm (MDTraj) to Angstrom.
        - Output SASA values are converted from Å² back to nm².

    Example:
        >>> import mdtraj as md
        >>> from freesasa_zig.mdtraj import compute_sasa
        >>>
        >>> traj = md.load('trajectory.xtc', top='topology.pdb')
        >>>
        >>> # Per-atom SASA
        >>> sasa_atom = compute_sasa(traj, mode='atom')
        >>> print(sasa_atom.shape)  # (n_frames, n_atoms)
        >>>
        >>> # Per-residue SASA
        >>> sasa_residue = compute_sasa(traj, mode='residue')
        >>> print(sasa_residue.shape)  # (n_frames, n_residues)
    """
    # Convert coordinates from nm to Angstrom
    coords_angstrom = np.ascontiguousarray(traj.xyz * 10.0, dtype=np.float32)

    # Get radii from topology (in Angstrom)
    radii = _get_radii_from_topology(traj.topology)

    # Calculate SASA using batch API
    if algorithm == "sr":
        result = calculate_sasa_batch(
            coords_angstrom,
            radii,
            algorithm="sr",
            n_points=n_points,
            probe_radius=probe_radius,
            n_threads=n_threads,
        )
    elif algorithm == "lr":
        result = calculate_sasa_batch(
            coords_angstrom,
            radii,
            algorithm="lr",
            n_slices=n_slices,
            probe_radius=probe_radius,
            n_threads=n_threads,
        )
    else:
        msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
        raise ValueError(msg)

    # Convert from Å² to nm² (1 nm² = 100 Å²)
    sasa_nm2 = result.atom_areas / 100.0

    # Aggregate by residue if requested
    if mode == "residue":
        n_frames = sasa_nm2.shape[0]
        n_residues = traj.n_residues

        # Create residue mapping
        atom_to_residue = np.array(
            [atom.residue.index for atom in traj.topology.atoms],
            dtype=np.int32,
        )

        # Aggregate using numpy bincount
        sasa_residue = np.zeros((n_frames, n_residues), dtype=np.float32)
        for frame_idx in range(n_frames):
            sasa_residue[frame_idx] = np.bincount(
                atom_to_residue,
                weights=sasa_nm2[frame_idx],
                minlength=n_residues,
            )

        return sasa_residue

    return sasa_nm2


# Alias for compatibility
shrake_rupley = compute_sasa

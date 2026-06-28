"""MDAnalysis integration for zsasa.

This module provides a SASAAnalysis class that integrates with MDAnalysis,
offering a high-performance alternative to mdakit-sasa.

Example:
    >>> import MDAnalysis as mda
    >>> from zsasa.mdanalysis import SASAAnalysis
    >>>
    >>> u = mda.Universe("topology.pdb", "trajectory.xtc")
    >>> sasa = SASAAnalysis(u, select="protein")
    >>> sasa.run()
    >>>
    >>> print(f"Mean SASA: {sasa.results.mean_total_area:.2f} Å²")
    >>> print(f"Per-frame: {sasa.results.total_area}")
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from zsasa.sasa import calculate_sasa_batch

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe

logger = logging.getLogger(__name__)

# Import MDAnalysis vdwradii table for consistency
# Keys are uppercase element symbols (e.g., "C", "N", "CA" for calcium)
try:
    from MDAnalysis.guesser.tables import vdwradii as _MDA_VDW_RADII
except ImportError:
    # Fallback if MDAnalysis structure changes
    _MDA_VDW_RADII: dict[str, float] = {
        "H": 1.1,
        "C": 1.7,
        "N": 1.55,
        "O": 1.52,
        "S": 1.8,
        "P": 1.8,
    }

# Default radius for unknown elements (conservative estimate, larger than most common elements)
_DEFAULT_RADIUS = 2.0


def _get_element(atom) -> str:  # noqa: ANN001
    """Get element symbol from MDAnalysis atom.

    Tries multiple methods with fallbacks.
    """
    # Try element attribute first
    try:
        if hasattr(atom, "element") and atom.element:
            return str(atom.element).capitalize()
    except (AttributeError, TypeError):
        pass

    # Try type attribute (first character)
    try:
        if hasattr(atom, "type") and atom.type:
            return str(atom.type)[0].upper()
    except (AttributeError, TypeError, IndexError):
        pass

    # Try name attribute (first character)
    try:
        if hasattr(atom, "name") and atom.name:
            # Handle names like "CA", "CB" -> "C"
            name = str(atom.name).strip()
            if name:
                return name[0].upper()
    except (AttributeError, TypeError, IndexError):
        pass

    return "C"  # Default to carbon


def _get_radius(atom) -> float:  # noqa: ANN001
    """Get van der Waals radius for an atom in Angstrom.

    Uses MDAnalysis vdwradii table for consistency with MDAnalysis ecosystem.
    """
    element = _get_element(atom)
    # MDAnalysis uses uppercase keys (e.g., "C", "N", "CA" for calcium)
    return _MDA_VDW_RADII.get(element.upper(), _DEFAULT_RADIUS)


def _get_radii_from_atomgroup(atomgroup: AtomGroup) -> NDArray[np.float32]:
    """Extract atomic radii from MDAnalysis AtomGroup.

    Args:
        atomgroup: MDAnalysis AtomGroup object.

    Returns:
        Array of atomic radii in Angstroms.
    """
    radii = np.array([_get_radius(atom) for atom in atomgroup], dtype=np.float32)
    return radii


class SASAAnalysis:
    """Solvent Accessible Surface Area analysis for MDAnalysis.

    This class computes SASA for MDAnalysis Universe or AtomGroup objects,
    providing a high-performance alternative to mdakit-sasa.

    Parameters
    ----------
    universe_or_atomgroup : Universe or AtomGroup
        MDAnalysis Universe or AtomGroup to analyze.
    select : str, optional
        Atom selection string (default: "all").

    Attributes
    ----------
    atomgroup : AtomGroup
        The atoms being analyzed.
    results : Results
        Results object containing SASA data after run().

    Example
    -------
    >>> import MDAnalysis as mda
    >>> from zsasa.mdanalysis import SASAAnalysis
    >>>
    >>> u = mda.Universe("protein.pdb", "trajectory.xtc")
    >>> sasa = SASAAnalysis(u, select="protein")
    >>> sasa.run(start=0, stop=100, step=10)
    >>>
    >>> print(sasa.results.total_area)       # Per-frame total SASA
    >>> print(sasa.results.residue_area)     # Per-residue SASA
    >>> print(sasa.results.atom_area)        # Per-atom SASA
    >>> print(sasa.results.mean_total_area)  # Mean total SASA
    """

    def __init__(
        self,
        universe_or_atomgroup: Universe | AtomGroup,
        select: str = "all",
    ) -> None:
        """Initialize SASAAnalysis."""
        # Get universe and atomgroup
        if hasattr(universe_or_atomgroup, "universe"):
            self.universe = universe_or_atomgroup.universe
            if hasattr(universe_or_atomgroup, "select_atoms"):
                self.atomgroup = universe_or_atomgroup.select_atoms(select)
            else:
                # Already an AtomGroup
                self.atomgroup = universe_or_atomgroup
        else:
            self.universe = universe_or_atomgroup
            self.atomgroup = universe_or_atomgroup.select_atoms(select)

        self._trajectory = self.universe.trajectory

        # Pre-compute radii (reused across all frames)
        self._radii = _get_radii_from_atomgroup(self.atomgroup)

        # Build atom-to-residue mapping. MDAnalysis resindices belong to the
        # whole Universe, so selections can be nonzero/discontiguous. Remap to
        # dense 0..N-1 columns for residue aggregation.
        raw_atom_to_residue = np.array(
            [atom.resindex for atom in self.atomgroup],
            dtype=np.int64,
        )
        _, dense_atom_to_residue = np.unique(raw_atom_to_residue, return_inverse=True)
        self._atom_to_residue = dense_atom_to_residue.astype(np.int32, copy=False)
        self._n_residues = int(self._atom_to_residue.max()) + 1 if self._atom_to_residue.size else 0

        # Results container
        self.results = _Results()

        # Frame info (set after run)
        self.n_frames: int = 0
        self.times: NDArray[np.float64] | None = None
        self.frames: NDArray[np.int64] | None = None

    def run(
        self,
        start: int = 0,
        stop: int | None = None,
        step: int = 1,
        *,
        probe_radius: float = 1.4,
        n_points: int = 960,
        algorithm: Literal["sr", "lr"] = "sr",
        n_slices: int = 20,
        n_threads: int = 0,
        chunk_size: int | None = None,
        store_atom_areas: bool = True,
        use_bitmask: bool = False,
        bitmask_correction: bool = False,
        bitmask_correction_coeff: float | None = None,
    ) -> SASAAnalysis:
        """Run the SASA analysis.

        Parameters
        ----------
        start : int, optional
            First frame to analyze (default: 0).
        stop : int, optional
            Last frame to analyze (default: None, meaning last frame).
        step : int, optional
            Step between frames (default: 1).
        probe_radius : float, optional
            Probe radius in Angstroms (default: 1.4).
        n_points : int, optional
            Number of points per atom for SR algorithm (default: 960).
        algorithm : {"sr", "lr"}, optional
            Algorithm: "sr" (Shrake-Rupley) or "lr" (Lee-Richards).
            Default: "sr".
        n_slices : int, optional
            Number of slices per atom for LR algorithm (default: 20).
        n_threads : int, optional
            Number of threads (0 = auto-detect). Default: 0.
        chunk_size : int, optional
            Number of frames to process per native batch. None processes all
            selected frames in one batch.
        store_atom_areas : bool, optional
            Store per-atom SASA in results.atom_area. Set False to retain only
            totals and residue aggregates.
        use_bitmask : bool, optional
            Use bitmask LUT optimization for SR algorithm.
            Supports n_points 1..1024. Default: False.
        bitmask_correction : bool, optional
            Apply experimental bitmask exposed-fraction correction.
            Requires use_bitmask=True. Default: False.
        bitmask_correction_coeff : float, optional
            Override the experimental correction coefficient.
            None uses the library default.

        Returns
        -------
        SASAAnalysis
            Self, for method chaining.
        """
        # Determine frame range
        if stop is None:
            stop = len(self._trajectory)

        frame_indices = list(range(start, stop, step))
        self.n_frames = len(frame_indices)

        if self.n_frames == 0:
            msg = "No frames to analyze"
            raise ValueError(msg)
        if chunk_size is None:
            chunk_size = self.n_frames
        if chunk_size <= 0:
            msg = f"chunk_size must be positive, got {chunk_size}"
            raise ValueError(msg)

        n_atoms = len(self.atomgroup)
        times = np.zeros(self.n_frames, dtype=np.float64)
        frames = np.zeros(self.n_frames, dtype=np.int64)
        atom_chunks: list[NDArray[np.float32]] = []
        residue_chunks: list[NDArray[np.float32]] = []
        total_chunks: list[NDArray[np.float32]] = []

        def calculate_chunk(coords: NDArray[np.float32]) -> NDArray[np.float32]:
            if algorithm == "sr":
                return calculate_sasa_batch(
                    coords,
                    self._radii,
                    algorithm="sr",
                    n_points=n_points,
                    probe_radius=probe_radius,
                    n_threads=n_threads,
                    use_bitmask=use_bitmask,
                    bitmask_correction=bitmask_correction,
                    bitmask_correction_coeff=bitmask_correction_coeff,
                ).atom_areas
            if algorithm == "lr":
                return calculate_sasa_batch(
                    coords,
                    self._radii,
                    algorithm="lr",
                    n_slices=n_slices,
                    probe_radius=probe_radius,
                    n_threads=n_threads,
                    use_bitmask=use_bitmask,
                    bitmask_correction=bitmask_correction,
                    bitmask_correction_coeff=bitmask_correction_coeff,
                ).atom_areas
            msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
            raise ValueError(msg)

        for chunk_start in range(0, self.n_frames, chunk_size):
            chunk_indices = frame_indices[chunk_start : chunk_start + chunk_size]
            coords = np.zeros((len(chunk_indices), n_atoms, 3), dtype=np.float32)

            for offset, frame_idx in enumerate(chunk_indices):
                self._trajectory[frame_idx]
                coords[offset] = self.atomgroup.positions.astype(np.float32)
                times[chunk_start + offset] = self._trajectory.time
                frames[chunk_start + offset] = frame_idx

            atom_areas = calculate_chunk(coords)
            if store_atom_areas:
                atom_chunks.append(atom_areas)

            residue_areas = np.zeros((atom_areas.shape[0], self._n_residues), dtype=np.float32)
            np.add.at(
                residue_areas,
                (np.arange(atom_areas.shape[0])[:, np.newaxis], self._atom_to_residue),
                atom_areas,
            )
            residue_chunks.append(residue_areas)
            total_chunks.append(atom_areas.sum(axis=1))

        self.times = times
        self.frames = frames

        self.results.atom_area = np.concatenate(atom_chunks) if store_atom_areas else None
        self.results.residue_area = np.concatenate(residue_chunks)
        self.results.total_area = np.concatenate(total_chunks)
        self.results.mean_total_area = float(self.results.total_area.mean())

        return self


class _Results:
    """Container for SASA analysis results."""

    def __init__(self) -> None:
        self.atom_area: NDArray[np.float32] | None = None
        self.residue_area: NDArray[np.float32] | None = None
        self.total_area: NDArray[np.float32] | None = None
        self.mean_total_area: float | None = None

    def __getitem__(self, key: str) -> NDArray[np.float32] | float | None:
        """Allow dict-like access for compatibility."""
        return getattr(self, key, None)


# Convenience function (alternative to class-based API)
def compute_sasa(
    universe_or_atomgroup: Universe | AtomGroup,
    *,
    select: str = "all",
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    probe_radius: float = 1.4,
    n_points: int = 960,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    chunk_size: int | None = None,
    mode: Literal["atom", "residue", "total"] = "atom",
    use_bitmask: bool = False,
    bitmask_correction: bool = False,
    bitmask_correction_coeff: float | None = None,
) -> NDArray[np.float32]:
    """Compute SASA for MDAnalysis Universe or AtomGroup.

    This is a convenience function that wraps SASAAnalysis for simple use cases.

    Parameters
    ----------
    universe_or_atomgroup : Universe or AtomGroup
        MDAnalysis Universe or AtomGroup to analyze.
    select : str, optional
        Atom selection string (default: "all").
    start : int, optional
        First frame to analyze (default: 0).
    stop : int, optional
        Last frame to analyze (default: None).
    step : int, optional
        Step between frames (default: 1).
    probe_radius : float, optional
        Probe radius in Angstroms (default: 1.4).
    n_points : int, optional
        Number of points per atom for SR algorithm (default: 960).
    algorithm : {"sr", "lr"}, optional
        Algorithm to use (default: "sr").
    n_slices : int, optional
        Number of slices for LR algorithm (default: 20).
    n_threads : int, optional
        Number of threads (0 = auto). Default: 0.
    chunk_size : int, optional
        Number of frames to process per native batch. None processes all
        selected frames in one batch.
    mode : {"atom", "residue", "total"}, optional
        Output mode (default: "atom").
    use_bitmask : bool, optional
        Use bitmask LUT optimization for SR algorithm.
        Supports n_points 1..1024. Default: False.
    bitmask_correction : bool, optional
        Apply experimental bitmask exposed-fraction correction.
        Requires use_bitmask=True. Default: False.
    bitmask_correction_coeff : float, optional
        Override the experimental correction coefficient.
        None uses the library default.

    Returns
    -------
    numpy.ndarray
        SASA values in Å².
        Shape depends on mode:
        - mode="atom": (n_frames, n_atoms)
        - mode="residue": (n_frames, n_residues)
        - mode="total": (n_frames,)
    """
    analysis = SASAAnalysis(universe_or_atomgroup, select=select)
    analysis.run(
        start=start,
        stop=stop,
        step=step,
        probe_radius=probe_radius,
        n_points=n_points,
        algorithm=algorithm,
        n_slices=n_slices,
        n_threads=n_threads,
        chunk_size=chunk_size,
        store_atom_areas=mode == "atom",
        use_bitmask=use_bitmask,
        bitmask_correction=bitmask_correction,
        bitmask_correction_coeff=bitmask_correction_coeff,
    )

    if mode == "total":
        return analysis.results.total_area
    if mode == "residue":
        return analysis.results.residue_area
    return analysis.results.atom_area

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
import re
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from zsasa._ffi import _validate_frame_selection
from zsasa.sasa import calculate_sasa_batch

if TYPE_CHECKING:
    from MDAnalysis.core.groups import AtomGroup
    from MDAnalysis.core.universe import Universe

logger = logging.getLogger(__name__)

# Import the MDAnalysis vdwradii table for consistency with the MDAnalysis ecosystem.
# Keys are uppercase element symbols (e.g., "C", "N", "CA" for calcium).
# The table lives in MDAnalysis.guesser.tables from MDAnalysis 2.8 on; earlier
# releases (2.0 to 2.7, which pyproject.toml allows) only have it in
# MDAnalysis.topology.tables, which 2.8 deprecates.
try:
    from MDAnalysis.guesser.tables import vdwradii as _MDA_VDW_RADII
except ImportError:
    try:
        from MDAnalysis.topology.tables import vdwradii as _MDA_VDW_RADII
    except ImportError:
        # MDAnalysis is not installed (the module is importable without it); use
        # the values MDAnalysis has for the elements this module can infer.
        _MDA_VDW_RADII: dict[str, float] = {
            "H": 1.1,
            "C": 1.7,
            "N": 1.55,
            "O": 1.52,
            "F": 1.47,
            "NA": 2.27,
            "MG": 1.73,
            "P": 1.8,
            "S": 1.8,
            "CL": 1.75,
            "K": 2.75,
            "CA": 2.31,
            "ZN": 1.39,
            "SE": 1.9,
            "BR": 1.85,
            "I": 1.98,
        }

# Default radius for unknown elements (conservative estimate, larger than most common elements)
_DEFAULT_RADIUS = 2.0

# Two-letter symbols _get_element may return when there is no element attribute
# (ions, halogens and metals of biomolecular systems).
_TWO_LETTER_ELEMENTS = frozenset(
    {
        "LI",
        "NA",
        "MG",
        "AL",
        "SI",
        "CL",
        "CA",
        "MN",
        "FE",
        "CO",
        "NI",
        "CU",
        "ZN",
        "SE",
        "BR",
        "RB",
        "SR",
        "CD",
        "CS",
        "BA",
        "HG",
    }
)

# Two-letter symbols that no protein, nucleic acid or common ligand atom name
# starts with, so a name or type that begins with one is that element. The others
# are also ordinary atom names (CA alpha carbon, CD/CE/CG side-chain carbons, NA
# nitrogen A of a heme, HG/HE/HD hydrogens, SE selenium or sulfur, CO carbonyl) and
# are only taken as the element when the residue name says the atom is an ion.
_UNAMBIGUOUS_TWO_LETTER = frozenset({"CL", "BR", "FE", "ZN", "MG", "MN", "NI", "CU", "LI", "AL"})

# CHARMM names for single-atom ions, which do not start with their element symbol.
_ION_ALIASES = {
    "SOD": "NA",
    "POT": "K",
    "CLA": "CL",
    "CAL": "CA",
    "CES": "CS",
    "LIT": "LI",
    "RUB": "RB",
    "BAR": "BA",
}

_LEADING_LETTERS = re.compile(r"[^A-Za-z]*([A-Za-z]+)")


def _letters(text: object) -> str:
    """Return the first run of letters in ``text``, upper-cased, skipping leading numbers.

    ``"1HB"`` gives ``"HB"``, ``"ZN2+"`` gives ``"ZN"``, ``"O5'"`` gives ``"O"``.
    """
    if text is None:
        return ""
    match = _LEADING_LETTERS.match(str(text).strip())
    return match.group(1).upper() if match else ""


def _element_from_label(label: object, resname: str) -> str:
    """Infer an upper-case element symbol from an atom type or atom name.

    The first letter is the element unless the label starts with a two-letter
    symbol that is more than an ordinary atom name (see ``_infer_element``).
    """
    letters = _letters(label)
    if not letters:
        return ""
    if letters in _ION_ALIASES:
        return _ION_ALIASES[letters]

    two = letters[:2]
    if (
        len(letters) >= 2
        and two in _TWO_LETTER_ELEMENTS
        and (two == resname or two in _UNAMBIGUOUS_TWO_LETTER)
    ):
        return two
    return letters[0]


def _infer_element(atom) -> str:  # noqa: ANN001
    """Infer the element of an atom that has no ``element`` attribute.

    Taking the first letter is not enough: ``CL``, ``FE``, ``ZN`` and ``MG`` must
    not become carbon, fluorine and so on, while the alpha carbon ``CA`` of a
    protein must stay carbon. Which reading is right depends on what the topology
    provides. The atom type and the atom name are each read as follows:

    1. CHARMM ion names (``SOD``, ``CLA``, ``CAL`` ...) give their element.
    2. A two-letter symbol that is never an ordinary atom name (``CL``, ``BR``, ``FE``,
       ``ZN``, ``MG``, ``MN``, ``NI``, ``CU``, ``LI``, ``AL``) is that element.
    3. A two-letter symbol that is also an ordinary atom name (``CA``, ``CD``, ``NA``,
       ``HG``, ``SE``, ``CO`` ...) is that element only when the residue name is the
       same symbol (``CA`` in a residue ``CA``, ``NA`` in ``NA+``); inside ``ALA``,
       ``GLN`` or ``HEM`` it is the first letter (C, C, N).
    4. Otherwise the first letter. Force-field types such as ``CT1``, ``NH1`` or
       ``HGA2`` are not element symbols and end here.

    A two-letter reading from either the type or the name wins, because the type is
    often guessed from the name by the topology parser and can be wrong for an ion
    (MDAnalysis types the ion ``CA`` in residue ``CA`` as carbon). Otherwise the
    type decides, then the name. The mass is not used: when a topology has no
    masses MDAnalysis guesses them from these same types, so a mass cannot confirm
    them.

    Returns an upper-case symbol, or "" when the atom has neither type nor name.
    """
    try:
        resname = _letters(atom.resname)
    except (AttributeError, TypeError):
        resname = ""

    symbols = []
    for attribute in ("type", "name"):
        try:
            label = getattr(atom, attribute)
        except (AttributeError, TypeError):
            continue
        symbol = _element_from_label(label, resname)
        if symbol:
            symbols.append(symbol)

    for symbol in symbols:
        if len(symbol) == 2:
            return symbol
    return symbols[0] if symbols else ""


def _get_element(atom) -> str:  # noqa: ANN001
    """Get element symbol from MDAnalysis atom.

    Uses the ``element`` attribute when the topology has one. Otherwise the element
    is inferred from the atom type, name and residue name (see ``_infer_element``);
    carbon is the last resort.
    """
    # Try element attribute first
    try:
        if hasattr(atom, "element") and atom.element:
            return str(atom.element).capitalize()
    except (AttributeError, TypeError):
        pass

    symbol = _infer_element(atom)
    return symbol.capitalize() if symbol else "C"  # Default to carbon


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
            First frame to analyze (default: 0; must not be negative).
        stop : int, optional
            Stop before this frame (default: None, meaning run through the
            last frame; must not be negative).
        step : int, optional
            Step between frames (default: 1; must be at least 1, otherwise
            ValueError is raised before any frame is read).
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
        _validate_frame_selection(start, stop, step)

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
        First frame to analyze (default: 0; must not be negative).
    stop : int, optional
        Stop before this frame (default: None, meaning run through the last
        frame; must not be negative).
    step : int, optional
        Step between frames (default: 1; must be at least 1).
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

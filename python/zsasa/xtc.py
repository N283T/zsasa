"""Native XTC trajectory reader using zsasa's high-performance Zig implementation.

This module provides a standalone XTC reader that doesn't require MDTraj or MDAnalysis,
using zsasa's native Zig XTC implementation for potentially better performance.

Example:
    >>> from zsasa.xtc import XtcReader, compute_sasa_trajectory
    >>>
    >>> # Low-level reader API
    >>> with XtcReader("trajectory.xtc") as reader:
    ...     for frame in reader:
    ...         print(f"Step {frame.step}, {frame.natoms} atoms")
    >>>
    >>> # High-level SASA calculation
    >>> from zsasa.xtc import compute_sasa_trajectory
    >>> result = compute_sasa_trajectory(
    ...     "trajectory.xtc",
    ...     "topology.pdb",  # for radii
    ... )
    >>> print(result.total_areas)  # Per-frame total SASA

Note:
    Coordinates are returned in nanometers (nm), matching GROMACS convention.
    For SASA calculation, coordinates are automatically converted to Angstroms.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from zsasa._ffi import _get_lib
from zsasa.sasa import calculate_sasa_batch

if TYPE_CHECKING:
    from cffi import FFI


# Error codes from C API
_ZSASA_OK = 0
_ZSASA_XTC_END_OF_FILE = 1
_ZSASA_ERROR_INVALID_INPUT = -1
_ZSASA_ERROR_OUT_OF_MEMORY = -2


@dataclass
class XtcFrame:
    """A single frame from an XTC trajectory.

    Attributes:
        step: Simulation step number.
        time: Simulation time in ps.
        coords: Atom coordinates as (n_atoms, 3) array in nanometers.
        box: Simulation box as 3x3 matrix in nanometers.
        precision: Compression precision.
    """

    step: int
    time: float
    coords: NDArray[np.float32]
    box: NDArray[np.float32]
    precision: float

    @property
    def natoms(self) -> int:
        """Number of atoms."""
        return self.coords.shape[0]


class XtcReader:
    """Native XTC trajectory reader.

    This reader uses zsasa's high-performance Zig XTC implementation,
    providing a standalone alternative to MDTraj/MDAnalysis for reading
    XTC files.

    Parameters
    ----------
    path : str or Path
        Path to XTC trajectory file.

    Attributes
    ----------
    natoms : int
        Number of atoms in the trajectory.

    Example
    -------
    >>> with XtcReader("trajectory.xtc") as reader:
    ...     print(f"Trajectory has {reader.natoms} atoms")
    ...     for frame in reader:
    ...         print(f"Step {frame.step}: {frame.coords.shape}")

    >>> # Or without context manager (remember to close!)
    >>> reader = XtcReader("trajectory.xtc")
    >>> frame = reader.read_frame()
    >>> reader.close()
    """

    def __init__(self, path: str | Path) -> None:
        """Open an XTC file for reading."""
        self._ffi: FFI
        self._lib: object
        self._ffi, self._lib = _get_lib()

        self._path = str(path)
        self._handle = None
        self._natoms = 0
        self._closed = False

        # Open the file
        natoms_out = self._ffi.new("int*")
        error_code = self._ffi.new("int*")

        self._handle = self._lib.zsasa_xtc_open(
            self._path.encode("utf-8"),
            natoms_out,
            error_code,
        )

        if self._handle == self._ffi.NULL:
            if error_code[0] == _ZSASA_ERROR_INVALID_INPUT:
                msg = f"Cannot open XTC file: {self._path}"
                raise FileNotFoundError(msg)
            elif error_code[0] == _ZSASA_ERROR_OUT_OF_MEMORY:
                msg = "Out of memory opening XTC file"
                raise MemoryError(msg)
            else:
                msg = f"Error opening XTC file: {error_code[0]}"
                raise RuntimeError(msg)

        self._natoms = natoms_out[0]

        # Pre-allocate buffers for reading frames
        self._coords_buffer = np.zeros(self._natoms * 3, dtype=np.float32)
        self._box_buffer = np.zeros(9, dtype=np.float32)

    @property
    def natoms(self) -> int:
        """Number of atoms in the trajectory."""
        return self._natoms

    def read_frame(self) -> XtcFrame | None:
        """Read the next frame from the trajectory.

        Returns
        -------
        XtcFrame or None
            The next frame, or None if end of file is reached.

        Raises
        ------
        RuntimeError
            If an error occurs during reading.
        """
        if self._closed:
            msg = "XtcReader is closed"
            raise RuntimeError(msg)

        step = self._ffi.new("int*")
        time = self._ffi.new("float*")
        precision = self._ffi.new("float*")

        coords_ptr = self._ffi.cast("float*", self._coords_buffer.ctypes.data)
        box_ptr = self._ffi.cast("float*", self._box_buffer.ctypes.data)

        result = self._lib.zsasa_xtc_read_frame(
            self._handle,
            coords_ptr,
            step,
            time,
            box_ptr,
            precision,
        )

        if result == _ZSASA_XTC_END_OF_FILE:
            return None
        elif result != _ZSASA_OK:
            msg = f"Error reading XTC frame: {result}"
            raise RuntimeError(msg)

        # Copy buffers to new arrays (so they're independent of the reader)
        coords = self._coords_buffer.copy().reshape(self._natoms, 3)
        box = self._box_buffer.copy().reshape(3, 3)

        return XtcFrame(
            step=step[0],
            time=time[0],
            coords=coords,
            box=box,
            precision=precision[0],
        )

    def __iter__(self) -> Iterator[XtcFrame]:
        """Iterate over all frames in the trajectory."""
        while True:
            frame = self.read_frame()
            if frame is None:
                break
            yield frame

    def __enter__(self) -> XtcReader:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the XTC file and release resources."""
        if not self._closed and self._handle is not None:
            self._lib.zsasa_xtc_close(self._handle)
            self._handle = None
            self._closed = True

    def __del__(self) -> None:
        """Destructor - ensure file is closed."""
        self.close()


@dataclass
class TrajectorySasaResult:
    """Result of trajectory SASA calculation.

    Attributes:
        atom_areas: Per-atom SASA for all frames, shape (n_frames, n_atoms) in Å².
        steps: Step numbers for each frame.
        times: Time values for each frame in ps.
    """

    atom_areas: NDArray[np.float32]
    steps: NDArray[np.int32]
    times: NDArray[np.float32]

    @property
    def n_frames(self) -> int:
        """Number of frames."""
        return self.atom_areas.shape[0]

    @property
    def n_atoms(self) -> int:
        """Number of atoms."""
        return self.atom_areas.shape[1]

    @property
    def total_areas(self) -> NDArray[np.float32]:
        """Total SASA per frame, shape (n_frames,)."""
        return self.atom_areas.sum(axis=1)

    def __repr__(self) -> str:
        return f"TrajectorySasaResult(n_frames={self.n_frames}, n_atoms={self.n_atoms})"


def compute_sasa_trajectory(
    xtc_path: str | Path,
    radii: NDArray[np.floating] | list[float],
    *,
    probe_radius: float = 1.4,
    n_points: int = 100,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    use_bitmask: bool = False,
) -> TrajectorySasaResult:
    """Compute SASA for an XTC trajectory using native Zig reader.

    This function reads an XTC file directly using zsasa's native Zig XTC reader,
    without requiring MDTraj or MDAnalysis.

    Parameters
    ----------
    xtc_path : str or Path
        Path to XTC trajectory file.
    radii : array-like
        Atomic radii in Angstroms, shape (n_atoms,).
        Must match the number of atoms in the trajectory.
    probe_radius : float, optional
        Water probe radius in Angstroms. Default: 1.4.
    n_points : int, optional
        Number of test points per atom (SR algorithm). Default: 100.
    algorithm : {"sr", "lr"}, optional
        Algorithm: "sr" (Shrake-Rupley) or "lr" (Lee-Richards). Default: "sr".
    n_slices : int, optional
        Number of slices per atom (LR algorithm). Default: 20.
    n_threads : int, optional
        Number of threads (0 = auto-detect). Default: 0.
    start : int, optional
        First frame to process (0-indexed). Default: 0.
    stop : int, optional
        Stop before this frame (exclusive). Default: None (all frames).
    step : int, optional
        Process every Nth frame. Default: 1.
    use_bitmask : bool, optional
        Use bitmask LUT optimization for SR algorithm.
        Supports n_points 1..1024. Default: False.

    Returns
    -------
    TrajectorySasaResult
        Result containing per-atom SASA values for all frames.

    Note
    ----
    XTC coordinates are in nanometers. This function automatically converts
    them to Angstroms (x10) for SASA calculation. Output SASA values are in Å².

    Example
    -------
    >>> import numpy as np
    >>> from zsasa.xtc import compute_sasa_trajectory
    >>>
    >>> # Define radii for each atom (e.g., from a topology file)
    >>> radii = np.array([1.7, 1.55, 1.52, ...])  # Carbon, Nitrogen, Oxygen, etc.
    >>>
    >>> result = compute_sasa_trajectory("trajectory.xtc", radii)
    >>> print(f"Processed {result.n_frames} frames")
    >>> print(f"Total SASA: {result.total_areas}")
    """
    radii = np.asarray(radii, dtype=np.float32)

    # Read all frames first (to know total count and validate radii)
    frames: list[XtcFrame] = []
    steps: list[int] = []
    times: list[float] = []

    with XtcReader(xtc_path) as reader:
        # Validate radii length
        if len(radii) != reader.natoms:
            msg = f"radii length ({len(radii)}) doesn't match trajectory atoms ({reader.natoms})"
            raise ValueError(msg)

        frame_idx = 0
        for frame in reader:
            # Check if this frame should be processed
            if frame_idx < start:
                frame_idx += 1
                continue
            if stop is not None and frame_idx >= stop:
                break
            if (frame_idx - start) % step != 0:
                frame_idx += 1
                continue

            frames.append(frame)
            steps.append(frame.step)
            times.append(frame.time)
            frame_idx += 1

    if not frames:
        msg = "No frames to process"
        raise ValueError(msg)

    # Stack coordinates: (n_frames, n_atoms, 3) and convert nm -> Angstrom
    coords = np.stack([f.coords for f in frames], axis=0)
    coords_angstrom = coords * 10.0  # nm -> Angstrom

    # Calculate SASA using batch API
    result = calculate_sasa_batch(
        coords_angstrom,
        radii,
        algorithm=algorithm,
        n_points=n_points,
        n_slices=n_slices,
        probe_radius=probe_radius,
        n_threads=n_threads,
        use_bitmask=use_bitmask,
    )

    return TrajectorySasaResult(
        atom_areas=result.atom_areas,
        steps=np.array(steps, dtype=np.int32),
        times=np.array(times, dtype=np.float32),
    )

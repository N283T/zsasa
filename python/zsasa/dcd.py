"""Native DCD trajectory reader using zsasa's high-performance Zig implementation.

This module provides a standalone DCD reader that doesn't require MDTraj or MDAnalysis,
using zsasa's native Zig DCD implementation.

DCD is a binary trajectory format used by NAMD and CHARMM. Unlike XTC, DCD
stores coordinates in Angstroms (not nanometers).

Example:
    >>> from zsasa.dcd import DcdReader, compute_sasa_trajectory
    >>>
    >>> # Low-level reader API
    >>> with DcdReader("trajectory.dcd") as reader:
    ...     for frame in reader:
    ...         print(f"Step {frame.step}, {frame.natoms} atoms")
    >>>
    >>> # High-level SASA calculation
    >>> result = compute_sasa_trajectory("trajectory.dcd", radii)
    >>> print(result.total_areas)  # Per-frame total SASA

Note:
    DCD coordinates are in Angstroms, matching SASA calculation convention.
    No unit conversion is needed (unlike XTC which requires nm -> Angstrom).
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

from zsasa.core import _get_lib, calculate_sasa_batch
from zsasa.xtc import TrajectorySasaResult

if TYPE_CHECKING:
    from cffi import FFI


# Error codes from C API
_ZSASA_OK = 0
_ZSASA_DCD_END_OF_FILE = 2
_ZSASA_ERROR_INVALID_INPUT = -1
_ZSASA_ERROR_OUT_OF_MEMORY = -2


@dataclass
class DcdFrame:
    """A single frame from a DCD trajectory.

    Attributes:
        step: Simulation step number.
        time: Simulation time.
        coords: Atom coordinates as (n_atoms, 3) array in Angstroms.
        unitcell: Unit cell parameters (6 doubles), or None if not present.
    """

    step: int
    time: float
    coords: NDArray[np.float32]
    unitcell: NDArray[np.float64] | None

    @property
    def natoms(self) -> int:
        """Number of atoms."""
        return self.coords.shape[0]


class DcdReader:
    """Native DCD trajectory reader.

    This reader uses zsasa's high-performance Zig DCD implementation,
    providing a standalone alternative to MDTraj/MDAnalysis for reading
    DCD files.

    Parameters
    ----------
    path : str or Path
        Path to DCD trajectory file.

    Attributes
    ----------
    natoms : int
        Number of atoms in the trajectory.

    Example
    -------
    >>> with DcdReader("trajectory.dcd") as reader:
    ...     print(f"Trajectory has {reader.natoms} atoms")
    ...     for frame in reader:
    ...         print(f"Step {frame.step}: {frame.coords.shape}")

    >>> # Or without context manager (remember to close!)
    >>> reader = DcdReader("trajectory.dcd")
    >>> frame = reader.read_frame()
    >>> reader.close()
    """

    def __init__(self, path: str | Path) -> None:
        """Open a DCD file for reading."""
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

        self._handle = self._lib.zsasa_dcd_open(
            self._path.encode("utf-8"),
            natoms_out,
            error_code,
        )

        if self._handle == self._ffi.NULL:
            if error_code[0] == _ZSASA_ERROR_INVALID_INPUT:
                msg = f"Cannot open DCD file: {self._path}"
                raise FileNotFoundError(msg)
            elif error_code[0] == _ZSASA_ERROR_OUT_OF_MEMORY:
                msg = "Out of memory opening DCD file"
                raise MemoryError(msg)
            else:
                msg = f"Error opening DCD file: {error_code[0]}"
                raise RuntimeError(msg)

        self._natoms = natoms_out[0]

        # Pre-allocate buffers for reading frames
        self._coords_buffer = np.zeros(self._natoms * 3, dtype=np.float32)
        self._unitcell_buffer = np.zeros(6, dtype=np.float64)

    @property
    def natoms(self) -> int:
        """Number of atoms in the trajectory."""
        return self._natoms

    def read_frame(self) -> DcdFrame | None:
        """Read the next frame from the trajectory.

        Returns
        -------
        DcdFrame or None
            The next frame, or None if end of file is reached.

        Raises
        ------
        RuntimeError
            If an error occurs during reading.
        """
        if self._closed:
            msg = "DcdReader is closed"
            raise RuntimeError(msg)

        step = self._ffi.new("int*")
        time = self._ffi.new("float*")

        coords_ptr = self._ffi.cast("float*", self._coords_buffer.ctypes.data)
        unitcell_ptr = self._ffi.cast("double*", self._unitcell_buffer.ctypes.data)

        result = self._lib.zsasa_dcd_read_frame(
            self._handle,
            coords_ptr,
            step,
            time,
            unitcell_ptr,
        )

        if result == _ZSASA_DCD_END_OF_FILE:
            return None
        elif result != _ZSASA_OK:
            error_messages = {
                _ZSASA_ERROR_INVALID_INPUT: "invalid or corrupt DCD frame data",
                _ZSASA_ERROR_OUT_OF_MEMORY: "out of memory",
            }
            detail = error_messages.get(result, f"error code {result}")
            msg = f"Error reading DCD frame: {detail}"
            raise RuntimeError(msg)

        # Copy buffers to new arrays
        coords = self._coords_buffer.copy().reshape(self._natoms, 3)

        # Check if unitcell has non-zero values
        unitcell = self._unitcell_buffer.copy()
        has_unitcell = np.any(unitcell != 0.0)

        return DcdFrame(
            step=step[0],
            time=time[0],
            coords=coords,
            unitcell=unitcell if has_unitcell else None,
        )

    def __iter__(self) -> Iterator[DcdFrame]:
        """Iterate over all frames in the trajectory."""
        while True:
            frame = self.read_frame()
            if frame is None:
                break
            yield frame

    def __enter__(self) -> DcdReader:
        """Context manager entry."""
        return self

    def __exit__(self, exc_type: object, exc_val: object, exc_tb: object) -> None:
        """Context manager exit."""
        self.close()

    def close(self) -> None:
        """Close the DCD file and release resources."""
        if not self._closed and self._handle is not None:
            self._lib.zsasa_dcd_close(self._handle)
            self._handle = None
            self._closed = True

    def __del__(self) -> None:
        """Destructor - ensure file is closed."""
        self.close()


def compute_sasa_trajectory(
    dcd_path: str | Path,
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
    """Compute SASA for a DCD trajectory using native Zig reader.

    This function reads a DCD file directly using zsasa's native Zig DCD reader,
    without requiring MDTraj or MDAnalysis.

    Parameters
    ----------
    dcd_path : str or Path
        Path to DCD trajectory file.
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
        Only supports n_points of 64, 128, or 256. Default: False.

    Returns
    -------
    TrajectorySasaResult
        Result containing per-atom SASA values for all frames.

    Note
    ----
    DCD coordinates are in Angstroms. No unit conversion is needed
    (unlike XTC which stores coordinates in nanometers).
    Output SASA values are in Å².

    Example
    -------
    >>> import numpy as np
    >>> from zsasa.dcd import compute_sasa_trajectory
    >>>
    >>> radii = np.array([1.7, 1.55, 1.52, ...])
    >>> result = compute_sasa_trajectory("trajectory.dcd", radii)
    >>> print(f"Processed {result.n_frames} frames")
    >>> print(f"Total SASA: {result.total_areas}")
    """
    radii = np.asarray(radii, dtype=np.float32)

    # Read all frames first
    frames: list[DcdFrame] = []
    steps: list[int] = []
    times: list[float] = []

    with DcdReader(dcd_path) as reader:
        # Validate radii length
        if len(radii) != reader.natoms:
            msg = f"radii length ({len(radii)}) doesn't match trajectory atoms ({reader.natoms})"
            raise ValueError(msg)

        frame_idx = 0
        for frame in reader:
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

    # Stack coordinates: (n_frames, n_atoms, 3) — already in Angstroms
    coords_angstrom = np.stack([f.coords for f in frames], axis=0)

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

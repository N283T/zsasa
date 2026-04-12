"""SASA calculation functions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from zsasa._ffi import (
    _BITMASK_MAX_N_POINTS,
    _BITMASK_MIN_N_POINTS,
    ZSASA_ERROR_CALCULATION,
    ZSASA_ERROR_INVALID_INPUT,
    ZSASA_ERROR_OUT_OF_MEMORY,
    ZSASA_ERROR_UNSUPPORTED_N_POINTS,
    ZSASA_OK,
    _get_lib,
    _validate_bitmask_params,
)


@dataclass
class SasaResult:
    """Result of SASA calculation.

    Attributes:
        total_area: Total solvent accessible surface area in Å².
        atom_areas: Per-atom SASA values in Å².
    """

    total_area: float
    atom_areas: NDArray[np.float64]


def calculate_sasa(
    coords: NDArray[np.float64],
    radii: NDArray[np.float64],
    *,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    use_bitmask: bool = False,
) -> SasaResult:
    """Calculate Solvent Accessible Surface Area (SASA).

    Args:
        coords: Atom coordinates as (N, 3) array.
        radii: Atom radii as (N,) array.
        algorithm: Algorithm to use: "sr" (Shrake-Rupley) or "lr" (Lee-Richards).
        n_points: Number of test points per atom (for SR algorithm). Default: 100.
        n_slices: Number of slices per atom (for LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads to use. 0 = auto-detect. Default: 0.
        use_bitmask: Use bitmask LUT optimization for SR algorithm.
            Supports n_points 1..1024. Default: False.

    Returns:
        SasaResult containing total_area and per-atom atom_areas.

    Raises:
        ValueError: If input arrays have invalid shapes or calculation fails.

    Example:
        >>> import numpy as np
        >>> from zsasa import calculate_sasa
        >>>
        >>> # Two atoms
        >>> coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        >>> radii = np.array([1.5, 1.5])
        >>> result = calculate_sasa(coords, radii)
        >>> print(f"Total: {result.total_area:.2f}")
    """
    ffi, lib = _get_lib()

    # Validate parameters
    if n_points <= 0:
        msg = f"n_points must be positive, got {n_points}"
        raise ValueError(msg)
    if n_slices <= 0:
        msg = f"n_slices must be positive, got {n_slices}"
        raise ValueError(msg)
    if probe_radius <= 0:
        msg = f"probe_radius must be positive, got {probe_radius}"
        raise ValueError(msg)
    if n_threads < 0:
        msg = f"n_threads must be non-negative, got {n_threads}"
        raise ValueError(msg)

    # Validate bitmask constraints
    if use_bitmask:
        use_bitmask = _validate_bitmask_params(algorithm, n_points)

    # Validate and convert inputs
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    radii = np.ascontiguousarray(radii, dtype=np.float64)

    if coords.ndim != 2 or coords.shape[1] != 3:
        msg = f"coords must be (N, 3) array, got shape {coords.shape}"
        raise ValueError(msg)

    n_atoms = coords.shape[0]
    if radii.shape != (n_atoms,):
        msg = f"radii must be ({n_atoms},) array, got shape {radii.shape}"
        raise ValueError(msg)

    if np.any(radii < 0):
        msg = "All radii must be non-negative"
        raise ValueError(msg)

    # Extract x, y, z as contiguous arrays
    x = np.ascontiguousarray(coords[:, 0])
    y = np.ascontiguousarray(coords[:, 1])
    z = np.ascontiguousarray(coords[:, 2])

    # Allocate output arrays
    atom_areas = np.zeros(n_atoms, dtype=np.float64)
    total_area = ffi.new("double*")

    # Get cffi pointers from numpy arrays
    x_ptr = ffi.cast("double*", x.ctypes.data)
    y_ptr = ffi.cast("double*", y.ctypes.data)
    z_ptr = ffi.cast("double*", z.ctypes.data)
    radii_ptr = ffi.cast("double*", radii.ctypes.data)
    areas_ptr = ffi.cast("double*", atom_areas.ctypes.data)

    # Call the appropriate function
    if use_bitmask:
        result = lib.zsasa_calc_sr_bitmask(
            x_ptr,
            y_ptr,
            z_ptr,
            radii_ptr,
            n_atoms,
            n_points,
            probe_radius,
            n_threads,
            areas_ptr,
            total_area,
        )
    elif algorithm == "sr":
        result = lib.zsasa_calc_sr(
            x_ptr,
            y_ptr,
            z_ptr,
            radii_ptr,
            n_atoms,
            n_points,
            probe_radius,
            n_threads,
            areas_ptr,
            total_area,
        )
    elif algorithm == "lr":
        result = lib.zsasa_calc_lr(
            x_ptr,
            y_ptr,
            z_ptr,
            radii_ptr,
            n_atoms,
            n_slices,
            probe_radius,
            n_threads,
            areas_ptr,
            total_area,
        )
    else:
        msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
        raise ValueError(msg)

    # Check for errors
    if result == ZSASA_ERROR_INVALID_INPUT:
        msg = "Invalid input to SASA calculation"
        raise ValueError(msg)
    elif result == ZSASA_ERROR_OUT_OF_MEMORY:
        msg = "Out of memory during SASA calculation"
        raise MemoryError(msg)
    elif result == ZSASA_ERROR_CALCULATION:
        msg = "Error during SASA calculation"
        raise RuntimeError(msg)
    elif result == ZSASA_ERROR_UNSUPPORTED_N_POINTS:
        msg = (
            f"Unsupported n_points for bitmask: {n_points}. "
            f"Must be {_BITMASK_MIN_N_POINTS}..{_BITMASK_MAX_N_POINTS}"
        )
        raise ValueError(msg)
    elif result != ZSASA_OK:
        msg = f"Unknown error code: {result}"
        raise RuntimeError(msg)

    return SasaResult(
        total_area=total_area[0],
        atom_areas=atom_areas,
    )


@dataclass
class BatchSasaResult:
    """Result of batch SASA calculation for multiple frames.

    Attributes:
        atom_areas: Per-atom SASA values for all frames, shape (n_frames, n_atoms).
                    Values are in Angstrom² (Å²).
    """

    atom_areas: NDArray[np.float32]

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
        return f"BatchSasaResult(n_frames={self.n_frames}, n_atoms={self.n_atoms})"


def calculate_sasa_batch(
    coordinates: NDArray[np.floating],
    radii: NDArray[np.floating],
    *,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    precision: Literal["f64", "f32"] = "f64",
    use_bitmask: bool = False,
) -> BatchSasaResult:
    """Calculate SASA for multiple frames (batch processing).

    Optimized for MD trajectory analysis where the same atoms are processed
    across multiple frames. Parallelizes across frames for maximum performance.

    Args:
        coordinates: Atom coordinates as (n_frames, n_atoms, 3) array in Angstroms.
        radii: Atom radii as (n_atoms,) array in Angstroms.
        algorithm: Algorithm to use: "sr" (Shrake-Rupley) or "lr" (Lee-Richards).
        n_points: Number of test points per atom (for SR algorithm). Default: 100.
        n_slices: Number of slices per atom (for LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads to use. 0 = auto-detect. Default: 0.
        precision: Internal calculation precision: "f64" (default, higher precision)
            or "f32" (matches RustSASA/mdsasa-bolt for comparison). Default: "f64".
        use_bitmask: Use bitmask LUT optimization for SR algorithm.
            Supports n_points 1..1024. Default: False.

    Returns:
        BatchSasaResult containing per-atom SASA for all frames.

    Raises:
        ValueError: If input arrays have invalid shapes or calculation fails.

    Example:
        >>> import numpy as np
        >>> from zsasa import calculate_sasa_batch
        >>>
        >>> # 10 frames, 100 atoms
        >>> coords = np.random.randn(10, 100, 3).astype(np.float32)
        >>> radii = np.full(100, 1.5, dtype=np.float32)
        >>> result = calculate_sasa_batch(coords, radii)
        >>> print(f"Shape: {result.atom_areas.shape}")  # (10, 100)
        >>> print(f"Total SASA per frame: {result.total_areas}")
    """
    ffi, lib = _get_lib()

    # Validate parameters
    if n_points <= 0:
        msg = f"n_points must be positive, got {n_points}"
        raise ValueError(msg)
    if n_slices <= 0:
        msg = f"n_slices must be positive, got {n_slices}"
        raise ValueError(msg)
    if probe_radius <= 0:
        msg = f"probe_radius must be positive, got {probe_radius}"
        raise ValueError(msg)
    if n_threads < 0:
        msg = f"n_threads must be non-negative, got {n_threads}"
        raise ValueError(msg)

    # Validate bitmask constraints
    if use_bitmask:
        use_bitmask = _validate_bitmask_params(algorithm, n_points)

    # Validate and convert inputs
    coordinates = np.ascontiguousarray(coordinates, dtype=np.float32)
    radii = np.ascontiguousarray(radii, dtype=np.float32)

    if coordinates.ndim != 3 or coordinates.shape[2] != 3:
        msg = f"coordinates must be (n_frames, n_atoms, 3) array, got shape {coordinates.shape}"
        raise ValueError(msg)

    n_frames = coordinates.shape[0]
    n_atoms = coordinates.shape[1]

    if radii.shape != (n_atoms,):
        msg = f"radii must be ({n_atoms},) array, got shape {radii.shape}"
        raise ValueError(msg)

    if np.any(radii < 0):
        msg = "All radii must be non-negative"
        raise ValueError(msg)

    # Allocate output array
    atom_areas = np.zeros((n_frames, n_atoms), dtype=np.float32)

    # Get cffi pointers from numpy arrays
    coords_ptr = ffi.cast("float*", coordinates.ctypes.data)
    radii_ptr = ffi.cast("float*", radii.ctypes.data)
    areas_ptr = ffi.cast("float*", atom_areas.ctypes.data)

    # Call the appropriate batch function based on algorithm, precision, and bitmask
    if use_bitmask:
        if precision == "f32":
            result = lib.zsasa_calc_sr_batch_bitmask_f32(
                coords_ptr,
                n_frames,
                n_atoms,
                radii_ptr,
                n_points,
                probe_radius,
                n_threads,
                areas_ptr,
            )
        else:
            result = lib.zsasa_calc_sr_batch_bitmask(
                coords_ptr,
                n_frames,
                n_atoms,
                radii_ptr,
                n_points,
                probe_radius,
                n_threads,
                areas_ptr,
            )
    elif precision == "f64":
        # Default: f32 I/O with f64 internal precision
        if algorithm == "sr":
            result = lib.zsasa_calc_sr_batch(
                coords_ptr,
                n_frames,
                n_atoms,
                radii_ptr,
                n_points,
                probe_radius,
                n_threads,
                areas_ptr,
            )
        elif algorithm == "lr":
            result = lib.zsasa_calc_lr_batch(
                coords_ptr,
                n_frames,
                n_atoms,
                radii_ptr,
                n_slices,
                probe_radius,
                n_threads,
                areas_ptr,
            )
        else:
            msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
            raise ValueError(msg)
    elif precision == "f32":
        # Pure f32 precision (matches RustSASA/mdsasa-bolt)
        if algorithm == "sr":
            result = lib.zsasa_calc_sr_batch_f32(
                coords_ptr,
                n_frames,
                n_atoms,
                radii_ptr,
                n_points,
                probe_radius,
                n_threads,
                areas_ptr,
            )
        elif algorithm == "lr":
            result = lib.zsasa_calc_lr_batch_f32(
                coords_ptr,
                n_frames,
                n_atoms,
                radii_ptr,
                n_slices,
                probe_radius,
                n_threads,
                areas_ptr,
            )
        else:
            msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
            raise ValueError(msg)
    else:
        msg = f"Unknown precision: {precision}. Use 'f64' or 'f32'."
        raise ValueError(msg)

    # Check for errors
    if result == ZSASA_ERROR_INVALID_INPUT:
        msg = "Invalid input to batch SASA calculation"
        raise ValueError(msg)
    elif result == ZSASA_ERROR_OUT_OF_MEMORY:
        msg = "Out of memory during batch SASA calculation"
        raise MemoryError(msg)
    elif result == ZSASA_ERROR_CALCULATION:
        msg = "Error during batch SASA calculation"
        raise RuntimeError(msg)
    elif result == ZSASA_ERROR_UNSUPPORTED_N_POINTS:
        msg = (
            f"Unsupported n_points for bitmask: {n_points}. "
            f"Must be {_BITMASK_MIN_N_POINTS}..{_BITMASK_MAX_N_POINTS}"
        )
        raise ValueError(msg)
    elif result != ZSASA_OK:
        msg = f"Unknown error code: {result}"
        raise RuntimeError(msg)

    return BatchSasaResult(atom_areas=atom_areas)

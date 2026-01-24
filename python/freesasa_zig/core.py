"""Core ctypes bindings for freesasa-zig library."""

from __future__ import annotations

import ctypes
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

# Error codes from C API
FREESASA_OK = 0
FREESASA_ERROR_INVALID_INPUT = -1
FREESASA_ERROR_OUT_OF_MEMORY = -2
FREESASA_ERROR_CALCULATION = -3


def _find_library() -> Path:
    """Find the freesasa_zig shared library."""
    # Check environment variable first
    if lib_path := os.environ.get("FREESASA_ZIG_LIB"):
        return Path(lib_path)

    # Platform-specific library names
    if sys.platform == "darwin":
        lib_name = "libfreesasa_zig.dylib"
    elif sys.platform == "win32":
        lib_name = "freesasa_zig.dll"
    else:
        lib_name = "libfreesasa_zig.so"

    # Search paths
    search_paths = [
        # Relative to this file (development)
        Path(__file__).parent.parent.parent / "zig-out" / "lib" / lib_name,
        # System paths
        Path("/usr/local/lib") / lib_name,
        Path("/usr/lib") / lib_name,
        # Current directory
        Path.cwd() / lib_name,
        Path.cwd() / "zig-out" / "lib" / lib_name,
    ]

    for path in search_paths:
        if path.exists():
            return path

    msg = (
        f"Could not find {lib_name}. "
        f"Set FREESASA_ZIG_LIB environment variable or build the library first."
    )
    raise FileNotFoundError(msg)


def _load_library() -> ctypes.CDLL:
    """Load the freesasa_zig shared library."""
    lib_path = _find_library()
    lib = ctypes.CDLL(str(lib_path))

    # Define function signatures

    # freesasa_version() -> const char*
    lib.freesasa_version.argtypes = []
    lib.freesasa_version.restype = ctypes.c_char_p

    # freesasa_calc_sr(x, y, z, radii, n_atoms, n_points, probe_radius, n_threads,
    #                  atom_areas, total_area) -> int
    lib.freesasa_calc_sr.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.POINTER(ctypes.c_double),  # z
        ctypes.POINTER(ctypes.c_double),  # radii
        ctypes.c_size_t,  # n_atoms
        ctypes.c_uint32,  # n_points
        ctypes.c_double,  # probe_radius
        ctypes.c_size_t,  # n_threads
        ctypes.POINTER(ctypes.c_double),  # atom_areas (output)
        ctypes.POINTER(ctypes.c_double),  # total_area (output)
    ]
    lib.freesasa_calc_sr.restype = ctypes.c_int

    # freesasa_calc_lr(x, y, z, radii, n_atoms, n_slices, probe_radius, n_threads,
    #                  atom_areas, total_area) -> int
    lib.freesasa_calc_lr.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.POINTER(ctypes.c_double),  # z
        ctypes.POINTER(ctypes.c_double),  # radii
        ctypes.c_size_t,  # n_atoms
        ctypes.c_uint32,  # n_slices
        ctypes.c_double,  # probe_radius
        ctypes.c_size_t,  # n_threads
        ctypes.POINTER(ctypes.c_double),  # atom_areas (output)
        ctypes.POINTER(ctypes.c_double),  # total_area (output)
    ]
    lib.freesasa_calc_lr.restype = ctypes.c_int

    return lib


# Global library instance (lazy loaded)
_lib: ctypes.CDLL | None = None


def _get_lib() -> ctypes.CDLL:
    """Get or load the library."""
    global _lib
    if _lib is None:
        _lib = _load_library()
    return _lib


def get_version() -> str:
    """Get the library version string."""
    lib = _get_lib()
    return lib.freesasa_version().decode("utf-8")


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
) -> SasaResult:
    """Calculate Solvent Accessible Surface Area (SASA).

    Args:
        coords: Atom coordinates as (N, 3) array or separate x, y, z arrays.
        radii: Atom radii as (N,) array.
        algorithm: Algorithm to use: "sr" (Shrake-Rupley) or "lr" (Lee-Richards).
        n_points: Number of test points per atom (for SR algorithm). Default: 100.
        n_slices: Number of slices per atom (for LR algorithm). Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads to use. 0 = auto-detect. Default: 0.

    Returns:
        SasaResult containing total_area and per-atom atom_areas.

    Raises:
        ValueError: If input arrays have invalid shapes or calculation fails.

    Example:
        >>> import numpy as np
        >>> from freesasa_zig import calculate_sasa
        >>>
        >>> # Two atoms
        >>> coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        >>> radii = np.array([1.5, 1.5])
        >>> result = calculate_sasa(coords, radii)
        >>> print(f"Total: {result.total_area:.2f}")
    """
    lib = _get_lib()

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

    # Extract x, y, z
    x = np.ascontiguousarray(coords[:, 0])
    y = np.ascontiguousarray(coords[:, 1])
    z = np.ascontiguousarray(coords[:, 2])

    # Allocate output arrays
    atom_areas = np.zeros(n_atoms, dtype=np.float64)
    total_area = ctypes.c_double(0.0)

    # Get pointers
    x_ptr = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    y_ptr = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    z_ptr = z.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    radii_ptr = radii.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    areas_ptr = atom_areas.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    # Call the appropriate function
    if algorithm == "sr":
        result = lib.freesasa_calc_sr(
            x_ptr,
            y_ptr,
            z_ptr,
            radii_ptr,
            n_atoms,
            n_points,
            probe_radius,
            n_threads,
            areas_ptr,
            ctypes.byref(total_area),
        )
    elif algorithm == "lr":
        result = lib.freesasa_calc_lr(
            x_ptr,
            y_ptr,
            z_ptr,
            radii_ptr,
            n_atoms,
            n_slices,
            probe_radius,
            n_threads,
            areas_ptr,
            ctypes.byref(total_area),
        )
    else:
        msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
        raise ValueError(msg)

    # Check for errors
    if result == FREESASA_ERROR_INVALID_INPUT:
        msg = "Invalid input to SASA calculation"
        raise ValueError(msg)
    elif result == FREESASA_ERROR_OUT_OF_MEMORY:
        msg = "Out of memory during SASA calculation"
        raise MemoryError(msg)
    elif result == FREESASA_ERROR_CALCULATION:
        msg = "Error during SASA calculation"
        raise RuntimeError(msg)
    elif result != FREESASA_OK:
        msg = f"Unknown error code: {result}"
        raise RuntimeError(msg)

    return SasaResult(
        total_area=total_area.value,
        atom_areas=atom_areas,
    )

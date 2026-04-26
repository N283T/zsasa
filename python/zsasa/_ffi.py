"""FFI bindings and library loading for zsasa."""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path
from typing import Any

from cffi import FFI

# Error codes from C API
ZSASA_OK = 0
ZSASA_ERROR_INVALID_INPUT = -1
ZSASA_ERROR_OUT_OF_MEMORY = -2
ZSASA_ERROR_CALCULATION = -3
ZSASA_ERROR_FILE_IO = -4
ZSASA_ERROR_UNSUPPORTED_N_POINTS = -5

# Valid n_points range for bitmask algorithm
_BITMASK_MIN_N_POINTS = 1
_BITMASK_MAX_N_POINTS = 1024


def _validate_bitmask_params(algorithm: str, n_points: int) -> bool:
    """Validate parameters for bitmask mode.

    Returns True if bitmask should be used, False if falling back to standard SR.
    Raises ValueError if algorithm is incompatible (not SR).
    """
    if algorithm != "sr":
        msg = "use_bitmask=True only supports algorithm='sr' (Shrake-Rupley)"
        raise ValueError(msg)
    if not (_BITMASK_MIN_N_POINTS <= n_points <= _BITMASK_MAX_N_POINTS):
        warnings.warn(
            f"use_bitmask=True requires n_points in "
            f"{_BITMASK_MIN_N_POINTS}..{_BITMASK_MAX_N_POINTS}, got {n_points}. "
            f"Falling back to standard Shrake-Rupley.",
            stacklevel=3,
        )
        return False
    return True


# Algorithm constants
ZSASA_ALGORITHM_SR = 0
ZSASA_ALGORITHM_LR = 1

# Classifier types
ZSASA_CLASSIFIER_NACCESS = 0
ZSASA_CLASSIFIER_PROTOR = 1
ZSASA_CLASSIFIER_OONS = 2
ZSASA_CLASSIFIER_CCD = 3

# Atom classes
ZSASA_ATOM_CLASS_POLAR = 0
ZSASA_ATOM_CLASS_APOLAR = 1
ZSASA_ATOM_CLASS_UNKNOWN = 2

# C API definitions for cffi
_CDEF = """
    // Version
    const char* zsasa_version(void);

    // SASA calculation
    int zsasa_calc_sr(
        const double* x, const double* y, const double* z, const double* radii,
        size_t n_atoms, uint32_t n_points, double probe_radius, size_t n_threads,
        double* atom_areas, double* total_area
    );

    int zsasa_calc_lr(
        const double* x, const double* y, const double* z, const double* radii,
        size_t n_atoms, uint32_t n_slices, double probe_radius, size_t n_threads,
        double* atom_areas, double* total_area
    );

    // Batch SASA calculation (for MD trajectories)
    int zsasa_calc_sr_batch(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    int zsasa_calc_lr_batch(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_slices, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    // Batch SASA calculation (pure f32 precision for RustSASA compatibility)
    int zsasa_calc_sr_batch_f32(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    int zsasa_calc_lr_batch_f32(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_slices, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    // Bitmask Shrake-Rupley (single frame, f64 internal)
    int zsasa_calc_sr_bitmask(
        const double* x, const double* y, const double* z, const double* radii,
        size_t n_atoms, uint32_t n_points, double probe_radius, size_t n_threads,
        double* atom_areas, double* total_area
    );

    // Bitmask batch SASA calculation (f64 internal precision)
    int zsasa_calc_sr_batch_bitmask(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    // Bitmask batch SASA calculation (f32 internal precision)
    int zsasa_calc_sr_batch_bitmask_f32(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    // Classifier functions
    double zsasa_classifier_get_radius(
        int classifier_type, const char* residue, const char* atom);
    int zsasa_classifier_get_class(
        int classifier_type, const char* residue, const char* atom);
    double zsasa_guess_radius(const char* element);
    double zsasa_guess_radius_from_atom_name(const char* atom_name);

    int zsasa_classify_atoms(
        int classifier_type,
        const char** residues, const char** atoms, size_t n_atoms,
        double* radii_out, int* classes_out
    );

    // RSA functions
    double zsasa_get_max_sasa(const char* residue_name);
    double zsasa_calculate_rsa(double sasa, const char* residue_name);
    int zsasa_calculate_rsa_batch(
        const double* sasas, const char** residue_names, size_t n_residues,
        double* rsa_out
    );

    // XTC trajectory reader
    void* zsasa_xtc_open(const char* path, int* natoms_out, int* error_code);
    void zsasa_xtc_close(void* handle);
    int zsasa_xtc_read_frame(
        void* handle,
        float* coords_out,
        int* step_out,
        float* time_out,
        float* box_out,
        float* precision_out
    );
    int zsasa_xtc_get_natoms(void* handle);

    // DCD trajectory reader
    void* zsasa_dcd_open(const char* path, int* natoms_out, int* error_code);
    void zsasa_dcd_close(void* handle);
    int zsasa_dcd_read_frame(
        void* handle,
        float* coords_out,
        int* step_out,
        float* time_out,
        double* unitcell_out
    );
    int zsasa_dcd_get_natoms(void* handle);

    // Batch directory processing
    void* zsasa_batch_dir_process(
        const char* input_dir, const char* output_dir,
        int algorithm, uint32_t n_points, double probe_radius,
        size_t n_threads, int classifier_type,
        int include_hydrogens, int include_hetatm,
        int* error_code
    );
    size_t zsasa_batch_dir_get_total_files(void* handle);
    size_t zsasa_batch_dir_get_successful(void* handle);
    size_t zsasa_batch_dir_get_failed(void* handle);
    const char* zsasa_batch_dir_get_filename(void* handle, size_t index);
    size_t zsasa_batch_dir_get_n_atoms(void* handle, size_t index);
    double zsasa_batch_dir_get_total_sasa(void* handle, size_t index);
    int zsasa_batch_dir_get_status(void* handle, size_t index);
    void zsasa_batch_dir_free(void* handle);
"""


def _find_library() -> Path:
    """Find the zsasa shared library."""
    # Check environment variable first
    if lib_path := os.environ.get("ZSASA_LIB"):
        return Path(lib_path)

    # Platform-specific library names
    if sys.platform == "darwin":
        lib_name = "libzsasa.dylib"
    elif sys.platform == "win32":
        lib_name = "zsasa.dll"
    else:
        lib_name = "libzsasa.so"

    # Search paths
    search_paths = [
        # Bundled in package (wheel installation)
        Path(__file__).parent / lib_name,
        # Relative to this file (development: python/zsasa -> zig-out/lib)
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
        f"Please install with: pip install zsasa "
        f"(requires Zig 0.16.0+ to be installed)"
    )
    raise FileNotFoundError(msg)


def _load_library() -> tuple[FFI, Any]:
    """Load the zsasa shared library using cffi."""
    ffi = FFI()
    ffi.cdef(_CDEF)
    lib_path = _find_library()
    lib = ffi.dlopen(str(lib_path))
    return ffi, lib


# Global library instance (lazy loaded)
_ffi: FFI | None = None
_lib: Any = None


def _get_lib() -> tuple[FFI, Any]:
    """Get or load the library."""
    global _ffi, _lib
    if _ffi is None:
        _ffi, _lib = _load_library()
    return _ffi, _lib


def get_version() -> str:
    """Get the library version string."""
    ffi, lib = _get_lib()
    return ffi.string(lib.zsasa_version()).decode("utf-8")

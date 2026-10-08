"""FFI bindings and library loading for zsasa."""

from __future__ import annotations

import filecmp
import os
import sys
import warnings
from pathlib import Path
from typing import Any, NoReturn

from cffi import FFI

# Error codes from C API
ZSASA_OK = 0
ZSASA_ERROR_INVALID_INPUT = -1
ZSASA_ERROR_OUT_OF_MEMORY = -2
ZSASA_ERROR_CALCULATION = -3
ZSASA_ERROR_FILE_IO = -4
ZSASA_ERROR_UNSUPPORTED_N_POINTS = -5
ZSASA_ERROR_OUTPUT_NAME_COLLISION = -6
ZSASA_ERROR_INVALID_FORMAT = -7
ZSASA_ERROR_OUTPUT_DIR = -8

# Version of the C ABI this package's cdef (below) was written for. It must equal
# the value zsasa_abi_version() returns in the loaded library. Bump it together with
# ABI_VERSION in src/c_api.zig whenever an existing exported signature or struct
# layout changes (not for new exports or new error codes).
_EXPECTED_ABI_VERSION = 1

# Valid n_points range for bitmask algorithm
_BITMASK_MIN_N_POINTS = 1
_BITMASK_MAX_N_POINTS = 1024


def _validate_bitmask_params(algorithm: str, n_points: int, *, strict: bool = False) -> bool:
    """Validate parameters for bitmask mode.

    Returns True if bitmask should be used, False if falling back to standard SR.
    Raises ValueError if algorithm is incompatible (not SR), or if strict=True
    and n_points is outside the bitmask range.
    """
    if algorithm != "sr":
        msg = "use_bitmask=True only supports algorithm='sr' (Shrake-Rupley)"
        raise ValueError(msg)
    if not (_BITMASK_MIN_N_POINTS <= n_points <= _BITMASK_MAX_N_POINTS):
        if strict:
            msg = (
                f"bitmask_correction=True requires n_points in "
                f"{_BITMASK_MIN_N_POINTS}..{_BITMASK_MAX_N_POINTS}, got {n_points}"
            )
            raise ValueError(msg)
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
    int zsasa_abi_version(void);

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

    int zsasa_calc_sr_bitmask_corrected(
        const double* x, const double* y, const double* z, const double* radii,
        size_t n_atoms, uint32_t n_points, double probe_radius, size_t n_threads,
        double correction_coeff, double* atom_areas, double* total_area
    );

    // Bitmask batch SASA calculation (f64 internal precision)
    int zsasa_calc_sr_batch_bitmask(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    int zsasa_calc_sr_batch_bitmask_corrected(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, double correction_coeff, float* atom_areas
    );

    // Bitmask batch SASA calculation (f32 internal precision)
    int zsasa_calc_sr_batch_bitmask_f32(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, float* atom_areas
    );

    int zsasa_calc_sr_batch_bitmask_f32_corrected(
        const float* coordinates, size_t n_frames, size_t n_atoms,
        const float* radii, uint32_t n_points, float probe_radius,
        size_t n_threads, double correction_coeff, float* atom_areas
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


def _library_file_name() -> str:
    """Return the platform-specific file name of the shared library."""
    if sys.platform == "darwin":
        return "libzsasa.dylib"
    if sys.platform == "win32":
        return "zsasa.dll"
    return "libzsasa.so"


def _choose_checkout_library(bundled: Path, built: Path) -> Path:
    """Choose between the library bundled in the package and a ``zig-out`` build.

    A development checkout can have both: the package directory holds a copy made
    by the build hook (untracked), ``zig-out`` holds the latest ``zig build``.
    The bundled copy goes stale when the sources change without a reinstall, and
    loading it instead of the fresh build would run old code with no warning. The
    newer file wins; a warning names the one that was skipped when the two differ.
    """
    if filecmp.cmp(bundled, built, shallow=False):
        return bundled  # the build hook's copy of this very build

    if built.stat().st_mtime > bundled.stat().st_mtime:
        chosen, skipped = built, bundled
    else:
        chosen, skipped = bundled, built
    warnings.warn(
        f"Found two different zsasa libraries: {bundled} (bundled in the package) and "
        f"{built} (build output). Loading {chosen}, which is newer; skipping {skipped}. "
        "Rebuild with 'zig build -Doptimize=ReleaseFast', reinstall the package, or delete "
        "the stale file to silence this warning. Set ZSASA_LIB to the library you want to "
        "load to override the choice.",
        stacklevel=2,
    )
    return chosen


def _find_library() -> Path:
    """Find the zsasa shared library.

    Order: the ``ZSASA_LIB`` environment variable, then the package directory and the
    ``zig-out`` directory of the checkout this file belongs to (the newer wins when both
    exist, see ``_choose_checkout_library``), then the system library directories. The
    current directory is deliberately not searched: loading a shared library from
    wherever the process happens to run would let a stray file run code in it.
    """
    # Check environment variable first
    if lib_path := os.environ.get("ZSASA_LIB"):
        return Path(lib_path)

    lib_name = _library_file_name()
    package_dir = Path(__file__).parent
    checkout_root = package_dir.parent.parent

    # Bundled in the package (wheel installation, or a copy made by the build hook)
    bundled = package_dir / lib_name
    # Development checkout: python/zsasa -> zig-out/lib or zig-out/bin (Windows DLLs)
    built = next(
        (
            path
            for path in (
                checkout_root / "zig-out" / "lib" / lib_name,
                checkout_root / "zig-out" / "bin" / lib_name,
            )
            if path.exists()
        ),
        None,
    )

    if bundled.exists() and built is not None:
        return _choose_checkout_library(bundled, built)
    if bundled.exists():
        return bundled
    if built is not None:
        return built

    for path in (Path("/usr/local/lib") / lib_name, Path("/usr/lib") / lib_name):
        if path.exists():
            return path

    msg = (
        f"Could not find {lib_name}. "
        f"Please install with: pip install zsasa "
        f"(requires Zig 0.16.0+ to be installed)"
    )
    raise FileNotFoundError(msg)


def _check_abi_version(lib: Any, lib_path: Path) -> None:
    """Refuse a library whose C ABI differs from the one this package declares.

    The signatures in ``_CDEF`` are written by hand. Calling a function through a
    declaration that no longer matches the library is undefined behavior (wrong
    arguments, corrupted memory), so a mismatch is an error, not a warning.
    """
    rebuild = (
        "Rebuild it with 'zig build -Doptimize=ReleaseFast' from the same checkout as this "
        "package, reinstall zsasa, or point ZSASA_LIB at a matching library."
    )
    try:
        found = lib.zsasa_abi_version()
    except AttributeError as e:
        msg = (
            f"The zsasa library {lib_path} does not export zsasa_abi_version(): it was built "
            f"from an older version than this Python package (expects ABI version "
            f"{_EXPECTED_ABI_VERSION}). {rebuild}"
        )
        raise ImportError(msg) from e
    if found != _EXPECTED_ABI_VERSION:
        msg = (
            f"The zsasa library {lib_path} has ABI version {found}, but this Python package "
            f"expects ABI version {_EXPECTED_ABI_VERSION}. {rebuild}"
        )
        raise ImportError(msg)


def _load_library() -> tuple[FFI, Any]:
    """Load the zsasa shared library using cffi and check its ABI version."""
    ffi = FFI()
    ffi.cdef(_CDEF)
    lib_path = _find_library()
    lib = ffi.dlopen(str(lib_path))
    _check_abi_version(lib, lib_path)
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


def _validate_frame_selection(start: int, stop: int | None, step: int) -> None:
    """Reject a ``start``/``stop``/``step`` frame selection that cannot select frames.

    Called before a trajectory is opened, so that a bad argument fails at once with a
    message about the argument instead of as a ``ZeroDivisionError`` after the file
    has been opened, or as a silent empty selection.
    """
    if step < 1:
        msg = f"step must be a positive integer (1 selects every frame), got {step}"
        raise ValueError(msg)
    if start < 0:
        msg = f"start must be non-negative, got {start}"
        raise ValueError(msg)
    if stop is not None and stop < 0:
        msg = f"stop must be non-negative or None, got {stop}"
        raise ValueError(msg)


def _raise_trajectory_open_error(kind: str, path: str, code: int) -> NoReturn:
    """Raise the exception that matches why a trajectory file could not be opened.

    ``kind`` is the format name ("XTC" or "DCD"), ``code`` the error code from the
    C open function. The C library reports a missing file and a file that is not
    valid for the format as different codes; the file system tells a directory or an
    unreadable file apart from a malformed one.
    """
    if code == ZSASA_ERROR_OUT_OF_MEMORY:
        msg = f"Out of memory opening {kind} file: {path}"
        raise MemoryError(msg)
    if code in (ZSASA_ERROR_INVALID_INPUT, ZSASA_ERROR_INVALID_FORMAT):
        target = Path(path)
        if not target.exists():
            msg = f"{kind} file not found: {path}"
            raise FileNotFoundError(msg)
        if target.is_dir():
            msg = f"{kind} path is a directory, not a file: {path}"
            raise IsADirectoryError(msg)
        if not os.access(target, os.R_OK):
            msg = f"Permission denied reading {kind} file: {path}"
            raise PermissionError(msg)
        if code == ZSASA_ERROR_INVALID_FORMAT:
            msg = f"{path} is not a valid {kind} file: it is empty, truncated, or in another format"
            raise ValueError(msg)
        msg = f"Cannot open {kind} file: {path}"
        raise OSError(msg)
    msg = f"Error opening {kind} file {path}: error code {code}"
    raise RuntimeError(msg)


def _raise_trajectory_read_error(kind: str, code: int) -> NoReturn:
    """Raise the exception for a failed ``zsasa_*_read_frame`` call."""
    if code == ZSASA_ERROR_OUT_OF_MEMORY:
        msg = f"Out of memory reading {kind} frame"
        raise MemoryError(msg)
    if code == ZSASA_ERROR_INVALID_FORMAT:
        msg = f"Error reading {kind} frame: the file is corrupt or truncated"
        raise RuntimeError(msg)
    msg = f"Error reading {kind} frame: error code {code}"
    raise RuntimeError(msg)

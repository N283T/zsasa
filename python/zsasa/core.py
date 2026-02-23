"""Core cffi bindings for zsasa library."""

from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Any, Literal

import numpy as np
from cffi import FFI
from numpy.typing import NDArray

# Error codes from C API
ZSASA_OK = 0
ZSASA_ERROR_INVALID_INPUT = -1
ZSASA_ERROR_OUT_OF_MEMORY = -2
ZSASA_ERROR_CALCULATION = -3
ZSASA_ERROR_FILE_IO = -4

# Algorithm constants
ZSASA_ALGORITHM_SR = 0
ZSASA_ALGORITHM_LR = 1

# Classifier types
ZSASA_CLASSIFIER_NACCESS = 0
ZSASA_CLASSIFIER_PROTOR = 1
ZSASA_CLASSIFIER_OONS = 2

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
        f"(requires Zig 0.15.2+ to be installed)"
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
        coords: Atom coordinates as (N, 3) array.
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
    if algorithm == "sr":
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
    elif result != ZSASA_OK:
        msg = f"Unknown error code: {result}"
        raise RuntimeError(msg)

    return SasaResult(
        total_area=total_area[0],
        atom_areas=atom_areas,
    )


# =============================================================================
# Classifier Types and Classes
# =============================================================================


class ClassifierType(IntEnum):
    """Available classifier types for atom radius assignment.

    Attributes:
        NACCESS: NACCESS-compatible radii (default, most commonly used).
        PROTOR: ProtOr radii based on hybridization state.
        OONS: Ooi, Oobatake, Nemethy, Scheraga radii (older FreeSASA default).
    """

    NACCESS = ZSASA_CLASSIFIER_NACCESS
    PROTOR = ZSASA_CLASSIFIER_PROTOR
    OONS = ZSASA_CLASSIFIER_OONS


class AtomClass(IntEnum):
    """Atom polarity classes.

    Attributes:
        POLAR: Polar atoms (e.g., N, O).
        APOLAR: Apolar/hydrophobic atoms (e.g., C).
        UNKNOWN: Unknown classification.
    """

    POLAR = ZSASA_ATOM_CLASS_POLAR
    APOLAR = ZSASA_ATOM_CLASS_APOLAR
    UNKNOWN = ZSASA_ATOM_CLASS_UNKNOWN


# =============================================================================
# Classifier Functions
# =============================================================================


def get_radius(
    residue: str,
    atom: str,
    classifier_type: ClassifierType = ClassifierType.NACCESS,
) -> float | None:
    """Get van der Waals radius for an atom using the specified classifier.

    Args:
        residue: Residue name (e.g., "ALA", "GLY").
        atom: Atom name (e.g., "CA", "CB").
        classifier_type: Classifier to use. Default: ClassifierType.NACCESS.

    Returns:
        Radius in Angstroms, or None if atom is not found in classifier.

    Example:
        >>> from zsasa import get_radius, ClassifierType
        >>> get_radius("ALA", "CA")
        1.87
        >>> get_radius("ALA", "XX")  # Unknown atom
        None
    """
    _, lib = _get_lib()
    radius = lib.zsasa_classifier_get_radius(
        classifier_type,
        residue.encode("utf-8"),
        atom.encode("utf-8"),
    )
    if np.isnan(radius):
        return None
    return radius


def get_atom_class(
    residue: str,
    atom: str,
    classifier_type: ClassifierType = ClassifierType.NACCESS,
) -> int:
    """Get atom polarity class using the specified classifier.

    Args:
        residue: Residue name (e.g., "ALA", "GLY").
        atom: Atom name (e.g., "CA", "CB").
        classifier_type: Classifier to use. Default: ClassifierType.NACCESS.

    Returns:
        AtomClass constant (POLAR, APOLAR, or UNKNOWN).

    Example:
        >>> from zsasa import get_atom_class, AtomClass
        >>> get_atom_class("ALA", "CA") == AtomClass.APOLAR
        True
        >>> get_atom_class("ALA", "O") == AtomClass.POLAR
        True
    """
    _, lib = _get_lib()
    return lib.zsasa_classifier_get_class(
        classifier_type,
        residue.encode("utf-8"),
        atom.encode("utf-8"),
    )


def guess_radius(element: str) -> float | None:
    """Guess van der Waals radius from element symbol.

    Args:
        element: Element symbol (e.g., "C", "N", "FE").
                 Case-insensitive, whitespace is trimmed.

    Returns:
        Radius in Angstroms, or None if element is not recognized.

    Example:
        >>> from zsasa import guess_radius
        >>> guess_radius("C")
        1.7
        >>> guess_radius("FE")
        1.26
        >>> guess_radius("XX")  # Unknown
        None
    """
    _, lib = _get_lib()
    radius = lib.zsasa_guess_radius(element.encode("utf-8"))
    if np.isnan(radius):
        return None
    return radius


def guess_radius_from_atom_name(atom_name: str) -> float | None:
    """Guess van der Waals radius from PDB-style atom name.

    Extracts element symbol from atom name following PDB conventions:
    - Leading space indicates single-char element (e.g., " CA " = Carbon alpha)
    - No leading space may indicate 2-char element (e.g., "FE  " = Iron)

    Args:
        atom_name: PDB-style atom name (e.g., " CA ", "FE  ").

    Returns:
        Radius in Angstroms, or None if element cannot be determined.

    Example:
        >>> from zsasa import guess_radius_from_atom_name
        >>> guess_radius_from_atom_name(" CA ")  # Carbon alpha
        1.7
        >>> guess_radius_from_atom_name("FE  ")  # Iron
        1.26
    """
    _, lib = _get_lib()
    radius = lib.zsasa_guess_radius_from_atom_name(atom_name.encode("utf-8"))
    if np.isnan(radius):
        return None
    return radius


@dataclass
class ClassificationResult:
    """Result of batch atom classification.

    Attributes:
        radii: Per-atom van der Waals radii in Angstroms.
               NaN values indicate atoms not found in the classifier.
               Use np.isnan(result.radii) to find unknown atoms.
        classes: Per-atom polarity classes (AtomClass constants).

    Example:
        >>> result = classify_atoms(residues, atoms)
        >>> unknown_mask = np.isnan(result.radii)
        >>> if unknown_mask.any():
        ...     print(f"Found {unknown_mask.sum()} unknown atoms")
    """

    radii: NDArray[np.float64]
    classes: NDArray[np.int32]

    def __repr__(self) -> str:
        return f"ClassificationResult(n_atoms={len(self.radii)})"


def classify_atoms(
    residues: list[str],
    atoms: list[str],
    classifier_type: ClassifierType = ClassifierType.NACCESS,
    *,
    include_classes: bool = True,
) -> ClassificationResult:
    """Classify multiple atoms at once (batch operation).

    This is more efficient than calling get_radius for each atom individually.

    Args:
        residues: List of residue names.
        atoms: List of atom names (must be same length as residues).
        classifier_type: Classifier to use. Default: ClassifierType.NACCESS.
        include_classes: Whether to compute atom classes. Default: True.

    Returns:
        ClassificationResult with radii and classes arrays.
        Unknown atoms have NaN radius and UNKNOWN class.

    Raises:
        ValueError: If residues and atoms have different lengths.

    Example:
        >>> from zsasa import classify_atoms
        >>> result = classify_atoms(
        ...     ["ALA", "ALA", "GLY"],
        ...     ["CA", "O", "N"],
        ... )
        >>> result.radii
        array([1.87, 1.4 , 1.65])
    """
    ffi, lib = _get_lib()

    if len(residues) != len(atoms):
        msg = f"residues and atoms must have same length: {len(residues)} != {len(atoms)}"
        raise ValueError(msg)

    n_atoms = len(residues)
    if n_atoms == 0:
        return ClassificationResult(
            radii=np.array([], dtype=np.float64),
            classes=np.array([], dtype=np.int32),
        )

    # Encode strings and create cffi arrays
    # Keep references to prevent garbage collection
    residues_bytes = [ffi.new("char[]", r.encode("utf-8")) for r in residues]
    atoms_bytes = [ffi.new("char[]", a.encode("utf-8")) for a in atoms]

    residues_arr = ffi.new("char*[]", residues_bytes)
    atoms_arr = ffi.new("char*[]", atoms_bytes)

    # Allocate output arrays
    radii = np.zeros(n_atoms, dtype=np.float64)
    classes = np.zeros(n_atoms, dtype=np.int32) if include_classes else None

    radii_ptr = ffi.cast("double*", radii.ctypes.data)
    classes_ptr = ffi.cast("int*", classes.ctypes.data) if include_classes else ffi.NULL

    result = lib.zsasa_classify_atoms(
        classifier_type,
        residues_arr,
        atoms_arr,
        n_atoms,
        radii_ptr,
        classes_ptr,
    )

    if result == ZSASA_ERROR_INVALID_INPUT:
        msg = f"Invalid classifier type: {classifier_type}"
        raise ValueError(msg)
    elif result != ZSASA_OK:
        msg = f"Classification error: {result}"
        raise RuntimeError(msg)

    if not include_classes:
        classes = np.full(n_atoms, AtomClass.UNKNOWN, dtype=np.int32)

    return ClassificationResult(radii=radii, classes=classes)


# =============================================================================
# RSA (Relative Solvent Accessibility) Functions
# =============================================================================

# Maximum SASA values for standard amino acids (Tien et al. 2013)
MAX_SASA = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLN": 225.0,
    "GLU": 223.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}


def get_max_sasa(residue_name: str) -> float | None:
    """Get maximum SASA value for a standard amino acid.

    Values from Tien et al. (2013) "Maximum allowed solvent accessibilities
    of residues in proteins".

    Args:
        residue_name: 3-letter residue code (e.g., "ALA", "GLY").

    Returns:
        Maximum SASA in Angstroms², or None if residue is not a standard amino acid.

    Example:
        >>> from zsasa import get_max_sasa
        >>> get_max_sasa("ALA")
        129.0
        >>> get_max_sasa("TRP")
        285.0
        >>> get_max_sasa("HOH")  # Water - not a standard amino acid
        None
    """
    _, lib = _get_lib()
    max_sasa = lib.zsasa_get_max_sasa(residue_name.encode("utf-8"))
    if np.isnan(max_sasa):
        return None
    return max_sasa


def calculate_rsa(sasa: float, residue_name: str) -> float | None:
    """Calculate RSA (Relative Solvent Accessibility) for a single residue.

    RSA = SASA / MaxSASA

    Args:
        sasa: Observed SASA value in Angstroms².
        residue_name: 3-letter residue code (e.g., "ALA", "GLY").

    Returns:
        RSA value (typically 0.0-1.0), or None if residue is not a standard amino acid.
        Note: RSA > 1.0 is possible for exposed terminal residues.

    Example:
        >>> from zsasa import calculate_rsa
        >>> calculate_rsa(64.5, "ALA")  # 64.5 / 129.0 = 0.5
        0.5
        >>> calculate_rsa(150.0, "GLY")  # RSA > 1.0 is possible
        1.4423076923076923
    """
    _, lib = _get_lib()
    rsa = lib.zsasa_calculate_rsa(sasa, residue_name.encode("utf-8"))
    if np.isnan(rsa):
        return None
    return rsa


def calculate_rsa_batch(
    sasas: NDArray[np.float64] | list[float],
    residue_names: list[str],
) -> NDArray[np.float64]:
    """Calculate RSA for multiple residues at once (batch operation).

    This is more efficient than calling calculate_rsa for each residue individually.

    Args:
        sasas: Array of SASA values in Angstroms².
        residue_names: List of 3-letter residue codes (must be same length as sasas).

    Returns:
        Array of RSA values. NaN values indicate non-standard amino acids.

    Raises:
        ValueError: If sasas and residue_names have different lengths.

    Example:
        >>> import numpy as np
        >>> from zsasa import calculate_rsa_batch
        >>> sasas = np.array([64.5, 52.0, 100.0])
        >>> residues = ["ALA", "GLY", "HOH"]  # HOH is not standard
        >>> rsa = calculate_rsa_batch(sasas, residues)
        >>> rsa
        array([0.5       , 0.5       ,        nan])
    """
    ffi, lib = _get_lib()

    sasas = np.ascontiguousarray(sasas, dtype=np.float64)
    n_residues = len(sasas)

    if len(residue_names) != n_residues:
        msg = f"sasas and residue_names must have same length: {n_residues} != {len(residue_names)}"
        raise ValueError(msg)

    if n_residues == 0:
        return np.array([], dtype=np.float64)

    # Encode strings and create cffi array
    # Keep references to prevent garbage collection
    residues_bytes = [ffi.new("char[]", r.encode("utf-8")) for r in residue_names]
    residues_arr = ffi.new("char*[]", residues_bytes)

    # Allocate output array
    rsa_out = np.zeros(n_residues, dtype=np.float64)

    sasas_ptr = ffi.cast("double*", sasas.ctypes.data)
    rsa_ptr = ffi.cast("double*", rsa_out.ctypes.data)

    result = lib.zsasa_calculate_rsa_batch(sasas_ptr, residues_arr, n_residues, rsa_ptr)
    if result != ZSASA_OK:
        msg = f"RSA batch calculation failed with error code: {result}"
        raise RuntimeError(msg)

    return rsa_out


# =============================================================================
# Batch SASA Calculation (for MD trajectories)
# =============================================================================


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

    # Call the appropriate batch function based on algorithm and precision
    if precision == "f64":
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
    elif result != ZSASA_OK:
        msg = f"Unknown error code: {result}"
        raise RuntimeError(msg)

    return BatchSasaResult(atom_areas=atom_areas)


# =============================================================================
# Batch Directory Processing
# =============================================================================


@dataclass
class BatchDirResult:
    """Result of directory batch processing.

    Attributes:
        total_files: Number of supported structure files found in the directory.
        successful: Number of successfully processed files.
        failed: Number of files that failed processing.
        filenames: List of filenames (file name only, without directory path).
        n_atoms: Per-file atom counts.
        total_sasa: Per-file total SASA in Angstroms² (NaN for failed files).
        status: Per-file status (1=ok, 0=failed).
    """

    total_files: int
    successful: int
    failed: int
    filenames: list[str]
    n_atoms: list[int]
    total_sasa: list[float]
    status: list[int]

    def __post_init__(self) -> None:
        n = len(self.filenames)
        if len(self.n_atoms) != n or len(self.total_sasa) != n or len(self.status) != n:
            msg = (
                f"All per-file lists must have the same length as filenames ({n}), "
                f"got n_atoms={len(self.n_atoms)}, total_sasa={len(self.total_sasa)}, "
                f"status={len(self.status)}"
            )
            raise ValueError(msg)
        if self.total_files != n:
            msg = f"total_files ({self.total_files}) != len(filenames) ({n})"
            raise ValueError(msg)
        if self.successful + self.failed != self.total_files:
            msg = (
                f"successful ({self.successful}) + failed ({self.failed}) "
                f"!= total_files ({self.total_files})"
            )
            raise ValueError(msg)

    def __repr__(self) -> str:
        return (
            f"BatchDirResult(total_files={self.total_files}, "
            f"successful={self.successful}, failed={self.failed})"
        )


def process_directory(
    input_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    classifier: ClassifierType | None = ClassifierType.PROTOR,
    include_hydrogens: bool = False,
    include_hetatm: bool = False,
) -> BatchDirResult:
    """Process all supported structure files in a directory for SASA calculation.

    Supported formats: PDB (.pdb), mmCIF (.cif, .mmcif), PDB/ENT (.ent),
    JSON (.json), and their gzip-compressed variants (.gz).

    Args:
        input_dir: Path to directory containing structure files.
        output_dir: Optional path for per-file output. None = no file output.
        algorithm: Algorithm to use: "sr" (Shrake-Rupley) or "lr" (Lee-Richards).
        n_points: Number of test points per atom (SR only; ignored for LR).
            Default: 100.
        n_slices: Number of slices per atom (LR only; ignored for SR).
            Default: 20.
        probe_radius: Water probe radius in Angstroms. Default: 1.4.
        n_threads: Number of threads to use. 0 = auto-detect. Default: 0.
        classifier: Classifier for radius assignment. None = use input radii.
            Default: ClassifierType.PROTOR.
        include_hydrogens: Whether to include hydrogen atoms. Default: False.
        include_hetatm: Whether to include HETATM records. Default: False.

    Returns:
        BatchDirResult with per-file details.

    Raises:
        ValueError: If input parameters are invalid.
        FileNotFoundError: If the input directory does not exist.
        MemoryError: If out of memory.
        RuntimeError: For other processing errors.

    Example:
        >>> from zsasa import process_directory
        >>> result = process_directory("path/to/pdbs/")
        >>> print(f"Processed {result.successful}/{result.total_files} files")
    """
    ffi, lib = _get_lib()

    # Validate algorithm
    if algorithm == "sr":
        algo_int = ZSASA_ALGORITHM_SR
        n_points_val = n_points
    elif algorithm == "lr":
        algo_int = ZSASA_ALGORITHM_LR
        n_points_val = n_slices
    else:
        msg = f"Unknown algorithm: {algorithm}. Use 'sr' or 'lr'."
        raise ValueError(msg)

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

    # Map classifier
    classifier_int = int(classifier) if classifier is not None else -1

    # Convert paths
    input_dir_bytes = str(input_dir).encode("utf-8")
    output_dir_bytes = str(output_dir).encode("utf-8") if output_dir is not None else ffi.NULL

    error_code = ffi.new("int*")

    handle = lib.zsasa_batch_dir_process(
        input_dir_bytes,
        output_dir_bytes,
        algo_int,
        n_points_val,
        probe_radius,
        n_threads,
        classifier_int,
        int(include_hydrogens),
        int(include_hetatm),
        error_code,
    )

    if handle == ffi.NULL:
        ec = error_code[0]
        if ec == ZSASA_ERROR_INVALID_INPUT:
            msg = "Invalid input parameters for directory batch processing"
            raise ValueError(msg)
        elif ec == ZSASA_ERROR_OUT_OF_MEMORY:
            msg = "Out of memory during directory batch processing"
            raise MemoryError(msg)
        elif ec == ZSASA_ERROR_CALCULATION:
            msg = "SASA calculation failed during directory batch processing"
            raise RuntimeError(msg)
        elif ec == ZSASA_ERROR_FILE_IO:
            msg = f"Directory not found or not readable: {input_dir}"
            raise FileNotFoundError(msg)
        else:
            msg = f"Directory batch processing failed with error code: {ec}"
            raise RuntimeError(msg)

    try:
        total_files = lib.zsasa_batch_dir_get_total_files(handle)
        successful = lib.zsasa_batch_dir_get_successful(handle)
        failed = lib.zsasa_batch_dir_get_failed(handle)

        filenames: list[str] = []
        n_atoms_list: list[int] = []
        total_sasa_list: list[float] = []
        status_list: list[int] = []

        for i in range(total_files):
            fname_ptr = lib.zsasa_batch_dir_get_filename(handle, i)
            if fname_ptr == ffi.NULL:
                msg = (
                    f"Internal error: zsasa_batch_dir_get_filename returned NULL "
                    f"for index {i} (total_files={total_files})"
                )
                raise RuntimeError(msg)
            filenames.append(ffi.string(fname_ptr).decode("utf-8"))
            n_atoms_list.append(lib.zsasa_batch_dir_get_n_atoms(handle, i))
            sasa = lib.zsasa_batch_dir_get_total_sasa(handle, i)
            total_sasa_list.append(float("nan") if math.isnan(sasa) else sasa)
            st = lib.zsasa_batch_dir_get_status(handle, i)
            if st not in (0, 1):
                msg = f"Internal error: unexpected status {st} for file index {i} (expected 0 or 1)"
                raise RuntimeError(msg)
            status_list.append(st)

        return BatchDirResult(
            total_files=total_files,
            successful=successful,
            failed=failed,
            filenames=filenames,
            n_atoms=n_atoms_list,
            total_sasa=total_sasa_list,
            status=status_list,
        )
    finally:
        lib.zsasa_batch_dir_free(handle)

"""RSA (Relative Solvent Accessibility) functions."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from zsasa._ffi import ZSASA_OK, _get_lib

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

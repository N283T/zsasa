"""Batch directory processing for SASA calculation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from zsasa._ffi import (
    ZSASA_ALGORITHM_LR,
    ZSASA_ALGORITHM_SR,
    ZSASA_ERROR_CALCULATION,
    ZSASA_ERROR_FILE_IO,
    ZSASA_ERROR_INVALID_INPUT,
    ZSASA_ERROR_OUT_OF_MEMORY,
    _get_lib,
)
from zsasa.classifier import ClassifierType


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

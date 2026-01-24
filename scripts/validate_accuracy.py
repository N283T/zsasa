#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = []
# ///
"""Validate Zig SASA implementation against FreeSASA reference values.

Compares the Zig implementation output with FreeSASA reference data
and reports accuracy metrics.

Usage:
    ./validate_accuracy.py [--algorithm=sr|lr] [--tolerance=2.0]

Examples:
    ./validate_accuracy.py                    # Validate with SR algorithm
    ./validate_accuracy.py --algorithm=lr     # Validate with LR algorithm
    ./validate_accuracy.py --tolerance=1.0    # Require <1% difference
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ValidationResult:
    """Result of validating one structure."""

    pdb_id: str
    n_atoms: int
    reference_sasa: float
    zig_sasa: float
    difference_percent: float
    passed: bool
    time_ms: float | None = None


def run_zig_sasa(
    input_path: Path,
    algorithm: str = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    classifier: str | None = None,
) -> tuple[float, float]:
    """Run Zig SASA implementation and return (total_area, time_ms)."""
    zig_binary = Path(__file__).parent.parent / "zig-out" / "bin" / "freesasa_zig"

    if not zig_binary.exists():
        raise FileNotFoundError(f"Zig binary not found: {zig_binary}")

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
        output_path = Path(f.name)

    try:
        cmd = [
            str(zig_binary),
            f"--algorithm={algorithm}",
        ]
        if algorithm == "sr":
            cmd.append(f"--n-points={n_points}")
        else:
            cmd.append(f"--n-slices={n_slices}")

        if classifier:
            cmd.append(f"--classifier={classifier}")

        cmd.extend([str(input_path), str(output_path)])

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            raise RuntimeError(f"Zig failed: {result.stderr}")

        # Parse time from stderr (e.g., "Total area: 18908.90 Å²")
        time_ms = None
        for line in result.stderr.split("\n"):
            if "ms" in line.lower():
                # Try to extract time
                match = re.search(r"(\d+(?:\.\d+)?)\s*ms", line, re.IGNORECASE)
                if match:
                    time_ms = float(match.group(1))

        with open(output_path) as f:
            data = json.load(f)

        return data["total_area"], time_ms

    finally:
        output_path.unlink(missing_ok=True)


def validate_structure(
    pdb_id: str,
    inputs_dir: Path,
    references_dir: Path,
    algorithm: str,
    tolerance: float,
    classifier: str | None = None,
) -> ValidationResult | None:
    """Validate one structure against reference."""
    input_path = inputs_dir / f"{pdb_id}.json"
    reference_path = references_dir / f"{pdb_id}_n100_p1.4.json"

    if not input_path.exists() or not reference_path.exists():
        return None

    with open(reference_path) as f:
        ref = json.load(f)

    reference_sasa = ref["total_area"]
    n_atoms = ref["n_atoms"]

    try:
        zig_sasa, time_ms = run_zig_sasa(input_path, algorithm, classifier=classifier)
    except Exception as e:
        print(f"  {pdb_id}: ERROR - {e}")
        return None

    if reference_sasa == 0:
        diff_percent = 0.0 if zig_sasa == 0 else float("inf")
    else:
        diff_percent = abs(zig_sasa - reference_sasa) / reference_sasa * 100
    passed = diff_percent <= tolerance

    return ValidationResult(
        pdb_id=pdb_id,
        n_atoms=n_atoms,
        reference_sasa=reference_sasa,
        zig_sasa=zig_sasa,
        difference_percent=diff_percent,
        passed=passed,
        time_ms=time_ms,
    )


def main() -> int:
    # Parse arguments
    algorithm = "sr"
    tolerance = 2.0
    classifier = "protor"  # Default: use ProtOr for fair comparison with FreeSASA

    for arg in sys.argv[1:]:
        if arg.startswith("--algorithm="):
            algorithm = arg.split("=")[1]
        elif arg.startswith("--tolerance="):
            tolerance = float(arg.split("=")[1])
        elif arg.startswith("--classifier="):
            classifier = arg.split("=")[1]
        elif arg == "--no-classifier":
            classifier = None

    # Setup paths
    base_dir = Path(__file__).parent.parent / "benchmarks"
    inputs_dir = base_dir / "inputs"
    references_dir = base_dir / "references"

    structures = ["1crn", "1ubq", "1a0q", "3hhb", "1aon", "4v6x"]

    print("=" * 70)
    clf_str = classifier if classifier else "none (element-based)"
    print(
        f"SASA Validation (algorithm={algorithm}, classifier={clf_str}, tolerance={tolerance}%)"
    )
    print("=" * 70)
    print(
        f"{'PDB':<8} {'Atoms':>8} {'FreeSASA':>12} {'Zig':>12} {'Diff%':>8} {'Status':<8}"
    )
    print("-" * 70)

    results = []
    for pdb_id in structures:
        result = validate_structure(
            pdb_id, inputs_dir, references_dir, algorithm, tolerance, classifier
        )
        if result:
            results.append(result)
            status = "PASS" if result.passed else "FAIL"
            print(
                f"{result.pdb_id:<8} {result.n_atoms:>8} "
                f"{result.reference_sasa:>12.2f} {result.zig_sasa:>12.2f} "
                f"{result.difference_percent:>7.3f}% {status:<8}"
            )
        else:
            print(f"{pdb_id:<8} {'SKIPPED':>8}")

    print("-" * 70)

    # Summary
    if results:
        passed = sum(1 for r in results if r.passed)
        total = len(results)
        avg_diff = sum(r.difference_percent for r in results) / len(results)

        print(f"\nSummary: {passed}/{total} passed (avg diff: {avg_diff:.3f}%)")

        if passed == total:
            print("All validations PASSED!")
            return 0
        else:
            print("Some validations FAILED!")
            return 1

    return 1


if __name__ == "__main__":
    sys.exit(main())

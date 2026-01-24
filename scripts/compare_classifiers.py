#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
# ]
# ///
"""Print FreeSASA default classifier radii for reference.

This script prints the atom radii returned by FreeSASA's Python bindings
for all standard amino acid atoms. FreeSASA Python uses ProtOr radii by default.

To compare with our Zig NACCESS implementation, run `zig build test` and
check the test assertions in src/classifier_naccess.zig.

Usage:
    ./compare_classifiers.py

Note: The radii shown are ProtOr values, not NACCESS. For NACCESS comparison,
refer to FreeSASA's share/naccess.config file.
"""

from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass


# Standard amino acid atoms to test
AMINO_ACIDS = {
    "ALA": ["N", "CA", "C", "O", "CB"],
    "ARG": ["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"],
    "ASP": ["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"],
    "CYS": ["N", "CA", "C", "O", "CB", "SG"],
    "GLN": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"],
    "GLY": ["N", "CA", "C", "O"],
    "HIS": ["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"],
    "LEU": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"],
    "LYS": ["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"],
    "MET": ["N", "CA", "C", "O", "CB", "CG", "SD", "CE"],
    "PHE": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["N", "CA", "C", "O", "CB", "CG", "CD"],
    "SER": ["N", "CA", "C", "O", "CB", "OG"],
    "THR": ["N", "CA", "C", "O", "CB", "OG1", "CG2"],
    "TRP": [
        "N",
        "CA",
        "C",
        "O",
        "CB",
        "CG",
        "CD1",
        "CD2",
        "NE1",
        "CE2",
        "CE3",
        "CZ2",
        "CZ3",
        "CH2",
    ],
    "TYR": ["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["N", "CA", "C", "O", "CB", "CG1", "CG2"],
}


@dataclass
class ComparisonResult:
    residue: str
    atom: str
    freesasa_radius: float | None
    zig_radius: float | None
    match: bool
    diff: float | None


def get_freesasa_radius(residue: str, atom: str) -> float | None:
    """Get radius from FreeSASA's Python bindings.

    Note: FreeSASA Python module uses ProtOr radii by default,
    regardless of any classifier parameter. To compare with NACCESS,
    refer to FreeSASA's naccess.config file directly.
    """
    import freesasa

    try:
        clf = freesasa.Classifier()
        return clf.radius(residue, atom)
    except Exception:
        return None


def compare_classifiers() -> list[ComparisonResult]:
    """Compare radii between FreeSASA (ProtOr default) and Zig for all amino acid atoms."""
    results = []

    for residue, atoms in AMINO_ACIDS.items():
        for atom in atoms:
            fs_radius = get_freesasa_radius(residue, atom)

            results.append(
                ComparisonResult(
                    residue=residue,
                    atom=atom,
                    freesasa_radius=fs_radius,
                    zig_radius=None,  # Compare with Zig test output manually
                    match=True,  # Assume match for now
                    diff=None,
                )
            )

    return results


def print_freesasa_radii() -> None:
    """Print all FreeSASA radii for verification.

    Note: FreeSASA Python uses ProtOr radii by default.
    For NACCESS comparison, refer to FreeSASA's naccess.config file.
    """
    print("FreeSASA Default (ProtOr) Classifier Radii")
    print("=" * 60)
    print(f"{'Residue':<8} {'Atom':<6} {'Radius (Å)':<12}")
    print("-" * 60)

    for residue, atoms in sorted(AMINO_ACIDS.items()):
        for atom in atoms:
            radius = get_freesasa_radius(residue, atom)
            if radius is not None:
                print(f"{residue:<8} {atom:<6} {radius:<12.2f}")
            else:
                print(f"{residue:<8} {atom:<6} {'N/A':<12}")
        print()


def main() -> int:
    print_freesasa_radii()

    # Summary
    print("\nSummary:")
    print("Note: FreeSASA Python uses ProtOr radii by default")
    print(f"Amino acids tested: {len(AMINO_ACIDS)}")
    total_atoms = sum(len(atoms) for atoms in AMINO_ACIDS.values())
    print(f"Total atoms: {total_atoms}")
    print("\nTo compare with Zig NACCESS implementation:")
    print("  Run: zig build test")
    print("  Check: src/classifier_naccess.zig test assertions")

    return 0


if __name__ == "__main__":
    sys.exit(main())

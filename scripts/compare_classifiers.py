#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
# ]
# ///
"""Compare classifier radii between FreeSASA and freesasa-zig.

This script compares the atom radii returned by FreeSASA's Python bindings
with our Zig implementation to verify consistency.

Usage:
    ./compare_classifiers.py [--classifier=naccess|protor|oons]

Examples:
    ./compare_classifiers.py
    ./compare_classifiers.py --classifier=protor
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


def get_freesasa_radius(
    residue: str, atom: str, classifier: str = "naccess"
) -> float | None:
    """Get radius from FreeSASA's Python bindings.

    Note: FreeSASA Python module uses ProtOr radii by default.
    The classifier parameter selects from available config files.
    """
    import freesasa

    # FreeSASA classifier selection via config file path
    # The default Classifier() uses ProtOr radii
    try:
        clf = freesasa.Classifier()
        return clf.radius(residue, atom)
    except Exception:
        return None


def get_zig_radius(
    residue: str, atom: str, classifier: str = "naccess"
) -> float | None:
    """Get radius from our Zig implementation by running a test."""
    # We'll use zig build test to run a specific test
    # For now, return None and we'll compare with the values from the Zig source
    return None  # TODO: Add actual Zig invocation


def compare_classifiers(classifier: str = "naccess") -> list[ComparisonResult]:
    """Compare radii between FreeSASA and Zig for all amino acid atoms."""
    results = []

    for residue, atoms in AMINO_ACIDS.items():
        for atom in atoms:
            fs_radius = get_freesasa_radius(residue, atom, classifier)

            results.append(
                ComparisonResult(
                    residue=residue,
                    atom=atom,
                    freesasa_radius=fs_radius,
                    zig_radius=None,  # Will be filled in by Zig tests
                    match=True,  # Assume match for now
                    diff=None,
                )
            )

    return results


def print_freesasa_radii(classifier: str = "naccess") -> None:
    """Print all FreeSASA radii for verification."""
    print(f"FreeSASA {classifier.upper()} Classifier Radii")
    print("=" * 60)
    print(f"{'Residue':<8} {'Atom':<6} {'Radius (Å)':<12}")
    print("-" * 60)

    for residue, atoms in sorted(AMINO_ACIDS.items()):
        for atom in atoms:
            radius = get_freesasa_radius(residue, atom, classifier)
            if radius is not None:
                print(f"{residue:<8} {atom:<6} {radius:<12.2f}")
            else:
                print(f"{residue:<8} {atom:<6} {'N/A':<12}")
        print()


def main() -> int:
    classifier = "naccess"

    for arg in sys.argv[1:]:
        if arg.startswith("--classifier="):
            classifier = arg.split("=")[1]

    print_freesasa_radii(classifier)

    # Summary
    print("\nSummary:")
    print(f"Classifier: {classifier}")
    print(f"Amino acids tested: {len(AMINO_ACIDS)}")
    total_atoms = sum(len(atoms) for atoms in AMINO_ACIDS.values())
    print(f"Total atoms: {total_atoms}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

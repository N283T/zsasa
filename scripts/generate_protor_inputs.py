#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
# ]
# ///
"""Generate input JSON files with ProtOr radii pre-applied.

This creates inputs that can be used directly with the C benchmark
without needing a classifier.

Usage:
    ./generate_protor_inputs.py
"""

from __future__ import annotations

import json
from pathlib import Path

import freesasa
import gemmi


def get_protor_radius(residue: str, atom_name: str) -> float:
    """Get ProtOr radius for an atom using FreeSASA's classifier."""
    # FreeSASA ProtOr classifier
    classifier = freesasa.Classifier.getStandardClassifier("protor")

    # Try to get radius from classifier
    radius = classifier.radius(residue, atom_name)

    # If not found, fall back to element-based
    if radius <= 0:
        # Common element radii
        element = atom_name[0] if atom_name else "C"
        fallback = {
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "S": 1.80,
            "H": 1.20,
            "P": 1.80,
        }
        radius = fallback.get(element, 1.70)

    return radius


def convert_input_to_protor(input_path: Path, output_path: Path) -> None:
    """Convert input JSON to use ProtOr radii."""
    with open(input_path) as f:
        data = json.load(f)

    n_atoms = len(data["x"])
    new_radii = []

    for i in range(n_atoms):
        residue = data["residue"][i]
        atom_name = data["atom_name"][i]
        radius = get_protor_radius(residue, atom_name)
        new_radii.append(radius)

    data["r"] = new_radii

    with open(output_path, "w") as f:
        json.dump(data, f)

    print(f"  Converted {input_path.name} -> {output_path.name} ({n_atoms} atoms)")


def main() -> int:
    base_dir = Path(__file__).parent.parent / "benchmarks"
    inputs_dir = base_dir / "inputs"
    protor_dir = base_dir / "inputs_protor"

    protor_dir.mkdir(exist_ok=True)

    print("Generating ProtOr-radii inputs...")

    for input_file in sorted(inputs_dir.glob("*.json")):
        output_file = protor_dir / input_file.name
        convert_input_to_protor(input_file, output_file)

    print(f"\nDone! ProtOr inputs saved to: {protor_dir}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

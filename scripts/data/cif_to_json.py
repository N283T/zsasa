#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "gemmi>=0.7.4",
#     "pydantic>=2.0",
# ]
# ///
"""Convert coordinate files (mmCIF/PDB) to SASA input JSON format.

Usage:
    ./cif_to_input_json.py <input_file> [output_file]

Examples:
    ./cif_to_input_json.py examples/1A0Q.cif.gz
    ./cif_to_input_json.py 1abc.pdb output.json
"""

from __future__ import annotations

import sys
from pathlib import Path

import gemmi
from pydantic import BaseModel


class AtomCoordinates(BaseModel):
    """Atom coordinates and radii for SASA calculation."""

    x: list[float]
    y: list[float]
    z: list[float]
    r: list[float]
    residue: list[str] | None = None
    atom_name: list[str] | None = None
    element: list[int] | None = None  # Atomic numbers (e.g., 6=C, 7=N, 8=O)

    @property
    def n_atoms(self) -> int:
        return len(self.x)


def extract_atoms(path: str | Path) -> AtomCoordinates:
    """Extract atom coordinates and van der Waals radii from a coordinate file.

    Uses gemmi's built-in vdw_r for radii (based on J. Phys. Chem. A 2009).

    Args:
        path: Path to mmCIF or PDB file (gzipped OK)

    Returns:
        AtomCoordinates with x, y, z coordinates, r radii, residue/atom names, and elements
    """
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    rs: list[float] = []
    residues: list[str] = []
    atom_names: list[str] = []
    elements: list[int] = []

    # Use first model only (common for X-ray structures)
    model = st[0]

    for cra in model.all():
        xs.append(cra.atom.pos.x)
        ys.append(cra.atom.pos.y)
        zs.append(cra.atom.pos.z)
        rs.append(cra.atom.element.vdw_r)
        residues.append(cra.residue.name)
        atom_names.append(cra.atom.name)
        elements.append(cra.atom.element.atomic_number)

    return AtomCoordinates(
        x=xs, y=ys, z=zs, r=rs, residue=residues, atom_name=atom_names, element=elements
    )


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 1

    input_path = Path(sys.argv[1])
    if not input_path.exists():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        return 1

    # Default output: same name with .json extension
    if len(sys.argv) >= 3:
        output_path = Path(sys.argv[2])
    else:
        # Remove .gz and .cif/.pdb extensions, add .json
        stem = input_path.name
        for suffix in [".gz", ".cif", ".pdb", ".ent"]:
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
        output_path = input_path.parent / f"{stem}.json"

    result = extract_atoms(input_path)
    print(f"Extracted {result.n_atoms} atoms from {input_path.name}")

    output_path.write_text(result.model_dump_json())

    print(f"Output written to {output_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
#     "pydantic>=2.0",
# ]
# ///
"""Calculate reference SASA values using FreeSASA.

Converts mmCIF/PDB to a clean PDB (no hydrogens, no altlocs, no ligands),
then calculates SASA using FreeSASA's Shrake-Rupley algorithm.

Usage:
    ./scripts/data/calc_reference.py <input_file> [output_file]

Examples:
    ./scripts/data/calc_reference.py examples/1A0Q.cif.gz
    ./scripts/data/calc_reference.py 1abc.pdb reference.json
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import freesasa
import gemmi
from pydantic import BaseModel


class SASAResult(BaseModel):
    """Reference SASA values from FreeSASA."""

    total_area: float
    atom_areas: list[float]

    @property
    def n_atoms(self) -> int:
        return len(self.atom_areas)


def prepare_structure(path: str | Path) -> str:
    """Load structure with gemmi, clean it, and write to temp PDB.

    Returns path to temporary PDB file.
    """
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    # Write to temp PDB for FreeSASA
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    st.write_pdb(tmp.name)
    return tmp.name


def calc_sasa(pdb_path: str) -> SASAResult:
    """Calculate SASA using FreeSASA."""
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)

    atom_areas = [result.atomArea(i) for i in range(result.nAtoms())]

    return SASAResult(
        total_area=result.totalArea(),
        atom_areas=atom_areas,
    )


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 1

    input_path = Path(sys.argv[1])
    if not input_path.exists():
        print(f"Error: File not found: {input_path}", file=sys.stderr)
        return 1

    # Default output: same name with _sasa.json extension
    if len(sys.argv) >= 3:
        output_path = Path(sys.argv[2])
    else:
        stem = input_path.name
        for suffix in [".gz", ".cif", ".pdb", ".ent"]:
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
        output_path = input_path.parent / f"{stem}_sasa.json"

    # Prepare and calculate
    pdb_path = prepare_structure(input_path)
    result = calc_sasa(pdb_path)

    # Cleanup temp file
    Path(pdb_path).unlink()

    print(f"Calculated SASA for {result.n_atoms} atoms")
    print(f"Total area: {result.total_area:.2f} Å²")

    output_path.write_text(result.model_dump_json())
    print(f"Output written to {output_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())

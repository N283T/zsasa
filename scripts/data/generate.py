#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
#     "httpx>=0.27",
#     "pydantic>=2.0",
# ]
# ///
"""Generate benchmark dataset: download structures, create inputs, calculate references.

Downloads PDB structures from RCSB, converts to input JSON format,
and calculates FreeSASA reference SASA values.

Usage:
    ./scripts/data/generate.py [--structures-only] [--references-only]

Examples:
    ./scripts/data/generate.py                    # Full generation
    ./scripts/data/generate.py --structures-only  # Download structures only
    ./scripts/data/generate.py --references-only  # Generate references only
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import freesasa
import gemmi
import httpx

# Benchmark structures
STRUCTURES = [
    ("1crn", "tiny", "Crambin - small protein"),
    ("1ubq", "small", "Ubiquitin - small protein"),
    ("1a0q", "medium", "Lipid transfer protein"),
    ("3hhb", "medium", "Hemoglobin - multi-chain"),
    ("1aon", "large", "GroEL-GroES complex"),
    ("4v6x", "xlarge", "Ribosome"),
]

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"


@dataclass
class AtomData:
    """Atom data for SASA calculation."""

    x: list[float]
    y: list[float]
    z: list[float]
    r: list[float]
    residue: list[str]
    atom_name: list[str]
    element: list[int]

    def to_dict(self) -> dict:
        return {
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "r": self.r,
            "residue": self.residue,
            "atom_name": self.atom_name,
            "element": self.element,
        }


def download_structure(pdb_id: str, output_dir: Path) -> Path:
    """Download structure from RCSB PDB."""
    output_path = output_dir / f"{pdb_id}.cif.gz"

    if output_path.exists():
        print(f"  {pdb_id}: Already exists")
        return output_path

    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    print(f"  {pdb_id}: Downloading from RCSB...")

    with httpx.Client(timeout=60.0) as client:
        response = client.get(url)
        response.raise_for_status()

    output_path.write_bytes(response.content)
    print(f"  {pdb_id}: Downloaded ({len(response.content) / 1024:.1f} KB)")
    return output_path


def extract_atoms(structure_path: Path) -> AtomData:
    """Extract atom coordinates and radii from structure file."""
    st = gemmi.read_structure(str(structure_path))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    xs, ys, zs, rs = [], [], [], []
    residues, atom_names, elements = [], [], []

    model = st[0]
    for cra in model.all():
        xs.append(cra.atom.pos.x)
        ys.append(cra.atom.pos.y)
        zs.append(cra.atom.pos.z)
        rs.append(cra.atom.element.vdw_r)
        residues.append(cra.residue.name)
        atom_names.append(cra.atom.name)
        elements.append(cra.atom.element.atomic_number)

    return AtomData(
        x=xs, y=ys, z=zs, r=rs, residue=residues, atom_name=atom_names, element=elements
    )


def generate_input_json(structure_path: Path, output_path: Path) -> int:
    """Generate input JSON from structure file. Returns atom count."""
    atoms = extract_atoms(structure_path)
    with open(output_path, "w") as f:
        json.dump(atoms.to_dict(), f)
    return len(atoms.x)


def calculate_freesasa_reference(
    structure_path: Path,
    output_path: Path,
    n_points: int = 100,
    probe_radius: float = 1.4,
) -> dict:
    """Calculate FreeSASA reference SASA values."""
    # Load structure with gemmi, write clean PDB for FreeSASA
    st = gemmi.read_structure(str(structure_path))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    # Write temporary PDB
    tmp_pdb = output_path.parent / f"{output_path.stem}_tmp.pdb"
    st.write_pdb(str(tmp_pdb))

    try:
        # Calculate SASA with FreeSASA
        params = freesasa.Parameters()
        params.setAlgorithm(freesasa.ShrakeRupley)
        params.setNPoints(n_points)
        params.setProbeRadius(probe_radius)

        structure = freesasa.Structure(str(tmp_pdb))
        result = freesasa.calc(structure, params)

        # Extract per-atom areas
        atom_areas = []
        for i in range(structure.nAtoms()):
            atom_areas.append(result.atomArea(i))

        reference = {
            "structure": structure_path.stem.split(".")[0].upper(),
            "parameters": {
                "algorithm": "shrake-rupley",
                "n_points": n_points,
                "probe_radius": probe_radius,
            },
            "freesasa_version": "2.x",
            "n_atoms": len(atom_areas),
            "total_area": result.totalArea(),
            "atom_areas": atom_areas,
        }

        with open(output_path, "w") as f:
            json.dump(reference, f, indent=2)

        return reference

    finally:
        tmp_pdb.unlink(missing_ok=True)


def main() -> int:
    # Parse arguments
    structures_only = "--structures-only" in sys.argv
    references_only = "--references-only" in sys.argv

    # Setup directories
    base_dir = Path(__file__).parent.parent / "benchmarks"
    structures_dir = base_dir / "structures"
    inputs_dir = base_dir / "inputs"
    references_dir = base_dir / "references"

    for d in [structures_dir, inputs_dir, references_dir]:
        d.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Benchmark Data Generator")
    print("=" * 60)

    # Step 1: Download structures
    if not references_only:
        print("\n[1/3] Downloading structures...")
        for pdb_id, category, description in STRUCTURES:
            try:
                download_structure(pdb_id, structures_dir)
            except Exception as e:
                print(f"  {pdb_id}: ERROR - {e}")

    # Step 2: Generate input JSONs
    if not references_only:
        print("\n[2/3] Generating input JSONs...")
        for pdb_id, category, description in STRUCTURES:
            structure_path = structures_dir / f"{pdb_id}.cif.gz"
            input_path = inputs_dir / f"{pdb_id}.json"

            if not structure_path.exists():
                print(f"  {pdb_id}: Skipping (structure not found)")
                continue

            if input_path.exists() and not structures_only:
                with open(input_path) as f:
                    data = json.load(f)
                    n_atoms = len(data.get("x", []))
                print(f"  {pdb_id}: Already exists ({n_atoms} atoms)")
                continue

            try:
                n_atoms = generate_input_json(structure_path, input_path)
                print(f"  {pdb_id}: Generated ({n_atoms} atoms)")
            except Exception as e:
                print(f"  {pdb_id}: ERROR - {e}")

    if structures_only:
        print("\n[Done] Structures and inputs generated.")
        return 0

    # Step 3: Calculate FreeSASA references
    print("\n[3/3] Calculating FreeSASA references...")
    for pdb_id, category, description in STRUCTURES:
        input_path = inputs_dir / f"{pdb_id}.json"
        structure_path = structures_dir / f"{pdb_id}.cif.gz"
        reference_path = references_dir / f"{pdb_id}_n100_p1.4.json"

        if not structure_path.exists():
            print(f"  {pdb_id}: Skipping (structure not found)")
            continue

        if reference_path.exists():
            with open(reference_path) as f:
                ref = json.load(f)
            print(f"  {pdb_id}: Already exists (total: {ref['total_area']:.2f} Å²)")
            continue

        try:
            print(f"  {pdb_id}: Calculating...", end=" ", flush=True)
            ref = calculate_freesasa_reference(structure_path, reference_path)
            print(f"done (total: {ref['total_area']:.2f} Å²)")
        except Exception as e:
            print(f"ERROR - {e}")

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"{'PDB':<8} {'Category':<8} {'Atoms':>10} {'SASA (Å²)':>12}")
    print("-" * 40)

    for pdb_id, category, description in STRUCTURES:
        input_path = inputs_dir / f"{pdb_id}.json"
        reference_path = references_dir / f"{pdb_id}_n100_p1.4.json"

        n_atoms = "-"
        sasa = "-"

        if input_path.exists():
            with open(input_path) as f:
                data = json.load(f)
                n_atoms = len(data.get("x", []))

        if reference_path.exists():
            with open(reference_path) as f:
                ref = json.load(f)
                sasa = f"{ref['total_area']:.2f}"

        print(f"{pdb_id:<8} {category:<8} {n_atoms:>10} {sasa:>12}")

    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

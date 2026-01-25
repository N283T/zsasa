#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "freesasa>=2.0",
#     "gemmi>=0.7.4",
#     "rich>=13.0",
#     "typer>=0.9",
# ]
# ///
"""Generate benchmark JSON inputs from CIF files with ProtOr radii.

Creates minimal JSON files (x, y, z, r only) for fair benchmarking.
- Protein atoms only (HETATM excluded)
- Hydrogens excluded
- ProtOr radii pre-applied

Usage:
    ./scripts/data/generate_benchmark_json.py <cif_file> [output.json]
    ./scripts/data/generate_benchmark_json.py --all  # Generate all benchmark structures
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import freesasa
import gemmi
import typer
from rich.console import Console
from rich.table import Table

console = Console()

# ProtOr classifier (cached)
_protor_classifier: freesasa.Classifier | None = None


def get_protor_classifier() -> freesasa.Classifier:
    """Get cached ProtOr classifier."""
    global _protor_classifier
    if _protor_classifier is None:
        _protor_classifier = freesasa.Classifier.getStandardClassifier("protor")
    return _protor_classifier


def get_protor_radius(residue: str, atom_name: str) -> float:
    """Get ProtOr radius for an atom."""
    classifier = get_protor_classifier()
    radius = classifier.radius(residue, atom_name)

    # Fallback to element-based if not found
    if radius <= 0:
        element = atom_name.strip()[0] if atom_name.strip() else "C"
        fallback = {
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "S": 1.80,
            "P": 1.80,
        }
        radius = fallback.get(element, 1.70)

    return radius


def cif_to_benchmark_json(cif_path: Path, output_path: Path) -> int:
    """Convert CIF to minimal benchmark JSON with ProtOr radii.

    Returns number of atoms extracted.
    """
    # Read structure
    st = gemmi.read_structure(str(cif_path))

    # Clean up
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    rs: list[float] = []

    # Use first model only
    model = st[0]

    for cra in model.all():
        residue_name = cra.residue.name
        atom_name = cra.atom.name

        # Get ProtOr radius
        radius = get_protor_radius(residue_name, atom_name)

        xs.append(round(cra.atom.pos.x, 3))
        ys.append(round(cra.atom.pos.y, 3))
        zs.append(round(cra.atom.pos.z, 3))
        rs.append(round(radius, 2))

    # Minimal JSON (no residue/atom_name - radii are pre-computed)
    data = {"x": xs, "y": ys, "z": zs, "r": rs}

    with open(output_path, "w") as f:
        json.dump(data, f, separators=(",", ":"))  # Compact format

    return len(xs)


# Benchmark structures
BENCHMARK_STRUCTURES = [
    ("1crn", "Crambin"),
    ("1ubq", "Ubiquitin"),
    ("1a0q", "Lipid transfer protein"),
    ("3hhb", "Hemoglobin"),
    ("1aon", "GroEL-GroES"),
    ("4v6x", "Ribosome"),
]


def download_cif(pdb_id: str, output_dir: Path) -> Path:
    """Download CIF from RCSB if not present."""
    import urllib.request

    cif_path = output_dir / f"{pdb_id}.cif.gz"
    if cif_path.exists():
        return cif_path

    url = f"https://files.rcsb.org/download/{pdb_id}.cif.gz"
    console.print(f"  Downloading {pdb_id}...", style="dim")
    urllib.request.urlretrieve(url, cif_path)
    return cif_path


app = typer.Typer(help=__doc__)


@app.command()
def generate(
    cif_file: Annotated[
        Path | None,
        typer.Argument(help="Input CIF file (or use --all for benchmark set)"),
    ] = None,
    output: Annotated[
        Path | None,
        typer.Argument(help="Output JSON file"),
    ] = None,
    all_benchmarks: Annotated[
        bool,
        typer.Option("--all", "-a", help="Generate all benchmark structures"),
    ] = False,
) -> None:
    """Generate benchmark JSON from CIF file(s)."""

    if all_benchmarks:
        generate_all_benchmarks()
        return

    if cif_file is None:
        console.print("[red]Error: Provide CIF file or use --all[/red]")
        raise typer.Exit(1)

    if not cif_file.exists():
        console.print(f"[red]Error: File not found: {cif_file}[/red]")
        raise typer.Exit(1)

    # Default output path
    if output is None:
        stem = cif_file.name
        for suffix in [".gz", ".cif", ".pdb", ".ent"]:
            if stem.lower().endswith(suffix):
                stem = stem[: -len(suffix)]
        output = cif_file.parent / f"{stem}.json"

    n_atoms = cif_to_benchmark_json(cif_file, output)
    console.print(f"[green]Generated {output.name}[/green] ({n_atoms} atoms)")


def generate_all_benchmarks() -> None:
    """Generate JSON inputs for all benchmark structures."""
    base_dir = Path(__file__).parent.parent.parent / "benchmarks"
    structures_dir = base_dir / "structures"
    json_dir = base_dir / "inputs_json"

    structures_dir.mkdir(parents=True, exist_ok=True)
    json_dir.mkdir(parents=True, exist_ok=True)

    console.print("[bold]Generating benchmark JSON inputs...[/bold]\n")

    table = Table(title="Benchmark Inputs")
    table.add_column("PDB", style="cyan")
    table.add_column("Description")
    table.add_column("Atoms", justify="right")
    table.add_column("Size", justify="right")

    for pdb_id, description in BENCHMARK_STRUCTURES:
        # Download CIF if needed
        cif_path = download_cif(pdb_id, structures_dir)

        # Generate JSON
        json_path = json_dir / f"{pdb_id}.json"
        n_atoms = cif_to_benchmark_json(cif_path, json_path)

        # File size
        size_kb = json_path.stat().st_size / 1024

        table.add_row(
            pdb_id.upper(),
            description,
            f"{n_atoms:,}",
            f"{size_kb:.1f} KB",
        )

    console.print(table)
    console.print(f"\n[green]Done![/green] JSON inputs saved to: {json_dir}")


if __name__ == "__main__":
    app()

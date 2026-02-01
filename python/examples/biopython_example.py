#!/usr/bin/env python3
"""BioPython integration example.

This example demonstrates SASA calculation using BioPython for
structure file parsing.

BioPython is a comprehensive library for computational biology.
Use this integration if you're already using BioPython for:
- Structure manipulation
- Sequence analysis
- PDB/mmCIF file handling

Requirements:
    pip install zsasa[biopython]

Examples:
    # Run all examples
    python biopython_example.py

    # Use in your own code
    >>> from zsasa.integrations.biopython import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.pdb")
    >>> print(f"Total: {result.total_area:.1f} A^2")
"""

from pathlib import Path
from typing import Optional

from zsasa import aggregate_from_result
from zsasa.integrations.biopython import (
    calculate_sasa_from_model,
    calculate_sasa_from_structure,
    extract_atoms_from_model,
)

# Find example structure files
EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def get_structure_file() -> Optional[Path]:
    """Find an available example structure file."""
    pdb_file = EXAMPLES_DIR / "1ubq.pdb"
    cif_file = EXAMPLES_DIR / "1ubq.cif"

    if pdb_file.exists():
        return pdb_file
    elif cif_file.exists():
        return cif_file
    return None


def example_basic_sasa() -> None:
    """Calculate SASA from a structure file.

    The simplest usage: pass a file path to calculate_sasa_from_structure().
    Supports both PDB and mmCIF formats.
    """
    print("\n1. Basic SASA calculation from file")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    result = calculate_sasa_from_structure(structure_file)

    print(f"Total SASA:  {result.total_area:>8.1f} Å²")
    print(f"Polar SASA:  {result.polar_area:>8.1f} Å²")
    print(f"Apolar SASA: {result.apolar_area:>8.1f} Å²")
    print(f"Number of atoms: {len(result.atom_areas)}")


def example_biopython_structure() -> None:
    """Work with BioPython Structure objects.

    If you already have a Bio.PDB.Structure object from your own
    loading or manipulation, you can pass it directly.

    BioPython structure hierarchy:
    Structure -> Model -> Chain -> Residue -> Atom
    """
    print("\n2. Using BioPython Structure object")
    print("-" * 30)

    from Bio.PDB import PDBParser

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    # Load with BioPython's parser
    parser = PDBParser(QUIET=True)  # QUIET suppresses warnings
    structure = parser.get_structure("1ubq", str(structure_file))

    print(f"Structure ID: {structure.id}")
    print(f"Number of models: {len(list(structure.get_models()))}")

    # Calculate SASA for the first model
    model = structure[0]
    result = calculate_sasa_from_model(model)

    print(f"Model 0 SASA: {result.total_area:.1f} Å²")


def example_atom_inspection() -> None:
    """Extract and inspect atom data from a model.

    This shows what atoms are included in the SASA calculation
    and their properties.
    """
    print("\n3. Extract and inspect atom data")
    print("-" * 30)

    from Bio.PDB import PDBParser

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1ubq", str(structure_file))
    model = structure[0]

    # Extract protein atoms only (no HETATM)
    atom_data = extract_atoms_from_model(model, include_hetatm=False)

    print(f"Number of atoms (protein only): {len(atom_data)}")
    print()

    # Show first 5 atoms
    print("First 5 atoms:")
    print(f"{'Chain':<6} {'Res':<6} {'ResID':>6} {'Atom':<6}")
    print("-" * 28)

    for i in range(min(5, len(atom_data))):
        print(
            f"{atom_data.chain_ids[i]:<6} "
            f"{atom_data.residue_names[i]:<6} "
            f"{atom_data.residue_ids[i]:>6} "
            f"{atom_data.atom_names[i]:<6}"
        )


def example_per_residue_analysis() -> None:
    """Analyze per-residue SASA with RSA.

    Find exposed and buried residues based on their
    Relative Solvent Accessibility (RSA).
    """
    print("\n4. Per-residue SASA with RSA")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    result = calculate_sasa_from_structure(structure_file)
    residue_results = aggregate_from_result(result)

    # Categorize residues by exposure
    exposed = []  # RSA > 50%
    buried = []   # RSA < 10%

    for res in residue_results:
        if res.rsa is not None:
            if res.rsa > 0.5:
                exposed.append(res)
            elif res.rsa < 0.1:
                buried.append(res)

    print(f"Highly exposed residues (RSA > 50%): {len(exposed)}")
    for res in exposed[:5]:  # Show first 5
        print(f"  {res.chain_id}:{res.residue_name}{res.residue_id} RSA={res.rsa:.1%}")

    print()
    print(f"Buried residues (RSA < 10%): {len(buried)}")
    for res in buried[:5]:  # Show first 5
        print(f"  {res.chain_id}:{res.residue_name}{res.residue_id} RSA={res.rsa:.1%}")


def example_chain_analysis() -> None:
    """Analyze SASA by chain.

    Useful for multi-chain structures to see contribution
    of each chain to the total surface.
    """
    print("\n5. Chain-specific analysis")
    print("-" * 30)

    from Bio.PDB import PDBParser

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1ubq", str(structure_file))
    model = structure[0]

    # Iterate over chains
    print(f"{'Chain':<8} {'Residues':>10} {'Atoms':>8}")
    print("-" * 30)

    for chain in model:
        chain_id = chain.get_id()
        n_residues = len(list(chain.get_residues()))
        n_atoms = len(list(chain.get_atoms()))
        print(f"{chain_id:<8} {n_residues:>10} {n_atoms:>8}")

    print()
    print("Note: For multi-chain structures, you can calculate")
    print("SASA per chain by extracting atoms from specific chains")


def example_mmcif_support() -> None:
    """Load mmCIF format files.

    BioPython supports mmCIF through the MMCIFParser.
    mmCIF is the newer, more comprehensive format used by PDB.
    """
    print("\n6. mmCIF file support")
    print("-" * 30)

    pdb_file = EXAMPLES_DIR / "1ubq.pdb"
    cif_file = EXAMPLES_DIR / "1ubq.cif"

    if pdb_file.exists():
        result_pdb = calculate_sasa_from_structure(pdb_file)
        print(f"PDB format:   {result_pdb.total_area:.1f} Å²")

    if cif_file.exists():
        result_cif = calculate_sasa_from_structure(cif_file)
        print(f"mmCIF format: {result_cif.total_area:.1f} Å²")

        if pdb_file.exists():
            diff = abs(result_pdb.total_area - result_cif.total_area)
            print(f"Difference:   {diff:.1f} Å²")

    if not pdb_file.exists() and not cif_file.exists():
        print("Example files not found")


def example_filter_comparison() -> None:
    """Compare SASA with different atom filters.

    HETATM records include:
    - Water molecules (HOH)
    - Ligands
    - Metal ions
    - Other non-standard residues
    """
    print("\n7. HETATM comparison")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    result_with = calculate_sasa_from_structure(
        structure_file,
        include_hetatm=True,
    )

    result_without = calculate_sasa_from_structure(
        structure_file,
        include_hetatm=False,
    )

    print(f"With HETATM:    {len(result_with.atom_areas):4d} atoms, {result_with.total_area:.1f} Å²")
    print(f"Without HETATM: {len(result_without.atom_areas):4d} atoms, {result_without.total_area:.1f} Å²")
    print()
    print("Tip: Exclude HETATM for protein-only analysis")


def main() -> None:
    """Run all BioPython integration examples."""
    print("BioPython Integration Examples")
    print("=" * 50)

    structure_file = get_structure_file()
    if not structure_file:
        print(f"Example file not found: {EXAMPLES_DIR / '1ubq.pdb'}")
        return

    print(f"Using structure file: {structure_file.name}")

    # Run each example
    example_basic_sasa()
    example_biopython_structure()
    example_atom_inspection()
    example_per_residue_analysis()
    example_chain_analysis()
    example_mmcif_support()
    example_filter_comparison()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

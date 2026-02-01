#!/usr/bin/env python3
"""BioPython integration example.

This example demonstrates SASA calculation using BioPython for
structure file parsing.

Requirements:
    pip install zsasa[biopython]
"""

from pathlib import Path

from zsasa import aggregate_from_result
from zsasa.integrations.biopython import (
    calculate_sasa_from_model,
    calculate_sasa_from_structure,
    extract_atoms_from_model,
)


def main() -> None:
    print("BioPython Integration Examples")
    print("=" * 50)

    # Find example structure file
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    pdb_file = examples_dir / "1ubq.pdb"

    if not pdb_file.exists():
        print(f"Example file not found: {pdb_file}")
        return

    print(f"Using structure file: {pdb_file.name}")

    # Example 1: Simple SASA calculation from file
    print("\n1. Basic SASA calculation from file")
    print("-" * 30)

    result = calculate_sasa_from_structure(pdb_file)

    print(f"Total SASA:  {result.total_area:.1f} A^2")
    print(f"Polar SASA:  {result.polar_area:.1f} A^2")
    print(f"Apolar SASA: {result.apolar_area:.1f} A^2")
    print(f"Number of atoms: {len(result.atom_areas)}")

    # Example 2: Using BioPython Structure object
    print("\n2. Using BioPython Structure object")
    print("-" * 30)

    from Bio.PDB import PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("1ubq", str(pdb_file))

    print(f"Structure ID: {structure.id}")
    print(f"Number of models: {len(list(structure.get_models()))}")

    model = structure[0]
    result = calculate_sasa_from_model(model)

    print(f"Model 0 SASA: {result.total_area:.1f} A^2")

    # Example 3: Extract atoms and inspect
    print("\n3. Extract and inspect atom data")
    print("-" * 30)

    atom_data = extract_atoms_from_model(model, include_hetatm=False)

    print(f"Number of atoms (protein only): {len(atom_data)}")

    # Show first 5 atoms
    print("\nFirst 5 atoms:")
    print(f"{'Chain':<6} {'Res':<6} {'ResID':<6} {'Atom':<6}")
    print("-" * 28)
    for i in range(min(5, len(atom_data))):
        print(
            f"{atom_data.chain_ids[i]:<6} "
            f"{atom_data.residue_names[i]:<6} "
            f"{atom_data.residue_ids[i]:<6} "
            f"{atom_data.atom_names[i]:<6}"
        )

    # Example 4: Per-residue analysis
    print("\n4. Per-residue SASA with RSA")
    print("-" * 30)

    residue_results = aggregate_from_result(result)

    # Find most exposed and most buried residues
    exposed = []
    buried = []

    for res in residue_results:
        if res.rsa is not None:
            if res.rsa > 0.5:
                exposed.append(res)
            elif res.rsa < 0.1:
                buried.append(res)

    print(f"Highly exposed residues (RSA > 50%): {len(exposed)}")
    for res in exposed[:5]:
        print(f"  {res.chain_id}:{res.residue_name}{res.residue_id} RSA={res.rsa:.1%}")

    print(f"\nBuried residues (RSA < 10%): {len(buried)}")
    for res in buried[:5]:
        print(f"  {res.chain_id}:{res.residue_name}{res.residue_id} RSA={res.rsa:.1%}")

    # Example 5: Chain-specific analysis
    print("\n5. Chain-specific analysis")
    print("-" * 30)

    for chain in model:
        chain_id = chain.get_id()
        n_residues = len(list(chain.get_residues()))
        print(f"Chain {chain_id}: {n_residues} residues")

    # Example 6: mmCIF file support
    print("\n6. mmCIF file support")
    print("-" * 30)

    cif_file = examples_dir / "1ubq.cif"
    if cif_file.exists():
        result_cif = calculate_sasa_from_structure(cif_file)
        print(f"mmCIF SASA: {result_cif.total_area:.1f} A^2")
        print(f"PDB SASA:   {result.total_area:.1f} A^2")
        print(f"Difference: {abs(result_cif.total_area - result.total_area):.1f} A^2")
    else:
        print(f"mmCIF file not found: {cif_file}")


if __name__ == "__main__":
    main()

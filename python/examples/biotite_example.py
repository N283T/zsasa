#!/usr/bin/env python3
"""Biotite integration example.

This example demonstrates SASA calculation using Biotite for
structure file parsing. Also works with AtomWorks (built on Biotite).

Requirements:
    pip install zsasa[biotite]

For AtomWorks:
    pip install atomworks
"""

from pathlib import Path

from zsasa import aggregate_from_result
from zsasa.integrations.biotite import (
    calculate_sasa_from_atom_array,
    calculate_sasa_from_structure,
    extract_atoms_from_atom_array,
)


def main() -> None:
    print("Biotite Integration Examples")
    print("=" * 50)

    # Find example structure file
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    pdb_file = examples_dir / "1ubq.pdb"
    cif_file = examples_dir / "1ubq.cif"

    if not pdb_file.exists() and not cif_file.exists():
        print(f"Example files not found in {examples_dir}")
        return

    structure_file = cif_file if cif_file.exists() else pdb_file
    print(f"Using structure file: {structure_file.name}")

    # Example 1: Simple SASA calculation from file
    print("\n1. Basic SASA calculation from file")
    print("-" * 30)

    result = calculate_sasa_from_structure(structure_file)

    print(f"Total SASA:  {result.total_area:.1f} A^2")
    print(f"Polar SASA:  {result.polar_area:.1f} A^2")
    print(f"Apolar SASA: {result.apolar_area:.1f} A^2")
    print(f"Number of atoms: {len(result.atom_areas)}")

    # Example 2: Using Biotite AtomArray
    print("\n2. Using Biotite AtomArray")
    print("-" * 30)

    import biotite.structure.io as strucio

    atom_array = strucio.load_structure(str(structure_file))

    # Handle AtomArrayStack (multiple models)
    if hasattr(atom_array, "stack_depth"):
        print(f"Structure has {atom_array.stack_depth()} models")
        atom_array = atom_array[0]  # Use first model

    print(f"Number of atoms: {len(atom_array)}")
    print(f"Chains: {set(atom_array.chain_id)}")

    result = calculate_sasa_from_atom_array(atom_array)
    print(f"SASA: {result.total_area:.1f} A^2")

    # Example 3: Extract and inspect atom data
    print("\n3. Extract and inspect atom data")
    print("-" * 30)

    atom_data = extract_atoms_from_atom_array(atom_array, include_hetatm=False)

    print(f"Number of atoms (protein only): {len(atom_data)}")

    # Show first 5 atoms
    print("\nFirst 5 atoms:")
    print(f"{'Chain':<6} {'Res':<6} {'ResID':<6} {'Atom':<6} {'Element':<8}")
    print("-" * 36)
    for i in range(min(5, len(atom_data))):
        print(
            f"{atom_data.chain_ids[i]:<6} "
            f"{atom_data.residue_names[i]:<6} "
            f"{atom_data.residue_ids[i]:<6} "
            f"{atom_data.atom_names[i]:<6} "
            f"{atom_data.elements[i]:<8}"
        )

    # Example 4: Biotite filtering capabilities
    print("\n4. Biotite filtering capabilities")
    print("-" * 30)

    import biotite.structure as struc

    # Filter protein atoms only
    protein_mask = struc.filter_amino_acids(atom_array)
    protein_atoms = atom_array[protein_mask]
    print(f"Protein atoms: {len(protein_atoms)}")

    # Filter backbone atoms
    backbone_mask = struc.filter_backbone(atom_array)
    backbone_atoms = atom_array[backbone_mask]
    print(f"Backbone atoms: {len(backbone_atoms)}")

    # Calculate SASA for backbone only
    result_backbone = calculate_sasa_from_atom_array(backbone_atoms)
    print(f"Backbone SASA: {result_backbone.total_area:.1f} A^2")

    # Example 5: Per-residue analysis
    print("\n5. Per-residue SASA analysis")
    print("-" * 30)

    residue_results = aggregate_from_result(result)

    # Group by secondary structure tendency
    # (simplified: just show statistics)
    total_residues = len([r for r in residue_results if r.rsa is not None])
    exposed = len([r for r in residue_results if r.rsa is not None and r.rsa > 0.25])
    buried = len([r for r in residue_results if r.rsa is not None and r.rsa < 0.25])

    print(f"Total residues with RSA: {total_residues}")
    print(f"Exposed (RSA > 25%): {exposed} ({exposed/total_residues:.1%})")
    print(f"Buried (RSA < 25%):  {buried} ({buried/total_residues:.1%})")

    # Example 6: AtomWorks compatibility
    print("\n6. AtomWorks compatibility note")
    print("-" * 30)
    print("AtomWorks is built on Biotite, so you can use:")
    print()
    print("  from atomworks.io.utils.io_utils import load_any")
    print("  from zsasa.integrations.biotite import calculate_sasa_from_atom_array")
    print()
    print("  atom_array = load_any('protein.cif')")
    print("  result = calculate_sasa_from_atom_array(atom_array)")

    # Example 7: Compare with/without HETATM
    print("\n7. HETATM comparison")
    print("-" * 30)

    result_with = calculate_sasa_from_structure(structure_file, include_hetatm=True)
    result_without = calculate_sasa_from_structure(structure_file, include_hetatm=False)

    print(f"With HETATM:    {len(result_with.atom_areas):4d} atoms, {result_with.total_area:.1f} A^2")
    print(f"Without HETATM: {len(result_without.atom_areas):4d} atoms, {result_without.total_area:.1f} A^2")


if __name__ == "__main__":
    main()

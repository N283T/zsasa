#!/usr/bin/env python3
"""Gemmi integration example.

This example demonstrates SASA calculation using gemmi for
structure file parsing (mmCIF, PDB).

Requirements:
    pip install zsasa[gemmi]
"""

from pathlib import Path

from zsasa import aggregate_from_result
from zsasa.integrations.gemmi import (
    calculate_sasa_from_model,
    calculate_sasa_from_structure,
    extract_atoms_from_model,
)


def main() -> None:
    print("Gemmi Integration Examples")
    print("=" * 50)

    # Find example structure file
    examples_dir = Path(__file__).parent.parent.parent / "examples"
    pdb_file = examples_dir / "1ubq.pdb"
    cif_file = examples_dir / "1ubq.cif"

    if not pdb_file.exists() and not cif_file.exists():
        print(f"Example files not found in {examples_dir}")
        print("Please ensure 1ubq.pdb or 1ubq.cif exists.")
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

    # Example 2: Using gemmi Structure object
    print("\n2. Using gemmi Structure object")
    print("-" * 30)

    import gemmi

    structure = gemmi.read_structure(str(structure_file))
    print(f"Structure name: {structure.name}")
    print(f"Number of models: {len(structure)}")

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
    print(f"{'Chain':<6} {'Res':<6} {'Atom':<6} {'Element':<8}")
    print("-" * 28)
    for i in range(min(5, len(atom_data))):
        print(
            f"{atom_data.chain_ids[i]:<6} "
            f"{atom_data.residue_names[i]:<6} "
            f"{atom_data.atom_names[i]:<6} "
            f"{atom_data.elements[i]:<8}"
        )

    # Example 4: Per-residue SASA with RSA
    print("\n4. Per-residue SASA with RSA")
    print("-" * 30)

    residue_results = aggregate_from_result(result)

    # Show first 10 residues with RSA
    print(f"{'Chain':<6} {'Residue':<10} {'SASA':>8} {'RSA':>8}")
    print("-" * 34)

    count = 0
    for res in residue_results:
        if res.rsa is not None and count < 10:
            print(
                f"{res.chain_id:<6} "
                f"{res.residue_name}{res.residue_id:<6} "
                f"{res.total_area:>8.1f} "
                f"{res.rsa:>7.1%}"
            )
            count += 1

    # Example 5: Different algorithms
    print("\n5. Algorithm comparison")
    print("-" * 30)

    result_sr = calculate_sasa_from_structure(
        structure_file,
        algorithm="sr",
        n_points=100,
    )

    result_lr = calculate_sasa_from_structure(
        structure_file,
        algorithm="lr",
        n_slices=20,
    )

    print(f"Shrake-Rupley (100 points): {result_sr.total_area:.1f} A^2")
    print(f"Lee-Richards (20 slices):   {result_lr.total_area:.1f} A^2")
    print(f"Difference: {abs(result_sr.total_area - result_lr.total_area):.1f} A^2")

    # Example 6: Filter options
    print("\n6. Filter options (HETATM, hydrogens)")
    print("-" * 30)

    # With HETATM (default)
    result_with_het = calculate_sasa_from_structure(
        structure_file,
        include_hetatm=True,
    )

    # Without HETATM
    result_no_het = calculate_sasa_from_structure(
        structure_file,
        include_hetatm=False,
    )

    print(f"With HETATM:    {len(result_with_het.atom_areas)} atoms, {result_with_het.total_area:.1f} A^2")
    print(f"Without HETATM: {len(result_no_het.atom_areas)} atoms, {result_no_het.total_area:.1f} A^2")

    # Example 7: Different classifiers
    print("\n7. Classifier comparison")
    print("-" * 30)

    from zsasa import ClassifierType

    for classifier in [ClassifierType.NACCESS, ClassifierType.PROTOR, ClassifierType.OONS]:
        result = calculate_sasa_from_structure(
            structure_file,
            classifier=classifier,
        )
        print(f"{classifier.name:<10}: Total={result.total_area:.1f}, Polar={result.polar_area:.1f}, Apolar={result.apolar_area:.1f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Gemmi integration example.

This example demonstrates SASA calculation using gemmi for
structure file parsing (mmCIF, PDB).

Gemmi is a fast C++ library for macromolecular crystallography.
It's recommended for production use due to:
- Fast file parsing (especially for mmCIF)
- Support for compressed files (.gz)
- Comprehensive PDB/mmCIF support

Requirements:
    pip install zsasa[gemmi]

Examples:
    # Run all examples
    python gemmi_example.py

    # Use in your own code
    >>> from zsasa.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> print(f"Total: {result.total_area:.1f} A^2")
"""

from pathlib import Path
from typing import Optional

from zsasa import ClassifierType, aggregate_from_result
from zsasa.integrations.gemmi import (
    calculate_sasa_from_model,
    calculate_sasa_from_structure,
    extract_atoms_from_model,
)

# Find example structure files
EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def get_structure_file() -> Optional[Path]:
    """Find an available example structure file.

    Returns mmCIF file preferentially over PDB (more complete metadata).
    """
    cif_file = EXAMPLES_DIR / "1ubq.cif"
    pdb_file = EXAMPLES_DIR / "1ubq.pdb"

    if cif_file.exists():
        return cif_file
    elif pdb_file.exists():
        return pdb_file
    return None


def example_basic_sasa() -> None:
    """Calculate SASA from a structure file.

    The simplest way to calculate SASA: just pass a file path.
    Supports PDB, mmCIF, and compressed (.gz) formats.

    Returns a StructureResult with:
    - total_area: Total SASA in Å²
    - polar_area: Polar atom SASA
    - apolar_area: Apolar (hydrophobic) atom SASA
    - atom_areas: Per-atom SASA values
    """
    print("\n1. Basic SASA calculation from file")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    print(f"Using: {structure_file.name}")

    # Single function call for SASA calculation
    result = calculate_sasa_from_structure(structure_file)

    print()
    print(f"Total SASA:  {result.total_area:>8.1f} Å²")
    print(f"Polar SASA:  {result.polar_area:>8.1f} Å²")
    print(f"Apolar SASA: {result.apolar_area:>8.1f} Å²")
    print(f"Number of atoms: {len(result.atom_areas)}")


def example_gemmi_structure() -> None:
    """Work with gemmi Structure objects directly.

    If you already have a gemmi.Structure (e.g., from custom loading),
    you can calculate SASA from a specific model.

    This is useful for:
    - Multi-model structures (NMR ensembles)
    - Custom structure manipulation before SASA calculation
    """
    print("\n2. Using gemmi Structure object")
    print("-" * 30)

    import gemmi

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    # Load structure with gemmi
    structure = gemmi.read_structure(str(structure_file))

    print(f"Structure name: {structure.name}")
    print(f"Number of models: {len(structure)}")

    # For NMR structures with multiple models, you can iterate:
    # for model in structure:
    #     result = calculate_sasa_from_model(model)

    # Calculate SASA for first model
    model = structure[0]
    result = calculate_sasa_from_model(model)

    print(f"Model 0 SASA: {result.total_area:.1f} Å²")


def example_atom_inspection() -> None:
    """Extract and inspect atom data.

    The extract_atoms_from_model() function provides direct access
    to the atom data used for SASA calculation.

    This is useful for:
    - Debugging atom selection issues
    - Custom analysis workflows
    - Understanding what atoms are included
    """
    print("\n3. Extract and inspect atom data")
    print("-" * 30)

    import gemmi

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    structure = gemmi.read_structure(str(structure_file))
    model = structure[0]

    # Extract protein atoms only (no HETATM)
    atom_data = extract_atoms_from_model(model, include_hetatm=False)

    print(f"Number of atoms (protein only): {len(atom_data)}")
    print()

    # Atom data includes: chain_ids, residue_names, residue_ids,
    #                     atom_names, elements, coordinates
    print("First 5 atoms:")
    print(f"{'Chain':<6} {'Res':<6} {'ResID':>6} {'Atom':<6} {'Element':<8}")
    print("-" * 38)

    for i in range(min(5, len(atom_data))):
        print(
            f"{atom_data.chain_ids[i]:<6} "
            f"{atom_data.residue_names[i]:<6} "
            f"{atom_data.residue_ids[i]:>6} "
            f"{atom_data.atom_names[i]:<6} "
            f"{atom_data.elements[i]:<8}"
        )


def example_per_residue_analysis() -> None:
    """Aggregate SASA to per-residue values with RSA.

    aggregate_from_result() converts per-atom SASA to per-residue
    values and calculates RSA (Relative Solvent Accessibility).

    RSA interpretation:
    - > 25%: Surface residue (exposed)
    - < 25%: Core residue (buried)
    - ~5%: Typically in protein-protein interfaces
    """
    print("\n4. Per-residue SASA with RSA")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    result = calculate_sasa_from_structure(structure_file)

    # Aggregate atom SASA to residue SASA
    residue_results = aggregate_from_result(result)

    print(f"{'Chain':<6} {'Residue':<10} {'SASA':>8} {'RSA':>8} {'State':<8}")
    print("-" * 44)

    count = 0
    for res in residue_results:
        if res.rsa is not None and count < 10:
            state = "Exposed" if res.rsa > 0.25 else "Buried"
            print(
                f"{res.chain_id:<6} "
                f"{res.residue_name}{res.residue_id:<6} "
                f"{res.total_area:>8.1f} "
                f"{res.rsa:>7.1%} "
                f"{state:<8}"
            )
            count += 1

    # Summary
    total = len([r for r in residue_results if r.rsa is not None])
    exposed = len([r for r in residue_results if r.rsa is not None and r.rsa > 0.25])
    print()
    print(f"Total residues: {total}")
    print(f"Exposed (RSA > 25%): {exposed} ({exposed/total:.1%})")
    print(f"Buried (RSA < 25%): {total - exposed} ({(total-exposed)/total:.1%})")


def example_algorithm_comparison() -> None:
    """Compare Shrake-Rupley and Lee-Richards algorithms.

    Both algorithms should give similar results for proteins.
    Lee-Richards is generally faster for large structures.
    """
    print("\n5. Algorithm comparison")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    # Shrake-Rupley: point-based method
    # More points = more accurate but slower
    result_sr = calculate_sasa_from_structure(
        structure_file,
        algorithm="sr",   # Shrake-Rupley
        n_points=100,     # Test points per atom
    )

    # Lee-Richards: slice-based method
    # More slices = more accurate but slower
    result_lr = calculate_sasa_from_structure(
        structure_file,
        algorithm="lr",   # Lee-Richards
        n_slices=20,      # Slices per atom
    )

    print(f"Shrake-Rupley (100 points): {result_sr.total_area:.1f} Å²")
    print(f"Lee-Richards (20 slices):   {result_lr.total_area:.1f} Å²")
    print()
    diff = abs(result_sr.total_area - result_lr.total_area)
    print(f"Absolute difference: {diff:.1f} Å²")
    print(f"Relative difference: {diff / result_sr.total_area:.2%}")


def example_filter_options() -> None:
    """Control which atoms are included in SASA calculation.

    Key filter options:
    - include_hetatm: Include non-standard residues (ligands, water)
    - Skip hydrogens: Hydrogen atoms are always excluded

    For protein analysis, typically exclude HETATM to focus on
    the protein surface only.
    """
    print("\n6. Filter options (HETATM)")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    # Include HETATM (ligands, water, ions)
    result_with_het = calculate_sasa_from_structure(
        structure_file,
        include_hetatm=True,
    )

    # Exclude HETATM (protein only)
    result_no_het = calculate_sasa_from_structure(
        structure_file,
        include_hetatm=False,
    )

    print(f"With HETATM:    {len(result_with_het.atom_areas):4d} atoms, {result_with_het.total_area:.1f} Å²")
    print(f"Without HETATM: {len(result_no_het.atom_areas):4d} atoms, {result_no_het.total_area:.1f} Å²")
    print()
    print("Note: HETATM includes water, ligands, and ions")
    print("For most analyses, exclude HETATM to focus on protein surface")


def example_classifier_comparison() -> None:
    """Compare different atomic radii classifiers.

    zsasa includes three classifier databases:
    - NACCESS: Most commonly used, Hubbard & Thornton (1993)
    - PROTOR: Alternative parameterization
    - OONS: Ooi, Oobatake, Némethy & Scheraga (1987)

    Results will vary slightly due to different radii definitions.
    """
    print("\n7. Classifier comparison")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    print(f"{'Classifier':<10} {'Total':>10} {'Polar':>10} {'Apolar':>10}")
    print("-" * 44)

    for classifier in [ClassifierType.NACCESS, ClassifierType.PROTOR, ClassifierType.OONS]:
        result = calculate_sasa_from_structure(
            structure_file,
            classifier=classifier,
        )
        print(
            f"{classifier.name:<10} "
            f"{result.total_area:>10.1f} "
            f"{result.polar_area:>10.1f} "
            f"{result.apolar_area:>10.1f}"
        )

    print()
    print("Note: NACCESS is the most widely used in published literature")


def example_compressed_files() -> None:
    """Load compressed structure files.

    Gemmi natively supports gzip-compressed files.
    Simply use .gz extension (e.g., 3hhb.cif.gz).
    """
    print("\n8. Compressed file support")
    print("-" * 30)

    # Check for compressed example file
    compressed_file = EXAMPLES_DIR / "3hhb.cif.gz"

    if compressed_file.exists():
        result = calculate_sasa_from_structure(compressed_file)
        print(f"Loaded compressed file: {compressed_file.name}")
        print(f"Total SASA: {result.total_area:.1f} Å²")
        print(f"Number of atoms: {len(result.atom_areas)}")
    else:
        print(f"Compressed example not found: {compressed_file.name}")
        print()
        print("To use compressed files, simply pass the .gz path:")
        print("  result = calculate_sasa_from_structure('protein.cif.gz')")


def main() -> None:
    """Run all gemmi integration examples."""
    print("Gemmi Integration Examples")
    print("=" * 50)

    structure_file = get_structure_file()
    if not structure_file:
        print(f"Example files not found in {EXAMPLES_DIR}")
        print("Please ensure 1ubq.pdb or 1ubq.cif exists.")
        return

    print(f"Using structure file: {structure_file.name}")

    # Run each example
    example_basic_sasa()
    example_gemmi_structure()
    example_atom_inspection()
    example_per_residue_analysis()
    example_algorithm_comparison()
    example_filter_options()
    example_classifier_comparison()
    example_compressed_files()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

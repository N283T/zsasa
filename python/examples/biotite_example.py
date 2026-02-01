#!/usr/bin/env python3
"""Biotite integration example.

This example demonstrates SASA calculation using Biotite for
structure file parsing. Also works with AtomWorks (built on Biotite).

Biotite is a Python library for computational biology with:
- Fast structure file parsing (Cython-optimized)
- NumPy-based data structures
- Comprehensive filtering capabilities

Requirements:
    pip install zsasa[biotite]

For AtomWorks:
    pip install atomworks

Examples:
    # Run all examples
    python biotite_example.py

    # Use in your own code
    >>> from zsasa.integrations.biotite import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
    >>> print(f"Total: {result.total_area:.1f} A^2")
"""

from pathlib import Path
from typing import Optional

from zsasa import aggregate_from_result
from zsasa.integrations.biotite import (
    calculate_sasa_from_atom_array,
    calculate_sasa_from_structure,
    extract_atoms_from_atom_array,
)

# Find example structure files
EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def get_structure_file() -> Optional[Path]:
    """Find an available example structure file.

    Prefers mmCIF over PDB for better metadata handling.
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

    The simplest usage: pass a file path.
    Supports PDB, mmCIF, and other formats that Biotite supports.
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


def example_atom_array() -> None:
    """Work with Biotite AtomArray objects.

    Biotite stores structure data in AtomArray objects.
    This is useful when you've already loaded the structure
    or want to manipulate it before SASA calculation.
    """
    print("\n2. Using Biotite AtomArray")
    print("-" * 30)

    import biotite.structure.io as strucio

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    # Load structure with Biotite
    atom_array = strucio.load_structure(str(structure_file))

    # Handle AtomArrayStack (multiple models, e.g., NMR structures)
    if hasattr(atom_array, "stack_depth"):
        print(f"Structure has {atom_array.stack_depth()} models")
        atom_array = atom_array[0]  # Use first model

    print(f"Number of atoms: {len(atom_array)}")
    print(f"Chains: {sorted(set(atom_array.chain_id))}")

    # Calculate SASA
    result = calculate_sasa_from_atom_array(atom_array)
    print(f"SASA: {result.total_area:.1f} Å²")


def example_atom_inspection() -> None:
    """Extract and inspect atom data.

    See what atoms are being used for SASA calculation.
    """
    print("\n3. Extract and inspect atom data")
    print("-" * 30)

    import biotite.structure.io as strucio

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    atom_array = strucio.load_structure(str(structure_file))
    if hasattr(atom_array, "stack_depth"):
        atom_array = atom_array[0]

    # Extract protein atoms only (no HETATM)
    atom_data = extract_atoms_from_atom_array(atom_array, include_hetatm=False)

    print(f"Number of atoms (protein only): {len(atom_data)}")
    print()

    # Show first 5 atoms
    print("First 5 atoms:")
    print(f"{'Chain':<6} {'Res':<6} {'ResID':>6} {'Atom':<6} {'Element':<8}")
    print("-" * 40)

    for i in range(min(5, len(atom_data))):
        print(
            f"{atom_data.chain_ids[i]:<6} "
            f"{atom_data.residue_names[i]:<6} "
            f"{atom_data.residue_ids[i]:>6} "
            f"{atom_data.atom_names[i]:<6} "
            f"{atom_data.elements[i]:<8}"
        )


def example_biotite_filtering() -> None:
    """Use Biotite's powerful filtering capabilities.

    Biotite provides many filter functions for atom selection:
    - filter_amino_acids: Protein residues only
    - filter_backbone: Backbone atoms only
    - filter_nucleotides: DNA/RNA residues
    - Custom boolean masks
    """
    print("\n4. Biotite filtering capabilities")
    print("-" * 30)

    import biotite.structure as struc
    import biotite.structure.io as strucio

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    atom_array = strucio.load_structure(str(structure_file))
    if hasattr(atom_array, "stack_depth"):
        atom_array = atom_array[0]

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
    print(f"Backbone SASA: {result_backbone.total_area:.1f} Å²")

    # Compare to full protein
    result_protein = calculate_sasa_from_atom_array(protein_atoms)
    print(f"Full protein SASA: {result_protein.total_area:.1f} Å²")


def example_per_residue_analysis() -> None:
    """Analyze SASA per residue with RSA.

    Calculate burial statistics for the protein.
    """
    print("\n5. Per-residue SASA analysis")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    result = calculate_sasa_from_structure(structure_file)
    residue_results = aggregate_from_result(result)

    # Calculate statistics
    residues_with_rsa = [r for r in residue_results if r.rsa is not None]
    total_residues = len(residues_with_rsa)
    exposed = len([r for r in residues_with_rsa if r.rsa > 0.25])
    buried = len([r for r in residues_with_rsa if r.rsa < 0.25])

    print(f"Total residues with RSA: {total_residues}")
    print(f"Exposed (RSA > 25%): {exposed} ({exposed/total_residues:.1%})")
    print(f"Buried (RSA < 25%):  {buried} ({buried/total_residues:.1%})")


def example_atomworks_compatibility() -> None:
    """AtomWorks compatibility.

    AtomWorks is built on Biotite, so you can use zsasa's
    Biotite integration directly with AtomWorks structures.
    """
    print("\n6. AtomWorks compatibility note")
    print("-" * 30)

    print("AtomWorks is built on Biotite, so you can use:")
    print()
    print("  from atomworks.io.utils.io_utils import load_any")
    print("  from zsasa.integrations.biotite import calculate_sasa_from_atom_array")
    print()
    print("  atom_array = load_any('protein.cif')")
    print("  result = calculate_sasa_from_atom_array(atom_array)")
    print()
    print("The same integration works for both libraries.")


def example_hetatm_comparison() -> None:
    """Compare SASA with and without HETATM.

    HETATM includes non-protein atoms:
    - Water molecules
    - Ligands
    - Metal ions
    """
    print("\n7. HETATM comparison")
    print("-" * 30)

    structure_file = get_structure_file()
    if not structure_file:
        print("Example files not found")
        return

    result_with = calculate_sasa_from_structure(structure_file, include_hetatm=True)
    result_without = calculate_sasa_from_structure(structure_file, include_hetatm=False)

    print(f"With HETATM:    {len(result_with.atom_areas):4d} atoms, {result_with.total_area:.1f} Å²")
    print(f"Without HETATM: {len(result_without.atom_areas):4d} atoms, {result_without.total_area:.1f} Å²")


def main() -> None:
    """Run all Biotite integration examples."""
    print("Biotite Integration Examples")
    print("=" * 50)

    structure_file = get_structure_file()
    if not structure_file:
        print(f"Example files not found in {EXAMPLES_DIR}")
        return

    print(f"Using structure file: {structure_file.name}")

    # Run each example
    example_basic_sasa()
    example_atom_array()
    example_atom_inspection()
    example_biotite_filtering()
    example_per_residue_analysis()
    example_atomworks_compatibility()
    example_hetatm_comparison()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

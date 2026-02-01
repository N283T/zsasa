#!/usr/bin/env python3
"""Atom classification and radius lookup example.

This example demonstrates how to use built-in classifiers
to get atomic radii and polar/apolar classifications.

zsasa includes three classifier databases:
- NACCESS: Most commonly used, good for proteins
- PROTOR: Alternative with slightly different values
- OONS: Ooi et al. parameterization

Requirements:
    pip install zsasa

Examples:
    # Run all examples
    python classifier.py

    # Use in your own code
    >>> from zsasa import classify_atoms, ClassifierType
    >>> residues = ["ALA", "ALA"]
    >>> atoms = ["CA", "O"]
    >>> result = classify_atoms(residues, atoms, ClassifierType.NACCESS)
    >>> print(result.radii)  # [1.87, 1.40]
"""

import numpy as np

from zsasa import (
    AtomClass,
    ClassifierType,
    calculate_rsa,
    calculate_sasa,
    classify_atoms,
    get_atom_class,
    get_max_sasa,
    get_radius,
    guess_radius,
    guess_radius_from_atom_name,
)


def single_atom_radius() -> None:
    """Look up radius for a specific atom type.

    The get_radius() function returns the van der Waals radius
    for a given residue/atom combination from the specified classifier.

    Different classifiers have slightly different parameterizations
    based on different training data and methodologies.
    """
    print("\n1. Single atom radius lookup")
    print("-" * 30)

    # Compare radii across classifiers for the same atom
    print("Alanine CA (alpha carbon) radius:")
    print()
    for classifier in [ClassifierType.NACCESS, ClassifierType.PROTOR, ClassifierType.OONS]:
        radius = get_radius("ALA", "CA", classifier)
        print(f"  {classifier.name:8s}: {radius:.2f} Å")

    print()
    print("Note: Small differences reflect different parameterization methods")


def polarity_classification() -> None:
    """Classify atoms as polar or apolar.

    Polar atoms (N, O, etc.) can form hydrogen bonds with water.
    Apolar atoms (C, S) are hydrophobic.

    This classification is useful for:
    - Calculating buried hydrophobic surface
    - Predicting protein-protein interfaces
    - Analyzing solvation effects
    """
    print("\n2. Atom polarity classification")
    print("-" * 30)

    # Test various atom types
    test_atoms = [
        ("ALA", "CA", "Alpha carbon - backbone"),
        ("ALA", "N", "Amide nitrogen - backbone"),
        ("ALA", "O", "Carbonyl oxygen - backbone"),
        ("ALA", "CB", "Beta carbon - sidechain"),
        ("SER", "OG", "Hydroxyl oxygen - polar sidechain"),
        ("PHE", "CZ", "Aromatic carbon - apolar sidechain"),
    ]

    print(f"{'Residue':<8} {'Atom':<6} {'Class':<8} {'Description'}")
    print("-" * 55)

    for res, atom, desc in test_atoms:
        atom_class = get_atom_class(res, atom, ClassifierType.NACCESS)
        class_name = {
            AtomClass.POLAR: "POLAR",
            AtomClass.APOLAR: "APOLAR",
            AtomClass.UNKNOWN: "UNKNOWN",
        }[atom_class]
        print(f"{res:<8} {atom:<6} {class_name:<8} {desc}")


def guess_from_element() -> None:
    """Guess atomic radius from element symbol.

    When residue/atom names are unavailable (e.g., for ligands),
    you can guess the radius from the element symbol alone.

    This uses standard van der Waals radii from:
    - Bondi (1964) for common elements
    - Cambridge Structural Database for metals
    """
    print("\n3. Guess radius from element symbol")
    print("-" * 30)

    elements = [
        ("C", "Carbon"),
        ("N", "Nitrogen"),
        ("O", "Oxygen"),
        ("S", "Sulfur"),
        ("P", "Phosphorus"),
        ("H", "Hydrogen"),
        ("FE", "Iron"),
        ("ZN", "Zinc"),
        ("MG", "Magnesium"),
        ("CA", "Calcium"),
    ]

    print(f"{'Element':<10} {'Radius':>8}  {'Name'}")
    print("-" * 35)

    for elem, name in elements:
        radius = guess_radius(elem)
        if not np.isnan(radius):
            print(f"{elem:<10} {radius:>8.2f} Å  {name}")
        else:
            print(f"{elem:<10} {'Unknown':>8}  {name}")


def guess_from_atom_name() -> None:
    """Guess radius from PDB atom name.

    PDB atom names encode element information:
    - First 1-2 characters often indicate element
    - "CA" = alpha carbon, "FE" = iron, etc.

    This is useful when you have PDB atom names but
    not full residue context.
    """
    print("\n4. Guess radius from PDB atom name")
    print("-" * 30)

    atom_names = [
        ("CA", "Alpha carbon"),
        ("CB", "Beta carbon"),
        ("N", "Nitrogen"),
        ("O", "Oxygen"),
        ("OG", "Gamma oxygen (Ser)"),
        ("FE", "Iron"),
        ("ZN", "Zinc"),
        ("MG1", "Magnesium (alternate)"),
    ]

    print(f"{'Atom Name':<12} {'Radius':>8}  {'Description'}")
    print("-" * 40)

    for name, desc in atom_names:
        radius = guess_radius_from_atom_name(name)
        if not np.isnan(radius):
            print(f"{name:<12} {radius:>8.2f} Å  {desc}")
        else:
            print(f"{name:<12} {'Unknown':>8}  {desc}")


def batch_classification() -> None:
    """Classify multiple atoms at once.

    For efficiency, use classify_atoms() to process many atoms
    in a single call. This returns:
    - radii: numpy array of van der Waals radii
    - classes: numpy array of polarity classes (0=polar, 1=apolar, 2=unknown)
    """
    print("\n5. Batch atom classification")
    print("-" * 30)

    # Example peptide backbone atoms
    residues = ["ALA", "ALA", "ALA", "ALA", "GLY", "GLY"]
    atoms = ["N", "CA", "C", "O", "N", "CA"]

    # Classify all at once
    result = classify_atoms(residues, atoms, ClassifierType.NACCESS)

    print(f"{'Residue':<10} {'Atom':<8} {'Radius':>10} {'Class':<10}")
    print("-" * 40)

    class_names = {0: "POLAR", 1: "APOLAR", 2: "UNKNOWN"}

    for i, (res, atom) in enumerate(zip(residues, atoms)):
        radius = result.radii[i]
        cls = result.classes[i]
        print(f"{res:<10} {atom:<8} {radius:>10.2f} {class_names[cls]:<10}")

    # Summary statistics
    n_polar = np.sum(result.classes == AtomClass.POLAR)
    n_apolar = np.sum(result.classes == AtomClass.APOLAR)
    print()
    print(f"Total: {len(atoms)} atoms")
    print(f"Polar: {n_polar}, Apolar: {n_apolar}")


def max_sasa_lookup() -> None:
    """Look up maximum SASA for residues.

    Max SASA values represent the theoretical maximum solvent-accessible
    surface area for each residue type in a Gly-X-Gly tripeptide context.

    These values are used to calculate RSA (Relative Solvent Accessibility).
    """
    print("\n6. Max SASA lookup for residues")
    print("-" * 30)

    # Standard amino acids by size
    residues = [
        ("GLY", "Glycine (smallest)"),
        ("ALA", "Alanine"),
        ("SER", "Serine (polar)"),
        ("PHE", "Phenylalanine (aromatic)"),
        ("TRP", "Tryptophan (largest)"),
        ("GLU", "Glutamate (charged)"),
    ]

    print(f"{'Residue':<6} {'Max SASA':>12}  {'Description'}")
    print("-" * 45)

    for res, desc in residues:
        max_sasa = get_max_sasa(res)
        if not np.isnan(max_sasa):
            print(f"{res:<6} {max_sasa:>12.1f} Å²  {desc}")
        else:
            print(f"{res:<6} {'Unknown':>12}  {desc}")


def rsa_calculation() -> None:
    """Calculate Relative Solvent Accessibility (RSA).

    RSA = observed_SASA / max_SASA × 100%

    RSA values indicate burial state:
    - RSA > 25%: Exposed (surface residue)
    - RSA < 25%: Buried (core residue)

    Typical cutoffs vary by application:
    - Structural biology: 25%
    - Interface detection: 5-10%
    """
    print("\n7. RSA (Relative Solvent Accessibility) calculation")
    print("-" * 30)

    # Example: Alanine residue with different burial states
    test_cases = [
        (107.95, "Fully exposed (max)"),
        (80.0, "Partially exposed"),
        (50.0, "Partially buried"),
        (20.0, "Mostly buried"),
        (5.0, "Nearly completely buried"),
    ]

    max_sasa = get_max_sasa("ALA")
    print(f"Alanine max SASA: {max_sasa:.1f} Å²")
    print()
    print(f"{'Observed SASA':>15} {'RSA':>10} {'State'}")
    print("-" * 50)

    for observed, description in test_cases:
        rsa = calculate_rsa(observed, "ALA")
        state = "Exposed" if rsa > 0.25 else "Buried"
        print(f"{observed:>15.1f} {rsa:>10.1%} {description}")


def complete_workflow() -> None:
    """Complete workflow: classify atoms, calculate SASA, analyze polarity.

    This demonstrates a realistic workflow:
    1. Define atom types (residue + atom names)
    2. Get radii and classifications
    3. Calculate SASA
    4. Analyze polar vs apolar surface
    """
    print("\n8. Complete workflow: SASA with polar/apolar breakdown")
    print("-" * 30)

    # Small peptide: Ala-Gly-Ser backbone + sidechains
    residues = [
        "ALA", "ALA", "ALA", "ALA",  # Ala backbone + CB
        "GLY", "GLY", "GLY",          # Gly backbone (no sidechain)
        "SER", "SER", "SER", "SER", "SER",  # Ser backbone + OG
    ]
    atoms = [
        "N", "CA", "C", "O",     # Ala
        "N", "CA", "C",          # Gly
        "N", "CA", "C", "O", "OG",  # Ser
    ]

    # Step 1: Classify atoms to get radii
    print("Step 1: Classify atoms")
    classification = classify_atoms(residues, atoms, ClassifierType.NACCESS)
    print(f"  Classified {len(atoms)} atoms")

    # Step 2: Create coordinates (linear chain for demo)
    print("Step 2: Generate coordinates")
    np.random.seed(42)
    n_atoms = len(residues)
    coords = np.zeros((n_atoms, 3))
    for i in range(n_atoms):
        coords[i] = [i * 1.5, np.random.uniform(-0.5, 0.5), np.random.uniform(-0.5, 0.5)]
    print(f"  Generated {n_atoms} atom positions")

    # Step 3: Calculate SASA
    print("Step 3: Calculate SASA")
    sasa_result = calculate_sasa(coords, classification.radii)
    print(f"  Total SASA: {sasa_result.total_area:.1f} Å²")

    # Step 4: Aggregate by polarity
    print("Step 4: Analyze by polarity")
    polar_mask = classification.classes == AtomClass.POLAR
    apolar_mask = classification.classes == AtomClass.APOLAR

    polar_sasa = np.sum(sasa_result.atom_areas[polar_mask])
    apolar_sasa = np.sum(sasa_result.atom_areas[apolar_mask])
    total = sasa_result.total_area

    print()
    print("Results:")
    print(f"  Total SASA:  {total:>8.1f} Å² (100%)")
    print(f"  Polar SASA:  {polar_sasa:>8.1f} Å² ({polar_sasa/total:>5.1%})")
    print(f"  Apolar SASA: {apolar_sasa:>8.1f} Å² ({apolar_sasa/total:>5.1%})")


def main() -> None:
    """Run all classification examples."""
    print("Atom Classification Examples")
    print("=" * 50)

    # Run each example
    single_atom_radius()
    polarity_classification()
    guess_from_element()
    guess_from_atom_name()
    batch_classification()
    max_sasa_lookup()
    rsa_calculation()
    complete_workflow()

    print()
    print("=" * 50)
    print("All examples completed!")


if __name__ == "__main__":
    main()

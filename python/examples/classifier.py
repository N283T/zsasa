#!/usr/bin/env python3
"""Atom classification and radius lookup example.

This example demonstrates how to use built-in classifiers
to get atomic radii and polar/apolar classifications.

Requirements:
    pip install zsasa
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


def main() -> None:
    print("Atom Classification Examples")
    print("=" * 50)

    # Example 1: Get radius for specific atom
    print("\n1. Single atom radius lookup")
    print("-" * 30)

    for classifier in [ClassifierType.NACCESS, ClassifierType.PROTOR, ClassifierType.OONS]:
        radius = get_radius("ALA", "CA", classifier)
        print(f"{classifier.name}: ALA CA = {radius:.2f} A")

    # Example 2: Atom polarity classification
    print("\n2. Atom polarity classification")
    print("-" * 30)

    test_atoms = [
        ("ALA", "CA"),   # Carbon - apolar
        ("ALA", "N"),    # Nitrogen - polar
        ("ALA", "O"),    # Oxygen - polar
        ("ALA", "CB"),   # Carbon - apolar
        ("SER", "OG"),   # Hydroxyl oxygen - polar
        ("PHE", "CZ"),   # Aromatic carbon - apolar
    ]

    for res, atom in test_atoms:
        atom_class = get_atom_class(res, atom, ClassifierType.NACCESS)
        class_name = {
            AtomClass.POLAR: "POLAR",
            AtomClass.APOLAR: "APOLAR",
            AtomClass.UNKNOWN: "UNKNOWN",
        }[atom_class]
        print(f"{res:3s} {atom:4s} -> {class_name}")

    # Example 3: Guess radius from element
    print("\n3. Guess radius from element symbol")
    print("-" * 30)

    elements = ["C", "N", "O", "S", "P", "H", "FE", "ZN", "MG", "CA"]
    for elem in elements:
        radius = guess_radius(elem)
        if not np.isnan(radius):
            print(f"{elem:2s} -> {radius:.2f} A")
        else:
            print(f"{elem:2s} -> Unknown")

    # Example 4: Guess radius from PDB atom name
    print("\n4. Guess radius from PDB atom name")
    print("-" * 30)

    atom_names = ["CA", "CB", "N", "O", "OG", "FE", "ZN", "MG1"]
    for name in atom_names:
        radius = guess_radius_from_atom_name(name)
        if not np.isnan(radius):
            print(f"{name:4s} -> {radius:.2f} A")
        else:
            print(f"{name:4s} -> Unknown")

    # Example 5: Batch classification
    print("\n5. Batch atom classification")
    print("-" * 30)

    residues = ["ALA", "ALA", "ALA", "ALA", "GLY", "GLY"]
    atoms = ["N", "CA", "C", "O", "N", "CA"]

    result = classify_atoms(residues, atoms, ClassifierType.NACCESS)

    print(f"{'Residue':<8} {'Atom':<6} {'Radius':>8} {'Class':<8}")
    print("-" * 32)
    for i, (res, atom) in enumerate(zip(residues, atoms)):
        radius = result.radii[i]
        cls = result.classes[i]
        class_name = {0: "POLAR", 1: "APOLAR", 2: "UNKNOWN"}[cls]
        print(f"{res:<8} {atom:<6} {radius:>8.2f} {class_name:<8}")

    # Example 6: Max SASA and RSA calculation
    print("\n6. Max SASA lookup and RSA calculation")
    print("-" * 30)

    residues_3letter = ["ALA", "GLY", "SER", "TRP", "GLU"]
    for res in residues_3letter:
        max_sasa = get_max_sasa(res)
        if not np.isnan(max_sasa):
            print(f"{res}: Max SASA = {max_sasa:.1f} A^2")
        else:
            print(f"{res}: Max SASA = Unknown")

    # Calculate RSA for a single residue
    print("\n7. RSA (Relative Solvent Accessibility) calculation")
    print("-" * 30)

    # Example: Alanine residue with some burial
    observed_sasa = 50.0  # A^2
    rsa = calculate_rsa(observed_sasa, "ALA")
    print(f"Observed SASA: {observed_sasa:.1f} A^2")
    print(f"Max SASA (ALA): {get_max_sasa('ALA'):.1f} A^2")
    print(f"RSA: {rsa:.1%}")

    # Example 8: Complete workflow - classify, calculate SASA, get polar/apolar
    print("\n8. Complete workflow: SASA with polar/apolar breakdown")
    print("-" * 30)

    # Small peptide: Ala-Gly-Ser backbone
    residues = ["ALA", "ALA", "ALA", "ALA", "GLY", "GLY", "GLY", "SER", "SER", "SER", "SER", "SER"]
    atoms = ["N", "CA", "C", "O", "N", "CA", "C", "N", "CA", "C", "O", "OG"]

    # Get radii and classes
    classification = classify_atoms(residues, atoms, ClassifierType.NACCESS)

    # Create some example coordinates (spread out for visibility)
    np.random.seed(42)
    n_atoms = len(residues)
    coords = np.zeros((n_atoms, 3))
    for i in range(n_atoms):
        coords[i] = [i * 1.5, np.random.uniform(-0.5, 0.5), np.random.uniform(-0.5, 0.5)]

    # Calculate SASA
    sasa_result = calculate_sasa(coords, classification.radii)

    # Aggregate by polarity
    polar_mask = classification.classes == AtomClass.POLAR
    apolar_mask = classification.classes == AtomClass.APOLAR

    polar_sasa = np.sum(sasa_result.atom_areas[polar_mask])
    apolar_sasa = np.sum(sasa_result.atom_areas[apolar_mask])

    print(f"Total SASA:  {sasa_result.total_area:.1f} A^2")
    print(f"Polar SASA:  {polar_sasa:.1f} A^2 ({polar_sasa/sasa_result.total_area:.1%})")
    print(f"Apolar SASA: {apolar_sasa:.1f} A^2 ({apolar_sasa/sasa_result.total_area:.1%})")


if __name__ == "__main__":
    main()

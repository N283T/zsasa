"""Tests for shared integration radius assignment policy."""

from __future__ import annotations

import numpy as np
import pytest

from zsasa.classifier import ClassifierType
from zsasa.integrations._types import AtomData, classify_atom_data


def _atom_data(*, element: str) -> AtomData:
    return AtomData(
        coords=np.array([[0.0, 0.0, 0.0]], dtype=np.float64),
        residue_names=["UNK"],
        atom_names=["ZZ1"],
        chain_ids=["A"],
        residue_ids=[42],
        elements=[element],
    )


def test_unknown_classifier_radius_falls_back_to_element_radius() -> None:
    """Unknown residue/atom pairs should use element-derived radii when possible."""
    classification = classify_atom_data(_atom_data(element="C"), ClassifierType.CCD)

    assert classification.radii[0] == pytest.approx(1.70)


@pytest.mark.parametrize("classifier", [ClassifierType.NACCESS, ClassifierType.OONS])
def test_element_decides_for_atoms_outside_naccess_and_oons_tables(
    classifier: ClassifierType,
) -> None:
    """The element, not the atom name, sets the radius of hydrogens and ligand atoms."""
    # (residue, atom, element, expected radius)
    atoms = [
        ("SER", "HG", "H", 1.10),  # not mercury
        ("VAL", "HG21", "H", 1.10),
        ("HEM", "NA", "N", 1.55),  # not sodium
        ("ATP", "PB", "P", 1.80),  # not lead
        ("PCA", "CD", "C", 1.70),  # not cadmium
        ("LIG", "CL1", "CL", 1.75),  # not carbon
        ("CUA", "CU1", "Cu", 1.40),  # not carbon
        ("HG", "HG", "HG", 1.55),  # mercury ion
        ("CD", "CD", "CD", 1.58),  # cadmium ion
    ]
    atom_data = AtomData(
        coords=np.zeros((len(atoms), 3), dtype=np.float64),
        residue_names=[a[0] for a in atoms],
        atom_names=[a[1] for a in atoms],
        chain_ids=["A"] * len(atoms),
        residue_ids=list(range(1, len(atoms) + 1)),
        elements=[a[2] for a in atoms],
    )

    classification = classify_atom_data(atom_data, classifier)

    np.testing.assert_allclose(classification.radii, [a[3] for a in atoms])


@pytest.mark.parametrize(
    ("classifier", "expected"),
    [(ClassifierType.NACCESS, [1.87, 1.76]), (ClassifierType.OONS, [2.00, 1.75])],
)
def test_element_does_not_replace_table_radii(
    classifier: ClassifierType, expected: list[float]
) -> None:
    """Atoms listed in the classifier tables keep their radii."""
    atom_data = AtomData(
        coords=np.zeros((2, 3), dtype=np.float64),
        residue_names=["ARG", "PHE"],
        atom_names=["CD", "CD1"],
        chain_ids=["A", "A"],
        residue_ids=[1, 2],
        elements=["C", "C"],
    )

    classification = classify_atom_data(atom_data, classifier)

    np.testing.assert_allclose(classification.radii, expected)


def test_name_guess_is_kept_when_element_is_missing() -> None:
    """Without an element, the NACCESS guess from the names is the last resort."""
    atom_data = AtomData(
        coords=np.zeros((2, 3), dtype=np.float64),
        residue_names=["SER", "ZN"],
        atom_names=["HG", "ZN"],
        chain_ids=["A", "A"],
        residue_ids=[1, 2],
        elements=["", ""],
    )

    classification = classify_atom_data(atom_data, ClassifierType.NACCESS)

    np.testing.assert_allclose(classification.radii, [1.10, 1.39])


def test_unknown_radius_without_element_fails_with_atom_identifier() -> None:
    """Unknown radii should identify the atom clearly instead of leaking NaN to CFFI."""
    with pytest.raises(ValueError, match=r"chain A residue 42 UNK atom ZZ1"):
        classify_atom_data(_atom_data(element=""), ClassifierType.CCD)

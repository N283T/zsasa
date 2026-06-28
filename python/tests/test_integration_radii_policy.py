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


def test_unknown_radius_without_element_fails_with_atom_identifier() -> None:
    """Unknown radii should identify the atom clearly instead of leaking NaN to CFFI."""
    with pytest.raises(ValueError, match=r"chain A residue 42 UNK atom ZZ1"):
        classify_atom_data(_atom_data(element=""), ClassifierType.CCD)

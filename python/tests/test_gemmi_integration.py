"""Tests for gemmi integration.

These tests require gemmi to be installed.
Run with: pip install freesasa-zig[gemmi] pytest
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

# Skip all tests if gemmi is not installed
gemmi = pytest.importorskip("gemmi")

from freesasa_zig import AtomClass, ClassifierType  # noqa: E402
from freesasa_zig.integrations.gemmi import (  # noqa: E402
    AtomData,
    SasaResultWithAtoms,
    calculate_sasa_from_model,
    calculate_sasa_from_structure,
    extract_atoms_from_model,
)

# =============================================================================
# Test Fixtures
# =============================================================================


@pytest.fixture
def simple_structure() -> gemmi.Structure:
    """Create a simple structure with a few atoms for testing."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")

    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.seqid = gemmi.SeqId("1")

    # Add backbone atoms
    for name, element, pos in [
        ("N", "N", (0.0, 0.0, 0.0)),
        ("CA", "C", (1.5, 0.0, 0.0)),
        ("C", "C", (2.5, 1.2, 0.0)),
        ("O", "O", (2.3, 2.4, 0.0)),
        ("CB", "C", (1.8, -1.5, 0.0)),
    ]:
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element(element)
        atom.pos = gemmi.Position(*pos)
        residue.add_atom(atom)

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)

    return structure


@pytest.fixture
def structure_with_hetatm() -> gemmi.Structure:
    """Create a structure with HETATM records (water)."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")

    # Protein chain
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.seqid = gemmi.SeqId("1")
    residue.het_flag = "A"  # ATOM record

    atom = gemmi.Atom()
    atom.name = "CA"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(0.0, 0.0, 0.0)
    residue.add_atom(atom)
    chain.add_residue(residue)

    # Water
    water = gemmi.Residue()
    water.name = "HOH"
    water.seqid = gemmi.SeqId("100")
    water.het_flag = "H"  # HETATM record

    water_atom = gemmi.Atom()
    water_atom.name = "O"
    water_atom.element = gemmi.Element("O")
    water_atom.pos = gemmi.Position(5.0, 0.0, 0.0)
    water.add_atom(water_atom)
    chain.add_residue(water)

    model.add_chain(chain)
    structure.add_model(model)

    return structure


@pytest.fixture
def structure_with_hydrogens() -> gemmi.Structure:
    """Create a structure with hydrogen atoms."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")

    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "ALA"
    residue.seqid = gemmi.SeqId("1")

    # Heavy atom
    ca = gemmi.Atom()
    ca.name = "CA"
    ca.element = gemmi.Element("C")
    ca.pos = gemmi.Position(0.0, 0.0, 0.0)
    residue.add_atom(ca)

    # Hydrogen
    ha = gemmi.Atom()
    ha.name = "HA"
    ha.element = gemmi.Element("H")
    ha.pos = gemmi.Position(1.0, 0.0, 0.0)
    residue.add_atom(ha)

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)

    return structure


# =============================================================================
# Tests for extract_atoms_from_model
# =============================================================================


class TestExtractAtoms:
    """Tests for extract_atoms_from_model function."""

    def test_extract_basic(self, simple_structure):
        """Should extract all atoms from a simple structure."""
        atoms = extract_atoms_from_model(simple_structure[0])

        assert isinstance(atoms, AtomData)
        assert len(atoms) == 5
        assert atoms.coords.shape == (5, 3)
        assert atoms.residue_names == ["ALA"] * 5
        assert set(atoms.atom_names) == {"N", "CA", "C", "O", "CB"}
        assert atoms.chain_ids == ["A"] * 5

    def test_extract_hetatm_included(self, structure_with_hetatm):
        """Should include HETATM by default."""
        atoms = extract_atoms_from_model(structure_with_hetatm[0])
        assert len(atoms) == 2  # CA + water O

    def test_extract_hetatm_excluded(self, structure_with_hetatm):
        """Should exclude HETATM when requested."""
        atoms = extract_atoms_from_model(structure_with_hetatm[0], include_hetatm=False)
        assert len(atoms) == 1  # Only CA

    def test_extract_hydrogens_excluded(self, structure_with_hydrogens):
        """Should exclude hydrogens by default."""
        atoms = extract_atoms_from_model(structure_with_hydrogens[0])
        assert len(atoms) == 1  # Only CA, no HA

    def test_extract_hydrogens_included(self, structure_with_hydrogens):
        """Should include hydrogens when requested."""
        atoms = extract_atoms_from_model(structure_with_hydrogens[0], include_hydrogens=True)
        assert len(atoms) == 2  # CA + HA

    def test_atom_data_repr(self, simple_structure):
        """AtomData should have a clean repr."""
        atoms = extract_atoms_from_model(simple_structure[0])
        assert "n_atoms=5" in repr(atoms)


# =============================================================================
# Tests for calculate_sasa_from_model
# =============================================================================


class TestCalculateSasaFromModel:
    """Tests for calculate_sasa_from_model function."""

    def test_basic_calculation(self, simple_structure):
        """Should calculate SASA from a model."""
        result = calculate_sasa_from_model(simple_structure[0])

        assert isinstance(result, SasaResultWithAtoms)
        assert result.total_area > 0
        assert len(result.atom_areas) == 5
        assert result.polar_area >= 0
        assert result.apolar_area >= 0
        assert result.polar_area + result.apolar_area <= result.total_area * 1.01  # Allow 1% error

    def test_result_has_atom_data(self, simple_structure):
        """Result should include atom metadata."""
        result = calculate_sasa_from_model(simple_structure[0])

        assert result.atom_data is not None
        assert len(result.atom_data) == 5
        assert result.atom_classes is not None
        assert len(result.atom_classes) == 5

    def test_polar_apolar_classification(self, simple_structure):
        """Should correctly classify polar and apolar atoms."""
        result = calculate_sasa_from_model(simple_structure[0])

        # N and O should be polar, C atoms should be apolar
        polar_count = np.sum(result.atom_classes == AtomClass.POLAR)
        apolar_count = np.sum(result.atom_classes == AtomClass.APOLAR)

        assert polar_count >= 2  # At least N and O
        assert apolar_count >= 2  # At least CA, C, CB

    def test_different_classifiers(self, simple_structure):
        """Should work with different classifiers."""
        result_naccess = calculate_sasa_from_model(
            simple_structure[0], classifier=ClassifierType.NACCESS
        )
        result_protor = calculate_sasa_from_model(
            simple_structure[0], classifier=ClassifierType.PROTOR
        )

        # Results should be similar but not identical
        assert result_naccess.total_area > 0
        assert result_protor.total_area > 0

    def test_different_algorithms(self, simple_structure):
        """Should work with both SR and LR algorithms."""
        result_sr = calculate_sasa_from_model(simple_structure[0], algorithm="sr")
        result_lr = calculate_sasa_from_model(simple_structure[0], algorithm="lr")

        # Results should be similar (within 10%)
        diff = abs(result_sr.total_area - result_lr.total_area)
        assert diff / result_sr.total_area < 0.1

    def test_result_repr(self, simple_structure):
        """SasaResultWithAtoms should have informative repr."""
        result = calculate_sasa_from_model(simple_structure[0])
        repr_str = repr(result)

        assert "total=" in repr_str
        assert "polar=" in repr_str
        assert "apolar=" in repr_str
        assert "n_atoms=" in repr_str


# =============================================================================
# Tests for calculate_sasa_from_structure
# =============================================================================


class TestCalculateSasaFromStructure:
    """Tests for calculate_sasa_from_structure function."""

    def test_from_gemmi_structure(self, simple_structure):
        """Should accept a gemmi Structure object."""
        result = calculate_sasa_from_structure(simple_structure)

        assert isinstance(result, SasaResultWithAtoms)
        assert result.total_area > 0

    def test_model_index(self, simple_structure):
        """Should use specified model index."""
        result = calculate_sasa_from_structure(simple_structure, model_index=0)
        assert result.total_area > 0

    def test_invalid_model_index(self, simple_structure):
        """Should raise error for invalid model index."""
        with pytest.raises(IndexError, match="out of range"):
            calculate_sasa_from_structure(simple_structure, model_index=99)


# =============================================================================
# Tests with real structure files
# =============================================================================


class TestWithRealFiles:
    """Tests using real structure files from the test directory."""

    @pytest.fixture
    def test_cif_path(self) -> Path | None:
        """Get path to test CIF file if it exists."""
        # Look for test files in common locations
        possible_paths = [
            Path(__file__).parent.parent.parent / "benchmarks" / "dataset" / "1aon.json",
            Path(__file__).parent.parent.parent / "test" / "1ubq.cif",
        ]
        for path in possible_paths:
            if path.exists():
                return path
        return None

    def test_from_file_path(self, test_cif_path, tmp_path):
        """Should load structure from file path."""
        if test_cif_path is None:
            pytest.skip("No test structure file available")

        # For now, skip this test as we don't have a CIF file in the test directory
        pytest.skip("No CIF test file available")


# =============================================================================
# Edge cases
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases."""

    def test_empty_model(self):
        """Should handle empty model."""
        structure = gemmi.Structure()
        model = gemmi.Model("1")
        structure.add_model(model)

        result = calculate_sasa_from_model(model)

        assert result.total_area == 0.0
        assert len(result.atom_areas) == 0
        assert result.polar_area == 0.0
        assert result.apolar_area == 0.0

    def test_single_atom(self):
        """Should handle single atom."""
        structure = gemmi.Structure()
        model = gemmi.Model("1")
        chain = gemmi.Chain("A")
        residue = gemmi.Residue()
        residue.name = "ALA"
        residue.seqid = gemmi.SeqId("1")

        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(0.0, 0.0, 0.0)
        residue.add_atom(atom)

        chain.add_residue(residue)
        model.add_chain(chain)
        structure.add_model(model)

        result = calculate_sasa_from_model(model)

        assert result.total_area > 0
        assert len(result.atom_areas) == 1

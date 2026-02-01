"""Tests for Biotite integration.

These tests require Biotite to be installed.
Run with: pip install zsasa[biotite] pytest
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip all tests if Biotite is not installed
biotite = pytest.importorskip("biotite")
import biotite.structure as struc  # noqa: E402

from zsasa import AtomClass, ClassifierType  # noqa: E402
from zsasa.integrations._types import AtomData, SasaResultWithAtoms  # noqa: E402
from zsasa.integrations.biotite import (  # noqa: E402
    calculate_sasa_from_atom_array,
    calculate_sasa_from_structure,
    extract_atoms_from_atom_array,
)

# =============================================================================
# Helper functions for creating test structures
# =============================================================================


def _create_atom_array(
    atoms: list[tuple[str, str, str, int, tuple[float, float, float], str, bool]],
) -> struc.AtomArray:
    """Create a Biotite AtomArray from a list of atom specifications.

    Args:
        atoms: List of (chain_id, res_name, atom_name, res_id, coord, element, hetero) tuples.

    Returns:
        AtomArray with the specified atoms.
    """
    n_atoms = len(atoms)
    atom_array = struc.AtomArray(n_atoms)

    for i, (chain_id, res_name, atom_name, res_id, coord, element, hetero) in enumerate(atoms):
        atom_array.chain_id[i] = chain_id
        atom_array.res_name[i] = res_name
        atom_array.atom_name[i] = atom_name
        atom_array.res_id[i] = res_id
        atom_array.coord[i] = coord
        atom_array.element[i] = element
        atom_array.hetero[i] = hetero

    return atom_array


# =============================================================================
# Test Fixtures
# =============================================================================


@pytest.fixture
def simple_structure() -> struc.AtomArray:
    """Create a simple structure with a few atoms for testing."""
    atoms = [
        ("A", "ALA", "N", 1, (0.0, 0.0, 0.0), "N", False),
        ("A", "ALA", "CA", 1, (1.5, 0.0, 0.0), "C", False),
        ("A", "ALA", "C", 1, (2.5, 1.2, 0.0), "C", False),
        ("A", "ALA", "O", 1, (2.3, 2.4, 0.0), "O", False),
        ("A", "ALA", "CB", 1, (1.8, -1.5, 0.0), "C", False),
    ]
    return _create_atom_array(atoms)


@pytest.fixture
def structure_with_hetatm() -> struc.AtomArray:
    """Create a structure with HETATM records (water)."""
    atoms = [
        ("A", "ALA", "CA", 1, (0.0, 0.0, 0.0), "C", False),
        ("A", "HOH", "O", 100, (5.0, 0.0, 0.0), "O", True),
    ]
    return _create_atom_array(atoms)


@pytest.fixture
def structure_with_hydrogens() -> struc.AtomArray:
    """Create a structure with hydrogen atoms."""
    atoms = [
        ("A", "ALA", "CA", 1, (0.0, 0.0, 0.0), "C", False),
        ("A", "ALA", "HA", 1, (1.0, 0.0, 0.0), "H", False),
    ]
    return _create_atom_array(atoms)


@pytest.fixture
def multi_chain_structure() -> struc.AtomArray:
    """Create a structure with multiple chains."""
    atoms = [
        ("A", "ALA", "CA", 1, (0.0, 0.0, 0.0), "C", False),
        ("A", "ALA", "N", 1, (1.5, 0.0, 0.0), "N", False),
        ("B", "ALA", "CA", 1, (20.0, 0.0, 0.0), "C", False),
        ("B", "ALA", "N", 1, (21.5, 0.0, 0.0), "N", False),
    ]
    return _create_atom_array(atoms)


# =============================================================================
# Tests for extract_atoms_from_atom_array
# =============================================================================


class TestExtractAtoms:
    """Tests for extract_atoms_from_atom_array function."""

    def test_extract_basic(self, simple_structure):
        """Should extract all atoms from a simple structure."""
        atoms = extract_atoms_from_atom_array(simple_structure)

        assert isinstance(atoms, AtomData)
        assert len(atoms) == 5
        assert atoms.coords.shape == (5, 3)
        assert atoms.residue_names == ["ALA"] * 5
        assert set(atoms.atom_names) == {"N", "CA", "C", "O", "CB"}
        assert atoms.chain_ids == ["A"] * 5

    def test_extract_hetatm_included(self, structure_with_hetatm):
        """Should include HETATM by default."""
        atoms = extract_atoms_from_atom_array(structure_with_hetatm)
        assert len(atoms) == 2  # CA + water O

    def test_extract_hetatm_excluded(self, structure_with_hetatm):
        """Should exclude HETATM when requested."""
        atoms = extract_atoms_from_atom_array(structure_with_hetatm, include_hetatm=False)
        assert len(atoms) == 1  # Only CA

    def test_extract_hydrogens_excluded(self, structure_with_hydrogens):
        """Should exclude hydrogens by default."""
        atoms = extract_atoms_from_atom_array(structure_with_hydrogens)
        assert len(atoms) == 1  # Only CA, no HA

    def test_extract_hydrogens_included(self, structure_with_hydrogens):
        """Should include hydrogens when requested."""
        atoms = extract_atoms_from_atom_array(structure_with_hydrogens, include_hydrogens=True)
        assert len(atoms) == 2  # CA + HA

    def test_atom_data_repr(self, simple_structure):
        """AtomData should have a clean repr."""
        atoms = extract_atoms_from_atom_array(simple_structure)
        assert "n_atoms=5" in repr(atoms)


# =============================================================================
# Tests for calculate_sasa_from_atom_array
# =============================================================================


class TestCalculateSasaFromAtomArray:
    """Tests for calculate_sasa_from_atom_array function."""

    def test_basic_calculation(self, simple_structure):
        """Should calculate SASA from an atom array."""
        result = calculate_sasa_from_atom_array(simple_structure)

        assert isinstance(result, SasaResultWithAtoms)
        assert result.total_area > 0
        assert len(result.atom_areas) == 5
        assert result.polar_area >= 0
        assert result.apolar_area >= 0
        assert result.polar_area + result.apolar_area <= result.total_area * 1.01  # Allow 1% error

    def test_result_has_atom_data(self, simple_structure):
        """Result should include atom metadata."""
        result = calculate_sasa_from_atom_array(simple_structure)

        assert result.atom_data is not None
        assert len(result.atom_data) == 5
        assert result.atom_classes is not None
        assert len(result.atom_classes) == 5

    def test_polar_apolar_classification(self, simple_structure):
        """Should correctly classify polar and apolar atoms."""
        result = calculate_sasa_from_atom_array(simple_structure)

        # N and O should be polar, C atoms should be apolar
        polar_count = np.sum(result.atom_classes == AtomClass.POLAR)
        apolar_count = np.sum(result.atom_classes == AtomClass.APOLAR)

        assert polar_count >= 2  # At least N and O
        assert apolar_count >= 2  # At least CA, C, CB

    def test_different_classifiers(self, simple_structure):
        """Should work with different classifiers."""
        result_naccess = calculate_sasa_from_atom_array(
            simple_structure, classifier=ClassifierType.NACCESS
        )
        result_protor = calculate_sasa_from_atom_array(
            simple_structure, classifier=ClassifierType.PROTOR
        )

        # Results should be similar but not identical
        assert result_naccess.total_area > 0
        assert result_protor.total_area > 0

    def test_different_algorithms(self, simple_structure):
        """Should work with both SR and LR algorithms."""
        result_sr = calculate_sasa_from_atom_array(simple_structure, algorithm="sr")
        result_lr = calculate_sasa_from_atom_array(simple_structure, algorithm="lr")

        # Results should be similar (within 10%)
        diff = abs(result_sr.total_area - result_lr.total_area)
        assert diff / result_sr.total_area < 0.1

    def test_result_repr(self, simple_structure):
        """SasaResultWithAtoms should have informative repr."""
        result = calculate_sasa_from_atom_array(simple_structure)
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

    def test_from_biotite_atom_array(self, simple_structure):
        """Should accept a Biotite AtomArray object."""
        result = calculate_sasa_from_structure(simple_structure)

        assert isinstance(result, SasaResultWithAtoms)
        assert result.total_area > 0

    def test_model_index_single_model(self, simple_structure):
        """Should use model_index=0 for single AtomArray."""
        result = calculate_sasa_from_structure(simple_structure, model_index=0)
        assert result.total_area > 0

    def test_invalid_model_index_single_model(self, simple_structure):
        """Should raise error for invalid model index on single AtomArray."""
        with pytest.raises(IndexError, match="out of range"):
            calculate_sasa_from_structure(simple_structure, model_index=1)

    def test_negative_model_index(self, simple_structure):
        """Should raise error for negative model index."""
        with pytest.raises(IndexError, match="out of range"):
            calculate_sasa_from_structure(simple_structure, model_index=-1)

    def test_atom_array_stack(self, simple_structure):
        """Should handle AtomArrayStack (multiple models)."""
        # Create a stack with 2 models
        stack = struc.stack([simple_structure, simple_structure])

        result = calculate_sasa_from_structure(stack, model_index=0)
        assert result.total_area > 0

        result1 = calculate_sasa_from_structure(stack, model_index=1)
        assert result1.total_area > 0

    def test_invalid_model_index_stack(self, simple_structure):
        """Should raise error for invalid model index on AtomArrayStack."""
        stack = struc.stack([simple_structure, simple_structure])

        with pytest.raises(IndexError, match="out of range"):
            calculate_sasa_from_structure(stack, model_index=99)


# =============================================================================
# Tests with real structure files
# =============================================================================


class TestWithRealFiles:
    """Tests using real structure files from the test directory."""

    def test_from_file_path(self, simple_structure, tmp_path):
        """Should load structure from file path."""
        import biotite.structure.io.pdb as pdb

        # Write structure to temp file
        pdb_path = tmp_path / "test.pdb"
        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(simple_structure)
        pdb_file.write(str(pdb_path))

        result = calculate_sasa_from_structure(pdb_path)
        assert result.total_area > 0
        assert len(result.atom_areas) == 5

    def test_from_string_path(self, simple_structure, tmp_path):
        """Should accept string path as well as Path object."""
        import biotite.structure.io.pdb as pdb

        pdb_path = tmp_path / "test.pdb"
        pdb_file = pdb.PDBFile()
        pdb_file.set_structure(simple_structure)
        pdb_file.write(str(pdb_path))

        # Pass as string instead of Path
        result = calculate_sasa_from_structure(str(pdb_path))
        assert result.total_area > 0

    def test_file_not_found(self):
        """Should raise FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError, match="Structure file not found"):
            calculate_sasa_from_structure("/nonexistent/path/to/file.pdb")


# =============================================================================
# Edge cases
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases."""

    def test_empty_atom_array(self):
        """Should handle empty atom array."""
        atom_array = struc.AtomArray(0)

        result = calculate_sasa_from_atom_array(atom_array)

        assert result.total_area == 0.0
        assert len(result.atom_areas) == 0
        assert result.polar_area == 0.0
        assert result.apolar_area == 0.0

    def test_single_atom(self):
        """Should handle single atom."""
        atoms = [("A", "ALA", "CA", 1, (0.0, 0.0, 0.0), "C", False)]
        atom_array = _create_atom_array(atoms)

        result = calculate_sasa_from_atom_array(atom_array)

        assert result.total_area > 0
        assert len(result.atom_areas) == 1

    def test_multi_chain(self, multi_chain_structure):
        """Should handle multiple chains correctly."""
        result = calculate_sasa_from_atom_array(multi_chain_structure)

        assert result.total_area > 0
        # Check that we have atoms from both chains
        chain_ids = set(result.atom_data.chain_ids)
        assert "A" in chain_ids
        assert "B" in chain_ids

    def test_deuterium_filtered(self):
        """Should filter deuterium atoms as hydrogen."""
        atoms = [
            ("A", "ALA", "CA", 1, (0.0, 0.0, 0.0), "C", False),
            ("A", "ALA", "DA", 1, (1.0, 0.0, 0.0), "D", False),  # Deuterium
        ]
        atom_array = _create_atom_array(atoms)

        atoms_extracted = extract_atoms_from_atom_array(atom_array)
        assert len(atoms_extracted) == 1  # Only CA, no DA


# =============================================================================
# Integration with analysis module
# =============================================================================


class TestAnalysisIntegration:
    """Tests for integration with the analysis module."""

    def test_aggregate_from_result(self, simple_structure):
        """Should work with aggregate_from_result."""
        from zsasa.analysis import aggregate_from_result

        result = calculate_sasa_from_atom_array(simple_structure)
        residues = aggregate_from_result(result)

        assert len(residues) == 1
        assert residues[0].residue_name == "ALA"
        assert residues[0].n_atoms == 5
        assert residues[0].rsa is not None

    def test_multi_residue_aggregation(self):
        """Should correctly aggregate multiple residues."""
        from zsasa.analysis import aggregate_from_result

        atoms = [
            ("A", "ALA", "CA", 1, (0.0, 0.0, 0.0), "C", False),
            ("A", "ALA", "N", 1, (1.5, 0.0, 0.0), "N", False),
            ("A", "GLY", "CA", 2, (5.0, 0.0, 0.0), "C", False),
            ("A", "GLY", "N", 2, (6.5, 0.0, 0.0), "N", False),
        ]
        atom_array = _create_atom_array(atoms)

        result = calculate_sasa_from_atom_array(atom_array)
        residues = aggregate_from_result(result)

        assert len(residues) == 2
        assert residues[0].residue_name == "ALA"
        assert residues[0].residue_id == 1
        assert residues[1].residue_name == "GLY"
        assert residues[1].residue_id == 2

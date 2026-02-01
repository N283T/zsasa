"""Tests for BioPython integration.

These tests require BioPython to be installed.
Run with: pip install zsasa[biopython] pytest
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip all tests if BioPython is not installed
Bio = pytest.importorskip("Bio")
from Bio.PDB import StructureBuilder  # noqa: E402

from zsasa import AtomClass, ClassifierType  # noqa: E402
from zsasa.integrations._types import AtomData, SasaResultWithAtoms  # noqa: E402
from zsasa.integrations.biopython import (  # noqa: E402
    calculate_sasa_from_model,
    calculate_sasa_from_structure,
    extract_atoms_from_model,
)

# =============================================================================
# Helper functions for creating test structures
# =============================================================================


def _create_structure(structure_id: str = "test"):
    """Create a StructureBuilder for building test structures."""
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure(structure_id)
    return sb


def _add_model(sb: StructureBuilder.StructureBuilder, model_id: int = 0):
    """Add a model to the structure."""
    sb.init_model(model_id)


def _add_chain(sb: StructureBuilder.StructureBuilder, chain_id: str):
    """Add a chain to the current model."""
    sb.init_chain(chain_id)


def _add_residue(
    sb: StructureBuilder.StructureBuilder,
    resname: str,
    resseq: int,
    hetflag: str = " ",
    icode: str = " ",
):
    """Add a residue to the current chain."""
    sb.init_residue(resname, hetflag, resseq, icode)


def _add_atom(
    sb: StructureBuilder.StructureBuilder,
    name: str,
    coord: tuple[float, float, float],
    element: str,
    serial_number: int = 1,
):
    """Add an atom to the current residue."""
    sb.init_atom(
        name=name,
        coord=np.array(coord, dtype=np.float32),
        b_factor=0.0,
        occupancy=1.0,
        altloc=" ",
        fullname=f" {name} " if len(name) < 4 else name,
        serial_number=serial_number,
        element=element,
    )


# =============================================================================
# Test Fixtures
# =============================================================================


@pytest.fixture
def simple_structure():
    """Create a simple structure with a few atoms for testing."""
    sb = _create_structure("simple")
    _add_model(sb)
    _add_chain(sb, "A")
    _add_residue(sb, "ALA", 1)

    # Add backbone atoms
    atoms = [
        ("N", (0.0, 0.0, 0.0), "N"),
        ("CA", (1.5, 0.0, 0.0), "C"),
        ("C", (2.5, 1.2, 0.0), "C"),
        ("O", (2.3, 2.4, 0.0), "O"),
        ("CB", (1.8, -1.5, 0.0), "C"),
    ]
    for i, (name, pos, elem) in enumerate(atoms, 1):
        _add_atom(sb, name, pos, elem, i)

    return sb.get_structure()


@pytest.fixture
def structure_with_hetatm():
    """Create a structure with HETATM records (water)."""
    sb = _create_structure("hetatm")
    _add_model(sb)
    _add_chain(sb, "A")

    # Protein residue
    _add_residue(sb, "ALA", 1, hetflag=" ")
    _add_atom(sb, "CA", (0.0, 0.0, 0.0), "C", 1)

    # Water (HETATM)
    _add_residue(sb, "HOH", 100, hetflag="W")
    _add_atom(sb, "O", (5.0, 0.0, 0.0), "O", 2)

    return sb.get_structure()


@pytest.fixture
def structure_with_hydrogens():
    """Create a structure with hydrogen atoms."""
    sb = _create_structure("hydrogens")
    _add_model(sb)
    _add_chain(sb, "A")
    _add_residue(sb, "ALA", 1)

    # Heavy atom
    _add_atom(sb, "CA", (0.0, 0.0, 0.0), "C", 1)

    # Hydrogen
    _add_atom(sb, "HA", (1.0, 0.0, 0.0), "H", 2)

    return sb.get_structure()


@pytest.fixture
def multi_chain_structure():
    """Create a structure with multiple chains."""
    sb = _create_structure("multichain")
    _add_model(sb)

    for chain_id in ["A", "B"]:
        _add_chain(sb, chain_id)
        _add_residue(sb, "ALA", 1)
        offset = 0.0 if chain_id == "A" else 20.0
        _add_atom(sb, "CA", (offset, 0.0, 0.0), "C", 1)
        _add_atom(sb, "N", (offset + 1.5, 0.0, 0.0), "N", 2)

    return sb.get_structure()


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

    def test_from_biopython_structure(self, simple_structure):
        """Should accept a BioPython Structure object."""
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

    def test_negative_model_index(self, simple_structure):
        """Should raise error for negative model index."""
        with pytest.raises(IndexError, match="out of range"):
            calculate_sasa_from_structure(simple_structure, model_index=-1)


# =============================================================================
# Tests with real structure files
# =============================================================================


class TestWithRealFiles:
    """Tests using real structure files from the test directory."""

    def test_from_file_path(self, simple_structure, tmp_path):
        """Should load structure from file path."""
        from Bio.PDB import PDBIO

        # Write structure to temp file
        pdb_path = tmp_path / "test.pdb"
        io = PDBIO()
        io.set_structure(simple_structure)
        io.save(str(pdb_path))

        result = calculate_sasa_from_structure(pdb_path)
        assert result.total_area > 0
        assert len(result.atom_areas) == 5

    def test_from_string_path(self, simple_structure, tmp_path):
        """Should accept string path as well as Path object."""
        from Bio.PDB import PDBIO

        pdb_path = tmp_path / "test.pdb"
        io = PDBIO()
        io.set_structure(simple_structure)
        io.save(str(pdb_path))

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

    def test_empty_model(self):
        """Should handle empty model."""
        sb = _create_structure("empty")
        _add_model(sb)
        structure = sb.get_structure()

        result = calculate_sasa_from_model(structure[0])

        assert result.total_area == 0.0
        assert len(result.atom_areas) == 0
        assert result.polar_area == 0.0
        assert result.apolar_area == 0.0

    def test_single_atom(self):
        """Should handle single atom."""
        sb = _create_structure("single")
        _add_model(sb)
        _add_chain(sb, "A")
        _add_residue(sb, "ALA", 1)
        _add_atom(sb, "CA", (0.0, 0.0, 0.0), "C", 1)

        structure = sb.get_structure()
        result = calculate_sasa_from_model(structure[0])

        assert result.total_area > 0
        assert len(result.atom_areas) == 1

    def test_multi_chain(self, multi_chain_structure):
        """Should handle multiple chains correctly."""
        result = calculate_sasa_from_model(multi_chain_structure[0])

        assert result.total_area > 0
        # Check that we have atoms from both chains
        chain_ids = set(result.atom_data.chain_ids)
        assert "A" in chain_ids
        assert "B" in chain_ids


# =============================================================================
# Integration with analysis module
# =============================================================================


class TestAnalysisIntegration:
    """Tests for integration with the analysis module."""

    def test_aggregate_from_result(self, simple_structure):
        """Should work with aggregate_from_result."""
        from zsasa.analysis import aggregate_from_result

        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        assert len(residues) == 1
        assert residues[0].residue_name == "ALA"
        assert residues[0].n_atoms == 5
        assert residues[0].rsa is not None

    def test_multi_residue_aggregation(self):
        """Should correctly aggregate multiple residues."""
        from zsasa.analysis import aggregate_from_result

        sb = _create_structure("multi_res")
        _add_model(sb)
        _add_chain(sb, "A")

        # Add two residues
        _add_residue(sb, "ALA", 1)
        _add_atom(sb, "CA", (0.0, 0.0, 0.0), "C", 1)
        _add_atom(sb, "N", (1.5, 0.0, 0.0), "N", 2)

        _add_residue(sb, "GLY", 2)
        _add_atom(sb, "CA", (5.0, 0.0, 0.0), "C", 3)
        _add_atom(sb, "N", (6.5, 0.0, 0.0), "N", 4)

        structure = sb.get_structure()
        result = calculate_sasa_from_model(structure[0])
        residues = aggregate_from_result(result)

        assert len(residues) == 2
        assert residues[0].residue_name == "ALA"
        assert residues[0].residue_id == 1
        assert residues[1].residue_name == "GLY"
        assert residues[1].residue_id == 2

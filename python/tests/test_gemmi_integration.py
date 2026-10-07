"""Tests for gemmi integration.

These tests require gemmi to be installed.
Run with: pip install zsasa[gemmi] pytest
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip all tests if gemmi is not installed
gemmi = pytest.importorskip("gemmi")

from zsasa import AtomClass, ClassifierType  # noqa: E402
from zsasa.integrations.gemmi import (  # noqa: E402
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
        """Should include HETATM when requested."""
        atoms = extract_atoms_from_model(structure_with_hetatm[0], include_hetatm=True)
        assert len(atoms) == 2  # CA + water O

    def test_extract_hetatm_excluded(self, structure_with_hetatm):
        """Should exclude HETATM by default."""
        atoms = extract_atoms_from_model(structure_with_hetatm[0])
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

    def test_extract_deuterium_excluded(self, structure_with_hydrogens):
        """Deuterium is an isotope of hydrogen and is excluded with it."""
        residue = structure_with_hydrogens[0]["A"][0]
        da = gemmi.Atom()
        da.name = "DA"
        da.element = gemmi.Element("D")
        da.pos = gemmi.Position(0.0, 1.0, 0.0)
        residue.add_atom(da)

        atoms = extract_atoms_from_model(structure_with_hydrogens[0])
        assert atoms.atom_names == ["CA"]

        with_hydrogens = extract_atoms_from_model(
            structure_with_hydrogens[0], include_hydrogens=True
        )
        assert with_hydrogens.atom_names == ["CA", "HA", "DA"]


# =============================================================================
# Tests for alternate locations
# =============================================================================

# Residue 2 is PRO as altloc A and SER as altloc B, with a shared N that has no
# altloc. Residue 3 is LEU as A and ILE as B at equal occupancy. The x
# coordinate is the atom serial number. The same file is a test fixture of the
# PDB parser in src/pdb_parser.zig.
MICROHETEROGENEITY_PDB = """\
ATOM      1  N   GLY A   1       1.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       2.000   0.000   0.000  1.00 10.00           C
ATOM      3  N   PRO A   2       3.000   0.000   0.000  1.00 10.00           N
ATOM      4  CA APRO A   2       4.000   0.000   0.000  0.40 10.00           C
ATOM      8  CA BSER A   2       8.000   0.000   0.000  0.60 10.00           C
ATOM      5  CB APRO A   2       5.000   0.000   0.000  0.40 10.00           C
ATOM      9  CB BSER A   2       9.000   0.000   0.000  0.60 10.00           C
ATOM      6  CG APRO A   2       6.000   0.000   0.000  0.40 10.00           C
ATOM     10  OG BSER A   2      10.000   0.000   0.000  0.60 10.00           O
ATOM      7  CD APRO A   2       7.000   0.000   0.000  0.40 10.00           C
ATOM     11  N  ALEU A   3      11.000   0.000   0.000  0.50 10.00           N
ATOM     14  N  BILE A   3      14.000   0.000   0.000  0.50 10.00           N
ATOM     12  CA ALEU A   3      12.000   0.000   0.000  0.50 10.00           C
ATOM     15  CA BILE A   3      15.000   0.000   0.000  0.50 10.00           C
ATOM     13  CD1ALEU A   3      13.000   0.000   0.000  0.50 10.00           C
ATOM     16  CG2BILE A   3      16.000   0.000   0.000  0.50 10.00           C
ATOM     17  CD1BILE A   3      17.000   0.000   0.000  0.50 10.00           C
END
"""

# CA of residue 1 has the alternates A and B, CB the alternates B and C without
# an A, and O a 0.50/0.50 tie whose first alternate is C. Water 101 has two
# alternates. The x coordinate is the atom serial number.
ALTLOC_PDB = """\
ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA AALA A   1       2.000   0.000   0.000  0.30 10.00           C
ATOM      3  CA BALA A   1       3.000   0.000   0.000  0.70 10.00           C
ATOM      4  CB BALA A   1       4.000   0.000   0.000  0.40 10.00           C
ATOM      5  CB CALA A   1       5.000   0.000   0.000  0.60 10.00           C
ATOM      6  O  CALA A   1       6.000   0.000   0.000  0.50 10.00           O
ATOM      7  O  BALA A   1       7.000   0.000   0.000  0.50 10.00           O
ATOM      8  DA AALA A   1       8.000   0.000   0.000  0.30 10.00           D
ATOM      9  DA BALA A   1       9.000   0.000   0.000  0.70 10.00           D
HETATM   10  O  AHOH A 101      10.000   0.000   0.000  0.50 10.00           O
HETATM   11  O  BHOH A 101      11.000   0.000   0.000  0.50 10.00           O
END
"""


class TestAlternateLocations:
    """One conformer per site, by the rules of ``--altloc=auto``."""

    def test_keeps_one_conformer_per_site(self):
        """A is preferred, then the highest occupancy, then the first alternate."""
        structure = gemmi.read_pdb_string(ALTLOC_PDB)
        atoms = extract_atoms_from_model(structure[0])

        assert atoms.coords[:, 0].tolist() == [1.0, 2.0, 5.0, 6.0]
        assert atoms.atom_names == ["N", "CA", "CB", "O"]

    def test_filters_apply_before_conformers_are_chosen(self):
        """HETATM and hydrogen alternates are resolved when they are included."""
        structure = gemmi.read_pdb_string(ALTLOC_PDB)
        atoms = extract_atoms_from_model(structure[0], include_hetatm=True, include_hydrogens=True)

        assert atoms.coords[:, 0].tolist() == [1.0, 2.0, 5.0, 6.0, 8.0, 10.0]
        assert atoms.residue_names == ["ALA"] * 5 + ["HOH"]

    def test_microheterogeneity_keeps_one_residue(self):
        """Alternates that are different residues do not both survive."""
        structure = gemmi.read_pdb_string(MICROHETEROGENEITY_PDB)
        atoms = extract_atoms_from_model(structure[0])

        # PRO (altloc A) at 2 with the shared N, and LEU (altloc A) at 3
        assert sorted(atoms.coords[:, 0].tolist()) == [
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            6.0,
            7.0,
            11.0,
            12.0,
            13.0,
        ]
        assert {
            (number, name)
            for number, name in zip(atoms.residue_ids, atoms.residue_names, strict=True)
        } == {(1, "GLY"), (2, "PRO"), (3, "LEU")}

    def test_all_conformers_were_counted_before(self):
        """The fixtures do hold more atoms than are kept."""
        structure = gemmi.read_pdb_string(MICROHETEROGENEITY_PDB)
        assert structure[0].count_atom_sites() == 17

    def test_calculation_uses_the_selected_conformer(self):
        """SASA is calculated for the kept atoms only."""
        structure = gemmi.read_pdb_string(MICROHETEROGENEITY_PDB)
        result = calculate_sasa_from_model(structure[0])

        assert len(result.atom_areas) == 10
        assert len(result.atom_data) == 10

    def test_structure_without_altlocs_is_unchanged(self, simple_structure):
        """Nothing is dropped from a structure without alternate locations."""
        atoms = extract_atoms_from_model(simple_structure[0])
        assert atoms.atom_names == ["N", "CA", "C", "O", "CB"]


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
        # Write structure to temp file
        cif_path = tmp_path / "test.cif"
        simple_structure.make_mmcif_document().write_file(str(cif_path))

        result = calculate_sasa_from_structure(cif_path)
        assert result.total_area > 0
        assert len(result.atom_areas) == 5

    def test_from_string_path(self, simple_structure, tmp_path):
        """Should accept string path as well as Path object."""
        cif_path = tmp_path / "test.cif"
        simple_structure.make_mmcif_document().write_file(str(cif_path))

        # Pass as string instead of Path
        result = calculate_sasa_from_structure(str(cif_path))
        assert result.total_area > 0

    def test_file_not_found(self):
        """Should raise FileNotFoundError for non-existent file."""
        with pytest.raises(FileNotFoundError, match="Structure file not found"):
            calculate_sasa_from_structure("/nonexistent/path/to/file.cif")


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

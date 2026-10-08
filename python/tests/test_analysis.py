"""Tests for per-residue aggregation analysis module."""

from __future__ import annotations

import numpy as np
import pytest

from zsasa import MAX_SASA, AtomClass
from zsasa.analysis import (
    ResidueResult,
    aggregate_by_residue,
    aggregate_from_result,
)

# =============================================================================
# Tests for aggregate_by_residue
# =============================================================================


class TestAggregateByResidue:
    """Tests for aggregate_by_residue function."""

    def test_single_residue(self):
        """Should aggregate all atoms in a single residue."""
        atom_areas = np.array([10.0, 20.0, 30.0])
        chain_ids = ["A", "A", "A"]
        residue_ids = [1, 1, 1]
        residue_names = ["ALA", "ALA", "ALA"]

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 1
        res = results[0]
        assert res.chain_id == "A"
        assert res.residue_id == 1
        assert res.residue_name == "ALA"
        assert res.total_area == 60.0
        assert res.n_atoms == 3

    def test_multiple_residues(self):
        """Should properly group atoms into separate residues."""
        atom_areas = np.array([10.0, 20.0, 15.0, 25.0])
        chain_ids = ["A", "A", "A", "A"]
        residue_ids = [1, 1, 2, 2]
        residue_names = ["ALA", "ALA", "GLY", "GLY"]

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 2
        assert results[0].chain_id == "A"
        assert results[0].residue_id == 1
        assert results[0].residue_name == "ALA"
        assert results[0].total_area == 30.0
        assert results[0].n_atoms == 2

        assert results[1].chain_id == "A"
        assert results[1].residue_id == 2
        assert results[1].residue_name == "GLY"
        assert results[1].total_area == 40.0
        assert results[1].n_atoms == 2

    def test_multi_chain(self):
        """Should distinguish same residue_id in different chains."""
        atom_areas = np.array([10.0, 20.0, 30.0, 40.0])
        chain_ids = ["A", "A", "B", "B"]
        residue_ids = [1, 1, 1, 1]  # Same residue_id in different chains
        residue_names = ["ALA", "ALA", "GLY", "GLY"]

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 2
        assert results[0].chain_id == "A"
        assert results[0].total_area == 30.0

        assert results[1].chain_id == "B"
        assert results[1].total_area == 70.0

    def test_polar_apolar_breakdown(self):
        """Should correctly calculate polar and apolar areas."""
        atom_areas = np.array([10.0, 20.0, 30.0, 5.0])
        chain_ids = ["A", "A", "A", "A"]
        residue_ids = [1, 1, 1, 1]
        residue_names = ["ALA", "ALA", "ALA", "ALA"]
        atom_classes = np.array(
            [
                AtomClass.POLAR,  # 10.0
                AtomClass.APOLAR,  # 20.0
                AtomClass.APOLAR,  # 30.0
                AtomClass.UNKNOWN,  # 5.0 - not counted in polar/apolar
            ],
            dtype=np.int32,
        )

        results = aggregate_by_residue(
            atom_areas, chain_ids, residue_ids, residue_names, atom_classes
        )

        assert len(results) == 1
        res = results[0]
        assert res.total_area == 65.0
        assert res.polar_area == 10.0
        assert res.apolar_area == 50.0
        # UNKNOWN area not included in polar or apolar
        assert res.polar_area + res.apolar_area == 60.0

    def test_rsa_standard_amino_acid(self):
        """Should calculate RSA for standard amino acids."""
        # ALA has MAX_SASA of 129.0
        atom_areas = np.array([64.5])  # Half of max
        chain_ids = ["A"]
        residue_ids = [1]
        residue_names = ["ALA"]

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 1
        assert results[0].rsa == pytest.approx(0.5, rel=1e-3)

    def test_rsa_non_standard_amino_acid(self):
        """Should return None RSA for non-standard amino acids."""
        atom_areas = np.array([100.0])
        chain_ids = ["A"]
        residue_ids = [1]
        residue_names = ["HOH"]  # Water - not a standard amino acid

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 1
        assert results[0].rsa is None

    def test_rsa_greater_than_one(self):
        """RSA can exceed 1.0 for exposed terminal residues."""
        # GLY has MAX_SASA of 104.0
        atom_areas = np.array([150.0])  # More than max
        chain_ids = ["A"]
        residue_ids = [1]
        residue_names = ["GLY"]

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 1
        assert results[0].rsa > 1.0
        assert results[0].rsa == pytest.approx(150.0 / 104.0, rel=1e-3)

    def test_empty_input(self):
        """Should handle empty input gracefully."""
        atom_areas = np.array([], dtype=np.float64)
        chain_ids: list[str] = []
        residue_ids: list[int] = []
        residue_names: list[str] = []

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert results == []

    def test_preserves_order(self):
        """Should preserve the order of residues as they appear."""
        atom_areas = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        chain_ids = ["A", "A", "B", "B", "A", "A"]
        residue_ids = [1, 1, 1, 1, 2, 2]
        residue_names = ["ALA", "ALA", "GLY", "GLY", "VAL", "VAL"]

        results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

        assert len(results) == 3
        # Order should be: A:1, B:1, A:2 (as they appear in input)
        assert results[0].chain_id == "A"
        assert results[0].residue_id == 1
        assert results[1].chain_id == "B"
        assert results[1].residue_id == 1
        assert results[2].chain_id == "A"
        assert results[2].residue_id == 2

    def test_all_standard_amino_acids_have_rsa(self):
        """All standard amino acids should produce valid RSA."""
        for aa in MAX_SASA:
            atom_areas = np.array([50.0])
            chain_ids = ["A"]
            residue_ids = [1]
            residue_names = [aa]

            results = aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

            assert len(results) == 1
            assert results[0].rsa is not None
            assert results[0].rsa == pytest.approx(50.0 / MAX_SASA[aa], rel=1e-6)

    def test_insertion_codes_keep_residues_apart(self):
        """Residues 10, 10A and 10B share a number and differ in the insertion code."""
        atom_areas = np.array([1.0, 2.0, 4.0, 8.0, 16.0])
        chain_ids = ["H"] * 5
        residue_ids = [10] * 5
        residue_names = ["GLY", "GLY", "SER", "THR", "THR"]
        atom_classes = np.array(
            [AtomClass.POLAR, AtomClass.APOLAR, AtomClass.POLAR, AtomClass.POLAR, AtomClass.APOLAR],
            dtype=np.int32,
        )
        # Structure libraries write a blank for "no insertion code"
        insertion_codes = ["", " ", "A", "B", "B"]

        results = aggregate_by_residue(
            atom_areas,
            chain_ids,
            residue_ids,
            residue_names,
            atom_classes,
            insertion_codes=insertion_codes,
        )

        assert [(r.residue_name, r.residue_id, r.insertion_code, r.n_atoms) for r in results] == [
            ("GLY", 10, "", 2),
            ("SER", 10, "A", 1),
            ("THR", 10, "B", 2),
        ]
        assert [r.total_area for r in results] == [3.0, 4.0, 24.0]
        assert [r.polar_area for r in results] == [1.0, 4.0, 8.0]
        assert [r.apolar_area for r in results] == [2.0, 0.0, 16.0]
        assert results[1].rsa == pytest.approx(4.0 / MAX_SASA["SER"])

    def test_without_insertion_codes_residues_have_none(self):
        """Callers that pass no insertion codes get the grouping by chain and number."""
        atom_areas = np.array([1.0, 2.0])
        results = aggregate_by_residue(atom_areas, ["A", "A"], [1, 1], ["ALA", "ALA"])

        assert len(results) == 1
        assert results[0].insertion_code == ""
        assert results[0].n_atoms == 2

    def test_same_insertion_code_in_different_residues(self):
        """The insertion code is part of the key, next to chain and number."""
        atom_areas = np.array([1.0, 2.0, 4.0])
        results = aggregate_by_residue(
            atom_areas,
            ["A", "A", "B"],
            [1, 2, 1],
            ["ALA", "GLY", "SER"],
            insertion_codes=["A", "A", "A"],
        )

        assert [(r.chain_id, r.residue_id, r.insertion_code) for r in results] == [
            ("A", 1, "A"),
            ("A", 2, "A"),
            ("B", 1, "A"),
        ]


class TestAggregateByResidueValidation:
    """Tests for input validation."""

    def test_mismatched_chain_ids_length(self):
        """Should raise error for mismatched chain_ids length."""
        atom_areas = np.array([10.0, 20.0])
        chain_ids = ["A"]  # Wrong length
        residue_ids = [1, 1]
        residue_names = ["ALA", "ALA"]

        with pytest.raises(ValueError, match="chain_ids length"):
            aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

    def test_mismatched_residue_ids_length(self):
        """Should raise error for mismatched residue_ids length."""
        atom_areas = np.array([10.0, 20.0])
        chain_ids = ["A", "A"]
        residue_ids = [1]  # Wrong length
        residue_names = ["ALA", "ALA"]

        with pytest.raises(ValueError, match="residue_ids length"):
            aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

    def test_mismatched_residue_names_length(self):
        """Should raise error for mismatched residue_names length."""
        atom_areas = np.array([10.0, 20.0])
        chain_ids = ["A", "A"]
        residue_ids = [1, 1]
        residue_names = ["ALA"]  # Wrong length

        with pytest.raises(ValueError, match="residue_names length"):
            aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names)

    def test_mismatched_atom_classes_length(self):
        """Should raise error for mismatched atom_classes length."""
        atom_areas = np.array([10.0, 20.0])
        chain_ids = ["A", "A"]
        residue_ids = [1, 1]
        residue_names = ["ALA", "ALA"]
        atom_classes = np.array([AtomClass.POLAR], dtype=np.int32)  # Wrong length

        with pytest.raises(ValueError, match="atom_classes length"):
            aggregate_by_residue(atom_areas, chain_ids, residue_ids, residue_names, atom_classes)

    def test_mismatched_insertion_codes_length(self):
        """Should raise error for mismatched insertion_codes length."""
        atom_areas = np.array([10.0, 20.0])
        chain_ids = ["A", "A"]
        residue_ids = [1, 1]
        residue_names = ["ALA", "ALA"]

        with pytest.raises(ValueError, match="insertion_codes length"):
            aggregate_by_residue(
                atom_areas, chain_ids, residue_ids, residue_names, insertion_codes=["A"]
            )


class TestResidueResultRepr:
    """Tests for ResidueResult repr."""

    def test_repr_with_rsa(self):
        """Should show RSA value in repr."""
        res = ResidueResult(
            chain_id="A",
            residue_id=1,
            residue_name="ALA",
            total_area=64.5,
            polar_area=20.0,
            apolar_area=44.5,
            rsa=0.5,
            n_atoms=5,
        )
        repr_str = repr(res)
        assert "A:ALA1" in repr_str
        assert "total=64.5" in repr_str
        assert "rsa=0.500" in repr_str
        assert "n_atoms=5" in repr_str

    def test_repr_without_rsa(self):
        """Should show None for RSA when not available."""
        res = ResidueResult(
            chain_id="A",
            residue_id=100,
            residue_name="HOH",
            total_area=50.0,
            polar_area=50.0,
            apolar_area=0.0,
            rsa=None,
            n_atoms=1,
        )
        repr_str = repr(res)
        assert "rsa=None" in repr_str

    def test_repr_with_insertion_code(self):
        """Should show the insertion code after the residue number."""
        res = ResidueResult(
            chain_id="H",
            residue_id=10,
            residue_name="SER",
            total_area=30.0,
            polar_area=10.0,
            apolar_area=20.0,
            rsa=0.19,
            n_atoms=6,
            insertion_code="A",
        )
        assert "H:SER10A," in repr(res)

    def test_positional_construction_without_insertion_code(self):
        """The insertion code is the last field and optional."""
        res = ResidueResult("A", 1, "ALA", 64.5, 20.0, 44.5, 0.5, 5)
        assert res.insertion_code == ""
        assert res.n_atoms == 5


class TestAtomDataInsertionCodes:
    """AtomData carries optional per-atom insertion codes."""

    def test_positional_construction_without_insertion_codes(self):
        from zsasa.integrations._types import AtomData

        atoms = AtomData(np.zeros((1, 3)), ["ALA"], ["CA"], ["A"], [1], ["C"])
        assert atoms.insertion_codes is None
        assert len(atoms) == 1

    def test_insertion_codes_field(self):
        from zsasa.integrations._types import AtomData

        atoms = AtomData(
            coords=np.zeros((2, 3)),
            residue_names=["GLY", "SER"],
            atom_names=["CA", "CA"],
            chain_ids=["H", "H"],
            residue_ids=[10, 10],
            elements=["C", "C"],
            insertion_codes=["", "A"],
        )
        assert atoms.insertion_codes == ["", "A"]


# =============================================================================
# Tests for aggregate_from_result (integration with gemmi)
# =============================================================================

# Skip these tests if gemmi is not installed
gemmi = pytest.importorskip("gemmi")

from zsasa.integrations.gemmi import (  # noqa: E402
    calculate_sasa_from_model,
)


@pytest.fixture
def simple_structure() -> gemmi.Structure:
    """Create a simple structure with multiple residues."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")

    # Residue 1: ALA
    res1 = gemmi.Residue()
    res1.name = "ALA"
    res1.seqid = gemmi.SeqId("1")
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
        res1.add_atom(atom)
    chain.add_residue(res1)

    # Residue 2: GLY
    res2 = gemmi.Residue()
    res2.name = "GLY"
    res2.seqid = gemmi.SeqId("2")
    for name, element, pos in [
        ("N", "N", (5.0, 0.0, 0.0)),
        ("CA", "C", (6.5, 0.0, 0.0)),
        ("C", "C", (7.5, 1.2, 0.0)),
        ("O", "O", (7.3, 2.4, 0.0)),
    ]:
        atom = gemmi.Atom()
        atom.name = name
        atom.element = gemmi.Element(element)
        atom.pos = gemmi.Position(*pos)
        res2.add_atom(atom)
    chain.add_residue(res2)

    model.add_chain(chain)
    structure.add_model(model)

    return structure


@pytest.fixture
def multi_chain_structure() -> gemmi.Structure:
    """Create a structure with multiple chains."""
    structure = gemmi.Structure()
    model = gemmi.Model("1")

    for chain_name in ["A", "B"]:
        chain = gemmi.Chain(chain_name)
        res = gemmi.Residue()
        res.name = "ALA"
        res.seqid = gemmi.SeqId("1")

        offset = 0.0 if chain_name == "A" else 20.0
        for name, element, pos in [
            ("N", "N", (0.0 + offset, 0.0, 0.0)),
            ("CA", "C", (1.5 + offset, 0.0, 0.0)),
            ("C", "C", (2.5 + offset, 1.2, 0.0)),
            ("O", "O", (2.3 + offset, 2.4, 0.0)),
        ]:
            atom = gemmi.Atom()
            atom.name = name
            atom.element = gemmi.Element(element)
            atom.pos = gemmi.Position(*pos)
            res.add_atom(atom)
        chain.add_residue(res)
        model.add_chain(chain)

    structure.add_model(model)
    return structure


class TestAggregateFromResult:
    """Tests for aggregate_from_result convenience function."""

    def test_basic_aggregation(self, simple_structure):
        """Should aggregate atoms from SasaResultWithAtoms."""
        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        assert len(residues) == 2
        assert residues[0].chain_id == "A"
        assert residues[0].residue_name == "ALA"
        assert residues[0].residue_id == 1
        assert residues[0].n_atoms == 5

        assert residues[1].chain_id == "A"
        assert residues[1].residue_name == "GLY"
        assert residues[1].residue_id == 2
        assert residues[1].n_atoms == 4

    def test_areas_sum_correctly(self, simple_structure):
        """Residue areas should sum to total atom areas."""
        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        total_from_residues = sum(r.total_area for r in residues)
        total_from_atoms = float(np.sum(result.atom_areas))

        assert total_from_residues == pytest.approx(total_from_atoms, rel=1e-6)

    def test_polar_apolar_passed_through(self, simple_structure):
        """Should include polar/apolar areas from atom classification."""
        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        # Both residues should have some polar and apolar area
        for res in residues:
            assert res.polar_area >= 0
            assert res.apolar_area >= 0
            # polar + apolar should not exceed total (UNKNOWN atoms exist)
            assert res.polar_area + res.apolar_area <= res.total_area * 1.01

    def test_rsa_calculated(self, simple_structure):
        """RSA should be calculated for standard amino acids."""
        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        for res in residues:
            assert res.rsa is not None
            assert res.rsa >= 0

    def test_multi_chain(self, multi_chain_structure):
        """Should distinguish residues in different chains."""
        result = calculate_sasa_from_model(multi_chain_structure[0])
        residues = aggregate_from_result(result)

        assert len(residues) == 2
        chain_ids = [r.chain_id for r in residues]
        assert "A" in chain_ids
        assert "B" in chain_ids

    def test_empty_result(self):
        """Should handle empty structure."""
        structure = gemmi.Structure()
        model = gemmi.Model("1")
        structure.add_model(model)

        result = calculate_sasa_from_model(model)
        residues = aggregate_from_result(result)

        assert residues == []


class TestRealWorldUsage:
    """Tests demonstrating real-world usage patterns."""

    def test_find_buried_residues(self, simple_structure):
        """Demonstrate finding buried residues (RSA < 0.25)."""
        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        buried = [r for r in residues if r.rsa is not None and r.rsa < 0.25]
        exposed = [r for r in residues if r.rsa is not None and r.rsa >= 0.25]

        # Just verify we can filter
        assert len(buried) + len(exposed) == len(residues)

    def test_sum_polar_apolar_areas(self, simple_structure):
        """Demonstrate summing polar/apolar areas across structure."""
        result = calculate_sasa_from_model(simple_structure[0])
        residues = aggregate_from_result(result)

        total_polar = sum(r.polar_area for r in residues)
        total_apolar = sum(r.apolar_area for r in residues)

        # Should be consistent with SasaResultWithAtoms totals (within tolerance)
        assert total_polar == pytest.approx(result.polar_area, rel=0.01)
        assert total_apolar == pytest.approx(result.apolar_area, rel=0.01)

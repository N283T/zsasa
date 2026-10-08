"""Element inference of the MDAnalysis integration when a topology has no element column.

These tests use plain stand-in atoms, so they run without MDAnalysis installed.
"""

from __future__ import annotations

import importlib
import sys
import types
from dataclasses import dataclass

import pytest

from zsasa import mdanalysis as zmda


@dataclass
class FakeAtom:
    """The attributes of an MDAnalysis atom the element inference reads."""

    name: str
    type: str | None = None
    resname: str | None = None
    mass: float | None = None


@dataclass
class FakeAtomWithElement(FakeAtom):
    element: str = ""


class TestTwoLetterElements:
    """Ions and halogens are recognized from their type or name."""

    @pytest.mark.parametrize(
        ("atom", "expected"),
        [
            # Types as the topology parsers guess them (GRO, PDB without element column)
            (FakeAtom("CL", "CL", "CL", 35.45), "Cl"),
            (FakeAtom("FE", "FE", "HEM", 55.847), "Fe"),
            (FakeAtom("ZN", "ZN", "ZN", 65.37), "Zn"),
            (FakeAtom("MG", "MG", "MG", 24.305), "Mg"),
            (FakeAtom("NA", "NA", "NA", 22.98977), "Na"),
            (FakeAtom("BR", "BR", "BR", 79.904), "Br"),
            (FakeAtom("CA", "CA", "CA", 40.08), "Ca"),
            (FakeAtom("CA", "C0", "CA", 40.08), "Ca"),  # AMBER Ca2+, type C0
            (FakeAtom("HG", "HG", "HG", 200.59), "Hg"),
            # Decorated names: numbers and charges
            (FakeAtom("CL1", "CL1", "LIG"), "Cl"),
            (FakeAtom("ZN2+", "ZN2+", "ZN2"), "Zn"),
            (FakeAtom("Na+", "Na+", "Na+"), "Na"),
            (FakeAtom("BR2", None, "LIG"), "Br"),
            # No type at all: the name decides
            (FakeAtom("FE", None, "HEM"), "Fe"),
            (FakeAtom("CL", None, "CL"), "Cl"),
            # CHARMM ion names
            (FakeAtom("SOD", "SOD", "SOD", 22.99), "Na"),
            (FakeAtom("CLA", "CLA", "CLA", 35.45), "Cl"),
            (FakeAtom("POT", "POT", "POT", 39.098), "K"),
            (FakeAtom("CAL", "CAL", "CAL", 40.08), "Ca"),
            # CHARMM force-field types for halogens
            (FakeAtom("CL1", "CLGA1", "LIG", 35.45), "Cl"),
            (FakeAtom("BR1", "BRGR1", "LIG", 79.9), "Br"),
        ],
    )
    def test_ions_and_halogens(self, atom: FakeAtom, expected: str) -> None:
        assert zmda._get_element(atom) == expected

    def test_first_letter_regression(self) -> None:
        """The first character alone (the old behavior) read these as C, F, Z... and N."""
        assert zmda._get_element(FakeAtom("CL", "CL", "CL", 35.45)) != "C"
        assert zmda._get_element(FakeAtom("FE", "FE", "HEM", 55.8)) != "F"
        assert zmda._get_element(FakeAtom("ZN", "ZN", "ZN", 65.4)) != "Z"
        assert zmda._get_element(FakeAtom("NA", "NA", "NA", 22.99)) != "N"

    def test_radius_of_inferred_ions(self) -> None:
        """The inferred element selects the table radius instead of the 2.0 default."""
        assert zmda._get_radius(FakeAtom("CL", "CL", "CL", 35.45)) == pytest.approx(1.75)
        assert zmda._get_radius(FakeAtom("ZN", "ZN", "ZN", 65.37)) == pytest.approx(1.39)
        assert zmda._get_radius(FakeAtom("NA", "NA", "NA", 22.99)) == pytest.approx(2.27)
        assert zmda._get_radius(FakeAtom("CA", "CA", "CA", 40.08)) == pytest.approx(2.31)


class TestOrganicAtomNames:
    """Names that look like an element symbol but are ordinary atom names stay organic."""

    @pytest.mark.parametrize(
        ("atom", "expected"),
        [
            # The alpha carbon of a protein is carbon, with or without a guessed type
            (FakeAtom("CA", "C", "ALA", 12.011), "C"),
            (FakeAtom("CA", "CA", "ALA", 12.011), "C"),
            (FakeAtom("CA", None, "GLY"), "C"),
            (FakeAtom("CA", "CA", "ALA"), "C"),
            # Side-chain and base names that begin with a symbol
            (FakeAtom("CD", "CD", "GLN"), "C"),
            (FakeAtom("CD1", "CD1", "LEU"), "C"),
            (FakeAtom("CE", "CE", "LYS"), "C"),
            (FakeAtom("CO", "CO", "LIG"), "C"),
            (FakeAtom("NA", "NA", "HEM", 22.99), "N"),
            (FakeAtom("NA", None, "HEM"), "N"),
            (FakeAtom("NE", "NE", "ARG"), "N"),
            (FakeAtom("HG", "HG", "SER"), "H"),
            (FakeAtom("HG21", "HG21", "VAL"), "H"),
            (FakeAtom("HE1", None, "HIS"), "H"),
            (FakeAtom("1HB", None, "ALA"), "H"),
            (FakeAtom("SE", "SE", "LIG"), "S"),
            (FakeAtom("O5'", None, "DA"), "O"),
            # Force-field types are not elements
            (FakeAtom("CZ", "CA", "PHE", 12.011), "C"),
            (FakeAtom("CB", "CT1", "ALA", 12.011), "C"),
            (FakeAtom("N", "NH1", "ALA", 14.007), "N"),
            (FakeAtom("OW", "OT", "TIP3", 15.999), "O"),
            (FakeAtom("HW1", "HT", "TIP3", 1.008), "H"),
            (FakeAtom("HA", "HGA1", "ALA", 1.008), "H"),
        ],
    )
    def test_organic_names(self, atom: FakeAtom, expected: str) -> None:
        assert zmda._get_element(atom) == expected

    def test_ion_is_found_even_when_the_parser_guessed_a_carbon_type(self) -> None:
        """MDAnalysis types the ion CA in residue CA as carbon (type "C", mass 12.011)."""
        assert zmda._get_element(FakeAtom("CA", "C", "CA", 12.011)) == "Ca"
        assert zmda._get_element(FakeAtom("NA", "N", "NA", 14.007)) == "Na"
        # The same names in a protein or a heme stay organic.
        assert zmda._get_element(FakeAtom("CA", "C", "ALA", 12.011)) == "C"
        assert zmda._get_element(FakeAtom("NA", "N", "HEM", 14.007)) == "N"

    def test_residue_name_decides_between_ion_and_atom(self) -> None:
        assert zmda._get_element(FakeAtom("CA", "CA", "CA")) == "Ca"
        assert zmda._get_element(FakeAtom("CA", "CA", "ALA")) == "C"
        assert zmda._get_element(FakeAtom("CD", "CD", "CD")) == "Cd"
        assert zmda._get_element(FakeAtom("CD", "CD", "GLN")) == "C"

    def test_protein_radii_are_unchanged(self) -> None:
        carbon = zmda._get_radius(FakeAtom("CA", "C", "ALA", 12.011))
        assert carbon == pytest.approx(1.70)
        assert zmda._get_radius(FakeAtom("N", "N", "ALA", 14.007)) == pytest.approx(1.55)
        assert zmda._get_radius(FakeAtom("HG", "H", "SER", 1.008)) == pytest.approx(1.1)


class TestElementAttribute:
    def test_element_attribute_has_priority(self) -> None:
        atom = FakeAtomWithElement("CA", "CA", "CA", 40.08, element="C")
        assert zmda._get_element(atom) == "C"

    def test_element_attribute_is_normalized(self) -> None:
        atom = FakeAtomWithElement("X", "X", "X", element="CL")
        assert zmda._get_element(atom) == "Cl"

    def test_empty_element_falls_back_to_inference(self) -> None:
        atom = FakeAtomWithElement("ZN", "ZN", "ZN", 65.38, element="")
        assert zmda._get_element(atom) == "Zn"


class TestLastResort:
    def test_nothing_known_is_carbon(self) -> None:
        assert zmda._get_element(FakeAtom("", None, None)) == "C"

    def test_missing_attributes_are_tolerated(self) -> None:
        assert zmda._get_element(types.SimpleNamespace(name="FE")) == "Fe"
        assert zmda._get_element(types.SimpleNamespace(type="CL")) == "Cl"
        assert zmda._get_element(types.SimpleNamespace()) == "C"


class TestVdwTableImport:
    """The radii table comes from guesser.tables (MDAnalysis 2.8+) or topology.tables (older)."""

    @pytest.fixture
    def restore_module(self):  # noqa: ANN201
        yield
        importlib.reload(zmda)

    @staticmethod
    def _fake_tables(marker: dict[str, float]) -> types.ModuleType:
        module = types.ModuleType("fake_tables")
        module.vdwradii = marker
        return module

    def test_old_module_path_is_used_before_mdanalysis_2_8(
        self, monkeypatch: pytest.MonkeyPatch, restore_module: None
    ) -> None:
        old_table = {"C": 9.0}
        monkeypatch.setitem(sys.modules, "MDAnalysis.guesser.tables", None)  # not importable
        monkeypatch.setitem(sys.modules, "MDAnalysis.topology.tables", self._fake_tables(old_table))

        reloaded = importlib.reload(zmda)

        assert reloaded._MDA_VDW_RADII is old_table

    def test_new_module_path_is_preferred(
        self, monkeypatch: pytest.MonkeyPatch, restore_module: None
    ) -> None:
        new_table = {"C": 8.0}
        monkeypatch.setitem(sys.modules, "MDAnalysis.guesser.tables", self._fake_tables(new_table))
        monkeypatch.setitem(
            sys.modules, "MDAnalysis.topology.tables", self._fake_tables({"C": 9.0})
        )

        reloaded = importlib.reload(zmda)

        assert reloaded._MDA_VDW_RADII is new_table

    def test_builtin_table_covers_the_inferred_elements(
        self, monkeypatch: pytest.MonkeyPatch, restore_module: None
    ) -> None:
        """Without MDAnalysis, the elements the inference returns still get a radius."""
        monkeypatch.setitem(sys.modules, "MDAnalysis.guesser.tables", None)
        monkeypatch.setitem(sys.modules, "MDAnalysis.topology.tables", None)

        reloaded = importlib.reload(zmda)

        for symbol in ("H", "C", "N", "O", "S", "P", "CL", "NA", "MG", "ZN", "CA", "BR"):
            assert symbol in reloaded._MDA_VDW_RADII

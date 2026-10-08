"""Integration tests using real PDB/mmCIF structure files.

These tests verify that the full pipeline works correctly with real-world
structures from the examples/ directory.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

# Get the examples directory path
EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def _skip_if_missing(path: Path) -> None:
    """Skip test if file doesn't exist."""
    if not path.exists():
        pytest.skip(f"Test file not found: {path}")


def _assert_matches_cli(result, path: Path, tmp_path: Path) -> None:
    """The integration must give what `zsasa calc` gives for the same file.

    The defaults of the integrations (CCD classifier, 100 points, probe radius 1.4,
    no HETATM, no hydrogens) are spelled out here so that the comparison does not
    depend on the CLI defaults.
    """
    output = tmp_path / "cli.json"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "zsasa.cli",
            "calc",
            str(path),
            str(output),
            "-q",
            "--algorithm=sr",
            "--classifier=ccd",
            "--n-points=100",
            "--probe-radius=1.4",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    expected = json.loads(output.read_text())

    assert len(result.atom_areas) == len(expected["atom_areas"])
    assert result.total_area == pytest.approx(expected["total_area"], rel=1e-9)
    np.testing.assert_allclose(result.atom_areas, expected["atom_areas"], rtol=1e-9, atol=1e-9)


class TestGemmiRealStructures:
    """Test gemmi integration with real structure files."""

    @pytest.fixture(autouse=True)
    def check_gemmi(self):
        """Skip if gemmi is not installed."""
        pytest.importorskip("gemmi")

    def test_1crn_pdb(self, tmp_path):
        """Test SASA calculation for crambin (1CRN) from PDB."""
        from zsasa.integrations.gemmi import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        # Crambin has 327 heavy atoms (no hydrogens, no HETATM)
        assert len(result.atom_areas) == 327
        _assert_matches_cli(result, pdb_path, tmp_path)
        assert result.polar_area > 0
        assert result.apolar_area > 0
        assert result.polar_area + result.apolar_area == pytest.approx(result.total_area)

    def test_1ubq_cif(self, tmp_path):
        """Test SASA calculation for ubiquitin (1UBQ) from mmCIF."""
        from zsasa.integrations.gemmi import calculate_sasa_from_structure

        cif_path = EXAMPLES_DIR / "1ubq.cif"
        _skip_if_missing(cif_path)

        result = calculate_sasa_from_structure(cif_path)

        # Ubiquitin has 602 heavy atoms (no hydrogens, no HETATM)
        assert len(result.atom_areas) == 602
        _assert_matches_cli(result, cif_path, tmp_path)

    def test_1crn_residue_aggregation(self):
        """Test per-residue aggregation for crambin."""
        from zsasa.analysis import aggregate_from_result
        from zsasa.integrations.gemmi import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)
        residues = aggregate_from_result(result)

        # Crambin has 46 residues
        assert len(residues) >= 40
        assert len(residues) <= 50

        # Check RSA values are reasonable
        for res in residues:
            if res.rsa is not None:
                assert 0.0 <= res.rsa <= 2.0  # RSA can exceed 1.0 for exposed residues


class TestBioPythonRealStructures:
    """Test BioPython integration with real structure files."""

    @pytest.fixture(autouse=True)
    def check_biopython(self):
        """Skip if BioPython is not installed."""
        pytest.importorskip("Bio")

    def test_1crn_pdb(self, tmp_path):
        """Test SASA calculation for crambin (1CRN) from PDB."""
        from zsasa.integrations.biopython import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert len(result.atom_areas) == 327
        _assert_matches_cli(result, pdb_path, tmp_path)

    def test_1ubq_pdb(self, tmp_path):
        """Test SASA calculation for ubiquitin (1UBQ) from PDB."""
        from zsasa.integrations.biopython import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1ubq.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert len(result.atom_areas) == 602
        _assert_matches_cli(result, pdb_path, tmp_path)

    def test_1crn_residue_aggregation(self):
        """Test per-residue aggregation for crambin."""
        from zsasa.analysis import aggregate_from_result
        from zsasa.integrations.biopython import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)
        residues = aggregate_from_result(result)

        assert len(residues) >= 40
        assert len(residues) <= 50


class TestBiotiteRealStructures:
    """Test Biotite integration with real structure files."""

    @pytest.fixture(autouse=True)
    def check_biotite(self):
        """Skip if Biotite is not installed."""
        pytest.importorskip("biotite")

    def test_1crn_pdb(self, tmp_path):
        """Test SASA calculation for crambin (1CRN) from PDB."""
        from zsasa.integrations.biotite import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert len(result.atom_areas) == 327
        _assert_matches_cli(result, pdb_path, tmp_path)

    def test_1ubq_pdb(self, tmp_path):
        """Test SASA calculation for ubiquitin (1UBQ) from PDB."""
        from zsasa.integrations.biotite import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1ubq.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert len(result.atom_areas) == 602
        _assert_matches_cli(result, pdb_path, tmp_path)

    def test_1crn_residue_aggregation(self):
        """Test per-residue aggregation for crambin."""
        from zsasa.analysis import aggregate_from_result
        from zsasa.integrations.biotite import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)
        residues = aggregate_from_result(result)

        assert len(residues) >= 40
        assert len(residues) <= 50


class TestCrossLibraryConsistency:
    """Test that different libraries produce consistent results."""

    @pytest.fixture(autouse=True)
    def check_all_libraries(self):
        """Skip if any library is not installed."""
        pytest.importorskip("gemmi")
        pytest.importorskip("Bio")
        pytest.importorskip("biotite")

    def test_1crn_consistency(self):
        """All libraries should produce similar SASA for 1CRN."""
        from zsasa.integrations.biopython import (
            calculate_sasa_from_structure as bp_calc,
        )
        from zsasa.integrations.biotite import (
            calculate_sasa_from_structure as bt_calc,
        )
        from zsasa.integrations.gemmi import (
            calculate_sasa_from_structure as gm_calc,
        )

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result_gm = gm_calc(pdb_path)
        result_bp = bp_calc(pdb_path)
        result_bt = bt_calc(pdb_path)

        # Results should be within 5% of each other
        # (differences due to hydrogen handling, occupancy, etc.)
        mean_area = (result_gm.total_area + result_bp.total_area + result_bt.total_area) / 3

        assert abs(result_gm.total_area - mean_area) / mean_area < 0.05
        assert abs(result_bp.total_area - mean_area) / mean_area < 0.05
        assert abs(result_bt.total_area - mean_area) / mean_area < 0.05


class TestInsertionCodes:
    """Residues that differ only in the insertion code stay separate residues."""

    @staticmethod
    def ubq_with_insertion_code(tmp_path: Path) -> Path:
        """Write 1UBQ with residue 2 renumbered to 1A (so 1 and 1A share a number)."""
        source = EXAMPLES_DIR / "1ubq.pdb"
        _skip_if_missing(source)
        lines = []
        for line in source.read_text().splitlines(keepends=True):
            if line.startswith("ATOM") and line[21] == "A" and line[22:27] == "   2 ":
                line = line[:22] + "   1A" + line[27:]
            lines.append(line)
        path = tmp_path / "1ubq_1A.pdb"
        path.write_text("".join(lines))
        return path

    @pytest.mark.parametrize(
        ("package", "module"),
        [("gemmi", "gemmi"), ("Bio", "biopython"), ("biotite", "biotite")],
    )
    def test_residues_1_and_1a_are_two_residues(self, tmp_path, package, module):
        import importlib

        pytest.importorskip(package)
        from zsasa.analysis import aggregate_from_result

        integration = importlib.import_module(f"zsasa.integrations.{module}")
        result = integration.calculate_sasa_from_structure(self.ubq_with_insertion_code(tmp_path))

        # One insertion code per atom: "" for a residue without one
        codes = result.atom_data.insertion_codes
        assert codes is not None
        assert len(codes) == len(result.atom_areas)
        assert set(codes) == {"", "A"}
        assert codes.count("A") == 9

        residues = aggregate_from_result(result)

        # 76 residues as in the unmodified file, not 75 with MET 1 holding 17 atoms
        assert len(residues) == 76
        met, gln, ile = residues[:3]
        assert (met.residue_name, met.residue_id, met.insertion_code, met.n_atoms) == (
            "MET",
            1,
            "",
            8,
        )
        assert (gln.residue_name, gln.residue_id, gln.insertion_code, gln.n_atoms) == (
            "GLN",
            1,
            "A",
            9,
        )
        assert (ile.residue_name, ile.residue_id, ile.insertion_code, ile.n_atoms) == (
            "ILE",
            3,
            "",
            8,
        )
        assert sum(res.total_area for res in residues) == pytest.approx(result.total_area)

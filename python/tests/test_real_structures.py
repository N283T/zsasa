"""Integration tests using real PDB/mmCIF structure files.

These tests verify that the full pipeline works correctly with real-world
structures from the examples/ directory.
"""

from pathlib import Path

import pytest

# Get the examples directory path
EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def _skip_if_missing(path: Path) -> None:
    """Skip test if file doesn't exist."""
    if not path.exists():
        pytest.skip(f"Test file not found: {path}")


class TestGemmiRealStructures:
    """Test gemmi integration with real structure files."""

    @pytest.fixture(autouse=True)
    def check_gemmi(self):
        """Skip if gemmi is not installed."""
        pytest.importorskip("gemmi")

    def test_1crn_pdb(self):
        """Test SASA calculation for crambin (1CRN) from PDB."""
        from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        # Crambin: ~327 atoms, total SASA should be around 2400-2600 A^2
        assert result.total_area > 2000
        assert result.total_area < 3000
        assert len(result.atom_areas) > 300
        assert result.polar_area > 0
        assert result.apolar_area > 0

    def test_1ubq_cif(self):
        """Test SASA calculation for ubiquitin (1UBQ) from mmCIF."""
        from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure

        cif_path = EXAMPLES_DIR / "1ubq.cif"
        _skip_if_missing(cif_path)

        result = calculate_sasa_from_structure(cif_path)

        # Ubiquitin: ~602 atoms, total SASA should be around 4500-5500 A^2
        assert result.total_area > 4000
        assert result.total_area < 6000
        assert len(result.atom_areas) > 500

    def test_1crn_residue_aggregation(self):
        """Test per-residue aggregation for crambin."""
        from freesasa_zig.analysis import aggregate_from_result
        from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure

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

    def test_1crn_pdb(self):
        """Test SASA calculation for crambin (1CRN) from PDB."""
        from freesasa_zig.integrations.biopython import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert result.total_area > 2000
        assert result.total_area < 3000
        assert len(result.atom_areas) > 300

    def test_1ubq_pdb(self):
        """Test SASA calculation for ubiquitin (1UBQ) from PDB."""
        from freesasa_zig.integrations.biopython import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1ubq.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert result.total_area > 4000
        assert result.total_area < 6000
        assert len(result.atom_areas) > 500

    def test_1crn_residue_aggregation(self):
        """Test per-residue aggregation for crambin."""
        from freesasa_zig.analysis import aggregate_from_result
        from freesasa_zig.integrations.biopython import calculate_sasa_from_structure

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

    def test_1crn_pdb(self):
        """Test SASA calculation for crambin (1CRN) from PDB."""
        from freesasa_zig.integrations.biotite import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1crn.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert result.total_area > 2000
        assert result.total_area < 3000
        assert len(result.atom_areas) > 300

    def test_1ubq_pdb(self):
        """Test SASA calculation for ubiquitin (1UBQ) from PDB."""
        from freesasa_zig.integrations.biotite import calculate_sasa_from_structure

        pdb_path = EXAMPLES_DIR / "1ubq.pdb"
        _skip_if_missing(pdb_path)

        result = calculate_sasa_from_structure(pdb_path)

        assert result.total_area > 4000
        assert result.total_area < 6000
        assert len(result.atom_areas) > 500

    def test_1crn_residue_aggregation(self):
        """Test per-residue aggregation for crambin."""
        from freesasa_zig.analysis import aggregate_from_result
        from freesasa_zig.integrations.biotite import calculate_sasa_from_structure

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
        from freesasa_zig.integrations.biopython import (
            calculate_sasa_from_structure as bp_calc,
        )
        from freesasa_zig.integrations.biotite import (
            calculate_sasa_from_structure as bt_calc,
        )
        from freesasa_zig.integrations.gemmi import (
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

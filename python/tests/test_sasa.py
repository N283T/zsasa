"""Tests for freesasa-zig Python bindings."""

import numpy as np
import pytest

from freesasa_zig import (
    AtomClass,
    ClassificationResult,
    ClassifierType,
    SasaResult,
    calculate_sasa,
    classify_atoms,
    get_atom_class,
    get_radius,
    get_version,
    guess_radius,
    guess_radius_from_atom_name,
)


class TestVersion:
    """Tests for version function."""

    def test_get_version(self):
        version = get_version()
        assert isinstance(version, str)
        assert len(version) > 0
        # Check version format (e.g., "0.0.5")
        parts = version.split(".")
        assert len(parts) >= 2


class TestSingleAtom:
    """Tests for single isolated atom."""

    def test_single_atom_sr(self):
        """Single atom should have full spherical surface area."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        result = calculate_sasa(coords, radii)

        assert isinstance(result, SasaResult)
        # Expected: 4π * (1.5 + 1.4)² ≈ 105.68 Å²
        expected = 4 * np.pi * (1.5 + 1.4) ** 2
        assert abs(result.total_area - expected) < 1.0
        assert len(result.atom_areas) == 1
        assert result.atom_areas[0] == pytest.approx(result.total_area)

    def test_single_atom_lr(self):
        """Single atom with Lee-Richards algorithm."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        result = calculate_sasa(coords, radii, algorithm="lr")

        expected = 4 * np.pi * (1.5 + 1.4) ** 2
        assert abs(result.total_area - expected) < 1.0


class TestTwoAtoms:
    """Tests for two atom systems."""

    def test_overlapping_atoms(self):
        """Two overlapping atoms should have less total area than two isolated."""
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        radii = np.array([1.5, 1.5])

        result = calculate_sasa(coords, radii)

        # Total should be less than 2 * single atom area due to overlap
        single_area = 4 * np.pi * (1.5 + 1.4) ** 2
        assert result.total_area < 2 * single_area
        assert len(result.atom_areas) == 2

    def test_distant_atoms(self):
        """Two distant atoms should have nearly full spherical areas."""
        coords = np.array([[0.0, 0.0, 0.0], [100.0, 0.0, 0.0]])
        radii = np.array([1.5, 1.5])

        result = calculate_sasa(coords, radii)

        # Total should be approximately 2 * single atom area
        single_area = 4 * np.pi * (1.5 + 1.4) ** 2
        assert abs(result.total_area - 2 * single_area) < 1.0


class TestParameters:
    """Tests for parameter variations."""

    def test_n_points_affect_result(self):
        """More test points should give more accurate results."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        result_low = calculate_sasa(coords, radii, n_points=20)
        result_high = calculate_sasa(coords, radii, n_points=500)

        expected = 4 * np.pi * (1.5 + 1.4) ** 2
        # Higher n_points should be closer to expected
        error_low = abs(result_low.total_area - expected)
        error_high = abs(result_high.total_area - expected)
        assert error_high < error_low or error_high < 0.1

    def test_probe_radius(self):
        """Larger probe radius should give larger area."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        result_small = calculate_sasa(coords, radii, probe_radius=1.0)
        result_large = calculate_sasa(coords, radii, probe_radius=2.0)

        assert result_large.total_area > result_small.total_area

    def test_threading(self):
        """Results should be consistent regardless of thread count."""
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        radii = np.array([1.5, 1.5])

        result_1 = calculate_sasa(coords, radii, n_threads=1)
        result_auto = calculate_sasa(coords, radii, n_threads=0)

        assert result_1.total_area == pytest.approx(result_auto.total_area)


class TestAlgorithmComparison:
    """Tests comparing SR and LR algorithms."""

    def test_sr_lr_similar_results(self):
        """SR and LR should give similar results for same structure."""
        coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
        radii = np.array([1.5, 1.5])

        result_sr = calculate_sasa(coords, radii, algorithm="sr", n_points=100)
        result_lr = calculate_sasa(coords, radii, algorithm="lr", n_slices=20)

        # Results should be within 5%
        diff_percent = abs(result_sr.total_area - result_lr.total_area) / result_sr.total_area * 100
        assert diff_percent < 5.0


class TestInputValidation:
    """Tests for input validation."""

    def test_invalid_coords_shape(self):
        """Should raise error for wrong coordinate shape."""
        coords = np.array([0.0, 0.0, 0.0])  # 1D instead of 2D
        radii = np.array([1.5])

        with pytest.raises(ValueError, match="coords must be"):
            calculate_sasa(coords, radii)

    def test_mismatched_sizes(self):
        """Should raise error for mismatched array sizes."""
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        radii = np.array([1.5])  # Only 1 radius for 2 atoms

        with pytest.raises(ValueError, match="radii must be"):
            calculate_sasa(coords, radii)

    def test_invalid_algorithm(self):
        """Should raise error for unknown algorithm."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        with pytest.raises(ValueError, match="Unknown algorithm"):
            calculate_sasa(coords, radii, algorithm="invalid")

    def test_empty_input(self):
        """Should raise error for empty arrays."""
        coords = np.array([]).reshape(0, 3)
        radii = np.array([])

        with pytest.raises((ValueError, RuntimeError)):
            calculate_sasa(coords, radii)

    def test_negative_n_points(self):
        """Should raise error for negative n_points."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        with pytest.raises(ValueError, match="n_points must be positive"):
            calculate_sasa(coords, radii, n_points=-1)

    def test_zero_n_slices(self):
        """Should raise error for zero n_slices."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        with pytest.raises(ValueError, match="n_slices must be positive"):
            calculate_sasa(coords, radii, algorithm="lr", n_slices=0)

    def test_negative_probe_radius(self):
        """Should raise error for negative probe_radius."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([1.5])

        with pytest.raises(ValueError, match="probe_radius must be positive"):
            calculate_sasa(coords, radii, probe_radius=-1.0)

    def test_negative_radii(self):
        """Should raise error for negative radii."""
        coords = np.array([[0.0, 0.0, 0.0]])
        radii = np.array([-1.5])

        with pytest.raises(ValueError, match="non-negative"):
            calculate_sasa(coords, radii)


# =============================================================================
# Classifier Tests
# =============================================================================


class TestClassifierTypes:
    """Tests for ClassifierType enum."""

    def test_classifier_type_values(self):
        """Classifier types should have correct integer values."""
        assert ClassifierType.NACCESS == 0
        assert ClassifierType.PROTOR == 1
        assert ClassifierType.OONS == 2

    def test_classifier_type_is_int(self):
        """Classifier types should be usable as integers."""
        assert isinstance(ClassifierType.NACCESS, int)


class TestAtomClass:
    """Tests for AtomClass enum."""

    def test_atom_class_values(self):
        """Atom classes should have correct integer values."""
        assert AtomClass.POLAR == 0
        assert AtomClass.APOLAR == 1
        assert AtomClass.UNKNOWN == 2

    def test_atom_class_is_int(self):
        """Atom classes should be usable as integers."""
        assert isinstance(AtomClass.POLAR, int)


class TestGetRadius:
    """Tests for get_radius function."""

    def test_get_radius_naccess(self):
        """NACCESS classifier should return correct radii."""
        assert get_radius("ALA", "CA") == pytest.approx(1.87, abs=0.01)
        assert get_radius("ALA", "O") == pytest.approx(1.40, abs=0.01)
        assert get_radius("ALA", "N") == pytest.approx(1.65, abs=0.01)
        assert get_radius("ALA", "CB") == pytest.approx(1.87, abs=0.01)

    def test_get_radius_protor(self):
        """ProtoR classifier should return valid radii."""
        radius = get_radius("ALA", "CA", ClassifierType.PROTOR)
        assert radius is not None
        assert 1.0 < radius < 3.0

    def test_get_radius_oons(self):
        """OONS classifier should return valid radii."""
        radius = get_radius("ALA", "CA", ClassifierType.OONS)
        assert radius is not None
        assert 1.0 < radius < 3.0

    def test_get_radius_unknown_atom(self):
        """Unknown atom should return None."""
        assert get_radius("ALA", "XX") is None
        assert get_radius("XXX", "YY") is None

    def test_get_radius_any_fallback(self):
        """Should fall back to ANY residue for common atoms."""
        # Backbone atoms should work for any residue
        assert get_radius("UNK", "CA") == pytest.approx(1.87, abs=0.01)
        assert get_radius("UNK", "O") == pytest.approx(1.40, abs=0.01)


class TestGetAtomClass:
    """Tests for get_atom_class function."""

    def test_carbon_atoms_apolar(self):
        """Carbon atoms should be classified as apolar."""
        assert get_atom_class("ALA", "CA") == AtomClass.APOLAR
        assert get_atom_class("ALA", "CB") == AtomClass.APOLAR
        assert get_atom_class("ALA", "C") == AtomClass.APOLAR

    def test_nitrogen_oxygen_polar(self):
        """Nitrogen and oxygen should be classified as polar."""
        assert get_atom_class("ALA", "N") == AtomClass.POLAR
        assert get_atom_class("ALA", "O") == AtomClass.POLAR

    def test_unknown_atom_class(self):
        """Unknown atom should return UNKNOWN class."""
        assert get_atom_class("ALA", "XX") == AtomClass.UNKNOWN

    def test_different_classifiers(self):
        """Different classifiers should give consistent polarity."""
        for classifier in [ClassifierType.NACCESS, ClassifierType.PROTOR, ClassifierType.OONS]:
            assert get_atom_class("ALA", "CA", classifier) == AtomClass.APOLAR
            assert get_atom_class("ALA", "O", classifier) == AtomClass.POLAR


class TestGuessRadius:
    """Tests for guess_radius function."""

    def test_common_elements(self):
        """Common biological elements should return correct radii."""
        assert guess_radius("C") == pytest.approx(1.70, abs=0.01)
        assert guess_radius("N") == pytest.approx(1.55, abs=0.01)
        assert guess_radius("O") == pytest.approx(1.52, abs=0.01)
        assert guess_radius("S") == pytest.approx(1.80, abs=0.01)
        assert guess_radius("H") == pytest.approx(1.10, abs=0.01)

    def test_metals(self):
        """Metal elements should return correct radii."""
        assert guess_radius("FE") == pytest.approx(1.26, abs=0.01)
        assert guess_radius("ZN") == pytest.approx(1.39, abs=0.01)
        assert guess_radius("CA") == pytest.approx(2.31, abs=0.01)
        assert guess_radius("MG") == pytest.approx(1.73, abs=0.01)

    def test_case_insensitive(self):
        """Element lookup should be case-insensitive."""
        assert guess_radius("c") == guess_radius("C")
        assert guess_radius("fe") == guess_radius("FE")
        assert guess_radius("Fe") == guess_radius("FE")

    def test_unknown_element(self):
        """Unknown element should return None."""
        assert guess_radius("XX") is None
        assert guess_radius("") is None


class TestGuessRadiusFromAtomName:
    """Tests for guess_radius_from_atom_name function."""

    def test_standard_pdb_atoms(self):
        """Standard PDB atom names should be parsed correctly."""
        # Leading space = single-char element
        assert guess_radius_from_atom_name(" CA ") == pytest.approx(1.70, abs=0.01)  # Carbon
        assert guess_radius_from_atom_name(" N  ") == pytest.approx(1.55, abs=0.01)
        assert guess_radius_from_atom_name(" O  ") == pytest.approx(1.52, abs=0.01)

    def test_metal_atoms(self):
        """Metal atom names without leading space should be 2-char elements."""
        assert guess_radius_from_atom_name("FE  ") == pytest.approx(1.26, abs=0.01)
        assert guess_radius_from_atom_name("ZN  ") == pytest.approx(1.39, abs=0.01)
        assert guess_radius_from_atom_name("CA  ") == pytest.approx(2.31, abs=0.01)  # Calcium

    def test_ca_disambiguation(self):
        """CA with leading space is Carbon, without is Calcium."""
        ca_carbon = guess_radius_from_atom_name(" CA ")  # Carbon alpha
        ca_calcium = guess_radius_from_atom_name("CA  ")  # Calcium ion
        assert ca_carbon == pytest.approx(1.70, abs=0.01)
        assert ca_calcium == pytest.approx(2.31, abs=0.01)
        assert ca_carbon != ca_calcium


class TestClassifyAtoms:
    """Tests for classify_atoms batch function."""

    def test_basic_classification(self):
        """Basic batch classification should work."""
        result = classify_atoms(["ALA", "ALA", "GLY"], ["CA", "O", "N"])

        assert isinstance(result, ClassificationResult)
        assert len(result.radii) == 3
        assert len(result.classes) == 3

        assert result.radii[0] == pytest.approx(1.87, abs=0.01)
        assert result.radii[1] == pytest.approx(1.40, abs=0.01)
        assert result.radii[2] == pytest.approx(1.65, abs=0.01)

        assert result.classes[0] == AtomClass.APOLAR
        assert result.classes[1] == AtomClass.POLAR
        assert result.classes[2] == AtomClass.POLAR

    def test_empty_input(self):
        """Empty input should return empty result."""
        result = classify_atoms([], [])
        assert len(result.radii) == 0
        assert len(result.classes) == 0

    def test_mismatched_lengths(self):
        """Mismatched lengths should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            classify_atoms(["ALA"], ["CA", "O"])

    def test_unknown_atoms_nan(self):
        """Unknown atoms should have NaN radii."""
        result = classify_atoms(["ALA", "XXX"], ["CA", "YY"])

        assert result.radii[0] == pytest.approx(1.87, abs=0.01)
        assert np.isnan(result.radii[1])
        assert result.classes[1] == AtomClass.UNKNOWN

    def test_without_classes(self):
        """Should work without computing classes."""
        result = classify_atoms(["ALA", "ALA"], ["CA", "O"], include_classes=False)

        assert result.radii[0] == pytest.approx(1.87, abs=0.01)
        assert result.radii[1] == pytest.approx(1.40, abs=0.01)
        # Classes should all be UNKNOWN when not computed
        assert all(c == AtomClass.UNKNOWN for c in result.classes)

    def test_different_classifiers(self):
        """Different classifiers should give different radii."""
        result_naccess = classify_atoms(["ALA"], ["CA"], ClassifierType.NACCESS)
        result_protor = classify_atoms(["ALA"], ["CA"], ClassifierType.PROTOR)
        result_oons = classify_atoms(["ALA"], ["CA"], ClassifierType.OONS)

        # All should return valid radii
        assert not np.isnan(result_naccess.radii[0])
        assert not np.isnan(result_protor.radii[0])
        assert not np.isnan(result_oons.radii[0])

    def test_classification_result_repr(self):
        """ClassificationResult should have a clean repr."""
        result = classify_atoms(["ALA", "ALA", "GLY"], ["CA", "O", "N"])
        repr_str = repr(result)
        assert "n_atoms=3" in repr_str


# =============================================================================
# RSA Tests
# =============================================================================


class TestGetMaxSasa:
    """Tests for get_max_sasa function."""

    def test_standard_amino_acids(self):
        """Should return correct values for standard amino acids."""
        from freesasa_zig import MAX_SASA, get_max_sasa

        assert get_max_sasa("ALA") == pytest.approx(129.0, abs=0.01)
        assert get_max_sasa("GLY") == pytest.approx(104.0, abs=0.01)
        assert get_max_sasa("TRP") == pytest.approx(285.0, abs=0.01)
        assert get_max_sasa("ARG") == pytest.approx(274.0, abs=0.01)

        # Check consistency with MAX_SASA dict
        for residue, expected in MAX_SASA.items():
            assert get_max_sasa(residue) == pytest.approx(expected, abs=0.01)

    def test_unknown_residue(self):
        """Unknown residue should return None."""
        from freesasa_zig import get_max_sasa

        assert get_max_sasa("XXX") is None
        assert get_max_sasa("HOH") is None
        assert get_max_sasa("") is None


class TestCalculateRsa:
    """Tests for calculate_rsa function."""

    def test_basic_rsa(self):
        """Should calculate RSA correctly."""
        from freesasa_zig import calculate_rsa

        # ALA: RSA = 64.5 / 129.0 = 0.5
        assert calculate_rsa(64.5, "ALA") == pytest.approx(0.5, abs=0.001)

        # GLY: RSA = 52.0 / 104.0 = 0.5
        assert calculate_rsa(52.0, "GLY") == pytest.approx(0.5, abs=0.001)

        # TRP: RSA = 142.5 / 285.0 = 0.5
        assert calculate_rsa(142.5, "TRP") == pytest.approx(0.5, abs=0.001)

    def test_rsa_greater_than_one(self):
        """RSA > 1.0 is possible for exposed terminal residues."""
        from freesasa_zig import calculate_rsa

        rsa = calculate_rsa(150.0, "GLY")  # 150 / 104 ≈ 1.44
        assert rsa is not None
        assert rsa > 1.0
        assert rsa == pytest.approx(150.0 / 104.0, abs=0.001)

    def test_rsa_zero(self):
        """Zero SASA should give zero RSA."""
        from freesasa_zig import calculate_rsa

        assert calculate_rsa(0.0, "ALA") == pytest.approx(0.0, abs=0.001)

    def test_unknown_residue(self):
        """Unknown residue should return None."""
        from freesasa_zig import calculate_rsa

        assert calculate_rsa(100.0, "XXX") is None
        assert calculate_rsa(100.0, "HOH") is None


class TestCalculateRsaBatch:
    """Tests for calculate_rsa_batch function."""

    def test_basic_batch(self):
        """Should calculate RSA for multiple residues."""
        from freesasa_zig import calculate_rsa_batch

        sasas = np.array([64.5, 52.0, 142.5])
        residues = ["ALA", "GLY", "TRP"]
        rsa = calculate_rsa_batch(sasas, residues)

        assert len(rsa) == 3
        assert rsa[0] == pytest.approx(0.5, abs=0.001)  # ALA: 64.5/129
        assert rsa[1] == pytest.approx(0.5, abs=0.001)  # GLY: 52/104
        assert rsa[2] == pytest.approx(0.5, abs=0.001)  # TRP: 142.5/285

    def test_batch_with_unknown(self):
        """Unknown residues should have NaN RSA."""
        from freesasa_zig import calculate_rsa_batch

        sasas = np.array([64.5, 100.0, 52.0])
        residues = ["ALA", "HOH", "GLY"]
        rsa = calculate_rsa_batch(sasas, residues)

        assert rsa[0] == pytest.approx(0.5, abs=0.001)  # ALA
        assert np.isnan(rsa[1])  # HOH - unknown
        assert rsa[2] == pytest.approx(0.5, abs=0.001)  # GLY

    def test_batch_empty(self):
        """Empty input should return empty array."""
        from freesasa_zig import calculate_rsa_batch

        rsa = calculate_rsa_batch(np.array([]), [])
        assert len(rsa) == 0

    def test_batch_from_list(self):
        """Should accept list as input."""
        from freesasa_zig import calculate_rsa_batch

        sasas = [64.5, 52.0]
        residues = ["ALA", "GLY"]
        rsa = calculate_rsa_batch(sasas, residues)

        assert len(rsa) == 2
        assert rsa[0] == pytest.approx(0.5, abs=0.001)

    def test_batch_mismatched_lengths(self):
        """Mismatched lengths should raise ValueError."""
        from freesasa_zig import calculate_rsa_batch

        with pytest.raises(ValueError, match="same length"):
            calculate_rsa_batch(np.array([64.5]), ["ALA", "GLY"])

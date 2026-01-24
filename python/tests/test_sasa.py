"""Tests for freesasa-zig Python bindings."""

import numpy as np
import pytest

from freesasa_zig import SasaResult, calculate_sasa, get_version


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

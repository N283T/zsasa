"""Tests for batch SASA calculation."""

from __future__ import annotations

import numpy as np
import pytest

from zsasa import BatchSasaResult, calculate_sasa_batch


class TestBatchSasaResult:
    """Tests for BatchSasaResult dataclass."""

    def test_properties(self) -> None:
        """Test n_frames, n_atoms, total_areas properties."""
        atom_areas = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32)
        result = BatchSasaResult(atom_areas=atom_areas)

        assert result.n_frames == 2
        assert result.n_atoms == 3
        np.testing.assert_array_almost_equal(result.total_areas, [6.0, 15.0])

    def test_repr(self) -> None:
        """Test string representation."""
        atom_areas = np.array([[1.0, 2.0]], dtype=np.float32)
        result = BatchSasaResult(atom_areas=atom_areas)

        repr_str = repr(result)
        assert "BatchSasaResult" in repr_str
        assert "n_frames=1" in repr_str
        assert "n_atoms=2" in repr_str


class TestCalculateSasaBatch:
    """Tests for calculate_sasa_batch function."""

    def test_single_frame(self) -> None:
        """Test batch with a single frame."""
        coords = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        result = calculate_sasa_batch(coords, radii)

        assert result.n_frames == 1
        assert result.n_atoms == 1
        assert result.atom_areas.shape == (1, 1)
        # Single isolated atom: area = 4 * pi * (r + probe)^2
        expected_area = 4 * np.pi * (1.5 + 1.4) ** 2
        assert result.atom_areas[0, 0] == pytest.approx(expected_area, rel=0.05)

    def test_multiple_frames(self) -> None:
        """Test batch with multiple frames."""
        n_frames = 10
        n_atoms = 50
        # Create random coordinates with atoms spread apart
        coords = np.random.randn(n_frames, n_atoms, 3).astype(np.float32) * 5.0
        radii = np.full(n_atoms, 1.5, dtype=np.float32)

        result = calculate_sasa_batch(coords, radii)

        assert result.atom_areas.shape == (n_frames, n_atoms)
        assert result.n_frames == n_frames
        assert result.n_atoms == n_atoms
        # All areas should be positive
        assert np.all(result.atom_areas >= 0)
        # Total areas should be reasonable
        assert np.all(result.total_areas > 0)

    def test_overlapping_atoms(self) -> None:
        """Test that overlapping atoms reduce SASA."""
        # Two atoms at same position
        coords_overlap = np.array([[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]], dtype=np.float32)
        # Two atoms far apart
        coords_far = np.array([[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]], dtype=np.float32)
        radii = np.array([1.5, 1.5], dtype=np.float32)

        result_overlap = calculate_sasa_batch(coords_overlap, radii)
        result_far = calculate_sasa_batch(coords_far, radii)

        # Overlapping atoms should have less total SASA
        assert result_overlap.total_areas[0] < result_far.total_areas[0]

    def test_threading_consistency(self) -> None:
        """Results should match regardless of thread count."""
        n_frames = 5
        n_atoms = 20
        coords = np.random.randn(n_frames, n_atoms, 3).astype(np.float32) * 3.0
        radii = np.full(n_atoms, 1.5, dtype=np.float32)

        result_1thread = calculate_sasa_batch(coords, radii, n_threads=1)
        result_4threads = calculate_sasa_batch(coords, radii, n_threads=4)

        np.testing.assert_array_almost_equal(
            result_1thread.atom_areas, result_4threads.atom_areas, decimal=4
        )

    def test_lr_algorithm(self) -> None:
        """Test Lee-Richards algorithm in batch mode."""
        coords = np.array(
            [[[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]], [[1.0, 0.0, 0.0], [4.0, 0.0, 0.0]]],
            dtype=np.float32,
        )
        radii = np.array([1.5, 1.5], dtype=np.float32)

        result = calculate_sasa_batch(coords, radii, algorithm="lr", n_slices=20)

        assert result.atom_areas.shape == (2, 2)
        assert np.all(result.atom_areas > 0)

    def test_sr_lr_similar_results(self) -> None:
        """SR and LR algorithms should give similar results."""
        # Use fixed seed for reproducibility
        rng = np.random.default_rng(42)
        coords = rng.standard_normal((3, 10, 3)).astype(np.float32) * 3.0
        radii = np.full(10, 1.5, dtype=np.float32)

        result_sr = calculate_sasa_batch(coords, radii, algorithm="sr", n_points=500)
        result_lr = calculate_sasa_batch(coords, radii, algorithm="lr", n_slices=50)

        # Should be within 2% of each other
        np.testing.assert_allclose(result_sr.total_areas, result_lr.total_areas, rtol=0.02)

    def test_probe_radius(self) -> None:
        """Different probe radii should affect results."""
        coords = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        result_small = calculate_sasa_batch(coords, radii, probe_radius=1.0)
        result_large = calculate_sasa_batch(coords, radii, probe_radius=2.0)

        # Larger probe = larger effective radius = larger area
        assert result_large.atom_areas[0, 0] > result_small.atom_areas[0, 0]


class TestCalculateSasaBatchPerformance:
    """Performance sanity checks for batch API."""

    def test_batch_faster_than_sequential(self) -> None:
        """Batch should be faster than sequential single-frame calls."""
        import time

        from zsasa import calculate_sasa

        n_frames = 20
        n_atoms = 100
        rng = np.random.default_rng(42)
        coords = rng.standard_normal((n_frames, n_atoms, 3)).astype(np.float32) * 3.0
        radii = np.full(n_atoms, 1.5, dtype=np.float32)

        # Time batch (with multiple threads)
        start = time.perf_counter()
        calculate_sasa_batch(coords, radii, n_threads=4)
        batch_time = time.perf_counter() - start

        # Time sequential
        start = time.perf_counter()
        for i in range(n_frames):
            calculate_sasa(coords[i], radii)
        seq_time = time.perf_counter() - start

        # Batch should be faster (at least not slower)
        # Note: On small inputs, overhead may make batch similar speed
        assert batch_time < seq_time * 1.5, (
            f"Batch ({batch_time:.3f}s) should not be much slower than sequential ({seq_time:.3f}s)"
        )


class TestCalculateSasaBatchValidation:
    """Tests for input validation in calculate_sasa_batch."""

    def test_invalid_coords_shape_2d(self) -> None:
        """Should reject 2D coordinates."""
        coords = np.array([[0.0, 0.0, 0.0]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        with pytest.raises(ValueError, match=r"n_frames, n_atoms, 3"):
            calculate_sasa_batch(coords, radii)

    def test_invalid_coords_shape_wrong_last_dim(self) -> None:
        """Should reject coordinates without 3 components."""
        coords = np.array([[[0.0, 0.0]]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        with pytest.raises(ValueError, match="3"):
            calculate_sasa_batch(coords, radii)

    def test_mismatched_atom_count(self) -> None:
        """Should reject mismatched atom counts."""
        coords = np.array([[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)  # Only 1 radius for 2 atoms

        with pytest.raises(ValueError, match=r"radii must be"):
            calculate_sasa_batch(coords, radii)

    def test_invalid_algorithm(self) -> None:
        """Should reject invalid algorithm."""
        coords = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
        radii = np.array([1.5], dtype=np.float32)

        with pytest.raises(ValueError, match="algorithm"):
            calculate_sasa_batch(coords, radii, algorithm="invalid")

    def test_empty_coordinates(self) -> None:
        """Should raise error for empty coordinates."""
        coords = np.zeros((0, 0, 3), dtype=np.float32)
        radii = np.array([], dtype=np.float32)

        with pytest.raises(ValueError, match="Invalid input"):
            calculate_sasa_batch(coords, radii)

    def test_negative_radii(self) -> None:
        """Should reject negative radii."""
        coords = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
        radii = np.array([-1.5], dtype=np.float32)

        with pytest.raises(ValueError, match="negative|positive"):
            calculate_sasa_batch(coords, radii)

"""Tests for native XTC reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from zsasa.xtc import XtcReader, TrajectorySasaResult, compute_sasa_trajectory


# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"
XTC_FILE = TEST_DATA_DIR / "1l2y.xtc"


class TestXtcReader:
    """Test XtcReader class."""

    @pytest.fixture
    def reader(self) -> XtcReader:
        """Create XtcReader for test file."""
        return XtcReader(XTC_FILE)

    def test_open_close(self) -> None:
        """Test opening and closing XTC file."""
        reader = XtcReader(XTC_FILE)
        assert reader.natoms == 304
        reader.close()

    def test_context_manager(self) -> None:
        """Test using XtcReader as context manager."""
        with XtcReader(XTC_FILE) as reader:
            assert reader.natoms == 304

    def test_file_not_found(self) -> None:
        """Test opening non-existent file raises error."""
        with pytest.raises(FileNotFoundError):
            XtcReader("nonexistent.xtc")

    def test_natoms(self, reader: XtcReader) -> None:
        """Test natoms property."""
        assert reader.natoms == 304
        reader.close()

    def test_read_first_frame(self, reader: XtcReader) -> None:
        """Test reading first frame."""
        frame = reader.read_frame()
        reader.close()

        assert frame is not None
        assert frame.step == 1
        assert frame.coords.shape == (304, 3)
        assert frame.box.shape == (3, 3)

        # Check first atom coordinates (in nm, from C xdrfile reference)
        np.testing.assert_allclose(frame.coords[0], [-0.8901, 0.4127, -0.0555], atol=0.001)

    def test_read_all_frames(self, reader: XtcReader) -> None:
        """Test reading all frames via iteration."""
        frame_count = 0
        for frame in reader:
            frame_count += 1
            assert frame.coords.shape == (304, 3)

        # 1l2y.xtc has 38 frames
        assert frame_count == 38

    def test_read_frame_returns_none_at_eof(self, reader: XtcReader) -> None:
        """Test that read_frame returns None at end of file."""
        # Skip all frames
        for _ in reader:
            pass

        # Next read should return None
        frame = reader.read_frame()
        reader.close()
        assert frame is None

    def test_frame_properties(self) -> None:
        """Test XtcFrame properties."""
        with XtcReader(XTC_FILE) as reader:
            frame = reader.read_frame()
            assert frame is not None

            assert frame.natoms == 304
            assert frame.step == 1
            assert frame.time >= 0
            assert frame.precision > 0


class TestComputeSasaTrajectory:
    """Test compute_sasa_trajectory function."""

    @pytest.fixture
    def radii(self) -> np.ndarray:
        """Create radii array for 1l2y (304 atoms)."""
        # Use uniform radii for simplicity (1.7 Å = typical carbon)
        return np.full(304, 1.7, dtype=np.float32)

    def test_basic(self, radii: np.ndarray) -> None:
        """Test basic SASA calculation."""
        result = compute_sasa_trajectory(XTC_FILE, radii)

        assert isinstance(result, TrajectorySasaResult)
        assert result.n_frames == 38
        assert result.n_atoms == 304
        assert result.atom_areas.shape == (38, 304)

    def test_total_areas(self, radii: np.ndarray) -> None:
        """Test total_areas property."""
        result = compute_sasa_trajectory(XTC_FILE, radii)

        assert result.total_areas.shape == (38,)
        # Total SASA should be positive
        assert np.all(result.total_areas > 0)

    def test_steps_and_times(self, radii: np.ndarray) -> None:
        """Test steps and times arrays."""
        result = compute_sasa_trajectory(XTC_FILE, radii)

        assert len(result.steps) == 38
        assert len(result.times) == 38

    def test_frame_range(self, radii: np.ndarray) -> None:
        """Test processing subset of frames."""
        # Process only frames 5-15
        result = compute_sasa_trajectory(XTC_FILE, radii, start=5, stop=15)

        assert result.n_frames == 10

    def test_frame_step(self, radii: np.ndarray) -> None:
        """Test processing every Nth frame."""
        # Process every 5th frame
        result = compute_sasa_trajectory(XTC_FILE, radii, step=5)

        # 38 frames, every 5th = frames 0, 5, 10, 15, 20, 25, 30, 35 = 8 frames
        assert result.n_frames == 8

    def test_algorithms(self, radii: np.ndarray) -> None:
        """Test different algorithms give similar results."""
        result_sr = compute_sasa_trajectory(XTC_FILE, radii, algorithm="sr", n_points=500)
        result_lr = compute_sasa_trajectory(XTC_FILE, radii, algorithm="lr", n_slices=50)

        # Results should be similar (within 5%)
        np.testing.assert_allclose(
            result_sr.total_areas,
            result_lr.total_areas,
            rtol=0.05,
        )

    def test_radii_length_mismatch(self) -> None:
        """Test that mismatched radii length raises error."""
        wrong_radii = np.full(100, 1.7, dtype=np.float32)  # Wrong length

        with pytest.raises(ValueError, match="radii length"):
            compute_sasa_trajectory(XTC_FILE, wrong_radii)

    def test_threading(self, radii: np.ndarray) -> None:
        """Test that threading doesn't affect results."""
        result_1 = compute_sasa_trajectory(XTC_FILE, radii, n_threads=1)
        result_4 = compute_sasa_trajectory(XTC_FILE, radii, n_threads=4)

        np.testing.assert_array_almost_equal(
            result_1.atom_areas,
            result_4.atom_areas,
            decimal=4,
        )


class TestXtcReaderClosed:
    """Test behavior with closed reader."""

    def test_read_after_close_raises(self) -> None:
        """Test that reading after close raises error."""
        reader = XtcReader(XTC_FILE)
        reader.close()

        with pytest.raises(RuntimeError, match="closed"):
            reader.read_frame()

    def test_double_close_safe(self) -> None:
        """Test that closing twice is safe."""
        reader = XtcReader(XTC_FILE)
        reader.close()
        reader.close()  # Should not raise


class TestXtcCoordinateUnits:
    """Test coordinate unit handling."""

    def test_coords_in_nanometers(self) -> None:
        """Test that coordinates are in nanometers."""
        with XtcReader(XTC_FILE) as reader:
            frame = reader.read_frame()
            assert frame is not None

            # 1l2y is a small protein, coordinates should be in nm range
            # (roughly -1 to +1.5 nm for this protein)
            assert np.all(frame.coords > -5)  # No coordinate > 5 nm away
            assert np.all(frame.coords < 5)

    def test_sasa_result_in_angstrom_squared(self) -> None:
        """Test that SASA output is in Å²."""
        radii = np.full(304, 1.7, dtype=np.float32)
        result = compute_sasa_trajectory(XTC_FILE, radii)

        # Total SASA for a small protein should be ~5000-10000 Å²
        # If it were nm², it would be ~50-100
        assert np.all(result.total_areas > 100)  # Must be in Å²
        assert np.all(result.total_areas < 50000)  # Reasonable upper bound


class TestRepr:
    """Test string representations."""

    def test_trajectory_sasa_result_repr(self) -> None:
        """Test TrajectorySasaResult repr."""
        radii = np.full(304, 1.7, dtype=np.float32)
        result = compute_sasa_trajectory(XTC_FILE, radii, stop=5)

        repr_str = repr(result)
        assert "TrajectorySasaResult" in repr_str
        assert "n_frames=5" in repr_str
        assert "n_atoms=304" in repr_str

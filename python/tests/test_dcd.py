"""Tests for native DCD reader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from zsasa.dcd import DcdFrame, DcdReader, compute_sasa_trajectory
from zsasa.xtc import TrajectorySasaResult

# Path to test data
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"
DCD_FILE = TEST_DATA_DIR / "1l2y.dcd"
XTC_FILE = TEST_DATA_DIR / "1l2y.xtc"


def _xtc_enabled() -> bool:
    """Return True if XTC support was compiled in."""
    from zsasa.xtc import XtcReader

    try:
        XtcReader("/nonexistent/path.xtc")
    except RuntimeError as e:
        if "Error opening XTC file: -6" in str(e):
            return False
        return True
    except (FileNotFoundError, OSError):
        return True
    return True


_skip_if_xtc_disabled = pytest.mark.skipif(
    not _xtc_enabled(),
    reason="XTC support disabled (rebuild with -Dxtc=true)",
)


@pytest.fixture(autouse=True)
def _skip_if_no_dcd() -> None:
    """Skip tests if DCD test file is not available."""
    if not DCD_FILE.exists():
        pytest.skip(f"DCD test file not found: {DCD_FILE}")


class TestDcdReader:
    """Test DcdReader class."""

    @pytest.fixture
    def reader(self) -> DcdReader:
        """Create DcdReader for test file."""
        return DcdReader(DCD_FILE)

    def test_open_close(self) -> None:
        """Test opening and closing DCD file."""
        reader = DcdReader(DCD_FILE)
        assert reader.natoms == 304
        reader.close()

    def test_context_manager(self) -> None:
        """Test using DcdReader as context manager."""
        with DcdReader(DCD_FILE) as reader:
            assert reader.natoms == 304

    def test_file_not_found(self) -> None:
        """Test opening non-existent file raises error."""
        with pytest.raises(FileNotFoundError):
            DcdReader("nonexistent.dcd")

    def test_natoms(self, reader: DcdReader) -> None:
        """Test natoms property."""
        assert reader.natoms == 304
        reader.close()

    def test_read_first_frame(self, reader: DcdReader) -> None:
        """Test reading first frame."""
        frame = reader.read_frame()
        reader.close()

        assert frame is not None
        assert frame.coords.shape == (304, 3)

        # DCD coordinates are in Angstroms
        # XTC reference: atom[0] = [-0.8901, 0.4127, -0.0555] nm
        # In Angstroms: [-8.901, 4.127, -0.555]
        np.testing.assert_allclose(frame.coords[0], [-8.901, 4.127, -0.555], atol=0.05)

    def test_read_all_frames(self, reader: DcdReader) -> None:
        """Test reading all frames via iteration."""
        frame_count = 0
        for frame in reader:
            frame_count += 1
            assert frame.coords.shape == (304, 3)

        # 1l2y.dcd has 38 frames (same as XTC)
        assert frame_count == 38

    def test_read_frame_returns_none_at_eof(self, reader: DcdReader) -> None:
        """Test that read_frame returns None at end of file."""
        for _ in reader:
            pass

        frame = reader.read_frame()
        reader.close()
        assert frame is None

    def test_frame_properties(self) -> None:
        """Test DcdFrame properties."""
        with DcdReader(DCD_FILE) as reader:
            frame = reader.read_frame()
            assert frame is not None

            assert frame.natoms == 304
            assert isinstance(frame.step, int)
            assert isinstance(frame.time, float)


class TestDcdCoordinateUnits:
    """Test coordinate unit handling."""

    def test_coords_in_angstroms(self) -> None:
        """Test that coordinates are in Angstroms."""
        with DcdReader(DCD_FILE) as reader:
            frame = reader.read_frame()
            assert frame is not None

            # 1l2y is a small protein, coordinates in Angstroms
            # should be roughly -10 to +15 A range
            assert np.all(frame.coords > -50)
            assert np.all(frame.coords < 50)

    @_skip_if_xtc_disabled
    def test_dcd_matches_xtc_coordinates(self) -> None:
        """Test that DCD coordinates match XTC coordinates (after unit conversion)."""
        from zsasa.xtc import XtcReader

        if not XTC_FILE.exists():
            pytest.skip("XTC test file not found")

        with DcdReader(DCD_FILE) as dcd_reader, XtcReader(XTC_FILE) as xtc_reader:
            dcd_frame = dcd_reader.read_frame()
            xtc_frame = xtc_reader.read_frame()

            assert dcd_frame is not None
            assert xtc_frame is not None

            # XTC coords in nm, DCD in Angstroms
            xtc_angstrom = xtc_frame.coords * 10.0

            np.testing.assert_allclose(dcd_frame.coords, xtc_angstrom, atol=0.05)


class TestComputeSasaTrajectory:
    """Test compute_sasa_trajectory function."""

    @pytest.fixture
    def radii(self) -> np.ndarray:
        """Create radii array for 1l2y (304 atoms)."""
        return np.full(304, 1.7, dtype=np.float32)

    def test_basic(self, radii: np.ndarray) -> None:
        """Test basic SASA calculation."""
        result = compute_sasa_trajectory(DCD_FILE, radii)

        assert isinstance(result, TrajectorySasaResult)
        assert result.n_frames == 38
        assert result.n_atoms == 304
        assert result.atom_areas.shape == (38, 304)

    def test_total_areas(self, radii: np.ndarray) -> None:
        """Test total_areas property."""
        result = compute_sasa_trajectory(DCD_FILE, radii)

        assert result.total_areas.shape == (38,)
        assert np.all(result.total_areas > 0)

    @_skip_if_xtc_disabled
    def test_sasa_matches_xtc(self, radii: np.ndarray) -> None:
        """Test that DCD SASA matches XTC SASA."""
        from zsasa.xtc import compute_sasa_trajectory as compute_xtc

        if not XTC_FILE.exists():
            pytest.skip("XTC test file not found")

        dcd_result = compute_sasa_trajectory(DCD_FILE, radii)
        xtc_result = compute_xtc(XTC_FILE, radii)

        # Total SASA should match within tolerance
        np.testing.assert_allclose(
            dcd_result.total_areas,
            xtc_result.total_areas,
            rtol=0.01,  # 1% tolerance
        )

    def test_frame_range(self, radii: np.ndarray) -> None:
        """Test processing subset of frames."""
        result = compute_sasa_trajectory(DCD_FILE, radii, start=5, stop=15)
        assert result.n_frames == 10

    def test_frame_step(self, radii: np.ndarray) -> None:
        """Test processing every Nth frame."""
        result = compute_sasa_trajectory(DCD_FILE, radii, step=5)
        # 38 frames, every 5th = 8 frames
        assert result.n_frames == 8

    def test_radii_length_mismatch(self) -> None:
        """Test that mismatched radii length raises error."""
        wrong_radii = np.full(100, 1.7, dtype=np.float32)
        with pytest.raises(ValueError, match="radii length"):
            compute_sasa_trajectory(DCD_FILE, wrong_radii)

    def test_sasa_in_angstrom_squared(self, radii: np.ndarray) -> None:
        """Test that SASA output is in Å²."""
        result = compute_sasa_trajectory(DCD_FILE, radii)

        # Total SASA for a small protein should be ~5000-10000 Å²
        assert np.all(result.total_areas > 100)
        assert np.all(result.total_areas < 50000)


class TestDcdReaderClosed:
    """Test behavior with closed reader."""

    def test_read_after_close_raises(self) -> None:
        """Test that reading after close raises error."""
        reader = DcdReader(DCD_FILE)
        reader.close()

        with pytest.raises(RuntimeError, match="closed"):
            reader.read_frame()

    def test_double_close_safe(self) -> None:
        """Test that closing twice is safe."""
        reader = DcdReader(DCD_FILE)
        reader.close()
        reader.close()  # Should not raise


class TestRepr:
    """Test string representations."""

    def test_dcd_frame_repr(self) -> None:
        """Test DcdFrame has expected fields."""
        with DcdReader(DCD_FILE) as reader:
            frame = reader.read_frame()
            assert frame is not None
            assert isinstance(frame, DcdFrame)
            assert hasattr(frame, "step")
            assert hasattr(frame, "time")
            assert hasattr(frame, "coords")
            assert hasattr(frame, "unitcell")

"""Tests for directory batch processing via process_directory()."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from zsasa import BatchDirResult, ClassifierType, process_directory

# test_data/ lives at the project root
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"


class TestProcessDirectory:
    """Tests for process_directory()."""

    def test_process_directory_with_test_data(self) -> None:
        """Process test_data with SR algorithm and verify result fields."""
        result = process_directory(TEST_DATA_DIR, algorithm="sr")

        assert isinstance(result, BatchDirResult)
        assert result.total_files > 0
        assert result.successful > 0
        assert result.failed == 0
        assert len(result.filenames) == result.total_files
        assert len(result.n_atoms) == result.total_files
        assert len(result.total_sasa) == result.total_files
        assert len(result.status) == result.total_files

    def test_process_directory_lr(self) -> None:
        """Process test_data with LR algorithm."""
        result = process_directory(TEST_DATA_DIR, algorithm="lr")

        assert result.successful > 0
        for i in range(result.total_files):
            if result.status[i] == 1:
                assert result.total_sasa[i] > 0.0
                assert result.n_atoms[i] > 0

    def test_process_directory_nonexistent(self) -> None:
        """FileNotFoundError for a non-existent path."""
        with pytest.raises(FileNotFoundError):
            process_directory("/nonexistent/path/that/does/not/exist")

    def test_process_directory_classifier_none(self) -> None:
        """classifier=None uses input radii (classifier_type=-1)."""
        result = process_directory(TEST_DATA_DIR, classifier=None)

        assert result.total_files > 0
        assert result.successful > 0

    def test_process_directory_result_properties(self) -> None:
        """Verify all BatchDirResult dataclass fields are populated correctly."""
        result = process_directory(TEST_DATA_DIR)

        assert result.total_files == result.successful + result.failed

        for i in range(result.total_files):
            assert isinstance(result.filenames[i], str)
            assert len(result.filenames[i]) > 0
            assert isinstance(result.n_atoms[i], int)
            assert isinstance(result.status[i], int)
            assert result.status[i] in (0, 1)

            if result.status[i] == 1:
                assert result.n_atoms[i] > 0
                assert result.total_sasa[i] > 0.0
                assert not math.isnan(result.total_sasa[i])

    def test_process_directory_output_dir(self, tmp_path: Path) -> None:
        """Process with an output directory writes per-file results."""
        result = process_directory(TEST_DATA_DIR, output_dir=tmp_path)

        assert result.successful > 0
        output_files = list(tmp_path.iterdir())
        assert len(output_files) > 0

    def test_process_directory_invalid_algorithm(self) -> None:
        """ValueError for an invalid algorithm string."""
        with pytest.raises(ValueError, match="Unknown algorithm"):
            process_directory(TEST_DATA_DIR, algorithm="invalid")

    def test_process_directory_invalid_probe_radius(self) -> None:
        """ValueError for non-positive probe_radius."""
        with pytest.raises(ValueError, match="probe_radius must be positive"):
            process_directory(TEST_DATA_DIR, probe_radius=-1.0)

    def test_process_directory_zero_probe_radius(self) -> None:
        """probe_radius=0 should be rejected (boundary of <= 0 check)."""
        with pytest.raises(ValueError, match="probe_radius must be positive"):
            process_directory(TEST_DATA_DIR, probe_radius=0.0)

    def test_process_directory_string_path(self) -> None:
        """Accept string paths (not just Path objects)."""
        result = process_directory(str(TEST_DATA_DIR))

        assert result.total_files > 0
        assert result.successful > 0

    def test_process_directory_naccess_classifier(self) -> None:
        """Process with NACCESS classifier."""
        result = process_directory(
            TEST_DATA_DIR, classifier=ClassifierType.NACCESS
        )

        assert result.successful > 0

    def test_process_directory_empty_dir(self, tmp_path: Path) -> None:
        """Empty directory returns result with zero files."""
        result = process_directory(tmp_path)

        assert result.total_files == 0
        assert result.successful == 0
        assert result.failed == 0
        assert result.filenames == []

    def test_process_directory_include_hydrogens(self) -> None:
        """include_hydrogens=True should be accepted without error."""
        result = process_directory(TEST_DATA_DIR, include_hydrogens=True)

        assert result.successful > 0

    def test_process_directory_include_hetatm(self) -> None:
        """include_hetatm=True should be accepted without error."""
        result = process_directory(TEST_DATA_DIR, include_hetatm=True)

        assert result.successful > 0

    def test_process_directory_explicit_threads(self) -> None:
        """Explicit thread count should produce valid results."""
        result = process_directory(TEST_DATA_DIR, n_threads=1)

        assert result.successful > 0

    def test_process_directory_repr(self) -> None:
        """BatchDirResult repr shows summary counts."""
        result = process_directory(TEST_DATA_DIR)

        repr_str = repr(result)
        assert "BatchDirResult" in repr_str
        assert f"total_files={result.total_files}" in repr_str
        assert f"successful={result.successful}" in repr_str

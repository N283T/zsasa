"""Tests for directory batch processing via process_directory()."""

from __future__ import annotations

import math
import shutil
from pathlib import Path

import pytest

from zsasa import BatchDirResult, ClassifierType, process_directory

# test_data/ lives at the project root
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"

# Tiny structures for tests that must stay fast (e.g. Lee-Richards).
SMALL_ALA_PDB = (
    "ATOM      1  N   ALA A   1      -0.966   0.493   1.500  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       0.257   0.418   0.692  1.00  0.00           C\n"
    "ATOM      3  C   ALA A   1      -0.094   0.017  -0.716  1.00  0.00           C\n"
    "ATOM      4  O   ALA A   1      -1.056  -0.682  -0.923  1.00  0.00           O\n"
    "ATOM      5  CB  ALA A   1       1.204  -0.620   1.296  1.00  0.00           C\n"
    "END\n"
)
SMALL_GLY_PDB = (
    "ATOM      1  N   GLY A   1      10.000  10.000  10.000  1.00  0.00           N\n"
    "ATOM      2  CA  GLY A   1      11.450  10.000  10.000  1.00  0.00           C\n"
    "ATOM      3  C   GLY A   1      11.980  11.420  10.000  1.00  0.00           C\n"
    "ATOM      4  O   GLY A   1      11.230  12.390  10.000  1.00  0.00           O\n"
    "END\n"
)


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

    def test_process_directory_lr(self, tmp_path: Path) -> None:
        """LR on a small directory matches `zsasa calc --algorithm=lr` and differs from SR.

        test_data/1l2y.pdb (38 NMR models superimposed into 11,552 atoms) takes
        seconds with Lee-Richards, so this uses two tiny structures instead.
        """
        (tmp_path / "ala.pdb").write_text(SMALL_ALA_PDB)
        (tmp_path / "gly.ent").write_text(SMALL_GLY_PDB)

        result = process_directory(tmp_path, algorithm="lr", classifier=ClassifierType.NACCESS)

        assert result.total_files == 2
        assert result.successful == 2
        assert result.failed == 0
        n_atoms = dict(zip(result.filenames, result.n_atoms, strict=True))
        lr_area = dict(zip(result.filenames, result.total_sasa, strict=True))
        assert n_atoms == {"ala.pdb": 5, "gly.ent": 4}
        # Reference values from `zsasa calc --algorithm=lr --classifier=naccess`
        # (exact arc angles, the default; --lr-trig=fast gives 215.14412396821933
        # for ala.pdb, the value of zsasa 0.9.1).
        assert lr_area["ala.pdb"] == pytest.approx(215.07786201196504, abs=1e-6)
        assert lr_area["gly.ent"] == pytest.approx(188.2086929717762, abs=1e-6)

        # Shrake-Rupley gives a different, close area, so the algorithm is honored.
        sr = process_directory(tmp_path, algorithm="sr", classifier=ClassifierType.NACCESS)
        sr_area = dict(zip(sr.filenames, sr.total_sasa, strict=True))
        assert sr_area["ala.pdb"] == pytest.approx(215.8249020274959, abs=1e-6)
        assert sr_area["gly.ent"] == pytest.approx(187.0015182803052, abs=1e-6)
        for name, area in lr_area.items():
            assert abs(area - sr_area[name]) > 0.1

    def test_process_directory_nonexistent(self) -> None:
        """FileNotFoundError for a non-existent path."""
        with pytest.raises(FileNotFoundError):
            process_directory("/nonexistent/path/that/does/not/exist")

    def test_process_directory_input_is_a_file(self, tmp_path: Path) -> None:
        """A file given as the input directory is not a missing directory."""
        not_a_dir = tmp_path / "ala.pdb"
        not_a_dir.write_text(SMALL_ALA_PDB)

        with pytest.raises(NotADirectoryError, match="Input path is not a directory"):
            process_directory(not_a_dir)

    def test_process_directory_unreadable_output_dir_blames_the_output(
        self, tmp_path: Path
    ) -> None:
        """An uncreatable output directory is not reported as a bad input directory.

        The input directory exists and is fine; the output path lies below a
        regular file, so it can never be created.
        """
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        (input_dir / "ala.pdb").write_text(SMALL_ALA_PDB)
        blocker = tmp_path / "blocker"
        blocker.write_text("a regular file")
        output_dir = blocker / "sub"

        for n_threads in (1, 4):
            with pytest.raises(
                NotADirectoryError, match="Cannot create output directory"
            ) as excinfo:
                process_directory(input_dir, output_dir=output_dir, n_threads=n_threads)
            message = str(excinfo.value)
            assert str(output_dir) in message
            assert str(blocker) in message
            assert "Input" not in message

    def test_process_directory_output_dir_is_a_file(self, tmp_path: Path) -> None:
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        (input_dir / "ala.pdb").write_text(SMALL_ALA_PDB)
        output_path = tmp_path / "out.txt"
        output_path.write_text("a regular file")

        with pytest.raises(FileExistsError, match="not a directory"):
            process_directory(input_dir, output_dir=output_path)

    def test_process_directory_missing_input_with_valid_output(self, tmp_path: Path) -> None:
        """The missing directory is the input one, so the input error is raised."""
        with pytest.raises(FileNotFoundError, match="Input directory not found"):
            process_directory(tmp_path / "missing", output_dir=tmp_path / "out")
        assert not (tmp_path / "out").exists()

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

    def test_process_directory_output_name_collision(self, tmp_path: Path) -> None:
        """Inputs sharing a stem are rejected instead of overwriting each other."""
        input_dir = tmp_path / "in"
        output_dir = tmp_path / "out"
        input_dir.mkdir()
        shutil.copy(TEST_DATA_DIR / "1l2y.pdb", input_dir / "1l2y.pdb")
        shutil.copy(TEST_DATA_DIR / "1l2y.pdb", input_dir / "1l2y.ent")

        with pytest.raises(ValueError, match="same output file name"):
            process_directory(input_dir, output_dir=output_dir)
        assert not output_dir.exists()

        result = process_directory(input_dir)
        assert result.total_files == 2

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
        result = process_directory(TEST_DATA_DIR, classifier=ClassifierType.NACCESS)

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

    def test_process_directory_parallel_matches_single_thread(self) -> None:
        """Parallel directory processing should match single-thread results."""
        result_1 = process_directory(TEST_DATA_DIR, n_threads=1)
        result_4 = process_directory(TEST_DATA_DIR, n_threads=4)

        assert result_4.total_files == result_1.total_files
        assert result_4.successful == result_1.successful
        assert result_4.failed == result_1.failed

        by_name_1 = {
            filename: (n_atoms, total_sasa, status)
            for filename, n_atoms, total_sasa, status in zip(
                result_1.filenames,
                result_1.n_atoms,
                result_1.total_sasa,
                result_1.status,
                strict=True,
            )
        }
        by_name_4 = {
            filename: (n_atoms, total_sasa, status)
            for filename, n_atoms, total_sasa, status in zip(
                result_4.filenames,
                result_4.n_atoms,
                result_4.total_sasa,
                result_4.status,
                strict=True,
            )
        }
        assert by_name_4.keys() == by_name_1.keys()
        for filename, (n_atoms_1, total_sasa_1, status_1) in by_name_1.items():
            n_atoms_4, total_sasa_4, status_4 = by_name_4[filename]
            assert n_atoms_4 == n_atoms_1
            assert status_4 == status_1
            if status_1 == 1:
                assert total_sasa_4 == pytest.approx(total_sasa_1)

    def test_process_directory_repr(self) -> None:
        """BatchDirResult repr shows summary counts."""
        result = process_directory(TEST_DATA_DIR)

        repr_str = repr(result)
        assert "BatchDirResult" in repr_str
        assert f"total_files={result.total_files}" in repr_str
        assert f"successful={result.successful}" in repr_str

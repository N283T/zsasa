"""The exceptions of the native trajectory readers match the cause of the failure."""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from zsasa import dcd, xtc

TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"
XTC_FILE = TEST_DATA_DIR / "1l2y.xtc"
DCD_FILE = TEST_DATA_DIR / "1l2y.dcd"
PDB_FILE = TEST_DATA_DIR / "1l2y.pdb"

READERS = {
    "xtc": (xtc.XtcReader, XTC_FILE, DCD_FILE),
    "dcd": (dcd.DcdReader, DCD_FILE, XTC_FILE),
}
COMPUTE = {
    "xtc": (xtc.compute_sasa_trajectory, xtc.compute_sasa_trajectory_summary, XTC_FILE),
    "dcd": (dcd.compute_sasa_trajectory, dcd.compute_sasa_trajectory_summary, DCD_FILE),
}
RADII = np.full(304, 1.7, dtype=np.float32)


@pytest.fixture(params=sorted(READERS))
def reader_case(request: pytest.FixtureRequest):  # noqa: ANN201
    """(format name, reader class, a file of that format, a file of the other format)."""
    name = request.param
    reader_class, own, other = READERS[name]
    return name.upper(), reader_class, own, other


class TestOpenErrors:
    def test_missing_file_is_file_not_found(self, reader_case, tmp_path: Path) -> None:  # noqa: ANN001
        kind, reader_class, _, _ = reader_case
        missing = tmp_path / "missing.bin"

        with pytest.raises(FileNotFoundError, match=f"{kind} file not found"):
            reader_class(missing)

    def test_other_trajectory_format_is_value_error(self, reader_case) -> None:  # noqa: ANN001
        """An XTC reader on a DCD file (and the reverse): the file exists."""
        kind, reader_class, _, other = reader_case
        assert other.exists()

        with pytest.raises(ValueError, match=f"not a valid {kind} file") as excinfo:
            reader_class(other)

        assert str(other) in str(excinfo.value)

    def test_text_file_is_value_error(self, reader_case) -> None:  # noqa: ANN001
        kind, reader_class, _, _ = reader_case

        with pytest.raises(ValueError, match=f"not a valid {kind} file"):
            reader_class(PDB_FILE)

    def test_empty_file_is_value_error(self, reader_case, tmp_path: Path) -> None:  # noqa: ANN001
        kind, reader_class, _, _ = reader_case
        empty = tmp_path / "empty.bin"
        empty.write_bytes(b"")

        with pytest.raises(ValueError, match=f"not a valid {kind} file"):
            reader_class(empty)

    def test_truncated_header_is_value_error(self, reader_case, tmp_path: Path) -> None:  # noqa: ANN001
        kind, reader_class, own, _ = reader_case
        truncated = tmp_path / "truncated.bin"
        truncated.write_bytes(own.read_bytes()[:4])

        with pytest.raises(ValueError, match=f"not a valid {kind} file"):
            reader_class(truncated)

    def test_directory_is_is_a_directory_error(self, reader_case, tmp_path: Path) -> None:  # noqa: ANN001
        _, reader_class, _, _ = reader_case

        with pytest.raises(IsADirectoryError):
            reader_class(tmp_path)

    @pytest.mark.skipif(os.name == "nt" or os.geteuid() == 0, reason="needs POSIX permissions")
    def test_unreadable_file_is_permission_error(self, reader_case, tmp_path: Path) -> None:  # noqa: ANN001
        _, reader_class, own, _ = reader_case
        locked = tmp_path / "locked.bin"
        locked.write_bytes(own.read_bytes())
        locked.chmod(0)
        try:
            with pytest.raises(PermissionError):
                reader_class(locked)
        finally:
            locked.chmod(0o600)

    def test_valid_file_still_opens(self, reader_case) -> None:  # noqa: ANN001
        _, reader_class, own, _ = reader_case

        with reader_class(own) as reader:
            assert reader.natoms == 304


class TestReadErrors:
    def test_truncated_xtc_frame_is_runtime_error(self, tmp_path: Path) -> None:
        truncated = tmp_path / "truncated.xtc"
        truncated.write_bytes(XTC_FILE.read_bytes()[:200])

        with (
            xtc.XtcReader(truncated) as reader,
            pytest.raises(RuntimeError, match="corrupt or truncated"),
        ):
            list(reader)


@pytest.mark.parametrize("fmt", sorted(COMPUTE))
class TestFrameSelection:
    """step, start and stop are checked before the file is opened."""

    @pytest.mark.parametrize("function_index", [0, 1], ids=["trajectory", "summary"])
    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            ({"step": 0}, "step must be a positive integer.*got 0"),
            ({"step": -1}, "step must be a positive integer.*got -1"),
            ({"start": -5}, "start must be non-negative, got -5"),
            ({"stop": -1}, "stop must be non-negative"),
        ],
    )
    def test_bad_selection_is_value_error(
        self,
        fmt: str,
        function_index: int,
        kwargs: dict,
        message: str,
    ) -> None:
        function = COMPUTE[fmt][function_index]

        with pytest.raises(ValueError, match=message):
            function(COMPUTE[fmt][2], RADII, **kwargs)

    @pytest.mark.parametrize("function_index", [0, 1], ids=["trajectory", "summary"])
    def test_bad_selection_wins_over_a_missing_file(
        self,
        fmt: str,
        function_index: int,
        tmp_path: Path,
    ) -> None:
        """No work is done for a bad argument: the file is not even opened."""
        function = COMPUTE[fmt][function_index]

        with pytest.raises(ValueError, match="step must be a positive integer"):
            function(tmp_path / "missing.bin", RADII, step=0)

    def test_valid_selection_still_works(self, fmt: str) -> None:
        function, _, path = COMPUTE[fmt]

        result = function(path, RADII, start=2, stop=10, step=3, n_points=20)

        assert result.n_frames == 3

    def test_missing_file_is_file_not_found(self, fmt: str, tmp_path: Path) -> None:
        function = COMPUTE[fmt][0]

        with pytest.raises(FileNotFoundError):
            function(tmp_path / "missing.bin", RADII)

"""Thread-stress tests for Python FFI entry points."""

from __future__ import annotations

import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pytest

from zsasa import calculate_sasa, process_directory
from zsasa.xtc import XtcReader

TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"
PDB_FILE = TEST_DATA_DIR / "1l2y.pdb"
XTC_FILE = TEST_DATA_DIR / "1l2y.xtc"


def test_calculate_sasa_is_stable_from_multiple_python_threads() -> None:
    """Concurrent calculate_sasa calls should not crash or drift."""
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 3.0],
        ],
        dtype=np.float64,
    )
    radii = np.full(4, 1.7, dtype=np.float64)
    baseline = calculate_sasa(coords, radii, n_points=64, n_threads=1)

    def worker(_: int) -> tuple[float, np.ndarray]:
        last = baseline
        for _ in range(20):
            last = calculate_sasa(coords, radii, n_points=64, n_threads=2)
        return last.total_area, last.atom_areas.copy()

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(worker, range(4)))

    for total_area, atom_areas in results:
        assert total_area == pytest.approx(baseline.total_area)
        np.testing.assert_allclose(atom_areas, baseline.atom_areas)


def test_xtc_readers_can_open_and_read_concurrently() -> None:
    """Independent XtcReader handles should work from multiple Python threads."""
    if not XTC_FILE.exists():
        pytest.skip("XTC fixture not available")

    def worker(_: int) -> list[int]:
        steps: list[int] = []
        with XtcReader(XTC_FILE) as reader:
            for _ in range(3):
                frame = reader.read_frame()
                assert frame is not None
                steps.append(frame.step)
        return steps

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(worker, range(4)))

    assert results == [[1, 2, 3]] * 4


def test_process_directory_can_run_concurrently(tmp_path: Path) -> None:
    """Independent process_directory calls should be safe in parallel."""
    input_a = tmp_path / "a"
    input_b = tmp_path / "b"
    input_a.mkdir()
    input_b.mkdir()
    shutil.copy(PDB_FILE, input_a / "1l2y-a.pdb")
    shutil.copy(PDB_FILE, input_b / "1l2y-b.pdb")

    def worker(path: Path) -> tuple[int, int, int, float]:
        result = process_directory(path, n_threads=2)
        return result.total_files, result.successful, result.failed, result.total_sasa[0]

    with ThreadPoolExecutor(max_workers=2) as executor:
        result_a, result_b = list(executor.map(worker, [input_a, input_b]))

    assert result_a[:3] == (1, 1, 0)
    assert result_b[:3] == (1, 1, 0)
    assert result_a[3] == pytest.approx(result_b[3])

"""Thread-stress tests for Python FFI entry points."""

from __future__ import annotations

import shutil
import sys
import threading
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


_SIGPIPE_PROBE_LOCK = threading.Lock()


def _sigpipe_handler() -> int:
    """Address of the process's SIGPIPE handler (0 = default, 1 = ignored), via libc.

    Reads the handler by swapping it for SIG_IGN and back, so the swaps are serialized.
    """
    import ctypes

    libc = ctypes.CDLL(None)
    libc.signal.restype = ctypes.c_void_p
    libc.signal.argtypes = [ctypes.c_int, ctypes.c_void_p]
    sigpipe, sig_ignore = 13, 1
    with _SIGPIPE_PROBE_LOCK:
        current = libc.signal(sigpipe, sig_ignore) or 0
        libc.signal(sigpipe, current)
    return current


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX signal handlers")
def test_concurrent_process_directory_keeps_one_sigpipe_handler(tmp_path: Path) -> None:
    """Overlapping calls must not restore each other's handlers.

    The native batch code runs on a process-wide thread runtime that replaces the
    SIGPIPE handler once. A runtime created and torn down per call put the handler
    back to whatever it found, which with overlapping calls was another call's
    replacement, so the handler depended on the order in which calls finished.
    """
    input_dir = tmp_path / "in"
    input_dir.mkdir()
    shutil.copy(PDB_FILE, input_dir / "1l2y.pdb")

    process_directory(input_dir, n_threads=1)
    installed = _sigpipe_handler()
    assert installed not in (0, 1), "the batch runtime should have installed its handler"

    def worker(_: int) -> int:
        process_directory(input_dir, n_threads=2)
        return _sigpipe_handler()

    with ThreadPoolExecutor(max_workers=6) as executor:
        handlers = list(executor.map(worker, range(12)))

    assert set(handlers) == {installed}
    assert _sigpipe_handler() == installed

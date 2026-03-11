"""Tests for the CLI entry point."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"


def run_zsasa(*args: str) -> subprocess.CompletedProcess[str]:
    """Run zsasa CLI via the Python entry point."""
    return subprocess.run(
        [sys.executable, "-m", "zsasa.cli", *args],
        capture_output=True,
        text=True,
        timeout=30,
    )


def get_output(result: subprocess.CompletedProcess[str]) -> str:
    """Get combined stdout+stderr output (zsasa writes help/version to stderr)."""
    return result.stdout + result.stderr


class TestCLIEntryPoint:
    """Test that the CLI binary is bundled and executable."""

    def test_help(self):
        result = run_zsasa("--help")
        assert result.returncode == 0
        output = get_output(result)
        assert "USAGE" in output
        assert "calc" in output

    def test_version(self):
        result = run_zsasa("--version")
        assert result.returncode == 0
        assert "zsasa" in get_output(result)

    def test_calc_help(self):
        result = run_zsasa("calc", "--help")
        assert result.returncode == 0
        assert "SASA" in get_output(result)

    def test_calc_structure(self, tmp_path):
        input_file = EXAMPLES_DIR / "1ubq.cif"
        if not input_file.exists():
            pytest.skip("Example structure not available")

        output_file = tmp_path / "output.json"
        result = run_zsasa("calc", str(input_file), str(output_file))
        assert result.returncode == 0
        assert output_file.exists()
        assert output_file.stat().st_size > 0

    def test_binary_exists(self):
        from zsasa.cli import _find_binary

        binary = _find_binary()
        assert Path(binary).exists()

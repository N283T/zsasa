"""Tests for the CLI entry point."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).parent.parent.parent / "examples"
TEST_DATA_DIR = Path(__file__).parent.parent.parent / "test_data"


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


class TestSdfClassification:
    """An SDF/MOL molecule is classified from its own bond topology."""

    @staticmethod
    def calc_csv(tmp_path: Path, input_file: Path, *args: str) -> tuple[list[list[str]], str]:
        """Run `calc --format=csv`; return the CSV rows and the progress output."""
        output_file = tmp_path / f"{input_file.stem}.csv"
        result = run_zsasa("calc", "--format=csv", *args, str(input_file), str(output_file))
        assert result.returncode == 0, result.stderr
        rows = [line.split(",") for line in output_file.read_text().splitlines()]
        return rows, result.stderr

    @pytest.mark.parametrize("fixture", ["ethanol_v2000.sdf", "ethanol_v3000.sdf"])
    @pytest.mark.parametrize(
        ("args", "summary"),
        [
            ((), "Classifier 'CCD': 3 atoms classified, 0 fallback"),
            (("--include-hydrogens",), "Classifier 'CCD': 3 atoms classified, 6 fallback"),
        ],
    )
    def test_blank_title_gives_the_radii_of_the_titled_molecule(
        self, tmp_path: Path, fixture: str, args: tuple[str, ...], summary: str
    ):
        titled_file = TEST_DATA_DIR / fixture
        lines = titled_file.read_text().split("\n")
        assert lines[0] == "ethanol"
        blank_file = tmp_path / f"blank_{fixture}"
        blank_file.write_text("\n".join(["", *lines[1:]]))

        titled, titled_log = self.calc_csv(tmp_path, titled_file, *args)
        blank, blank_log = self.calc_csv(tmp_path, blank_file, *args)

        # Both are classified from the bond table, not by element
        assert summary in titled_log
        assert summary in blank_log

        # chain,residue,resnum,atom_name,x,y,z,radius,area: only the residue
        # name (the title) differs
        assert titled[0][1] == "residue"
        assert {row[1] for row in titled[1:-1]} == {"ethan"}
        assert {row[1] for row in blank[1:-1]} == {""}
        assert [row[:1] + row[2:] for row in titled] == [row[:1] + row[2:] for row in blank]
        assert [row[7] for row in blank[1:4]] == ["1.880", "1.880", "1.460"]

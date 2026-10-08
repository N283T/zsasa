"""Tests for the CLI entry point."""

from __future__ import annotations

import csv
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


RICH_CSV_HEADER = [
    "chain",
    "residue",
    "resnum",
    "insertion_code",
    "atom_name",
    "x",
    "y",
    "z",
    "radius",
    "area",
]


def calc_csv(tmp_path: Path, input_file: Path, *args: str) -> tuple[list[list[str]], str]:
    """Run `calc --format=csv`; return the CSV rows and the progress output."""
    output_file = tmp_path / f"{input_file.stem}.csv"
    result = run_zsasa("calc", "--format=csv", *args, str(input_file), str(output_file))
    assert result.returncode == 0, result.stderr
    with output_file.open(newline="") as f:
        rows = list(csv.reader(f))
    return rows, result.stderr


class TestSdfClassification:
    """An SDF/MOL molecule is classified from its own bond topology."""

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

        titled, titled_log = calc_csv(tmp_path, titled_file, *args)
        blank, blank_log = calc_csv(tmp_path, blank_file, *args)

        # Both are classified from the bond table, not by element
        assert summary in titled_log
        assert summary in blank_log

        # Only the residue name (the title) differs
        assert titled[0] == RICH_CSV_HEADER
        residue, radius = RICH_CSV_HEADER.index("residue"), RICH_CSV_HEADER.index("radius")
        assert {row[residue] for row in titled[1:-1]} == {"ethan"}
        assert {row[residue] for row in blank[1:-1]} == {""}
        assert [row[:residue] + row[residue + 1 :] for row in titled] == [
            row[:residue] + row[residue + 1 :] for row in blank
        ]
        assert [row[radius] for row in blank[1:4]] == ["1.880", "1.880", "1.460"]


# Residues 10, 10A and 10B of chain H (antibody numbering), far apart
INSERTION_CODE_PDB = """\
ATOM      1  N   GLY H  10       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY H  10       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  N   SER H  10A     20.000   0.000   0.000  1.00 20.00           N
ATOM      4  CA  SER H  10A     21.458   0.000   0.000  1.00 20.00           C
ATOM      5  N   THR H  10B     40.000   0.000   0.000  1.00 20.00           N
END
"""


class TestCsvOutput:
    """`--format=csv` for structure input."""

    def test_calc_has_an_insertion_code_column_after_resnum(self, tmp_path: Path):
        input_file = tmp_path / "insertion.pdb"
        input_file.write_text(INSERTION_CODE_PDB)

        rows, _ = calc_csv(tmp_path, input_file)

        assert rows[0] == RICH_CSV_HEADER
        assert [row[:5] for row in rows[1:-1]] == [
            ["H", "GLY", "10", "", "N"],
            ["H", "GLY", "10", "", "CA"],
            ["H", "SER", "10", "A", "N"],
            ["H", "SER", "10", "A", "CA"],
            ["H", "THR", "10", "B", "N"],
        ]
        assert all(len(row) == len(RICH_CSV_HEADER) for row in rows)
        # The total row leaves every column but the area empty
        assert rows[-1][:-1] == [""] * (len(RICH_CSV_HEADER) - 1)
        assert float(rows[-1][-1]) == pytest.approx(sum(float(row[-1]) for row in rows[1:-1]))

    def test_batch_writes_the_basic_csv_without_residue_columns(self, tmp_path: Path):
        """`batch` per-file CSV is `atom_index,area` for every input.

        It has no residue columns, so the insertion code column of the rich
        CSV does not apply to it.
        """
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        (input_dir / "insertion.pdb").write_text(INSERTION_CODE_PDB)
        output_dir = tmp_path / "out"

        result = run_zsasa("batch", "--format=csv", "--quiet", str(input_dir), str(output_dir))
        assert result.returncode == 0, result.stderr

        with (output_dir / "insertion.csv").open(newline="") as f:
            rows = list(csv.reader(f))
        assert rows[0] == ["atom_index", "area"]
        assert [row[0] for row in rows[1:]] == ["0", "1", "2", "3", "4", "total"]

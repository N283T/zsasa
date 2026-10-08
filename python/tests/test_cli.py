"""Tests for the CLI entry point."""

from __future__ import annotations

import csv
import json
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

    def test_fields_with_csv_metacharacters_survive_a_csv_parser(self, tmp_path: Path):
        """RFC 4180 quoting: a chain ID of `,` or `"` stays one field."""
        input_file = tmp_path / "metachars.pdb"
        input_file.write_text(
            "ATOM      1  N   GLY ,   1       0.000   0.000   0.000  1.00 20.00           N\n"
            'ATOM      2  N   GLY "   2      20.000   0.000   0.000  1.00 20.00           N\n'
            "ATOM      3  N   GLY A   3      40.000   0.000   0.000  1.00 20.00           N\n"
            "END\n"
        )

        rows, _ = calc_csv(tmp_path, input_file)

        assert rows[0] == RICH_CSV_HEADER
        assert all(len(row) == len(RICH_CSV_HEADER) for row in rows)
        assert [row[:5] for row in rows[1:-1]] == [
            [",", "GLY", "1", "", "N"],
            ['"', "GLY", "2", "", "N"],
            ["A", "GLY", "3", "", "N"],
        ]
        # A field that needs no quoting is written as before
        lines = (tmp_path / "metachars.csv").read_text().splitlines()
        assert lines[1].startswith('",",GLY,1,,N,')
        assert lines[2].startswith('"""",GLY,2,,N,')
        assert lines[3].startswith("A,GLY,3,,N,40.000,0.000,0.000,")

    def test_sdf_title_with_csv_metacharacters_survives_a_csv_parser(self, tmp_path: Path):
        """The residue name of an SDF molecule is the start of its title."""
        lines = (TEST_DATA_DIR / "ethanol_v2000.sdf").read_text().split("\n")
        input_file = tmp_path / "title.sdf"
        input_file.write_text("\n".join(['a,"b', *lines[1:]]))

        rows, _ = calc_csv(tmp_path, input_file)

        assert all(len(row) == len(RICH_CSV_HEADER) for row in rows)
        assert {row[1] for row in rows[1:-1]} == {'a,"b'}


def pdb_atom(serial: int, atom: str, residue: str, chain: str, number: str, x: float) -> str:
    """One ATOM record; `number` is the residue number with its insertion code."""
    resseq, icode = (number[:-1], number[-1]) if number[-1].isalpha() else (number, " ")
    element = atom[0]
    return (
        f"ATOM  {serial:5d}  {atom:<3s} {residue:>3s} {chain}{resseq:>4s}{icode}   "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00          {element:>2s}\n"
    )


# (chain, residue name, residue number, insertion code, atom count) per residue
ResidueKey = tuple[str, str, int, str, int]

RESIDUE_CASES: dict[str, tuple[str, list[ResidueKey]]] = {
    "insertion_codes": (
        INSERTION_CODE_PDB,
        [("H", "GLY", 10, "", 2), ("H", "SER", 10, "A", 2), ("H", "THR", 10, "B", 1)],
    ),
    # Same chain, number and insertion code, different residue names
    "same_number_different_names": (
        pdb_atom(1, "N", "GLY", "A", "10", 0.0)
        + pdb_atom(2, "CA", "GLY", "A", "10", 1.458)
        + pdb_atom(3, "N", "LYS", "A", "10", 20.0)
        + "END\n",
        [("A", "GLY", 10, "", 2), ("A", "LYS", 10, "", 1)],
    ),
    # ALA A 1 is interrupted by a residue of chain B: one entry per run
    "non_contiguous_residue": (
        pdb_atom(1, "N", "ALA", "A", "1", 0.0)
        + pdb_atom(2, "N", "GLY", "B", "2", 20.0)
        + pdb_atom(3, "CB", "ALA", "A", "1", 40.0)
        + "END\n",
        [("A", "ALA", 1, "", 1), ("B", "GLY", 2, "", 1), ("A", "ALA", 1, "", 1)],
    ),
    # All models are read superimposed by default: one entry per residue and model
    "two_models": (
        "MODEL        1\n"
        + pdb_atom(1, "N", "MET", "A", "1", 0.0)
        + pdb_atom(2, "CA", "MET", "A", "1", 1.458)
        + pdb_atom(3, "N", "GLY", "A", "2", 20.0)
        + "ENDMDL\nMODEL        2\n"
        + pdb_atom(1, "N", "MET", "A", "1", 0.5)
        + pdb_atom(2, "CA", "MET", "A", "1", 1.958)
        + pdb_atom(3, "N", "GLY", "A", "2", 20.5)
        + "ENDMDL\nEND\n",
        [
            ("A", "MET", 1, "", 2),
            ("A", "GLY", 2, "", 1),
            ("A", "MET", 1, "", 2),
            ("A", "GLY", 2, "", 1),
        ],
    ),
}


class TestResidueOutputsAgree:
    """`--per-residue`, `--format=rsa` and the JSONL residue map share one residue identity."""

    @staticmethod
    def per_residue_table(stderr: str) -> list[tuple[ResidueKey, float]]:
        """Rows of the `--per-residue` table: `Chain  Res    Num       SASA  Atoms`."""
        lines = stderr.splitlines()
        start = lines.index("Per-residue SASA:") + 3
        rows = []
        for line in lines[start:]:
            if not line.strip() or line.startswith("Output written"):
                break
            number = line[11:17].strip()
            icode = number[-1] if number[-1].isalpha() else ""
            key = (
                line[0:5].strip(),
                line[6:10].strip(),
                int(number.removesuffix(icode) if icode else number),
                icode,
                int(line[29:35]),
            )
            rows.append((key, float(line[18:28])))
        return rows

    @staticmethod
    def rsa_rows(text: str) -> list[tuple[tuple[str, str, int, str], float]]:
        """`RES` rows by the NACCESS fixed columns (as Bio.PDB.NACCESS slices them)."""
        rows = []
        for line in text.splitlines():
            if line.startswith("RES"):
                key = (line[8].strip(), line[4:7].strip(), int(line[9:13]), line[13].strip())
                rows.append((key, float(line[15:22])))
        return rows

    @pytest.mark.parametrize("case", list(RESIDUE_CASES))
    def test_same_residues_atom_counts_and_areas(self, tmp_path: Path, case: str):
        pdb_text, expected = RESIDUE_CASES[case]
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        input_file = input_dir / f"{case}.pdb"
        input_file.write_text(pdb_text)

        # --per-residue table on stderr
        result = run_zsasa("calc", "--per-residue", str(input_file), str(tmp_path / "out.json"))
        assert result.returncode == 0, result.stderr
        table = self.per_residue_table(result.stderr)

        # RSA file
        rsa_file = tmp_path / "out.rsa"
        result = run_zsasa("calc", "--quiet", "--format=rsa", str(input_file), str(rsa_file))
        assert result.returncode == 0, result.stderr
        assert "Warning" not in result.stderr
        rsa = self.rsa_rows(rsa_file.read_text())

        # JSONL residue map
        jsonl_file = tmp_path / "out.jsonl"
        result = run_zsasa(
            "batch",
            "--quiet",
            "--format=jsonl",
            "--residue-map",
            "-o",
            str(jsonl_file),
            str(input_dir),
        )
        assert result.returncode == 0, result.stderr
        row = json.loads(jsonl_file.read_text())
        residue_map = list(
            zip(
                row["residue_chain"],
                row["residue_name"],
                row["residue_number"],
                row["residue_insertion_code"],
                row["residue_atom_count"],
                strict=True,
            )
        )

        assert [key for key, _ in table] == expected
        assert residue_map == expected
        assert [key for key, _ in rsa] == [key[:4] for key in expected]

        # Areas: the table and the RSA file print two decimals
        for (_, table_area), (_, rsa_area), map_area in zip(
            table, rsa, row["residue_sasa"], strict=True
        ):
            assert table_area == pytest.approx(map_area, abs=0.0051)
            assert rsa_area == pytest.approx(map_area, abs=0.0051)
        assert sum(row["residue_sasa"]) == pytest.approx(row["total_area"])
        # Atom ranges of the residue map cover every atom exactly once, in order
        starts, counts = row["residue_atom_start"], row["residue_atom_count"]
        assert starts == [sum(counts[:i]) for i in range(len(counts))]
        assert sum(counts) == len(row["atom_areas"])


class TestPolarPartition:
    """The polar/non-polar split by atom follows the classes of the classifier."""

    @pytest.mark.parametrize("classifier", ["naccess", "oons", "ccd"])
    def test_rsa_totals_and_polar_summary_are_the_sums_by_class(
        self, tmp_path: Path, classifier: str
    ):
        from zsasa import AtomClass, ClassifierType, classify_atoms

        input_file = EXAMPLES_DIR / "1ubq.pdb"
        # The atoms that calc reads by default: ATOM records without hydrogens
        atom_lines = [
            line
            for line in input_file.read_text().splitlines()
            if line.startswith("ATOM") and line[76:78].strip() != "H"
        ]
        residues = [line[17:20].strip() for line in atom_lines]
        atom_names = [line[12:16].strip() for line in atom_lines]

        json_file = tmp_path / "out.json"
        result = run_zsasa(
            "calc", f"--classifier={classifier}", "--polar", str(input_file), str(json_file)
        )
        assert result.returncode == 0, result.stderr
        areas = json.loads(json_file.read_text())["atom_areas"]
        assert len(areas) == len(atom_lines)

        # Classes from the same classifier through the C API
        classes = classify_atoms(residues, atom_names, ClassifierType[classifier.upper()]).classes
        assert AtomClass.UNKNOWN not in set(classes)
        polar = sum(
            area for area, cls in zip(areas, classes, strict=True) if cls == AtomClass.POLAR
        )
        apolar = sum(
            area for area, cls in zip(areas, classes, strict=True) if cls == AtomClass.APOLAR
        )

        # TOTAL row of the RSA file: non-polar in columns 51-60, polar in columns 64-73
        rsa_file = tmp_path / "out.rsa"
        result_rsa = run_zsasa(
            "calc",
            "--quiet",
            "--format=rsa",
            f"--classifier={classifier}",
            str(input_file),
            str(rsa_file),
        )
        assert result_rsa.returncode == 0, result_rsa.stderr
        total_row = next(
            line for line in rsa_file.read_text().splitlines() if line.startswith("TOTAL")
        )
        assert float(total_row[50:60]) == pytest.approx(apolar, abs=0.051)
        assert float(total_row[63:73]) == pytest.approx(polar, abs=0.051)

        # Atom summary of --polar on stderr
        lines = result.stderr.splitlines()
        start = next(
            i
            for i, line in enumerate(lines)
            if line.startswith("Polar/Nonpolar SASA by atom class")
        )
        assert float(lines[start + 1].split()[1]) == pytest.approx(polar, abs=0.0051)
        assert float(lines[start + 2].split()[1]) == pytest.approx(apolar, abs=0.0051)
        assert lines[start + 1].endswith(
            f"- {sum(cls == AtomClass.POLAR for cls in classes)} atoms"
        )
        assert lines[start + 2].endswith(
            f"- {sum(cls == AtomClass.APOLAR for cls in classes)} atoms"
        )

    def test_classifiers_give_different_partitions(self, tmp_path: Path):
        """NACCESS classes sulfur as apolar, OONS classes carbonyl carbon as polar."""
        totals = {}
        for classifier in ("naccess", "oons", "ccd"):
            rsa_file = tmp_path / f"{classifier}.rsa"
            result = run_zsasa(
                "calc",
                "--quiet",
                "--format=rsa",
                f"--classifier={classifier}",
                str(EXAMPLES_DIR / "1ubq.pdb"),
                str(rsa_file),
            )
            assert result.returncode == 0, result.stderr
            total_row = next(
                line for line in rsa_file.read_text().splitlines() if line.startswith("TOTAL")
            )
            totals[classifier] = (float(total_row[50:60]), float(total_row[63:73]))

        # 1ubq, 100 test points: non-polar and polar area by class
        assert totals["naccess"] == pytest.approx((2469.4, 2353.9), abs=0.11)
        assert totals["oons"] == pytest.approx((2542.6, 2236.9), abs=0.11)
        assert totals["ccd"] == pytest.approx((2318.9, 2515.8), abs=0.11)

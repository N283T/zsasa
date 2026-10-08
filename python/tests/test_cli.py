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


class TestCLIEntryPoint:
    """Test that the CLI binary is found and executable."""

    def test_help(self):
        result = run_zsasa("--help")
        assert result.returncode == 0
        assert "USAGE" in result.stdout
        assert "calc" in result.stdout
        assert result.stderr == ""

    def test_version(self):
        from zsasa import get_version

        result = run_zsasa("--version")
        assert result.returncode == 0
        assert result.stdout == f"zsasa {get_version()}\n"
        assert result.stderr == ""

    @pytest.mark.parametrize("flag", ["--help", "-h"])
    @pytest.mark.parametrize("command", ["calc", "batch", "traj", "compile-dict"])
    def test_command_help_is_written_to_stdout(self, command: str, flag: str):
        result = run_zsasa(command, flag)
        assert result.returncode == 0
        assert f"{command}" in result.stdout
        assert "SASA" in result.stdout or "ZSDC" in result.stdout
        assert result.stderr == ""

    @pytest.mark.parametrize("args", [(), ("no-such-command",)])
    def test_usage_after_an_error_is_written_to_stderr(self, args: tuple[str, ...]):
        result = run_zsasa(*args)
        assert result.returncode == 1
        assert "USAGE" in result.stderr
        assert result.stdout == ""

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


# Two chains, one residue each.
TWO_CHAIN_PDB = (
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
    "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C\n"
    "ATOM      3  N   GLY B   1       5.000   0.000   0.000  1.00 20.00           N\n"
    "ATOM      4  CA  GLY B   1       6.500   0.000   0.000  1.00 20.00           C\n"
    "END\n"
)

# What every runner says about the two unreadable inputs of `write_inputs`.
TWO_FAILURES = (
    "  bad1.pdb: read/parse failed: NoAtomsFound\n  bad2.cif: read/parse failed: NoAtomSiteLoop\n"
)

# The failure report lists this many inputs and counts the rest.
MAX_LISTED_FAILURES = 20

THREADS = ["1", "4"]

# How a workflow is routed: every file parsed once and shared by the jobs, or
# one batch run per job. A job whose `auth_chain` differs from the shared
# setting selects the second runner for the whole workflow.
WORKFLOW_RUNNERS = {"file_first": "", "job_first": "auth_chain = true\n"}


def write_inputs(input_dir: Path, *, n_good: int = 2, n_extra_bad: int = 0) -> Path:
    """Fill a directory with readable structures and unreadable files."""
    input_dir.mkdir()
    for i in range(1, n_good + 1):
        (input_dir / f"good{i}.pdb").write_text(TWO_CHAIN_PDB)
    (input_dir / "bad1.pdb").write_text("not a structure\n")
    (input_dir / "bad2.cif").write_text("data_bad\nloop_\n_cell.length_a\n1.0\n")
    for i in range(n_extra_bad):
        (input_dir / f"worse{i:02d}.pdb").write_text("not a structure\n")
    return input_dir


def write_workflow(
    path: Path,
    runner: str,
    input_dir: Path,
    output_dir: Path | None,
    *,
    fmt: str = "jsonl",
    extra: str = "",
    first_job: str = '[[jobs]]\nname = "chain_a"\nchains = ["A"]\n\n',
) -> Path:
    """Write a workflow with the jobs `chain_a` and `everything` for a runner.

    Without an output directory the workflow has the one job `everything`.
    """
    output = f'dir = "{output_dir.as_posix()}"\n' if output_dir else ""
    path.write_text(
        "version = 1\n"
        'kind = "workflow"\n\n'
        f'[input]\ndir = "{input_dir.as_posix()}"\n\n'
        f'[output]\n{output}format = "{fmt}"\n\n'
        '[calculation]\nn_points = 8\n\n[classifier]\ntype = "naccess"\n\n'
        f"{extra}"
        f"{first_job if output_dir else ''}"
        f'[[jobs]]\nname = "everything"\n{WORKFLOW_RUNNERS[runner]}'
    )
    return path


def jsonl_rows(path: Path) -> list[str]:
    """Rows of a JSONL file in a fixed order (workers write as they finish)."""
    return sorted(path.read_text().splitlines())


class TestBatchFailureReporting:
    """Inputs that fail are reported on stderr, with or without `--quiet`."""

    @pytest.mark.parametrize("threads", THREADS)
    def test_quiet_reports_failed_inputs_of_per_file_output(self, tmp_path: Path, threads: str):
        input_dir = write_inputs(tmp_path / "in")
        output_dir = tmp_path / "out"

        result = run_zsasa("batch", "-q", f"--threads={threads}", str(input_dir), str(output_dir))

        # A failed input leaves no output file: the report is its only trace.
        assert sorted(p.name for p in output_dir.iterdir()) == ["good1.json", "good2.json"]
        assert result.stderr == "2 of 4 inputs failed:\n" + TWO_FAILURES
        assert result.stdout == ""
        # Partial failure is reported, not fatal.
        assert result.returncode == 0

    @pytest.mark.parametrize("threads", THREADS)
    def test_quiet_reports_failed_inputs_of_jsonl_output(self, tmp_path: Path, threads: str):
        input_dir = write_inputs(tmp_path / "in")
        output_file = tmp_path / "out.jsonl"

        result = run_zsasa(
            "batch",
            "-q",
            f"--threads={threads}",
            "--format=jsonl",
            str(input_dir),
            str(output_file),
        )

        assert result.returncode == 0
        assert result.stderr == "2 of 4 inputs failed:\n" + TWO_FAILURES
        rows = jsonl_rows(output_file)
        assert [row[:16] for row in rows] == ['{"status":"err",'] * 2 + ['{"status":"ok","'] * 2

    @pytest.mark.parametrize("threads", THREADS)
    @pytest.mark.parametrize("output_args", [("out",), ("--format=jsonl", "out.jsonl"), ()])
    def test_quiet_prints_nothing_when_every_input_succeeds(
        self, tmp_path: Path, threads: str, output_args: tuple[str, ...]
    ):
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        for name in ("a.pdb", "b.pdb", "c.pdb"):
            (input_dir / name).write_text(TWO_CHAIN_PDB)
        output_args = tuple(str(tmp_path / a) if not a.startswith("-") else a for a in output_args)

        result = run_zsasa("batch", "-q", f"--threads={threads}", str(input_dir), *output_args)

        assert result.returncode == 0
        assert result.stderr == ""
        assert result.stdout == ""

    @pytest.mark.parametrize("threads", THREADS)
    def test_every_input_failing_is_reported_too(self, tmp_path: Path, threads: str):
        input_dir = write_inputs(tmp_path / "in", n_good=0)

        result = run_zsasa(
            "batch", "-q", f"--threads={threads}", str(input_dir), str(tmp_path / "o")
        )

        assert result.stderr == "2 of 2 inputs failed:\n" + TWO_FAILURES
        # The exit status of a run without a single result is unchanged (#432).
        assert result.returncode == 0

    @pytest.mark.parametrize("threads", THREADS)
    def test_summary_lists_each_failed_input_once(self, tmp_path: Path, threads: str):
        input_dir = write_inputs(tmp_path / "in")

        result = run_zsasa("batch", f"--threads={threads}", str(input_dir), str(tmp_path / "out"))

        assert result.returncode == 0
        assert "  Failed:          2\n" in result.stderr
        # The report closes the summary instead of being printed next to it.
        assert result.stderr.endswith("\n\n2 of 4 inputs failed:\n" + TWO_FAILURES)
        assert result.stderr.count("bad1.pdb") == 1
        assert result.stderr.count("bad2.cif") == 1

    @pytest.mark.parametrize("threads", THREADS)
    @pytest.mark.parametrize("quiet", [("-q",), ()])
    def test_more_failures_than_the_report_lists(
        self, tmp_path: Path, threads: str, quiet: tuple[str, ...]
    ):
        n_failed = MAX_LISTED_FAILURES + 5
        input_dir = write_inputs(tmp_path / "in", n_extra_bad=n_failed - 2)
        listed = TWO_FAILURES + "".join(
            f"  worse{i:02d}.pdb: read/parse failed: NoAtomsFound\n"
            for i in range(MAX_LISTED_FAILURES - 2)
        )
        head = f"{n_failed} of {n_failed + 2} inputs failed:\n" + listed

        # Per-file output: the rest is counted.
        per_file = run_zsasa(
            "batch", *quiet, f"--threads={threads}", str(input_dir), str(tmp_path / "out")
        )
        assert per_file.returncode == 0
        assert per_file.stderr.endswith(head + "  ... and 5 more\n")
        assert ("Batch Results:" in per_file.stderr) == (not quiet)
        if quiet:
            assert per_file.stderr == head + "  ... and 5 more\n"

        # JSONL output: the last line says where every failure is.
        output_file = tmp_path / "out.jsonl"
        jsonl = run_zsasa(
            "batch",
            *quiet,
            f"--threads={threads}",
            "--format=jsonl",
            str(input_dir),
            "-o",
            str(output_file),
        )
        assert jsonl.stderr.endswith(
            head + f'  ... and 5 more (every failure is a "status":"err" row in {output_file})\n'
        )
        assert sum('"status":"err"' in row for row in jsonl_rows(output_file)) == n_failed

        # JSONL on standard output
        stdout = run_zsasa(
            "batch", *quiet, f"--threads={threads}", "--format=jsonl", str(input_dir)
        )
        assert stdout.stderr.endswith(
            head + '  ... and 5 more (every failure is a "status":"err" row in the JSONL output)\n'
        )
        assert stdout.stdout.count('"status":"err"') == n_failed

    @pytest.mark.parametrize("threads", THREADS)
    def test_empty_directory_leaves_an_empty_jsonl_file(self, tmp_path: Path, threads: str):
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        output_file = tmp_path / "out.jsonl"
        output_file.write_text('{"status":"ok","filename":"stale.pdb"}\n')

        result = run_zsasa(
            "batch",
            "-q",
            f"--threads={threads}",
            "--format=jsonl",
            str(input_dir),
            str(output_file),
        )

        assert result.returncode == 0
        assert result.stderr == ""
        assert output_file.read_text() == ""

    def test_empty_chain_list_is_rejected(self, tmp_path: Path):
        input_dir = write_inputs(tmp_path / "in")

        for value in (",", "", " , "):
            result = run_zsasa(
                "batch", "-q", f"--chain={value}", str(input_dir), str(tmp_path / "o")
            )
            assert result.returncode == 1
            assert "Error: --chain needs at least one chain ID" in result.stderr
            assert not (tmp_path / "o").exists()


class TestWorkflowFailureReporting:
    """A workflow reports failures the same way whichever runner it is routed to."""

    @staticmethod
    def run_workflow(workflow: Path, threads: str, *args: str) -> subprocess.CompletedProcess[str]:
        return run_zsasa("batch", "-q", f"--threads={threads}", *args, "--workflow", str(workflow))

    @pytest.mark.parametrize("threads", THREADS)
    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_reports_failed_inputs_per_job(self, tmp_path: Path, runner: str, threads: str):
        input_dir = write_inputs(tmp_path / "in")
        output_dir = tmp_path / "out"
        workflow = write_workflow(tmp_path / "wf.toml", runner, input_dir, output_dir)

        result = self.run_workflow(workflow, threads)

        assert result.returncode == 0
        assert result.stderr == (
            "Workflow complete: 4 successful, 4 failed\n"
            "Job 'chain_a': 2 of 4 inputs failed:\n" + TWO_FAILURES + "Job 'everything': "
            "2 of 4 inputs failed:\n" + TWO_FAILURES
        )
        row = '{{"status":"err","filename":"{}","error":"read/parse failed: {}"}}'
        for job in ("chain_a", "everything"):
            rows = jsonl_rows(output_dir / f"{job}.jsonl")
            assert rows[:2] == [
                row.format("bad1.pdb", "NoAtomsFound"),
                row.format("bad2.cif", "NoAtomSiteLoop"),
            ]
            assert len(rows) == 4
            assert rows[2].startswith('{"status":"ok","filename":"good1.pdb",')
            assert rows[3].startswith('{"status":"ok","filename":"good2.pdb",')

    @pytest.mark.parametrize("threads", THREADS)
    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_quiet_prints_only_the_totals_without_failures(
        self, tmp_path: Path, runner: str, threads: str
    ):
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        (input_dir / "a.pdb").write_text(TWO_CHAIN_PDB)
        (input_dir / "b.pdb").write_text(TWO_CHAIN_PDB)
        workflow = write_workflow(tmp_path / "wf.toml", runner, input_dir, tmp_path / "out")

        result = self.run_workflow(workflow, threads)

        assert result.returncode == 0
        assert result.stderr == "Workflow complete: 4 successful, 0 failed\n"

    @pytest.mark.parametrize("threads", THREADS)
    @pytest.mark.parametrize(
        ("fmt", "has_output_dir", "extra"),
        [
            ("jsonl", True, ""),
            ("jsonl", True, "[output.jsonl]\natom_areas = false\n\n"),
            ("json", True, ""),
            ("jsonl", False, ""),
            # Used to write nothing at all in the job-first runner.
            ("jsonl", False, "[output.jsonl]\natom_areas = false\n\n"),
        ],
    )
    def test_runners_agree_on_rows_stdout_and_exit_status(
        self, tmp_path: Path, threads: str, fmt: str, has_output_dir: bool, extra: str
    ):
        input_dir = write_inputs(tmp_path / "in")
        outcomes = {}
        for runner in WORKFLOW_RUNNERS:
            output_dir = tmp_path / f"out-{runner}" if has_output_dir else None
            workflow = write_workflow(
                tmp_path / f"{runner}.toml", runner, input_dir, output_dir, fmt=fmt, extra=extra
            )
            result = self.run_workflow(workflow, threads)
            files = {}
            if output_dir:
                for path in sorted(output_dir.rglob("*")):
                    if path.is_file():
                        key = path.relative_to(output_dir).as_posix()
                        files[key] = (
                            jsonl_rows(path) if path.suffix == ".jsonl" else path.read_text()
                        )
            outcomes[runner] = {
                "returncode": result.returncode,
                "stdout": sorted(result.stdout.splitlines()),
                "stderr": result.stderr,
                "files": files,
            }

        file_first, job_first = outcomes["file_first"], outcomes["job_first"]
        assert file_first == job_first
        assert file_first["returncode"] == 0
        n_jobs = 2 if has_output_dir else 1
        assert file_first["stderr"].startswith(
            f"Workflow complete: {2 * n_jobs} successful, {2 * n_jobs} failed\n"
        )

        # One row per input and job, wherever the rows go.
        if fmt == "json":
            assert sorted(file_first["files"]) == [
                "chain_a/good1.json",
                "chain_a/good2.json",
                "everything/good1.json",
                "everything/good2.json",
            ]
            assert file_first["stdout"] == []
        elif has_output_dir:
            assert sorted(file_first["files"]) == ["chain_a.jsonl", "everything.jsonl"]
            assert all(len(rows) == 4 for rows in file_first["files"].values())
            assert file_first["stdout"] == []
        else:
            assert len(file_first["stdout"]) == 4
            assert sum('"status":"err"' in row for row in file_first["stdout"]) == 2
        if fmt == "jsonl":
            rows = file_first["stdout"] or file_first["files"]["everything.jsonl"]
            assert all(('"atom_areas"' in row) == (not extra) for row in rows if '"ok"' in row)

    @pytest.mark.parametrize("threads", THREADS)
    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_more_failures_than_the_report_lists(self, tmp_path: Path, runner: str, threads: str):
        n_failed = MAX_LISTED_FAILURES + 5
        input_dir = write_inputs(tmp_path / "in", n_extra_bad=n_failed - 2)
        output_dir = tmp_path / "out"
        workflow = write_workflow(tmp_path / "wf.toml", runner, input_dir, output_dir)

        result = self.run_workflow(workflow, threads)

        assert result.returncode == 0
        lines = result.stderr.splitlines()
        assert lines[0] == f"Workflow complete: 4 successful, {2 * n_failed} failed"
        # Per job: a heading, the listed inputs and the line that counts the rest.
        assert len(lines) == 1 + 2 * (MAX_LISTED_FAILURES + 2)
        for job, start in (("chain_a", 1), ("everything", MAX_LISTED_FAILURES + 3)):
            assert lines[start] == f"Job '{job}': {n_failed} of {n_failed + 2} inputs failed:"
            assert lines[start + 1] == "  bad1.pdb: read/parse failed: NoAtomsFound"
            last_line = lines[start + MAX_LISTED_FAILURES + 1]
            assert last_line.startswith(
                '  ... and 5 more (every failure is a "status":"err" row in '
            )
            assert last_line.endswith(f"{job}.jsonl)")

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_missing_input_directory_is_an_error(self, tmp_path: Path, runner: str):
        workflow = write_workflow(
            tmp_path / "wf.toml", runner, tmp_path / "no-such-dir", tmp_path / "out"
        )

        result = self.run_workflow(workflow, "2")

        assert result.returncode == 1
        assert "FileNotFound" in result.stderr
        if runner == "job_first":
            # The job and the cause, apart from the per-input counts
            assert (
                "Error running workflow job 'chain_a': cannot read input directory "
                f"'{(tmp_path / 'no-such-dir').as_posix()}': FileNotFound\n"
            ) in result.stderr
            assert "Workflow complete: 0 successful, 0 failed\n" in result.stderr
            assert "2 of 2 jobs failed: chain_a, everything\n" in result.stderr
            assert result.stderr.endswith("Error: WorkflowJobFailed\n")

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_job_that_cannot_write_its_output_is_an_error(self, tmp_path: Path, runner: str):
        input_dir = write_inputs(tmp_path / "in")
        output_dir = tmp_path / "out"
        # A directory where the JSONL file of the second job belongs
        (output_dir / "everything.jsonl").mkdir(parents=True)
        workflow = write_workflow(tmp_path / "wf.toml", runner, input_dir, output_dir)

        result = self.run_workflow(workflow, "2")

        assert result.returncode == 1
        if runner == "job_first":
            assert "Error running workflow job 'everything': cannot create JSONL output" in (
                result.stderr
            )
            # The other job still ran; its failed inputs are not failed jobs.
            assert "Workflow complete: 2 successful, 2 failed\n" in result.stderr
            assert "Job 'chain_a': 2 of 4 inputs failed:\n" + TWO_FAILURES in result.stderr
            assert "1 of 2 jobs failed: everything\n" in result.stderr
            assert len(jsonl_rows(output_dir / "chain_a.jsonl")) == 4

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_colliding_output_names_are_an_error(self, tmp_path: Path, runner: str):
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        (input_dir / "same.pdb").write_text(TWO_CHAIN_PDB)
        (input_dir / "same.ent").write_text(TWO_CHAIN_PDB)
        output_dir = tmp_path / "out"
        workflow = write_workflow(tmp_path / "wf.toml", runner, input_dir, output_dir, fmt="json")

        result = self.run_workflow(workflow, "2")

        assert result.returncode == 1
        assert "same.json <- same.ent, same.pdb" in result.stderr
        assert not output_dir.exists()

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_empty_chains_array_is_rejected(self, tmp_path: Path, runner: str):
        input_dir = write_inputs(tmp_path / "in")
        output_dir = tmp_path / "out"
        workflow = write_workflow(
            tmp_path / "wf.toml",
            runner,
            input_dir,
            output_dir,
            first_job='[[jobs]]\nname = "nothing"\nchains = []\n\n',
        )

        result = self.run_workflow(workflow, "2")

        assert result.returncode == 1
        assert f"Error reading workflow file '{workflow}': EmptyJobChains\n" in result.stderr
        assert "a [[jobs]] entry has an empty chains array" in result.stderr
        assert not output_dir.exists()


BSA_ANALYSIS = '[analysis]\ntype = "bsa"\npartner_a = ["A"]\npartner_b = ["B"]\n'


def write_key_workflow(
    path: Path,
    input_dir: Path,
    output_dir: Path,
    *,
    input_extra: str = "",
    calc_extra: str = "",
    body: str = '[[jobs]]\nname = "everything"\n',
) -> Path:
    """Write a JSONL batch workflow with extra lines in `[input]` and `[calculation]`."""
    path.write_text(
        "version = 1\n"
        'kind = "workflow"\n\n'
        f'[input]\ndir = "{input_dir.as_posix()}"\n{input_extra}\n'
        f'[output]\ndir = "{output_dir.as_posix()}"\nformat = "jsonl"\n\n'
        f"[calculation]\nn_points = 8\n{calc_extra}\n"
        '[classifier]\ntype = "naccess"\n\n'
        f"{body}"
    )
    return path


class TestWorkflowIgnoredKeys:
    """Nothing a workflow says is dropped silently: an error or a warning."""

    @staticmethod
    def structure_dir(tmp_path: Path) -> Path:
        input_dir = tmp_path / "in"
        input_dir.mkdir()
        (input_dir / "two.pdb").write_text(TWO_CHAIN_PDB)
        return input_dir

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    @pytest.mark.parametrize(
        ("input_extra", "calc_extra", "key"),
        [
            ("model = 2\n", "", "[input] model"),
            ('mol = "1"\n', "", "[input] mol"),
            ('path = "two.pdb"\n', "", "[input] path"),
            ("", "rsa = true\n", "[calculation] rsa = true"),
            ("", "per_residue = true\n", "[calculation] per_residue = true"),
            ("", "polar = true\n", "[calculation] polar = true"),
            ("", "validate_only = true\n", "[calculation] validate_only = true"),
        ],
    )
    def test_batch_rejects_keys_it_cannot_honor(
        self, tmp_path: Path, runner: str, input_extra: str, calc_extra: str, key: str
    ):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            input_extra=input_extra,
            calc_extra=calc_extra,
            body=f'[[jobs]]\nname = "everything"\n{WORKFLOW_RUNNERS[runner]}',
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 1
        assert f"Error: {key}" in result.stderr
        assert "zsasa calc --workflow" in result.stderr
        assert not output_dir.exists()

    def test_batch_analysis_rejects_input_chain(self, tmp_path: Path):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            input_extra='chain = "A"\n',
            body=BSA_ANALYSIS,
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 1
        assert "Error: [input] chain is not supported by an [analysis] workflow" in result.stderr
        assert "partner_a and partner_b" in result.stderr
        assert not output_dir.exists()

    @pytest.mark.parametrize("flag", ["--residue-map", "--format=csv"])
    def test_batch_analysis_rejects_options_it_would_override(self, tmp_path: Path, flag: str):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            body=BSA_ANALYSIS,
        )

        result = run_zsasa("batch", "-q", flag, "--workflow", str(workflow))

        assert result.returncode == 1
        assert f"Error: {flag.split('=')[0]} " in result.stderr
        assert not output_dir.exists()

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    def test_input_chain_is_the_default_chains_of_jobs_without_any(
        self, tmp_path: Path, runner: str
    ):
        output_dir = tmp_path / "out"
        body = (
            '[[jobs]]\nname = "default"\n\n'
            f'[[jobs]]\nname = "explicit"\nchains = ["B"]\n{WORKFLOW_RUNNERS[runner]}\n'
            '[[jobs]]\nname = "both"\nchains = ["A", "B"]\n'
        )
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            input_extra='chain = "B"\n',
            body=body,
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 0, result.stderr
        assert "Warning" not in result.stderr
        default_rows = (output_dir / "default.jsonl").read_text()
        assert default_rows == (output_dir / "explicit.jsonl").read_text()
        assert default_rows != (output_dir / "both.jsonl").read_text()

    def test_input_chain_that_no_job_uses_is_rejected(self, tmp_path: Path):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            input_extra='chain = "B"\n',
            body='[[jobs]]\nname = "own"\nchains = ["A"]\n',
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 1
        assert "Error: [input] chain is used by no job" in result.stderr
        assert not output_dir.exists()

    def test_analysis_with_jobs_is_a_manifest_error(self, tmp_path: Path):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            body=BSA_ANALYSIS + '\n[[jobs]]\nname = "never"\n',
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 1
        assert f"Error reading workflow file '{workflow}': AnalysisWithJobs\n" in result.stderr
        assert "[analysis] and [[jobs]] cannot be combined" in result.stderr
        assert not output_dir.exists()

    @pytest.mark.parametrize(
        "jobs",
        [
            '[[jobs]]\n[[jobs]]\nname = "second"\n',
            '[[jobs]]\nname = "first"\n[[jobs]]\n',
            "[[jobs]]\n",
        ],
    )
    def test_jobs_table_without_a_name_is_rejected(self, tmp_path: Path, jobs: str):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml", self.structure_dir(tmp_path), output_dir, body=jobs
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 1
        assert f"Error reading workflow file '{workflow}': MissingJobName\n" in result.stderr
        assert not output_dir.exists()

    @pytest.mark.parametrize("runner", WORKFLOW_RUNNERS)
    @pytest.mark.parametrize("how", ["manifest", "option"])
    def test_timing_in_a_workflow_with_jobs_warns_and_runs(
        self, tmp_path: Path, runner: str, how: str
    ):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            calc_extra="timing = true\n" if how == "manifest" else "",
            body=f'[[jobs]]\nname = "everything"\n{WORKFLOW_RUNNERS[runner]}',
        )
        extra = ["--timing"] if how == "option" else []

        result = run_zsasa("batch", "-q", *extra, "--workflow", str(workflow))

        # -q keeps the warning, as it keeps the failure reports
        assert result.returncode == 0, result.stderr
        timing = "--timing" if how == "option" else "[calculation] timing = true"
        assert result.stderr.startswith(
            f"Warning: {timing} has no effect on a workflow with [[jobs]]"
        )
        assert (output_dir / "everything.jsonl").exists()

    def test_timing_in_an_analysis_workflow_is_honored(self, tmp_path: Path):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml",
            self.structure_dir(tmp_path),
            output_dir,
            calc_extra="timing = true\n",
            body=BSA_ANALYSIS,
        )

        result = run_zsasa("batch", "--workflow", str(workflow))

        assert result.returncode == 0, result.stderr
        assert "Warning" not in result.stderr
        assert "BSA analysis SASA time:" in result.stderr

    def test_profile_stages_in_a_workflow_warns(self, tmp_path: Path):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml", self.structure_dir(tmp_path), output_dir
        )

        result = run_zsasa("batch", "-q", "--profile-stages", "--workflow", str(workflow))

        assert result.returncode == 0, result.stderr
        assert result.stderr.startswith("Warning: --profile-stages has no effect on a workflow")
        assert (output_dir / "everything.jsonl").exists()

    def test_batch_output_path_warns_and_runs(self, tmp_path: Path):
        output_dir = tmp_path / "out"
        workflow = write_key_workflow(
            tmp_path / "wf.toml", self.structure_dir(tmp_path), output_dir
        )
        workflow.write_text(
            workflow.read_text().replace("[output]\n", '[output]\npath = "ignored.json"\n')
        )

        result = run_zsasa("batch", "-q", "--workflow", str(workflow))

        assert result.returncode == 0, result.stderr
        assert result.stderr.startswith(
            "Warning: [output] path is not read by 'zsasa batch --workflow'"
        )
        assert (output_dir / "everything.jsonl").exists()

    @staticmethod
    def write_calc_workflow(
        path: Path,
        structure: Path,
        result_file: Path,
        *,
        input_extra: str = "",
        output_extra: str = "",
        tail: str = "",
    ) -> Path:
        """Write a `calc` workflow with extra lines in `[input]` and `[output]`."""
        path.write_text(
            "version = 1\n"
            'kind = "workflow"\n\n'
            f'[input]\npath = "{structure.as_posix()}"\n{input_extra}\n'
            f'[output]\npath = "{result_file.as_posix()}"\nformat = "json"\n{output_extra}\n'
            "[calculation]\nn_points = 8\nquiet = true\n\n"
            '[classifier]\ntype = "naccess"\n\n'
            f"{tail}"
        )
        return path

    @pytest.mark.parametrize(
        ("input_extra", "tail", "message"),
        [
            (
                "",
                BSA_ANALYSIS,
                "Error: [analysis] is read only by 'zsasa batch --workflow'",
            ),
            (
                "",
                '[[jobs]]\nname = "all"\n',
                "Error: [[jobs]] is read only by 'zsasa batch --workflow'",
            ),
            (
                'dir = "structures"\n',
                "",
                "Error: [input] dir is read only by 'zsasa batch --workflow'",
            ),
        ],
    )
    def test_calc_rejects_batch_only_structure(
        self, tmp_path: Path, input_extra: str, tail: str, message: str
    ):
        structure = tmp_path / "two.pdb"
        structure.write_text(TWO_CHAIN_PDB)
        result_file = tmp_path / "result.json"
        workflow = self.write_calc_workflow(
            tmp_path / "wf.toml", structure, result_file, input_extra=input_extra, tail=tail
        )

        result = run_zsasa("calc", "--workflow", str(workflow))

        assert result.returncode == 1
        assert message in result.stderr
        assert "run it with 'zsasa batch --workflow'" in result.stderr
        assert not result_file.exists()

    def test_calc_warns_about_batch_output_keys_and_runs(self, tmp_path: Path):
        structure = tmp_path / "two.pdb"
        structure.write_text(TWO_CHAIN_PDB)
        result_file = tmp_path / "result.json"
        workflow = self.write_calc_workflow(
            tmp_path / "wf.toml",
            structure,
            result_file,
            output_extra='dir = "results"\n',
            tail="[output.jsonl]\ndecimals = 3\n",
        )

        result = run_zsasa("calc", "-q", "--workflow", str(workflow))

        assert result.returncode == 0, result.stderr
        assert "Warning: [output] dir is read only by 'zsasa batch --workflow'" in result.stderr
        assert "Warning: [output.jsonl] is read only by 'zsasa batch --workflow'" in result.stderr
        assert result_file.exists()
        assert not (tmp_path / "results").exists()


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

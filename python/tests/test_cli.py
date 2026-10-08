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

#!/usr/bin/env python3
import contextlib
import importlib.util
import io
import shutil
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str):
    path = ROOT.joinpath("scripts", name)
    spec = importlib.util.spec_from_file_location(name.replace(".py", ""), path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CheckVersionsTests(unittest.TestCase):
    def setUp(self):
        self.cv = load_script("check_versions.py")
        self.tmp = Path(tempfile.mkdtemp(prefix="zsasa-check-versions-"))
        self.addCleanup(lambda: shutil.rmtree(self.tmp, ignore_errors=True))
        self.write_repo(self.tmp, "0.10.0")

    @staticmethod
    def write_repo(root: Path, version: str, **overrides: str):
        versions = {
            "build.zig": version,
            "build.zig.zon": version,
            "src/c_api.zig": version,
            "python/pyproject.toml": version,
        }
        versions.update(overrides)
        root.joinpath("src").mkdir(exist_ok=True)
        root.joinpath("python").mkdir(exist_ok=True)
        root.joinpath("build.zig").write_text(
            f'const version = "{versions["build.zig"]}";\n'
        )
        # A dependency hash that contains a version-like string must not be picked up.
        root.joinpath("build.zig.zon").write_text(
            textwrap.dedent(f"""\
                .{{
                    .name = .zsasa,
                    .version = "{versions["build.zig.zon"]}",
                    .dependencies = .{{
                        .example = .{{ .hash = "example-0.1.0-abcdef" }},
                    }},
                }}
                """)
        )
        root.joinpath("src", "c_api.zig").write_text(
            f'const VERSION = "{versions["src/c_api.zig"]}";\n'
        )
        root.joinpath("python", "pyproject.toml").write_text(
            textwrap.dedent(f"""\
                [project]
                name = "zsasa"
                version = "{versions["python/pyproject.toml"]}"
                requires-python = ">=3.11"

                [tool.ty.environment]
                python-version = "3.11"
                """)
        )

    def test_matching_versions_pass_without_and_with_tag(self):
        self.assertEqual(self.cv.check(self.tmp), [])
        self.assertEqual(self.cv.check(self.tmp, "v0.10.0"), [])
        self.assertEqual(self.cv.read_versions(self.tmp)["build.zig.zon"], "0.10.0")

    def test_files_that_disagree_are_reported(self):
        self.write_repo(self.tmp, "0.10.0", **{"src/c_api.zig": "0.9.1"})
        problems = self.cv.check(self.tmp)
        self.assertEqual(len(problems), 1)
        self.assertIn("src/c_api.zig=0.9.1", problems[0])
        self.assertIn("build.zig=0.10.0", problems[0])

    def test_tag_mismatch_names_every_file_that_differs(self):
        problems = self.cv.check(self.tmp, "v0.10.1")
        self.assertEqual(len(problems), 4)
        self.assertTrue(all("0.10.0" in p and "0.10.1" in p for p in problems))

    def test_single_stale_file_is_caught_against_the_tag(self):
        self.write_repo(self.tmp, "0.10.0", **{"python/pyproject.toml": "0.9.1"})
        problems = self.cv.check(self.tmp, "v0.10.0")
        self.assertEqual(len(problems), 1)
        self.assertIn("python/pyproject.toml", problems[0])

    def test_tag_without_v_prefix_is_accepted_and_malformed_tag_rejected(self):
        self.assertEqual(self.cv.check(self.tmp, "0.10.0"), [])
        with self.assertRaises(ValueError):
            self.cv.check(self.tmp, "v0.10")

    def test_missing_file_and_missing_version_are_errors(self):
        self.tmp.joinpath("src", "c_api.zig").write_text("// no version here\n")
        with self.assertRaisesRegex(
            RuntimeError, "src/c_api.zig: could not find a version"
        ):
            self.cv.read_versions(self.tmp)
        self.tmp.joinpath("src", "c_api.zig").unlink()
        with self.assertRaisesRegex(RuntimeError, "src/c_api.zig: file not found"):
            self.cv.read_versions(self.tmp)

    def run_main(self, *argv: str) -> tuple[int, str, str]:
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = self.cv.main(["--root", str(self.tmp), *argv])
        return code, out.getvalue(), err.getvalue()

    def test_main_exit_codes(self):
        code, out, _ = self.run_main("--tag", "v0.10.0")
        self.assertEqual(code, 0)
        self.assertIn("OK: all versions agree and tag v0.10.0", out)

        code, _, err = self.run_main("--tag", "v0.11.0")
        self.assertEqual(code, 1)
        self.assertIn("does not match tag version 0.11.0", err)

        code, _, err = self.run_main("--tag", "not-a-tag")
        self.assertEqual(code, 1)
        self.assertIn("version must be", err)

    def test_repository_versions_agree(self):
        # The real checkout must be self-consistent; this is the same check CI runs.
        self.assertEqual(self.cv.check(ROOT), [])


if __name__ == "__main__":
    unittest.main()

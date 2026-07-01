#!/usr/bin/env python3
import importlib.util
import shutil
import subprocess
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


class ReleaseBumpTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="zsasa-release-test-"))
        self.addCleanup(lambda: shutil.rmtree(self.tmp, ignore_errors=True))
        self.write_minimal_repo(self.tmp)
        subprocess.run(["git", "init"], cwd=self.tmp, check=True, stdout=subprocess.DEVNULL)
        subprocess.run(["git", "config", "user.email", "test@example.invalid"], cwd=self.tmp, check=True)
        subprocess.run(["git", "config", "user.name", "Test User"], cwd=self.tmp, check=True)
        subprocess.run(["git", "add", "."], cwd=self.tmp, check=True)
        subprocess.run(["git", "commit", "-m", "initial"], cwd=self.tmp, check=True, stdout=subprocess.DEVNULL)

    def write_minimal_repo(self, root: Path):
        root.joinpath("src").mkdir()
        root.joinpath("python").mkdir()
        root.joinpath("packaging", "conda-forge").mkdir(parents=True)
        root.joinpath("website", "docs").mkdir(parents=True)
        root.joinpath("build.zig").write_text('const version = "0.7.1";\n')
        root.joinpath("build.zig.zon").write_text('    .version = "0.7.1",\n')
        root.joinpath("flake.nix").write_text('version = "0.7.1";\n')
        root.joinpath("python", "pyproject.toml").write_text('[project]\nname = "zsasa"\nversion = "0.7.1"\n')
        root.joinpath("python", "uv.lock").write_text('[[package]]\nname = "zsasa"\nversion = "0.7.1"\n')
        root.joinpath("packaging", "conda-forge", "meta.yaml").write_text('{% set version = "0.7.1" %}\n')
        root.joinpath("src", "c_api.zig").write_text('const VERSION = "0.7.1";\n')
        root.joinpath("CITATION.cff").write_text('version: "0.7.1"\ndate-released: "2026-06-29"\n')
        root.joinpath("CHANGELOG.md").write_text(textwrap.dedent("""\
            # Changelog

            ## [Unreleased]

            ### Changed

            - Something changed. (#1)

            ## [0.7.1] - 2026-06-29

            ### Fixed

            - Previous fix.

            [Unreleased]: https://github.com/N283T/zsasa/compare/v0.7.1...HEAD
            [0.7.1]: https://github.com/N283T/zsasa/compare/v0.7.0...v0.7.1
            """))
        root.joinpath("website", "docs", "changelog.md").write_text(textwrap.dedent("""\
            ---
            sidebar_position: 10
            ---

            # Changelog

            ## Unreleased

            ### Changed

            - Something changed. (#1)

            ## [v0.7.1](https://github.com/N283T/zsasa/releases/tag/v0.7.1) — 2026-06-29

            ### Fixed

            - Previous fix.
            """))

    def test_bump_updates_known_version_files_and_promotes_unreleased_notes(self):
        bump = load_script("release_bump.py")
        result = bump.run(self.tmp, "0.8.0", release_date="2026-07-01", check_clean=True)

        self.assertEqual(result.version, "0.8.0")
        self.assertEqual(result.tag, "v0.8.0")
        for rel in [
            "build.zig",
            "build.zig.zon",
            "flake.nix",
            "python/pyproject.toml",
            "python/uv.lock",
            "packaging/conda-forge/meta.yaml",
            "src/c_api.zig",
            "CITATION.cff",
        ]:
            self.assertIn("0.8.0", self.tmp.joinpath(rel).read_text(), rel)
            self.assertNotIn("0.7.1", self.tmp.joinpath(rel).read_text(), rel)

        changelog = self.tmp.joinpath("CHANGELOG.md").read_text()
        self.assertIn("## [Unreleased]\n\n## [0.8.0] - 2026-07-01", changelog)
        self.assertIn("- Something changed. (#1)", changelog)
        self.assertIn("[Unreleased]: https://github.com/N283T/zsasa/compare/v0.8.0...HEAD", changelog)
        self.assertIn("[0.8.0]: https://github.com/N283T/zsasa/compare/v0.7.1...v0.8.0", changelog)

        website = self.tmp.joinpath("website/docs/changelog.md").read_text()
        self.assertIn("## Unreleased\n\n## [v0.8.0](https://github.com/N283T/zsasa/releases/tag/v0.8.0) — 2026-07-01", website)
        self.assertIn("- Something changed. (#1)", website)

    def test_bump_rejects_dirty_tree_when_requested(self):
        bump = load_script("release_bump.py")
        self.tmp.joinpath("README.md").write_text("dirty\n")
        with self.assertRaisesRegex(RuntimeError, "working tree is not clean"):
            bump.run(self.tmp, "0.8.0", release_date="2026-07-01", check_clean=True)


class ReleaseTagTests(unittest.TestCase):
    def test_normalize_version_accepts_prefixed_and_unprefixed(self):
        tag = load_script("release_tag.py")
        self.assertEqual(tag.normalize_version("0.8.0"), ("0.8.0", "v0.8.0"))
        self.assertEqual(tag.normalize_version("v0.8.0"), ("0.8.0", "v0.8.0"))


if __name__ == "__main__":
    unittest.main()

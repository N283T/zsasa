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
        subprocess.run(
            ["git", "config", "user.email", "test@example.invalid"],
            cwd=self.tmp,
            check=True,
        )
        subprocess.run(["git", "config", "user.name", "Test User"], cwd=self.tmp, check=True)
        subprocess.run(["git", "add", "."], cwd=self.tmp, check=True)
        subprocess.run(
            ["git", "commit", "-m", "initial"],
            cwd=self.tmp,
            check=True,
            stdout=subprocess.DEVNULL,
        )

    def write_minimal_repo(self, root: Path):
        root.joinpath("src").mkdir()
        root.joinpath("python").mkdir()
        root.joinpath("website", "docs").mkdir(parents=True)
        root.joinpath("build.zig").write_text('const version = "0.7.1";\n')
        root.joinpath("build.zig.zon").write_text(
            textwrap.dedent("""\
                .version = "0.7.1",
                .dependencies = .{
                    .example = .{ .hash = "example-0.7.1-abcdef" },
                },
            """)
        )
        root.joinpath("flake.nix").write_text('version = "0.7.1";\n')
        root.joinpath("python", "pyproject.toml").write_text('[project]\nname = "zsasa"\nversion = "0.7.1"\n')
        root.joinpath("python", "uv.lock").write_text('[[package]]\nname = "zsasa"\nversion = "0.7.1"\n')
        root.joinpath("src", "c_api.zig").write_text('const VERSION = "0.7.1";\n')
        root.joinpath("CITATION.cff").write_text('version: "0.7.1"\ndate-released: "2026-06-29"\n')
        root.joinpath("CHANGELOG.md").write_text(
            textwrap.dedent("""\
            # Changelog

            ## [Unreleased]

            ### Changed

            - Something changed. (#1)

            ## [0.7.1] - 2026-06-29

            ### Fixed

            - Previous fix.

            [Unreleased]: https://github.com/N283T/zsasa/compare/v0.7.1...HEAD
            [0.7.1]: https://github.com/N283T/zsasa/compare/v0.7.0...v0.7.1
            """)
        )
        root.joinpath("website", "docs", "changelog.md").write_text(
            textwrap.dedent("""\
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
            """)
        )

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
            "src/c_api.zig",
            "CITATION.cff",
        ]:
            self.assertIn("0.8.0", self.tmp.joinpath(rel).read_text(), rel)
            if rel != "build.zig.zon":
                self.assertNotIn("0.7.1", self.tmp.joinpath(rel).read_text(), rel)

        self.assertIn("example-0.7.1-abcdef", self.tmp.joinpath("build.zig.zon").read_text())

        changelog = self.tmp.joinpath("CHANGELOG.md").read_text()
        self.assertIn("## [Unreleased]\n\n## [0.8.0] - 2026-07-01", changelog)
        self.assertIn("- Something changed. (#1)", changelog)
        self.assertIn(
            "[Unreleased]: https://github.com/N283T/zsasa/compare/v0.8.0...HEAD",
            changelog,
        )
        self.assertIn("[0.8.0]: https://github.com/N283T/zsasa/compare/v0.7.1...v0.8.0", changelog)

        website = self.tmp.joinpath("website/docs/changelog.md").read_text()
        self.assertIn(
            "## Unreleased\n\n## [v0.8.0](https://github.com/N283T/zsasa/releases/tag/v0.8.0) — 2026-07-01",
            website,
        )
        self.assertIn("- Something changed. (#1)", website)

    def test_bump_rejects_dirty_tree_when_requested(self):
        bump = load_script("release_bump.py")
        self.tmp.joinpath("README.md").write_text("dirty\n")
        with self.assertRaisesRegex(RuntimeError, "working tree is not clean"):
            bump.run(self.tmp, "0.8.0", release_date="2026-07-01", check_clean=True)


ZON = """\
.{
    .name = .example,
    // The package version changes with every release; the dependencies do not.
    .version = "VERSION",
    .dependencies = .{
        // a comment with a "quote
        .dep = .{
            .url = "git+https://example.invalid/dep.git?ref=main#abc",
            .hash = "dep-0.1.0-AAAA",
        },
        .other = .{ .path = "../other" },
    },
    .paths = .{ "build.zig" },
}
"""

FLAKE = """\
          zig = zig-overlay.packages.${system}."ZIG";
          zigDeps = pkgs.runCommand "deps"
            {
              outputHashMode = "recursive";
              outputHash = "sha256-BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB=";
            }
            "";
"""


class CheckNixDepsHashTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="zsasa-nix-test-"))
        self.addCleanup(lambda: shutil.rmtree(self.tmp, ignore_errors=True))
        self.script = load_script("check_nix_deps_hash.py")
        self.write(zon=ZON.replace("VERSION", "0.1.0"), flake=FLAKE.replace("ZIG", "0.16.0"))

    def write(self, *, zon: str | None = None, flake: str | None = None):
        if zon is not None:
            self.tmp.joinpath("build.zig.zon").write_text(zon)
        if flake is not None:
            self.tmp.joinpath("flake.nix").write_text(flake)

    def fingerprint(self) -> str:
        return self.script.compute_fingerprint(
            self.tmp.joinpath("build.zig.zon").read_text(),
            self.tmp.joinpath("flake.nix").read_text(),
        )

    def fake_build(self, real_hash: str):
        seen: list[str] = []

        def build(root: Path) -> str:
            seen.append(root.joinpath("flake.nix").read_text())
            specified = f"         specified: {self.script.FAKE_HASH}"
            return f"error: hash mismatch in fixed-output derivation\n{specified}\n            got:    {real_hash}\n"

        return build, seen

    def test_fingerprint_follows_the_dependencies_and_the_zig_version_only(self):
        zon = ZON.replace("VERSION", "0.1.0")
        base = self.fingerprint()
        # A new package version and edited comments leave the dependencies as they were.
        self.write(zon=ZON.replace("VERSION", "0.2.0").replace("The package version", "Edited comment"))
        self.assertEqual(self.fingerprint(), base)
        for old, new in [
            ("dep-0.1.0-AAAA", "dep-0.2.0-BBBB"),
            ("?ref=main#abc", "?ref=main#def"),
            ("../other", "../elsewhere"),
        ]:
            self.write(zon=zon.replace(old, new))
            self.assertNotEqual(self.fingerprint(), base, new)
        self.write(zon=zon, flake=FLAKE.replace("ZIG", "0.17.0"))
        self.assertNotEqual(self.fingerprint(), base)

    def test_dependency_urls_with_double_slash_are_not_cut_as_comments(self):
        lines = self.script.dependency_lines(ZON.replace("VERSION", "0.1.0"))
        self.assertIn('dep.url="git+https://example.invalid/dep.git?ref=main#abc"', lines)
        self.assertIn('other.path="../other"', lines)
        self.assertEqual(self.script.dependency_lines(".{ .name = .x }"), [])

    def test_check_needs_a_matching_fingerprint(self):
        self.assertFalse(self.script.check(self.tmp))  # nothing recorded yet
        flake = self.tmp.joinpath("flake.nix")
        flake.write_text(self.script.set_fingerprint(flake.read_text(), self.fingerprint()))
        self.assertTrue(self.script.check(self.tmp))
        self.write(zon=ZON.replace("VERSION", "0.1.0").replace("dep-0.1.0-AAAA", "dep-0.2.0-BBBB"))
        self.assertFalse(self.script.check(self.tmp))

    def test_refresh_records_the_hash_and_the_fingerprint(self):
        real = "sha256-" + "C" * 43 + "="
        build, seen = self.fake_build(real)
        self.assertEqual(self.script.refresh(self.tmp, build), real)

        self.assertIn(f'outputHash = "{self.script.FAKE_HASH}"', seen[0])  # built with the fake hash
        flake = self.tmp.joinpath("flake.nix").read_text()
        self.assertIn(f'outputHash = "{real}";', flake)
        self.assertIn(
            f"              # zig-deps-fingerprint: {self.fingerprint()}\n              outputHash",
            flake,
        )
        self.assertTrue(self.script.check(self.tmp))
        # A second refresh updates the marker in place instead of adding another.
        self.script.refresh(self.tmp, build)
        self.assertEqual(self.tmp.joinpath("flake.nix").read_text().count("zig-deps-fingerprint"), 1)

    def test_refresh_restores_flake_when_nix_gives_no_hash(self):
        original = self.tmp.joinpath("flake.nix").read_text()
        with self.assertRaisesRegex(RuntimeError, "did not report a hash mismatch"):
            self.script.refresh(self.tmp, lambda root: "error: attribute 'zsasa' missing\n")
        self.assertEqual(self.tmp.joinpath("flake.nix").read_text(), original)

        def interrupted(root: Path) -> str:
            raise KeyboardInterrupt

        with self.assertRaises(KeyboardInterrupt):
            self.script.refresh(self.tmp, interrupted)
        self.assertEqual(self.tmp.joinpath("flake.nix").read_text(), original)

    def test_the_flake_in_this_repository_is_up_to_date_with_build_zig_zon(self):
        self.assertTrue(
            self.script.check(ROOT),
            "build.zig.zon dependencies changed since flake.nix outputHash was computed: run scripts/check_nix_deps_hash.py --refresh",
        )


if __name__ == "__main__":
    unittest.main()

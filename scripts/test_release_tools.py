#!/usr/bin/env python3
import hashlib
import importlib.util
import json
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CONDA_RECIPE = """\
{{% set name = "zsasa" %}}
{{% set version = "{version}" %}}

source:
  - url: https://github.com/N283T/zsasa/releases/download/v{{{{ version }}}}/zsasa-{{{{ version }}}}-linux-x86_64  # [linux and x86_64]
    sha256: {sha}  # [linux and x86_64]
    fn: zsasa  # [linux and x86_64]
  - url: https://github.com/N283T/zsasa/releases/download/v{{{{ version }}}}/zsasa-{{{{ version }}}}-macos-aarch64  # [osx and arm64]
    sha256: {sha}  # [osx and arm64]
    fn: zsasa  # [osx and arm64]
  - url: https://raw.githubusercontent.com/N283T/zsasa/v{{{{ version }}}}/LICENSE
    sha256: {sha}
    fn: LICENSE

build:
  number: 3
"""


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
        root.joinpath("packaging", "conda-forge").mkdir(parents=True)
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
        root.joinpath("packaging", "conda-forge", "meta.yaml").write_text(CONDA_RECIPE.format(version="0.7.1", sha="a" * 64))
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
            "packaging/conda-forge/meta.yaml",
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

    def test_bump_marks_conda_checksums_pending_and_resets_build_number(self):
        bump = load_script("release_bump.py")
        bump.run(self.tmp, "0.8.0", release_date="2026-07-01", check_clean=True)

        recipe = self.tmp.joinpath("packaging", "conda-forge", "meta.yaml").read_text()
        self.assertEqual(recipe.count(bump.PENDING_CHECKSUM), 3)
        self.assertNotIn("a" * 64, recipe)
        # The selector comments after the checksums survive.
        self.assertIn(f"sha256: {bump.PENDING_CHECKSUM}  # [linux and x86_64]", recipe)
        self.assertIn("number: 0", recipe)

    def test_bump_rejects_dirty_tree_when_requested(self):
        bump = load_script("release_bump.py")
        self.tmp.joinpath("README.md").write_text("dirty\n")
        with self.assertRaisesRegex(RuntimeError, "working tree is not clean"):
            bump.run(self.tmp, "0.8.0", release_date="2026-07-01", check_clean=True)


LICENSE_TEXT = b"MIT License\n"
TARGETS = (
    "linux-x86_64",
    "linux-aarch64",
    "macos-x86_64",
    "macos-aarch64",
    "windows-x86_64.exe",
)
PACKAGING_FILES = ("packaging/conda-forge/meta.yaml",)


def fake_sha(name: str) -> str:
    return hashlib.sha256(name.encode()).hexdigest()


def fake_release(version: str, *, api_digests: bool = True, tamper: dict[str, str] | None = None):
    """Return (fetch, expected) imitating the files of a published release; no network."""
    names = [f"zsasa-{version}-{target}" for target in TARGETS]
    expected = {name: fake_sha(name) for name in names}
    sums = "".join(f"{expected[name]}  {name}\n" for name in names)
    api = {
        "assets": [
            {
                "name": name,
                "digest": f"sha256:{expected[name]}" if api_digests else None,
            }
            for name in names
        ]
    }
    for name, digest in (tamper or {}).items():
        api["assets"][names.index(name)]["digest"] = f"sha256:{digest}"
    files = {
        f"https://github.com/N283T/zsasa/releases/download/v{version}/SHA256SUMS": sums.encode(),
        f"https://api.github.com/repos/N283T/zsasa/releases/tags/v{version}": json.dumps(api).encode(),
        f"https://raw.githubusercontent.com/N283T/zsasa/v{version}/LICENSE": LICENSE_TEXT,
    }

    def fetch(url: str) -> bytes:
        try:
            return files[url]
        except KeyError:
            raise OSError(f"404 {url}") from None

    return fetch, expected


class UpdatePackagingChecksumsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="zsasa-packaging-test-"))
        self.addCleanup(lambda: shutil.rmtree(self.tmp, ignore_errors=True))
        self.script = load_script("update_packaging_checksums.py")
        self.tmp.joinpath("build.zig").write_text('const version = "0.8.0";\n')
        self.tmp.joinpath("packaging", "conda-forge").mkdir(parents=True)
        self.write_packaging("0.7.1", "b" * 64)

    def write_packaging(self, version: str, sha: str):
        self.tmp.joinpath("packaging", "conda-forge", "meta.yaml").write_text(CONDA_RECIPE.format(version=version, sha=sha))

    def read(self, rel: str) -> str:
        return self.tmp.joinpath(rel).read_text()

    def test_update_rewrites_the_conda_recipe(self):
        fetch, expected = fake_release("0.8.0")
        stale = self.script.run(self.tmp, None, fetch=fetch)

        self.assertEqual(sorted(stale), sorted(PACKAGING_FILES))
        recipe = self.read("packaging/conda-forge/meta.yaml")
        self.assertIn('{% set version = "0.8.0" %}', recipe)
        self.assertIn(
            f"sha256: {expected['zsasa-0.8.0-linux-x86_64']}  # [linux and x86_64]",
            recipe,
        )
        self.assertIn(
            f"sha256: {expected['zsasa-0.8.0-macos-aarch64']}  # [osx and arm64]",
            recipe,
        )
        self.assertIn(
            f"sha256: {hashlib.sha256(LICENSE_TEXT).hexdigest()}\n    fn: LICENSE",
            recipe,
        )
        self.assertNotIn("b" * 64, recipe)
        # Assets the recipe does not use are ignored.
        self.assertNotIn(expected["zsasa-0.8.0-windows-x86_64.exe"], recipe)


    def test_update_is_idempotent_and_check_agrees(self):
        fetch, _ = fake_release("0.8.0")
        self.script.run(self.tmp, "0.8.0", fetch=fetch)
        self.assertEqual(self.script.run(self.tmp, "0.8.0", fetch=fetch), [])
        self.assertEqual(self.script.run(self.tmp, "v0.8.0", check=True, fetch=fetch), [])

    def test_check_reports_stale_files_without_writing(self):
        fetch, _ = fake_release("0.8.0")
        before = {rel: self.read(rel) for rel in PACKAGING_FILES}
        stale = self.script.run(self.tmp, "0.8.0", check=True, fetch=fetch)
        self.assertEqual(sorted(stale), sorted(before))
        for rel, text in before.items():
            self.assertEqual(self.read(rel), text, rel)

    def test_check_detects_one_wrong_checksum_at_the_right_version(self):
        fetch, expected = fake_release("0.8.0")
        self.script.run(self.tmp, "0.8.0", fetch=fetch)
        path = self.tmp.joinpath("packaging", "conda-forge", "meta.yaml")
        path.write_text(path.read_text().replace(expected["zsasa-0.8.0-macos-aarch64"], "c" * 64))
        self.assertEqual(
            self.script.run(self.tmp, "0.8.0", check=True, fetch=fetch),
            ["packaging/conda-forge/meta.yaml"],
        )

    def test_cross_check_rejects_a_digest_that_disagrees(self):
        fetch, _ = fake_release("0.8.0", tamper={"zsasa-0.8.0-linux-x86_64": "d" * 64})
        before = self.read("packaging/conda-forge/meta.yaml")
        with self.assertRaisesRegex(RuntimeError, "asset digest"):
            self.script.run(self.tmp, "0.8.0", fetch=fetch)
        self.assertEqual(self.read("packaging/conda-forge/meta.yaml"), before)

    def test_cross_check_requires_digests_unless_disabled(self):
        fetch, expected = fake_release("0.8.0", api_digests=False)
        with self.assertRaisesRegex(RuntimeError, "cannot be cross-checked"):
            self.script.run(self.tmp, "0.8.0", fetch=fetch)
        self.script.run(self.tmp, "0.8.0", fetch=fetch, cross_check=False)
        self.assertIn(expected["zsasa-0.8.0-linux-x86_64"], self.read("packaging/conda-forge/meta.yaml"))

    def test_missing_asset_checksum_fails_before_any_file_is_written(self):
        fetch, expected = fake_release("0.8.0")
        sums = "".join(f"{digest}  {name}\n" for name, digest in expected.items() if "macos-aarch64" not in name)

        def partial(url: str) -> bytes:
            return sums.encode() if url.endswith("/SHA256SUMS") else fetch(url)

        before = {rel: self.read(rel) for rel in PACKAGING_FILES}
        with self.assertRaisesRegex(RuntimeError, "no checksum for zsasa-0.8.0-macos-aarch64"):
            self.script.run(self.tmp, "0.8.0", fetch=partial, cross_check=False)
        for rel, text in before.items():
            self.assertEqual(self.read(rel), text, rel)

    def test_unpublished_release_is_reported(self):
        def offline(url: str) -> bytes:
            raise OSError("HTTP Error 404: Not Found")

        with self.assertRaisesRegex(RuntimeError, "is the release published"):
            self.script.run(self.tmp, "0.8.0", fetch=offline)

    def test_parse_sha256sums(self):
        digest = "e" * 64
        self.assertEqual(
            self.script.parse_sha256sums(f"{digest}  a\n\n{digest.upper()} *b c\n"),
            {"a": digest, "b c": digest},
        )
        with self.assertRaisesRegex(RuntimeError, "line 1"):
            self.script.parse_sha256sums("not a checksum line\n")
        with self.assertRaisesRegex(RuntimeError, "twice"):
            self.script.parse_sha256sums(f"{digest}  a\n{'f' * 64}  a\n")
        with self.assertRaisesRegex(RuntimeError, "empty"):
            self.script.parse_sha256sums("\n")

    def test_bump_then_update_goes_from_pending_to_published_checksums(self):
        bump = load_script("release_bump.py")
        self.write_packaging("0.7.1", "b" * 64)
        # Minimal files the bump touches besides the packaging ones.
        for rel, text in {
            "build.zig": 'const version = "0.7.1";\n',
            "build.zig.zon": '.version = "0.7.1",\n',
            "flake.nix": 'version = "0.7.1";\n',
            "python/pyproject.toml": 'version = "0.7.1"\n',
            "python/uv.lock": 'name = "zsasa"\nversion = "0.7.1"\n',
            "src/c_api.zig": 'const VERSION = "0.7.1";\n',
            "CITATION.cff": 'version: "0.7.1"\ndate-released: "2026-06-29"\n',
            "CHANGELOG.md": "## [Unreleased]\n\n- A change.\n\n[Unreleased]: https://github.com/N283T/zsasa/compare/v0.7.1...HEAD\n",
            "website/docs/changelog.md": "## Unreleased\n\n- A change.\n",
        }.items():
            path = self.tmp.joinpath(rel)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        bump.run(self.tmp, "0.8.0", release_date="2026-07-01", check_clean=False)

        # The recipe cannot pass for the published release until the checksums are filled in.
        fetch, _ = fake_release("0.8.0")
        self.assertEqual(
            self.script.run(self.tmp, "0.8.0", check=True, fetch=fetch),
            list(PACKAGING_FILES),
        )

        self.script.run(self.tmp, "0.8.0", fetch=fetch)
        self.assertNotIn(bump.PENDING_CHECKSUM, self.read("packaging/conda-forge/meta.yaml"))
        self.assertEqual(self.script.run(self.tmp, "0.8.0", check=True, fetch=fetch), [])

    def test_the_real_packaging_files_are_understood(self):
        for rel in PACKAGING_FILES:
            self.tmp.joinpath(rel).write_text(ROOT.joinpath(rel).read_text())
        version = self.script.read_current_version(ROOT)
        fetch, expected = fake_release(version)
        self.script.run(self.tmp, version, fetch=fetch)
        recipe = self.read("packaging/conda-forge/meta.yaml")
        for target in (
            "linux-x86_64",
            "linux-aarch64",
            "macos-x86_64",
            "macos-aarch64",
        ):
            self.assertIn(expected[f"zsasa-{version}-{target}"], recipe, target)


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

#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Bump zsasa release metadata and promote changelog notes."""

from __future__ import annotations

import argparse
import datetime as dt
import re
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

REPO = "https://github.com/N283T/zsasa"
VERSION_RE = re.compile(r"^v?(\d+\.\d+\.\d+)$")


class BumpResult(NamedTuple):
    version: str
    tag: str
    previous_version: str
    changed_files: tuple[str, ...]
    stale_refs: tuple[str, ...]


def normalize_version(raw: str) -> tuple[str, str]:
    match = VERSION_RE.match(raw.strip())
    if not match:
        raise ValueError(f"version must be X.Y.Z or vX.Y.Z, got {raw!r}")
    version = match.group(1)
    return version, f"v{version}"


def run_cmd(args: list[str], cwd: Path, *, capture: bool = False) -> str:
    result = subprocess.run(
        args,
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )
    return result.stdout if capture else ""


def require_clean_tree(root: Path) -> None:
    status = run_cmd(["git", "status", "--porcelain"], root, capture=True)
    if status.strip():
        raise RuntimeError("working tree is not clean; commit/stash changes before bumping a release")


def read_current_version(root: Path) -> str:
    text = root.joinpath("build.zig").read_text()
    match = re.search(r'const version = "(\d+\.\d+\.\d+)";', text)
    if not match:
        raise RuntimeError("could not determine current version from build.zig")
    return match.group(1)


def replace_once(path: Path, old: str, new: str) -> bool:
    text = path.read_text()
    if old not in text:
        raise RuntimeError(f"{path}: expected text not found: {old!r}")
    path.write_text(text.replace(old, new, 1))
    return old != new


def bump_fixed_version_files(root: Path, old: str, new: str, release_date: str) -> list[str]:
    replacements = {
        "build.zig": [(f'const version = "{old}";', f'const version = "{new}";')],
        "build.zig.zon": [(f'.version = "{old}",', f'.version = "{new}",')],
        "flake.nix": [(f'version = "{old}";', f'version = "{new}";')],
        "python/pyproject.toml": [(f'version = "{old}"', f'version = "{new}"')],
        "python/uv.lock": [(f'name = "zsasa"\nversion = "{old}"', f'name = "zsasa"\nversion = "{new}"')],
        "packaging/conda-forge/meta.yaml": [(f'{{% set version = "{old}" %}}', f'{{% set version = "{new}" %}}')],
        "src/c_api.zig": [(f'const VERSION = "{old}";', f'const VERSION = "{new}";')],
        "CITATION.cff": [(f'version: "{old}"', f'version: "{new}"'), (r'date-released: "', r'date-released: "')],
    }
    changed: list[str] = []
    for rel, pairs in replacements.items():
        path = root.joinpath(rel)
        text = path.read_text()
        original = text
        if rel == "CITATION.cff":
            text = text.replace(f'version: "{old}"', f'version: "{new}"', 1)
            text = re.sub(r'date-released: "\d{4}-\d{2}-\d{2}"', f'date-released: "{release_date}"', text, count=1)
            if text == original:
                raise RuntimeError(f"{path}: citation version/date replacement did not change file")
            path.write_text(text)
        else:
            for old_text, new_text in pairs:
                if old_text not in text:
                    raise RuntimeError(f"{path}: expected text not found: {old_text!r}")
                text = text.replace(old_text, new_text, 1)
            path.write_text(text)
        if path.read_text() != original:
            changed.append(rel)
    return changed


def find_next_heading(lines: list[str], start: int, prefix: str) -> int:
    for idx in range(start + 1, len(lines)):
        if lines[idx].startswith(prefix):
            return idx
    return len(lines)


def promote_changelog(root: Path, version: str, tag: str, old: str, release_date: str, *, allow_empty_notes: bool) -> list[str]:
    changed: list[str] = []
    path = root.joinpath("CHANGELOG.md")
    text = path.read_text()
    if f"## [{version}]" in text:
        raise RuntimeError(f"CHANGELOG.md already contains section for {version}")
    lines = text.splitlines(keepends=True)
    try:
        start = next(i for i, line in enumerate(lines) if line.rstrip() == "## [Unreleased]")
    except StopIteration as exc:
        raise RuntimeError("CHANGELOG.md missing ## [Unreleased] section") from exc
    end = find_next_heading(lines, start, "## [")
    notes = "".join(lines[start + 1 : end]).strip()
    if not notes and not allow_empty_notes:
        raise RuntimeError("CHANGELOG.md Unreleased section is empty; add release notes or pass --allow-empty-notes")
    new_section = f"## [Unreleased]\n\n## [{version}] - {release_date}\n\n{notes}\n\n"
    updated = "".join(lines[:start]) + new_section + "".join(lines[end:])
    unreleased_old = f"[Unreleased]: {REPO}/compare/v{old}...HEAD"
    unreleased_new = f"[Unreleased]: {REPO}/compare/{tag}...HEAD"
    version_link = f"[{version}]: {REPO}/compare/v{old}...{tag}"
    if unreleased_old in updated:
        updated = updated.replace(unreleased_old, f"{unreleased_new}\n{version_link}", 1)
    elif "[Unreleased]:" in updated:
        updated = re.sub(r"^\[Unreleased\]: .*$", f"{unreleased_new}\n{version_link}", updated, count=1, flags=re.MULTILINE)
    else:
        updated = updated.rstrip() + f"\n\n{unreleased_new}\n{version_link}\n"
    path.write_text(updated)
    changed.append("CHANGELOG.md")

    website_path = root.joinpath("website", "docs", "changelog.md")
    website = website_path.read_text()
    if f"## [v{version}]" in website:
        raise RuntimeError(f"website/docs/changelog.md already contains section for v{version}")
    website_lines = website.splitlines(keepends=True)
    try:
        w_start = next(i for i, line in enumerate(website_lines) if line.rstrip() == "## Unreleased")
    except StopIteration as exc:
        raise RuntimeError("website/docs/changelog.md missing ## Unreleased section") from exc
    w_end = find_next_heading(website_lines, w_start, "## [")
    website_notes = "".join(website_lines[w_start + 1 : w_end]).strip() or notes
    if not website_notes and not allow_empty_notes:
        raise RuntimeError("website changelog Unreleased section is empty; add release notes or pass --allow-empty-notes")
    website_section = f"## Unreleased\n\n## [v{version}]({REPO}/releases/tag/{tag}) — {release_date}\n\n{website_notes}\n\n"
    website_path.write_text("".join(website_lines[:w_start]) + website_section + "".join(website_lines[w_end:]))
    changed.append("website/docs/changelog.md")
    return changed


def collect_stale_refs(root: Path, old: str) -> tuple[str, ...]:
    active_paths = [
        "build.zig",
        "build.zig.zon",
        "flake.nix",
        "python/pyproject.toml",
        "python/uv.lock",
        "packaging/conda-forge/meta.yaml",
        "src/c_api.zig",
        "CITATION.cff",
    ]
    stale: list[str] = []
    pattern = re.compile(rf"\b{re.escape(old)}\b")
    for rel in active_paths:
        path = root.joinpath(rel)
        for number, line in enumerate(path.read_text().splitlines(), start=1):
            if pattern.search(line):
                stale.append(f"{rel}:{number}:{line}")
    return tuple(stale)


def run(
    root: Path,
    version_arg: str,
    *,
    release_date: str | None = None,
    check_clean: bool = True,
    allow_empty_notes: bool = False,
) -> BumpResult:
    root = root.resolve()
    version, tag = normalize_version(version_arg)
    release_date = release_date or dt.date.today().isoformat()
    if check_clean:
        require_clean_tree(root)
    old = read_current_version(root)
    if old == version:
        raise RuntimeError(f"release version is already {version}")
    changed = bump_fixed_version_files(root, old, version, release_date)
    changed.extend(promote_changelog(root, version, tag, old, release_date, allow_empty_notes=allow_empty_notes))
    stale = collect_stale_refs(root, old)
    if stale:
        raise RuntimeError("stale active version references remain:\n" + "\n".join(stale))
    return BumpResult(version=version, tag=tag, previous_version=old, changed_files=tuple(changed), stale_refs=stale)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="Release version, X.Y.Z or vX.Y.Z")
    parser.add_argument("--date", dest="release_date", help="Release date, YYYY-MM-DD (default: today)")
    parser.add_argument("--allow-dirty", action="store_true", help="Allow running with a dirty git tree")
    parser.add_argument("--allow-empty-notes", action="store_true", help="Allow empty Unreleased changelog notes")
    args = parser.parse_args(argv)
    try:
        result = run(
            Path.cwd(),
            args.version,
            release_date=args.release_date,
            check_clean=not args.allow_dirty,
            allow_empty_notes=args.allow_empty_notes,
        )
    except Exception as exc:  # noqa: BLE001 - CLI should print concise errors.
        print(f"release-bump: {exc}", file=sys.stderr)
        return 1
    print(f"Bumped zsasa {result.previous_version} -> {result.version} ({result.tag})")
    print("Changed files:")
    for rel in result.changed_files:
        print(f"  {rel}")
    print("Note: packaging/aur/PKGBUILD is intentionally not bumped before release asset checksums exist.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Check that the version in every build-time source file agrees.

    python3 scripts/check_versions.py              # the files must agree with each other
    python3 scripts/check_versions.py --tag v0.10.0  # ...and with the release tag

The checked files are the ones whose version ends up in a published artifact:

- ``build.zig``           (CLI ``--version``)
- ``build.zig.zon``       (Zig package version)
- ``src/c_api.zig``       (``zsasa_version()``, Python ``get_version()``)
- ``python/pyproject.toml`` (wheel and sdist metadata)

CI runs it without ``--tag`` on every change; the publish workflow runs it with the
pushed tag before anything is built, so a tag that does not match the sources is
rejected before any artifact is published. Exits 1 on any mismatch, missing file or
unparsable version.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# release_bump.py is the single source of truth for how a release version is
# spelled and for how build.zig is parsed; import it rather than copying it.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from release_bump import normalize_version, read_current_version  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

# (path relative to the repository root, regex with one group for the version).
# The patterns mirror the ones release_bump.py rewrites. build.zig is read with
# release_bump.read_current_version instead of a pattern here.
PATTERNS: tuple[tuple[str, str], ...] = (
    ("build.zig.zon", r'^\s*\.version = "([^"]+)",$'),
    ("src/c_api.zig", r'^const VERSION = "([^"]+)";$'),
    ("python/pyproject.toml", r'^version = "([^"]+)"$'),
)


def read_versions(root: Path) -> dict[str, str]:
    """Return {relative path: version} for every checked file.

    Raises RuntimeError if a file is missing or contains no version.
    """
    versions = {"build.zig": read_current_version(root)}
    for rel, pattern in PATTERNS:
        path = root.joinpath(rel)
        if not path.is_file():
            raise RuntimeError(f"{rel}: file not found")
        match = re.search(pattern, path.read_text(), flags=re.MULTILINE)
        if not match:
            raise RuntimeError(f"{rel}: could not find a version")
        versions[rel] = match.group(1)
    return versions


def check(root: Path, tag: str | None = None) -> list[str]:
    """Return a list of problems; an empty list means everything agrees."""
    versions = read_versions(root)
    problems: list[str] = []
    expected: str | None = None
    if tag is not None:
        expected, _ = normalize_version(tag)
        problems.extend(
            f"{rel}: version {found} does not match tag version {expected}"
            for rel, found in versions.items()
            if found != expected
        )
    elif len(set(versions.values())) > 1:
        problems.append(
            "versions disagree: "
            + ", ".join(f"{rel}={found}" for rel, found in versions.items())
        )
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", help="Release tag (vX.Y.Z) every file must match")
    parser.add_argument(
        "--root",
        type=Path,
        default=ROOT,
        help="Repository root (default: this checkout)",
    )
    args = parser.parse_args(argv)
    try:
        problems = check(args.root, args.tag)
        versions = read_versions(args.root)
    except (RuntimeError, ValueError) as exc:
        print(f"check-versions: {exc}", file=sys.stderr)
        return 1
    for rel, found in versions.items():
        print(f"{rel}: {found}")
    if problems:
        for problem in problems:
            print(f"check-versions: {problem}", file=sys.stderr)
        return 1
    suffix = f" and tag {args.tag}" if args.tag else ""
    print(f"OK: all versions agree{suffix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

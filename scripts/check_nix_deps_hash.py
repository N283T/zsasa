#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Check (or refresh) the Nix fixed-output hash of the Zig dependencies in flake.nix.

``flake.nix`` fetches the dependencies of ``build.zig.zon`` in a fixed-output
derivation whose ``outputHash`` has to be updated by hand whenever those
dependencies (or the Zig version) change. Nothing fails until somebody builds
the flake, so ``flake.nix`` also records a fingerprint of the inputs the hash
was computed for:

    # zig-deps-fingerprint: <sha256>
    outputHash = "sha256-...";

Without options this script recomputes the fingerprint from ``build.zig.zon``
(its ``.dependencies`` block, not the package version) and the Zig version
named in ``flake.nix``, and exits 1 when it differs from the recorded one. That
needs no Nix and no network. The fingerprint cannot prove the hash is right; it
only proves the hash was refreshed after the dependencies last changed.

``--refresh`` computes the hash with ``nix build`` (it builds with a fake hash
and reads the ``got:`` value from the mismatch error) and records the hash and
the fingerprint.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

FAKE_HASH = "sha256-AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA="
FLAKE_ATTR = ".#zsasa"

_HASH_RE = re.compile(r'(outputHash = ")(sha256-[A-Za-z0-9+/]+=*)(";)')
_FINGERPRINT_RE = re.compile(r"^(\s*# zig-deps-fingerprint: )([0-9a-f]{64})$", re.MULTILINE)
_ZIG_VERSION_RE = re.compile(r'zig-overlay\.packages\.\$\{system\}\."([^"]+)"')

# Runs `nix build` for the flake in `root` and returns its combined output. Replaced in tests.
NixBuild = Callable[[Path], str]


def strip_comments(zon: str) -> str:
    """Remove ``//`` comments from Zon source, leaving string literals (URLs contain ``//``) alone."""
    out: list[str] = []
    in_string = False
    i = 0
    while i < len(zon):
        c = zon[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < len(zon):
                out.append(zon[i + 1])
                i += 1
            elif c == '"':
                in_string = False
        elif c == '"':
            in_string = True
            out.append(c)
        elif zon.startswith("//", i):
            while i < len(zon) and zon[i] != "\n":
                i += 1
            continue
        else:
            out.append(c)
        i += 1
    return "".join(out)


def matching_brace(text: str, open_index: int) -> int:
    """Return the index of the ``}`` closing the ``{`` at ``open_index`` (strings are skipped)."""
    depth = 0
    in_string = False
    i = open_index
    while i < len(text):
        c = text[i]
        if in_string:
            if c == "\\":
                i += 1
            elif c == '"':
                in_string = False
        elif c == '"':
            in_string = True
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise RuntimeError("unbalanced braces in build.zig.zon")


def dependency_lines(zon: str) -> list[str]:
    """Return one canonical line per dependency field of ``.dependencies``, sorted."""
    zon = strip_comments(zon)
    start = re.search(r"\.dependencies\s*=\s*\.\{", zon)
    if not start:
        return []  # no dependencies: nothing to fetch
    block_open = start.end() - 1
    block = zon[block_open + 1 : matching_brace(zon, block_open)]
    lines: list[str] = []
    for entry in re.finditer(r'\.(@"[^"]+"|\w+)\s*=\s*\.\{', block):
        name = entry.group(1)
        body = block[entry.end() : matching_brace(block, entry.end() - 1)]
        for field in re.finditer(r'\.(url|hash|path|lazy)\s*=\s*("(?:[^"\\]|\\.)*"|true|false)', body):
            lines.append(f"{name}.{field.group(1)}={field.group(2)}")
    return sorted(lines)


def zig_version(flake: str) -> str:
    match = _ZIG_VERSION_RE.search(flake)
    if not match:
        raise RuntimeError('flake.nix: could not find the Zig version (zig-overlay.packages.${system}."X.Y.Z")')
    return match.group(1)


def compute_fingerprint(zon: str, flake: str) -> str:
    canonical = "\n".join([f"zig={zig_version(flake)}", *dependency_lines(zon)]) + "\n"
    return hashlib.sha256(canonical.encode()).hexdigest()


def recorded_fingerprint(flake: str) -> str | None:
    match = _FINGERPRINT_RE.search(flake)
    return match.group(2) if match else None


def set_output_hash(flake: str, new_hash: str) -> str:
    if len(_HASH_RE.findall(flake)) != 1:
        raise RuntimeError('flake.nix: expected exactly one outputHash = "sha256-...";')
    return _HASH_RE.sub(lambda m: f"{m.group(1)}{new_hash}{m.group(3)}", flake, count=1)


def set_fingerprint(flake: str, fingerprint: str) -> str:
    if _FINGERPRINT_RE.search(flake):
        return _FINGERPRINT_RE.sub(lambda m: f"{m.group(1)}{fingerprint}", flake, count=1)
    # First use: put the marker on the line above outputHash, with the same indentation.
    match = re.search(r"^([ \t]*)outputHash = ", flake, re.MULTILINE)
    if not match:
        raise RuntimeError("flake.nix: no outputHash line to attach the fingerprint to")
    marker = f"{match.group(1)}# zig-deps-fingerprint: {fingerprint}\n"
    return flake[: match.start()] + marker + flake[match.start() :]


def nix_build(root: Path) -> str:
    result = subprocess.run(
        ["nix", "build", FLAKE_ATTR, "--no-link"],
        cwd=root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return result.stdout


def check(root: Path) -> bool:
    """Return True when the recorded fingerprint matches build.zig.zon and flake.nix."""
    flake = root.joinpath("flake.nix").read_text()
    zon = root.joinpath("build.zig.zon").read_text()
    return recorded_fingerprint(flake) == compute_fingerprint(zon, flake)


def refresh(root: Path, build: NixBuild = nix_build) -> str:
    """Recompute ``outputHash`` with ``nix build`` and record it with the fingerprint; return the hash."""
    flake_path = root.joinpath("flake.nix")
    original = flake_path.read_text()
    fingerprint = compute_fingerprint(root.joinpath("build.zig.zon").read_text(), original)
    flake_path.write_text(set_output_hash(original, FAKE_HASH))
    try:
        output = build(root)
        match = re.search(r"got:\s+(sha256-[A-Za-z0-9+/]+=*)", output)
        if not match:
            tail = "\n".join(output.splitlines()[-15:])
            raise RuntimeError(f"nix build did not report a hash mismatch, so the hash is unknown. Output:\n{tail}")
    except BaseException:
        flake_path.write_text(original)
        raise
    flake_path.write_text(set_fingerprint(set_output_hash(original, match.group(1)), fingerprint))
    return match.group(1)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--refresh",
        action="store_true",
        help="Recompute outputHash with `nix build` and record it",
    )
    args = parser.parse_args(argv)
    root = Path.cwd()
    try:
        if args.refresh:
            new_hash = refresh(root)
            print(f"flake.nix: outputHash = {new_hash}; fingerprint recorded.")
            print("Verify with: nix build && ./result/bin/zsasa --version")
            return 0
        if check(root):
            print("flake.nix: the Zig dependency hash is up to date with build.zig.zon.")
            return 0
    except Exception as exc:  # noqa: BLE001 - CLI should print concise errors.
        print(f"check-nix-deps-hash: {exc}", file=sys.stderr)
        return 2
    print(
        "flake.nix: the Zig dependency hash is stale: the dependencies in build.zig.zon (or the Zig version in flake.nix)\n"
        "changed since outputHash was computed. Run: scripts/check_nix_deps_hash.py --refresh",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Fill the conda-forge recipe and the AUR package files with the checksums of a published release.

The checksums of the release binaries only exist once the publish workflow has
run, so this is a post-release step: ``scripts/release_bump.py`` writes the new
version into ``packaging/conda-forge/meta.yaml`` and marks its checksums as
pending, and this script replaces them with the published values. It also brings
``packaging/aur/PKGBUILD`` and ``packaging/aur/.SRCINFO`` to the same version.

The checksums come from the ``SHA256SUMS`` asset of the release and are
cross-checked against the digests GitHub records for the uploaded assets.
``--check`` changes nothing and exits non-zero when the files are stale.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import NamedTuple

REPO_SLUG = "N283T/zsasa"
RELEASE_DOWNLOAD_URL = f"https://github.com/{REPO_SLUG}/releases/download"
RELEASE_API_URL = f"https://api.github.com/repos/{REPO_SLUG}/releases/tags"
RAW_URL = f"https://raw.githubusercontent.com/{REPO_SLUG}"

VERSION_RE = re.compile(r"^v?(\d+\.\d+\.\d+)$")

CONDA_RECIPE = "packaging/conda-forge/meta.yaml"
AUR_PKGBUILD = "packaging/aur/PKGBUILD"
AUR_SRCINFO = "packaging/aur/.SRCINFO"

# Fetches a URL and returns the response body. Replaced in tests.
Fetch = Callable[[str], bytes]


class ReleaseChecksums(NamedTuple):
    assets: dict[str, str]  # asset name -> sha256
    license: str  # sha256 of LICENSE at the release tag


def normalize_version(raw: str) -> str:
    match = VERSION_RE.match(raw.strip())
    if not match:
        raise ValueError(f"version must be X.Y.Z or vX.Y.Z, got {raw!r}")
    return match.group(1)


def read_current_version(root: Path) -> str:
    text = root.joinpath("build.zig").read_text()
    match = re.search(r'const version = "(\d+\.\d+\.\d+)";', text)
    if not match:
        raise RuntimeError("could not determine current version from build.zig")
    return match.group(1)


# ---------------------------------------------------------------------------
# Published checksums
# ---------------------------------------------------------------------------
def http_fetch(url: str) -> bytes:
    request = urllib.request.Request(
        url, headers={"User-Agent": "zsasa-update-packaging-checksums"}
    )
    token = os.environ.get("GITHUB_TOKEN")
    if token and url.startswith("https://api.github.com/"):
        request.add_header("Authorization", f"Bearer {token}")
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def parse_sha256sums(text: str) -> dict[str, str]:
    """Parse ``sha256sum`` output (``<hex>  <name>`` or ``<hex> *<name>``)."""
    sums: dict[str, str] = {}
    for number, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        match = re.fullmatch(r"([0-9a-fA-F]{64}) [ *](.+)", line.rstrip())
        if not match:
            raise RuntimeError(
                f"SHA256SUMS line {number} is not '<sha256>  <name>': {line!r}"
            )
        digest, name = match.group(1).lower(), match.group(2)
        if sums.setdefault(name, digest) != digest:
            raise RuntimeError(
                f"SHA256SUMS lists {name} twice with different checksums"
            )
    if not sums:
        raise RuntimeError("SHA256SUMS is empty")
    return sums


def parse_api_digests(body: bytes) -> dict[str, str]:
    """Return the ``sha256:`` digests GitHub records for the assets of a release."""
    digests: dict[str, str] = {}
    for asset in json.loads(body).get("assets", []):
        digest = asset.get("digest") or ""
        if digest.startswith("sha256:"):
            digests[asset["name"]] = digest.removeprefix("sha256:").lower()
    return digests


def load_release_checksums(
    version: str, fetch: Fetch, *, cross_check: bool = True
) -> ReleaseChecksums:
    tag = f"v{version}"
    try:
        sums = parse_sha256sums(
            fetch(f"{RELEASE_DOWNLOAD_URL}/{tag}/SHA256SUMS").decode()
        )
        license_sha = hashlib.sha256(fetch(f"{RAW_URL}/{tag}/LICENSE")).hexdigest()
        digests = (
            parse_api_digests(fetch(f"{RELEASE_API_URL}/{tag}")) if cross_check else {}
        )
    except OSError as exc:  # urllib.error.URLError and HTTPError are OSErrors
        raise RuntimeError(
            f"could not fetch the published checksums of {tag} ({exc}); is the release published?"
        ) from exc
    if cross_check:
        for name, digest in sums.items():
            if name not in digests:
                raise RuntimeError(
                    f"GitHub records no sha256 digest for {name}, so SHA256SUMS cannot be cross-checked; "
                    "pass --no-cross-check to trust SHA256SUMS alone"
                )
            if digests[name] != digest:
                raise RuntimeError(
                    f"{name}: SHA256SUMS says {digest} but GitHub's asset digest is {digests[name]}"
                )
    return ReleaseChecksums(assets=sums, license=license_sha)


def lookup(checksums: ReleaseChecksums, key: str) -> str:
    """Return the checksum of a release asset, or of LICENSE."""
    if key == "LICENSE":
        return checksums.license
    try:
        return checksums.assets[key]
    except KeyError:
        known = ", ".join(sorted(checksums.assets))
        raise RuntimeError(
            f"the release has no checksum for {key} (SHA256SUMS lists: {known})"
        ) from None


# ---------------------------------------------------------------------------
# conda-forge recipe
# ---------------------------------------------------------------------------
_CONDA_VERSION_RE = re.compile(r'^\{% set version = "[^"]*" %\}$', re.MULTILINE)
_CONDA_URL_RE = re.compile(r"^\s*-?\s*url:\s*(.*?)\s*(?:#.*)?$")
_CONDA_SHA_RE = re.compile(r"^(\s*sha256:\s*)(\S+)(.*)$", re.DOTALL)


def _conda_source_key(url: str) -> str:
    """Name an entry of the recipe's ``source:`` list by the file it downloads."""
    if url.endswith("/LICENSE"):
        return "LICENSE"
    name = url.rsplit("/", 1)[-1]
    if "/releases/download/" not in url:
        raise RuntimeError(f"{CONDA_RECIPE}: unrecognized source url {url!r}")
    return name


def update_conda_recipe(text: str, version: str, checksums: ReleaseChecksums) -> str:
    if not _CONDA_VERSION_RE.search(text):
        raise RuntimeError(f'{CONDA_RECIPE}: no {{% set version = "..." %}} line')
    text = _CONDA_VERSION_RE.sub(f'{{% set version = "{version}" %}}', text, count=1)
    out: list[str] = []
    key: str | None = None
    entries = 0
    for line in text.splitlines(keepends=True):
        url_match = _CONDA_URL_RE.match(line)
        if url_match:
            key = _conda_source_key(
                url_match.group(1).replace("{{ version }}", version)
            )
        sha_match = _CONDA_SHA_RE.match(line)
        if sha_match:
            if key is None:
                raise RuntimeError(
                    f"{CONDA_RECIPE}: sha256 line without a preceding url: {line!r}"
                )
            line = f"{sha_match.group(1)}{lookup(checksums, key)}{sha_match.group(3)}"
            key = None
            entries += 1
        out.append(line)
    if entries == 0:
        raise RuntimeError(f"{CONDA_RECIPE}: no sha256 entries found")
    return "".join(out)


# ---------------------------------------------------------------------------
# AUR package
# ---------------------------------------------------------------------------
_PKGVER_RE = re.compile(r"^pkgver=(\S+)$", re.MULTILINE)
_PKGREL_RE = re.compile(r"^pkgrel=(\S+)$", re.MULTILINE)
_ARRAY_RE = r"^{name}=\((?P<body>[^)]*)\)"


def _pkgbuild_array(text: str, name: str) -> re.Match[str]:
    match = re.search(_ARRAY_RE.format(name=name), text, re.MULTILINE)
    if not match:
        raise RuntimeError(f"{AUR_PKGBUILD}: no {name}=(...) array")
    return match


def _quoted_items(body: str) -> list[str]:
    return re.findall(r"""["']([^"']*)["']""", body)


def pkgbuild_sources(text: str, version: str) -> list[str]:
    """Return the source URLs of a PKGBUILD with ``${pkgver}`` expanded."""
    items = _quoted_items(_pkgbuild_array(text, "source").group("body"))
    return [
        item.replace("${pkgver}", version).replace("$pkgver", version) for item in items
    ]


def update_pkgbuild(text: str, version: str, checksums: ReleaseChecksums) -> str:
    match = _PKGVER_RE.search(text)
    if not match or not _PKGREL_RE.search(text):
        raise RuntimeError(f"{AUR_PKGBUILD}: missing pkgver= or pkgrel= line")
    if match.group(1) != version:
        text = _PKGVER_RE.sub(f"pkgver={version}", text, count=1)
        text = _PKGREL_RE.sub(
            "pkgrel=1", text, count=1
        )  # a new upstream version restarts the package release
    sources = pkgbuild_sources(text, version)
    sums = [lookup(checksums, url.rsplit("/", 1)[-1]) for url in sources]
    sha_array = _pkgbuild_array(text, "sha256sums")
    body = ("\n" + " " * 12).join(f"'{digest}'" for digest in sums)
    return text[: sha_array.start("body")] + body + text[sha_array.end("body") :]


def update_srcinfo(srcinfo: str, pkgbuild: str, version: str) -> str:
    """Rewrite the version, release, source and checksum fields of ``.SRCINFO`` from the PKGBUILD."""
    pkgrel = _PKGREL_RE.search(pkgbuild)
    if not pkgrel:
        raise RuntimeError(f"{AUR_PKGBUILD}: missing pkgrel= line")
    sources = pkgbuild_sources(pkgbuild, version)
    sums = _quoted_items(_pkgbuild_array(pkgbuild, "sha256sums").group("body"))
    if len(sources) != len(sums):
        raise RuntimeError(
            f"{AUR_PKGBUILD}: {len(sources)} sources but {len(sums)} sha256sums"
        )
    replacements = {
        "pkgver": [version],
        "pkgrel": [pkgrel.group(1)],
        "source": sources,
        "sha256sums": sums,
    }
    seen = dict.fromkeys(replacements, 0)
    out: list[str] = []
    for line in srcinfo.splitlines(keepends=True):
        field = re.match(r"^\t(\w+) = ", line)
        if field and field.group(1) in replacements:
            name = field.group(1)
            values = replacements[name]
            if seen[name] >= len(values):
                raise RuntimeError(
                    f"{AUR_SRCINFO}: more {name} lines than the PKGBUILD has"
                )
            line = f"\t{name} = {values[seen[name]]}\n"
            seen[name] += 1
        out.append(line)
    for name, values in replacements.items():
        if seen[name] != len(values):
            raise RuntimeError(
                f"{AUR_SRCINFO}: expected {len(values)} {name} line(s), found {seen[name]}"
            )
    return "".join(out)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def compute_updates(
    root: Path, version: str, checksums: ReleaseChecksums
) -> dict[str, str]:
    """Return the new content of each packaging file (changed or not)."""
    recipe = update_conda_recipe(
        root.joinpath(CONDA_RECIPE).read_text(), version, checksums
    )
    pkgbuild = update_pkgbuild(
        root.joinpath(AUR_PKGBUILD).read_text(), version, checksums
    )
    srcinfo = update_srcinfo(root.joinpath(AUR_SRCINFO).read_text(), pkgbuild, version)
    return {CONDA_RECIPE: recipe, AUR_PKGBUILD: pkgbuild, AUR_SRCINFO: srcinfo}


def run(
    root: Path,
    version_arg: str | None,
    *,
    check: bool = False,
    fetch: Fetch = http_fetch,
    cross_check: bool = True,
) -> list[str]:
    """Update (or, with ``check``, only compare) the packaging files; return the files that were or are stale."""
    root = root.resolve()
    version = (
        normalize_version(version_arg) if version_arg else read_current_version(root)
    )
    checksums = load_release_checksums(version, fetch, cross_check=cross_check)
    stale: list[str] = []
    for rel, content in compute_updates(root, version, checksums).items():
        path = root.joinpath(rel)
        if path.read_text() != content:
            stale.append(rel)
            if not check:
                path.write_text(content)
    return stale


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "version",
        nargs="?",
        help="Release version, X.Y.Z or vX.Y.Z (default: the version in build.zig)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Change nothing; exit 1 if the packaging files are stale",
    )
    parser.add_argument(
        "--no-cross-check",
        action="store_true",
        help="Trust SHA256SUMS without comparing it to the asset digests GitHub records",
    )
    args = parser.parse_args(argv)
    try:
        version = (
            normalize_version(args.version)
            if args.version
            else read_current_version(Path.cwd())
        )
        stale = run(
            Path.cwd(), version, check=args.check, cross_check=not args.no_cross_check
        )
    except Exception as exc:  # noqa: BLE001 - CLI should print concise errors.
        print(f"update-packaging-checksums: {exc}", file=sys.stderr)
        return 2
    if args.check:
        if stale:
            print(f"Packaging files are stale for v{version}:", file=sys.stderr)
            for rel in stale:
                print(f"  {rel}", file=sys.stderr)
            print(
                f"Run: scripts/update_packaging_checksums.py {version}", file=sys.stderr
            )
            return 1
        print(f"Packaging files match the published checksums of v{version}.")
        return 0
    if stale:
        print(f"Updated packaging files for v{version}:")
        for rel in stale:
            print(f"  {rel}")
    else:
        print(f"Packaging files already match the published checksums of v{version}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

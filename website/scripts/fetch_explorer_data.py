# /// script
# requires-python = ">=3.11"
# ///
"""Fetch the Parquet tables of the SASA explorer into ``website/static/data/explorer``.

The tables are built in ``zsasa-benchmarks`` (``scripts/build_explorer_tables.py``)
and archived on Zenodo; ``website/data/explorer.json`` lists each file with its
download URL, size and SHA-256. The site builder copies ``static/`` to the site,
so the explorer page reads the files from ``data/explorer/`` next to it.

    uv run website/scripts/fetch_explorer_data.py
    uv run website/scripts/fetch_explorer_data.py --from ../zsasa-benchmarks/workflows/sasa-explorer/tables

``--from`` copies the files from a local directory instead of downloading them,
and checks them against the same SHA-256 values. Files that are already in place
with the right checksum are kept.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import urllib.request
from pathlib import Path

WEB = Path(__file__).resolve().parents[1]
MANIFEST = WEB / "data" / "explorer.json"
OUT = WEB / "static" / "data" / "explorer"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--from", dest="source", type=Path, help="local directory holding the files"
    )
    args = parser.parse_args(argv)

    manifest = json.loads(MANIFEST.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    for entry in manifest["files"]:
        target = OUT / entry["name"]
        if target.exists() and sha256(target) == entry["sha256"]:
            print(f"{entry['name']}: up to date")
            continue
        partial = target.with_suffix(target.suffix + ".part")
        if args.source:
            shutil.copyfile(args.source / entry["name"], partial)
        else:
            print(
                f"{entry['name']}: downloading {entry['size'] / 1e6:.1f} MB from {entry['url']}"
            )
            with (
                urllib.request.urlopen(entry["url"]) as response,
                partial.open("wb") as f,
            ):
                shutil.copyfileobj(response, f)
        actual = sha256(partial)
        if actual != entry["sha256"]:
            partial.unlink()
            print(
                f"error: {entry['name']}: SHA-256 {actual} does not match {entry['sha256']}",
                file=sys.stderr,
            )
            return 1
        partial.replace(target)
        print(f"{entry['name']}: ok")
    return 0


if __name__ == "__main__":
    sys.exit(main())

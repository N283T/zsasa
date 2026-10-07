#!/usr/bin/env python3
"""Check that every Zig test runs in exactly one test artifact.

build.zig builds three test artifacts (the zsasa module, the CLI executable and
the C library). Zig runs the tests of every file reachable from each root, so
without filters most tests would run two or three times. build.zig therefore
filters the module and executable artifacts down to the tests no other artifact
reaches. This script fails when that partition drifts, for example after a new
file is added that only one root reaches:

* a test that exists but runs in no artifact (silently lost coverage), or
* a test that runs in more than one artifact (wasted time).

It builds the three test executables twice, with and without the filters
(`zig build test-bins`, `-Dall-tests=true`), asks each one for its test list
over the `--listen=-` protocol that the build runner uses, and compares the sets.
Only the standard library is needed. Run it from anywhere:

    python3 scripts/check_test_partition.py
"""

from __future__ import annotations

import argparse
import struct
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
ARTIFACTS = ("mod-tests", "exe-tests", "lib-tests")

# Message tags of std.zig.Client.Message and std.zig.Server.Message.
CLIENT_EXIT = 0
CLIENT_QUERY_TEST_METADATA = 4
SERVER_TEST_METADATA = 3


def build_test_bins(zig: str, prefix: Path, *, all_tests: bool) -> Path:
    cmd = [zig, "build", "test-bins", "--prefix", str(prefix)]
    if all_tests:
        cmd.append("-Dall-tests=true")
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)
    return prefix / "test-bin"


def read_exact(stream, n: int) -> bytes:
    data = stream.read(n)
    if len(data) != n:
        raise RuntimeError("test runner closed the pipe before sending test metadata")
    return data


def list_tests(binary: Path) -> list[str]:
    """Return the names of the tests compiled into a test executable."""
    proc = subprocess.Popen(
        [str(binary), "--listen=-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        cwd=REPO_ROOT,
    )
    assert proc.stdin is not None and proc.stdout is not None
    try:
        proc.stdin.write(struct.pack("<II", CLIENT_QUERY_TEST_METADATA, 0))
        proc.stdin.flush()
        while True:
            tag, length = struct.unpack("<II", read_exact(proc.stdout, 8))
            body = read_exact(proc.stdout, length)
            if tag != SERVER_TEST_METADATA:
                continue  # e.g. the zig_version greeting
            string_bytes_len, tests_len = struct.unpack_from("<II", body, 0)
            name_offsets = struct.unpack_from(f"<{tests_len}I", body, 8)
            strings_start = (
                8 + 8 * tests_len
            )  # skip the name and expected-panic index arrays
            strings = body[strings_start : strings_start + string_bytes_len]
            return [strings[o : strings.index(b"\0", o)].decode() for o in name_offsets]
    finally:
        try:
            proc.stdin.write(struct.pack("<II", CLIENT_EXIT, 0))
            proc.stdin.flush()
        except BrokenPipeError:
            pass
        proc.stdin.close()
        proc.stdout.close()
        proc.wait()


def list_artifacts(bin_dir: Path) -> dict[str, list[str]]:
    return {name: list_tests(bin_dir / name) for name in ARTIFACTS}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--zig", default="zig", help="zig executable (default: zig)")
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="zsasa-test-partition-") as tmp:
        tmp_path = Path(tmp)
        unfiltered = list_artifacts(
            build_test_bins(args.zig, tmp_path / "all", all_tests=True)
        )
        filtered = list_artifacts(
            build_test_bins(args.zig, tmp_path / "filtered", all_tests=False)
        )

    distinct = {name for names in unfiltered.values() for name in names}
    runs: dict[str, list[str]] = defaultdict(list)
    for artifact, names in filtered.items():
        for name in names:
            runs[name].append(artifact)

    print(f"{'artifact':<10} {'unfiltered':>10} {'filtered':>10}")
    for artifact in ARTIFACTS:
        print(
            f"{artifact:<10} {len(unfiltered[artifact]):>10} {len(filtered[artifact]):>10}"
        )
    executions = sum(len(names) for names in filtered.values())
    print(
        f"{'total':<10} {sum(len(n) for n in unfiltered.values()):>10} {executions:>10}"
        f"  ({len(distinct)} distinct tests)"
    )

    lost = sorted(distinct - runs.keys())
    repeated = sorted((name, arts) for name, arts in runs.items() if len(arts) > 1)
    unknown = sorted(runs.keys() - distinct)
    if not (lost or repeated or unknown):
        print("OK: every test runs in exactly one artifact")
        return 0

    print(
        "\nFAILED: the test filters in build.zig no longer partition the tests.",
        file=sys.stderr,
    )
    if lost:
        print(f"\n{len(lost)} test(s) run in no artifact:", file=sys.stderr)
        for name in lost:
            print(f"  {name}", file=sys.stderr)
    if repeated:
        print(
            f"\n{len(repeated)} test(s) run in more than one artifact:", file=sys.stderr
        )
        for name, arts in repeated:
            print(f"  {name}  [{', '.join(arts)}]", file=sys.stderr)
    if unknown:
        print(
            f"\n{len(unknown)} filtered test(s) missing from the unfiltered list:",
            file=sys.stderr,
        )
        for name in unknown:
            print(f"  {name}", file=sys.stderr)
    print(
        "\nAdjust the `.filters` of the module and executable test artifacts in build.zig so that "
        "each test runs in exactly one artifact (see CONTRIBUTING.md).",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())

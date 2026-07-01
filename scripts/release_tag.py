#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Merge a zsasa release PR when requested, then create and push the release tag."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

VERSION_RE = re.compile(r"^v?(\d+\.\d+\.\d+)$")
PASSING_CONCLUSIONS = {"SUCCESS", "SKIPPED", "NEUTRAL"}


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
        raise RuntimeError("working tree is not clean; commit/stash changes before tagging")


def load_pr(root: Path, pr: str) -> dict:
    raw = run_cmd(
        [
            "gh",
            "pr",
            "view",
            pr,
            "--json",
            "state,mergeable,isDraft,reviewDecision,statusCheckRollup,headRefName,baseRefName,url,mergedAt",
        ],
        root,
        capture=True,
    )
    return json.loads(raw)


def require_pr_checks_pass(pr_data: dict) -> None:
    if pr_data.get("isDraft"):
        raise RuntimeError("release PR is still draft")
    if pr_data.get("reviewDecision") == "CHANGES_REQUESTED":
        raise RuntimeError("release PR has requested changes")
    checks = pr_data.get("statusCheckRollup") or []
    bad: list[str] = []
    pending: list[str] = []
    for check in checks:
        name = check.get("name", "<unnamed>")
        status = check.get("status")
        conclusion = check.get("conclusion")
        if status != "COMPLETED":
            pending.append(f"{name}: {status}")
        elif conclusion not in PASSING_CONCLUSIONS:
            bad.append(f"{name}: {conclusion}")
    if pending:
        raise RuntimeError("release PR checks are still pending:\n" + "\n".join(pending))
    if bad:
        raise RuntimeError("release PR checks are not passing:\n" + "\n".join(bad))


def merge_release_pr(root: Path, pr: str, tag: str) -> None:
    pr_data = load_pr(root, pr)
    if pr_data.get("state") == "MERGED":
        return
    if pr_data.get("state") != "OPEN":
        raise RuntimeError(f"release PR #{pr} is not open or merged: {pr_data.get('state')}")
    if pr_data.get("mergeable") != "MERGEABLE":
        raise RuntimeError(f"release PR #{pr} is not mergeable: {pr_data.get('mergeable')}")
    require_pr_checks_pass(pr_data)
    run_cmd(["gh", "pr", "merge", pr, "--squash", "--delete-branch", "--subject", f"release: {tag}", "--body", ""], root)


def read_file_from_ref(root: Path, ref: str, rel: str) -> str:
    return run_cmd(["git", "show", f"{ref}:{rel}"], root, capture=True)


def verify_release_files(root: Path, ref: str, version: str, tag: str) -> None:
    required = {
        "build.zig": f'const version = "{version}";',
        "build.zig.zon": f'.version = "{version}",',
        "src/c_api.zig": f'const VERSION = "{version}";',
        "python/pyproject.toml": f'version = "{version}"',
        "CHANGELOG.md": f"## [{version}]",
        "website/docs/changelog.md": f"## [{tag}]",
    }
    missing: list[str] = []
    for rel, needle in required.items():
        if needle not in read_file_from_ref(root, ref, rel):
            missing.append(f"{rel}: missing {needle!r}")
    if missing:
        raise RuntimeError("release files do not match requested version:\n" + "\n".join(missing))


def tag_exists(root: Path, tag: str, remote: str) -> bool:
    local = subprocess.run(["git", "rev-parse", "-q", "--verify", f"refs/tags/{tag}"], cwd=root, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if local.returncode == 0:
        return True
    remote_result = subprocess.run(["git", "ls-remote", "--exit-code", "--tags", remote, f"refs/tags/{tag}"], cwd=root, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return remote_result.returncode == 0


def run(root: Path, version_arg: str, *, pr: str | None = None, merge: bool = False, remote: str = "origin") -> str:
    root = root.resolve()
    version, tag = normalize_version(version_arg)
    require_clean_tree(root)
    if pr is not None:
        pr_data = load_pr(root, pr)
        if pr_data.get("state") == "OPEN":
            if not merge:
                raise RuntimeError(f"release PR #{pr} is open; pass --merge to merge before tagging")
            merge_release_pr(root, pr, tag)
        elif pr_data.get("state") != "MERGED":
            raise RuntimeError(f"release PR #{pr} is not merged: {pr_data.get('state')}")
    run_cmd(["git", "fetch", "--prune", "--tags", remote, "main"], root)
    target_ref = "FETCH_HEAD"
    verify_release_files(root, target_ref, version, tag)
    if tag_exists(root, tag, remote):
        raise RuntimeError(f"tag {tag} already exists locally or on {remote}")
    run_cmd(["git", "tag", "-a", tag, target_ref, "-m", f"Release {tag}"], root)
    run_cmd(["git", "push", remote, tag], root)
    return tag


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="Release version, X.Y.Z or vX.Y.Z")
    parser.add_argument("--pr", help="Release PR number to verify or merge before tagging")
    parser.add_argument("--merge", action="store_true", help="Squash-merge the release PR if it is still open and checks pass")
    parser.add_argument("--remote", default="origin", help="Git remote for tag push (default: origin)")
    args = parser.parse_args(argv)
    try:
        tag = run(Path.cwd(), args.version, pr=args.pr, merge=args.merge, remote=args.remote)
    except Exception as exc:  # noqa: BLE001 - CLI should print concise errors.
        print(f"release-tag: {exc}", file=sys.stderr)
        return 1
    print(f"Pushed {tag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

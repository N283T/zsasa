---
name: zsasa-release
description: Use when preparing, merging, tagging, or troubleshooting a zsasa project release, including requests like release, version bump, changelog, create release PR, merge release PR, tag push, or publish vX.Y.Z.
---

# zsasa Release

Use this instead of the generic release flow when working in the `N283T/zsasa` repository.

## Principles

- Work from a non-`main` branch for release PR preparation.
- Keep release edits mechanical and scoped.
- Treat tag push as the publish trigger; do not tag before the release PR is merged to `main`.
- Do not bump `packaging/aur/PKGBUILD` or `.SRCINFO` in the release PR. Their checksums depend on release assets generated after tag push. `release_bump.py` leaves them alone and marks the checksums in `packaging/conda-forge/meta.yaml` as `PENDING-...`; Phase 3 fills them in.
- Use concrete dates in release metadata and changelogs.

## Phase 1: create the release PR

1. Start from clean, up-to-date `main`.
2. Create `release/vX.Y.Z`.
3. Run:

   ```bash
   ./scripts/release_bump.py X.Y.Z
   ```

   This updates the fixed zsasa version surfaces and promotes `CHANGELOG.md` / `website/docs/changelog.md` Unreleased notes into the release section.

4. Inspect for stale versions:

   ```bash
   rg -n --glob '!CHANGELOG.md' --glob '!website/docs/changelog.md' --glob '!*lock*' '\bOLD_VERSION\b' .
   ```

5. Run the release checks that match touched surfaces:

   ```bash
   zig fmt --check src/
   zig build test
   zig build -Doptimize=ReleaseFast
   ./zig-out/bin/zsasa --version
   mkdir -p /tmp/zsasa-check
   ./zig-out/bin/zsasa calc examples/1ubq.pdb /tmp/zsasa-check/output.json
   ```

   Also run `scripts/check_nix_deps_hash.py` (and `nix build` where Nix is available), the release tooling tests (`uv run --no-project --python 3.12 --with pytest python -m pytest scripts/ -q`), Python checks when Python/C ABI files changed, and website build when website docs changed.

6. Commit as `release: vX.Y.Z`, push, and open the PR.

## Phase 2: merge and tag

When the user asks to proceed after CI, do not stop at “PR is open.” Verify checks and merge/tag in one flow:

```bash
./scripts/release_tag.py X.Y.Z --pr PR_NUMBER --merge
```

The script verifies CI/check state, squash-merges the release PR when open, fetches `origin/main`, verifies release files from that commit, creates an annotated `vX.Y.Z` tag on it, and pushes the tag.

## Phase 3: packaging checksums after the release is published

Once the publish workflow has attached the binaries and `SHA256SUMS` to the GitHub release, from up-to-date `main` on a new branch:

```bash
./scripts/update_packaging_checksums.py X.Y.Z
```

This fills the conda recipe checksums and brings `packaging/aur/PKGBUILD` and `.SRCINFO` to `X.Y.Z` (the checksums are cross-checked against the asset digests GitHub records). Commit, open a PR, then submit the AUR files and update the conda-forge feedstock by hand. `./scripts/update_packaging_checksums.py X.Y.Z --check` exits non-zero while any of them is stale.

## Common mistakes

- Forgetting `python/uv.lock` or `src/c_api.zig` version bumps.
- Editing website changelog but not root `CHANGELOG.md`, or vice versa.
- Treating skipped deploy checks as failures; skipped deploy is normal on PRs.
- Tagging a release branch instead of updated `main`.
- Updating AUR `pkgver` before the asset checksum exists.
- Changing `build.zig.zon` dependencies without refreshing the Nix hash (`scripts/check_nix_deps_hash.py --refresh`).

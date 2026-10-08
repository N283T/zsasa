# AGENTS.md

Instructions for AI coding agents working on `zsasa`, a Solvent Accessible Surface Area calculator written in Zig with Python bindings. `CONTRIBUTING.md` has the developer setup and explains how the tests are organized; read it before changing tests or the build.

## Layout

- `src/` — Zig library, CLI, parsers, algorithms, C ABI (`c_api.zig`) and their tests.
- `python/` — Python package (cffi bindings over the C ABI) and its tests.
- `website/` — documentation site: Markdown in `website/docs/`, built by `website/build.py`.
- `scripts/` — release and maintenance scripts, with tests.
- `packaging/`, `flake.nix`, `Dockerfile`, `install.sh` — distribution.
- `benchmarks/`, `examples/`, `test_data/` — benchmark scripts and fixtures.

## Rules

- Do not commit to `main`. Work on a `feature/`, `fix/`, `docs/` or `release/` branch.
- Do not merge pull requests, create tags or publish anything without explicit approval from the user.
- A change in behavior comes with its tests, its documentation (`website/docs/`, `README.md`, `python/README.md`) and a `CHANGELOG.md` entry under `[Unreleased]`, mirrored in `website/docs/changelog.md`.
- Do not change the CLI, the JSON/CSV/RSA output, the C ABI or the Python API unless the task asks for it. Mark a change of default results or of an output schema as such in the changelog.
- Do not delete tracked fixtures or docs, and do not edit `website/data/benchmarks/*.json` by hand (regenerate with `website/scripts/export_benchmarks.py`). Benchmark and accuracy figures in the docs come from the paper; do not change them on your own.
- Parsers, units, atom order, classifier radii and numerical precision are high-risk: compare against existing fixtures and add a test with a tolerance.

## Checks

Run the narrowest set that covers the change, and say which checks were skipped and why.

```bash
# Zig
zig fmt --check src/
zig build test                            # prints nothing from the tests when it passes
python3 scripts/check_test_partition.py   # every test runs in exactly one artifact
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa calc examples/1ubq.pdb /tmp/zsasa-check/output.json

# Python (when python/ or the C ABI changed)
cd python && ruff format . && ruff check . && pytest tests/ -v

# Website (when website/ or the docs changed)
uv run website/build.py                   # fails on broken links and anchors

# Release tooling (scripts/, install.sh, flake.nix, build.zig.zon, packaging/)
python3 -m unittest discover -s scripts -p 'test_*.py'
python3 scripts/check_versions.py
python3 scripts/check_nix_deps_hash.py
sh -n install.sh
```

## Things that are easy to get wrong

- **Zig 0.16.** `std.Io` is passed explicitly. Release builds are ReleaseFast, where illegal behavior is silent corruption, not a panic: validate lengths and enum values read from files.
- **Python tests and the library.** A stale `python/zsasa/libzsasa.*` copy can shadow a fresh build. Set `ZSASA_LIB=$PWD/zig-out/lib/libzsasa.dylib` (or `.so`) when testing a local build, and `PYTHONPATH=<worktree>/python` in a git worktree. The integration tests need `pip install -e ".[all,dev]"`; without the optional packages they are skipped.
- **C ABI.** Keep `src/c_api.zig`, the cdef in `python/zsasa/_ffi.py` and the Python error mapping in step. Bump `ABI_VERSION` (and `_EXPECTED_ABI_VERSION`) when an existing exported signature or struct layout changes; adding an export does not need it.
- **Test partition.** A test file that only the module root or the executable root reaches needs a filter in `build.zig`, or its tests run nowhere.
- **Stderr in tests.** A test that runs a code path printing with `std.debug.print` starts with `var muted = test_support.muteStderr(); defer muted.restore();`.
- **Windows.** CI builds for Windows but runs no tests there; they first run in the publish workflow. Do not query the size of a file opened write-only (`AccessDenied` on Windows).
- **Nix.** When the dependencies in `build.zig.zon` or the Zig version change, run `scripts/check_nix_deps_hash.py --refresh` and `nix build`.
- **Zig version bump.** Also update the tarball checksums in `Dockerfile` and `python/pyproject.toml` (cibuildwheel `before-all`) and the version in `.github/workflows/`.

## Release

A pushed `vX.Y.Z` tag publishes to PyPI, GitHub Releases, GHCR, Homebrew and Scoop, and cannot be undone. Merge and tag only after the user says so.

1. From an up-to-date `main`: `git switch -c release/vX.Y.Z`, then `./scripts/release_bump.py X.Y.Z`. It bumps every version file, promotes the `[Unreleased]` notes in both changelogs, and marks the conda checksums `PENDING-...`. `packaging/aur/` stays at the previous release.
2. If defaults, output formats or accepted input changed, add an "Upgrade notes" block at the top of the new changelog section.
3. Run the checks above, plus `python3 scripts/check_versions.py --tag vX.Y.Z`. Commit as `release: vX.Y.Z`, push and open the pull request.
4. Rehearse the publish workflow on the release branch. It builds every wheel, CLI binary and the Docker image and publishes nothing:

   ```bash
   gh workflow run publish.yml --ref release/vX.Y.Z -f release_tag=vX.Y.Z -f target=none -f jobs=all
   ```

5. After approval, with CI and the rehearsal green:

   ```bash
   gh pr merge <PR> --squash --subject "release: vX.Y.Z"
   git switch main && git pull --ff-only
   python3 scripts/check_versions.py --tag vX.Y.Z
   git tag -a vX.Y.Z -m "Release vX.Y.Z" && git push origin vX.Y.Z
   ```

6. When the publish run has finished, on a new branch: `./scripts/update_packaging_checksums.py X.Y.Z`, commit and open a pull request. The AUR repository (`zsasa-bin`) and the conda-forge feedstock are updated by hand from those files.

A failed publish run can be repeated for one part with `workflow_dispatch` (`target=pypi` and the job to repeat); `target=testpypi` uploads to TestPyPI only.

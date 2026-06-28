# Review Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the prioritized whole-repository review backlog for zsasa: safety gates, parser correctness, Python reliability/streaming, documentation, CI, release, and packaging polish.

**Architecture:** Keep changes small and localized. Add regression tests before behavior changes, prefer shared validation helpers at public boundaries, and avoid schema changes unless explicitly covered by tests/docs. Split work into independent slices so subagents can operate on disjoint files and the controller can integrate and verify.

**Tech Stack:** Zig 0.16 (`zig build test`, `zig fmt`), Python 3.11+ with NumPy/CFFI/Pytest/Ruff, GitHub Actions, Docusaurus docs, shell release tooling.

---

### Task 1: Zig safety gates and build coverage

**Files:** `build.zig`, `src/format_detect.zig`, `src/c_api.zig`, `src/shrake_rupley.zig`, `src/lee_richards.zig`, `src/types.zig` if needed.

- [ ] Add a failing build/test coverage check proving `src/c_api.zig` tests are not run by `zig build test`, then include the library module test in `build.zig` so C ABI tests run.
- [ ] Add tests for uppercase compressed suffixes such as `MODEL.PDB.GZ` and `MODEL.CIF.ZST`, then make format detection case-insensitive for base and compression suffixes.
- [ ] Add C ABI/Python-facing regression tests for NaN/Inf coordinates, radii, and probe radius returning invalid-input errors rather than producing NaN output.
- [ ] Implement shared finite/range validation at C ABI and algorithm public boundaries without changing valid calculations.
- [ ] Run focused Zig tests for touched files.

### Task 2: Parser correctness

**Files:** `src/pdb_parser.zig`, `src/mmcif_parser.zig`, `src/bcif_parser.zig`, `src/calc.zig`, `src/json_writer.zig`, `src/analysis.zig`, relevant parser tests/fixtures.

- [ ] Add regression tests documenting the chosen model default semantics across PDB/mmCIF/BCIF; default should be documented and consistent with CLI docs.
- [ ] Add altLoc regression tests where a later residue/site only has conformer `B`; it must not be dropped because an earlier site chose `A`.
- [ ] Implement per-site altLoc selection, preferring blank/`A`/highest occupancy when available and otherwise keeping the only available conformer.
- [ ] Add tests for default CCD vs explicit `--classifier=ccd` HETATM behavior in `calc`, then align with batch behavior or document any intentional difference.
- [ ] Add tests for long mmCIF/BCIF chain IDs that share the first four characters; grouping/output must use full chain IDs when present.
- [ ] Run focused parser and CLI tests.

### Task 3: Python reliability and trajectory memory

**Files:** `python/zsasa/sasa.py`, `python/zsasa/classifier.py`, `python/zsasa/integrations/*.py`, `python/zsasa/mdanalysis.py`, `python/zsasa/xtc.py`, `python/zsasa/dcd.py`, `python/zsasa/_ffi.py`, `python/tests/*`.

- [ ] Add tests that `calculate_sasa` and batch APIs reject NaN/Inf coords/radii/probe and out-of-range integer parameters before CFFI calls.
- [ ] Add integration tests for unknown atom radii policy; fail with clear atom identifiers or fall back through element-derived radii where elements are available.
- [ ] Add MDAnalysis selection tests with nonzero/discontiguous residue indices, then remap selected residues densely for aggregation.
- [ ] Add chunked/streaming trajectory processing APIs for XTC/DCD/MDAnalysis paths so large trajectories can compute totals/residue aggregates without retaining all atom areas.
- [ ] Add Windows editable library discovery for `zig-out/bin`.
- [ ] Run focused pytest and Ruff checks.

### Task 4: Distribution, CI, release, docs polish

**Files:** `.github/workflows/*.yml`, `python/pyproject.toml`, `python/hatch_build.py`, `packaging/conda-forge/meta.yaml`, `CITATION.cff`, `README.md`, `website/docs/**`, `python/examples/README.md`, `examples/README.md`, `.github/ISSUE_TEMPLATE/*.yml`.

- [ ] Add PR docs/site validation so website changes are built before merge.
- [ ] Harden manual publish workflow by requiring an explicit `vX.Y.Z` version/tag input and using it for assets/releases/tags.
- [ ] Fix publish job dependencies so each manual job choice is runnable and release gates are explicit.
- [ ] Ensure PyPI sdist/source installs include required Zig source/build files or add a validated source-install path; add CI smoke coverage.
- [ ] Fix conda recipe license/source gap.
- [ ] Sync docs for SDF/MOL/TRR/AMBER NetCDF, installer behavior, Python CCD defaults, stale links, fixture counts, issue templates, and `CITATION.cff` version metadata.
- [ ] Run `npm run build` for docs when docs/site changes are present.

### Task 5: Integration verification and publishing

**Files:** whole repository.

- [ ] Run `zig fmt --check src/`.
- [ ] Run `zig build test`.
- [ ] Run `zig build -Doptimize=ReleaseFast`.
- [ ] Run CLI smoke tests: `--help`, `--version`, and example calc.
- [ ] Run focused Python checks for changed Python surface.
- [ ] Run website build if docs/site changed.
- [ ] Review full diff for secrets, generated artifacts, and unintended tracked deletions.
- [ ] Commit with conventional messages, push `feat/review-hardening`, open PR or merge only after verification and user-approved merge policy.

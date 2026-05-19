# AGENTS.md

Repository-specific instructions for AI coding agents working on `zsasa`.

## Scope

These instructions apply to the entire repository. Follow them together with the user's global Codex instructions and prefer the more specific instruction when there is a conflict.

## Project Overview

`zsasa` is a high-performance Solvent Accessible Surface Area (SASA) calculator written primarily in Zig, with Python bindings and a Docusaurus documentation site.

Key areas:

- `src/` — Zig library, CLI, parsers, algorithms, C ABI, and tests.
- `python/` — Python package and tests for the Zig-backed bindings.
- `website/` — Docusaurus documentation site.
- `docs/` and `plans/` — design docs and implementation notes.
- `benchmarks/` — benchmark and validation scripts/data.
- `examples/` and `test_data/` — small fixtures used by docs, tests, and smoke checks.

## Development Workflow

- Keep changes small, focused, and reversible.
- Do not commit directly to `main`; create a feature/fix/docs branch before committing.
- Preserve existing public CLI, JSON/CSV output, C ABI, and Python API behavior unless the task explicitly changes them.
- Do not delete tracked fixtures, generated reference data, or docs without explicit user approval.
- Avoid committing local build artifacts and caches such as `zig-out/`, `.zig-cache/`, `result`, `website/node_modules/`, and Python caches.
- If changing behavior, update relevant docs, examples, and tests in the same change.

## Build and Test Commands

Core Zig checks:

```bash
zig fmt --check src/
zig build test
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa --version
./zig-out/bin/zsasa calc examples/1ubq.pdb /tmp/zsasa-output.json
```

Python package checks, when touching `python/` or the C ABI:

```bash
cd python
pip install -e ".[dev]"
ruff format .
ruff check .
pytest tests/ -v
```

Documentation site checks, when touching `website/` or documentation build plumbing:

```bash
cd website
npm ci
npm run build
```

Use the narrowest checks that cover the changed surface area. If a check is skipped, say why.

## Zig Guidelines

- Target Zig 0.16.0+ unless the project metadata changes.
- Run `zig fmt` on touched Zig files before finalizing.
- Follow existing naming and layout conventions:
  - `snake_case` for functions and variables.
  - `PascalCase` for types.
  - Small, focused functions with explicit error handling.
- Keep allocation ownership clear and documented at API boundaries.
- Preserve scalar fallbacks when adding SIMD or platform-specific optimizations.
- Benchmark performance-sensitive changes with `-Doptimize=ReleaseFast` when practical.
- Add or update tests for parser, algorithm, CLI, and output-format changes.

## Python Guidelines

- Python lives under `python/` and targets Python 3.11+.
- Use type hints for new or changed public functions.
- Keep the Python API aligned with the C ABI in `src/c_api.zig`.
- Prefer focused tests in `python/tests/` for binding/API changes.
- Use Ruff formatting/linting for Python changes.

## Documentation Guidelines

- Update `README.md`, `CONTRIBUTING.md`, `docs/`, `website/docs/`, or `python/README.md` when changing user-visible behavior.
- Keep examples runnable and consistent with the current CLI/API.
- Do not hand-edit generated autodoc outputs; update source docs/comments instead.

## Data, Scientific Correctness, and Compatibility

- Treat parser correctness, units, atom ordering, classifier behavior, and numerical precision as high-risk areas.
- For SASA algorithm changes, compare against existing examples/fixtures and include tolerance-aware tests where possible.
- Be careful with mmCIF/PDB/SDF edge cases, compressed inputs, trajectory formats, and residue/atom metadata preservation.
- Do not silently change output schemas. If schema changes are required, update docs, examples, and downstream Python handling.

## Release and Packaging Notes

- `CHANGELOG.md`, `build.zig.zon`, `python/pyproject.toml`, packaging metadata, and install scripts may need coordinated updates for releases.
- Do not publish packages, create tags, or merge PRs without explicit user approval.
- For Nix changes, verify the flake path touched and mention any follow-up commands the user should run.

# Phase 1: README Restructure (English)

## Goal

Simplify README.md to be concise, with links to docs for details.

## Current State

- README.md: ~450 lines, very detailed
- Contains full CLI options, input/output formats, algorithm details
- Should link to docs instead of duplicating content

## Target Structure

```markdown
# freesasa-zig

One-liner description

## Features (bullet points only)

## Quick Start
- Build: `zig build -Doptimize=ReleaseFast`
- Run: `freesasa_zig structure.cif output.json`

## Installation
- Zig CLI
- Python bindings (pip install)

## Basic Usage (3-5 examples max)

## Documentation
Links to docs/ for:
- [CLI Reference](docs/cli.md)
- [Python API](docs/python.md)
- [Algorithms](docs/algorithm.md)
- [Classifiers](docs/classifier.md)
- [Benchmarks](docs/benchmark/)

## Performance (brief summary + link)

## License

## References
```

## Tasks

- [x] **DONE** - Read current README.md
- [x] Verify all CLI options against `src/main.zig`
- [x] Verify Python examples against `python/freesasa_zig/core.py`
- [x] Create new concise README.md structure
- [x] Keep important options in table
- [x] Link to existing docs (cli-io.md, python/README.md)

## Files to Read

| File | Purpose |
|------|---------|
| `src/main.zig` | CLI options, help text |
| `python/freesasa_zig/core.py` | Python API |
| `python/freesasa_zig/__init__.py` | Public exports |

## Notes

- Keep all content in English
- Japanese version created in separate phase

---
- [x] **DONE** - Phase complete

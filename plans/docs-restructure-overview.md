# Documentation Restructure Plan

## Overview

Restructure and verify all documentation for v0.1.0 release.

## Principles

1. **README is concise** - Links to docs for details
2. **English first** - Create English docs, then translate
3. **Verify against code** - All examples must work
4. **Phase by phase** - Small, focused changes

## Phases

| Phase | Focus | Status |
|-------|-------|--------|
| [Phase 1](docs-restructure-phase1-readme.md) | README restructure | **Done** |
| [Phase 2](docs-restructure-phase2-cli-docs.md) | CLI documentation | **Done** |
| [Phase 3](docs-restructure-phase3-python-docs.md) | Python documentation | **Done** |
| [Phase 4](docs-restructure-phase4-existing-docs.md) | Review existing docs | Pending |
| [Phase 5](docs-restructure-phase5-changelog.md) | CHANGELOG review | Pending |
| [Phase 6](docs-restructure-phase6-japanese.md) | Japanese documentation | Pending |

## Dependencies

```
Phase 1 (README)
    ↓
Phase 2 (CLI docs) ←──┐
    ↓                 │
Phase 3 (Python docs) │
    ↓                 │
Phase 4 (Existing docs)
    ↓
Phase 5 (CHANGELOG)
    ↓
Phase 6 (Japanese)
```

## Key Files

### Source Code (for verification)

| File | Documents |
|------|-----------|
| `src/main.zig` | CLI options |
| `src/analysis.zig` | Analysis features |
| `src/classifier*.zig` | Classifiers |
| `python/freesasa_zig/core.py` | Python API |
| `python/freesasa_zig/*.py` | Integrations |

### Documentation Output

| File | Content |
|------|---------|
| `README.md` | Concise overview |
| `docs/cli.md` | CLI reference (NEW) |
| `docs/python.md` | Python API (NEW) |
| `docs/algorithm.md` | Algorithms (translate) |
| `docs/classifier.md` | Classifiers (translate) |
| `docs/benchmark/` | Benchmarks (translate) |

## Success Criteria

- [ ] README < 200 lines
- [ ] All CLI options documented with working examples
- [ ] All Python API documented with working examples
- [ ] All docs in English
- [ ] Japanese translations for key docs
- [ ] CHANGELOG complete and accurate
- [ ] No broken links

---
- [ ] **DONE** - All phases complete

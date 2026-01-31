# Phase 4: Review Existing Docs (English)

## Goal

Review and update existing docs to be accurate and in English.

## Depends On

- Phases 1-3 complete

## Current Docs (Japanese)

**Note**: README.md currently links to these Japanese docs. After translation, English readers can use the links properly.

| File | Content | Action | README Link |
|------|---------|--------|-------------|
| `docs/algorithm.md` | SASA algorithms | Translate, verify | Yes |
| `docs/architecture.md` | Code architecture | Translate, update | No |
| `docs/classifier.md` | Atom classifiers | Translate, verify | Yes |
| `docs/optimizations.md` | Optimization techniques | Translate, verify | Yes |
| `docs/cli-io.md` | CLI and I/O | Translate (or replace) | Yes |
| `docs/cpu-efficiency.md` | CPU analysis | Translate | No |
| `docs/ci.md` | CI/CD config | Translate, update | No |
| `docs/benchmark/` | Benchmark docs | Translate both | Yes |

## Tasks

For each doc:
- [ ] Read corresponding source code
- [ ] Translate to English
- [ ] Verify technical accuracy
- [ ] Update outdated sections
- [ ] Fix broken links

## Priority Order

1. `algorithm.md` - Core technical content
2. `classifier.md` - Important for users
3. `benchmark/methodology.md` - Reproducibility
4. `benchmark/results.md` - Performance claims
5. `architecture.md` - Developer reference
6. `optimizations.md` - Technical details
7. `cpu-efficiency.md` - Advanced analysis
8. `ci.md` - Infrastructure

## Output

- All docs in `docs/` translated to English
- `docs/cli-io.md` removed (replaced by `docs/cli.md`)

---
- [ ] **DONE** - Phase complete

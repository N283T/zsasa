# Phase 4: Review Existing Docs (English)

## Goal

Review and update existing docs to be accurate and in English.

## Depends On

- Phases 1-3 complete

## Current Docs (Japanese)

| File | Content | Action |
|------|---------|--------|
| `docs/algorithm.md` | SASA algorithms | Translate to English, verify accuracy |
| `docs/architecture.md` | Code architecture | Translate, update for current structure |
| `docs/classifier.md` | Atom classifiers | Translate, verify against code |
| `docs/optimizations.md` | Optimization techniques | Translate, verify |
| `docs/cli-io.md` | CLI and I/O | **Replace with cli.md** (Phase 2) |
| `docs/cpu-efficiency.md` | CPU analysis | Translate |
| `docs/ci.md` | CI/CD config | Translate, update for current CI |
| `docs/benchmark/` | Benchmark docs | Translate both files |

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

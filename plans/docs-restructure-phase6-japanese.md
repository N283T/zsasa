# Phase 6: Japanese Documentation

## Goal

Create Japanese versions of all documentation.

## Depends On

- Phases 1-5 complete (English docs finalized)

## Tasks

### README

- [ ] Create/Update `README.ja.md` from new `README.md`
- [ ] Ensure all links point to Japanese docs where available

### Docs

Create Japanese versions in `docs/ja/` or as `*.ja.md`:

- [ ] `docs/cli.ja.md` or `docs/ja/cli.md`
- [ ] `docs/python.ja.md` or `docs/ja/python.md`
- [ ] `docs/algorithm.ja.md` (may already exist)
- [ ] `docs/classifier.ja.md`
- [ ] `docs/benchmark/methodology.ja.md`
- [ ] `docs/benchmark/results.ja.md`
- [ ] Other docs as needed

## Directory Structure Options

### Option A: Suffix style
```
docs/
├── cli.md
├── cli.ja.md
├── python.md
├── python.ja.md
└── ...
```

### Option B: Directory style
```
docs/
├── en/
│   ├── cli.md
│   └── python.md
├── ja/
│   ├── cli.md
│   └── python.md
└── README.md (index)
```

## Decision

TBD - choose during implementation

## Notes

- Keep technical terms consistent (SASA, RSA, etc.)
- Use standard Japanese computing terminology
- Reference existing Japanese docs for style

---
- [ ] **DONE** - Phase complete

# Phase 2: CLI Documentation (English)

## Goal

Create comprehensive CLI documentation in docs/cli.md.

## Depends On

- Phase 1 (README restructure)

## Tasks

- [x] Read `src/main.zig` for all CLI options
- [x] Read `src/analysis.zig` for analysis features
- [x] Read `src/mmcif_parser.zig` for input handling
- [x] Create `docs/cli.md` with:
  - Full options table
  - Input formats (JSON, mmCIF)
  - Output formats (JSON, compact, CSV)
  - Examples for each feature
  - Error messages reference
- [x] Verify examples actually work

## Files to Read

| File | Purpose |
|------|---------|
| `src/main.zig` | All CLI options, help text, validation |
| `src/json_parser.zig` | JSON input format |
| `src/mmcif_parser.zig` | mmCIF input handling |
| `src/json_writer.zig` | Output formats |
| `src/analysis.zig` | --per-residue, --rsa, --polar |

## Output

- `docs/cli.md` - Comprehensive CLI reference

---
- [x] **DONE** - Phase complete

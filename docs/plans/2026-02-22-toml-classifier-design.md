# Custom Classifier from TOML Config - Design

**Issue:** #158
**Date:** 2026-02-22

## Goal

Allow users to define custom atom classifiers in TOML format, auto-detected by `.toml` file extension alongside the existing FreeSASA config format.

## TOML Schema

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }
C_CAR = { radius = 1.76, class = "apolar" }
N     = { radius = 1.65, class = "polar" }
O     = { radius = 1.40, class = "polar" }
S     = { radius = 1.85, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "N"
type = "N"

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"

[[atoms]]
residue = "ALA"
atom = "CB"
type = "C_ALI"
```

**Rationale:**
- `[types]` table: each key is a type name, value is inline table with radius + class
- `[[atoms]]` array of tables: natural for list of (residue, atom, type) mappings
- Maps 1:1 to existing FreeSASA format semantics
- Human-readable with TOML syntax highlighting in editors

## Architecture

```
--config=FILE
      |
      +--.toml extension? --> toml_classifier_parser.parseConfig()
      |                              |
      +--(else)           --> classifier_parser.parseConfig()  (existing)
                                     |
                          Both produce --> classifier.Classifier
```

### New Files

- `src/toml_parser.zig` - Minimal TOML subset parser (generic, reusable)
- `src/toml_classifier_parser.zig` - TOML-to-Classifier converter

### Modified Files

- `src/classifier_parser.zig` - Dispatch on `.toml` extension in `parseConfigFile`
- `docs/cli.md` - Document TOML format
- `CHANGELOG.md` - Add entry

### Unchanged

- `src/main.zig` - Already calls `parseConfigFile`, no changes needed

## TOML Subset Parser

Handles only what classifier configs require:

| Feature | Supported |
|---------|-----------|
| `key = "string"` | Yes |
| `key = float/int` | Yes |
| `[table]` | Yes |
| `{ k = v, ... }` inline tables | Yes |
| `[[array_of_tables]]` | Yes |
| `# comments` | Yes |
| Multiline strings | No |
| Literal strings (`'...'`) | No |
| Datetime | No |
| Nested tables beyond 1 level | No |
| Dotted keys | No |

Estimated size: ~200-400 lines.

## Format Detection

Auto-detect by file extension:
- `.toml` -> TOML parser
- Anything else (`.config`, `.cfg`, no extension) -> existing FreeSASA parser

No new CLI flags needed.

## Error Handling

Same error quality as existing FreeSASA parser:
- `UndefinedType` - atoms reference a type not in `[types]`
- `DuplicateType` - same type name defined twice
- `InvalidRadius` / `InvalidClass` - bad values in type definitions
- `MissingField` - atom entry missing `residue`, `atom`, or `type`
- Line numbers in error messages where possible

## Testing

- TOML parser unit tests (strings, floats, tables, array of tables, comments, edge cases)
- TOML classifier parser tests (valid config, error cases, ANY fallback)
- Round-trip test: equivalent FreeSASA and TOML configs produce identical `Classifier`
- Integration: `--config=test.toml` works end-to-end

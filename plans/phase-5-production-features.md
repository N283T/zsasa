# Phase 5: Production Features

## Goal

Add production-ready CLI features for better usability and configuration.

**Current**: Basic CLI with `--threads` option only
**Target**: Full-featured CLI with configurable parameters and multiple output formats

## Sub-Phases

### Phase 5.1: Extended CLI Options

**Goal**: Add commonly-used algorithm parameters as CLI options.

**Tasks**:
- [x] Add `--probe-radius=<float>` option (default: 1.4 Å)
- [x] Add `--n-points=<int>` option (default: 100)
- [x] Add `--quiet` flag (suppress progress output)
- [x] Add `--help` / `-h` flag
- [x] Add `--version` / `-V` flag
- [x] Update usage help text
- [x] Unit tests for argument parsing

**CLI Design**:
```bash
freesasa_zig [OPTIONS] <input.json> [output.json]

Options:
  --threads=N        Number of threads (default: auto-detect)
  --probe-radius=R   Probe radius in Å (default: 1.4)
  --n-points=N       Test points per atom (default: 100)
  --quiet            Suppress progress output
  -h, --help         Show help message
  -V, --version      Show version
```

**Files**:
| File | Action |
|------|--------|
| `src/main.zig` | MODIFY |
| `build.zig` | MODIFY (add version) |

---
- [x] **DONE** - Sub-phase 5.1 complete

---

### Phase 5.2: Output Format Options

**Goal**: Support multiple output formats for different use cases.

**Tasks**:
- [x] Add `--format=<json|csv|compact>` option
- [x] Implement CSV output (atom_index, area)
- [x] Implement compact JSON (no pretty-print)
- [x] Implement pretty-printed JSON (default)
- [x] Update json_writer.zig for format options

**Output Formats**:

**JSON (default)**:
```json
{
  "total_area": 1234.56,
  "atom_areas": [10.5, 20.3, ...]
}
```

**CSV**:
```csv
atom_index,area
0,10.5
1,20.3
...
total,1234.56
```

**Compact JSON**:
```json
{"total_area":1234.56,"atom_areas":[10.5,20.3,...]}
```

**Files**:
| File | Action |
|------|--------|
| `src/main.zig` | MODIFY |
| `src/json_writer.zig` | MODIFY |

---
- [x] **DONE** - Sub-phase 5.2 complete

---

### Phase 5.3: Input Validation & Error Messages

**Goal**: Provide clear, actionable error messages.

**Tasks**:
- [ ] Validate input JSON schema (required fields)
- [ ] Check array length consistency (x, y, z, r same length)
- [ ] Validate radius values (positive, reasonable range)
- [ ] Validate coordinates (finite, reasonable range)
- [ ] Provide line/column numbers for JSON parse errors
- [ ] Add `--validate` flag for dry-run validation

**Error Message Examples**:
```
Error: Input validation failed
  - 'x' array has 100 elements, but 'r' has 99
  - Fix: Ensure all coordinate arrays have the same length

Error: Invalid radius at index 42
  - Value: -1.5 (must be positive)
  - Fix: Check van der Waals radii data

Error: JSON parse error at line 15, column 8
  - Expected ',' or ']', found 'x'
```

**Files**:
| File | Action |
|------|--------|
| `src/json_parser.zig` | MODIFY |
| `src/main.zig` | MODIFY |

---
- [ ] **DONE** - Sub-phase 5.3 complete

---

### Phase 5.4: Documentation & Polish

**Goal**: Complete documentation and final polish.

**Tasks**:
- [ ] Update README with all CLI options
- [ ] Add examples for common use cases
- [ ] Add CHANGELOG.md
- [ ] Add man page or detailed --help output
- [ ] Final benchmark with all optimizations

**Files**:
| File | Action |
|------|--------|
| `README.md` | MODIFY |
| `CHANGELOG.md` | CREATE |

---
- [ ] **DONE** - Sub-phase 5.4 complete

---

## Architecture

```
CLI Input
    │
    ▼
┌─────────────────┐
│ Argument Parser │ ─── --help, --version (early exit)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Input Validator │ ─── --validate (dry-run exit)
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ SASA Calculator │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│ Output Writer   │ ─── json/csv/compact
└─────────────────┘
```

## Success Criteria

- [ ] All existing tests pass
- [ ] New CLI options work correctly
- [ ] Error messages are clear and actionable
- [ ] Documentation is complete and accurate
- [ ] Backward compatible (existing usage still works)

## Verification Commands

```bash
# Build and test
zig build test

# Test new options
./zig-out/bin/freesasa_zig --help
./zig-out/bin/freesasa_zig --version
./zig-out/bin/freesasa_zig --probe-radius=1.5 --n-points=200 examples/input_1a0q.json /tmp/out.json
./zig-out/bin/freesasa_zig --format=csv examples/input_1a0q.json /tmp/out.csv
./zig-out/bin/freesasa_zig --quiet examples/input_1a0q.json /tmp/out.json

# Validate input
./zig-out/bin/freesasa_zig --validate examples/input_1a0q.json
```

## Priority

Recommended implementation order:
1. **5.1** (CLI Options) - Most commonly requested
2. **5.3** (Validation) - Improves user experience
3. **5.2** (Output Formats) - Nice to have
4. **5.4** (Documentation) - Final polish

---
- [ ] **DONE** - Phase 5 complete

# Roadmap: FreeSASA Complete Backward Compatibility

## Goal

既存のFreeSASAユーザーがすぐにfreesasa-zigに乗り換えられるようにする。

## Current State

freesasa-zigは独自のCLI設計:
```bash
freesasa-zig input.cif output.json --algorithm=sr --n-points=100
```

FreeSASA Cは:
```bash
freesasa input.pdb --format=json --algorithm=sr -n 100
```

## Requirements

### Phase 1: CLI Compatibility

| FreeSASA Option | Current | Action |
|-----------------|---------|--------|
| `--format=json/csv/pdb` | JSON only | Add CSV, PDB output |
| `--algorithm=sr/lr` | ✅ Supported | - |
| `-n, --resolution` | `--n-points` | Add alias |
| `--probe-radius` | ✅ Supported | - |
| `--chain-groups` | ❌ | Add |
| `--separate-chains` | ❌ | Add |
| `--separate-models` | ❌ | Add |
| `--select` | ❌ | Add (selection syntax) |
| `--config-file` | `--classifier` | Add file path support |
| `--radii` | `--classifier` | Add alias |
| `--no-warnings` | ❌ | Add |
| `--no-log` | ❌ | Add |

### Phase 2: Output Format Compatibility

**JSON Output**
```json
// FreeSASA format
{
  "structure": {
    "name": "1crn",
    "atoms": [...],
    "residues": [...],
    "chains": [...]
  }
}

// freesasa-zig current
{
  "total_area": 1234.56,
  "atom_areas": [...]
}
```

→ `--format=freesasa-json` for FreeSASA-compatible output?

**CSV Output**
- Per-atom CSV
- Per-residue CSV
- Summary CSV

### Phase 3: Selection Syntax

FreeSASA selection syntax:
```
chain A
resn ALA
resi 1-10
name CA
```

→ Implement basic selection parser

### Phase 4: Configuration File

FreeSASA config file format:
```
types:
ALA C_ALI 2.00
...
atoms:
ALA CA C_ALI
...
```

→ Already supported via `--classifier` (need to verify format)

## Implementation Order

1. [ ] CLI aliases (`-n`, `--radii`)
2. [ ] Output formats (CSV, PDB)
3. [ ] Chain/model separation options
4. [ ] Selection syntax (basic)
5. [ ] FreeSASA-compatible JSON output
6. [ ] Documentation (migration guide)

## Non-Goals

- 100% identical output (floating-point differences OK)
- GUI support
- FreeSASA C library API compatibility

## Success Criteria

```bash
# Should work for most common use cases
freesasa-zig protein.pdb --format=json --algorithm=sr -n 100
freesasa-zig protein.cif --select "chain A" --format=csv
```

---
- [ ] **DONE** - Phase complete

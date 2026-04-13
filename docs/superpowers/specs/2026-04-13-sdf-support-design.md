# SDF File Support Design

## Overview

Add SDF (Structure-Data File) support to zsasa with two use cases:

1. **Direct input**: Calculate SASA for molecules in SDF files
2. **Supplementary bond info**: Provide bond topology for CCD-unregistered compounds via `--sdf` option

## Motivation

- SDF files contain bond information, enabling CCD-style hybridization-based radii derivation
- Compounds predicted by tools like Boltz lack CCD entries; `--sdf` eliminates the need to create custom definition files
- SDF is the most common small-molecule exchange format

## Use Cases

### Use Case 1: SDF as Input File

```bash
zsasa calc molecule.sdf                      # default: CCD classifier from bond info
zsasa calc molecule.sdf --classifier naccess  # explicit classifier
```

- Parse SDF, extract coordinates and bond topology
- Each molecule = 1 chain (A, B, C...), residue number = 1, residue name = molecule name (truncated to 5 chars)
- Atom names generated as element + index (C1, C2, O1...)
- Default classifier = CCD; bond info → hybridization analysis → radii
- All molecules in multi-molecule SDF are processed individually

### Use Case 2: `--sdf` Option for CCD Classifier

```bash
zsasa calc structure.pdb --classifier ccd --sdf ligand.sdf
zsasa calc structure.cif --sdf ligand1.sdf --sdf ligand2.sdf
```

- `--sdf` can be specified multiple times
- Each SDF molecule → `Component` (molecule name = comp_id)
- Registered via existing `addComponent` into CCD classifier
- Compatible with `--ccd` (both can be used together)
- SDF molecule name matched against PDB/mmCIF residue names

## Architecture

### New Module: `sdf_parser.zig`

Core SDF/MOL parser supporting V2000 and V3000 formats.

#### Data Structures

```zig
pub const SdfMolecule = struct {
    name: []const u8,           // header line 1 (trimmed)
    atoms: []const SdfAtom,
    bonds: []const SdfBond,
};

pub const SdfAtom = struct {
    x: f64,
    y: f64,
    z: f64,
    element: element.Element,
};

pub const SdfBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: hybridization.BondOrder,
};
```

#### Public API

```zig
/// Parse SDF/MOL file, return all molecules
pub fn parse(allocator: Allocator, source: []const u8) ![]const SdfMolecule;

/// Convert SdfMolecule to Component for CCD classifier
pub fn toComponent(molecule: SdfMolecule, allocator: Allocator) !hybridization.Component;

/// Convert multiple molecules to AtomInput (use case 1)
pub fn toAtomInput(molecules: []const SdfMolecule, allocator: Allocator) !types.AtomInput;
```

#### Format Detection

- Counts line (line 4) ending with `V2000` or `V3000` → auto-detect
- No version marker → default to V2000
- `$$$$` absent at EOF → treat as single MOL file

#### V2000 Parsing

```
MoleculeName          ← header line 1 = molecule name
  program info        ← line 2 (ignored)
Comment               ← line 3 (ignored)
 20 21  0  0  0  0  0  0  0  0999 V2000  ← counts line
    1.0000    2.0000    3.0000 C   0  0   ← atom block (x, y, z, element)
  1  2  1  0  0  0  0                     ← bond block (idx1, idx2, order)
M  END
$$$$
```

#### V3000 Parsing

```
MoleculeName
  program info
Comment
  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 20 21 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.0000 2.0000 3.0000 0        ← idx elem x y z charge
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2                             ← idx order atom1 atom2
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
```

#### Bond Order Mapping

| SDF value | `hybridization.BondOrder` |
|-----------|--------------------------|
| 1         | `.single`                |
| 2         | `.double`                |
| 3         | `.triple`                |
| 4         | `.aromatic`              |

### format_detect.zig Changes

```zig
// Add to supported_extensions
".sdf", ".sdf.gz",
".mol", ".mol.gz",

// Add to InputFormat enum
pub const InputFormat = enum { json, mmcif, pdb, sdf };
```

### calc.zig Changes

#### New CLI Option

```zig
// Add to CalcArgs
sdf_paths: []const []const u8 = &.{},  // --sdf=PATH (multiple allowed)
```

Add `--sdf` to `batch` and `traj` subcommands as well.

#### Use Case 1 Integration

- Detect SDF format via extension
- Parse with `sdf_parser.parse()`
- Convert to `AtomInput` via `sdf_parser.toAtomInput()`
- Default classifier = CCD
- For each molecule, create `Component` from bond info and `addComponent`

#### Use Case 2 Integration

- Parse each `--sdf` file
- Convert each molecule to `Component` via `sdf_parser.toComponent()`
- Register in CCD classifier alongside inline/external CCD data

#### CCD Classifier Lookup Priority

For non-standard residues:

1. Hardcoded table (standard amino acids)
2. SDF-derived Components (`--sdf`)
3. Inline CCD (from mmCIF file)
4. External CCD dictionary (`--ccd`)

SDF takes priority because user-specified definitions should override automatic sources.

## Error Handling

### SDF Parse Errors

| Case | Behavior |
|------|----------|
| Empty file | Error exit: `Error: SDF file contains no molecules` |
| Molecule with 0 atoms | Warning + skip to next molecule |
| Invalid coordinate | Error exit: `Error: Invalid coordinate at line N` |
| Out-of-range bond index | Error exit: `Error: Bond references non-existent atom at line N` |
| Missing `M  V30` prefix in V3000 | Error exit: `Error: Expected V3000 block at line N` |
| EOF without `$$$$` | Treat as final molecule (MOL compatibility) |

### Use Case 2 Matching Errors

| Case | Behavior |
|------|----------|
| SDF molecule name not in input structure | Warning: `Warning: SDF molecule 'XXX' not found in input structure` |
| Empty molecule name | Warning + skip |
| Name > 5 chars | Truncate to 5 chars + warning |
| Duplicate molecule names in SDF | Later definition overwrites + warning |

### Edge Cases

- **Chain ID limit (use case 1)**: 26 molecules max (A-Z). Molecules beyond 26 are skipped with warning.
- **`--include-hydrogens`**: Hydrogen atoms in SDF are skipped by default (consistent with PDB/mmCIF). Included when `--include-hydrogens` is specified.

## Test Strategy

### Unit Tests (sdf_parser.zig)

- **V2000 parse**: Minimal molecule (methane: 5 atoms, 4 bonds) — verify coordinates, elements, bond orders
- **V3000 parse**: Same molecule in V3000 — verify identical results
- **Multi-molecule**: 2-molecule SDF — verify `$$$$` splitting
- **MOL compatibility**: Single molecule without `$$$$` — verify correct parse
- **Bond order mapping**: All bond types (1-4) mapped correctly
- **Error cases**: Empty file, invalid coordinates, out-of-range bond indices

### Unit Tests (toComponent, toAtomInput)

- **toComponent**: Verify atom/bond count, comp_id correctness
- **toAtomInput**: Verify chain ID assignment (A, B, C...), residue names, coordinate copy

### Integration Tests

- **Use case 1**: `zsasa calc test.sdf` — SASA values in reasonable range
- **Use case 1 + classifier**: `zsasa calc test.sdf --classifier naccess` — runs without error
- **Use case 2**: `zsasa calc test.pdb --sdf ligand.sdf --classifier ccd` — ligand radii correctly set
- **Multiple `--sdf`**: Two SDF files specified simultaneously — both molecules registered

### Test Data

- Add small SDF files to `test_data/` (ethanol or similar, both V2000 and V3000 variants)

## Files to Create/Modify

| File | Action |
|------|--------|
| `src/sdf_parser.zig` | **Create** — SDF/MOL parser |
| `src/format_detect.zig` | **Modify** — add `.sdf`, `.mol` extensions and `sdf` format |
| `src/calc.zig` | **Modify** — add `--sdf` option, SDF input handling |
| `src/batch.zig` | **Modify** — add `--sdf` option |
| `src/traj.zig` | **Modify** — add `--sdf` option |
| `src/main.zig` | **Modify** — wire SDF format through CLI |
| `build.zig` | **Modify** — add `sdf_parser.zig` module |
| `test_data/*.sdf` | **Create** — test SDF files |

# Parsers

Three parsers convert structure files into `AtomInput` for SASA calculation.

PDB and mmCIF parsers:
- Return `AtomInput` with coordinates, radii, and metadata
- Assign default van der Waals radii based on element type during parsing
- Populate optional metadata fields (residue names, chain IDs, etc.)

The JSON parser expects radii to be provided in the input `r` array.

---

## PDB Parser

Import via `zsasa.pdb_parser`.

Parses fixed-width PDB format files (`.pdb`, `.ent`).

### PdbParser

```zig
pub const PdbParser = struct {
    allocator: Allocator,
    atom_only: bool = true,
    skip_hydrogens: bool = true,
    first_alt_loc_only: bool = true,
    model_num: ?u32 = null,
    chain_filter: ?[]const []const u8 = null,
};
```

**Configuration fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `atom_only` | `bool` | `true` | Include only ATOM records (exclude HETATM) |
| `skip_hydrogens` | `bool` | `true` | Exclude hydrogen atoms (H, D) |
| `first_alt_loc_only` | `bool` | `true` | Keep only first alternate location |
| `model_num` | `?u32` | `null` | Model number to extract (`null` = first model) |
| `chain_filter` | `?[]const []const u8` | `null` | Chain IDs to include (`null` = all chains) |

### init

```zig
pub fn init(allocator: Allocator) PdbParser
```

Creates a parser with default settings.

### parse

```zig
pub fn parse(self: *PdbParser, source: []const u8) !AtomInput
```

Parses PDB data from a string buffer.

### parseFile

```zig
pub fn parseFile(self: *PdbParser, path: []const u8) !AtomInput
```

Reads and parses a PDB file from disk.

### Errors

| Error | Description |
|-------|-------------|
| `InvalidCoordinate` | Non-numeric coordinate value |
| `FileReadError` | Cannot open or read file |
| `NoAtomsFound` | No atoms match the filter criteria |
| `LineTooShort` | PDB line shorter than required fields |

### Example

```zig
const zsasa = @import("zsasa");

var parser = zsasa.pdb_parser.PdbParser.init(allocator);

// Parse from file
var atoms = try parser.parseFile("protein.pdb");
defer atoms.deinit();

// Parse with custom settings
parser.skip_hydrogens = false;    // Include hydrogens
parser.atom_only = false;         // Include HETATM
parser.model_num = 2;             // Extract model 2

var atoms2 = try parser.parseFile("multi_model.pdb");
defer atoms2.deinit();
```

---

## mmCIF Parser

Import via `zsasa.mmcif_parser`.

Parses mmCIF/PDBx format files (`.cif`). Extracts data from `_atom_site` loop.

### MmcifParser

```zig
pub const MmcifParser = struct {
    allocator: Allocator,
    atom_only: bool = true,
    skip_hydrogens: bool = true,
    first_alt_loc_only: bool = true,
    model_num: ?u32 = null,
    chain_filter: ?[]const []const u8 = null,
    use_auth_chain: bool = false,
};
```

**Configuration fields:**

Same as `PdbParser`, plus:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `use_auth_chain` | `bool` | `false` | Use `auth_asym_id` instead of `label_asym_id` for chain IDs |

### init

```zig
pub fn init(allocator: Allocator) MmcifParser
```

### parse

```zig
pub fn parse(self: *MmcifParser, source: []const u8) !AtomInput
```

### parseFile

```zig
pub fn parseFile(self: *MmcifParser, path: []const u8) !AtomInput
```

### Errors

| Error | Description |
|-------|-------------|
| `NoAtomSiteLoop` | No `_atom_site` loop found |
| `MissingCoordinateField` | Missing `Cartn_x`, `Cartn_y`, or `Cartn_z` column |
| `InvalidCoordinate` | Non-numeric coordinate value |
| `FileReadError` | Cannot open or read file |

### Example

```zig
const zsasa = @import("zsasa");

var parser = zsasa.mmcif_parser.MmcifParser.init(allocator);
var atoms = try parser.parseFile("structure.cif");
defer atoms.deinit();

// Use auth chain IDs (PDB-style)
parser.use_auth_chain = true;
var atoms2 = try parser.parseFile("structure.cif");
defer atoms2.deinit();
```

---

## JSON Parser

Import via `zsasa.json_parser`.

Parses zsasa's JSON input format. Supports `.json` and `.json.gz` (gzip-compressed).

### parseAtomInput

```zig
pub fn parseAtomInput(
    allocator: Allocator,
    json_str: []const u8,
) !AtomInput
```

Parses JSON from a string buffer.

**Expected JSON format:**

```json
{
    "x": [1.0, 2.0, 3.0],
    "y": [4.0, 5.0, 6.0],
    "z": [7.0, 8.0, 9.0],
    "r": [1.5, 1.7, 1.5],
    "residue": ["ALA", "ALA", "GLY"],
    "atom_name": ["CA", "CB", "CA"]
}
```

The `x`, `y`, `z`, `r` fields are required. `residue`, `atom_name`, and `element` are optional.

### readAtomInputFromFile

```zig
pub fn readAtomInputFromFile(
    allocator: Allocator,
    path: []const u8,
) !AtomInput
```

Reads and parses a JSON file. Automatically decompresses `.json.gz` files.

### validateInput

```zig
pub fn validateInput(
    allocator: Allocator,
    input: AtomInput,
) !ValidationResult
```

Validates input data: checks for NaN/Inf coordinates, positive radii, and radius range (max 100 A).

**Returns:** `ValidationResult` with `valid: bool` and `errors: []ValidationError`.

```zig
var validation = try zsasa.json_parser.validateInput(allocator, atoms);
defer validation.deinit();

if (!validation.valid) {
    for (validation.errors) |err| {
        std.debug.print("Error: {s}\n", .{err.message});
    }
}
```

### checkDuplicateCoordinates

```zig
pub fn checkDuplicateCoordinates(
    allocator: Allocator,
    input: AtomInput,
) !usize
```

Returns the number of duplicate coordinate pairs found. Useful for detecting input issues.

### Example

```zig
const zsasa = @import("zsasa");

// From file (auto-detects .gz)
var atoms = try zsasa.json_parser.readAtomInputFromFile(allocator, "input.json.gz");
defer atoms.deinit();

// From string
const json = @embedFile("test_input.json");
var atoms2 = try zsasa.json_parser.parseAtomInput(allocator, json);
defer atoms2.deinit();
```

# Types

Core data types for SASA calculation. Import via `zsasa.types`.

```zig
const zsasa = @import("zsasa");
const types = zsasa.types;
```

---

## AtomInput

Input data for SASA calculation. Returned by all parsers.

```zig
pub const AtomInput = struct {
    x: []const f64,          // X coordinates
    y: []const f64,          // Y coordinates
    z: []const f64,          // Z coordinates
    r: []f64,                // Atomic radii (mutable for classifier)

    // Optional metadata (populated by parsers)
    residue: ?[]const FixedString5,     // Residue names ("ALA", "GLY")
    atom_name: ?[]const FixedString4,   // Atom names ("CA", "CB")
    element: ?[]const u8,               // Element atomic numbers
    chain_id: ?[]const FixedString4,    // Chain IDs ("A", "B")
    residue_num: ?[]const i32,          // Residue sequence numbers
    insertion_code: ?[]const FixedString4, // Insertion codes

    allocator: std.mem.Allocator,
};
```

All coordinate and radius arrays must have the same length (one entry per atom). Optional metadata arrays, when present, must also have the same length.

### Methods

#### atomCount

```zig
pub fn atomCount(self: AtomInput) usize
```

Returns the number of atoms (length of coordinate arrays).

#### hasClassificationInfo

```zig
pub fn hasClassificationInfo(self: AtomInput) bool
```

Returns `true` if both `residue` and `atom_name` are present. Required for classifier lookups.

#### hasElementInfo

```zig
pub fn hasElementInfo(self: AtomInput) bool
```

Returns `true` if `element` (atomic numbers) is present. Used for radius guessing fallback.

#### hasResidueInfo

```zig
pub fn hasResidueInfo(self: AtomInput) bool
```

Returns `true` if all per-residue fields are present (`residue`, `residue_num`, `chain_id`, `insertion_code`). Required for `analysis.aggregateByResidue()`.

#### deinit

```zig
pub fn deinit(self: *AtomInput) void
```

Frees all allocated memory. Always call via `defer`:

```zig
var atoms = try parser.parseFile("protein.pdb");
defer atoms.deinit();
```

---

## Config / Configf32

Configuration for the Shrake-Rupley algorithm.

```zig
pub const Config = ConfigGen(f64);
pub const Configf32 = ConfigGen(f32);
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `n_points` | `u32` | `100` | Test points per atom (higher = more accurate, slower) |
| `probe_radius` | `T` | `1.4` | Water probe radius in Angstroms |

**Example:**

```zig
// Default config
const config = zsasa.types.Config{};

// Custom config
const config = zsasa.types.Config{
    .n_points = 200,
    .probe_radius = 1.8,
};
```

---

## SasaResult / SasaResultf32

Output of SASA calculation.

```zig
pub const SasaResult = SasaResultGen(f64);
pub const SasaResultf32 = SasaResultGen(f32);
```

**Fields:**

| Field | Type | Description |
|-------|------|-------------|
| `total_area` | `T` | Total SASA in Angstroms squared |
| `atom_areas` | `[]T` | Per-atom SASA values |
| `allocator` | `Allocator` | Allocator used for `atom_areas` |

### Methods

#### deinit

```zig
pub fn deinit(self: *SasaResult) void
```

Frees `atom_areas`. Always call via `defer`.

#### toF64

```zig
pub fn toF64(self: SasaResultf32, allocator: Allocator) !SasaResult
```

Converts an f32 result to f64. The caller owns the returned result and must call `deinit()`.

```zig
var result_f32 = try zsasa.shrake_rupley.calculateSasaf32(allocator, atoms, .{});
defer result_f32.deinit();

var result_f64 = try result_f32.toF64(allocator);
defer result_f64.deinit();
```

---

## Vec3 / Vec3f32

3D vector type used internally.

```zig
pub const Vec3 = Vec3Gen(f64);
pub const Vec3f32 = Vec3Gen(f32);
```

**Fields:**

| Field | Type | Description |
|-------|------|-------------|
| `x` | `T` | X component |
| `y` | `T` | Y component |
| `z` | `T` | Z component |

### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `add` | `fn(self, other) Self` | Vector addition |
| `sub` | `fn(self, other) Self` | Vector subtraction |
| `scale` | `fn(self, scalar: T) Self` | Scalar multiplication |
| `length` | `fn(self) T` | Euclidean length |
| `dot` | `fn(self, other) T` | Dot product |
| `fromF64` | `fn(Vec3) Self` | Convert from f64 Vec3 |
| `toF64` | `fn(self) Vec3` | Convert to f64 Vec3 |

---

## Precision

Enum for selecting precision mode.

```zig
pub const Precision = enum { f32, f64 };
```

### fromString

```zig
pub fn fromString(s: []const u8) ?Precision
```

Parses `"f32"`, `"single"`, `"f64"`, or `"double"`. Returns `null` for unknown values.

---

## FixedString4

Fixed-size string for atom/chain names (max 4 characters). Zero-allocation alternative to heap-allocated strings.

```zig
pub const FixedString4 = struct {
    data: [4]u8,
    len: u8,
};
```

| Method | Description |
|--------|-------------|
| `fromSlice([]const u8)` | Create from a string slice (truncates at 4 chars) |
| `slice() []const u8` | Get contents as a slice |
| `eql(*const FixedString4) bool` | Compare with another FixedString4 |
| `eqlSlice([]const u8) bool` | Compare with a string slice |

## FixedString5

Same as `FixedString4` but max 5 characters. Used for residue names to support mmCIF 5-character `comp_id` values.

# Classifier

Atom radius and polarity classifier for SASA calculation. Import via `zsasa.classifier`.

The classifier assigns van der Waals radii and polarity classes (polar/apolar) to atoms based on residue name and atom name, following FreeSASA's lookup strategy:

1. Look up by `(residue_name, atom_name)` exact match
2. Fall back to `("ANY", atom_name)`
3. Fall back to element-based radius guessing

---

## ClassifierType

```zig
pub const ClassifierType = enum {
    naccess,  // NACCESS-compatible radii (Joao Rodrigues)
    protor,   // ProtOr hybridization-based radii (Tsai et al. 1999)
    oons,     // OONS radii (Ooi, Oobatake, Nemethy, Scheraga)
};
```

### Methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `name` | `fn(self) []const u8` | Display name (`"NACCESS"`, `"ProtOr"`, `"OONS"`) |
| `fromString` | `fn([]const u8) ?ClassifierType` | Parse from string (case-insensitive) |

```zig
const ct = zsasa.classifier.ClassifierType.fromString("protor").?;
std.debug.print("Using: {s}\n", .{ct.name()}); // "ProtOr"
```

---

## AtomClass

Atom polarity classification.

```zig
pub const AtomClass = enum {
    polar,
    apolar,
    unknown,
};
```

| Method | Signature | Description |
|--------|-----------|-------------|
| `fromString` | `fn([]const u8) AtomClass` | Parse from string |
| `toString` | `fn(self) []const u8` | Convert to string |

---

## AtomProperties

Combined radius and class for an atom type.

```zig
pub const AtomProperties = struct {
    radius: f64,       // van der Waals radius in Angstroms
    class: AtomClass,  // Polar/apolar classification
};
```

---

## Classifier

Dynamic classifier that can be loaded from config files or built programmatically.

### init

```zig
pub fn init(allocator: Allocator, name: []const u8) !Classifier
```

Creates an empty classifier with the given name.

### deinit

```zig
pub fn deinit(self: *Classifier) void
```

Frees all resources.

### addAtom

```zig
pub fn addAtom(
    self: *Classifier,
    residue: []const u8,
    atom: []const u8,
    radius: f64,
    class: AtomClass,
) !void
```

Adds an atom type entry. Use `"ANY"` as residue for wildcard backbone entries.

### getRadius

```zig
pub fn getRadius(
    self: *const Classifier,
    residue: []const u8,
    atom: []const u8,
) ?f64
```

Looks up radius for an atom. Tries exact `(residue, atom)` match first, then `("ANY", atom)` fallback. Returns `null` if not found.

### getClass

```zig
pub fn getClass(
    self: *const Classifier,
    residue: []const u8,
    atom: []const u8,
) AtomClass
```

Looks up polarity class. Same fallback logic as `getRadius`. Returns `.unknown` if not found.

### getProperties

```zig
pub fn getProperties(
    self: *const Classifier,
    residue: []const u8,
    atom: []const u8,
) ?AtomProperties
```

Returns both radius and class. Returns `null` if not found.

### count

```zig
pub fn count(self: *const Classifier) usize
```

Returns the number of atom type entries.

### Example

```zig
const zsasa = @import("zsasa");

var cls = try zsasa.classifier.Classifier.init(allocator, "custom");
defer cls.deinit();

try cls.addAtom("ALA", "CA", 1.87, .apolar);
try cls.addAtom("ANY", "N", 1.65, .polar);

const radius = cls.getRadius("ALA", "CA"); // 1.87
const class = cls.getClass("ALA", "N");    // .polar (via ANY fallback)
```

---

## Element-Based Radius Guessing

Fallback functions for atoms not found in any classifier.

### guessRadius

```zig
pub fn guessRadius(element: []const u8) ?f64
```

Returns van der Waals radius for an element symbol (e.g., `"C"` -> `1.70`, `"N"` -> `1.55`).

### guessRadiusFromAtomicNumber

```zig
pub fn guessRadiusFromAtomicNumber(atomic_number: u8) ?f64
```

Returns radius by atomic number (e.g., `6` -> `1.70` for Carbon).

### guessRadiusFromAtomName

```zig
pub fn guessRadiusFromAtomName(atom_name: []const u8) ?f64
```

Extracts element from PDB atom name and returns radius. Handles the PDB convention where the first character is the element for common atoms.

### extractElement

```zig
pub fn extractElement(atom_name: []const u8) []const u8
```

Extracts element symbol from a PDB-style atom name (e.g., `"CA"` -> `"C"`, `"FE"` -> `"FE"`).

### Common Radii

| Element | Radius (A) | Atomic Number |
|---------|-----------|---------------|
| H | 1.10 | 1 |
| C | 1.70 | 6 |
| N | 1.55 | 7 |
| O | 1.52 | 8 |
| P | 1.80 | 15 |
| S | 1.80 | 16 |
| Se | 1.90 | 34 |

Sources: Mantina et al. 2009 for common elements, gemmi for others.

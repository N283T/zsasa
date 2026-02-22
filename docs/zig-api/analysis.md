# Analysis

Per-residue SASA aggregation, RSA calculation, and polar/nonpolar classification. Import via `zsasa.analysis`.

---

## aggregateByResidue

```zig
pub fn aggregateByResidue(
    allocator: std.mem.Allocator,
    input: AtomInput,
    atom_areas: []const f64,
) !ResidueResult
```

Groups atom SASA values into per-residue totals. Atoms are grouped by `(chain_id, residue_num, insertion_code)`.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `allocator` | `Allocator` | Memory allocator |
| `input` | `AtomInput` | Atom data with residue metadata |
| `atom_areas` | `[]const f64` | Per-atom SASA values from `calculateSasa()` |

**Returns:** `ResidueResult`

**Errors:**
- `MissingChainInfo` - `input.chain_id` is `null`
- `MissingResidueInfo` - `input.residue` is `null`
- `MissingResidueNumInfo` - `input.residue_num` is `null`
- `MissingInsertionCodeInfo` - `input.insertion_code` is `null`
- `LengthMismatch` - `atom_areas.len != input.atomCount()`

**Requires:** `input.hasResidueInfo()` must be `true` (all per-residue metadata populated). PDB and mmCIF parsers always populate these fields.

**Example:**

```zig
const zsasa = @import("zsasa");

// Calculate SASA
var result = try zsasa.shrake_rupley.calculateSasaParallel(allocator, atoms, .{}, 0);
defer result.deinit();

// Aggregate by residue
var residues = try zsasa.analysis.aggregateByResidue(allocator, atoms, result.atom_areas);
defer residues.deinit();

for (residues.residues) |*res| {
    res.calculateRsa();  // Calculate RSA for each residue
    const rsa_str = if (res.rsa) |rsa| rsa else -1.0;
    std.debug.print("{s}:{s}{d} SASA={d:.2} RSA={d:.3}\n", .{
        res.chain_id.slice(),
        res.residue_name.slice(),
        res.residue_num,
        res.sasa,
        rsa_str,
    });
}
```

---

## ResidueResult

```zig
pub const ResidueResult = struct {
    residues: []ResidueSasa,
    allocator: std.mem.Allocator,
};
```

### deinit

```zig
pub fn deinit(self: *ResidueResult) void
```

Frees `residues` array.

---

## ResidueSasa

Per-residue SASA data.

```zig
pub const ResidueSasa = struct {
    chain_id: FixedString4,
    residue_name: FixedString5,
    residue_num: i32,
    insertion_code: FixedString4,
    sasa: f64,
    atom_count: usize,
    rsa: ?f64 = null,
};
```

| Field | Type | Description |
|-------|------|-------------|
| `chain_id` | `FixedString4` | Chain identifier (e.g., "A") |
| `residue_name` | `FixedString5` | Residue name (e.g., "ALA") |
| `residue_num` | `i32` | Residue sequence number |
| `insertion_code` | `FixedString4` | PDB insertion code |
| `sasa` | `f64` | Total residue SASA in Angstroms squared |
| `atom_count` | `usize` | Number of atoms in residue |
| `rsa` | `?f64` | Relative Solvent Accessibility (`null` until calculated) |

### calculateRsa

```zig
pub fn calculateRsa(self: *ResidueSasa) void
```

Computes RSA from SASA and the residue's theoretical maximum SASA (from `MaxSASA`). Sets `self.rsa` to `sasa / max_sasa`, or `null` if the residue type is unknown.

RSA values > 1.0 are possible for exposed terminal residues.

---

## calculatePolarSummary

```zig
pub fn calculatePolarSummary(residues: []const ResidueSasa) PolarSummary
```

Summarizes SASA by residue polarity. Classifies each residue as polar or nonpolar and sums their SASA values.

**Polar residues:** ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR, TYR
**Nonpolar residues:** ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, TRP, VAL

**Example:**

```zig
const summary = zsasa.analysis.calculatePolarSummary(residues.residues);
std.debug.print("Polar: {d:.2} A^2 ({d:.1}%)\n", .{
    summary.polar_sasa,
    summary.polarFraction() * 100,
});
std.debug.print("Nonpolar: {d:.2} A^2 ({d:.1}%)\n", .{
    summary.nonpolar_sasa,
    summary.nonpolarFraction() * 100,
});
```

---

## PolarSummary

```zig
pub const PolarSummary = struct {
    polar_sasa: f64,
    nonpolar_sasa: f64,
    unknown_sasa: f64,
    polar_residue_count: usize,
    nonpolar_residue_count: usize,
    unknown_residue_count: usize,
};
```

| Method | Returns | Description |
|--------|---------|-------------|
| `polarFraction()` | `f64` | Fraction of polar SASA (excluding unknown) |
| `nonpolarFraction()` | `f64` | Fraction of nonpolar SASA (excluding unknown) |

---

## MaxSASA

Theoretical maximum SASA values for standard amino acids (Tien et al. 2013).

```zig
pub const MaxSASA = struct {
    pub const ALA: f64 = 129.0;
    pub const ARG: f64 = 274.0;
    // ... (20 standard amino acids)
};
```

### get

```zig
pub fn get(residue_name: []const u8) ?f64
```

Returns the maximum SASA for a 3-letter residue code. Returns `null` for unknown residues.

```zig
const max = zsasa.analysis.MaxSASA.get("ALA"); // 129.0
const unknown = zsasa.analysis.MaxSASA.get("UNK"); // null
```

---

## ResidueClass

Residue polarity classification.

```zig
pub const ResidueClass = enum {
    polar,
    nonpolar,
    unknown,
};
```

### fromResidueName

```zig
pub fn fromResidueName(residue_name: []const u8) ResidueClass
```

Classifies a residue by its 3-letter code.

---

## Full Example

```zig
const zsasa = @import("zsasa");

// Parse and calculate SASA
var parser = zsasa.pdb_parser.PdbParser.init(allocator);
var atoms = try parser.parseFile("protein.pdb");
defer atoms.deinit();

var result = try zsasa.shrake_rupley.calculateSasaParallel(allocator, atoms, .{}, 0);
defer result.deinit();

// Per-residue analysis
var residues = try zsasa.analysis.aggregateByResidue(allocator, atoms, result.atom_areas);
defer residues.deinit();

// Calculate RSA for each residue
for (residues.residues) |*res| {
    res.calculateRsa();
}

// Polar summary
const summary = zsasa.analysis.calculatePolarSummary(residues.residues);
std.debug.print("Total SASA: {d:.2} A^2\n", .{result.total_area});
std.debug.print("Polar:    {d:.2} A^2 ({d:.1}%)\n", .{
    summary.polar_sasa, summary.polarFraction() * 100,
});
std.debug.print("Nonpolar: {d:.2} A^2 ({d:.1}%)\n", .{
    summary.nonpolar_sasa, summary.nonpolarFraction() * 100,
});
std.debug.print("Residues: {d}\n", .{residues.residues.len});
```

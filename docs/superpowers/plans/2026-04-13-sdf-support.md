# SDF File Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add SDF file support to zsasa for direct SASA calculation and supplementary bond topology for CCD-unregistered compounds.

**Architecture:** New `sdf_parser.zig` module handles V2000/V3000 parsing and conversion to existing `Component`/`AtomInput` types. Integration into `calc.zig`, `batch.zig`, and `traj.zig` via `--sdf` CLI option and SDF format detection. Reuses existing `hybridization.zig` bond analysis for radii derivation.

**Tech Stack:** Zig, existing `hybridization.zig`/`classifier_ccd.zig` infrastructure

---

## File Map

| File | Action | Responsibility |
|------|--------|----------------|
| `src/sdf_parser.zig` | **Create** | SDF V2000/V3000 parser, `toStoredComponent`, `toAtomInput` |
| `src/format_detect.zig` | **Modify** | Add `.sdf`, `.mol` extensions, `sdf` variant to `InputFormat` |
| `src/calc.zig` | **Modify** | Add `--sdf` option parsing, SDF input branch in `readInputFile`, SDF component loading |
| `src/batch.zig` | **Modify** | Add `--sdf` option parsing and SDF component loading |
| `src/traj.zig` | **Modify** | Add `--sdf` option parsing and SDF component loading |
| `src/root.zig` | **Modify** | Export `sdf_parser` in public API |
| `src/main.zig` | **Modify** | Add `sdf_parser` test import |
| `test_data/ethanol_v2000.sdf` | **Create** | V2000 test SDF (ethanol, 9 atoms, 8 bonds) |
| `test_data/ethanol_v3000.sdf` | **Create** | V3000 test SDF (same ethanol) |
| `test_data/two_molecules.sdf` | **Create** | Multi-molecule SDF (methane + water) |

---

### Task 1: Create Test Data Files

**Files:**
- Create: `test_data/ethanol_v2000.sdf`
- Create: `test_data/ethanol_v3000.sdf`
- Create: `test_data/two_molecules.sdf`

- [ ] **Step 1: Create V2000 ethanol SDF**

Ethanol (C2H5OH): 9 atoms, 8 bonds. Coordinates from a standard geometry.

```
test_data/ethanol_v2000.sdf
```

```sdf
ethanol
     RDKit          3D

  9  8  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5200    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0800    1.2100    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3900    0.9800   -0.2600 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3900   -0.5400    0.8700 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3900   -0.4400   -0.9200 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9100   -0.5400    0.8700 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9100    0.5400   -0.8700 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.0400    1.2100    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  3  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
$$$$
```

- [ ] **Step 2: Create V3000 ethanol SDF**

Same ethanol molecule in V3000 format.

```
test_data/ethanol_v3000.sdf
```

```sdf
ethanol
     RDKit          3D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 8 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.0000 0.0000 0.0000 0
M  V30 2 C 1.5200 0.0000 0.0000 0
M  V30 3 O 2.0800 1.2100 0.0000 0
M  V30 4 H -0.3900 0.9800 -0.2600 0
M  V30 5 H -0.3900 -0.5400 0.8700 0
M  V30 6 H -0.3900 -0.4400 -0.9200 0
M  V30 7 H 1.9100 -0.5400 0.8700 0
M  V30 8 H 1.9100 0.5400 -0.8700 0
M  V30 9 H 3.0400 1.2100 0.0000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 1 4
M  V30 3 1 1 5
M  V30 4 1 1 6
M  V30 5 1 2 3
M  V30 6 1 2 7
M  V30 7 1 2 8
M  V30 8 1 3 9
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
```

- [ ] **Step 3: Create multi-molecule SDF (methane + water)**

```
test_data/two_molecules.sdf
```

```sdf
methane
     RDKit          3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6300    0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6300   -0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6300    0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6300   -0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
$$$$
water
     RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1173 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.7572    0.0000   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7572    0.0000   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
$$$$
```

- [ ] **Step 4: Commit test data**

```bash
git add test_data/ethanol_v2000.sdf test_data/ethanol_v3000.sdf test_data/two_molecules.sdf
git commit -m "test: add SDF test data files (ethanol V2000/V3000, multi-molecule)"
```

---

### Task 2: SDF Parser — Core V2000 Parsing

**Files:**
- Create: `src/sdf_parser.zig`

- [ ] **Step 1: Write failing tests for V2000 parsing**

Add at the bottom of `src/sdf_parser.zig`:

```zig
// At the top of the file:
const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const hybridization = @import("hybridization.zig");
const ccd_parser = @import("ccd_parser.zig");
const types = @import("types.zig");

// ... (implementation stubs come in step 3)

test "parse V2000 single molecule — ethanol" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/ethanol_v2000.sdf");
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    const mol = molecules[0];
    try std.testing.expectEqualStrings("ethanol", mol.name);
    try std.testing.expectEqual(@as(usize, 9), mol.atoms.len);
    try std.testing.expectEqual(@as(usize, 8), mol.bonds.len);

    // First atom: C at (0, 0, 0)
    try std.testing.expectEqual(elem.Element.C, mol.atoms[0].element);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), mol.atoms[0].x, 0.001);

    // Second atom: C at (1.52, 0, 0)
    try std.testing.expectEqual(elem.Element.C, mol.atoms[1].element);
    try std.testing.expectApproxEqAbs(@as(f64, 1.52), mol.atoms[1].x, 0.001);

    // Third atom: O
    try std.testing.expectEqual(elem.Element.O, mol.atoms[2].element);

    // Fourth atom: H
    try std.testing.expectEqual(elem.Element.H, mol.atoms[3].element);

    // Bond 1: atom 0-1, single
    try std.testing.expectEqual(@as(u16, 0), mol.bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), mol.bonds[0].atom_idx_2);
    try std.testing.expectEqual(hybridization.BondOrder.single, mol.bonds[0].order);
}

test "parse V2000 multi-molecule SDF" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/two_molecules.sdf");
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 2), molecules.len);
    try std.testing.expectEqualStrings("methane", molecules[0].name);
    try std.testing.expectEqual(@as(usize, 5), molecules[0].atoms.len);
    try std.testing.expectEqual(@as(usize, 4), molecules[0].bonds.len);
    try std.testing.expectEqualStrings("water", molecules[1].name);
    try std.testing.expectEqual(@as(usize, 3), molecules[1].atoms.len);
    try std.testing.expectEqual(@as(usize, 2), molecules[1].bonds.len);
}

test "parse MOL (no $$$$ terminator)" {
    const allocator = std.testing.allocator;
    // Simple methane without $$$$ at end
    const source =
        \\methane
        \\     RDKit          3D
        \\
        \\  5  4  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.6300    0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.6300   -0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.6300    0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.6300   -0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  2  1  0
        \\  1  3  1  0
        \\  1  4  1  0
        \\  1  5  1  0
        \\M  END
    ;
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    try std.testing.expectEqualStrings("methane", molecules[0].name);
    try std.testing.expectEqual(@as(usize, 5), molecules[0].atoms.len);
}

test "parse empty SDF returns error" {
    const allocator = std.testing.allocator;
    const result = parse(allocator, "");
    try std.testing.expectError(error.EmptySdf, result);
}

test "parse SDF with bad bond index returns error" {
    const allocator = std.testing.allocator;
    const source =
        \\bad_bond
        \\     RDKit          3D
        \\
        \\  2  1  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  5  1  0
        \\M  END
    ;
    const result = parse(allocator, source);
    try std.testing.expectError(error.BondIndexOutOfRange, result);
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | head -30
```

Expected: Compilation error — `parse` and `freeMolecules` not defined yet.

- [ ] **Step 3: Implement V2000 parsing**

Add the data types and V2000 parser above the tests in `src/sdf_parser.zig`:

```zig
const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const hybridization = @import("hybridization.zig");
const types = @import("types.zig");

pub const SdfAtom = struct {
    x: f64,
    y: f64,
    z: f64,
    element: elem.Element,
};

pub const SdfBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: hybridization.BondOrder,
};

pub const SdfMolecule = struct {
    name: []const u8,
    atoms: []const SdfAtom,
    bonds: []const SdfBond,
};

pub const SdfError = error{
    EmptySdf,
    InvalidCountsLine,
    InvalidAtomLine,
    InvalidBondLine,
    InvalidCoordinate,
    BondIndexOutOfRange,
    InvalidV3000,
    OutOfMemory,
};

/// Parse an SDF or MOL file. Returns all molecules.
/// Caller owns the returned slice; free with `freeMolecules`.
pub fn parse(allocator: Allocator, source: []const u8) SdfError![]const SdfMolecule {
    if (source.len == 0) return error.EmptySdf;

    var molecules = std.ArrayListUnmanaged(SdfMolecule){};
    errdefer {
        for (molecules.items) |mol| {
            allocator.free(mol.atoms);
            allocator.free(mol.bonds);
        }
        molecules.deinit(allocator);
    }

    var lines = std.mem.splitScalar(u8, source, '\n');
    var has_content = false;

    while (true) {
        const mol = parseSingleMolecule(allocator, &lines) catch |err| {
            if (err == error.EmptySdf and has_content) break; // EOF after $$$$
            return err;
        };
        has_content = true;
        molecules.append(allocator, mol) catch return error.OutOfMemory;

        // Scan for $$$$ delimiter or EOF
        var found_delim = false;
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "$$$$")) {
                found_delim = true;
                break;
            }
        }
        if (!found_delim) break; // EOF without $$$$, MOL compat
    }

    if (molecules.items.len == 0) return error.EmptySdf;
    return molecules.toOwnedSlice(allocator) catch return error.OutOfMemory;
}

/// Free molecules returned by `parse`.
pub fn freeMolecules(allocator: Allocator, molecules: []const SdfMolecule) void {
    for (molecules) |mol| {
        allocator.free(mol.atoms);
        allocator.free(mol.bonds);
    }
    allocator.free(molecules);
}

/// Parse a single molecule from the line iterator.
fn parseSingleMolecule(
    allocator: Allocator,
    lines: *std.mem.SplitIterator(u8, .scalar),
) SdfError!SdfMolecule {
    // Line 1: molecule name
    const name_line = lines.next() orelse return error.EmptySdf;
    const name = std.mem.trim(u8, stripCr(name_line), " ");

    // Line 2: program/timestamp (skip)
    _ = lines.next() orelse return error.InvalidCountsLine;
    // Line 3: comment (skip)
    _ = lines.next() orelse return error.InvalidCountsLine;

    // Line 4: counts line
    const counts_line = lines.next() orelse return error.InvalidCountsLine;
    const trimmed_counts = stripCr(counts_line);

    // Detect V3000
    if (std.mem.indexOf(u8, trimmed_counts, "V3000") != null) {
        return parseV3000Body(allocator, name, lines);
    }

    // V2000: parse "aaabbblll..." fixed-width counts
    // Atoms: cols 0-2 (3 chars), Bonds: cols 3-5 (3 chars)
    if (trimmed_counts.len < 6) return error.InvalidCountsLine;

    const n_atoms = std.fmt.parseInt(u32, std.mem.trim(u8, trimmed_counts[0..3], " "), 10) catch
        return error.InvalidCountsLine;
    const n_bonds = std.fmt.parseInt(u32, std.mem.trim(u8, trimmed_counts[3..6], " "), 10) catch
        return error.InvalidCountsLine;

    // Parse atom block
    var atoms = std.ArrayListUnmanaged(SdfAtom){};
    atoms.ensureTotalCapacity(allocator, n_atoms) catch return error.OutOfMemory;
    errdefer atoms.deinit(allocator);

    for (0..n_atoms) |_| {
        const line = lines.next() orelse return error.InvalidAtomLine;
        const atom = parseV2000AtomLine(stripCr(line)) orelse return error.InvalidAtomLine;
        atoms.appendAssumeCapacity(atom);
    }

    // Parse bond block
    var bonds = std.ArrayListUnmanaged(SdfBond){};
    bonds.ensureTotalCapacity(allocator, n_bonds) catch return error.OutOfMemory;
    errdefer bonds.deinit(allocator);

    for (0..n_bonds) |_| {
        const line = lines.next() orelse return error.InvalidBondLine;
        const bond = parseV2000BondLine(stripCr(line), n_atoms) catch |err| return err;
        bonds.appendAssumeCapacity(bond);
    }

    return .{
        .name = name,
        .atoms = atoms.toOwnedSlice(allocator) catch return error.OutOfMemory,
        .bonds = bonds.toOwnedSlice(allocator) catch return error.OutOfMemory,
    };
}

/// Parse a V2000 atom line:
///   xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
/// Coords: cols 0-9, 10-19, 20-29. Element: cols 31-33.
fn parseV2000AtomLine(line: []const u8) ?SdfAtom {
    if (line.len < 34) return null;

    const x = std.fmt.parseFloat(f64, std.mem.trim(u8, line[0..10], " ")) catch return null;
    const y = std.fmt.parseFloat(f64, std.mem.trim(u8, line[10..20], " ")) catch return null;
    const z = std.fmt.parseFloat(f64, std.mem.trim(u8, line[20..30], " ")) catch return null;

    // Element symbol: cols 31-33 (sometimes 30-33 depending on alignment)
    const elem_str = std.mem.trim(u8, line[30..@min(line.len, 34)], " ");
    const element = elem.fromSymbol(elem_str);

    return .{ .x = x, .y = y, .z = z, .element = element };
}

/// Parse a V2000 bond line:
///   111222tttsssxxxrrrccc
/// Atom1: cols 0-2, Atom2: cols 3-5, Type: cols 6-8.
fn parseV2000BondLine(line: []const u8, n_atoms: u32) SdfError!SdfBond {
    if (line.len < 9) return error.InvalidBondLine;

    const a1 = std.fmt.parseInt(u16, std.mem.trim(u8, line[0..3], " "), 10) catch
        return error.InvalidBondLine;
    const a2 = std.fmt.parseInt(u16, std.mem.trim(u8, line[3..6], " "), 10) catch
        return error.InvalidBondLine;
    const bond_type = std.fmt.parseInt(u8, std.mem.trim(u8, line[6..9], " "), 10) catch
        return error.InvalidBondLine;

    // SDF uses 1-based indices
    if (a1 == 0 or a2 == 0 or a1 > n_atoms or a2 > n_atoms) return error.BondIndexOutOfRange;

    return .{
        .atom_idx_1 = a1 - 1, // Convert to 0-based
        .atom_idx_2 = a2 - 1,
        .order = sdfBondOrder(bond_type),
    };
}

/// Map SDF bond type integer to BondOrder.
fn sdfBondOrder(bond_type: u8) hybridization.BondOrder {
    return switch (bond_type) {
        1 => .single,
        2 => .double,
        3 => .triple,
        4 => .aromatic,
        else => .unknown,
    };
}

/// Strip trailing \r (for Windows line endings in SDF files).
fn stripCr(line: []const u8) []const u8 {
    if (line.len > 0 and line[line.len - 1] == '\r') {
        return line[0 .. line.len - 1];
    }
    return line;
}

// V3000 stub — will be implemented in Task 3
fn parseV3000Body(
    allocator: Allocator,
    name: []const u8,
    lines: *std.mem.SplitIterator(u8, .scalar),
) SdfError!SdfMolecule {
    _ = allocator;
    _ = name;
    _ = lines;
    return error.InvalidV3000;
}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -5
```

Expected: All V2000 tests pass. The V3000 test (if any) would fail — but we haven't added V3000 tests yet.

- [ ] **Step 5: Commit**

```bash
git add src/sdf_parser.zig
git commit -m "feat: add SDF V2000 parser with multi-molecule support"
```

---

### Task 3: SDF Parser — V3000 Parsing

**Files:**
- Modify: `src/sdf_parser.zig`

- [ ] **Step 1: Write failing test for V3000**

Add to tests in `src/sdf_parser.zig`:

```zig
test "parse V3000 single molecule — ethanol" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/ethanol_v3000.sdf");
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    const mol = molecules[0];
    try std.testing.expectEqualStrings("ethanol", mol.name);
    try std.testing.expectEqual(@as(usize, 9), mol.atoms.len);
    try std.testing.expectEqual(@as(usize, 8), mol.bonds.len);

    // Same coordinates as V2000
    try std.testing.expectEqual(elem.Element.C, mol.atoms[0].element);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), mol.atoms[0].x, 0.001);
    try std.testing.expectEqual(elem.Element.C, mol.atoms[1].element);
    try std.testing.expectApproxEqAbs(@as(f64, 1.52), mol.atoms[1].x, 0.001);
    try std.testing.expectEqual(elem.Element.O, mol.atoms[2].element);

    // Same bonds as V2000
    try std.testing.expectEqual(@as(u16, 0), mol.bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), mol.bonds[0].atom_idx_2);
    try std.testing.expectEqual(hybridization.BondOrder.single, mol.bonds[0].order);
}
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | grep -A2 "V3000"
```

Expected: FAIL — `parseV3000Body` returns `error.InvalidV3000`.

- [ ] **Step 3: Implement V3000 parser**

Replace the `parseV3000Body` stub in `src/sdf_parser.zig`:

```zig
/// Parse V3000 CTAB body (after counts line detected V3000).
/// Expects the iterator positioned after the counts line.
fn parseV3000Body(
    allocator: Allocator,
    name: []const u8,
    lines: *std.mem.SplitIterator(u8, .scalar),
) SdfError!SdfMolecule {
    var n_atoms: u32 = 0;
    var n_bonds: u32 = 0;
    var atoms = std.ArrayListUnmanaged(SdfAtom){};
    errdefer atoms.deinit(allocator);
    var bonds = std.ArrayListUnmanaged(SdfBond){};
    errdefer bonds.deinit(allocator);

    var in_atom_block = false;
    var in_bond_block = false;

    while (lines.next()) |raw_line| {
        const line = stripCr(raw_line);

        // V3000 lines start with "M  V30 "
        if (!std.mem.startsWith(u8, line, "M  V30 ")) {
            // "M  END" terminates the molecule
            if (std.mem.startsWith(u8, line, "M  END")) break;
            continue;
        }

        const content = line["M  V30 ".len..];

        if (std.mem.startsWith(u8, content, "COUNTS ")) {
            // Parse: COUNTS natoms nbonds ...
            var iter = std.mem.tokenizeScalar(u8, content["COUNTS ".len..], ' ');
            n_atoms = std.fmt.parseInt(u32, iter.next() orelse return error.InvalidV3000, 10) catch
                return error.InvalidV3000;
            n_bonds = std.fmt.parseInt(u32, iter.next() orelse return error.InvalidV3000, 10) catch
                return error.InvalidV3000;
            atoms.ensureTotalCapacity(allocator, n_atoms) catch return error.OutOfMemory;
            bonds.ensureTotalCapacity(allocator, n_bonds) catch return error.OutOfMemory;
        } else if (std.mem.eql(u8, content, "BEGIN ATOM")) {
            in_atom_block = true;
        } else if (std.mem.eql(u8, content, "END ATOM")) {
            in_atom_block = false;
        } else if (std.mem.eql(u8, content, "BEGIN BOND")) {
            in_bond_block = true;
        } else if (std.mem.eql(u8, content, "END BOND")) {
            in_bond_block = false;
        } else if (std.mem.eql(u8, content, "END CTAB")) {
            // done
        } else if (in_atom_block) {
            // Format: index element x y z charge [...]
            const atom = parseV3000AtomLine(content) orelse return error.InvalidAtomLine;
            atoms.appendAssumeCapacity(atom);
        } else if (in_bond_block) {
            // Format: index type atom1 atom2 [...]
            const bond = parseV3000BondLine(content, n_atoms) catch |err| return err;
            bonds.appendAssumeCapacity(bond);
        }
    }

    if (atoms.items.len == 0) return error.EmptySdf;

    return .{
        .name = name,
        .atoms = atoms.toOwnedSlice(allocator) catch return error.OutOfMemory,
        .bonds = bonds.toOwnedSlice(allocator) catch return error.OutOfMemory,
    };
}

/// Parse V3000 atom line: "index element x y z charge [...]"
fn parseV3000AtomLine(content: []const u8) ?SdfAtom {
    var iter = std.mem.tokenizeScalar(u8, content, ' ');
    _ = iter.next() orelse return null; // index (skip)
    const elem_str = iter.next() orelse return null; // element symbol
    const x_str = iter.next() orelse return null;
    const y_str = iter.next() orelse return null;
    const z_str = iter.next() orelse return null;

    const x = std.fmt.parseFloat(f64, x_str) catch return null;
    const y = std.fmt.parseFloat(f64, y_str) catch return null;
    const z = std.fmt.parseFloat(f64, z_str) catch return null;
    const element = elem.fromSymbol(elem_str);

    return .{ .x = x, .y = y, .z = z, .element = element };
}

/// Parse V3000 bond line: "index type atom1 atom2 [...]"
fn parseV3000BondLine(content: []const u8, n_atoms: u32) SdfError!SdfBond {
    var iter = std.mem.tokenizeScalar(u8, content, ' ');
    _ = iter.next() orelse return error.InvalidBondLine; // index (skip)
    const type_str = iter.next() orelse return error.InvalidBondLine;
    const a1_str = iter.next() orelse return error.InvalidBondLine;
    const a2_str = iter.next() orelse return error.InvalidBondLine;

    const bond_type = std.fmt.parseInt(u8, type_str, 10) catch return error.InvalidBondLine;
    const a1 = std.fmt.parseInt(u16, a1_str, 10) catch return error.InvalidBondLine;
    const a2 = std.fmt.parseInt(u16, a2_str, 10) catch return error.InvalidBondLine;

    if (a1 == 0 or a2 == 0 or a1 > n_atoms or a2 > n_atoms) return error.BondIndexOutOfRange;

    return .{
        .atom_idx_1 = a1 - 1,
        .atom_idx_2 = a2 - 1,
        .order = sdfBondOrder(bond_type),
    };
}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -5
```

Expected: All tests pass including V3000.

- [ ] **Step 5: Commit**

```bash
git add src/sdf_parser.zig
git commit -m "feat: add SDF V3000 parser support"
```

---

### Task 4: SDF Parser — `toStoredComponent` and `toAtomInput`

**Files:**
- Modify: `src/sdf_parser.zig`

- [ ] **Step 1: Write failing tests for toComponent and toAtomInput**

Add to tests in `src/sdf_parser.zig`:

```zig
test "toStoredComponent — ethanol molecule" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/ethanol_v2000.sdf");
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    var comp = try toStoredComponent(allocator, &molecules[0]);
    defer comp.deinit();

    // comp_id should be "ethan" (truncated to 5)
    const view = comp.view();
    try std.testing.expectEqualStrings("ethan", view.compIdSlice());
    try std.testing.expectEqual(@as(usize, 9), comp.atoms.len);
    try std.testing.expectEqual(@as(usize, 8), comp.bonds.len);

    // First atom type_symbol should be "C"
    try std.testing.expectEqualStrings("C", comp.atoms[0].typeSymbolSlice());
    // Third atom type_symbol should be "O"
    try std.testing.expectEqualStrings("O", comp.atoms[2].typeSymbolSlice());
}

test "toAtomInput — two molecules get separate chains" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/two_molecules.sdf");
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    var input = try toAtomInput(allocator, molecules, false);
    defer input.deinit();

    // methane(5) + water(3) = 8 total atoms
    try std.testing.expectEqual(@as(usize, 8), input.atomCount());

    // Chain IDs: methane=A, water=B
    const chains = input.chain_id.?;
    try std.testing.expectEqualStrings("A", chains[0].slice()); // methane atom
    try std.testing.expectEqualStrings("A", chains[4].slice()); // last methane atom
    try std.testing.expectEqualStrings("B", chains[5].slice()); // first water atom
    try std.testing.expectEqualStrings("B", chains[7].slice()); // last water atom

    // Residue names
    const residues = input.residue.?;
    try std.testing.expectEqualStrings("metha", residues[0].slice()); // truncated "methane"
    try std.testing.expectEqualStrings("water", residues[5].slice());
}

test "toAtomInput — skip hydrogens" {
    const allocator = std.testing.allocator;
    const source = @embedFile("../test_data/ethanol_v2000.sdf");
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    var input = try toAtomInput(allocator, molecules, true);
    defer input.deinit();

    // Ethanol has 3 heavy atoms (2C + 1O), 6 H skipped
    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | head -10
```

Expected: Compilation error — `toStoredComponent` and `toAtomInput` not defined.

- [ ] **Step 3: Implement toStoredComponent**

Add to `src/sdf_parser.zig` (before the tests):

```zig
/// Convert an SdfMolecule to a ccd_parser.StoredComponent for direct insertion
/// into a ComponentDict. The molecule name (truncated to 5 chars) becomes the comp_id.
/// Atom names are generated as element + per-element index (C1, C2, O1...).
/// Caller owns the returned StoredComponent; free with `.deinit()`.
pub fn toStoredComponent(
    allocator: Allocator,
    molecule: *const SdfMolecule,
) !ccd_parser.StoredComponent {
    const atoms = try allocator.alloc(hybridization.CompAtom, molecule.atoms.len);
    errdefer allocator.free(atoms);

    // Count per-element for generating atom names (C1, C2, O1...)
    var element_counts: [119]u16 = .{0} ** 119; // indexed by atomic number

    for (molecule.atoms, 0..) |sdf_atom, i| {
        const sym = sdf_atom.element.symbol();
        const atomic_num = sdf_atom.element.atomicNumber();
        element_counts[atomic_num] += 1;

        // Generate atom name: element + count (e.g., "C1", "O2")
        var name_buf: [4]u8 = .{ ' ', ' ', ' ', ' ' };
        var name_len: u3 = 0;

        // Copy element symbol
        for (sym) |c| {
            if (name_len < 4) {
                name_buf[name_len] = c;
                name_len += 1;
            }
        }

        // Append count as single digit if there's room
        if (name_len < 4) {
            const count = element_counts[atomic_num];
            if (count < 10) {
                name_buf[name_len] = '0' + @as(u8, @intCast(count));
                name_len += 1;
            } else if (count < 100 and name_len + 1 < 4) {
                name_buf[name_len] = '0' + @as(u8, @intCast(count / 10));
                name_len += 1;
                name_buf[name_len] = '0' + @as(u8, @intCast(count % 10));
                name_len += 1;
            }
        }

        var comp_atom = hybridization.CompAtom.init(name_buf[0..name_len], sym);
        comp_atom.aromatic = false;
        comp_atom.leaving = false;

        atoms[i] = comp_atom;
    }

    const bonds = try allocator.alloc(hybridization.CompBond, molecule.bonds.len);
    errdefer allocator.free(bonds);

    for (molecule.bonds, 0..) |sdf_bond, i| {
        bonds[i] = .{
            .atom_idx_1 = sdf_bond.atom_idx_1,
            .atom_idx_2 = sdf_bond.atom_idx_2,
            .order = sdf_bond.order,
            .aromatic = (sdf_bond.order == .aromatic),
        };
    }

    // comp_id from molecule name (truncated to 5 chars)
    var comp_id: [5]u8 = .{ 0, 0, 0, 0, 0 };
    const id_len: u3 = @intCast(@min(molecule.name.len, 5));
    for (molecule.name[0..id_len], 0..) |c, j| {
        comp_id[j] = c;
    }

    return .{
        .comp_id = comp_id,
        .comp_id_len = id_len,
        .atoms = atoms,
        .bonds = bonds,
        .allocator = allocator,
    };
}
```

- [ ] **Step 4: Implement toAtomInput**

Add to `src/sdf_parser.zig`:

```zig
/// Convert parsed SDF molecules to AtomInput for SASA calculation.
/// Each molecule becomes a separate chain (A, B, C...).
/// Residue name = molecule name (truncated to 5 chars), residue number = 1.
/// When `skip_hydrogens` is true, hydrogen atoms are excluded.
pub fn toAtomInput(
    allocator: Allocator,
    molecules: []const SdfMolecule,
    skip_hydrogens: bool,
) !types.AtomInput {
    // Count total atoms
    var total: usize = 0;
    for (molecules) |mol| {
        for (mol.atoms) |atom| {
            if (skip_hydrogens and atom.element == .H) continue;
            total += 1;
        }
    }

    if (total == 0) return error.NoAtomsFound;

    // Allocate all arrays
    const x = try allocator.alloc(f64, total);
    errdefer allocator.free(x);
    const y = try allocator.alloc(f64, total);
    errdefer allocator.free(y);
    const z = try allocator.alloc(f64, total);
    errdefer allocator.free(z);
    const r = try allocator.alloc(f64, total);
    errdefer allocator.free(r);
    const element_arr = try allocator.alloc(u8, total);
    errdefer allocator.free(element_arr);
    const atom_names = try allocator.alloc(types.FixedString4, total);
    errdefer allocator.free(atom_names);
    const residues = try allocator.alloc(types.FixedString5, total);
    errdefer allocator.free(residues);
    const chain_ids = try allocator.alloc(types.FixedString4, total);
    errdefer allocator.free(chain_ids);
    const residue_nums = try allocator.alloc(i32, total);
    errdefer allocator.free(residue_nums);
    const insertion_codes = try allocator.alloc(types.FixedString4, total);
    errdefer allocator.free(insertion_codes);

    var idx: usize = 0;
    const max_chains: usize = 26; // A-Z

    for (molecules, 0..) |mol, mol_idx| {
        if (mol_idx >= max_chains) break;

        const chain_letter: u8 = 'A' + @as(u8, @intCast(mol_idx));
        const chain_id = types.FixedString4.fromSlice(&.{chain_letter});
        const res_name = types.FixedString5.fromSlice(mol.name[0..@min(mol.name.len, 5)]);

        // Per-element counters for atom name generation
        var element_counts: [119]u16 = .{0} ** 119;

        for (mol.atoms) |atom| {
            if (skip_hydrogens and atom.element == .H) continue;

            x[idx] = atom.x;
            y[idx] = atom.y;
            z[idx] = atom.z;
            r[idx] = atom.element.vdwRadius(); // Default VdW; classifier will override
            element_arr[idx] = atom.element.atomicNumber();
            residues[idx] = res_name;
            chain_ids[idx] = chain_id;
            residue_nums[idx] = 1;
            insertion_codes[idx] = types.FixedString4.fromSlice("");

            // Generate atom name
            const sym = atom.element.symbol();
            const atomic_num = atom.element.atomicNumber();
            element_counts[atomic_num] += 1;
            var name_buf: [4]u8 = .{ ' ', ' ', ' ', ' ' };
            var name_len: usize = 0;
            for (sym) |c| {
                if (name_len < 4) {
                    name_buf[name_len] = c;
                    name_len += 1;
                }
            }
            const count = element_counts[atomic_num];
            if (count < 10 and name_len < 4) {
                name_buf[name_len] = '0' + @as(u8, @intCast(count));
                name_len += 1;
            } else if (count < 100 and name_len + 1 < 4) {
                name_buf[name_len] = '0' + @as(u8, @intCast(count / 10));
                name_len += 1;
                name_buf[name_len] = '0' + @as(u8, @intCast(count % 10));
                name_len += 1;
            }
            atom_names[idx] = types.FixedString4.fromSlice(name_buf[0..name_len]);

            idx += 1;
        }
    }

    return .{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .element = element_arr,
        .atom_name = atom_names,
        .residue = residues,
        .chain_id = chain_ids,
        .residue_num = residue_nums,
        .insertion_code = insertion_codes,
        .allocator = allocator,
    };
}
```

Note: also add `const NoAtomsFound = error.NoAtomsFound;` or add it to the `SdfError` set:

```zig
pub const SdfError = error{
    EmptySdf,
    InvalidCountsLine,
    InvalidAtomLine,
    InvalidBondLine,
    InvalidCoordinate,
    BondIndexOutOfRange,
    InvalidV3000,
    NoAtomsFound,
    OutOfMemory,
};
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -5
```

Expected: All tests pass.

- [ ] **Step 6: Commit**

```bash
git add src/sdf_parser.zig
git commit -m "feat: add toStoredComponent and toAtomInput for SDF molecules"
```

---

### Task 5: Format Detection — Add SDF/MOL Support

**Files:**
- Modify: `src/format_detect.zig`

- [ ] **Step 1: Write failing tests**

Add to the test blocks in `src/format_detect.zig`:

```zig
test "isSupportedFile accepts SDF/MOL files" {
    try std.testing.expect(isSupportedFile("molecule.sdf"));
    try std.testing.expect(isSupportedFile("molecule.sdf.gz"));
    try std.testing.expect(isSupportedFile("molecule.SDF"));
    try std.testing.expect(isSupportedFile("molecule.mol"));
    try std.testing.expect(isSupportedFile("molecule.mol.gz"));
    try std.testing.expect(isSupportedFile("molecule.MOL"));
}

test "detectInputFormat handles SDF/MOL files" {
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.SDF"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.MOL"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf.gz"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol.gz"));
}
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | head -10
```

Expected: Compilation error — `InputFormat.sdf` doesn't exist.

- [ ] **Step 3: Add SDF format to format_detect.zig**

In `src/format_detect.zig`, make these changes:

Add `sdf` to the `InputFormat` enum:
```zig
pub const InputFormat = enum {
    json,
    mmcif,
    pdb,
    sdf,
};
```

Add SDF/MOL extensions to `supported_extensions`:
```zig
pub const supported_extensions = [_][]const u8{
    ".json",
    ".json.gz",
    ".pdb",
    ".pdb.gz",
    ".cif",
    ".cif.gz",
    ".mmcif",
    ".mmcif.gz",
    ".ent",
    ".ent.gz",
    ".sdf",
    ".sdf.gz",
    ".mol",
    ".mol.gz",
};
```

Add uppercase SDF/MOL checks in `isSupportedFile`:
```zig
    if (std.mem.endsWith(u8, name, ".SDF")) return true;
    if (std.mem.endsWith(u8, name, ".MOL")) return true;
```

Add SDF/MOL detection in `detectInputFormat` (before the "Default to JSON" line):
```zig
    // SDF/MOL formats
    if (std.mem.endsWith(u8, base, ".sdf")) return .sdf;
    if (std.mem.endsWith(u8, base, ".SDF")) return .sdf;
    if (std.mem.endsWith(u8, base, ".mol")) return .sdf;
    if (std.mem.endsWith(u8, base, ".MOL")) return .sdf;
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -5
```

Expected: All tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/format_detect.zig
git commit -m "feat: add SDF/MOL file format detection"
```

---

### Task 6: Wire SDF Parser into root.zig and main.zig

**Files:**
- Modify: `src/root.zig`
- Modify: `src/main.zig`

- [ ] **Step 1: Add sdf_parser to root.zig**

Add after the existing `pub const ccd_binary` line in `src/root.zig`:

```zig
pub const sdf_parser = @import("sdf_parser.zig");
```

And add to the `test` block:
```zig
    _ = sdf_parser;
```

- [ ] **Step 2: Run zig build test to verify compilation**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -5
```

Expected: All tests pass (sdf_parser tests now discovered through root.zig).

- [ ] **Step 3: Commit**

```bash
git add src/root.zig
git commit -m "feat: export sdf_parser in public API"
```

---

### Task 7: calc.zig — SDF Input Support (Use Case 1)

**Files:**
- Modify: `src/calc.zig`

- [ ] **Step 1: Add sdf_parser import**

Add at the top of `src/calc.zig` alongside other imports:

```zig
const sdf_parser = @import("sdf_parser.zig");
```

- [ ] **Step 2: Add --sdf option to CalcArgs**

Add a new field to the `CalcArgs` struct (after `ccd_path`):

```zig
    sdf_paths: std.BoundedArray([]const u8, 16) = .{}, // --sdf=PATH (up to 16)
```

- [ ] **Step 3: Add --sdf argument parsing in parseArgs**

Add before the `--use-bitmask` else-if block (around line 395) in `parseArgs`:

```zig
        // --sdf=PATH or --sdf PATH (SDF file for CCD classifier)
        else if (std.mem.startsWith(u8, arg, "--sdf=")) {
            const value = arg["--sdf=".len..];
            result.sdf_paths.append(value) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--sdf")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --sdf\n", .{});
                std.process.exit(1);
            }
            result.sdf_paths.append(args[i]) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
        }
```

- [ ] **Step 4: Add SDF branch in readInputFile**

Modify the `readInputFile` function to handle `InputFormat.sdf`. Add a new branch:

```zig
fn readInputFile(allocator: std.mem.Allocator, path: []const u8, args: CalcArgs) !ReadResult {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => .{ .input = try json_parser.readAtomInputFromFile(allocator, path) },
        .mmcif => blk: {
            // ... (existing mmcif code unchanged)
        },
        .pdb => blk: {
            // ... (existing pdb code unchanged)
        },
        .sdf => blk: {
            // Read SDF file (plain or .gz)
            const source = if (std.mem.endsWith(u8, path, ".gz"))
                try gzip.readGzip(allocator, path)
            else file_blk: {
                const f = try std.fs.cwd().openFile(path, .{});
                defer f.close();
                break :file_blk try f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024);
            };
            defer allocator.free(source);

            const molecules = try sdf_parser.parse(allocator, source);
            defer sdf_parser.freeMolecules(allocator, molecules);

            break :blk .{
                .input = try sdf_parser.toAtomInput(allocator, molecules, !args.include_hydrogens),
            };
        },
    };
}
```

- [ ] **Step 5: Add SDF-specific default classifier behavior**

In the `run` function, after the existing `effective_classifier` logic (around line 879), modify the default classifier selection to also default to CCD for SDF input:

Change:
```zig
    const effective_classifier: ?ClassifierType = args.classifier_type orelse
        if (args.config_path == null and input_format != .json) .ccd else null;
```

This already covers SDF since `input_format` will be `.sdf` which is `!= .json`. No change needed.

- [ ] **Step 6: Load SDF components into CCD classifier**

In the `run` function, after external CCD loading (around line 941) and before `applyBuiltinClassifier`, add SDF component loading. Find the section that calls `applyBuiltinClassifier` and add SDF loading between ext_ccd setup and that call.

Add a new helper function and update the classifier flow:

```zig
/// Load SDF files and convert molecules to ComponentDict entries.
fn loadSdfComponents(
    allocator: std.mem.Allocator,
    sdf_paths: []const []const u8,
    quiet: bool,
) !?ccd_parser.ComponentDict {
    if (sdf_paths.len == 0) return null;

    var dict = ccd_parser.ComponentDict.init(allocator);
    errdefer dict.deinit();

    for (sdf_paths) |sdf_path| {
        const source = if (std.mem.endsWith(u8, sdf_path, ".gz"))
            gzip.readGzip(allocator, sdf_path) catch |err| {
                std.debug.print("Error reading SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            }
        else file_blk: {
            const f = std.fs.cwd().openFile(sdf_path, .{}) catch |err| {
                std.debug.print("Error opening SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            };
            defer f.close();
            break :file_blk f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024) catch |err| {
                std.debug.print("Error reading SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            };
        };
        defer allocator.free(source);

        const molecules = sdf_parser.parse(allocator, source) catch |err| {
            std.debug.print("Error parsing SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
            std.process.exit(1);
        };
        defer sdf_parser.freeMolecules(allocator, molecules);

        for (molecules) |mol| {
            if (mol.name.len == 0) {
                if (!quiet) std.debug.print("Warning: SDF molecule has no name, skipping\n", .{});
                continue;
            }
            const stored = sdf_parser.toStoredComponent(allocator, &mol) catch |err| {
                if (!quiet) std.debug.print("Warning: Could not convert SDF molecule '{s}': {s}\n", .{ mol.name, @errorName(err) });
                continue;
            };

            const comp_id_str = mol.name[0..@min(mol.name.len, 5)];
            const dict_key = allocator.dupe(u8, comp_id_str) catch {
                var s = stored;
                s.deinit();
                continue;
            };
            dict.owned_keys.append(allocator, dict_key) catch {
                allocator.free(dict_key);
                var s = stored;
                s.deinit();
                continue;
            };
            dict.components.put(dict_key, stored) catch continue;
        }

        if (!quiet) {
            std.debug.print("SDF: loaded from '{s}'\n", .{sdf_path});
        }
    }

    if (dict.components.count() == 0) {
        dict.deinit();
        return null;
    }
    return dict;
}
```

Then in the `run` function, after the `ext_ccd` loading block, add:

```zig
            // Load SDF-derived components
            var sdf_ccd: ?ccd_parser.ComponentDict = null;
            if (effective_args.sdf_paths.len > 0) {
                sdf_ccd = try loadSdfComponents(allocator, effective_args.sdf_paths.constSlice(), effective_args.quiet);
            }
            defer if (sdf_ccd) |*d| d.deinit();
```

Then update `applyBuiltinClassifier` call to pass SDF dict. This requires modifying the function signature to accept a third dict source. The cleanest approach is to update the priority order inside `applyBuiltinClassifier`.

Add a new parameter `sdf_ccd` to `applyBuiltinClassifier`:

```zig
fn applyBuiltinClassifier(
    input: *types.AtomInput,
    ct: ClassifierType,
    sdf_ccd: ?*const ccd_parser.ComponentDict,
    inline_ccd: ?*const ccd_parser.ComponentDict,
    external_ccd: ?*const ccd_parser.ComponentDict,
    quiet: bool,
) !void {
```

And update the `dicts` array to include all three sources in priority order:

```zig
            const dicts: [3]?*const ccd_parser.ComponentDict = .{ sdf_ccd, inline_ccd, external_ccd };
```

Update the call site:
```zig
            const sdf_ccd_ptr: ?*const ccd_parser.ComponentDict = if (sdf_ccd != null) &sdf_ccd.? else null;
            applyBuiltinClassifier(&input, ct, sdf_ccd_ptr, inline_ccd, ext_ccd_ptr, args.quiet) catch |err| {
```

- [ ] **Step 7: Add --sdf to help text**

In `printHelp`, add after the `--ccd=PATH` line:

```
        \\    --sdf=PATH         SDF file with bond topology for CCD classifier
        \\                       Can be specified multiple times for multiple ligands
```

- [ ] **Step 8: Verify build compiles**

```bash
cd /Users/nagaet/freesasa-zig && zig build 2>&1
```

Expected: Clean build with no errors.

- [ ] **Step 9: Commit**

```bash
git add src/calc.zig
git commit -m "feat: add SDF input support and --sdf option to calc subcommand"
```

---

### Task 8: batch.zig and traj.zig — Add --sdf Option

**Files:**
- Modify: `src/batch.zig`
- Modify: `src/traj.zig`

The implementation for `batch.zig` and `traj.zig` follows the same pattern as `calc.zig`. Both files already have `--ccd` option handling and CCD classifier integration. The `--sdf` option should be added alongside.

- [ ] **Step 1: Add --sdf to batch.zig**

Import `sdf_parser`:
```zig
const sdf_parser = @import("sdf_parser.zig");
```

Add `sdf_paths` field to `BatchArgs` (the batch args struct), following the same pattern as `CalcArgs`:
```zig
    sdf_paths: std.BoundedArray([]const u8, 16) = .{},
```

Add `--sdf` argument parsing in the batch `parseArgs` function (same pattern as calc).

Add SDF loading in the batch `run` function alongside existing CCD loading. Find where `ext_ccd` is loaded and add `sdf_ccd` loading after it. Update the `applyBuiltinClassifier` call (or equivalent) to pass `sdf_ccd`.

Add `--sdf` to the batch help text.

- [ ] **Step 2: Add --sdf to traj.zig**

Same changes as batch.zig:
- Import `sdf_parser`
- Add `sdf_paths` field to `TrajArgs`
- Add `--sdf` argument parsing
- Add SDF loading in `run` function
- Add `--sdf` to help text

- [ ] **Step 3: Update applyBuiltinClassifier signature in batch.zig and traj.zig**

Both `batch.zig` and `traj.zig` have their own inline CCD classifier flows (they don't share `calc.zig`'s `applyBuiltinClassifier`). Find the equivalent CCD classifier setup code and add SDF component loading there.

In `batch.zig` (around line 358) where `ccd_clf.?.addComponent` is called in the loop, add SDF dict to the `dicts` array:

```zig
const dicts: [3]?*const ccd_parser.ComponentDict = .{ sdf_ccd_ptr, inline_ccd, ext_ccd_ptr };
```

Same for `traj.zig` (around line 1095).

- [ ] **Step 4: Run build**

```bash
cd /Users/nagaet/freesasa-zig && zig build 2>&1
```

Expected: Clean build.

- [ ] **Step 5: Commit**

```bash
git add src/batch.zig src/traj.zig
git commit -m "feat: add --sdf option to batch and traj subcommands"
```

---

### Task 9: SDF Input for calc — CCD Component Auto-Registration

**Files:**
- Modify: `src/calc.zig`

When SDF is the input file (use case 1), the SDF bond topology should be automatically converted to Components and registered in the CCD classifier. This happens when no explicit `--sdf` is provided but the input is an SDF file.

- [ ] **Step 1: Add SDF auto-component registration**

In `calc.zig`'s `run` function, after `readInputFile` and before the classifier section, add logic to auto-register SDF molecules when the input is an SDF file:

```zig
    // For SDF input: auto-register molecules as CCD components
    var sdf_auto_ccd: ?ccd_parser.ComponentDict = null;
    if (input_format == .sdf) {
        // Re-read and parse the SDF file to extract components
        // (readInputFile already consumed the source for coordinates)
        const sdf_source = if (std.mem.endsWith(u8, input_path, ".gz"))
            try gzip.readGzip(allocator, input_path)
        else file_blk: {
            const f = try std.fs.cwd().openFile(input_path, .{});
            defer f.close();
            break :file_blk try f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024);
        };
        defer allocator.free(sdf_source);

        const molecules = sdf_parser.parse(allocator, sdf_source) catch |err| {
            if (!effective_args.quiet) {
                std.debug.print("Warning: Could not parse SDF for component extraction: {s}\n", .{@errorName(err)});
            }
            break; // skip auto-registration
        };
        defer sdf_parser.freeMolecules(allocator, molecules);

        var dict = ccd_parser.ComponentDict.init(allocator);
        for (molecules) |mol| {
            if (mol.name.len == 0) continue;
            const comp = sdf_parser.toComponent(allocator, &mol) catch continue;
            dict.putComponent(mol.name[0..@min(mol.name.len, 5)], comp) catch continue;
        }
        if (dict.components.count() > 0) {
            sdf_auto_ccd = dict;
        } else {
            dict.deinit();
        }
    }
    defer if (sdf_auto_ccd) |*d| d.deinit();
```

However, re-reading the file is wasteful. A better approach: modify `ReadResult` to optionally carry parsed SDF molecules, or store the component dict alongside the input.

Alternative cleaner approach — extend `ReadResult`:

Add a field to `ReadResult`:
```zig
const ReadResult = struct {
    input: types.AtomInput,
    mmcif: ?mmcif_parser.MmcifParser = null,
    sdf_ccd: ?ccd_parser.ComponentDict = null,

    fn deinitCcd(self: *ReadResult) void {
        if (self.mmcif) |*p| p.deinitCcd();
        if (self.sdf_ccd) |*d| d.deinit();
    }
};
```

Then in the `.sdf` branch of `readInputFile`, after parsing, build the component dict from the molecules before freeing them:

```zig
        .sdf => blk: {
            const source = if (std.mem.endsWith(u8, path, ".gz"))
                try gzip.readGzip(allocator, path)
            else file_blk: {
                const f = try std.fs.cwd().openFile(path, .{});
                defer f.close();
                break :file_blk try f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024);
            };
            defer allocator.free(source);

            const molecules = try sdf_parser.parse(allocator, source);
            defer sdf_parser.freeMolecules(allocator, molecules);

            const input = try sdf_parser.toAtomInput(allocator, molecules, !args.include_hydrogens);

            // Build CCD component dict from SDF bond topology
            var sdf_dict = ccd_parser.ComponentDict.init(allocator);
            for (molecules) |mol| {
                if (mol.name.len == 0) continue;
                const stored = sdf_parser.toStoredComponent(allocator, &mol) catch continue;
                const comp_id_str = mol.name[0..@min(mol.name.len, 5)];
                const dict_key = allocator.dupe(u8, comp_id_str) catch {
                    var s = stored;
                    s.deinit();
                    continue;
                };
                sdf_dict.owned_keys.append(allocator, dict_key) catch {
                    allocator.free(dict_key);
                    var s = stored;
                    s.deinit();
                    continue;
                };
                sdf_dict.components.put(dict_key, stored) catch continue;
            }
            const has_components = sdf_dict.components.count() > 0;

            break :blk .{
                .input = input,
                .sdf_ccd = if (has_components) sdf_dict else null_blk: {
                    sdf_dict.deinit();
                    break :null_blk null;
                },
            };
        },
```

Then in the classifier section, use `read_result.sdf_ccd` as another CCD source:

```zig
            const sdf_auto_ccd: ?*const ccd_parser.ComponentDict = if (read_result.sdf_ccd) |*d| d else null;
```

And pass it to `applyBuiltinClassifier` alongside the --sdf derived components. The SDF auto dict should have lower priority than explicit --sdf but higher than inline CCD.

- [ ] **Step 2: Run build and test**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -5
```

Expected: All tests pass.

- [ ] **Step 3: Commit**

```bash
git add src/calc.zig
git commit -m "feat: auto-register SDF bond topology as CCD components for SDF input"
```

---

### Task 10: Integration Test — End-to-End SDF Calculation

**Files:**
- No new files (test via CLI)

- [ ] **Step 1: Build the binary**

```bash
cd /Users/nagaet/freesasa-zig && zig build
```

- [ ] **Step 2: Test use case 1 — SDF as input**

```bash
./zig-out/bin/zsasa calc test_data/ethanol_v2000.sdf --include-hetatm --timing
```

Expected: SASA calculation completes. Output shows:
- Total SASA value (should be reasonable for ethanol, ~100-200 Å²)
- "Classifier 'ccd'" line (default classifier)

- [ ] **Step 3: Test use case 1 with explicit classifier**

```bash
./zig-out/bin/zsasa calc test_data/ethanol_v2000.sdf --classifier naccess --include-hetatm
```

Expected: Runs without error, uses NACCESS classifier.

- [ ] **Step 4: Test use case 1 with V3000**

```bash
./zig-out/bin/zsasa calc test_data/ethanol_v3000.sdf --include-hetatm
```

Expected: Same SASA result as V2000 (within floating-point tolerance).

- [ ] **Step 5: Test multi-molecule SDF**

```bash
./zig-out/bin/zsasa calc test_data/two_molecules.sdf --include-hetatm
```

Expected: Output shows two chains (A for methane, B for water).

- [ ] **Step 6: Test use case 2 — --sdf with PDB** (requires a test PDB with non-standard residue)

If there is no suitable test PDB with a matching residue name, this step can be validated by checking the CLI accepts the option:

```bash
./zig-out/bin/zsasa calc test_data/1ubq.pdb --sdf test_data/ethanol_v2000.sdf 2>&1
```

Expected: Runs. Warning message about "SDF molecule 'ethanol' not found in input structure" (since 1ubq doesn't have an "ethan" residue). This confirms the --sdf pipeline works even when no match is found.

- [ ] **Step 7: Run full test suite**

```bash
cd /Users/nagaet/freesasa-zig && zig build test 2>&1 | tail -10
```

Expected: All tests pass.

- [ ] **Step 8: Commit (if any fixes were needed)**

```bash
git add -A && git commit -m "fix: integration test fixes for SDF support"
```

Only if changes were needed. Skip if all tests passed on first try.

---

### Task 11: Documentation Update

**Files:**
- Modify: `docs/` (CLI documentation if it exists)

- [ ] **Step 1: Check for existing CLI docs**

```bash
ls docs/*.md 2>/dev/null
```

If CLI docs exist, update them to mention SDF support and the `--sdf` option. If not, skip this step — the `--help` text is sufficient.

- [ ] **Step 2: Commit doc changes**

```bash
git add docs/ && git commit -m "docs: add SDF file support documentation"
```

Only if documentation files were updated.

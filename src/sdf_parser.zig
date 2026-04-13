//! SDF/MOL Parser for extracting molecular structures.
//!
//! This module parses SDF (Structure-Data File) and MOL file formats,
//! extracting atom coordinates, elements, and bond connectivity.
//!
//! ## Supported Formats
//!
//! - V2000 MOL/SDF (fully supported)
//! - V3000 MOL/SDF (stub — returns `error.InvalidV3000`)
//!
//! ## SDF V2000 Record Format (Fixed Width)
//!
//! - Line 1:   Molecule name
//! - Line 2:   Program/timestamp line (skipped)
//! - Line 3:   Comment line (skipped)
//! - Line 4:   Counts line: cols 0-2 = atom count, cols 3-5 = bond count
//! - Atom block: cols 0-9 = x, 10-19 = y, 20-29 = z, 31-33 = element symbol
//! - Bond block: cols 0-2 = atom1 (1-based), 3-5 = atom2, 6-8 = bond type
//! - `M  END` terminates the connection table
//! - `$$$$` separates molecules in an SDF file
//!
//! ## Usage
//!
//! ```zig
//! const sdf = @import("sdf_parser.zig");
//! const molecules = try sdf.parse(allocator, sdf_data);
//! defer sdf.freeMolecules(allocator, molecules);
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const hybridization = @import("hybridization.zig");

// =============================================================================
// Types
// =============================================================================

/// A single atom parsed from an SDF/MOL file.
pub const SdfAtom = struct {
    x: f64,
    y: f64,
    z: f64,
    element: elem.Element,
};

/// A single bond parsed from an SDF/MOL file.
pub const SdfBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: hybridization.BondOrder,
};

/// A parsed molecule from an SDF/MOL file.
pub const SdfMolecule = struct {
    name: []const u8,
    atoms: []const SdfAtom,
    bonds: []const SdfBond,
};

/// Error types for SDF parsing.
pub const SdfError = error{
    /// Source data is empty.
    EmptySdf,
    /// Bond references an atom index out of range.
    BondIndexOutOfRange,
    /// Atom line is too short or malformed.
    InvalidAtomLine,
    /// Bond line is too short or malformed.
    InvalidBondLine,
    /// Counts line is too short or malformed.
    InvalidCountsLine,
    /// V3000 format detected (not yet supported).
    InvalidV3000,
    /// Float parsing failed.
    InvalidFloat,
    /// Integer parsing failed.
    InvalidInteger,
    /// Memory allocation failed.
    OutOfMemory,
};

// =============================================================================
// Public API
// =============================================================================

/// Parse an SDF or MOL file, returning all molecules found.
///
/// A single MOL file (no `$$$$` terminator) is treated as a one-molecule SDF.
/// Caller must free the result with `freeMolecules`.
pub fn parse(allocator: Allocator, source: []const u8) SdfError![]const SdfMolecule {
    if (source.len == 0) return error.EmptySdf;

    var molecules = std.ArrayListUnmanaged(SdfMolecule){};
    errdefer {
        for (molecules.items) |mol| {
            allocator.free(mol.atoms);
            allocator.free(mol.bonds);
            allocator.free(mol.name);
        }
        molecules.deinit(allocator);
    }

    var line_iter = std.mem.splitScalar(u8, source, '\n');
    var has_content = false;

    while (true) {
        const mol = (try parseSingleMolecule(allocator, &line_iter, &has_content)) orelse break;
        try molecules.append(allocator, mol.molecule);
        if (!mol.has_terminator) break;
        has_content = false;
    }

    if (molecules.items.len == 0) return error.EmptySdf;
    return try molecules.toOwnedSlice(allocator);
}

/// Free all molecules returned by `parse`.
pub fn freeMolecules(allocator: Allocator, molecules: []const SdfMolecule) void {
    for (molecules) |mol| {
        allocator.free(mol.atoms);
        allocator.free(mol.bonds);
        allocator.free(mol.name);
    }
    allocator.free(molecules);
}

// =============================================================================
// Internal Parsing Helpers
// =============================================================================

const MoleculeResult = struct {
    molecule: SdfMolecule,
    has_terminator: bool,
};

/// Parse a single molecule from the line iterator.
/// Returns `null` when there are no more molecules (EOF).
fn parseSingleMolecule(
    allocator: Allocator,
    line_iter: *std.mem.SplitIterator(u8, .scalar),
    has_content: *bool,
) SdfError!?MoleculeResult {
    // Find the first non-blank line (molecule name)
    while (true) {
        const first_line_raw = line_iter.next() orelse return null;
        const first_line = stripCr(first_line_raw);

        if (first_line.len == 0 and !has_content.*) continue;
        has_content.* = true;

        // Line 1 = molecule name
        const name = try allocator.dupe(u8, std.mem.trim(u8, first_line, " \t"));
        errdefer allocator.free(name);

        // Line 2 = program/timestamp (skip)
        _ = line_iter.next() orelse {
            allocator.free(name);
            return null;
        };
        // Line 3 = comment (skip)
        _ = line_iter.next() orelse {
            allocator.free(name);
            return null;
        };

        // Line 4 = counts line
        const counts_raw = line_iter.next() orelse {
            allocator.free(name);
            return null;
        };
        const counts_line = stripCr(counts_raw);

        // Check for V3000
        if (std.mem.indexOf(u8, counts_line, "V3000") != null) {
            return error.InvalidV3000;
        }

        const counts = parseCounts(counts_line) orelse return error.InvalidCountsLine;

        // Parse atom block
        var atom_list = std.ArrayListUnmanaged(SdfAtom){};
        errdefer atom_list.deinit(allocator);
        try atom_list.ensureTotalCapacity(allocator, counts.atom_count);

        for (0..counts.atom_count) |_| {
            const atom_raw = line_iter.next() orelse return error.InvalidAtomLine;
            const atom_line = stripCr(atom_raw);
            const atom = try parseAtomLine(atom_line);
            atom_list.appendAssumeCapacity(atom);
        }

        // Parse bond block
        var bond_list = std.ArrayListUnmanaged(SdfBond){};
        errdefer bond_list.deinit(allocator);
        try bond_list.ensureTotalCapacity(allocator, counts.bond_count);

        for (0..counts.bond_count) |_| {
            const bond_raw = line_iter.next() orelse return error.InvalidBondLine;
            const bond_line = stripCr(bond_raw);
            const bond = try parseBondLine(bond_line, counts.atom_count);
            bond_list.appendAssumeCapacity(bond);
        }

        // Skip remaining lines until $$$$ or EOF
        var found_terminator = false;
        while (line_iter.next()) |rest_raw| {
            const rest_line = stripCr(rest_raw);
            if (std.mem.startsWith(u8, rest_line, "$$$$")) {
                found_terminator = true;
                break;
            }
        }

        const atoms = try atom_list.toOwnedSlice(allocator);
        errdefer allocator.free(atoms);
        const bonds = try bond_list.toOwnedSlice(allocator);

        return .{
            .molecule = .{
                .name = name,
                .atoms = atoms,
                .bonds = bonds,
            },
            .has_terminator = found_terminator,
        };
    }
}

const Counts = struct {
    atom_count: u16,
    bond_count: u16,
};

/// Parse the V2000 counts line: first 3 chars = atom count, next 3 = bond count.
fn parseCounts(line: []const u8) ?Counts {
    if (line.len < 6) return null;
    const atom_count = parseFixedInt(u16, line[0..3]) orelse return null;
    const bond_count = parseFixedInt(u16, line[3..6]) orelse return null;
    return .{ .atom_count = atom_count, .bond_count = bond_count };
}

/// Parse a V2000 atom line.
/// Columns: 0-9=x, 10-19=y, 20-29=z, 31-33=element symbol.
fn parseAtomLine(line: []const u8) SdfError!SdfAtom {
    if (line.len < 34) return error.InvalidAtomLine;

    const x = parseFixedFloat(line[0..10]) orelse return error.InvalidFloat;
    const y = parseFixedFloat(line[10..20]) orelse return error.InvalidFloat;
    const z = parseFixedFloat(line[20..30]) orelse return error.InvalidFloat;

    // Element symbol at columns 31-33 (0-indexed), trimmed
    const element_str = std.mem.trim(u8, line[31..34], " ");
    const element = elem.fromSymbol(element_str);

    return .{ .x = x, .y = y, .z = z, .element = element };
}

/// Parse a V2000 bond line.
/// Columns: 0-2=atom1 (1-based), 3-5=atom2 (1-based), 6-8=bond type.
fn parseBondLine(line: []const u8, atom_count: u16) SdfError!SdfBond {
    if (line.len < 9) return error.InvalidBondLine;

    const idx1_raw = parseFixedInt(u16, line[0..3]) orelse return error.InvalidInteger;
    const idx2_raw = parseFixedInt(u16, line[3..6]) orelse return error.InvalidInteger;
    const bond_type = parseFixedInt(u8, line[6..9]) orelse return error.InvalidInteger;

    // Validate: 1-based indices must be within [1, atom_count]
    if (idx1_raw == 0 or idx1_raw > atom_count) return error.BondIndexOutOfRange;
    if (idx2_raw == 0 or idx2_raw > atom_count) return error.BondIndexOutOfRange;

    // Convert to 0-based
    return .{
        .atom_idx_1 = idx1_raw - 1,
        .atom_idx_2 = idx2_raw - 1,
        .order = sdfBondOrder(bond_type),
    };
}

/// Map an SDF bond-type integer to a BondOrder.
/// Reusable for both V2000 and (future) V3000 parsers.
fn sdfBondOrder(bond_type: u8) hybridization.BondOrder {
    return switch (bond_type) {
        1 => .single,
        2 => .double,
        3 => .triple,
        4 => .aromatic,
        else => .unknown,
    };
}

/// Parse a fixed-width integer field (trimming whitespace).
fn parseFixedInt(comptime T: type, field: []const u8) ?T {
    const trimmed = std.mem.trim(u8, field, " ");
    if (trimmed.len == 0) return null;
    return std.fmt.parseInt(T, trimmed, 10) catch null;
}

/// Parse a fixed-width float field (trimming whitespace).
fn parseFixedFloat(field: []const u8) ?f64 {
    const trimmed = std.mem.trim(u8, field, " ");
    if (trimmed.len == 0) return null;
    return std.fmt.parseFloat(f64, trimmed) catch null;
}

/// Strip trailing carriage return for Windows line endings.
fn stripCr(line: []const u8) []const u8 {
    if (line.len > 0 and line[line.len - 1] == '\r') {
        return line[0 .. line.len - 1];
    }
    return line;
}

// =============================================================================
// Tests
// =============================================================================

test "parse V2000 single molecule — ethanol" {
    const allocator = std.testing.allocator;
    const source =
        \\ethanol
        \\     zsasa   3D
        \\
        \\  9  8  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    1.5200    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    2.0800    1.2124    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.5200    0.9400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.5200   -0.5100    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.5200   -0.5100   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    1.8800   -0.5100    0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    1.8800   -0.5100   -0.8900 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    2.9200    1.2124    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  2  1  0  0  0  0
        \\  1  4  1  0  0  0  0
        \\  1  5  1  0  0  0  0
        \\  1  6  1  0  0  0  0
        \\  2  3  1  0  0  0  0
        \\  2  7  1  0  0  0  0
        \\  2  8  1  0  0  0  0
        \\  3  9  1  0  0  0  0
        \\M  END
        \\$$$$
    ;
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
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), mol.atoms[0].y, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), mol.atoms[0].z, 0.001);

    // Second atom: C at (1.52, 0, 0)
    try std.testing.expectEqual(elem.Element.C, mol.atoms[1].element);
    try std.testing.expectApproxEqAbs(@as(f64, 1.52), mol.atoms[1].x, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), mol.atoms[1].y, 0.001);

    // Third atom: O
    try std.testing.expectEqual(elem.Element.O, mol.atoms[2].element);

    // Bond 0: atoms 0-1, single
    try std.testing.expectEqual(@as(u16, 0), mol.bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), mol.bonds[0].atom_idx_2);
    try std.testing.expectEqual(hybridization.BondOrder.single, mol.bonds[0].order);
}

test "parse V2000 multi-molecule SDF" {
    const allocator = std.testing.allocator;
    const source =
        \\methane
        \\     zsasa   3D
        \\
        \\  5  4  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.6300    0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.6300   -0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.6300    0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.6300   -0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  2  1  0  0  0  0
        \\  1  3  1  0  0  0  0
        \\  1  4  1  0  0  0  0
        \\  1  5  1  0  0  0  0
        \\M  END
        \\$$$$
        \\water
        \\     zsasa   3D
        \\
        \\  3  2  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.7572    0.5858    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.7572    0.5858    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  2  1  0  0  0  0
        \\  1  3  1  0  0  0  0
        \\M  END
        \\$$$$
    ;
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 2), molecules.len);

    // Methane: 5 atoms, 4 bonds
    try std.testing.expectEqualStrings("methane", molecules[0].name);
    try std.testing.expectEqual(@as(usize, 5), molecules[0].atoms.len);
    try std.testing.expectEqual(@as(usize, 4), molecules[0].bonds.len);

    // Water: 3 atoms, 2 bonds
    try std.testing.expectEqualStrings("water", molecules[1].name);
    try std.testing.expectEqual(@as(usize, 3), molecules[1].atoms.len);
    try std.testing.expectEqual(@as(usize, 2), molecules[1].bonds.len);
}

test "parse MOL (no $$$$ terminator)" {
    const allocator = std.testing.allocator;
    const source =
        \\methane
        \\     zsasa   3D
        \\
        \\  5  4  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.6300    0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.6300   -0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   -0.6300    0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    0.6300   -0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  2  1  0  0  0  0
        \\  1  3  1  0  0  0  0
        \\  1  4  1  0  0  0  0
        \\  1  5  1  0  0  0  0
        \\M  END
    ;
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    try std.testing.expectEqualStrings("methane", molecules[0].name);
    try std.testing.expectEqual(@as(usize, 5), molecules[0].atoms.len);
    try std.testing.expectEqual(@as(usize, 4), molecules[0].bonds.len);
}

test "parse empty SDF returns error" {
    const allocator = std.testing.allocator;
    const result = parse(allocator, "");
    try std.testing.expectError(error.EmptySdf, result);
}

test "parse SDF with bad bond index returns error" {
    const allocator = std.testing.allocator;
    // Bond references atom 5, but only 2 atoms exist.
    const source =
        \\bad_bonds
        \\     zsasa   3D
        \\
        \\  2  1  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    1.5200    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\  1  5  1  0  0  0  0
        \\M  END
    ;
    const result = parse(allocator, source);
    try std.testing.expectError(error.BondIndexOutOfRange, result);
}

test "parse V3000 returns InvalidV3000 error" {
    const allocator = std.testing.allocator;
    const source =
        \\ethanol
        \\     zsasa   3D
        \\
        \\  0  0  0  0  0  0  0  0  0  0999 V3000
        \\M  V30 BEGIN CTAB
        \\M  V30 COUNTS 9 8 0 0 0
        \\M  V30 BEGIN ATOM
        \\M  V30 1 C 0.0000 0.0000 0.0000 0
        \\M  V30 END ATOM
        \\M  V30 END CTAB
        \\M  END
        \\$$$$
    ;
    const result = parse(allocator, source);
    try std.testing.expectError(error.InvalidV3000, result);
}

test "parse SDF with CRLF line endings" {
    const allocator = std.testing.allocator;
    const source = "methane\r\n     zsasa   3D\r\n\r\n  5  4  0  0  0  0  0  0  0  0999 V2000\r\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\r\n    0.6300    0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -0.6300   -0.6300    0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0\r\n   -0.6300    0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0\r\n    0.6300   -0.6300   -0.6300 H   0  0  0  0  0  0  0  0  0  0  0  0\r\n  1  2  1  0  0  0  0\r\n  1  3  1  0  0  0  0\r\n  1  4  1  0  0  0  0\r\n  1  5  1  0  0  0  0\r\nM  END\r\n$$$$\r\n";
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    try std.testing.expectEqualStrings("methane", molecules[0].name);
    try std.testing.expectEqual(@as(usize, 5), molecules[0].atoms.len);
}

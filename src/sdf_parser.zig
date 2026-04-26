//! SDF/MOL Parser for extracting molecular structures.
//!
//! This module parses SDF (Structure-Data File) and MOL file formats,
//! extracting atom coordinates, elements, and bond connectivity.
//!
//! ## Supported Formats
//!
//! - V2000 MOL/SDF (fully supported)
//! - V3000 MOL/SDF (fully supported)
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
const ccd_parser = @import("ccd_parser.zig");
const gzip = @import("gzip.zig");
const types = @import("types.zig");

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
    /// V3000 format is malformed.
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

    var molecules = std.ArrayListUnmanaged(SdfMolecule).empty;
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

        // Check for V3000 — parseV3000Body takes ownership of `name`
        // (handles cleanup on both success and error), so we must NOT
        // free `name` here on this path.
        if (std.mem.find(u8, counts_line, "V3000") != null) {
            return try parseV3000Body(allocator, name, line_iter);
        }

        // From this point, `name` is our responsibility on error paths.
        errdefer allocator.free(name);

        const counts = parseCounts(counts_line) orelse return error.InvalidCountsLine;

        // Parse atom block
        var atom_list = std.ArrayListUnmanaged(SdfAtom).empty;
        errdefer atom_list.deinit(allocator);
        try atom_list.ensureTotalCapacity(allocator, counts.atom_count);

        for (0..counts.atom_count) |_| {
            const atom_raw = line_iter.next() orelse return error.InvalidAtomLine;
            const atom_line = stripCr(atom_raw);
            const atom = try parseAtomLine(atom_line);
            atom_list.appendAssumeCapacity(atom);
        }

        // Parse bond block
        var bond_list = std.ArrayListUnmanaged(SdfBond).empty;
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

/// Parse a V3000 molecule body.
/// The line iterator is positioned just after the counts line (which contained "V3000").
/// We expect `M  V30 BEGIN CTAB`, `M  V30 COUNTS ...`, ATOM/BOND blocks, and `M  END`.
fn parseV3000Body(
    allocator: Allocator,
    name: []const u8,
    line_iter: *std.mem.SplitIterator(u8, .scalar),
) SdfError!MoleculeResult {
    errdefer allocator.free(name);

    // Helper: strip "M  V30 " prefix from a line and return the payload, or null.
    const stripV30 = struct {
        fn strip(line: []const u8) ?[]const u8 {
            const trimmed = std.mem.trimStart(u8, line, " ");
            if (std.mem.startsWith(u8, trimmed, "M  V30 ")) {
                return trimmed["M  V30 ".len..];
            }
            return null;
        }
    }.strip;

    // Read lines until we get COUNTS
    var atom_count: u16 = 0;
    var bond_count: u16 = 0;
    var found_counts = false;

    while (line_iter.next()) |raw_line| {
        const line = stripCr(raw_line);
        if (std.mem.startsWith(u8, std.mem.trimStart(u8, line, " "), "M  END")) break;

        if (stripV30(line)) |payload| {
            if (std.mem.startsWith(u8, payload, "COUNTS ")) {
                // Parse "COUNTS natoms nbonds ..."
                var tok = std.mem.tokenizeScalar(u8, payload, ' ');
                _ = tok.next(); // skip "COUNTS"
                const na = tok.next() orelse return error.InvalidCountsLine;
                const nb = tok.next() orelse return error.InvalidCountsLine;
                atom_count = std.fmt.parseInt(u16, na, 10) catch return error.InvalidInteger;
                bond_count = std.fmt.parseInt(u16, nb, 10) catch return error.InvalidInteger;
                found_counts = true;
            } else if (std.mem.startsWith(u8, payload, "BEGIN ATOM")) {
                break; // proceed to atom parsing
            }
        }
    }

    if (!found_counts) return error.InvalidCountsLine;

    // Parse ATOM block
    var atom_list = std.ArrayListUnmanaged(SdfAtom).empty;
    errdefer atom_list.deinit(allocator);
    try atom_list.ensureTotalCapacity(allocator, atom_count);

    while (line_iter.next()) |raw_line| {
        const line = stripCr(raw_line);
        if (stripV30(line)) |payload| {
            if (std.mem.startsWith(u8, payload, "END ATOM")) break;

            // "index element x y z charge [...]"
            var tok = std.mem.tokenizeScalar(u8, payload, ' ');
            _ = tok.next() orelse return error.InvalidAtomLine; // index (skip)
            const elem_str = tok.next() orelse return error.InvalidAtomLine;
            const x_str = tok.next() orelse return error.InvalidAtomLine;
            const y_str = tok.next() orelse return error.InvalidAtomLine;
            const z_str = tok.next() orelse return error.InvalidAtomLine;
            // charge and remaining fields are ignored

            const x = std.fmt.parseFloat(f64, x_str) catch return error.InvalidFloat;
            const y = std.fmt.parseFloat(f64, y_str) catch return error.InvalidFloat;
            const z = std.fmt.parseFloat(f64, z_str) catch return error.InvalidFloat;
            const element = elem.fromSymbol(elem_str);

            atom_list.appendAssumeCapacity(.{ .x = x, .y = y, .z = z, .element = element });
        }
    }

    // Look for BEGIN BOND (there may be lines between END ATOM and BEGIN BOND)
    var bond_list = std.ArrayListUnmanaged(SdfBond).empty;
    errdefer bond_list.deinit(allocator);
    try bond_list.ensureTotalCapacity(allocator, bond_count);

    // We may have already consumed "BEGIN ATOM"... now scan for "BEGIN BOND"
    var in_bond_block = false;
    while (line_iter.next()) |raw_line| {
        const line = stripCr(raw_line);
        if (std.mem.startsWith(u8, std.mem.trimStart(u8, line, " "), "M  END")) {
            // Reached M  END — no bond block or we're done
            break;
        }
        if (std.mem.startsWith(u8, std.mem.trimStart(u8, line, " "), "$$$$")) {
            // Terminator found before M  END
            const atoms = try atom_list.toOwnedSlice(allocator);
            errdefer allocator.free(atoms);
            const bonds = try bond_list.toOwnedSlice(allocator);

            return .{
                .molecule = .{ .name = name, .atoms = atoms, .bonds = bonds },
                .has_terminator = true,
            };
        }

        if (stripV30(line)) |payload| {
            if (std.mem.startsWith(u8, payload, "BEGIN BOND")) {
                in_bond_block = true;
                continue;
            }
            if (std.mem.startsWith(u8, payload, "END BOND")) {
                in_bond_block = false;
                continue;
            }
            if (std.mem.startsWith(u8, payload, "END CTAB")) {
                continue;
            }

            if (in_bond_block) {
                // "index bondtype atom1 atom2 [...]"
                var tok = std.mem.tokenizeScalar(u8, payload, ' ');
                _ = tok.next() orelse return error.InvalidBondLine; // index (skip)
                const bt_str = tok.next() orelse return error.InvalidBondLine;
                const a1_str = tok.next() orelse return error.InvalidBondLine;
                const a2_str = tok.next() orelse return error.InvalidBondLine;

                const bond_type = std.fmt.parseInt(u8, bt_str, 10) catch return error.InvalidInteger;
                const idx1_raw = std.fmt.parseInt(u16, a1_str, 10) catch return error.InvalidInteger;
                const idx2_raw = std.fmt.parseInt(u16, a2_str, 10) catch return error.InvalidInteger;

                if (idx1_raw == 0 or idx1_raw > atom_count) return error.BondIndexOutOfRange;
                if (idx2_raw == 0 or idx2_raw > atom_count) return error.BondIndexOutOfRange;

                bond_list.appendAssumeCapacity(.{
                    .atom_idx_1 = idx1_raw - 1,
                    .atom_idx_2 = idx2_raw - 1,
                    .order = sdfBondOrder(bond_type),
                });
            }
        }
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
        .molecule = .{ .name = name, .atoms = atoms, .bonds = bonds },
        .has_terminator = found_terminator,
    };
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
// Conversion Functions
// =============================================================================

/// Generate atom names like "C1", "C2", "O1" by appending a per-element counter
/// to the element symbol. Writes the result into a FixedString4-compatible [4]u8
/// and returns the length.
fn formatAtomName(symbol_str: []const u8, counter: u16, buf: *[4]u8) u3 {
    // Format counter as decimal (up to 3 digits for 4-char limit)
    var num_buf: [3]u8 = undefined;
    const num_str = std.fmt.bufPrint(&num_buf, "{d}", .{counter}) catch "";

    const total_len = @min(symbol_str.len + num_str.len, 4);
    buf.* = .{ 0, 0, 0, 0 };

    const sym_copy: usize = @min(symbol_str.len, 4);
    for (symbol_str[0..sym_copy], 0..) |c, i| {
        buf[i] = c;
    }
    const remaining = 4 - sym_copy;
    const num_copy = @min(num_str.len, remaining);
    for (num_str[0..num_copy], 0..) |c, i| {
        buf[sym_copy + i] = c;
    }

    return @intCast(total_len);
}

/// Convert an SdfMolecule to a StoredComponent for the CCD classifier.
///
/// - `comp_id` = molecule name truncated to 5 chars
/// - Atom names are generated as element symbol + per-element counter (C1, C2, O1...)
/// - Bond indices and orders are preserved from the SDF data
/// - Caller must call `.deinit()` on the returned StoredComponent.
pub fn toStoredComponent(allocator: Allocator, molecule: *const SdfMolecule) !ccd_parser.StoredComponent {
    // Build per-element counters for atom naming
    const atoms = try allocator.alloc(hybridization.CompAtom, molecule.atoms.len);
    errdefer allocator.free(atoms);

    // Element counter: keyed by element enum value
    var element_counts: [119]u16 = .{0} ** 119;

    for (molecule.atoms, 0..) |sdf_atom, i| {
        const sym = sdf_atom.element.symbol();
        const elem_idx = sdf_atom.element.atomicNumber();
        element_counts[elem_idx] += 1;

        var name_buf: [4]u8 = undefined;
        const name_len = formatAtomName(sym, element_counts[elem_idx], &name_buf);

        atoms[i] = hybridization.CompAtom{
            .atom_id = name_buf,
            .atom_id_len = name_len,
            .type_symbol = .{ 0, 0, 0, 0 },
            .type_symbol_len = 0,
            .aromatic = false,
            .leaving = false,
        };

        // Copy type_symbol from element symbol
        const ts_len: usize = @min(sym.len, 4);
        atoms[i].type_symbol_len = @intCast(ts_len);
        for (sym[0..ts_len], 0..) |c, j| {
            atoms[i].type_symbol[j] = c;
        }
    }

    const bonds = try allocator.alloc(hybridization.CompBond, molecule.bonds.len);
    errdefer allocator.free(bonds);

    for (molecule.bonds, 0..) |sdf_bond, i| {
        bonds[i] = .{
            .atom_idx_1 = sdf_bond.atom_idx_1,
            .atom_idx_2 = sdf_bond.atom_idx_2,
            .order = sdf_bond.order,
            .aromatic = sdf_bond.order == .aromatic,
        };
    }

    // Build comp_id from molecule name (truncated to 5, lowercased for consistency)
    var comp_id: [5]u8 = .{ 0, 0, 0, 0, 0 };
    const name_len: usize = @min(molecule.name.len, 5);
    for (molecule.name[0..name_len], 0..) |c, i| {
        comp_id[i] = c;
    }

    return .{
        .comp_id = comp_id,
        .comp_id_len = @intCast(name_len),
        .atoms = atoms,
        .bonds = bonds,
        .allocator = allocator,
    };
}

/// Convert SDF molecules into AtomInput for SASA calculation.
///
/// - Each molecule becomes one chain (A, B, C... up to Z, max 26)
/// - Residue name = molecule name truncated to 5 chars
/// - Atom names generated as element + per-element index (C1, C2, O1...)
/// - Radii default to element VdW radius (classifier will override later)
/// - When `skip_hydrogens` is true, H atoms are excluded
pub fn toAtomInput(allocator: Allocator, molecules: []const SdfMolecule, skip_hydrogens: bool) !types.AtomInput {
    // Limit to 26 chains (A-Z)
    const max_chains: usize = @min(molecules.len, 26);

    // First pass: count total atoms (only for molecules we will process)
    var total_atoms: usize = 0;
    for (molecules[0..max_chains]) |mol| {
        for (mol.atoms) |atom| {
            if (skip_hydrogens and atom.element == .H) continue;
            total_atoms += 1;
        }
    }

    // Allocate arrays
    const x = try allocator.alloc(f64, total_atoms);
    errdefer allocator.free(x);
    const y = try allocator.alloc(f64, total_atoms);
    errdefer allocator.free(y);
    const z = try allocator.alloc(f64, total_atoms);
    errdefer allocator.free(z);
    const r = try allocator.alloc(f64, total_atoms);
    errdefer allocator.free(r);
    const residue = try allocator.alloc(types.FixedString5, total_atoms);
    errdefer allocator.free(residue);
    const atom_name = try allocator.alloc(types.FixedString4, total_atoms);
    errdefer allocator.free(atom_name);
    const element_arr = try allocator.alloc(u8, total_atoms);
    errdefer allocator.free(element_arr);
    const chain_id = try allocator.alloc(types.FixedString4, total_atoms);
    errdefer allocator.free(chain_id);
    const residue_num = try allocator.alloc(i32, total_atoms);
    errdefer allocator.free(residue_num);
    const insertion_code = try allocator.alloc(types.FixedString4, total_atoms);
    errdefer allocator.free(insertion_code);

    var idx: usize = 0;

    for (molecules[0..max_chains], 0..) |mol, mol_idx| {
        const chain_letter: u8 = 'A' + @as(u8, @intCast(mol_idx));
        const chain = types.FixedString4.fromSlice(&[_]u8{chain_letter});
        const res_name = types.FixedString5.fromSlice(mol.name[0..@min(mol.name.len, 5)]);
        const empty_insertion = types.FixedString4.fromSlice("");

        // Per-element counters for atom naming (reset per molecule)
        var element_counts: [119]u16 = .{0} ** 119;

        for (mol.atoms) |atom| {
            if (skip_hydrogens and atom.element == .H) continue;

            const sym = atom.element.symbol();
            const elem_idx = atom.element.atomicNumber();
            element_counts[elem_idx] += 1;

            var name_buf: [4]u8 = undefined;
            const name_len = formatAtomName(sym, element_counts[elem_idx], &name_buf);

            x[idx] = atom.x;
            y[idx] = atom.y;
            z[idx] = atom.z;
            r[idx] = atom.element.vdwRadius();
            residue[idx] = res_name;
            atom_name[idx] = .{
                .data = name_buf,
                .len = name_len,
            };
            element_arr[idx] = atom.element.atomicNumber();
            chain_id[idx] = chain;
            residue_num[idx] = 1;
            insertion_code[idx] = empty_insertion;

            idx += 1;
        }
    }

    return .{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .residue = residue,
        .atom_name = atom_name,
        .element = element_arr,
        .chain_id = chain_id,
        .residue_num = residue_num,
        .insertion_code = insertion_code,
        .allocator = allocator,
    };
}

// =============================================================================
// SDF Path List and Component Loading
// =============================================================================

/// Fixed-capacity list for SDF paths (max 16).
pub const SdfPathList = struct {
    items: [max_sdf_paths][]const u8 = undefined,
    len: usize = 0,

    const max_sdf_paths = 16;

    pub fn append(self: *SdfPathList, value: []const u8) error{Overflow}!void {
        if (self.len >= max_sdf_paths) return error.Overflow;
        self.items[self.len] = value;
        self.len += 1;
    }

    pub fn constSlice(self: *const SdfPathList) []const []const u8 {
        return self.items[0..self.len];
    }
};

/// Load SDF files and build a ComponentDict from their bond topology.
///
/// Reads each SDF file (plain or gzip-compressed), parses molecules, and
/// converts them to StoredComponents. Duplicate molecule names (truncated
/// to 5 chars) are skipped to avoid leaking the first entry.
///
/// Returns `null` if no valid components were loaded.
pub fn loadSdfComponents(
    allocator: Allocator,
    io: std.Io,
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
            const f = std.Io.Dir.cwd().openFile(io, sdf_path, .{}) catch |err| {
                std.debug.print("Error opening SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            };
            defer f.close(io);
            var read_buf: [65536]u8 = undefined;
            var file_reader = f.reader(io, &read_buf);
            break :file_blk file_reader.interface.allocRemaining(allocator, .unlimited) catch |err| {
                std.debug.print("Error reading SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            };
        };
        defer allocator.free(source);

        const molecules = parse(allocator, source) catch |err| {
            std.debug.print("Error parsing SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
            std.process.exit(1);
        };
        defer freeMolecules(allocator, molecules);

        for (molecules) |mol| {
            if (mol.name.len == 0) {
                if (!quiet) std.debug.print("Warning: SDF molecule has no name, skipping\n", .{});
                continue;
            }
            const stored = toStoredComponent(allocator, &mol) catch |err| {
                if (!quiet) std.debug.print("Warning: Could not convert SDF molecule '{s}': {s}\n", .{ mol.name, @errorName(err) });
                continue;
            };
            const comp_id_str = mol.name[0..@min(mol.name.len, 5)];

            // Skip if already registered (avoid StoredComponent leak from duplicate names)
            if (dict.components.get(comp_id_str) != null) {
                var s = stored;
                s.deinit();
                continue;
            }

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
            dict.components.put(dict_key, stored) catch {
                // owned_keys already has dict_key; it will be freed by dict.deinit().
                // But we must free the stored component since it was not successfully
                // placed into the map.
                var s = stored;
                s.deinit();
                continue;
            };
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

test "parse V3000 single molecule — ethanol" {
    const allocator = std.testing.allocator;
    const source =
        \\ethanol
        \\     RDKit          3D
        \\
        \\  0  0  0  0  0  0  0  0  0  0999 V3000
        \\M  V30 BEGIN CTAB
        \\M  V30 COUNTS 9 8 0 0 0
        \\M  V30 BEGIN ATOM
        \\M  V30 1 C 0.0000 0.0000 0.0000 0
        \\M  V30 2 C 1.5200 0.0000 0.0000 0
        \\M  V30 3 O 2.0800 1.2100 0.0000 0
        \\M  V30 4 H -0.3900 0.9800 -0.2600 0
        \\M  V30 5 H -0.3900 -0.5400 0.8700 0
        \\M  V30 6 H -0.3900 -0.4400 -0.9200 0
        \\M  V30 7 H 1.9100 -0.5400 0.8700 0
        \\M  V30 8 H 1.9100 0.5400 -0.8700 0
        \\M  V30 9 H 3.0400 1.2100 0.0000 0
        \\M  V30 END ATOM
        \\M  V30 BEGIN BOND
        \\M  V30 1 1 1 2
        \\M  V30 2 1 1 4
        \\M  V30 3 1 1 5
        \\M  V30 4 1 1 6
        \\M  V30 5 1 2 3
        \\M  V30 6 1 2 7
        \\M  V30 7 1 2 8
        \\M  V30 8 1 3 9
        \\M  V30 END BOND
        \\M  V30 END CTAB
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

    // Same coordinates as V2000
    try std.testing.expectEqual(elem.Element.C, mol.atoms[0].element);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), mol.atoms[0].x, 0.001);
    try std.testing.expectEqual(elem.Element.C, mol.atoms[1].element);
    try std.testing.expectApproxEqAbs(@as(f64, 1.52), mol.atoms[1].x, 0.001);
    try std.testing.expectEqual(elem.Element.O, mol.atoms[2].element);

    // Bond: atom 0-1, single
    try std.testing.expectEqual(@as(u16, 0), mol.bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), mol.bonds[0].atom_idx_2);
    try std.testing.expectEqual(hybridization.BondOrder.single, mol.bonds[0].order);
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

test "toStoredComponent — ethanol molecule" {
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

    var stored = try toStoredComponent(allocator, &molecules[0]);
    defer stored.deinit();

    // comp_id = "ethan" (truncated to 5)
    const view = stored.view();
    try std.testing.expectEqualStrings("ethan", view.compIdSlice());

    // 9 atoms, 8 bonds
    try std.testing.expectEqual(@as(usize, 9), stored.atoms.len);
    try std.testing.expectEqual(@as(usize, 8), stored.bonds.len);

    // First atom: type_symbol = "C"
    try std.testing.expectEqualStrings("C", stored.atoms[0].typeSymbolSlice());
    // Third atom: type_symbol = "O"
    try std.testing.expectEqualStrings("O", stored.atoms[2].typeSymbolSlice());

    // Atom names: C1, C2, O1, H1, H2, ...
    try std.testing.expectEqualStrings("C1", stored.atoms[0].atomIdSlice());
    try std.testing.expectEqualStrings("C2", stored.atoms[1].atomIdSlice());
    try std.testing.expectEqualStrings("O1", stored.atoms[2].atomIdSlice());
    try std.testing.expectEqualStrings("H1", stored.atoms[3].atomIdSlice());
}

test "toAtomInput — two molecules get separate chains" {
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

    var input = try toAtomInput(allocator, molecules, false);
    defer input.deinit();

    // methane(5) + water(3) = 8 atoms total
    try std.testing.expectEqual(@as(usize, 8), input.atomCount());

    // Chain IDs: methane atoms = "A", water atoms = "B"
    const chains = input.chain_id.?;
    try std.testing.expectEqualStrings("A", chains[0].slice());
    try std.testing.expectEqualStrings("A", chains[4].slice());
    try std.testing.expectEqualStrings("B", chains[5].slice());
    try std.testing.expectEqualStrings("B", chains[7].slice());

    // Residue names: "metha" (truncated from "methane"), "water"
    const residues = input.residue.?;
    try std.testing.expectEqualStrings("metha", residues[0].slice());
    try std.testing.expectEqualStrings("water", residues[5].slice());

    // Residue numbers = 1
    const res_nums = input.residue_num.?;
    try std.testing.expectEqual(@as(i32, 1), res_nums[0]);
    try std.testing.expectEqual(@as(i32, 1), res_nums[5]);

    // Insertion codes are empty
    const ins_codes = input.insertion_code.?;
    try std.testing.expectEqualStrings("", ins_codes[0].slice());
}

test "toAtomInput — single molecule from multi-molecule SDF" {
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

    // Process only the first molecule (methane) as a slice of 1
    {
        var input = try toAtomInput(allocator, molecules[0..1], false);
        defer input.deinit();

        // methane has 5 atoms (1 C + 4 H)
        try std.testing.expectEqual(@as(usize, 5), input.atomCount());

        // All atoms should be chain "A"
        const chains = input.chain_id.?;
        try std.testing.expectEqualStrings("A", chains[0].slice());
        try std.testing.expectEqualStrings("A", chains[4].slice());

        // Residue name should be "metha" (truncated from "methane")
        try std.testing.expectEqualStrings("metha", input.residue.?[0].slice());
    }

    // Process only the second molecule (water) as a slice of 1
    {
        var input = try toAtomInput(allocator, molecules[1..2], false);
        defer input.deinit();

        // water has 3 atoms (1 O + 2 H)
        try std.testing.expectEqual(@as(usize, 3), input.atomCount());

        // All atoms should be chain "A" (not "B" — it's the first molecule in this slice)
        const chains = input.chain_id.?;
        try std.testing.expectEqualStrings("A", chains[0].slice());
        try std.testing.expectEqualStrings("A", chains[2].slice());

        // Residue name should be "water"
        try std.testing.expectEqualStrings("water", input.residue.?[0].slice());
    }

    // Process second molecule with skip_hydrogens
    {
        var input = try toAtomInput(allocator, molecules[1..2], true);
        defer input.deinit();

        // water without H: 1 atom (O only)
        try std.testing.expectEqual(@as(usize, 1), input.atomCount());
        try std.testing.expectEqualStrings("A", input.chain_id.?[0].slice());
    }
}

test "toAtomInput — skip hydrogens" {
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

    var input = try toAtomInput(allocator, molecules, true);
    defer input.deinit();

    // Ethanol has 3 heavy atoms (2C + 1O), 6H skipped
    try std.testing.expectEqual(@as(usize, 3), input.atomCount());

    // Elements should be C, C, O
    const elements = input.element.?;
    try std.testing.expectEqual(@as(u8, 6), elements[0]); // C
    try std.testing.expectEqual(@as(u8, 6), elements[1]); // C
    try std.testing.expectEqual(@as(u8, 8), elements[2]); // O

    // Atom names should be C1, C2, O1 (H counter not incremented)
    const names = input.atom_name.?;
    try std.testing.expectEqualStrings("C1", names[0].slice());
    try std.testing.expectEqualStrings("C2", names[1].slice());
    try std.testing.expectEqualStrings("O1", names[2].slice());
}

test "toAtomInput — max 26 chains" {
    const allocator = std.testing.allocator;

    // Build 28 single-atom molecules inline
    var source_buf: [28 * 256]u8 = undefined;
    var w = std.Io.Writer.fixed(&source_buf);
    for (0..28) |i| {
        w.print("mol{d:0>2}\n", .{i}) catch unreachable;
        w.writeAll("     zsasa   3D\n") catch unreachable;
        w.writeAll("\n") catch unreachable;
        w.writeAll("  1  0  0  0  0  0  0  0  0  0999 V2000\n") catch unreachable;
        w.writeAll("    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n") catch unreachable;
        w.writeAll("M  END\n") catch unreachable;
        w.writeAll("$$$$\n") catch unreachable;
    }
    const source = source_buf[0..w.end];

    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 28), molecules.len);

    var input = try toAtomInput(allocator, molecules, false);
    defer input.deinit();

    // Only 26 molecules should be processed (A-Z), not 28
    try std.testing.expectEqual(@as(usize, 26), input.atomCount());

    // Last chain should be 'Z'
    const chains = input.chain_id.?;
    try std.testing.expectEqualStrings("Z", chains[25].slice());
}

test "parse V3000 with bad bond index returns error" {
    const allocator = std.testing.allocator;
    const source =
        \\bad_v3k
        \\     RDKit          3D
        \\
        \\  0  0  0  0  0  0  0  0  0  0999 V3000
        \\M  V30 BEGIN CTAB
        \\M  V30 COUNTS 2 1 0 0 0
        \\M  V30 BEGIN ATOM
        \\M  V30 1 C 0.0000 0.0000 0.0000 0
        \\M  V30 2 C 1.5000 0.0000 0.0000 0
        \\M  V30 END ATOM
        \\M  V30 BEGIN BOND
        \\M  V30 1 1 1 5
        \\M  V30 END BOND
        \\M  V30 END CTAB
        \\M  END
        \\$$$$
    ;
    const result = parse(allocator, source);
    try std.testing.expectError(error.BondIndexOutOfRange, result);
}

test "sdfBondOrder maps all types correctly" {
    try std.testing.expectEqual(hybridization.BondOrder.single, sdfBondOrder(1));
    try std.testing.expectEqual(hybridization.BondOrder.double, sdfBondOrder(2));
    try std.testing.expectEqual(hybridization.BondOrder.triple, sdfBondOrder(3));
    try std.testing.expectEqual(hybridization.BondOrder.aromatic, sdfBondOrder(4));
    try std.testing.expectEqual(hybridization.BondOrder.unknown, sdfBondOrder(5));
    try std.testing.expectEqual(hybridization.BondOrder.unknown, sdfBondOrder(0));
    try std.testing.expectEqual(hybridization.BondOrder.unknown, sdfBondOrder(255));
}

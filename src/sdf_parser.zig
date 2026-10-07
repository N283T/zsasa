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
//! - Line 1:   Molecule name (may be blank)
//! - Line 2:   Program/timestamp line (skipped)
//! - Line 3:   Comment line (skipped)
//! - Line 4:   Counts line: cols 0-2 = atom count, cols 3-5 = bond count
//! - Atom block: cols 0-9 = x, 10-19 = y, 20-29 = z, 31-33 = element symbol
//! - Bond block: cols 0-2 = atom1 (1-based), 3-5 = atom2, 6-8 = bond type
//! - `M  END` terminates the connection table
//! - `$$$$` separates molecules in an SDF file; blank lines after it, before
//!   the next record or the end of the file, are skipped
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
const classifier = @import("classifier.zig");
const ccd_parser = @import("ccd_parser.zig");
const compressed = @import("compressed.zig");
const types = @import("types.zig");

// =============================================================================
// Types
// =============================================================================

/// A single atom parsed from an SDF/MOL file.
pub const SdfAtom = struct {
    x: f64,
    y: f64,
    z: f64,
    /// Deuterium (`D`) and tritium (`T`) are stored as `.H`.
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

    while (true) {
        const mol = (try parseSingleMolecule(allocator, &line_iter)) orelse break;
        try molecules.append(allocator, mol.molecule);
        if (!mol.has_terminator) break;
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
) SdfError!?MoleculeResult {
    // Blank lines may follow a `$$$$` separator or end the file, but a blank
    // line is also a legal title. A blank line is the title of a record when
    // the fourth line from it is a counts line, and is skipped otherwise.
    while (true) {
        const line = stripCr(line_iter.peek() orelse return null);
        if (line.len > 0 or recordStartsAt(line_iter.*)) break;
        _ = line_iter.next();
    }

    // The header block is exactly three lines, followed by the counts line.
    // Line 1 = molecule name
    const title_line = stripCr(line_iter.next() orelse return null);
    // Line 2 = program/timestamp (skip)
    _ = line_iter.next() orelse return null;
    // Line 3 = comment (skip)
    _ = line_iter.next() orelse return null;
    // Line 4 = counts line
    const counts_line = stripCr(line_iter.next() orelse return null);

    const name = try allocator.dupe(u8, std.mem.trim(u8, title_line, " \t"));

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

/// Whether a record starts at the position of `line_iter`, judged by its
/// fourth line, which must be the counts line. Only used for a blank line,
/// to tell a blank title from a stray blank line between records.
fn recordStartsAt(line_iter: std.mem.SplitIterator(u8, .scalar)) bool {
    var probe = line_iter;
    for (0..3) |_| _ = probe.next() orelse return false;
    const counts_line = stripCr(probe.next() orelse return false);

    if (std.mem.find(u8, counts_line, "V3000") != null) return true;
    const counts = parseCounts(counts_line) orelse return false;
    if (std.mem.find(u8, counts_line, "V2000") != null) return true;

    // Without a version tag a title or comment that starts with digits also
    // reads as two counts, so the line after it has to be an atom line.
    if (counts.atom_count == 0) return false;
    const atom_line = stripCr(probe.next() orelse return false);
    _ = parseAtomLine(atom_line) catch return false;
    return true;
}

/// One line of a V3000 connection table, as `V3000Reader` classifies it.
const V3000Line = union(enum) {
    /// Payload of an `M  V30 ` line (the text after the prefix), with its
    /// continuation lines joined. Valid until the next call to `next`.
    v30: []const u8,
    /// `M  END`: the end of the connection table.
    end,
    /// `$$$$`: the end of the record.
    terminator,
    /// Any other line.
    other,
};

/// Reads the lines of a V3000 connection table.
///
/// A V3000 line that ends in `-` continues on the next `M  V30 ` line: the
/// dash is dropped and the payload of the next line follows it directly.
const V3000Reader = struct {
    line_iter: *std.mem.SplitIterator(u8, .scalar),
    /// Holds a line joined from its continuation lines.
    joined: std.ArrayListUnmanaged(u8) = .empty,

    const prefix = "M  V30 ";

    fn deinit(self: *V3000Reader, allocator: Allocator) void {
        self.joined.deinit(allocator);
    }

    /// Returns the next line, or `null` at the end of the input.
    fn next(self: *V3000Reader, allocator: Allocator) SdfError!?V3000Line {
        const line = stripCr(self.line_iter.next() orelse return null);
        const trimmed = std.mem.trimStart(u8, line, " ");
        if (std.mem.startsWith(u8, trimmed, "M  END")) return .end;
        if (std.mem.startsWith(u8, trimmed, "$$$$")) return .terminator;
        if (!std.mem.startsWith(u8, trimmed, prefix)) return .other;

        var payload = trimmed[prefix.len..];
        if (!std.mem.endsWith(u8, payload, "-")) return .{ .v30 = payload };

        self.joined.clearRetainingCapacity();
        while (std.mem.endsWith(u8, payload, "-")) {
            try self.joined.appendSlice(allocator, payload[0 .. payload.len - 1]);
            // The line that continues it has to be a V3000 line as well
            const next_line = stripCr(self.line_iter.next() orelse return error.InvalidV3000);
            const next_trimmed = std.mem.trimStart(u8, next_line, " ");
            if (!std.mem.startsWith(u8, next_trimmed, prefix)) return error.InvalidV3000;
            payload = next_trimmed[prefix.len..];
        }
        try self.joined.appendSlice(allocator, payload);
        return .{ .v30 = self.joined.items };
    }
};

/// Parse a V3000 molecule body.
/// The line iterator is positioned just after the counts line (which contained "V3000").
/// We expect `M  V30 BEGIN CTAB`, `M  V30 COUNTS ...`, ATOM/BOND blocks, and `M  END`.
///
/// Atoms are read only between `BEGIN ATOM` and `END ATOM` and bonds only
/// between `BEGIN BOND` and `END BOND`; the connection table ends at `M  END`
/// (or at `$$$$`, for a record without one). A molecule may have no atom
/// block at all: `COUNTS 0 0` gives a molecule without atoms. Continued lines
/// are joined before they are read, in every block.
fn parseV3000Body(
    allocator: Allocator,
    name: []const u8,
    line_iter: *std.mem.SplitIterator(u8, .scalar),
) SdfError!MoleculeResult {
    errdefer allocator.free(name);

    var reader = V3000Reader{ .line_iter = line_iter };
    defer reader.deinit(allocator);

    // COUNTS only sizes the lists: the atoms and bonds actually present must
    // match it, so a body that disagrees with its header is an error instead
    // of a write past the reserved capacity.
    var counts: ?Counts = null;
    var atom_list = std.ArrayListUnmanaged(SdfAtom).empty;
    errdefer atom_list.deinit(allocator);
    var bond_list = std.ArrayListUnmanaged(SdfBond).empty;
    errdefer bond_list.deinit(allocator);

    var block: enum { none, atom, bond } = .none;
    var found_terminator = false;

    while (try reader.next(allocator)) |line| {
        const payload = switch (line) {
            .v30 => |text| text,
            .end => break,
            .terminator => {
                found_terminator = true;
                break;
            },
            .other => continue,
        };

        switch (block) {
            .none => if (std.mem.startsWith(u8, payload, "COUNTS ")) {
                // Parse "COUNTS natoms nbonds ..."
                var tok = std.mem.tokenizeScalar(u8, payload, ' ');
                _ = tok.next(); // skip "COUNTS"
                const na = tok.next() orelse return error.InvalidCountsLine;
                const nb = tok.next() orelse return error.InvalidCountsLine;
                const atom_count = std.fmt.parseInt(u16, na, 10) catch return error.InvalidInteger;
                const bond_count = std.fmt.parseInt(u16, nb, 10) catch return error.InvalidInteger;
                try atom_list.ensureTotalCapacity(allocator, atom_count);
                try bond_list.ensureTotalCapacity(allocator, bond_count);
                counts = .{ .atom_count = atom_count, .bond_count = bond_count };
            } else if (std.mem.startsWith(u8, payload, "BEGIN ATOM")) {
                if (counts == null) return error.InvalidCountsLine;
                block = .atom;
            } else if (std.mem.startsWith(u8, payload, "BEGIN BOND")) {
                if (counts == null) return error.InvalidCountsLine;
                block = .bond;
            },
            .atom => {
                if (std.mem.startsWith(u8, payload, "END ATOM")) {
                    block = .none;
                    continue;
                }
                if (atom_list.items.len >= counts.?.atom_count) return error.InvalidV3000;

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
                const element = elementFromSymbol(elem_str);

                atom_list.appendAssumeCapacity(.{ .x = x, .y = y, .z = z, .element = element });
            },
            .bond => {
                if (std.mem.startsWith(u8, payload, "END BOND")) {
                    block = .none;
                    continue;
                }
                if (bond_list.items.len >= counts.?.bond_count) return error.InvalidV3000;

                // "index bondtype atom1 atom2 [...]"
                var tok = std.mem.tokenizeScalar(u8, payload, ' ');
                _ = tok.next() orelse return error.InvalidBondLine; // index (skip)
                const bt_str = tok.next() orelse return error.InvalidBondLine;
                const a1_str = tok.next() orelse return error.InvalidBondLine;
                const a2_str = tok.next() orelse return error.InvalidBondLine;

                const bond_type = std.fmt.parseInt(u8, bt_str, 10) catch return error.InvalidInteger;
                const idx1_raw = std.fmt.parseInt(u16, a1_str, 10) catch return error.InvalidInteger;
                const idx2_raw = std.fmt.parseInt(u16, a2_str, 10) catch return error.InvalidInteger;

                // 1-based indices must refer to atoms that were read
                if (idx1_raw == 0 or idx1_raw > atom_list.items.len) return error.BondIndexOutOfRange;
                if (idx2_raw == 0 or idx2_raw > atom_list.items.len) return error.BondIndexOutOfRange;

                bond_list.appendAssumeCapacity(.{
                    .atom_idx_1 = idx1_raw - 1,
                    .atom_idx_2 = idx2_raw - 1,
                    .order = sdfBondOrder(bond_type),
                });
            },
        }
    }

    const declared = counts orelse return error.InvalidCountsLine;
    // A block that is still open was cut off by `M  END`, `$$$$` or the end
    // of the input.
    if (block != .none) return error.InvalidV3000;
    if (atom_list.items.len != declared.atom_count) return error.InvalidV3000;
    if (bond_list.items.len != declared.bond_count) return error.InvalidV3000;

    // Skip remaining lines until $$$$ or EOF
    while (!found_terminator) {
        const rest_line = stripCr(line_iter.next() orelse break);
        if (std.mem.startsWith(u8, rest_line, "$$$$")) found_terminator = true;
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
    const element = elementFromSymbol(element_str);

    return .{ .x = x, .y = y, .z = z, .element = element };
}

/// Element of an atom block symbol.
///
/// `D` (deuterium) and `T` (tritium) are atom symbols of their own in MOL
/// files. They are stored as hydrogen, so that they are excluded together
/// with the other hydrogens and count as hydrogens when radii are derived
/// from the bond table.
fn elementFromSymbol(symbol: []const u8) elem.Element {
    if (symbol.len == 1) {
        switch (std.ascii.toUpper(symbol[0])) {
            'D', 'T' => return .H,
            else => {},
        }
    }
    return elem.fromSymbol(symbol);
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

/// Generates the atom names of one molecule: the element symbol followed by
/// a per-element counter (`C1`, `C2`, `O1`, ...), unique within the molecule.
///
/// A name has at most four characters: that is the size of an atom name in
/// `types.AtomInput` and of `hybridization.CompAtom.atom_id`, and the
/// classifier looks atoms up by the first four characters of their name. The
/// counter therefore has three characters after a one-letter symbol and two
/// after a two-letter symbol. It is written
///
/// 1. in decimal while that fits: `C1` to `C999`, `Cl1` to `Cl99`;
/// 2. then, as in hybrid-36, as a base-36 number of the full width whose
///    first digit is a letter: `CA00` to `CZZZ` for carbons 1,000 to 34,695,
///    `ClA0` to `ClZZ` for chlorines 100 to 1,035. No element symbol has an
///    upper-case second letter, so these are not names of another element;
/// 3. past that, the name has no element symbol: it is a four-digit base-36
///    number whose first digit is a decimal digit (`0000`, `0001`, ...),
///    counted over all such atoms of the molecule. That is enough for every
///    atom of a molecule with up to 466,560 of them, more than the 65,535
///    atoms a bond can refer to.
const AtomNamer = struct {
    /// Atoms named so far, by atomic number.
    element_counts: [119]u16 = @splat(0),
    /// Atoms named without their element symbol (case 3).
    unprefixed_count: u32 = 0,

    const base36_digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    /// Names the next atom of the molecule.
    ///
    /// Every atom has to be passed, in file order, including atoms that are
    /// left out of the result: an atom then has the same name in
    /// `toAtomInput` and in `toStoredComponent`.
    fn next(self: *AtomNamer, element: elem.Element) types.FixedString4 {
        const sym = element.symbol();
        std.debug.assert(sym.len == 1 or sym.len == 2);

        const count = &self.element_counts[element.atomicNumber()];
        count.* +|= 1;

        var name = types.FixedString4{ .len = 4 };
        @memcpy(name.data[0..sym.len], sym);
        const counter = name.data[sym.len..];

        // Number of decimal and of base-36 counters that fit in the width
        const decimal_count: u32 = if (counter.len == 3) 1000 else 100;
        const base36_count: u32 = if (counter.len == 3) 36 * 36 * 36 else 36 * 36;
        const lettered_start = base36_count / 36 * 10; // "A00" or "A0"

        if (count.* < decimal_count) {
            const digits = std.fmt.bufPrint(counter, "{d}", .{count.*}) catch unreachable;
            name.len = @intCast(sym.len + digits.len);
            return name;
        }
        const lettered = lettered_start + (count.* - decimal_count);
        if (lettered < base36_count) {
            writeBase36(counter, lettered);
            return name;
        }
        writeBase36(&name.data, self.unprefixed_count);
        self.unprefixed_count += 1;
        return name;
    }

    /// Writes `value` in base 36, right-aligned and zero-padded to `buf.len`.
    fn writeBase36(buf: []u8, value: u32) void {
        var rest = value;
        var i = buf.len;
        while (i > 0) {
            i -= 1;
            buf[i] = base36_digits[rest % 36];
            rest /= 36;
        }
    }
};

/// Convert an SdfMolecule to a StoredComponent for the CCD classifier.
///
/// - `comp_id` = molecule name truncated to 5 chars
/// - Atom names are generated as element symbol + per-element counter
///   (C1, C2, O1...), see `AtomNamer`
/// - Bond indices and orders are preserved from the SDF data
/// - Caller must call `.deinit()` on the returned StoredComponent.
pub fn toStoredComponent(allocator: Allocator, molecule: *const SdfMolecule) !ccd_parser.StoredComponent {
    const atoms = try allocator.alloc(hybridization.CompAtom, molecule.atoms.len);
    errdefer allocator.free(atoms);

    var namer = AtomNamer{};

    for (molecule.atoms, 0..) |sdf_atom, i| {
        const sym = sdf_atom.element.symbol();
        const name = namer.next(sdf_atom.element);

        atoms[i] = hybridization.CompAtom{
            .atom_id = name.data,
            .atom_id_len = @intCast(name.len),
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
/// - Atom names generated as element + per-element index (C1, C2, O1...),
///   the same names as in `toStoredComponent`
/// - Radii default to element VdW radius (classifier will override later)
/// - When `skip_hydrogens` is true, H atoms (including D and T) are excluded
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

        // Atom names are unique per molecule. Skipped hydrogens are named
        // too, so that the other atoms keep the names of the component.
        var namer = AtomNamer{};

        for (mol.atoms) |atom| {
            const name = namer.next(atom.element);
            if (skip_hydrogens and atom.element == .H) continue;

            x[idx] = atom.x;
            y[idx] = atom.y;
            z[idx] = atom.z;
            r[idx] = atom.element.vdwRadius();
            residue[idx] = res_name;
            atom_name[idx] = name;
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

/// How many atoms `applyTopologyRadii` gave a bond-topology radius, and how
/// many an element-based fallback radius.
pub const TopologyRadiiCounts = struct {
    classified: usize = 0,
    fallback: usize = 0,
};

/// Give the atoms of an SDF/MOL molecule the CCD classifier's radii for the
/// molecule's own bond topology.
///
/// `input` is the result of `toAtomInput` for one molecule and `component`
/// the result of `toStoredComponent` for the same molecule. A molecule read
/// from an SDF/MOL file is its own component definition, so its atoms are
/// matched to `component` by their generated atom names alone and the
/// molecule title plays no part: a blank title, a title shared with another
/// molecule and a title that is also a residue name (`ALA`, `HOH`, `A`) all
/// give the radii of the bond table. This is what distinguishes it from the
/// lookup by residue name that structure files need.
///
/// An atom without a bond-topology radius (an element outside the ProtOr
/// table, or a hydrogen) gets the element-based fallback radius of the CCD
/// classifier, and keeps its radius if there is none.
pub fn applyTopologyRadii(input: *types.AtomInput, component: *const hybridization.Component) !TopologyRadiiCounts {
    const allocator = input.allocator;
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    const derived = try hybridization.deriveComponentProperties(allocator, component);
    defer allocator.free(derived);

    // Atom names are unique within a molecule and zero-padded (`AtomNamer`)
    var radius_by_name: std.AutoHashMapUnmanaged([4]u8, f64) = .empty;
    defer radius_by_name.deinit(allocator);
    try radius_by_name.ensureTotalCapacity(allocator, @intCast(derived.len));
    for (derived) |entry| radius_by_name.putAssumeCapacity(entry.atom_id, entry.props.radius);

    var counts = TopologyRadiiCounts{};
    for (input.r, atom_names, residues, 0..) |*radius, *atom_name, *residue, i| {
        if (radius_by_name.get(atom_name.data)) |derived_radius| {
            radius.* = derived_radius;
            counts.classified += 1;
        } else if (classifier.guessFallbackRadius(
            .ccd,
            if (input.element) |elements| elements[i] else null,
            residue.slice(),
            atom_name.slice(),
        )) |fallback_radius| {
            radius.* = fallback_radius;
            counts.fallback += 1;
        }
    }
    return counts;
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
/// Reads each SDF file (plain, gzip-compressed, or zstd-compressed), parses molecules, and
/// converts them to StoredComponents. Duplicate molecule names (truncated
/// to 5 chars) are skipped to avoid leaking the first entry.
///
/// The dictionary is for the residues of a structure file, which are matched
/// to its entries by residue name. A molecule without a title cannot be
/// matched to any residue and is skipped with a warning. (An SDF/MOL file
/// that is itself the input does not go through this dictionary, see
/// `applyTopologyRadii`.)
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
        const source = if (compressed.isCompressed(sdf_path))
            compressed.read(allocator, sdf_path) catch |err| {
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

        for (molecules, 1..) |mol, mol_number| {
            if (mol.name.len == 0) {
                if (!quiet) std.debug.print(
                    "Warning: molecule {d} of SDF file '{s}' has no title and cannot be matched to a residue name, skipping\n",
                    .{ mol_number, sdf_path },
                );
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
            dict.components.put(allocator, dict_key, stored) catch {
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

    // Radii derived from the bond table: the hydrogens listed in the SDF
    // make both carbons sp3 CHn (1.88), not hydrogen-free (1.61).
    const derived = try hybridization.deriveComponentProperties(allocator, &view);
    defer allocator.free(derived);
    try std.testing.expectEqual(@as(usize, 3), derived.len);
    try std.testing.expectEqual(@as(f64, 1.88), derived[0].props.radius); // C1
    try std.testing.expectEqual(@as(f64, 1.88), derived[1].props.radius); // C2
    try std.testing.expectEqual(@as(f64, 1.46), derived[2].props.radius); // O1
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

    // Atom names should be C1, C2, O1
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

/// Builds a V3000 molecule whose COUNTS line declares `counts` ("atoms bonds")
/// and whose body lists `atom_lines` atoms and `bond_lines` bonds; bond `i`
/// joins atoms `i` and `i + 1`. The bond block is left out when `bond_lines`
/// is 0.
fn buildV3000(allocator: Allocator, counts: []const u8, atom_lines: usize, bond_lines: usize) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    try writer.writeAll("mol\n     zsasa   3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");
    try writer.print("M  V30 BEGIN CTAB\nM  V30 COUNTS {s} 0 0 0\nM  V30 BEGIN ATOM\n", .{counts});
    for (0..atom_lines) |i| {
        try writer.print("M  V30 {d} C {d}.5000 0.0000 0.0000 0\n", .{ i + 1, i });
    }
    try writer.writeAll("M  V30 END ATOM\n");
    if (bond_lines > 0) {
        try writer.writeAll("M  V30 BEGIN BOND\n");
        for (0..bond_lines) |i| {
            try writer.print("M  V30 {d} 1 {d} {d}\n", .{ i + 1, i + 1, i + 2 });
        }
        try writer.writeAll("M  V30 END BOND\n");
    }
    try writer.writeAll("M  V30 END CTAB\nM  END\n$$$$\n");

    return aw.toOwnedSlice();
}

fn expectV3000Error(expected: SdfError, counts: []const u8, atom_lines: usize, bond_lines: usize) !void {
    const allocator = std.testing.allocator;
    const source = try buildV3000(allocator, counts, atom_lines, bond_lines);
    defer allocator.free(source);
    try std.testing.expectError(expected, parse(allocator, source));
}

test "parse V3000 rejects more atoms than COUNTS declares" {
    // Each of these used to write past the capacity reserved from COUNTS.
    try expectV3000Error(error.InvalidV3000, "0 0", 9, 8);
    try expectV3000Error(error.InvalidV3000, "2 0", 9, 0);
    try expectV3000Error(error.InvalidV3000, "8 8", 9, 8);
    try expectV3000Error(error.InvalidV3000, "1 0", 5000, 0);
}

test "parse V3000 rejects more bonds than COUNTS declares" {
    try expectV3000Error(error.InvalidV3000, "9 0", 9, 8);
    try expectV3000Error(error.InvalidV3000, "9 7", 9, 8);
    try expectV3000Error(error.InvalidV3000, "5000 1", 5000, 4999);
}

test "parse V3000 rejects fewer atoms or bonds than COUNTS declares" {
    try expectV3000Error(error.InvalidV3000, "9 8", 8, 7);
    try expectV3000Error(error.InvalidV3000, "9 8", 0, 0);
    try expectV3000Error(error.InvalidV3000, "9 8", 9, 7);
    try expectV3000Error(error.InvalidV3000, "9 8", 9, 0);
}

test "parse V3000 rejects a body cut off before the declared counts" {
    const allocator = std.testing.allocator;
    const source = try buildV3000(allocator, "9 8", 9, 8);
    defer allocator.free(source);

    // Input that ends inside the atom block, and inside the bond block
    const in_atoms = std.mem.find(u8, source, "M  V30 5 C").?;
    try std.testing.expectError(error.InvalidV3000, parse(allocator, source[0..in_atoms]));
    const in_bonds = std.mem.find(u8, source, "M  V30 5 1").?;
    try std.testing.expectError(error.InvalidV3000, parse(allocator, source[0..in_bonds]));
}

test "parse V3000 checks bond indices against the atoms that were read" {
    // The second bond names atom 3, and the molecule has two atoms.
    try expectV3000Error(error.BondIndexOutOfRange, "2 2", 2, 2);
}

test "parse V3000 keeps every atom and bond of a body that matches COUNTS" {
    const allocator = std.testing.allocator;
    const source = try buildV3000(allocator, "9 8", 9, 8);
    defer allocator.free(source);

    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    const mol = molecules[0];
    try std.testing.expectEqualStrings("mol", mol.name);
    try std.testing.expectEqual(@as(usize, 9), mol.atoms.len);
    try std.testing.expectEqual(@as(usize, 8), mol.bonds.len);
    for (mol.atoms, 0..) |atom, i| {
        try std.testing.expectEqual(elem.Element.C, atom.element);
        try std.testing.expectEqual(@as(f64, @floatFromInt(i)) + 0.5, atom.x);
        try std.testing.expectEqual(@as(f64, 0.0), atom.y);
        try std.testing.expectEqual(@as(f64, 0.0), atom.z);
    }
    for (mol.bonds, 0..) |bond, i| {
        try std.testing.expectEqual(@as(u16, @intCast(i)), bond.atom_idx_1);
        try std.testing.expectEqual(@as(u16, @intCast(i + 1)), bond.atom_idx_2);
        try std.testing.expectEqual(hybridization.BondOrder.single, bond.order);
    }
}

test "parse V3000 molecules without a bond block or M  END" {
    const allocator = std.testing.allocator;
    // The first molecule has no bond block and ends at $$$$ without M  END.
    const source =
        \\first
        \\     zsasa   3D
        \\
        \\  0  0  0  0  0  0  0  0  0  0999 V3000
        \\M  V30 BEGIN CTAB
        \\M  V30 COUNTS 2 0 0 0 0
        \\M  V30 BEGIN ATOM
        \\M  V30 1 O 0.0000 0.0000 0.0000 0
        \\M  V30 2 N 1.2000 0.0000 0.0000 0
        \\M  V30 END ATOM
        \\M  V30 END CTAB
        \\$$$$
        \\second
        \\     zsasa   3D
        \\
        \\  0  0  0  0  0  0  0  0  0  0999 V3000
        \\M  V30 BEGIN CTAB
        \\M  V30 COUNTS 2 1 0 0 0
        \\M  V30 BEGIN ATOM
        \\M  V30 1 C 0.0000 0.0000 0.0000 0
        \\M  V30 2 S 1.8000 0.0000 0.0000 0
        \\M  V30 END ATOM
        \\M  V30 BEGIN BOND
        \\M  V30 1 2 1 2
        \\M  V30 END BOND
        \\M  V30 END CTAB
        \\M  END
        \\$$$$
    ;
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    try std.testing.expectEqual(@as(usize, 2), molecules.len);

    try std.testing.expectEqualStrings("first", molecules[0].name);
    try std.testing.expectEqual(@as(usize, 2), molecules[0].atoms.len);
    try std.testing.expectEqual(elem.Element.O, molecules[0].atoms[0].element);
    try std.testing.expectEqual(elem.Element.N, molecules[0].atoms[1].element);
    try std.testing.expectApproxEqAbs(@as(f64, 1.2), molecules[0].atoms[1].x, 0.001);
    try std.testing.expectEqual(@as(usize, 0), molecules[0].bonds.len);

    try std.testing.expectEqualStrings("second", molecules[1].name);
    try std.testing.expectEqual(@as(usize, 2), molecules[1].atoms.len);
    try std.testing.expectEqual(elem.Element.S, molecules[1].atoms[1].element);
    try std.testing.expectApproxEqAbs(@as(f64, 1.8), molecules[1].atoms[1].x, 0.001);
    try std.testing.expectEqual(@as(usize, 1), molecules[1].bonds.len);
    try std.testing.expectEqual(@as(u16, 0), molecules[1].bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), molecules[1].bonds[0].atom_idx_2);
    try std.testing.expectEqual(hybridization.BondOrder.double, molecules[1].bonds[0].order);
}

test "parse V2000 stores exactly the atoms and bonds its counts line declares" {
    const allocator = std.testing.allocator;
    const header = "mol\n     zsasa   3D\n\n";
    const atom_line = "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n";
    const bond_line = "  1  2  1  0  0  0  0\n";
    const end = "M  END\n$$$$\n";

    // Fewer atom lines than declared: "M  END" is not an atom line
    try std.testing.expectError(error.InvalidAtomLine, parse(
        allocator,
        header ++ "  3  0  0  0  0  0  0  0  0  0999 V2000\n" ++ atom_line ++ atom_line ++ end,
    ));
    // Input that ends before the declared atoms, or before the declared bonds
    try std.testing.expectError(error.InvalidAtomLine, parse(
        allocator,
        header ++ "999  0  0  0  0  0  0  0  0  0999 V2000\n" ++ atom_line,
    ));
    try std.testing.expectError(error.InvalidBondLine, parse(
        allocator,
        header ++ "  2999  0  0  0  0  0  0  0  0999 V2000\n" ++ atom_line ++ atom_line ++ bond_line,
    ));
    // More atom lines than declared: the extra one is read as a bond line
    try std.testing.expectError(error.InvalidInteger, parse(
        allocator,
        header ++ "  1  1  0  0  0  0  0  0  0  0999 V2000\n" ++ atom_line ++ atom_line ++ end,
    ));

    // With no bonds declared, lines past the counts are skipped, not stored
    const molecules = try parse(
        allocator,
        header ++ "  1  0  0  0  0  0  0  0  0  0999 V2000\n" ++ atom_line ++ atom_line ++ bond_line ++ end,
    );
    defer freeMolecules(allocator, molecules);
    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    try std.testing.expectEqual(@as(usize, 1), molecules[0].atoms.len);
    try std.testing.expectEqual(@as(usize, 0), molecules[0].bonds.len);
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

// Header lines 2 and 3 (program/timestamp and comment) of the test records
const test_header_rest = "     zsasa   3D\n\n";

// Carbon monoxide as a V2000 body (counts line to `$$$$`), and hydrogen
// cyanide's heavy atoms as a V3000 body: two atoms and one bond each.
const test_v2000_body =
    "  2  1  0  0  0  0  0  0  0  0999 V2000\n" ++
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    1.1300    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  1  2  3  0  0  0  0\n" ++
    "M  END\n$$$$\n";
const test_v3000_body =
    "  0  0  0  0  0  0  0  0  0  0999 V3000\n" ++
    "M  V30 BEGIN CTAB\n" ++
    "M  V30 COUNTS 2 1 0 0 0\n" ++
    "M  V30 BEGIN ATOM\n" ++
    "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
    "M  V30 2 N 1.1600 0.0000 0.0000 0\n" ++
    "M  V30 END ATOM\n" ++
    "M  V30 BEGIN BOND\n" ++
    "M  V30 1 3 1 2\n" ++
    "M  V30 END BOND\n" ++
    "M  V30 END CTAB\n" ++
    "M  END\n$$$$\n";

const TestRecord = struct {
    name: []const u8,
    /// Element of the second atom: O for `test_v2000_body`, N for `test_v3000_body`
    second: elem.Element,
};

/// Parses `source` as it is and with CRLF line endings, and expects the given
/// records, each with the two atoms and one bond of the test bodies.
fn expectTestRecords(source: []const u8, expected: []const TestRecord) !void {
    const allocator = std.testing.allocator;
    const crlf = try std.mem.replaceOwned(u8, allocator, source, "\n", "\r\n");
    defer allocator.free(crlf);

    for ([_][]const u8{ source, crlf }) |text| {
        const molecules = try parse(allocator, text);
        defer freeMolecules(allocator, molecules);

        try std.testing.expectEqual(expected.len, molecules.len);
        for (expected, molecules) |want, mol| {
            try std.testing.expectEqualStrings(want.name, mol.name);
            try std.testing.expectEqual(@as(usize, 2), mol.atoms.len);
            try std.testing.expectEqual(elem.Element.C, mol.atoms[0].element);
            try std.testing.expectEqual(want.second, mol.atoms[1].element);
            try std.testing.expectEqual(@as(usize, 1), mol.bonds.len);
            try std.testing.expectEqual(hybridization.BondOrder.triple, mol.bonds[0].order);
        }
    }
}

test "parse accepts a blank title" {
    // RDKit writes a blank first line for a molecule without a name
    try expectTestRecords(
        "\n" ++ test_header_rest ++ test_v2000_body,
        &.{.{ .name = "", .second = .O }},
    );
    try expectTestRecords(
        "\n" ++ test_header_rest ++ test_v3000_body,
        &.{.{ .name = "", .second = .N }},
    );
    // A MOL file: no `$$$$` at the end
    try expectTestRecords(
        "\n" ++ test_header_rest ++ test_v2000_body[0 .. test_v2000_body.len - "$$$$\n".len],
        &.{.{ .name = "", .second = .O }},
    );
    // All three header lines blank
    try expectTestRecords(
        "\n\n\n" ++ test_v2000_body,
        &.{.{ .name = "", .second = .O }},
    );
}

test "parse accepts a blank title after a $$$$ separator" {
    // Only the second title is blank
    try expectTestRecords(
        "first\n" ++ test_header_rest ++ test_v2000_body ++
            "\n" ++ test_header_rest ++ test_v2000_body,
        &.{ .{ .name = "first", .second = .O }, .{ .name = "", .second = .O } },
    );
    try expectTestRecords(
        "first\n" ++ test_header_rest ++ test_v3000_body ++
            "\n" ++ test_header_rest ++ test_v3000_body,
        &.{ .{ .name = "first", .second = .N }, .{ .name = "", .second = .N } },
    );
    // Every title blank, V2000 and V3000 mixed, and a named record at the end
    try expectTestRecords(
        "\n" ++ test_header_rest ++ test_v3000_body ++
            "\n" ++ test_header_rest ++ test_v2000_body ++
            "\n" ++ test_header_rest ++ test_v3000_body ++
            "last\n" ++ test_header_rest ++ test_v2000_body,
        &.{
            .{ .name = "", .second = .N },
            .{ .name = "", .second = .O },
            .{ .name = "", .second = .N },
            .{ .name = "last", .second = .O },
        },
    );
}

test "parse skips blank lines between records and at the end of the file" {
    const first = "first\n" ++ test_header_rest ++ test_v2000_body;
    const second = "second\n" ++ test_header_rest ++ test_v3000_body;
    const both = [_]TestRecord{ .{ .name = "first", .second = .O }, .{ .name = "second", .second = .N } };

    // Trailing blank lines after the last `$$$$` are not another record
    try expectTestRecords(first ++ second ++ "\n", &both);
    try expectTestRecords(first ++ second ++ "\n\n\n\n\n", &both);
    // Stray blank lines before the first record and between records
    try expectTestRecords("\n" ++ first ++ "\n" ++ second, &both);
    try expectTestRecords("\n\n" ++ first ++ "\n\n" ++ second, &both);
    try expectTestRecords("\n\n\n\n" ++ first ++ "\n\n\n\n\n" ++ second ++ "\n\n", &both);

    // The fourth line from a stray blank line can be a title made of digits,
    // which reads as two counts but is not followed by an atom line.
    try expectTestRecords(
        first ++ "\n\n\n" ++ "5280343\n" ++ test_header_rest ++ test_v2000_body,
        &.{ .{ .name = "first", .second = .O }, .{ .name = "5280343", .second = .O } },
    );

    // A stray blank line before a record whose title is blank: the record
    // starts at the blank line that has the counts line three lines below it.
    try expectTestRecords(
        first ++ "\n" ++ "\n" ++ test_header_rest ++ test_v2000_body,
        &.{ .{ .name = "first", .second = .O }, .{ .name = "", .second = .O } },
    );
}

test "parse accepts a blank title above a counts line without a version tag" {
    const untagged = "  2  1  0  0  0  0  0  0  0  0  1\n" ++ test_v2000_body["  2  1  0  0  0  0  0  0  0  0999 V2000\n".len..];
    try expectTestRecords(
        "\n" ++ test_header_rest ++ untagged ++ "\n" ++ test_header_rest ++ untagged,
        &.{ .{ .name = "", .second = .O }, .{ .name = "", .second = .O } },
    );
}

test "parse still rejects a record without a counts line" {
    const allocator = std.testing.allocator;
    const atom_line = "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n";

    // Counts line missing: with a title, and with a blank title
    try std.testing.expectError(error.InvalidCountsLine, parse(
        allocator,
        "mol\n" ++ test_header_rest ++ atom_line ++ atom_line ++ "M  END\n$$$$\n",
    ));
    try std.testing.expectError(error.InvalidCountsLine, parse(
        allocator,
        "\n" ++ test_header_rest ++ atom_line ++ atom_line ++ "M  END\n$$$$\n",
    ));
    // Header one line short
    try std.testing.expectError(error.InvalidCountsLine, parse(
        allocator,
        "mol\n" ++ "     zsasa   3D\n" ++ test_v2000_body,
    ));
    // Only blank lines
    try std.testing.expectError(error.EmptySdf, parse(allocator, "\n\n\n\n\n\n"));
}

const test_ethanol_v2000 =
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
    \\
;

const test_ethanol_v3000 =
    \\ethanol
    \\     zsasa   3D
    \\
    \\  0  0  0  0  0  0  0  0  0  0999 V3000
    \\M  V30 BEGIN CTAB
    \\M  V30 COUNTS 9 8 0 0 0
    \\M  V30 BEGIN ATOM
    \\M  V30 1 C 0.0000 0.0000 0.0000 0
    \\M  V30 2 C 1.5200 0.0000 0.0000 0
    \\M  V30 3 O 2.0800 1.2124 0.0000 0
    \\M  V30 4 H -0.5200 0.9400 0.0000 0
    \\M  V30 5 H -0.5200 -0.5100 0.8900 0
    \\M  V30 6 H -0.5200 -0.5100 -0.8900 0
    \\M  V30 7 H 1.8800 -0.5100 0.8900 0
    \\M  V30 8 H 1.8800 -0.5100 -0.8900 0
    \\M  V30 9 H 2.9200 1.2124 0.0000 0
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
    \\
;

/// Parses ethanol with the first `n_renamed` of its hydrogens (the three on
/// C1 come first) written as `isotope`, and checks which atoms are kept and
/// the radii derived from the bond table.
fn expectHydrogenIsotope(ethanol: []const u8, comptime isotope: []const u8, n_renamed: usize) !void {
    const allocator = std.testing.allocator;

    // " H " matches the V2000 symbol column and the V3000 atom type
    var source = try allocator.dupe(u8, ethanol);
    defer allocator.free(source);
    for (0..n_renamed) |_| {
        const at = std.mem.find(u8, source, " H ").?;
        const renamed = try std.mem.concat(allocator, u8, &.{ source[0..at], " " ++ isotope ++ " ", source[at + 3 ..] });
        allocator.free(source);
        source = renamed;
    }

    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);
    try std.testing.expectEqual(@as(usize, 1), molecules.len);
    try std.testing.expectEqual(@as(usize, 9), molecules[0].atoms.len);
    for (molecules[0].atoms[3..]) |atom| try std.testing.expectEqual(elem.Element.H, atom.element);

    // Excluded with the hydrogens: C, C and O are left
    {
        var input = try toAtomInput(allocator, molecules, true);
        defer input.deinit();
        try std.testing.expectEqual(@as(usize, 3), input.atomCount());
        try std.testing.expectEqualSlices(u8, &.{ 6, 6, 8 }, input.element.?);
        const names = input.atom_name.?;
        try std.testing.expectEqualStrings("C1", names[0].slice());
        try std.testing.expectEqualStrings("C2", names[1].slice());
        try std.testing.expectEqualStrings("O1", names[2].slice());
    }

    // Included with the hydrogens, as hydrogen
    {
        var input = try toAtomInput(allocator, molecules, false);
        defer input.deinit();
        try std.testing.expectEqual(@as(usize, 9), input.atomCount());
        const names = input.atom_name.?;
        for (3..9) |i| {
            try std.testing.expectEqual(@as(u8, 1), input.element.?[i]);
            try std.testing.expectEqual(elem.Element.H.vdwRadius(), input.r[i]);
        }
        try std.testing.expectEqualStrings("H1", names[3].slice());
        try std.testing.expectEqualStrings("H6", names[8].slice());
    }

    // Counted as hydrogens of their carbon: C1 is C4H3 (1.88), not
    // hydrogen-free (1.61)
    var stored = try toStoredComponent(allocator, &molecules[0]);
    defer stored.deinit();
    const view = stored.view();
    const derived = try hybridization.deriveComponentProperties(allocator, &view);
    defer allocator.free(derived);
    try std.testing.expectEqual(@as(usize, 3), derived.len);
    try std.testing.expectEqualStrings("C1", derived[0].atomIdSlice());
    try std.testing.expectEqual(@as(f64, 1.88), derived[0].props.radius);
    try std.testing.expectEqual(@as(f64, 1.88), derived[1].props.radius);
    try std.testing.expectEqual(@as(f64, 1.46), derived[2].props.radius);
}

test "deuterium and tritium are hydrogen" {
    // One deuterium, and a CD3 group
    try expectHydrogenIsotope(test_ethanol_v2000, "D", 1);
    try expectHydrogenIsotope(test_ethanol_v2000, "D", 3);
    try expectHydrogenIsotope(test_ethanol_v3000, "D", 1);
    try expectHydrogenIsotope(test_ethanol_v3000, "D", 3);
    // Tritium
    try expectHydrogenIsotope(test_ethanol_v2000, "T", 3);
    try expectHydrogenIsotope(test_ethanol_v3000, "T", 3);
    // Every hydrogen replaced
    try expectHydrogenIsotope(test_ethanol_v2000, "D", 6);
    try expectHydrogenIsotope(test_ethanol_v3000, "D", 6);
}

test "elementFromSymbol maps D and T to hydrogen only as whole symbols" {
    try std.testing.expectEqual(elem.Element.H, elementFromSymbol("H"));
    try std.testing.expectEqual(elem.Element.H, elementFromSymbol("D"));
    try std.testing.expectEqual(elem.Element.H, elementFromSymbol("T"));
    try std.testing.expectEqual(elem.Element.H, elementFromSymbol("d"));
    // Elements whose symbols start with D or T keep their element
    try std.testing.expectEqual(elem.Element.Dy, elementFromSymbol("Dy"));
    try std.testing.expectEqual(elem.Element.Ti, elementFromSymbol("Ti"));
    try std.testing.expectEqual(elem.Element.Tc, elementFromSymbol("Tc"));
    try std.testing.expectEqual(elem.Element.C, elementFromSymbol("C"));
    try std.testing.expectEqual(elem.Element.Cl, elementFromSymbol("Cl"));
    try std.testing.expectEqual(elem.Element.X, elementFromSymbol("R#"));
    try std.testing.expectEqual(elem.Element.X, elementFromSymbol(""));
}

// V3000 records without atoms: with no atom block at all, and with an empty
// atom block.
const test_v3000_empty_body =
    "  0  0  0  0  0  0  0  0  0  0999 V3000\n" ++
    "M  V30 BEGIN CTAB\n" ++
    "M  V30 COUNTS 0 0 0 0 0\n" ++
    "M  V30 END CTAB\n" ++
    "M  END\n$$$$\n";
const test_v3000_empty_block_body =
    "  0  0  0  0  0  0  0  0  0  0999 V3000\n" ++
    "M  V30 BEGIN CTAB\n" ++
    "M  V30 COUNTS 0 0 0 0 0\n" ++
    "M  V30 BEGIN ATOM\n" ++
    "M  V30 END ATOM\n" ++
    "M  V30 END CTAB\n" ++
    "M  END\n$$$$\n";

/// Parses `source` as it is and with CRLF line endings, and expects molecules
/// with the given names and atom counts.
fn expectAtomCounts(source: []const u8, expected: []const struct { []const u8, usize }) !void {
    const allocator = std.testing.allocator;
    const crlf = try std.mem.replaceOwned(u8, allocator, source, "\n", "\r\n");
    defer allocator.free(crlf);

    for ([_][]const u8{ source, crlf }) |text| {
        const molecules = try parse(allocator, text);
        defer freeMolecules(allocator, molecules);

        try std.testing.expectEqual(expected.len, molecules.len);
        for (expected, molecules) |want, mol| {
            try std.testing.expectEqualStrings(want[0], mol.name);
            try std.testing.expectEqual(want[1], mol.atoms.len);
        }
    }
}

test "parse V3000 molecule without atoms does not read into the next record" {
    const first = "first\n" ++ test_header_rest ++ test_v3000_body;
    const last = "last\n" ++ test_header_rest ++ test_v3000_body;

    for ([_][]const u8{ test_v3000_empty_body, test_v3000_empty_block_body }) |empty_body| {
        const empty = try std.mem.concat(std.testing.allocator, u8, &.{ "empty\n", test_header_rest, empty_body });
        defer std.testing.allocator.free(empty);

        // Between two molecules, first, last, and alone
        const between = try std.mem.concat(std.testing.allocator, u8, &.{ first, empty, last });
        defer std.testing.allocator.free(between);
        try expectAtomCounts(between, &.{ .{ "first", 2 }, .{ "empty", 0 }, .{ "last", 2 } });

        const leading = try std.mem.concat(std.testing.allocator, u8, &.{ empty, first, last });
        defer std.testing.allocator.free(leading);
        try expectAtomCounts(leading, &.{ .{ "empty", 0 }, .{ "first", 2 }, .{ "last", 2 } });

        const trailing = try std.mem.concat(std.testing.allocator, u8, &.{ first, last, empty });
        defer std.testing.allocator.free(trailing);
        try expectAtomCounts(trailing, &.{ .{ "first", 2 }, .{ "last", 2 }, .{ "empty", 0 } });

        try expectAtomCounts(empty, &.{.{ "empty", 0 }});

        // Followed by a V2000 record, whose atom lines are not V3000 lines
        const before_v2000 = try std.mem.concat(std.testing.allocator, u8, &.{ empty, "last\n", test_header_rest, test_v2000_body });
        defer std.testing.allocator.free(before_v2000);
        try expectAtomCounts(before_v2000, &.{ .{ "empty", 0 }, .{ "last", 2 } });
    }

    // The molecules next to an empty one are complete
    const allocator = std.testing.allocator;
    const molecules = try parse(allocator, first ++ "empty\n" ++ test_header_rest ++ test_v3000_empty_body ++ last);
    defer freeMolecules(allocator, molecules);
    try std.testing.expectEqual(@as(usize, 3), molecules.len);
    try std.testing.expectEqual(@as(usize, 0), molecules[1].atoms.len);
    try std.testing.expectEqual(@as(usize, 0), molecules[1].bonds.len);
    for ([_]usize{ 0, 2 }) |i| {
        try std.testing.expectEqual(elem.Element.C, molecules[i].atoms[0].element);
        try std.testing.expectEqual(elem.Element.N, molecules[i].atoms[1].element);
        try std.testing.expectApproxEqAbs(@as(f64, 1.16), molecules[i].atoms[1].x, 1e-9);
        try std.testing.expectEqual(@as(usize, 1), molecules[i].bonds.len);
    }

    // An empty molecule gives no atoms, and the others keep theirs
    var empty_input = try toAtomInput(allocator, molecules[1..2], false);
    defer empty_input.deinit();
    try std.testing.expectEqual(@as(usize, 0), empty_input.atomCount());
    var all_input = try toAtomInput(allocator, molecules, false);
    defer all_input.deinit();
    try std.testing.expectEqual(@as(usize, 4), all_input.atomCount());
    try std.testing.expectEqualStrings("A", all_input.chain_id.?[1].slice());
    try std.testing.expectEqualStrings("C", all_input.chain_id.?[2].slice());

    var stored = try toStoredComponent(allocator, &molecules[1]);
    defer stored.deinit();
    try std.testing.expectEqual(@as(usize, 0), stored.atoms.len);
}

test "parse V3000 reads atoms and bonds only inside their blocks" {
    const allocator = std.testing.allocator;
    const header = "mol\n" ++ test_header_rest ++ "  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n";
    const atoms = "M  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 N 1.16 0 0 0\nM  V30 END ATOM\n";
    const bonds = "M  V30 BEGIN BOND\nM  V30 1 3 1 2\nM  V30 END BOND\n";
    const end = "M  V30 END CTAB\nM  END\n$$$$\n";

    // Lines of other blocks, before, between and after the atom and bond
    // blocks, are neither atoms nor bonds
    const sgroup = "M  V30 BEGIN SGROUP\nM  V30 1 SUP 0 ATOMS=(2 1 2) LABEL=CN\nM  V30 END SGROUP\n";
    const collection = "M  V30 BEGIN COLLECTION\nM  V30 MDLV30/STEABS ATOMS=(1 1)\nM  V30 END COLLECTION\n";
    try expectAtomCounts(
        header ++ "M  V30 COUNTS 2 1 1 0 0\n" ++ collection ++ atoms ++ sgroup ++ bonds ++ collection ++ end,
        &.{.{ "mol", 2 }},
    );

    // The bond block may come without an atom block only when it is empty
    try expectAtomCounts(
        header ++ "M  V30 COUNTS 0 0 0 0 0\nM  V30 BEGIN BOND\nM  V30 END BOND\n" ++ end,
        &.{.{ "mol", 0 }},
    );
    try std.testing.expectError(error.BondIndexOutOfRange, parse(
        allocator,
        header ++ "M  V30 COUNTS 0 1 0 0 0\n" ++ bonds ++ end,
    ));

    // An atom block that is not closed before `M  END`, `$$$$` or the end of
    // the input is an error, with and without a record after it
    const next = "next\n" ++ test_header_rest ++ test_v3000_body;
    const open_atoms = header ++ "M  V30 COUNTS 2 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0\nM  V30 2 N 1.16 0 0 0\n";
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atoms));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atoms ++ "M  END\n$$$$\n"));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atoms ++ "M  END\n$$$$\n" ++ next));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atoms ++ "$$$$\n" ++ next));
    // The same for a bond block
    const open_bonds = header ++ "M  V30 COUNTS 2 1 0 0 0\n" ++ atoms ++ "M  V30 BEGIN BOND\nM  V30 1 3 1 2\n";
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_bonds));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_bonds ++ "M  END\n$$$$\n" ++ next));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_bonds ++ "$$$$\n" ++ next));

    // Atoms of the next record do not make up for atoms that are missing
    try std.testing.expectError(error.InvalidV3000, parse(
        allocator,
        header ++ "M  V30 COUNTS 2 1 0 0 0\n" ++ end ++ next,
    ));

    // An atom or bond block needs a COUNTS line before it
    try std.testing.expectError(error.InvalidCountsLine, parse(allocator, header ++ atoms ++ bonds ++ end));
    try std.testing.expectError(error.InvalidCountsLine, parse(allocator, header ++ bonds ++ end));
    try std.testing.expectError(error.InvalidCountsLine, parse(allocator, header ++ end));
}

/// Parses a V3000 record made of `ctab` (the lines between `BEGIN CTAB` and
/// `END CTAB`), as it is and with CRLF line endings, and expects the atoms
/// and the bond of `test_v3000_body`.
fn expectV3000Cyanide(comptime ctab: []const u8) !void {
    const allocator = std.testing.allocator;
    const source = "mol\n" ++ test_header_rest ++ "  0  0  0  0  0  0  0  0  0  0999 V3000\n" ++
        "M  V30 BEGIN CTAB\n" ++ ctab ++ "M  V30 END CTAB\nM  END\n$$$$\n" ++
        "next\n" ++ test_header_rest ++ test_v2000_body;
    const crlf = try std.mem.replaceOwned(u8, allocator, source, "\n", "\r\n");
    defer allocator.free(crlf);

    for ([_][]const u8{ source, crlf }) |text| {
        const molecules = try parse(allocator, text);
        defer freeMolecules(allocator, molecules);

        try std.testing.expectEqual(@as(usize, 2), molecules.len);
        const mol = molecules[0];
        try std.testing.expectEqualStrings("mol", mol.name);
        try std.testing.expectEqual(@as(usize, 2), mol.atoms.len);
        try std.testing.expectEqual(elem.Element.C, mol.atoms[0].element);
        try std.testing.expectEqual(@as(f64, 0.0), mol.atoms[0].x);
        try std.testing.expectEqual(elem.Element.N, mol.atoms[1].element);
        try std.testing.expectEqual(@as(f64, 1.16), mol.atoms[1].x);
        try std.testing.expectEqual(@as(f64, 0.25), mol.atoms[1].y);
        try std.testing.expectEqual(@as(f64, -0.5), mol.atoms[1].z);
        try std.testing.expectEqual(@as(usize, 1), mol.bonds.len);
        try std.testing.expectEqual(@as(u16, 0), mol.bonds[0].atom_idx_1);
        try std.testing.expectEqual(@as(u16, 1), mol.bonds[0].atom_idx_2);
        try std.testing.expectEqual(hybridization.BondOrder.triple, mol.bonds[0].order);

        // The record after it starts where it should
        try std.testing.expectEqualStrings("next", molecules[1].name);
        try std.testing.expectEqual(@as(usize, 2), molecules[1].atoms.len);
    }
}

const test_v3000_counts = "M  V30 COUNTS 2 1 0 0 0\n";
const test_v3000_atoms =
    "M  V30 BEGIN ATOM\n" ++
    "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
    "M  V30 2 N 1.1600 0.2500 -0.5000 0\n" ++
    "M  V30 END ATOM\n";
const test_v3000_bonds = "M  V30 BEGIN BOND\nM  V30 1 3 1 2\nM  V30 END BOND\n";

test "parse V3000 joins continuation lines in the atom block" {
    // Not continued, for reference
    try expectV3000Cyanide(test_v3000_counts ++ test_v3000_atoms ++ test_v3000_bonds);

    // Continued between two fields
    try expectV3000Cyanide(test_v3000_counts ++
        "M  V30 BEGIN ATOM\n" ++
        "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
        "M  V30 2 N 1.1600 0.2500 -\n" ++
        "M  V30 -0.5000 0\n" ++
        "M  V30 END ATOM\n" ++ test_v3000_bonds);
    // Continued inside a field: the two parts are joined without a space
    try expectV3000Cyanide(test_v3000_counts ++
        "M  V30 BEGIN ATOM\n" ++
        "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
        "M  V30 2 N 1.16-\n" ++
        "M  V30 00 0.2500 -0.5000 0\n" ++
        "M  V30 END ATOM\n" ++ test_v3000_bonds);
    // Continued over several lines, in the first and in the last atom
    try expectV3000Cyanide(test_v3000_counts ++
        "M  V30 BEGIN ATOM\n" ++
        "M  V30 1 -\n" ++
        "M  V30 C -\n" ++
        "M  V30 0.0000 0.0000 -\n" ++
        "M  V30 0.0000 0\n" ++
        "M  V30 2 N 1.1600 0.2500 -0.5000 0 -\n" ++
        "M  V30 CFG=0 -\n" ++
        "M  V30 MASS=15\n" ++
        "M  V30 END ATOM\n" ++ test_v3000_bonds);
}

test "parse V3000 joins continuation lines in the bond block and elsewhere" {
    // Bond block
    try expectV3000Cyanide(test_v3000_counts ++ test_v3000_atoms ++
        "M  V30 BEGIN BOND\n" ++
        "M  V30 1 3 -\n" ++
        "M  V30 1 2\n" ++
        "M  V30 END BOND\n");
    try expectV3000Cyanide(test_v3000_counts ++ test_v3000_atoms ++
        "M  V30 BEGIN BOND\n" ++
        "M  V30 1 -\n" ++
        "M  V30 3 1 -\n" ++
        "M  V30 2 CFG=0\n" ++
        "M  V30 END BOND\n");
    // COUNTS line
    try expectV3000Cyanide("M  V30 COUNTS 2 -\n" ++ "M  V30 1 0 0 0\n" ++ test_v3000_atoms ++ test_v3000_bonds);
    // Block delimiters
    try expectV3000Cyanide(test_v3000_counts ++
        "M  V30 BEGIN -\n" ++
        "M  V30 ATOM\n" ++
        "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
        "M  V30 2 N 1.1600 0.2500 -0.5000 0\n" ++
        "M  V30 END AT-\n" ++
        "M  V30 OM\n" ++ test_v3000_bonds);
    // A block that is not read: its continuation line is not a new line,
    // whatever it looks like
    try expectV3000Cyanide(test_v3000_counts ++ test_v3000_atoms ++
        "M  V30 BEGIN SGROUP\n" ++
        "M  V30 1 SUP 0 ATOMS=(2 1 2) LABEL=-\n" ++
        "M  V30 BEGIN BOND\n" ++
        "M  V30 END SGROUP\n" ++ test_v3000_bonds);
}

test "parse V3000 does not read a continuation line as another atom or bond" {
    const allocator = std.testing.allocator;
    const header = "mol\n" ++ test_header_rest ++ "  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n";
    const end = "M  V30 END CTAB\nM  END\n$$$$\n";

    // The second line continues atom 2. It is not a third atom: the molecule
    // has the two atoms that COUNTS declares...
    const phantom_atom =
        "M  V30 BEGIN ATOM\n" ++
        "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
        "M  V30 2 N 1.1600 0.2500 -0.5000 0 -\n" ++
        "M  V30 3 Fe 9.0000 9.0000 9.0000 0\n" ++
        "M  V30 END ATOM\n";
    try expectV3000Cyanide(test_v3000_counts ++ phantom_atom ++ test_v3000_bonds);
    // ...and a COUNTS line that declares three is wrong
    try std.testing.expectError(error.InvalidV3000, parse(
        allocator,
        header ++ "M  V30 COUNTS 3 1 0 0 0\n" ++ phantom_atom ++ test_v3000_bonds ++ end,
    ));

    // The same for a bond
    const phantom_bond =
        "M  V30 BEGIN BOND\n" ++
        "M  V30 1 3 1 2 -\n" ++
        "M  V30 2 1 2 1\n" ++
        "M  V30 END BOND\n";
    try expectV3000Cyanide(test_v3000_counts ++ test_v3000_atoms ++ phantom_bond);
    try std.testing.expectError(error.InvalidV3000, parse(
        allocator,
        header ++ "M  V30 COUNTS 2 2 0 0 0\n" ++ test_v3000_atoms ++ phantom_bond ++ end,
    ));
}

test "parse V3000 rejects a continued line without a continuation" {
    const allocator = std.testing.allocator;
    const header = "mol\n" ++ test_header_rest ++ "  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\n";
    const open_atom = header ++ test_v3000_counts ++
        "M  V30 BEGIN ATOM\n" ++
        "M  V30 1 C 0.0000 0.0000 0.0000 0\n" ++
        "M  V30 2 N 1.1600 0.2500 -\n";

    // End of the input, with and without a final newline
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atom));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atom[0 .. open_atom.len - 1]));
    // A line that is not a V3000 line
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atom ++ "M  END\n$$$$\n"));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atom ++ "-0.5000 0\nM  V30 END ATOM\n"));
    try std.testing.expectError(error.InvalidV3000, parse(allocator, open_atom ++ "\nM  V30 -0.5000 0\nM  V30 END ATOM\n"));
    try std.testing.expectError(error.InvalidV3000, parse(
        allocator,
        open_atom ++ "$$$$\n" ++ "next\n" ++ test_header_rest ++ test_v3000_body,
    ));
}

const classifier_ccd = @import("classifier_ccd.zig");

/// Builds a V3000 molecule that is an unbranched chain of `n` atoms of
/// `element`. Bond `i` (1-based) joins atoms `i` and `i + 1` and is single,
/// except for the bonds listed in `orders` as `{ i, bond type }`.
fn buildV3000Chain(
    allocator: Allocator,
    name: []const u8,
    element: []const u8,
    n: usize,
    orders: []const struct { usize, u8 },
) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    try writer.print("{s}\n     zsasa   3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n", .{name});
    try writer.print("M  V30 BEGIN CTAB\nM  V30 COUNTS {d} {d} 0 0 0\nM  V30 BEGIN ATOM\n", .{ n, n - 1 });
    for (0..n) |i| {
        try writer.print("M  V30 {d} {s} {d}.5000 0.0000 0.0000 0\n", .{ i + 1, element, i });
    }
    try writer.writeAll("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for (1..n) |i| {
        var order: u8 = 1;
        for (orders) |entry| {
            if (entry[0] == i) order = entry[1];
        }
        try writer.print("M  V30 {d} {d} {d} {d}\n", .{ i, order, i, i + 1 });
    }
    try writer.writeAll("M  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n");

    return aw.toOwnedSlice();
}

/// Radii that the CCD classifier gives the atoms of `molecule` when it looks
/// them up by the names of `toAtomInput`, after checking that those names are
/// unique and are the names of `toStoredComponent`. Caller frees the result.
fn classifiedRadii(allocator: Allocator, molecule: *const SdfMolecule, names_out: *[]types.FixedString4) ![]?f64 {
    const mol_slice: []const SdfMolecule = @as([*]const SdfMolecule, @ptrCast(molecule))[0..1];
    var input = try toAtomInput(allocator, mol_slice, true);
    defer input.deinit();
    var stored = try toStoredComponent(allocator, molecule);
    defer stored.deinit();

    const names = input.atom_name.?;
    var seen = std.StringHashMapUnmanaged(void).empty;
    defer seen.deinit(allocator);
    var heavy: usize = 0;
    for (stored.atoms, molecule.atoms) |*comp_atom, sdf_atom| {
        try std.testing.expect(!seen.contains(comp_atom.atomIdSlice()));
        try seen.put(allocator, comp_atom.atomIdSlice(), {});
        if (sdf_atom.element == .H) continue;
        try std.testing.expectEqualStrings(comp_atom.atomIdSlice(), names[heavy].slice());
        heavy += 1;
    }
    try std.testing.expectEqual(heavy, names.len);

    var clf = classifier_ccd.CcdClassifier.init(allocator);
    defer clf.deinit();
    const view = stored.view();
    try clf.addComponent(&view);

    const radii = try allocator.alloc(?f64, names.len);
    errdefer allocator.free(radii);
    for (radii, names, input.residue.?) |*radius, *atom_name, *residue| {
        radius.* = clf.getRadius(residue.slice(), atom_name.slice());
    }
    names_out.* = try allocator.dupe(types.FixedString4, names);
    return radii;
}

test "atom names stay unique past 999 atoms of one element" {
    const allocator = std.testing.allocator;

    // 1,050 carbons: a triple bond between atoms 1000 and 1001 (sp, 1.61)
    // and a double bond between atoms 1049 and 1050 (sp2 with hydrogens,
    // 1.76). Every other carbon is sp3 with hydrogens (1.88).
    const source = try buildV3000Chain(allocator, "chain", "C", 1050, &.{ .{ 1000, 3 }, .{ 1049, 2 } });
    defer allocator.free(source);
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);
    try std.testing.expectEqual(@as(usize, 1050), molecules[0].atoms.len);

    var names: []types.FixedString4 = &.{};
    const radii = try classifiedRadii(allocator, &molecules[0], &names);
    defer allocator.free(radii);
    defer allocator.free(names);

    try std.testing.expectEqual(@as(usize, 1050), radii.len);
    for (radii, 1..) |radius, serial| {
        const expected: f64 = switch (serial) {
            1000, 1001 => 1.61,
            1049, 1050 => 1.76,
            else => 1.88,
        };
        try std.testing.expectEqual(@as(?f64, expected), radius);
    }

    try std.testing.expectEqualStrings("C1", names[0].slice());
    try std.testing.expectEqualStrings("C10", names[9].slice());
    try std.testing.expectEqualStrings("C999", names[998].slice());
    try std.testing.expectEqualStrings("CA00", names[999].slice());
    try std.testing.expectEqualStrings("CA01", names[1000].slice());
    try std.testing.expectEqualStrings("CA1E", names[1049].slice());
}

test "atom names stay unique past 99 atoms of a two-letter element" {
    const allocator = std.testing.allocator;

    // A carbon chain as above with 120 selenium atoms after it: the 100th
    // selenium is not named like the 10th, and no selenium like a carbon.
    const carbons = 1010;
    const total = carbons + 120;
    var aw = std.Io.Writer.Allocating.init(allocator);
    defer aw.deinit();
    const writer = &aw.writer;
    try writer.writeAll("mixed\n     zsasa   3D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\n");
    try writer.print("M  V30 BEGIN CTAB\nM  V30 COUNTS {d} {d} 0 0 0\nM  V30 BEGIN ATOM\n", .{ total, total - 1 });
    for (0..total) |i| {
        try writer.print("M  V30 {d} {s} {d}.5000 0.0000 0.0000 0\n", .{ i + 1, if (i < carbons) "C" else "Se", i });
    }
    try writer.writeAll("M  V30 END ATOM\nM  V30 BEGIN BOND\n");
    for (1..total) |i| {
        // A triple bond between carbons 1000 and 1001
        try writer.print("M  V30 {d} {d} {d} {d}\n", .{ i, @as(u8, if (i == 1000) 3 else 1), i, i + 1 });
    }
    try writer.writeAll("M  V30 END BOND\nM  V30 END CTAB\nM  END\n$$$$\n");

    const molecules = try parse(allocator, aw.written());
    defer freeMolecules(allocator, molecules);

    var names: []types.FixedString4 = &.{};
    const radii = try classifiedRadii(allocator, &molecules[0], &names);
    defer allocator.free(radii);
    defer allocator.free(names);

    try std.testing.expectEqual(@as(usize, total), radii.len);
    for (radii, 1..) |radius, serial| {
        const expected: f64 = switch (serial) {
            1000, 1001 => 1.61,
            // Including the last carbon, which is bonded to selenium
            1...999, 1002...carbons => 1.88,
            else => 1.90,
        };
        try std.testing.expectEqual(@as(?f64, expected), radius);
    }

    try std.testing.expectEqualStrings("Se1", names[carbons].slice());
    try std.testing.expectEqualStrings("Se10", names[carbons + 9].slice());
    try std.testing.expectEqualStrings("Se99", names[carbons + 98].slice());
    try std.testing.expectEqualStrings("SeA0", names[carbons + 99].slice());
    try std.testing.expectEqualStrings("SeA1", names[carbons + 100].slice());
    try std.testing.expectEqualStrings("SeAK", names[carbons + 119].slice());
}

test "atom names of a small molecule are the element symbol and a decimal counter" {
    const allocator = std.testing.allocator;
    const source =
        \\small
        \\     zsasa   3D
        \\
        \\ 12  0  0  0  0  0  0  0  0  0999 V2000
        \\    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    2.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
        \\    4.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\    6.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\    8.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        \\   10.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
        \\   12.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
        \\   14.0000    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0
        \\   16.0000    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
        \\   18.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
        \\   20.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        \\   22.0000    0.0000    0.0000 Fe  0  0  0  0  0  0  0  0  0  0  0  0
        \\M  END
        \\$$$$
    ;
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    const all = [_][]const u8{ "C1", "Cl1", "H1", "C2", "O1", "Br1", "Cl2", "X1", "H2", "N1", "C3", "Fe1" };
    const heavy = [_][]const u8{ "C1", "Cl1", "C2", "O1", "Br1", "Cl2", "X1", "N1", "C3", "Fe1" };

    var stored = try toStoredComponent(allocator, &molecules[0]);
    defer stored.deinit();
    try std.testing.expectEqual(all.len, stored.atoms.len);
    for (all, stored.atoms) |want, *atom| try std.testing.expectEqualStrings(want, atom.atomIdSlice());

    var with_h = try toAtomInput(allocator, molecules, false);
    defer with_h.deinit();
    try std.testing.expectEqual(all.len, with_h.atomCount());
    for (all, with_h.atom_name.?) |want, *name| try std.testing.expectEqualStrings(want, name.slice());

    var without_h = try toAtomInput(allocator, molecules, true);
    defer without_h.deinit();
    try std.testing.expectEqual(heavy.len, without_h.atomCount());
    for (heavy, without_h.atom_name.?) |want, *name| try std.testing.expectEqualStrings(want, name.slice());
}

test "AtomNamer switches to base 36 and then drops the element symbol" {
    var namer = AtomNamer{};
    var last = types.FixedString4{};

    // One-letter symbol: three characters for the counter
    for (1..34_697) |count| {
        last = namer.next(.C);
        switch (count) {
            1 => try std.testing.expectEqualStrings("C1", last.slice()),
            99 => try std.testing.expectEqualStrings("C99", last.slice()),
            100 => try std.testing.expectEqualStrings("C100", last.slice()),
            999 => try std.testing.expectEqualStrings("C999", last.slice()),
            1000 => try std.testing.expectEqualStrings("CA00", last.slice()),
            1035 => try std.testing.expectEqualStrings("CA0Z", last.slice()),
            1036 => try std.testing.expectEqualStrings("CA10", last.slice()),
            34_695 => try std.testing.expectEqualStrings("CZZZ", last.slice()),
            34_696 => try std.testing.expectEqualStrings("0000", last.slice()),
            else => {},
        }
    }

    // Two-letter symbol: two characters. Atoms without a symbol share one
    // counter, which the carbon above has started.
    for (1..1038) |count| {
        last = namer.next(.Cl);
        switch (count) {
            1 => try std.testing.expectEqualStrings("Cl1", last.slice()),
            99 => try std.testing.expectEqualStrings("Cl99", last.slice()),
            100 => try std.testing.expectEqualStrings("ClA0", last.slice()),
            135 => try std.testing.expectEqualStrings("ClAZ", last.slice()),
            136 => try std.testing.expectEqualStrings("ClB0", last.slice()),
            1035 => try std.testing.expectEqualStrings("ClZZ", last.slice()),
            1036 => try std.testing.expectEqualStrings("0001", last.slice()),
            1037 => try std.testing.expectEqualStrings("0002", last.slice()),
            else => {},
        }
    }
    last = namer.next(.C);
    try std.testing.expectEqualStrings("0003", last.slice());

    // Other elements are not affected
    last = namer.next(.Ca);
    try std.testing.expectEqualStrings("Ca1", last.slice());
    last = namer.next(.X);
    try std.testing.expectEqualStrings("X1", last.slice());
}

test "AtomNamer gives every atom of the largest molecule its own name" {
    const allocator = std.testing.allocator;

    // 65,535 atoms, the most a bond can refer to. Carbon and the two-letter
    // elements that start with C all run out of names with their symbol, and
    // so does hydrogen.
    const elements = [_]elem.Element{ .C, .Cl, .Ca, .Co, .Cu, .Cs, .Cd, .Cr, .H, .N, .X };
    const counts = [_]usize{ 36_000, 1100, 1100, 1100, 1100, 1100, 1100, 1100, 20_000, 1035, 800 };

    var seen = std.AutoHashMapUnmanaged([4]u8, void).empty;
    defer seen.deinit(allocator);
    try seen.ensureTotalCapacity(allocator, 65_535);

    var namer = AtomNamer{};
    var remaining = counts;
    var total: usize = 0;
    // Interleave the elements, as a file would
    var any = true;
    while (any) {
        any = false;
        for (elements, &remaining) |element, *left| {
            if (left.* == 0) continue;
            left.* -= 1;
            any = true;
            total += 1;

            const name = namer.next(element);
            try std.testing.expect(name.len >= 2 and name.len <= 4);
            // Unused bytes are zero, so the array identifies the name
            for (name.data[name.len..]) |c| try std.testing.expectEqual(@as(u8, 0), c);
            const entry = seen.getOrPutAssumeCapacity(name.data);
            try std.testing.expect(!entry.found_existing);
        }
    }
    try std.testing.expectEqual(@as(usize, 65_535), total);
}

// Heavy atoms of acetonitrile (methyl 1.88, nitrile carbon 1.61, N 1.64) and
// of acetaldehyde (methyl 1.88, carbonyl carbon 1.76, O 1.42), as V2000
// bodies. Both are three atoms named C1, C2 and N1 or O1, so the radius of C2
// tells which bond table was used.
const test_acetonitrile_body =
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n" ++
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    1.4600    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    2.6200    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  1  2  1  0  0  0  0\n" ++
    "  2  3  3  0  0  0  0\n" ++
    "M  END\n$$$$\n";
const test_acetaldehyde_body =
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n" ++
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "    2.1000    1.0500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  1  2  1  0  0  0  0\n" ++
    "  2  3  2  0  0  0  0\n" ++
    "M  END\n$$$$\n";

/// Converts one molecule and classifies it from its own bond topology.
/// Caller owns the returned input.
fn classifyOwnTopology(molecule: *const SdfMolecule, skip_hydrogens: bool, counts: *TopologyRadiiCounts) !types.AtomInput {
    const allocator = std.testing.allocator;
    const mol_slice: []const SdfMolecule = @as([*]const SdfMolecule, @ptrCast(molecule))[0..1];
    var input = try toAtomInput(allocator, mol_slice, skip_hydrogens);
    errdefer input.deinit();
    var stored = try toStoredComponent(allocator, molecule);
    defer stored.deinit();
    const view = stored.view();
    counts.* = try applyTopologyRadii(&input, &view);
    return input;
}

test "applyTopologyRadii gives a molecule the same radii whatever its title is" {
    const allocator = std.testing.allocator;

    // A title, no title, and titles that are residue names of the built-in
    // tables: amino acid, water, and nucleotides, which have atoms named C2
    const titles = [_][]const u8{ "ethanol", "", "   ", "ALA", "HOH", "A", "G", "DT" };
    const residues = [_][]const u8{ "ethan", "", "", "ALA", "HOH", "A", "G", "DT" };

    for ([_][]const u8{ test_ethanol_v2000, test_ethanol_v3000 }) |ethanol| {
        for (titles, residues) |title, residue| {
            const source = try std.mem.concat(allocator, u8, &.{ title, ethanol["ethanol".len..] });
            defer allocator.free(source);
            const molecules = try parse(allocator, source);
            defer freeMolecules(allocator, molecules);
            try std.testing.expectEqual(@as(usize, 1), molecules.len);

            var counts = TopologyRadiiCounts{};

            // Without hydrogens: every atom has a bond-topology radius
            {
                var input = try classifyOwnTopology(&molecules[0], true, &counts);
                defer input.deinit();
                try std.testing.expectEqualSlices(f64, &.{ 1.88, 1.88, 1.46 }, input.r);
                try std.testing.expectEqual(@as(usize, 3), counts.classified);
                try std.testing.expectEqual(@as(usize, 0), counts.fallback);
                // The residue name is the title, and is left alone
                for (input.residue.?) |*name| try std.testing.expectEqualStrings(residue, name.slice());
            }
            // With hydrogens: they get the element radius
            {
                var input = try classifyOwnTopology(&molecules[0], false, &counts);
                defer input.deinit();
                try std.testing.expectEqualSlices(f64, &.{ 1.88, 1.88, 1.46, 1.10, 1.10, 1.10, 1.10, 1.10, 1.10 }, input.r);
                try std.testing.expectEqual(@as(usize, 3), counts.classified);
                try std.testing.expectEqual(@as(usize, 6), counts.fallback);
            }
        }
    }
}

test "applyTopologyRadii uses the bond table of each molecule of a file" {
    const allocator = std.testing.allocator;
    const nitrile_radii = [_]f64{ 1.88, 1.61, 1.64 };
    const aldehyde_radii = [_]f64{ 1.88, 1.76, 1.42 };

    // Two molecules without a title, two with the same title, and one of
    // each: the titles never decide which bond table an atom is read from
    const title_pairs = [_][2][]const u8{ .{ "", "" }, .{ "lig", "lig" }, .{ "nitrile", "" }, .{ "", "aldehyde" } };
    for (title_pairs) |pair| {
        const source = try std.mem.concat(allocator, u8, &.{
            pair[0], "\n", test_header_rest, test_acetonitrile_body,
            pair[1], "\n", test_header_rest, test_acetaldehyde_body,
        });
        defer allocator.free(source);
        const molecules = try parse(allocator, source);
        defer freeMolecules(allocator, molecules);
        try std.testing.expectEqual(@as(usize, 2), molecules.len);

        for (molecules, [_][]const f64{ &nitrile_radii, &aldehyde_radii }) |*molecule, expected| {
            var counts = TopologyRadiiCounts{};
            var input = try classifyOwnTopology(molecule, true, &counts);
            defer input.deinit();
            try std.testing.expectEqualSlices(f64, expected, input.r);
            try std.testing.expectEqual(@as(usize, 3), counts.classified);
            try std.testing.expectEqual(@as(usize, 0), counts.fallback);
        }
    }
}

test "applyTopologyRadii falls back to the element where the bond table gives no radius" {
    const allocator = std.testing.allocator;
    // Chloromethane with an atom of an unknown element next to it
    const source = "\n" ++ test_header_rest ++
        "  3  1  0  0  0  0  0  0  0  0999 V2000\n" ++
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
        "    1.7800    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" ++
        "    9.0000    0.0000    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n" ++
        "  1  2  1  0  0  0  0\n" ++
        "M  END\n$$$$\n";
    const molecules = try parse(allocator, source);
    defer freeMolecules(allocator, molecules);

    var counts = TopologyRadiiCounts{};
    var input = try classifyOwnTopology(&molecules[0], true, &counts);
    defer input.deinit();

    // Carbon from the bond table, chlorine from its element, and the unknown
    // element keeps the radius it came with
    try std.testing.expectEqual(@as(f64, 1.88), input.r[0]);
    try std.testing.expectEqual(classifier.guessRadiusFromAtomicNumber(17).?, input.r[1]);
    try std.testing.expectEqual(elem.Element.X.vdwRadius(), input.r[2]);
    try std.testing.expectEqual(@as(usize, 1), counts.classified);
    try std.testing.expectEqual(@as(usize, 1), counts.fallback);
}

//! PDB Parser for extracting atom coordinates.
//!
//! This module provides a PDB format parser focused on extracting
//! atom coordinates for SASA calculation, following FreeSASA's approach.
//!
//! ## PDB Record Format (Fixed Width)
//!
//! ATOM/HETATM records (columns 1-indexed):
//! - 1-6:   Record name (ATOM/HETATM)
//! - 7-11:  Atom serial number
//! - 13-16: Atom name
//! - 17:    Alternate location indicator
//! - 18-20: Residue name
//! - 22:    Chain identifier
//! - 23-26: Residue sequence number
//! - 27:    Insertion code
//! - 31-38: X coordinate
//! - 39-46: Y coordinate
//! - 47-54: Z coordinate
//! - 55-60: Occupancy
//! - 61-66: Temperature factor
//! - 77-78: Element symbol
//!
//! ## Usage
//!
//! ```zig
//! const parser = @import("pdb_parser.zig");
//!
//! var pdb = parser.PdbParser.init(allocator);
//! const input = try pdb.parseFile("structure.pdb");
//! defer input.deinit();
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const input_io = @import("input_io.zig");
const mmap_reader = @import("mmap_reader.zig");
const compressed = @import("compressed.zig");
const types = @import("types.zig");
const AtomInput = types.AtomInput;
pub const InputIoMode = input_io.InputIoMode;

/// Error types for PDB parsing
pub const ParseError = error{
    /// Invalid coordinate value
    InvalidCoordinate,
    /// No atoms found in file
    NoAtomsFound,
    /// Line too short for required field
    LineTooShort,
};

/// PDB Parser
pub const PdbParser = struct {
    allocator: Allocator,
    /// Filter to include only ATOM records (exclude HETATM)
    /// Default: true (matches FreeSASA/RustSASA behavior)
    atom_only: bool = true,
    /// Skip hydrogen atoms
    /// Default: true (matches FreeSASA/RustSASA behavior)
    skip_hydrogens: bool = true,
    /// Filter to include only first alternate location
    first_alt_loc_only: bool = true,
    /// Model number to extract (null = all models)
    model_num: ?u32 = null,
    /// Read only the first model: stop at the first ENDMDL record, or at a
    /// second MODEL record. Unlike `model_num` this does not depend on how the
    /// models are numbered. Files without MODEL records are read in full.
    first_model_only: bool = false,
    /// Chain IDs to include (null = all chains)
    chain_filter: ?[]const []const u8 = null,
    /// Optional output for callers that keep hydrogens (`skip_hydrogens = false`)
    /// but need to know which atoms they are. When set, `parse` replaces the
    /// list contents with one flag per returned atom: true for the atoms that
    /// `skip_hydrogens` would drop (hydrogen and deuterium). The list is owned
    /// by the caller and grown with the parser's allocator.
    hydrogen_flags: ?*std.ArrayListUnmanaged(bool) = null,

    pub fn init(allocator: Allocator) PdbParser {
        return .{ .allocator = allocator };
    }

    /// Parse PDB from a string
    pub fn parse(self: *PdbParser, source: []const u8) !AtomInput {
        // Dynamic arrays for collecting atoms
        var x_list = std.ArrayListUnmanaged(f64).empty;
        defer x_list.deinit(self.allocator);
        var y_list = std.ArrayListUnmanaged(f64).empty;
        defer y_list.deinit(self.allocator);
        var z_list = std.ArrayListUnmanaged(f64).empty;
        defer z_list.deinit(self.allocator);
        var r_list = std.ArrayListUnmanaged(f64).empty;
        defer r_list.deinit(self.allocator);
        var element_list = std.ArrayListUnmanaged(u8).empty;
        defer element_list.deinit(self.allocator);
        var atom_name_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer atom_name_list.deinit(self.allocator);
        var residue_list = std.ArrayListUnmanaged(types.FixedString5).empty;
        defer residue_list.deinit(self.allocator);
        var chain_id_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer chain_id_list.deinit(self.allocator);
        var residue_num_list = std.ArrayListUnmanaged(i32).empty;
        defer residue_num_list.deinit(self.allocator);
        var insertion_code_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer insertion_code_list.deinit(self.allocator);
        var atom_records = std.ArrayListUnmanaged(AtomRecord).empty;
        defer atom_records.deinit(self.allocator);

        // Pre-allocate based on estimated atom count (PDB line ~80 chars)
        const estimated_atoms = source.len / 80;
        try x_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try y_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try z_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try r_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try element_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try atom_name_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try residue_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try chain_id_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try residue_num_list.ensureTotalCapacity(self.allocator, estimated_atoms);
        try insertion_code_list.ensureTotalCapacity(self.allocator, estimated_atoms);

        // Track model filtering. By default all models are included, matching
        // the calc CLI's documented `--model` default.
        var current_model: ?u32 = null;
        var in_target_model = true;
        var seen_model = false;
        const want_hydrogen_flags = self.hydrogen_flags != null;
        if (self.hydrogen_flags) |flags| flags.clearRetainingCapacity();

        // Only consulted for atoms without a usable element column
        var name_alignment = NameAlignment{ .source = source };

        // Parse line by line
        var lines = std.mem.splitScalar(u8, source, '\n');
        while (lines.next()) |line| {
            // Handle MODEL/ENDMDL records
            if (std.mem.startsWith(u8, line, "MODEL")) {
                if (self.first_model_only and seen_model) break;
                seen_model = true;
                current_model = parseModelNumber(line);
                if (self.model_num) |target| {
                    in_target_model = (current_model == target);
                } else {
                    in_target_model = true;
                }
                continue;
            }
            if (std.mem.startsWith(u8, line, "ENDMDL")) {
                if (self.first_model_only) break;
                continue;
            }

            if (!in_target_model) continue;

            // Check for ATOM/HETATM records
            const is_atom = std.mem.startsWith(u8, line, "ATOM  ");
            const is_hetatm = std.mem.startsWith(u8, line, "HETATM");

            if (!is_atom and !is_hetatm) continue;
            if (self.atom_only and is_hetatm) continue;

            // Parse atom record
            var atom = try self.parseAtomRecord(line, &name_alignment) orelse continue;
            atom.model_num = current_model;

            // Hydrogen filtering (also skip deuterium D, an isotope of H)
            if (self.skip_hydrogens or want_hydrogen_flags) {
                var is_hydrogen = atom.element == .H;
                // Check element column for deuterium (element symbol "D" maps to .X)
                if (!is_hydrogen and line.len >= 78) {
                    const elem_sym = std.mem.trim(u8, line[76..78], " ");
                    is_hydrogen = std.mem.eql(u8, elem_sym, "D");
                }
                if (self.skip_hydrogens and is_hydrogen) continue;
                atom.is_hydrogen = is_hydrogen;
            }

            // Chain filtering
            if (self.chain_filter) |chains| {
                var found = false;
                for (chains) |chain| {
                    if (std.mem.eql(u8, chain, atom.chain_id)) {
                        found = true;
                        break;
                    }
                }
                if (!found) continue;
            }

            try atom_records.append(self.allocator, atom);
        }

        for (atom_records.items, 0..) |atom, i| {
            if (!self.shouldKeepAltLoc(atom_records.items, i)) continue;

            try appendAtomRecord(
                self.allocator,
                atom,
                &x_list,
                &y_list,
                &z_list,
                &r_list,
                &element_list,
                &atom_name_list,
                &residue_list,
                &chain_id_list,
                &residue_num_list,
                &insertion_code_list,
            );
            if (self.hydrogen_flags) |flags| try flags.append(self.allocator, atom.is_hydrogen);
        }

        if (x_list.items.len == 0) {
            return ParseError.NoAtomsFound;
        }

        // Convert to owned slices
        return AtomInput{
            .x = try x_list.toOwnedSlice(self.allocator),
            .y = try y_list.toOwnedSlice(self.allocator),
            .z = try z_list.toOwnedSlice(self.allocator),
            .r = try r_list.toOwnedSlice(self.allocator),
            .element = try element_list.toOwnedSlice(self.allocator),
            .atom_name = try atom_name_list.toOwnedSlice(self.allocator),
            .residue = try residue_list.toOwnedSlice(self.allocator),
            .chain_id = try chain_id_list.toOwnedSlice(self.allocator),
            .residue_num = try residue_num_list.toOwnedSlice(self.allocator),
            .insertion_code = try insertion_code_list.toOwnedSlice(self.allocator),
            .allocator = self.allocator,
        };
    }

    /// Parse PDB from a file (handles plain, .gz, and .zst compressed)
    pub fn parseFile(self: *PdbParser, io: std.Io, path: []const u8) !AtomInput {
        return self.parseFileWithInputIo(io, path, .auto);
    }

    pub fn parseFileWithInputIo(self: *PdbParser, io: std.Io, path: []const u8, input_io_mode: InputIoMode) !AtomInput {
        if (compressed.isCompressed(path)) {
            const data = try compressed.read(self.allocator, path);
            defer self.allocator.free(data);
            return self.parse(data);
        }
        switch (input_io_mode.resolve(.mmap)) {
            .mmap => {
                const mapped = try mmap_reader.mmapFile(self.allocator, io, path);
                defer mapped.deinit();
                return self.parse(mapped.data);
            },
            .read => {
                const file = try std.Io.Dir.cwd().openFile(io, path, .{});
                defer file.close(io);
                var read_buf: [64 * 1024]u8 = undefined;
                var reader = file.reader(io, &read_buf);
                const data = try reader.interface.allocRemaining(self.allocator, .unlimited);
                defer self.allocator.free(data);
                return self.parse(data);
            },
            .auto => unreachable,
        }
    }

    /// Parsed atom data
    const AtomRecord = struct {
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        element: elem.Element,
        atom_name: []const u8,
        residue: []const u8,
        chain_id: []const u8,
        residue_num: i32,
        insertion_code: []const u8,
        alt_loc: u8,
        occupancy: f64,
        model_num: ?u32 = null,
        is_hydrogen: bool = false,
    };

    fn appendAtomRecord(
        allocator: Allocator,
        atom: AtomRecord,
        x_list: *std.ArrayListUnmanaged(f64),
        y_list: *std.ArrayListUnmanaged(f64),
        z_list: *std.ArrayListUnmanaged(f64),
        r_list: *std.ArrayListUnmanaged(f64),
        element_list: *std.ArrayListUnmanaged(u8),
        atom_name_list: *std.ArrayListUnmanaged(types.FixedString4),
        residue_list: *std.ArrayListUnmanaged(types.FixedString5),
        chain_id_list: *std.ArrayListUnmanaged(types.FixedString4),
        residue_num_list: *std.ArrayListUnmanaged(i32),
        insertion_code_list: *std.ArrayListUnmanaged(types.FixedString4),
    ) !void {
        try x_list.append(allocator, atom.x);
        try y_list.append(allocator, atom.y);
        try z_list.append(allocator, atom.z);
        try r_list.append(allocator, atom.radius);
        try element_list.append(allocator, atom.element.atomicNumber());
        try atom_name_list.append(allocator, types.FixedString4.fromSlice(atom.atom_name));
        try residue_list.append(allocator, types.FixedString5.fromSlice(atom.residue));
        try chain_id_list.append(allocator, types.FixedString4.fromSlice(atom.chain_id));
        try residue_num_list.append(allocator, atom.residue_num);
        try insertion_code_list.append(allocator, types.FixedString4.fromSlice(atom.insertion_code));
    }

    fn sameAltLocSite(a: AtomRecord, b: AtomRecord) bool {
        return a.model_num == b.model_num and
            a.residue_num == b.residue_num and
            std.mem.eql(u8, a.chain_id, b.chain_id) and
            std.mem.eql(u8, a.residue, b.residue) and
            std.mem.eql(u8, a.insertion_code, b.insertion_code) and
            std.mem.eql(u8, a.atom_name, b.atom_name);
    }

    fn shouldKeepAltLoc(self: *PdbParser, atoms: []const AtomRecord, index: usize) bool {
        if (!self.first_alt_loc_only) return true;

        const atom = atoms[index];
        if (atom.alt_loc == ' ') return true;

        var best_non_preferred: ?usize = null;
        for (atoms, 0..) |other, other_index| {
            if (!sameAltLocSite(atom, other)) continue;
            if (other.alt_loc == ' ') return false;
            if (other.alt_loc == 'A') return atom.alt_loc == 'A';
            if (best_non_preferred) |best_index| {
                if (other.occupancy > atoms[best_index].occupancy) {
                    best_non_preferred = other_index;
                }
            } else {
                best_non_preferred = other_index;
            }
        }
        return best_non_preferred == index;
    }

    /// Parse a single ATOM/HETATM record
    fn parseAtomRecord(self: *PdbParser, line: []const u8, name_alignment: *NameAlignment) !?AtomRecord {
        _ = self;

        // Minimum line length for coordinates (column 54)
        if (line.len < 54) return null;

        // Extract coordinates (columns 31-38, 39-46, 47-54, 0-indexed: 30-38, 38-46, 46-54)
        const x = parseCoordinate(line[30..38]) orelse return null;
        const y = parseCoordinate(line[38..46]) orelse return null;
        const z = parseCoordinate(line[46..54]) orelse return null;

        const atom_name_raw = if (line.len >= 16) line[12..16] else "    ";
        const atom_name = std.mem.trim(u8, atom_name_raw, " ");
        const residue_raw = if (line.len >= 20) line[17..20] else "   ";
        const residue = std.mem.trim(u8, residue_raw, " ");

        // Extract element (try columns 77-78 first, then infer from atom name)
        const element_field = if (line.len >= 78) line[76..78] else "";
        const element = parseElementField(element_field) orelse
            inferElementFromAtomName(atom_name_raw, residue, name_alignment.isColumnAligned());

        const radius = element.vdwRadius();

        // Extract other fields
        const alt_loc: u8 = if (line.len > 16) line[16] else ' ';
        const occupancy = if (line.len >= 60)
            std.fmt.parseFloat(f64, std.mem.trim(u8, line[54..60], " ")) catch 0.0
        else
            0.0;

        // Chain ID (column 22, 0-indexed 21) - return slice into line
        const chain_id: []const u8 = if (line.len > 21 and line[21] != ' ')
            line[21..22]
        else
            "";

        // Residue number (columns 23-26)
        const res_num_str = if (line.len >= 26) std.mem.trim(u8, line[22..26], " ") else "";
        const residue_num = std.fmt.parseInt(i32, res_num_str, 10) catch 0;

        // Insertion code (column 27, 0-indexed 26) - return slice into line
        const insertion_code: []const u8 = if (line.len > 26 and line[26] != ' ')
            line[26..27]
        else
            "";

        return AtomRecord{
            .x = x,
            .y = y,
            .z = z,
            .radius = radius,
            .element = element,
            .atom_name = atom_name,
            .residue = residue,
            .chain_id = chain_id,
            .residue_num = residue_num,
            .insertion_code = insertion_code,
            .alt_loc = alt_loc,
            .occupancy = occupancy,
        };
    }
};

/// Parse a coordinate value from a fixed-width field
/// Fast implementation avoiding std.fmt.parseFloat overhead
fn parseCoordinate(field: []const u8) ?f64 {
    const len = field.len;
    if (len == 0) return null;

    // Skip leading whitespace
    var start: usize = 0;
    while (start < len and field[start] == ' ') : (start += 1) {}
    if (start == len) return null;

    // Check for negative sign
    var negative = false;
    if (field[start] == '-') {
        negative = true;
        start += 1;
    } else if (field[start] == '+') {
        start += 1;
    }

    // Parse integer part with overflow detection
    var int_part: i64 = 0;
    var has_int_digits = false;
    while (start < len and field[start] >= '0' and field[start] <= '9') : (start += 1) {
        has_int_digits = true;
        const mul_result = @mulWithOverflow(int_part, 10);
        if (mul_result[1] != 0) return null; // Overflow
        const add_result = @addWithOverflow(mul_result[0], @as(i64, field[start] - '0'));
        if (add_result[1] != 0) return null; // Overflow
        int_part = add_result[0];
    }

    // Parse fractional part
    var frac: f64 = 0;
    var has_frac_digits = false;
    if (start < len and field[start] == '.') {
        start += 1;
        var mult: f64 = 0.1;
        while (start < len and field[start] >= '0' and field[start] <= '9') : (start += 1) {
            has_frac_digits = true;
            frac += @as(f64, @floatFromInt(field[start] - '0')) * mult;
            mult *= 0.1;
        }
    }

    // Reject sign-only input (e.g., "-" or "+")
    if (!has_int_digits and !has_frac_digits) return null;

    const result = @as(f64, @floatFromInt(int_part)) + frac;
    return if (negative) -result else result;
}

/// Parse MODEL record to get model number
fn parseModelNumber(line: []const u8) ?u32 {
    // MODEL record: columns 11-14 contain model serial number
    if (line.len < 14) return null;
    const num_str = std.mem.trim(u8, line[10..14], " ");
    return std.fmt.parseInt(u32, num_str, 10) catch null;
}

/// Read the element symbol field (columns 77-78).
///
/// Returns null when the field is blank or is not an element symbol, so the
/// caller infers the element from the atom name instead. Files written before
/// the element column existed can carry other text there (an ID code and line
/// number in columns 73-80).
fn parseElementField(field: []const u8) ?elem.Element {
    const symbol = std.mem.trim(u8, field, " ");
    if (elem.fromSymbolExact(symbol)) |element| return element;

    // "D" (deuterium) and "X" (unknown atom) are valid entries that have no
    // Element of their own
    if (symbol.len == 1) {
        const c = std.ascii.toUpper(symbol[0]);
        if (c == 'D' or c == 'X') return .X;
    }
    return null;
}

/// Element symbol for name-based inference. Transactinides never occur in
/// structures, while their symbols collide with common atom names (SG, NH,
/// DB, DS, HS, CN).
fn elementFromNamePrefix(prefix: []const u8) ?elem.Element {
    const element = elem.fromSymbolExact(prefix) orelse return null;
    return if (element.atomicNumber() >= elem.Element.Rf.atomicNumber()) null else element;
}

/// Elements whose atom names carry remoteness and branch suffixes in PDB
/// nomenclature ("CA", "CD1", "HG21", "NE2", "OG1", "PB", "SD").
/// Same rule as `classifier.extractElementInResidue`.
fn hasSuffixedAtomNames(first_char: u8) bool {
    return switch (std.ascii.toUpper(first_char)) {
        'H', 'C', 'N', 'O', 'P', 'S' => true,
        else => false,
    };
}

/// Whether the atom names of a PDB source follow the column rule, in which an
/// element symbol is right-justified in columns 13-14: " CA " is an alpha
/// carbon and "CA  " is calcium. Files that left-justify or center their names
/// break the rule, and their columns must not be read that way.
///
/// The whole source is scanned once, on first use, so files with an element
/// column never pay for it.
const NameAlignment = struct {
    source: []const u8,
    column_aligned: ?bool = null,

    fn isColumnAligned(self: *NameAlignment) bool {
        if (self.column_aligned == null) {
            self.column_aligned = atomNamesAreColumnAligned(self.source);
        }
        return self.column_aligned.?;
    }
};

/// A name of up to three characters that starts in column 13 must begin with
/// a two-letter element symbol. One that does not ("N   ", "CB  ", "OG1 ")
/// shows that the file does not follow the column rule.
fn atomNamesAreColumnAligned(source: []const u8) bool {
    var lines = std.mem.splitScalar(u8, source, '\n');
    while (lines.next()) |line| {
        if (line.len < 16) continue;
        if (!std.mem.startsWith(u8, line, "ATOM  ") and !std.mem.startsWith(u8, line, "HETATM")) continue;

        const name = line[12..16];
        // Names from column 14 and four-character names say nothing
        if (!std.ascii.isAlphabetic(name[0]) or name[3] != ' ') continue;
        if (elementFromNamePrefix(name[0..2]) == null) return false;
    }
    return true;
}

/// Infer element from the PDB atom name field (columns 13-16) for an atom
/// without a usable element column.
///
/// - A monatomic ion is a residue named after its atom: "CA" in "CA" is
///   calcium, "HG" in "HG" is mercury, "Na+" in "Na+" is sodium.
/// - In a column-aligned file (see `NameAlignment`) the columns decide:
///   " CA ", " NA " and "1HB " are one-letter elements, while "CA  ", "FE  "
///   and "CL1 " start with a two-letter element.
/// - Four-character names fill column 13 whatever their element ("HG21",
///   "HO5'"), and the names of a file that is not column-aligned carry no
///   column information. For both, a name starting with H, C, N, O, P or S is
///   that element, and any other name is a two-letter element if there is one.
fn inferElementFromAtomName(name_field: []const u8, residue: []const u8, column_aligned: bool) elem.Element {
    const trimmed = std.mem.trim(u8, name_field, " ");

    if (std.ascii.eqlIgnoreCase(trimmed, std.mem.trim(u8, residue, " "))) {
        // Charge or oxidation state suffix: "Na+", "Cl-", "FE2"
        const symbol = std.mem.trimEnd(u8, trimmed, "+-0123456789");
        if (elementFromNamePrefix(symbol)) |element| return element;
    }

    // Old-style hydrogen names carry a leading digit ("1HB ", "2HG1")
    const name = std.mem.trimStart(u8, trimmed, "0123456789");
    if (name.len == 0) return .X;

    const one_letter = elem.fromSymbolExact(name[0..1]) orelse .X;
    const two_letter = if (name.len >= 2) elementFromNamePrefix(name[0..2]) else null;

    const fills_all_columns = name_field.len >= 4 and name_field[3] != ' ';
    if (column_aligned and !fills_all_columns) {
        const starts_in_column_13 = std.ascii.isAlphabetic(name_field[0]);
        return if (starts_in_column_13) two_letter orelse one_letter else one_letter;
    }

    if (hasSuffixedAtomNames(name[0])) return one_letter;
    return two_letter orelse one_letter;
}

// Tests
test "parseCoordinate" {
    const testing = std.testing;

    // Basic cases
    try testing.expectEqual(@as(?f64, 11.104), parseCoordinate("  11.104"));
    try testing.expectEqual(@as(?f64, -6.504), parseCoordinate("  -6.504"));
    try testing.expectEqual(@as(?f64, 0.0), parseCoordinate("   0.000"));
    try testing.expectEqual(@as(?f64, null), parseCoordinate("        "));

    // Positive sign
    try testing.expectEqual(@as(?f64, 12.34), parseCoordinate("  +12.34"));

    // Sign-only input (should be null)
    try testing.expectEqual(@as(?f64, null), parseCoordinate("   -   "));
    try testing.expectEqual(@as(?f64, null), parseCoordinate("   +   "));

    // Decimal point only with digits
    try testing.expectEqual(@as(?f64, 0.5), parseCoordinate("     .5"));

    // Large numbers (PDB typical range)
    try testing.expect(parseCoordinate(" 9999.99") != null);
    try testing.expect(parseCoordinate("-9999.99") != null);
}

test "parseElementField" {
    const testing = std.testing;
    const E = elem.Element;

    try testing.expectEqual(@as(?E, .C), parseElementField(" C"));
    try testing.expectEqual(@as(?E, .C), parseElementField("C "));
    try testing.expectEqual(@as(?E, .Fe), parseElementField("FE"));
    try testing.expectEqual(@as(?E, .Fe), parseElementField("Fe"));
    try testing.expectEqual(@as(?E, .Hg), parseElementField("HG"));

    // Deuterium and unknown atoms keep their unknown element
    try testing.expectEqual(@as(?E, .X), parseElementField(" D"));
    try testing.expectEqual(@as(?E, .X), parseElementField(" X"));

    // Blank or not an element symbol: the atom name decides
    try testing.expectEqual(@as(?E, null), parseElementField(""));
    try testing.expectEqual(@as(?E, null), parseElementField("  "));
    try testing.expectEqual(@as(?E, null), parseElementField(" 1"));
    try testing.expectEqual(@as(?E, null), parseElementField("12"));
    try testing.expectEqual(@as(?E, null), parseElementField("C1"));
    try testing.expectEqual(@as(?E, null), parseElementField("1+"));
    try testing.expectEqual(@as(?E, null), parseElementField("QQ"));
}

test "inferElementFromAtomName column-aligned names" {
    const testing = std.testing;
    const E = elem.Element;

    // Names from column 14, and old-style hydrogens with a leading digit
    try testing.expectEqual(E.C, inferElementFromAtomName(" CA ", "ALA", true));
    try testing.expectEqual(E.N, inferElementFromAtomName(" N  ", "ALA", true));
    try testing.expectEqual(E.O, inferElementFromAtomName(" O  ", "ALA", true));
    try testing.expectEqual(E.C, inferElementFromAtomName(" CD1", "LEU", true));
    try testing.expectEqual(E.C, inferElementFromAtomName(" CD ", "PRO", true));
    try testing.expectEqual(E.H, inferElementFromAtomName(" HG ", "SER", true));
    try testing.expectEqual(E.N, inferElementFromAtomName(" NA ", "HEM", true));
    try testing.expectEqual(E.C, inferElementFromAtomName(" CAA", "HEM", true));
    try testing.expectEqual(E.P, inferElementFromAtomName(" PB ", "ATP", true));
    try testing.expectEqual(E.K, inferElementFromAtomName(" K  ", "K", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("1HB ", "ALA", true));

    // Two-letter elements start in column 13
    try testing.expectEqual(E.Fe, inferElementFromAtomName("FE  ", "HEM", true));
    try testing.expectEqual(E.Zn, inferElementFromAtomName("ZN  ", "ZN", true));
    try testing.expectEqual(E.Ca, inferElementFromAtomName("CA  ", "CA", true));
    try testing.expectEqual(E.Na, inferElementFromAtomName("NA  ", "NA", true));
    try testing.expectEqual(E.Cl, inferElementFromAtomName("CL  ", "CL", true));
    try testing.expectEqual(E.Mg, inferElementFromAtomName("MG  ", "MG", true));
    try testing.expectEqual(E.Mn, inferElementFromAtomName("MN  ", "MN", true));
    try testing.expectEqual(E.Cu, inferElementFromAtomName("CU  ", "CU", true));
    try testing.expectEqual(E.Cd, inferElementFromAtomName("CD  ", "CD", true));
    try testing.expectEqual(E.Br, inferElementFromAtomName("BR  ", "BR", true));
    try testing.expectEqual(E.Hg, inferElementFromAtomName("HG  ", "HG", true));
    try testing.expectEqual(E.Ho, inferElementFromAtomName("HO  ", "HO", true));
    // ... also inside a larger residue
    try testing.expectEqual(E.Se, inferElementFromAtomName("SE  ", "MSE", true));
    try testing.expectEqual(E.Hg, inferElementFromAtomName("HG  ", "MMC", true));
    try testing.expectEqual(E.Cu, inferElementFromAtomName("CU1 ", "CUA", true));
    try testing.expectEqual(E.Cl, inferElementFromAtomName("CL1 ", "LIG", true));
    try testing.expectEqual(E.Br, inferElementFromAtomName("BR1 ", "LIG", true));
    try testing.expectEqual(E.Na, inferElementFromAtomName("Na+ ", "Na+", true));

    // Four-character names fill column 13 whatever their element
    try testing.expectEqual(E.H, inferElementFromAtomName("HG21", "VAL", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("HG12", "ILE", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("HD11", "LEU", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("HE21", "GLN", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("HH11", "ARG", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("HO5'", "A", true));
    try testing.expectEqual(E.H, inferElementFromAtomName("2HG1", "VAL", true));
    try testing.expectEqual(E.C, inferElementFromAtomName("CA1B", "LIG", true));
    try testing.expectEqual(E.N, inferElementFromAtomName("NA1B", "LIG", true));
    try testing.expectEqual(E.Fe, inferElementFromAtomName("FE1A", "LIG", true));

    // An ion is recognized even where a name is misplaced
    try testing.expectEqual(E.Ca, inferElementFromAtomName(" CA ", "CA", true));
    try testing.expectEqual(E.Zn, inferElementFromAtomName(" ZN ", "ZN", true));

    // Unknown
    try testing.expectEqual(E.X, inferElementFromAtomName(" D  ", "ALA", true));
    try testing.expectEqual(E.X, inferElementFromAtomName("    ", "ALA", true));
    try testing.expectEqual(E.X, inferElementFromAtomName(" 12 ", "ALA", true));
}

test "inferElementFromAtomName names without column alignment" {
    const testing = std.testing;
    const E = elem.Element;

    // Left-justified names of ordinary atoms are not metals
    try testing.expectEqual(E.C, inferElementFromAtomName("CA  ", "ALA", false));
    try testing.expectEqual(E.N, inferElementFromAtomName("N   ", "ALA", false));
    try testing.expectEqual(E.C, inferElementFromAtomName("CD1 ", "LEU", false));
    try testing.expectEqual(E.C, inferElementFromAtomName("CE  ", "LYS", false));
    try testing.expectEqual(E.N, inferElementFromAtomName("NE2 ", "HIS", false));
    try testing.expectEqual(E.S, inferElementFromAtomName("SG  ", "CYS", false));
    try testing.expectEqual(E.H, inferElementFromAtomName("HG  ", "SER", false));
    try testing.expectEqual(E.H, inferElementFromAtomName("HE1 ", "HIS", false));
    try testing.expectEqual(E.N, inferElementFromAtomName("NA  ", "HEM", false));
    try testing.expectEqual(E.P, inferElementFromAtomName("PB  ", "ATP", false));

    // Ions and names that cannot be an organic atom
    try testing.expectEqual(E.Ca, inferElementFromAtomName("CA  ", "CA", false));
    try testing.expectEqual(E.Na, inferElementFromAtomName("NA  ", "NA", false));
    try testing.expectEqual(E.Hg, inferElementFromAtomName("HG  ", "HG", false));
    try testing.expectEqual(E.Zn, inferElementFromAtomName("ZN  ", "ZN", false));
    try testing.expectEqual(E.Fe, inferElementFromAtomName("FE  ", "HEM", false));
    try testing.expectEqual(E.Fe, inferElementFromAtomName(" FE ", "HEM", false));
    try testing.expectEqual(E.Br, inferElementFromAtomName("BR1 ", "LIG", false));
}

test "atomNamesAreColumnAligned" {
    const testing = std.testing;

    try testing.expect(atomNamesAreColumnAligned(
        \\ATOM      1  N   SER A   1       5.000   0.000   0.000  1.00 20.00
        \\ATOM      2  CA  SER A   1      10.000   0.000   0.000  1.00 20.00
        \\ATOM      3  OG  SER A   1      15.000   0.000   0.000  1.00 20.00
        \\ATOM      4 HG21 VAL A   2      20.000   0.000   0.000  1.00 20.00
        \\ATOM      5 1HB  ALA A   3      25.000   0.000   0.000  1.00 20.00
        \\HETATM    6 FE   HEM A   4      30.000   0.000   0.000  1.00 20.00
        \\HETATM    7 CL1  LIG A   5      35.000   0.000   0.000  1.00 20.00
        \\END
    ));

    // Left-justified names
    try testing.expect(!atomNamesAreColumnAligned(
        \\ATOM      1 N    SER A   1       5.000   0.000   0.000  1.00 20.00
        \\ATOM      2 CA   SER A   1      10.000   0.000   0.000  1.00 20.00
        \\END
    ));

    // Centered names: three-character names start in column 13
    try testing.expect(!atomNamesAreColumnAligned(
        \\ATOM      1  N   THR A   1       5.000   0.000   0.000  1.00 20.00
        \\ATOM      2  CA  THR A   1      10.000   0.000   0.000  1.00 20.00
        \\ATOM      3 OG1  THR A   1      15.000   0.000   0.000  1.00 20.00
        \\END
    ));

    // SG is a sulfur, not seaborgium
    try testing.expect(!atomNamesAreColumnAligned(
        \\ATOM      1 SG   CYS A   1       5.000   0.000   0.000  1.00 20.00
        \\END
    ));
}

/// Parse `pdb_content` with HETATM records and compare the element symbols.
fn expectParsedElements(pdb_content: []const u8, skip_hydrogens: bool, expected: []const []const u8) !void {
    var parser = PdbParser.init(std.testing.allocator);
    parser.atom_only = false;
    parser.skip_hydrogens = skip_hydrogens;
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try std.testing.expectEqual(expected.len, input.atomCount());
    for (expected, input.element.?, input.atom_name.?, input.r) |symbol, atomic_number, atom_name, radius| {
        const element = elem.fromAtomicNumber(atomic_number);
        std.testing.expectEqualStrings(symbol, element.symbol()) catch |err| {
            std.debug.print("atom name: '{s}'\n", .{atom_name.slice()});
            return err;
        };
        try std.testing.expectEqual(element.vdwRadius(), radius);
    }
}

test "PdbParser infers ions and ligand atoms without an element column" {
    const pdb_content =
        \\HETATM    1 NA    NA A   1       5.000   0.000   0.000  1.00 20.00
        \\HETATM    2 CL    CL A   2      10.000   0.000   0.000  1.00 20.00
        \\HETATM    3 ZN    ZN A   3      15.000   0.000   0.000  1.00 20.00
        \\HETATM    4 MG    MG A   4      20.000   0.000   0.000  1.00 20.00
        \\HETATM    5 MN    MN A   5      25.000   0.000   0.000  1.00 20.00
        \\HETATM    6 CU    CU A   6      30.000   0.000   0.000  1.00 20.00
        \\HETATM    7 CD    CD A   7      35.000   0.000   0.000  1.00 20.00
        \\HETATM    8 CA    CA A   8      40.000   0.000   0.000  1.00 20.00
        \\HETATM    9  K     K A   9      45.000   0.000   0.000  1.00 20.00
        \\HETATM   10 FE   HEM A  10      50.000   0.000   0.000  1.00 20.00
        \\HETATM   11  NA  HEM A  10      55.000   0.000   0.000  1.00 20.00
        \\HETATM   12  CAA HEM A  10      60.000   0.000   0.000  1.00 20.00
        \\HETATM   13 SE   MSE A  11      65.000   0.000   0.000  1.00 20.00
        \\HETATM   14  CA  MSE A  11      70.000   0.000   0.000  1.00 20.00
        \\HETATM   15  PB  ATP A  12      75.000   0.000   0.000  1.00 20.00
        \\HETATM   16 BR1  LIG A  13      80.000   0.000   0.000  1.00 20.00
        \\HETATM   17 CL1  LIG A  13      85.000   0.000   0.000  1.00 20.00
        \\HETATM   18  CD1 LIG A  13      90.000   0.000   0.000  1.00 20.00
        \\END
    ;
    try expectParsedElements(pdb_content, true, &.{
        "Na", "Cl", "Zn", "Mg", "Mn", "Cu", "Cd", "Ca", "K",
        "Fe", "N",  "C",  "Se", "C",  "P",  "Br", "Cl", "C",
    });
}

test "PdbParser hydrogen filter without an element column keeps mercury and holmium" {
    const pdb_content =
        \\ATOM      1  N   SER A   1       5.000   0.000   0.000  1.00 20.00
        \\ATOM      2  H   SER A   1      10.000   0.000   0.000  1.00 20.00
        \\ATOM      3  CA  SER A   1      15.000   0.000   0.000  1.00 20.00
        \\ATOM      4  HA  SER A   1      20.000   0.000   0.000  1.00 20.00
        \\ATOM      5  HG  SER A   1      25.000   0.000   0.000  1.00 20.00
        \\ATOM      6 HG21 VAL A   2      30.000   0.000   0.000  1.00 20.00
        \\ATOM      7 1HB  ALA A   3      35.000   0.000   0.000  1.00 20.00
        \\ATOM      8 2HG1 VAL A   2      40.000   0.000   0.000  1.00 20.00
        \\ATOM      9 HO5'   A A   4      45.000   0.000   0.000  1.00 20.00
        \\HETATM   10 HG    HG A   5      50.000   0.000   0.000  1.00 20.00
        \\HETATM   11 HO    HO A   6      55.000   0.000   0.000  1.00 20.00
        \\HETATM   12 HG   MMC A   7      60.000   0.000   0.000  1.00 20.00
        \\END
    ;
    // Default: the seven hydrogens are removed, the metals stay
    try expectParsedElements(pdb_content, true, &.{ "N", "C", "Hg", "Ho", "Hg" });
    try expectParsedElements(pdb_content, false, &.{
        "N", "H", "C", "H", "H", "H", "H", "H", "H", "Hg", "Ho", "Hg",
    });
}

test "PdbParser ignores an ID code and line number in the element columns" {
    // Columns 73-80 of files written before the element column existed
    const pdb_content =
        \\ATOM      1  N   SER A   1       5.000   0.000   0.000  1.00 20.00      1ABC 101
        \\ATOM      2  CA  SER A   1      10.000   0.000   0.000  1.00 20.00      1ABC 102
        \\ATOM      3  HA  SER A   1      15.000   0.000   0.000  1.00 20.00      1ABC 103
        \\ATOM      4  OG  SER A   1      20.000   0.000   0.000  1.00 20.00      1ABC 104
        \\ATOM      5  HG  SER A   1      25.000   0.000   0.000  1.00 20.00      1ABC 105
        \\HETATM    6 ZN    ZN A   2      30.000   0.000   0.000  1.00 20.00      1ABC1106
        \\HETATM    7 FE   HEM A   3      35.000   0.000   0.000  1.00 20.00      1ABC 107
        \\END
    ;
    try expectParsedElements(pdb_content, true, &.{ "N", "C", "O", "Zn", "Fe" });
    try expectParsedElements(pdb_content, false, &.{ "N", "C", "H", "O", "H", "Zn", "Fe" });
}

test "PdbParser left-justified names without an element column" {
    const pdb_content =
        \\ATOM      1 N    SER A   1       5.000   0.000   0.000  1.00 20.00
        \\ATOM      2 CA   SER A   1      10.000   0.000   0.000  1.00 20.00
        \\ATOM      3 HA   SER A   1      15.000   0.000   0.000  1.00 20.00
        \\ATOM      4 OG   SER A   1      20.000   0.000   0.000  1.00 20.00
        \\ATOM      5 HG   SER A   1      25.000   0.000   0.000  1.00 20.00
        \\ATOM      6 CD1  LEU A   2      30.000   0.000   0.000  1.00 20.00
        \\ATOM      7 NE2  HIS A   3      35.000   0.000   0.000  1.00 20.00
        \\ATOM      8 SG   CYS A   4      40.000   0.000   0.000  1.00 20.00
        \\HETATM    9 CA    CA A   5      45.000   0.000   0.000  1.00 20.00
        \\HETATM   10 FE   HEM A   6      50.000   0.000   0.000  1.00 20.00
        \\HETATM   11 NA   HEM A   6      55.000   0.000   0.000  1.00 20.00
        \\HETATM   12 NA    NA A   7      60.000   0.000   0.000  1.00 20.00
        \\END
    ;
    // Alpha carbon, gamma hydrogen and heme nitrogen, not calcium, mercury and sodium
    try expectParsedElements(pdb_content, false, &.{
        "N", "C", "H", "O", "H", "C", "N", "S", "Ca", "Fe", "N", "Na",
    });
    try expectParsedElements(pdb_content, true, &.{
        "N", "C", "O", "C", "N", "S", "Ca", "Fe", "N", "Na",
    });
}

test "PdbParser element column decides whatever the atom name" {
    // Names that would be read differently without the element column, in
    // both alignments, and the deuterium and unknown symbols
    const pdb_content =
        \\ATOM      1 CA   ALA A   1       5.000   0.000   0.000  1.00 20.00           C
        \\ATOM      2 HG   SER A   2      10.000   0.000   0.000  1.00 20.00           H
        \\HETATM    3  CA   CA A   3      15.000   0.000   0.000  1.00 20.00          CA
        \\HETATM    4  HG   HG A   4      20.000   0.000   0.000  1.00 20.00          HG
        \\HETATM    5 CL1  LIG A   5      25.000   0.000   0.000  1.00 20.00           C
        \\HETATM    6 FE   HEM A   6      30.000   0.000   0.000  1.00 20.00          Fe
        \\ATOM      7  DA  ALA A   1      35.000   0.000   0.000  1.00 20.00           D
        \\HETATM    8  C1  UNL A   7      40.000   0.000   0.000  1.00 20.00           X
        \\END
    ;
    try expectParsedElements(pdb_content, false, &.{ "C", "H", "Ca", "Hg", "C", "Fe", "X", "X" });
    // Hydrogen and deuterium are removed by default
    try expectParsedElements(pdb_content, true, &.{ "C", "Ca", "Hg", "C", "Fe", "X" });
}

test "PdbParser basic" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N
        \\ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C
        \\ATOM      3  C   ALA A   1      10.480   5.927  -4.153  1.00  7.65           C
        \\HETATM  100  O   HOH A 101       5.000   5.000   5.000  1.00 20.00           O
        \\END
    ;

    // Default: atom_only=true, skip_hydrogens=true (HETATM excluded)
    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 3), input.atomCount());
    try testing.expectApproxEqAbs(@as(f64, 11.104), input.x[0], 0.001);
    try testing.expectApproxEqAbs(@as(f64, 6.134), input.y[0], 0.001);
    try testing.expectApproxEqAbs(@as(f64, -6.504), input.z[0], 0.001);
}

test "PdbParser include HETATM" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N
        \\ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C
        \\ATOM      3  C   ALA A   1      10.480   5.927  -4.153  1.00  7.65           C
        \\HETATM  100  O   HOH A 101       5.000   5.000   5.000  1.00 20.00           O
        \\END
    ;

    var parser = PdbParser.init(allocator);
    parser.atom_only = false;
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 4), input.atomCount());
}

test "PdbParser default model selection includes all models" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\MODEL        1
        \\ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 10.00           C
        \\ENDMDL
        \\MODEL        2
        \\ATOM      2  CA  GLY B   2      11.000  21.000  31.000  1.00 10.00           C
        \\ENDMDL
        \\END
    ;

    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 2), input.atomCount());
    try testing.expectEqualStrings("A", input.chain_id.?[0].slice());
    try testing.expectEqualStrings("B", input.chain_id.?[1].slice());
}

test "PdbParser explicit model selection filters requested model" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\MODEL        1
        \\ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 10.00           C
        \\ENDMDL
        \\MODEL        2
        \\ATOM      2  CA  GLY B   2      11.000  21.000  31.000  1.00 10.00           C
        \\ENDMDL
        \\END
    ;

    var parser = PdbParser.init(allocator);
    parser.model_num = 2;
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 1), input.atomCount());
    try testing.expectEqualStrings("B", input.chain_id.?[0].slice());
}

test "PdbParser first_model_only stops after the first model whatever its number" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\MODEL        7
        \\ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 10.00           C
        \\HETATM    2  O   HOH A   2      12.000  20.000  30.000  1.00 10.00           O
        \\ENDMDL
        \\MODEL        8
        \\ATOM      1  CA  ALA A   1      11.000  21.000  31.000  1.00 10.00           C
        \\HETATM    2  O   HOH A   2      13.000  21.000  31.000  1.00 10.00           O
        \\ENDMDL
        \\END
    ;

    var parser = PdbParser.init(allocator);
    parser.atom_only = false;
    parser.first_model_only = true;
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 2), input.atomCount());
    try testing.expectEqual(@as(f64, 10.0), input.x[0]);
    try testing.expectEqual(@as(f64, 12.0), input.x[1]);

    // A second MODEL record ends the first model even without ENDMDL.
    const unterminated =
        \\MODEL        1
        \\ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 10.00           C
        \\MODEL        2
        \\ATOM      1  CA  ALA A   1      11.000  21.000  31.000  1.00 10.00           C
        \\END
    ;
    var unterminated_input = try parser.parse(unterminated);
    defer unterminated_input.deinit();
    try testing.expectEqual(@as(usize, 1), unterminated_input.atomCount());

    // Without MODEL records the whole file is one model.
    const single =
        \\ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00 10.00           C
        \\ATOM      2  CB  ALA A   1      11.000  21.000  31.000  1.00 10.00           C
        \\END
    ;
    var single_input = try parser.parse(single);
    defer single_input.deinit();
    try testing.expectEqual(@as(usize, 2), single_input.atomCount());
}

test "PdbParser hydrogen_flags marks the atoms skip_hydrogens would drop" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N
        \\ATOM      2  H   ALA A   1      11.500   6.900  -6.900  1.00 10.00           H
        \\ATOM      3  D   ALA A   1      12.000   7.000  -5.000  1.00 10.00           D
        \\ATOM      4  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C
        \\ATOM      5 1HB  ALA A   1      12.639   6.071  -5.147  1.00  9.13
        \\END
    ;

    var flags: std.ArrayListUnmanaged(bool) = .empty;
    defer flags.deinit(allocator);

    var parser = PdbParser.init(allocator);
    parser.skip_hydrogens = false;
    parser.hydrogen_flags = &flags;
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 5), input.atomCount());
    try testing.expectEqualSlices(bool, &.{ false, true, true, false, true }, flags.items);

    // The flags describe the returned atoms, so a second parse replaces them.
    parser.skip_hydrogens = true;
    var heavy = try parser.parse(pdb_content);
    defer heavy.deinit();
    try testing.expectEqual(@as(usize, 2), heavy.atomCount());
    try testing.expectEqualSlices(bool, &.{ false, false }, flags.items);
}

test "PdbParser altLoc selection is per atom site and keeps later B-only sites" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  CA AALA A   1      10.000  20.000  30.000  0.60 10.00           C
        \\ATOM      2  CA BALA A   1      12.000  22.000  32.000  0.40 10.00           C
        \\ATOM      3  CA BGLY A   2      14.000  24.000  34.000  0.50 10.00           C
        \\END
    ;

    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 2), input.atomCount());
    try testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try testing.expectApproxEqAbs(@as(f64, 14.0), input.x[1], 0.001);
    try testing.expectEqual(@as(i32, 2), input.residue_num.?[1]);
}

test "PdbParser altLoc selection falls back to highest occupancy without blank or A" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  CA BALA A   1      10.000  20.000  30.000  0.40 10.00           C
        \\ATOM      2  CA CALA A   1      12.000  22.000  32.000  0.70 10.00           C
        \\END
    ;

    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 1), input.atomCount());
    try testing.expectApproxEqAbs(@as(f64, 12.0), input.x[0], 0.001);
}

test "PdbParser altLoc selection is scoped by model" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\MODEL        1
        \\ATOM      1  CA AALA A   1      10.000  20.000  30.000  0.60 10.00           C
        \\ATOM      2  CA BALA A   1      12.000  22.000  32.000  0.40 10.00           C
        \\ENDMDL
        \\MODEL        2
        \\ATOM      3  CA BALA A   1      14.000  24.000  34.000  0.50 10.00           C
        \\ENDMDL
        \\END
    ;

    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 2), input.atomCount());
    try testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try testing.expectApproxEqAbs(@as(f64, 14.0), input.x[1], 0.001);
}

test "PdbParser atom_only filter (default)" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N
        \\HETATM  100  O   HOH A 101       5.000   5.000   5.000  1.00 20.00           O
        \\END
    ;

    // Default atom_only=true excludes HETATM
    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 1), input.atomCount());
}

test "PdbParser skip_hydrogens filter" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N
        \\ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C
        \\ATOM      3 1HB  ALA A   1      12.000   7.000  -5.000  1.00 10.00           H
        \\ATOM      4  H   ALA A   1      10.500   6.500  -7.000  1.00 12.00           H
        \\ATOM      5  O   ALA A   1      10.480   5.927  -4.153  1.00  7.65           O
        \\END
    ;

    // Default skip_hydrogens=true: should exclude H atoms
    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 3), input.atomCount()); // N, CA, O

    // Include hydrogens
    var parser2 = PdbParser.init(allocator);
    parser2.skip_hydrogens = false;
    var input2 = try parser2.parse(pdb_content);
    defer input2.deinit();

    try testing.expectEqual(@as(usize, 5), input2.atomCount()); // All atoms
}

test "PdbParser deuterium filter" {
    const testing = std.testing;
    const allocator = testing.allocator;

    const pdb_content =
        \\ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N
        \\ATOM      2  D   ALA A   1      12.000   7.000  -5.000  1.00 10.00           D
        \\ATOM      3  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C
        \\END
    ;

    // Default skip_hydrogens=true should also skip deuterium
    var parser = PdbParser.init(allocator);
    var input = try parser.parse(pdb_content);
    defer input.deinit();

    try testing.expectEqual(@as(usize, 2), input.atomCount()); // N, CA only
}

test "fuzz pdb parser" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, smith: *std.testing.Smith) !void {
            const input = smith.in orelse return;
            var parser = PdbParser.init(std.testing.allocator);
            var result = parser.parse(input) catch return;
            result.deinit();
        }
    }.testOne, .{
        .corpus = &.{
            "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N\n",
            "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N\nATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C\n",
            "HETATM  100  O   HOH A 101       5.000   5.000   5.000  1.00 20.00           O\n",
            "MODEL        1\nATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N\nENDMDL\n",
        },
    });
}

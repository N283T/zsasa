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
const mmap_reader = @import("mmap_reader.zig");
const types = @import("types.zig");
const AtomInput = types.AtomInput;

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
    /// Model number to extract (null = first model)
    model_num: ?u32 = null,
    /// Chain IDs to include (null = all chains)
    chain_filter: ?[]const []const u8 = null,

    pub fn init(allocator: Allocator) PdbParser {
        return .{ .allocator = allocator };
    }

    /// Parse PDB from a string
    pub fn parse(self: *PdbParser, source: []const u8) !AtomInput {
        // Dynamic arrays for collecting atoms
        var x_list = std.ArrayListUnmanaged(f64){};
        defer x_list.deinit(self.allocator);
        var y_list = std.ArrayListUnmanaged(f64){};
        defer y_list.deinit(self.allocator);
        var z_list = std.ArrayListUnmanaged(f64){};
        defer z_list.deinit(self.allocator);
        var r_list = std.ArrayListUnmanaged(f64){};
        defer r_list.deinit(self.allocator);
        var element_list = std.ArrayListUnmanaged(u8){};
        defer element_list.deinit(self.allocator);
        var atom_name_list = std.ArrayListUnmanaged(types.FixedString4){};
        defer atom_name_list.deinit(self.allocator);
        var residue_list = std.ArrayListUnmanaged(types.FixedString5){};
        defer residue_list.deinit(self.allocator);
        var chain_id_list = std.ArrayListUnmanaged(types.FixedString4){};
        defer chain_id_list.deinit(self.allocator);
        var residue_num_list = std.ArrayListUnmanaged(i32){};
        defer residue_num_list.deinit(self.allocator);
        var insertion_code_list = std.ArrayListUnmanaged(types.FixedString4){};
        defer insertion_code_list.deinit(self.allocator);

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

        // Track alt location for filtering
        var first_alt_loc: u8 = ' ';
        var current_model: ?u32 = null;
        var in_target_model = true;

        // Parse line by line
        var lines = std.mem.splitScalar(u8, source, '\n');
        while (lines.next()) |line| {
            // Handle MODEL/ENDMDL records
            if (std.mem.startsWith(u8, line, "MODEL")) {
                current_model = parseModelNumber(line);
                if (self.model_num) |target| {
                    in_target_model = (current_model == target);
                } else {
                    // No model specified: use first model only
                    in_target_model = (current_model == 1 or current_model == null);
                }
                continue;
            }
            if (std.mem.startsWith(u8, line, "ENDMDL")) {
                if (self.model_num == null and current_model != null) {
                    // First model complete, stop parsing
                    break;
                }
                continue;
            }

            if (!in_target_model) continue;

            // Check for ATOM/HETATM records
            const is_atom = std.mem.startsWith(u8, line, "ATOM  ");
            const is_hetatm = std.mem.startsWith(u8, line, "HETATM");

            if (!is_atom and !is_hetatm) continue;
            if (self.atom_only and is_hetatm) continue;

            // Parse atom record
            const atom = try self.parseAtomRecord(line) orelse continue;

            // Hydrogen filtering (also skip deuterium D, an isotope of H)
            if (self.skip_hydrogens) {
                if (atom.element == .H) continue;
                // Check element column for deuterium (element symbol "D" maps to .X)
                if (line.len >= 78) {
                    const elem_sym = std.mem.trim(u8, line[76..78], " ");
                    if (std.mem.eql(u8, elem_sym, "D")) continue;
                }
            }

            // Alt location filtering
            if (self.first_alt_loc_only) {
                if (atom.alt_loc != ' ') {
                    if (first_alt_loc == ' ') {
                        first_alt_loc = atom.alt_loc;
                    } else if (atom.alt_loc != first_alt_loc) {
                        continue;
                    }
                }
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

            // Add atom data
            try x_list.append(self.allocator, atom.x);
            try y_list.append(self.allocator, atom.y);
            try z_list.append(self.allocator, atom.z);
            try r_list.append(self.allocator, atom.radius);
            try element_list.append(self.allocator, atom.element.atomicNumber());

            // Use FixedString4 - no per-atom allocation needed
            try atom_name_list.append(self.allocator, types.FixedString4.fromSlice(atom.atom_name));
            try residue_list.append(self.allocator, types.FixedString5.fromSlice(atom.residue));
            try chain_id_list.append(self.allocator, types.FixedString4.fromSlice(atom.chain_id));
            try residue_num_list.append(self.allocator, atom.residue_num);
            try insertion_code_list.append(self.allocator, types.FixedString4.fromSlice(atom.insertion_code));
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

    /// Parse PDB from a file
    pub fn parseFile(self: *PdbParser, path: []const u8) !AtomInput {
        const mapped = try mmap_reader.mmapFile(path);
        defer mapped.deinit();
        return self.parse(mapped.data);
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
    };

    /// Parse a single ATOM/HETATM record
    fn parseAtomRecord(self: *PdbParser, line: []const u8) !?AtomRecord {
        _ = self;

        // Minimum line length for coordinates (column 54)
        if (line.len < 54) return null;

        // Extract coordinates (columns 31-38, 39-46, 47-54, 0-indexed: 30-38, 38-46, 46-54)
        const x = parseCoordinate(line[30..38]) orelse return null;
        const y = parseCoordinate(line[38..46]) orelse return null;
        const z = parseCoordinate(line[46..54]) orelse return null;

        // Extract element (try columns 77-78 first, then infer from atom name)
        const element_symbol = if (line.len >= 78)
            std.mem.trim(u8, line[76..78], " ")
        else
            "";

        const atom_name_raw = if (line.len >= 16) line[12..16] else "    ";
        const atom_name = std.mem.trim(u8, atom_name_raw, " ");

        const element = if (element_symbol.len > 0)
            elem.fromSymbol(element_symbol)
        else
            inferElementFromAtomName(atom_name_raw);

        const radius = element.vdwRadius();

        // Extract other fields
        const alt_loc: u8 = if (line.len > 16) line[16] else ' ';
        const residue_raw = if (line.len >= 20) line[17..20] else "   ";
        const residue = std.mem.trim(u8, residue_raw, " ");

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

/// Infer element from PDB atom name field (columns 13-16)
/// Following FreeSASA's approach:
/// - Position 13-14 (0-indexed 12-13): element symbol for standard atoms
/// - First letter after leading digit/space is typically the element
fn inferElementFromAtomName(atom_name: []const u8) elem.Element {
    if (atom_name.len < 2) return .X;

    // Check first character
    const first = atom_name[0];
    const second = atom_name[1];

    // Standard PDB: element symbol is at positions 13-14 (right-justified for 1-char elements)
    // " CA " -> C (alpha carbon)
    // " N  " -> N
    // "FE  " -> Fe (iron)
    // "1HB " -> H (hydrogen with digit prefix)

    if (first == ' ' or (first >= '0' and first <= '9')) {
        // Single-letter element at position 2
        return elem.fromSymbol(&[_]u8{second});
    }

    // Two-letter element (e.g., FE, CA for calcium in HETATM)
    // But be careful: CA in ATOM records is carbon-alpha, not calcium
    // For safety, just use first letter
    return elem.fromSymbol(&[_]u8{first});
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

test "inferElementFromAtomName" {
    const testing = std.testing;

    try testing.expectEqual(elem.Element.C, inferElementFromAtomName(" CA "));
    try testing.expectEqual(elem.Element.N, inferElementFromAtomName(" N  "));
    try testing.expectEqual(elem.Element.O, inferElementFromAtomName(" O  "));
    try testing.expectEqual(elem.Element.H, inferElementFromAtomName("1HB "));
    try testing.expectEqual(elem.Element.F, inferElementFromAtomName("FE  ")); // Will be F, not Fe
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
        fn testOne(_: void, input: []const u8) !void {
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

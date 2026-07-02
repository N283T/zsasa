const std = @import("std");
const compressed = @import("compressed.zig");
const elem = @import("element.zig");
const mmap_reader = @import("mmap_reader.zig");
const types = @import("types.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;

pub const ParseError = error{
    MissingSequence,
    NonStandardResidue,
    NoAtomRecord,
    MissingCoordinateMarker,
    InvalidCoordinate,
    AtomCountMismatch,
};

const ResidueSpec = struct {
    one_letter: u8,
    name: []const u8,
    atoms: []const []const u8,
};

const gly_atoms = [_][]const u8{ "N", "CA", "C", "O" };
const ala_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O" };
const val_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG1", "CG2" };
const leu_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2" };
const ile_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG1", "CG2", "CD1" };
const ser_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "OG" };
const thr_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG2", "OG1" };
const cys_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "SG" };
const met_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "SD", "CE" };
const pro_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD" };
const phe_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE1", "CE2", "CZ" };
const tyr_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE1", "CE2", "OH", "CZ" };
const trp_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD1", "CD2", "CE2", "CE3", "NE1", "CH2", "CZ2", "CZ3" };
const his_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD2", "ND1", "CE1", "NE2" };
const glu_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD", "OE1", "OE2" };
const asp_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "OD1", "OD2" };
const asn_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "ND2", "OD1" };
const gln_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD", "NE2", "OE1" };
const lys_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD", "CE", "NZ" };
const arg_atoms = [_][]const u8{ "N", "CA", "C", "CB", "O", "CG", "CD", "NE", "NH1", "NH2", "CZ" };

fn residueSpec(one_letter: u8) ?ResidueSpec {
    return switch (one_letter) {
        'A' => .{ .one_letter = 'A', .name = "ALA", .atoms = &ala_atoms },
        'C' => .{ .one_letter = 'C', .name = "CYS", .atoms = &cys_atoms },
        'D' => .{ .one_letter = 'D', .name = "ASP", .atoms = &asp_atoms },
        'E' => .{ .one_letter = 'E', .name = "GLU", .atoms = &glu_atoms },
        'F' => .{ .one_letter = 'F', .name = "PHE", .atoms = &phe_atoms },
        'G' => .{ .one_letter = 'G', .name = "GLY", .atoms = &gly_atoms },
        'H' => .{ .one_letter = 'H', .name = "HIS", .atoms = &his_atoms },
        'I' => .{ .one_letter = 'I', .name = "ILE", .atoms = &ile_atoms },
        'K' => .{ .one_letter = 'K', .name = "LYS", .atoms = &lys_atoms },
        'L' => .{ .one_letter = 'L', .name = "LEU", .atoms = &leu_atoms },
        'M' => .{ .one_letter = 'M', .name = "MET", .atoms = &met_atoms },
        'N' => .{ .one_letter = 'N', .name = "ASN", .atoms = &asn_atoms },
        'P' => .{ .one_letter = 'P', .name = "PRO", .atoms = &pro_atoms },
        'Q' => .{ .one_letter = 'Q', .name = "GLN", .atoms = &gln_atoms },
        'R' => .{ .one_letter = 'R', .name = "ARG", .atoms = &arg_atoms },
        'S' => .{ .one_letter = 'S', .name = "SER", .atoms = &ser_atoms },
        'T' => .{ .one_letter = 'T', .name = "THR", .atoms = &thr_atoms },
        'V' => .{ .one_letter = 'V', .name = "VAL", .atoms = &val_atoms },
        'W' => .{ .one_letter = 'W', .name = "TRP", .atoms = &trp_atoms },
        'Y' => .{ .one_letter = 'Y', .name = "TYR", .atoms = &tyr_atoms },
        else => null,
    };
}

pub fn parseFile(allocator: Allocator, io: std.Io, path: []const u8) !AtomInput {
    if (compressed.isCompressed(path)) {
        const data = try compressed.read(allocator, path);
        defer allocator.free(data);
        return parse(allocator, data);
    }

    const mapped = try mmap_reader.mmapFile(allocator, io, path);
    defer mapped.deinit();
    return parse(allocator, mapped.data);
}

pub fn parse(allocator: Allocator, source: []const u8) !AtomInput {
    const sequence = try extractSequence(allocator, source);
    defer allocator.free(sequence);
    if (sequence.len == 0) return ParseError.MissingSequence;

    const expected_atoms = try expectedAtomCount(sequence);
    var x = try allocator.alloc(f64, expected_atoms);
    errdefer allocator.free(x);
    var y = try allocator.alloc(f64, expected_atoms);
    errdefer allocator.free(y);
    var z = try allocator.alloc(f64, expected_atoms);
    errdefer allocator.free(z);
    var r = try allocator.alloc(f64, expected_atoms);
    errdefer allocator.free(r);
    var residue = try allocator.alloc(types.FixedString5, expected_atoms);
    errdefer allocator.free(residue);
    var atom_name = try allocator.alloc(types.FixedString4, expected_atoms);
    errdefer allocator.free(atom_name);
    var element = try allocator.alloc(u8, expected_atoms);
    errdefer allocator.free(element);
    var chain_id = try allocator.alloc(types.FixedString4, expected_atoms);
    errdefer allocator.free(chain_id);
    var residue_num = try allocator.alloc(i32, expected_atoms);
    errdefer allocator.free(residue_num);
    var insertion_code = try allocator.alloc(types.FixedString4, expected_atoms);
    errdefer allocator.free(insertion_code);

    var line_it = std.mem.splitScalar(u8, source, '\n');
    var atom_index: usize = 0;
    var cursor = AtomDescriptionCursor{ .sequence = sequence };
    while (line_it.next()) |raw_line| {
        const line = std.mem.trim(u8, raw_line, "\r");
        if (!std.mem.startsWith(u8, line, "ATOM ")) continue;
        if (atom_index >= expected_atoms) return ParseError.AtomCountMismatch;

        const xyz = try parseAtomLineCoordinates(line);
        x[atom_index] = xyz[0];
        y[atom_index] = xyz[1];
        z[atom_index] = xyz[2];

        const desc = try cursor.next();
        residue[atom_index] = types.FixedString5.fromSlice(desc.residue_name);
        atom_name[atom_index] = types.FixedString4.fromSlice(desc.atom_name);
        const e = elementFromAtomName(desc.atom_name);
        element[atom_index] = e.atomicNumber();
        r[atom_index] = e.vdwRadius();
        chain_id[atom_index] = types.FixedString4.fromSlice("A");
        residue_num[atom_index] = @intCast(desc.residue_index + 1);
        insertion_code[atom_index] = types.FixedString4.fromSlice("");

        atom_index += 1;
    }

    if (atom_index == 0) return ParseError.NoAtomRecord;
    if (atom_index != expected_atoms) return ParseError.AtomCountMismatch;

    return .{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .residue = residue,
        .atom_name = atom_name,
        .element = element,
        .chain_id = chain_id,
        .residue_num = residue_num,
        .insertion_code = insertion_code,
        .allocator = allocator,
    };
}

fn extractSequence(allocator: Allocator, source: []const u8) ![]u8 {
    const marker = "_struct_ref.pdbx_seq_one_letter_code";
    const marker_pos = std.mem.indexOf(u8, source, marker) orelse return ParseError.MissingSequence;
    var p: usize = marker_pos + marker.len;
    while (p < source.len and (source[p] == ' ' or source[p] == '\t')) p += 1;

    var raw: []const u8 = "";
    if (p < source.len and source[p] != '\n' and source[p] != '\r') {
        const line_end = std.mem.indexOfScalarPos(u8, source, p, '\n') orelse source.len;
        raw = source[p..line_end];
    } else {
        const newline = std.mem.indexOfScalarPos(u8, source, p, '\n') orelse return ParseError.MissingSequence;
        p = newline + 1;
        if (p >= source.len or source[p] != ';') return ParseError.MissingSequence;
        const content_start = p + 1;
        const content_end = std.mem.indexOf(u8, source[content_start..], "\n;") orelse return ParseError.MissingSequence;
        raw = source[content_start .. content_start + content_end];
    }

    var sequence = std.ArrayListUnmanaged(u8).empty;
    errdefer sequence.deinit(allocator);
    for (raw) |c| {
        if (std.ascii.isWhitespace(c) or c == ';' or c == '\'' or c == '"') continue;
        const aa = std.ascii.toUpper(c);
        if (residueSpec(aa) == null) return ParseError.NonStandardResidue;
        try sequence.append(allocator, aa);
    }
    return sequence.toOwnedSlice(allocator);
}

fn expectedAtomCount(sequence: []const u8) !usize {
    var count: usize = 1; // terminal OXT
    for (sequence) |aa| {
        const spec = residueSpec(aa) orelse return ParseError.NonStandardResidue;
        count += spec.atoms.len;
    }
    return count;
}

fn parseAtomLineCoordinates(line: []const u8) ![3]f64 {
    const q_index = std.mem.indexOfScalar(u8, line, '?') orelse return ParseError.MissingCoordinateMarker;
    var token_it = std.mem.tokenizeAny(u8, line[q_index + 1 ..], " \t\r");
    const x_token = token_it.next() orelse return ParseError.InvalidCoordinate;
    const y_token = token_it.next() orelse return ParseError.InvalidCoordinate;
    const z_token = token_it.next() orelse return ParseError.InvalidCoordinate;
    return .{
        std.fmt.parseFloat(f64, x_token) catch return ParseError.InvalidCoordinate,
        std.fmt.parseFloat(f64, y_token) catch return ParseError.InvalidCoordinate,
        std.fmt.parseFloat(f64, z_token) catch return ParseError.InvalidCoordinate,
    };
}

const AtomDescription = struct {
    residue_index: usize,
    residue_name: []const u8,
    atom_name: []const u8,
};

const AtomDescriptionCursor = struct {
    sequence: []const u8,
    residue_index: usize = 0,
    atom_index_in_residue: usize = 0,
    emitted_terminal_oxt: bool = false,

    fn next(self: *AtomDescriptionCursor) !AtomDescription {
        while (self.residue_index < self.sequence.len) {
            const spec = residueSpec(self.sequence[self.residue_index]) orelse return ParseError.NonStandardResidue;
            if (self.atom_index_in_residue < spec.atoms.len) {
                const desc = AtomDescription{
                    .residue_index = self.residue_index,
                    .residue_name = spec.name,
                    .atom_name = spec.atoms[self.atom_index_in_residue],
                };
                self.atom_index_in_residue += 1;
                return desc;
            }
            self.residue_index += 1;
            self.atom_index_in_residue = 0;
        }

        if (!self.emitted_terminal_oxt and self.sequence.len > 0) {
            self.emitted_terminal_oxt = true;
            const last = residueSpec(self.sequence[self.sequence.len - 1]) orelse return ParseError.NonStandardResidue;
            return .{
                .residue_index = self.sequence.len - 1,
                .residue_name = last.name,
                .atom_name = "OXT",
            };
        }

        return ParseError.AtomCountMismatch;
    }
};

fn elementFromAtomName(atom_name: []const u8) elem.Element {
    for (atom_name) |c| {
        const upper = std.ascii.toUpper(c);
        return switch (upper) {
            'H' => .H,
            'C' => .C,
            'N' => .N,
            'O' => .O,
            'P' => .P,
            'S' => .S,
            else => .X,
        };
    }
    return .X;
}

test "parse AlphaFold-like mmCIF atom records from sequence and ATOM rows" {
    const source =
        \\data_AF_TEST
        \\_struct_ref.pdbx_seq_one_letter_code
        \\;AG
        \\;
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_alt_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_ins_code
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\_atom_site.B_iso_or_equiv
        \\ATOM 1 N N . ALA A 1 1 ? 1.000 2.000 3.000 1.00 90.00
        \\ATOM 2 C CA . ALA A 1 1 ? 2.000 3.000 4.000 1.00 90.00
        \\ATOM 3 C C . ALA A 1 1 ? 3.000 4.000 5.000 1.00 90.00
        \\ATOM 4 C CB . ALA A 1 1 ? 4.000 5.000 6.000 1.00 90.00
        \\ATOM 5 O O . ALA A 1 1 ? 5.000 6.000 7.000 1.00 90.00
        \\ATOM 6 N N . GLY A 1 2 ? 6.000 7.000 8.000 1.00 80.00
        \\ATOM 7 C CA . GLY A 1 2 ? 7.000 8.000 9.000 1.00 80.00
        \\ATOM 8 C C . GLY A 1 2 ? 8.000 9.000 10.000 1.00 80.00
        \\ATOM 9 O O . GLY A 1 2 ? 9.000 10.000 11.000 1.00 80.00
        \\ATOM 10 O OXT . GLY A 1 2 ? 10.000 11.000 12.000 1.00 80.00
        \\#
        \\
    ;

    var input = try parse(std.testing.allocator, source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 10), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), input.x[0], 0.0001);
    try std.testing.expectApproxEqAbs(@as(f64, 12.0), input.z[9], 0.0001);
    try std.testing.expectEqualStrings("ALA", input.residue.?[0].slice());
    try std.testing.expectEqualStrings("CB", input.atom_name.?[3].slice());
    try std.testing.expectEqualStrings("GLY", input.residue.?[9].slice());
    try std.testing.expectEqualStrings("OXT", input.atom_name.?[9].slice());
    try std.testing.expectEqual(@as(i32, 2), input.residue_num.?[9]);
    try std.testing.expectEqual(@as(u8, 8), input.element.?[9]);
}

test "rejects AlphaFold-like mmCIF when ATOM row count does not match sequence" {
    const source =
        \\data_AF_BAD
        \\_struct_ref.pdbx_seq_one_letter_code AG
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_alt_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_ins_code
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N . ALA A 1 1 ? 1.000 2.000 3.000
        \\
    ;

    try std.testing.expectError(error.AtomCountMismatch, parse(std.testing.allocator, source));
}

const std = @import("std");
const compressed = @import("compressed.zig");
const elem = @import("element.zig");
const input_io = @import("input_io.zig");
const mmap_reader = @import("mmap_reader.zig");
const types = @import("types.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
pub const InputIoMode = input_io.InputIoMode;

pub const ParseOptions = struct {
    io_mode: InputIoMode = .mmap,
};

pub const ParseError = error{
    UnsupportedLayout,
    NoAtomRecord,
    InvalidCoordinate,
    InvalidResidueNumber,
};

pub fn parseFile(allocator: Allocator, io: std.Io, path: []const u8) !AtomInput {
    return parseFileWithOptions(allocator, io, path, .{});
}

pub fn parseFileWithOptions(allocator: Allocator, io: std.Io, path: []const u8, options: ParseOptions) !AtomInput {
    if (compressed.isCompressed(path)) {
        const data = try compressed.read(allocator, path);
        defer allocator.free(data);
        return parse(allocator, data);
    }

    switch (options.io_mode.resolve(.mmap)) {
        .mmap => {
            const mapped = try mmap_reader.mmapFile(allocator, io, path);
            defer mapped.deinit();
            return parse(allocator, mapped.data);
        },
        .read => {
            const file = try std.Io.Dir.cwd().openFile(io, path, .{});
            defer file.close(io);
            var read_buf: [64 * 1024]u8 = undefined;
            var reader = file.reader(io, &read_buf);
            const data = try reader.interface.allocRemaining(allocator, .unlimited);
            defer allocator.free(data);
            return parse(allocator, data);
        },
        .auto => unreachable,
    }
}

const AtomSiteLayout = struct {
    num_cols: usize = 0,
    data_start: usize = 0,
    group_pdb: ?usize = null,
    type_symbol: ?usize = null,
    atom_name: ?usize = null,
    residue: ?usize = null,
    chain_id: ?usize = null,
    residue_num: ?usize = null,
    insertion_code: ?usize = null,
    alt_loc: ?usize = null,
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    model_num: ?usize = null,

    fn hasRequiredFields(self: AtomSiteLayout) bool {
        return self.group_pdb != null and
            self.type_symbol != null and
            self.atom_name != null and
            self.residue != null and
            self.chain_id != null and
            self.residue_num != null and
            self.cartn_x != null and
            self.cartn_y != null and
            self.cartn_z != null;
    }
};

const AtomRow = struct {
    x: f64,
    y: f64,
    z: f64,
    radius: f64,
    residue: types.FixedString5,
    atom_name: types.FixedString4,
    element: u8,
    chain_id: types.FixedString4,
    residue_num: i32,
    insertion_code: types.FixedString4,
};

pub fn parse(allocator: Allocator, source: []const u8) !AtomInput {
    const layout = try findAtomSiteLayout(source);
    const atom_count = try countSupportedRows(source, layout);
    if (atom_count == 0) return ParseError.NoAtomRecord;

    var x = try allocator.alloc(f64, atom_count);
    errdefer allocator.free(x);
    var y = try allocator.alloc(f64, atom_count);
    errdefer allocator.free(y);
    var z = try allocator.alloc(f64, atom_count);
    errdefer allocator.free(z);
    var r = try allocator.alloc(f64, atom_count);
    errdefer allocator.free(r);
    var residue = try allocator.alloc(types.FixedString5, atom_count);
    errdefer allocator.free(residue);
    var atom_name = try allocator.alloc(types.FixedString4, atom_count);
    errdefer allocator.free(atom_name);
    var element = try allocator.alloc(u8, atom_count);
    errdefer allocator.free(element);
    var chain_id = try allocator.alloc(types.FixedString4, atom_count);
    errdefer allocator.free(chain_id);
    var residue_num = try allocator.alloc(i32, atom_count);
    errdefer allocator.free(residue_num);
    var insertion_code = try allocator.alloc(types.FixedString4, atom_count);
    errdefer allocator.free(insertion_code);

    var row_index: usize = 0;
    var line_start = layout.data_start;
    while (line_start < source.len and row_index < atom_count) {
        const raw_line_end = lineEnd(source, line_start);
        const line = std.mem.trim(u8, source[line_start..trimCarriageReturn(source, raw_line_end)], " \t");
        if (line.len != 0) {
            var fields_storage: [64][]const u8 = undefined;
            const fields = try splitAtomRow(line, layout.num_cols, &fields_storage);
            if (!std.mem.eql(u8, fields[layout.group_pdb.?], "ATOM")) break;

            const row = try parseAtomRow(fields, layout);
            x[row_index] = row.x;
            y[row_index] = row.y;
            z[row_index] = row.z;
            r[row_index] = row.radius;
            residue[row_index] = row.residue;
            atom_name[row_index] = row.atom_name;
            element[row_index] = row.element;
            chain_id[row_index] = row.chain_id;
            residue_num[row_index] = row.residue_num;
            insertion_code[row_index] = row.insertion_code;
            row_index += 1;
        }

        if (raw_line_end >= source.len) break;
        line_start = raw_line_end + 1;
    }

    if (row_index != atom_count) return ParseError.UnsupportedLayout;

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

fn findAtomSiteLayout(source: []const u8) !AtomSiteLayout {
    var line_start: usize = 0;
    var in_loop = false;
    var layout = AtomSiteLayout{};

    while (line_start < source.len) {
        const raw_line_end = lineEnd(source, line_start);
        const line = std.mem.trim(u8, source[line_start..trimCarriageReturn(source, raw_line_end)], " \t");

        if (std.mem.eql(u8, line, "loop_")) {
            in_loop = true;
            layout = .{};
        } else if (in_loop and std.mem.startsWith(u8, line, "_atom_site.")) {
            mapAtomSiteColumn(&layout, line["_atom_site.".len..], layout.num_cols);
            layout.num_cols += 1;
        } else if (in_loop and layout.num_cols > 0 and line.len != 0) {
            layout.data_start = line_start;
            if (layout.num_cols > 64 or !layout.hasRequiredFields() or layout.group_pdb.? != 0) {
                return ParseError.UnsupportedLayout;
            }
            return layout;
        } else if (in_loop and line.len != 0 and !std.mem.startsWith(u8, line, "#")) {
            in_loop = false;
        }

        if (raw_line_end >= source.len) break;
        line_start = raw_line_end + 1;
    }

    return ParseError.UnsupportedLayout;
}

fn mapAtomSiteColumn(layout: *AtomSiteLayout, field: []const u8, index: usize) void {
    if (std.ascii.eqlIgnoreCase(field, "group_PDB")) {
        layout.group_pdb = index;
    } else if (std.ascii.eqlIgnoreCase(field, "type_symbol")) {
        layout.type_symbol = index;
    } else if (std.ascii.eqlIgnoreCase(field, "label_atom_id")) {
        layout.atom_name = index;
    } else if (std.ascii.eqlIgnoreCase(field, "label_comp_id")) {
        layout.residue = index;
    } else if (std.ascii.eqlIgnoreCase(field, "label_asym_id")) {
        layout.chain_id = index;
    } else if (std.ascii.eqlIgnoreCase(field, "label_seq_id")) {
        layout.residue_num = index;
    } else if (std.ascii.eqlIgnoreCase(field, "pdbx_PDB_ins_code")) {
        layout.insertion_code = index;
    } else if (std.ascii.eqlIgnoreCase(field, "label_alt_id")) {
        layout.alt_loc = index;
    } else if (std.ascii.eqlIgnoreCase(field, "Cartn_x")) {
        layout.cartn_x = index;
    } else if (std.ascii.eqlIgnoreCase(field, "Cartn_y")) {
        layout.cartn_y = index;
    } else if (std.ascii.eqlIgnoreCase(field, "Cartn_z")) {
        layout.cartn_z = index;
    } else if (std.ascii.eqlIgnoreCase(field, "pdbx_PDB_model_num")) {
        layout.model_num = index;
    }
}

fn countSupportedRows(source: []const u8, layout: AtomSiteLayout) !usize {
    var count: usize = 0;
    var line_start = layout.data_start;
    var finished = false;

    while (line_start < source.len) {
        const raw_line_end = lineEnd(source, line_start);
        const line = std.mem.trim(u8, source[line_start..trimCarriageReturn(source, raw_line_end)], " \t");

        if (line.len != 0) {
            if (finished) {
                if (std.mem.startsWith(u8, line, "ATOM ")) return ParseError.UnsupportedLayout;
            } else {
                if (!std.mem.startsWith(u8, line, "ATOM ")) {
                    if (count == 0) return ParseError.UnsupportedLayout;
                    finished = true;
                } else {
                    count += 1;
                }
            }
        }

        if (raw_line_end >= source.len) break;
        line_start = raw_line_end + 1;
    }

    return count;
}

fn splitAtomRow(line: []const u8, expected_fields: usize, storage: *[64][]const u8) ![]const []const u8 {
    var tokenizer = std.mem.tokenizeAny(u8, line, " \t");
    var count: usize = 0;
    while (tokenizer.next()) |field| {
        if (count >= storage.len) return ParseError.UnsupportedLayout;
        if (field[0] == '\'' or field[0] == '"' or field[0] == ';') return ParseError.UnsupportedLayout;
        storage[count] = field;
        count += 1;
    }
    if (count != expected_fields) return ParseError.UnsupportedLayout;
    return storage[0..count];
}

fn parseAtomRow(fields: []const []const u8, layout: AtomSiteLayout) !AtomRow {
    const alt_loc = if (layout.alt_loc) |index| fields[index] else ".";
    if (!isNullValue(alt_loc)) return ParseError.UnsupportedLayout;

    if (layout.model_num) |index| {
        const model = fields[index];
        if (!isNullValue(model) and !std.mem.eql(u8, model, "1")) {
            return ParseError.UnsupportedLayout;
        }
    }

    const residue_value = fields[layout.residue.?];
    const atom_name_value = fields[layout.atom_name.?];
    const chain_value = fields[layout.chain_id.?];
    if (isNullValue(residue_value) or residue_value.len > 5 or
        isNullValue(atom_name_value) or atom_name_value.len > 4 or
        isNullValue(chain_value) or chain_value.len > 4)
    {
        return ParseError.UnsupportedLayout;
    }

    const element_value = fields[layout.type_symbol.?];
    if (isNullValue(element_value)) return ParseError.UnsupportedLayout;
    const element_enum = elem.fromSymbol(element_value);
    if (element_enum == .X or element_enum == .H) return ParseError.UnsupportedLayout;

    const residue_number_value = fields[layout.residue_num.?];
    const residue_number: i32 = if (isNullValue(residue_number_value))
        0
    else
        std.fmt.parseInt(i32, residue_number_value, 10) catch return ParseError.InvalidResidueNumber;

    const insertion_value = if (layout.insertion_code) |index| fields[index] else "?";
    if (!isNullValue(insertion_value) and insertion_value.len > 4) return ParseError.UnsupportedLayout;

    return .{
        .x = try parseCoordinate(fields[layout.cartn_x.?]),
        .y = try parseCoordinate(fields[layout.cartn_y.?]),
        .z = try parseCoordinate(fields[layout.cartn_z.?]),
        .radius = element_enum.vdwRadius(),
        .residue = types.FixedString5.fromSlice(residue_value),
        .atom_name = types.FixedString4.fromSlice(atom_name_value),
        .element = element_enum.atomicNumber(),
        .chain_id = types.FixedString4.fromSlice(chain_value),
        .residue_num = residue_number,
        .insertion_code = types.FixedString4.fromSlice(if (isNullValue(insertion_value)) "" else insertion_value),
    };
}

fn isNullValue(value: []const u8) bool {
    return std.mem.eql(u8, value, ".") or std.mem.eql(u8, value, "?");
}

fn lineEnd(source: []const u8, line_start: usize) usize {
    return std.mem.indexOfScalarPos(u8, source, line_start, '\n') orelse source.len;
}

fn trimCarriageReturn(source: []const u8, raw_line_end: usize) usize {
    if (raw_line_end > 0 and source[raw_line_end - 1] == '\r') return raw_line_end - 1;
    return raw_line_end;
}

fn parseCoordinate(field: []const u8) !f64 {
    if (std.mem.indexOfAny(u8, field, "eE") != null) {
        return std.fmt.parseFloat(f64, field) catch ParseError.InvalidCoordinate;
    }
    return parseFixedFloat(field) orelse ParseError.InvalidCoordinate;
}

fn parseFixedFloat(field: []const u8) ?f64 {
    if (field.len == 0) return null;
    var i: usize = 0;
    var negative = false;
    if (field[i] == '-') {
        negative = true;
        i += 1;
    } else if (field[i] == '+') {
        i += 1;
    }

    var value: f64 = 0.0;
    var has_digits = false;
    while (i < field.len and std.ascii.isDigit(field[i])) : (i += 1) {
        has_digits = true;
        value = value * 10.0 + @as(f64, @floatFromInt(field[i] - '0'));
    }

    if (i < field.len and field[i] == '.') {
        i += 1;
        var factor: f64 = 0.1;
        while (i < field.len and std.ascii.isDigit(field[i])) : (i += 1) {
            has_digits = true;
            value += @as(f64, @floatFromInt(field[i] - '0')) * factor;
            factor *= 0.1;
        }
    }

    if (!has_digits or i != field.len) return null;
    return if (negative) -value else value;
}

test "parse AlphaFold-like mmCIF metadata and coordinates from ATOM rows" {
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
        \\ATOM 1    N N   . ALA A 1 1   ? 1.000  2.000  3.000  1.00 90.00
        \\ATOM 2    C CA  . ALA A 1 1   ? 2.000  3.000  4.000  1.00 90.00
        \\ATOM 3    C C   . ALA A 1 1   ? 3.000  4.000  5.000  1.00 90.00
        \\ATOM 4    C CB  . ALA A 1 1   ? 4.000  5.000  6.000  1.00 90.00
        \\ATOM 5    O O   . ALA A 1 1   ? 5.000  6.000  7.000  1.00 90.00
        \\ATOM 6    N N   . GLY A 1 2   ? 6.000  7.000  8.000  1.00 80.00
        \\ATOM 7    C CA  . GLY A 1 2   ? 7.000  8.000  9.000  1.00 80.00
        \\ATOM 8    C C   . GLY A 1 2   ? 8.000  9.000  10.000 1.00 80.00
        \\ATOM 9    O O   . GLY A 1 2   ? 9.000  10.000 11.000 1.00 80.00
        \\ATOM 10   O OXT . GLY A 1 2   ? 10.000 11.000 12.000 1.00 80.00
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

test "rejects AF model mmCIF when ATOM rows are not consecutive" {
    const allocator = std.testing.allocator;
    const cif =
        "_struct_ref.pdbx_seq_one_letter_code A\n" ++
        "ATOM 1    N N   . ALA A 1 1   ? 1.000  2.000  3.000  1.00 90.00\n" ++
        "# non-ATOM row between atom records forces generic parser fallback\n" ++
        "ATOM 2    C CA  . ALA A 1 1   ? 2.000  3.000  4.000  1.00 90.00\n" ++
        "ATOM 3    C C   . ALA A 1 1   ? 3.000  4.000  5.000  1.00 90.00\n" ++
        "ATOM 4    C CB  . ALA A 1 1   ? 4.000  5.000  6.000  1.00 90.00\n" ++
        "ATOM 5    O O   . ALA A 1 1   ? 5.000  6.000  7.000  1.00 90.00\n" ++
        "ATOM 6    O OXT . ALA A 1 1   ? 6.000  7.000  8.000  1.00 90.00\n";

    try std.testing.expectError(ParseError.UnsupportedLayout, parse(allocator, cif));
}

test "preserves multiple chains from AF model atom rows" {
    const source =
        \\data_AF_MULTIMER
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
        \\_atom_site.pdbx_PDB_model_num
        \\ATOM 1 N N  . GLY A 1 1 ? 1.0 2.0 3.0 1
        \\ATOM 2 C CA . GLY A 1 1 ? 2.0 3.0 4.0 1
        \\ATOM 3 N N  . ALA B 2 1 ? 4.0 5.0 6.0 1
        \\ATOM 4 C CA . ALA B 2 1 ? 5.0 6.0 7.0 1
        \\#
    ;

    var input = try parse(std.testing.allocator, source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 4), input.atomCount());
    try std.testing.expectEqualStrings("A", input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("B", input.chain_id.?[2].slice());
    try std.testing.expectEqualStrings("ALA", input.residue.?[2].slice());
    try std.testing.expectEqualStrings("CA", input.atom_name.?[3].slice());
    try std.testing.expectEqual(@as(i32, 1), input.residue_num.?[3]);
}

test "rejects alternate locations as unsupported fast layout" {
    const source =
        \\data_AF_ALTLOC
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_alt_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_ins_code
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 C CA A GLY A 1 ? 1.0 2.0 3.0
        \\#
    ;

    try std.testing.expectError(ParseError.UnsupportedLayout, parse(std.testing.allocator, source));
}

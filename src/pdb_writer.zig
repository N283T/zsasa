//! Write per-atom SASA back into the PDB file it was calculated from.
//!
//! `writeAnnotated` copies the source lines, with these changes (the PDB
//! output of FreeSASA follows the same convention):
//!
//! - The ATOM/HETATM records of the atoms in the calculation get the radius
//!   in the occupancy column (55-60) and the area in the temperature factor
//!   column (61-66). The other ATOM/HETATM records are dropped, together
//!   with their ANISOU, SIGATM and SIGUIJ records, the MODEL blocks and TER
//!   records left without atoms, and the CONECT records that name a dropped
//!   atom. MASTER records are dropped because their counts no longer hold.
//! - `REMARK 999 ZSASA` lines describing the columns and the calculation are
//!   added after the HEADER record (or at the top). Earlier `REMARK 999
//!   ZSASA` lines are removed, so running zsasa on its own output gives the
//!   same file.
//!
//! The atoms are matched to lines by `AtomInput.source_row`, which the PDB
//! parser sets when `record_source_rows` is true.

const std = @import("std");
const types = @import("types.zig");
const pdb_parser = @import("pdb_parser.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;

const remark_prefix = "REMARK 999 ZSASA";

pub const Error = error{
    /// `AtomInput.source_row` is not set
    MissingSourceRows,
    /// The atom arrays, the areas or the source rows differ in length
    LengthMismatch,
    /// The source rows are not strictly increasing or do not all point at
    /// ATOM/HETATM records of the source
    SourceRowMismatch,
    /// An area or radius is NaN or infinite
    NonFiniteValue,
    /// An area or radius does not fit the six columns of its field
    ValueTooWide,
};

pub const Options = struct {
    /// Algorithm name, e.g. `shrake_rupley`
    algorithm: []const u8 = "unknown",
    classifier: []const u8 = "unknown",
    probe_radius: f64 = 1.4,
    /// Name of the resolution parameter (`n_points`, `n_slices`), or null to
    /// leave it out
    detail_item: ?[]const u8 = null,
    detail_count: u32 = 0,
};

/// Write `value` in six columns: two decimals, or one when it would not fit
fn writeField(out: *std.Io.Writer, value: f64) !void {
    if (value < 999.995) return out.print("{d:6.2}", .{value});
    if (value < 9999.95) return out.print("{d:6.1}", .{value});
    return Error.ValueTooWide;
}

/// Column offset of the standard columns of an ATOM/HETATM record: 1 when an
/// overflowing serial number pushed them one column to the right
fn recordOffset(line: []const u8, kind: pdb_parser.AtomRecordKind) usize {
    return @intFromPtr(kind.fields.ptr) - @intFromPtr(line.ptr);
}

/// Serial number field of an ATOM/HETATM record, trimmed
fn serialOf(line: []const u8, offset: usize) []const u8 {
    const end = @min(line.len, 11 + offset);
    if (end <= 6) return "";
    return std.mem.trim(u8, line[6..end], " ");
}

fn isRecord(line: []const u8, name: []const u8) bool {
    return std.mem.startsWith(u8, line, name);
}

/// Write `line` with the radius and area in place of occupancy and
/// temperature factor
fn writeAtomLine(out: *std.Io.Writer, line: []const u8, offset: usize, radius: f64, area: f64) !void {
    const has_cr = line.len > 0 and line[line.len - 1] == '\r';
    const body = if (has_cr) line[0 .. line.len - 1] else line;
    const occ_start = offset + 54;
    const b_end = offset + 66;
    if (body.len >= occ_start) {
        try out.writeAll(body[0..occ_start]);
    } else {
        try out.writeAll(body);
        try out.splatByteAll(' ', occ_start - body.len);
    }
    try writeField(out, radius);
    try writeField(out, area);
    if (body.len > b_end) try out.writeAll(body[b_end..]);
    if (has_cr) try out.writeAll("\r");
}

fn writeRemarks(out: *std.Io.Writer, total_area: f64, options: Options) !void {
    try out.writeAll(remark_prefix ++ ": OCCUPANCY = RADIUS (A), B-FACTOR = SASA (A^2)\n");
    var buf: [256]u8 = undefined;
    const params = if (options.detail_item) |item|
        try std.fmt.bufPrint(&buf, remark_prefix ++ ": ALGORITHM {s}, PROBE RADIUS {d:.2}, {s} {d}", .{ options.algorithm, options.probe_radius, item, options.detail_count })
    else
        try std.fmt.bufPrint(&buf, remark_prefix ++ ": ALGORITHM {s}, PROBE RADIUS {d:.2}", .{ options.algorithm, options.probe_radius });
    try out.print("{s}\n", .{params[0..@min(params.len, 80)]});
    const classifier = try std.fmt.bufPrint(&buf, remark_prefix ++ ": CLASSIFIER {s}", .{options.classifier[0..@min(options.classifier.len, 128)]});
    try out.print("{s}\n", .{classifier[0..@min(classifier.len, 80)]});
    try out.print(remark_prefix ++ ": TOTAL AREA {d:.2}\n", .{total_area});
}

/// Write `source` to `out` with the areas and radii of `input` in its
/// ATOM/HETATM records; see the module documentation.
pub fn writeAnnotated(
    allocator: Allocator,
    out: *std.Io.Writer,
    source: []const u8,
    input: AtomInput,
    atom_areas: []const f64,
    options: Options,
) !void {
    const rows = input.source_row orelse return Error.MissingSourceRows;
    const n = input.atomCount();
    if (rows.len != n or atom_areas.len != n or input.r.len != n) return Error.LengthMismatch;
    for (1..n) |i| {
        if (rows[i] <= rows[i - 1]) return Error.SourceRowMismatch;
    }
    var total_area: f64 = 0;
    for (atom_areas, input.r) |area, radius| {
        if (!std.math.isFinite(area) or !std.math.isFinite(radius)) return Error.NonFiniteValue;
        total_area += area;
    }

    // Serial numbers of the atoms that are written, for CONECT records
    var kept_serials = std.StringHashMapUnmanaged(void).empty;
    defer kept_serials.deinit(allocator);

    // A MODEL block is held back until ENDMDL shows whether it has atoms
    var block = std.Io.Writer.Allocating.init(allocator);
    defer block.deinit();
    var in_model = false;
    var model_has_atoms = false;

    var next_atom: usize = 0;
    var last_atom_kept = false;
    var atoms_since_ter = false;
    var remarks_written = false;

    var lines = std.mem.splitScalar(u8, source, '\n');
    var index: usize = 0;
    while (lines.next()) |line| : (index += 1) {
        const is_last = lines.peek() == null;
        if (is_last and line.len == 0) break;
        const w: *std.Io.Writer = if (in_model) &block.writer else out;

        if (!remarks_written and !isRecord(line, "HEADER")) {
            try writeRemarks(out, total_area, options);
            remarks_written = true;
        }

        var keep = true;
        if (pdb_parser.atomRecordKind(line)) |kind| {
            const offset = recordOffset(line, kind);
            if (next_atom < n and rows[next_atom] == index) {
                try writeAtomLine(w, line, offset, input.r[next_atom], atom_areas[next_atom]);
                if (!is_last) try w.writeAll("\n");
                try kept_serials.put(allocator, serialOf(line, offset), {});
                next_atom += 1;
                last_atom_kept = true;
                atoms_since_ter = true;
                model_has_atoms = true;
                continue;
            }
            keep = false;
            last_atom_kept = false;
        } else if (isRecord(line, "ANISOU") or isRecord(line, "SIGATM") or isRecord(line, "SIGUIJ")) {
            keep = last_atom_kept;
        } else if (isRecord(line, "TER")) {
            keep = atoms_since_ter;
            atoms_since_ter = false;
        } else if (isRecord(line, "MODEL")) {
            in_model = true;
            model_has_atoms = false;
            atoms_since_ter = false;
            block.clearRetainingCapacity();
            try block.writer.writeAll(line);
            if (!is_last) try block.writer.writeAll("\n");
            continue;
        } else if (isRecord(line, "ENDMDL")) {
            if (in_model and model_has_atoms) {
                try out.writeAll(block.written());
                try out.writeAll(line);
                if (!is_last) try out.writeAll("\n");
            }
            in_model = false;
            block.clearRetainingCapacity();
            continue;
        } else if (isRecord(line, "CONECT")) {
            keep = conectAtomsKept(line, kept_serials);
        } else if (isRecord(line, "MASTER") or isRecord(line, remark_prefix)) {
            keep = false;
        }

        if (keep) {
            try w.writeAll(line);
            if (!is_last) try w.writeAll("\n");
        }
    }
    // A MODEL block without ENDMDL at the end of the file
    if (in_model and model_has_atoms) try out.writeAll(block.written());
    if (!remarks_written) try writeRemarks(out, total_area, options);
    if (next_atom != n) return Error.SourceRowMismatch;
}

/// Every serial number a CONECT record names belongs to a written atom
fn conectAtomsKept(line: []const u8, kept: std.StringHashMapUnmanaged(void)) bool {
    var start: usize = 6;
    while (start < line.len) : (start += 5) {
        const field = std.mem.trim(u8, line[start..@min(start + 5, line.len)], " \r");
        if (field.len > 0 and !kept.contains(field)) return false;
    }
    return true;
}

/// `writeAnnotated` into a new string
pub fn annotate(
    allocator: Allocator,
    source: []const u8,
    input: AtomInput,
    atom_areas: []const f64,
    options: Options,
) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    try writeAnnotated(allocator, &aw.writer, source, input, atom_areas, options);
    return aw.toOwnedSlice();
}

// =============================================================================
// Tests
// =============================================================================

const test_source =
    \\HEADER    TEST                                    01-JAN-26   XXXX
    \\REMARK   2 RESOLUTION. 1.00 ANGSTROMS.
    \\ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N
    \\ANISOU    1  N   GLY A   1     1000   1000   1000      0      0      0       N
    \\ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00 20.00           C
    \\TER       3      GLY A   1
    \\ATOM      4  CA  GLY B   1      20.000   0.000   0.000  1.00 20.00           C
    \\TER       5      GLY B   1
    \\HETATM    6  O   HOH A 101      40.000   0.000   0.000  1.00 20.00           O
    \\ANISOU    6  O   HOH A 101     1000   1000   1000      0      0      0       O
    \\CONECT    1    2
    \\CONECT    2    6
    \\MASTER        0    0    0    0    0    0    0    0    4    2    0    0
    \\END
    \\
;

/// Parse `source` with source rows recorded, give every atom the area
/// `10 * (index + 1)` and radius 1.5
const Parsed = struct {
    input: AtomInput,
    areas: []f64,

    fn init(source: []const u8, configure: anytype) !Parsed {
        var parser = pdb_parser.PdbParser.init(std.testing.allocator);
        parser.record_source_rows = true;
        parser.quiet = true;
        configure(&parser);
        var input = try parser.parse(source);
        errdefer input.deinit();
        const areas = try std.testing.allocator.alloc(f64, input.atomCount());
        for (areas, 0..) |*a, i| a.* = 10.0 * @as(f64, @floatFromInt(i + 1));
        for (input.r) |*r| r.* = 1.5;
        return .{ .input = input, .areas = areas };
    }

    fn deinit(self: *Parsed) void {
        self.input.deinit();
        std.testing.allocator.free(self.areas);
    }

    fn annotate(self: Parsed, source: []const u8) ![]u8 {
        return testAnnotate(source, self.input, self.areas);
    }
};

fn testAnnotate(source: []const u8, input: AtomInput, areas: []const f64) ![]u8 {
    return annotate(std.testing.allocator, source, input, areas, .{
        .algorithm = "shrake_rupley",
        .classifier = "CCD",
        .detail_item = "n_points",
        .detail_count = 100,
    });
}

fn defaults(_: *pdb_parser.PdbParser) void {}

fn chainA(p: *pdb_parser.PdbParser) void {
    p.chain_filter = &.{"A"};
}

test "annotate writes radius and area and drops the atoms left out" {
    var parsed = try Parsed.init(test_source, chainA);
    defer parsed.deinit();
    const out = try parsed.annotate(test_source);
    defer std.testing.allocator.free(out);

    const expected =
        \\HEADER    TEST                                    01-JAN-26   XXXX
        \\REMARK 999 ZSASA: OCCUPANCY = RADIUS (A), B-FACTOR = SASA (A^2)
        \\REMARK 999 ZSASA: ALGORITHM shrake_rupley, PROBE RADIUS 1.40, n_points 100
        \\REMARK 999 ZSASA: CLASSIFIER CCD
        \\REMARK 999 ZSASA: TOTAL AREA 30.00
        \\REMARK   2 RESOLUTION. 1.00 ANGSTROMS.
        \\ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.50 10.00           N
        \\ANISOU    1  N   GLY A   1     1000   1000   1000      0      0      0       N
        \\ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.50 20.00           C
        \\TER       3      GLY A   1
        \\CONECT    1    2
        \\END
        \\
    ;
    try std.testing.expectEqualStrings(expected, out);
}

test "annotate on its own output gives the same file" {
    var parsed = try Parsed.init(test_source, struct {
        fn f(p: *pdb_parser.PdbParser) void {
            p.atom_only = false;
        }
    }.f);
    defer parsed.deinit();
    const first = try parsed.annotate(test_source);
    defer std.testing.allocator.free(first);
    try std.testing.expect(std.mem.indexOf(u8, first, "HOH A 101      40.000   0.000   0.000  1.50 40.00           O\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, first, "CONECT    2    6\n") != null);

    var reparsed = try Parsed.init(first, struct {
        fn f(p: *pdb_parser.PdbParser) void {
            p.atom_only = false;
        }
    }.f);
    defer reparsed.deinit();
    const second = try reparsed.annotate(first);
    defer std.testing.allocator.free(second);
    try std.testing.expectEqualStrings(first, second);
}

test "annotate drops models without atoms" {
    const source =
        \\MODEL        1
        \\ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C
        \\TER       2      GLY A   1
        \\ENDMDL
        \\MODEL        2
        \\ATOM      1  CA  GLY A   1       0.500   0.000   0.000  1.00 20.00           C
        \\TER       2      GLY A   1
        \\ENDMDL
        \\END
    ;
    var parsed = try Parsed.init(source, struct {
        fn f(p: *pdb_parser.PdbParser) void {
            p.model_num = 2;
        }
    }.f);
    defer parsed.deinit();
    const out = try parsed.annotate(source);
    defer std.testing.allocator.free(out);

    const expected =
        \\REMARK 999 ZSASA: OCCUPANCY = RADIUS (A), B-FACTOR = SASA (A^2)
        \\REMARK 999 ZSASA: ALGORITHM shrake_rupley, PROBE RADIUS 1.40, n_points 100
        \\REMARK 999 ZSASA: CLASSIFIER CCD
        \\REMARK 999 ZSASA: TOTAL AREA 10.00
        \\MODEL        2
        \\ATOM      1  CA  GLY A   1       0.500   0.000   0.000  1.50 10.00           C
        \\TER       2      GLY A   1
        \\ENDMDL
        \\END
    ;
    try std.testing.expectEqualStrings(expected, out);
}

test "annotate pads short records, keeps CRLF and follows shifted columns" {
    const source = "ATOM      1  CA  GLY A   1       0.000   0.000   0.000\r\n" ++
        "ATOM  100000  CA  GLY A   2       5.000   0.000   0.000  1.00 20.00           C\r\n";
    var parsed = try Parsed.init(source, defaults);
    defer parsed.deinit();
    try std.testing.expectEqual(@as(usize, 2), parsed.input.atomCount());
    parsed.areas[1] = 1234.5;
    const out = try parsed.annotate(source);
    defer std.testing.allocator.free(out);

    try std.testing.expect(std.mem.indexOf(u8, out, "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.50 10.00\r\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "ATOM  100000  CA  GLY A   2       5.000   0.000   0.000  1.501234.5           C\r\n") != null);
}

test "annotate rejects values that do not fit and atoms that do not match" {
    var parsed = try Parsed.init(test_source, defaults);
    defer parsed.deinit();
    const allocator = std.testing.allocator;

    parsed.areas[0] = 10000;
    try std.testing.expectError(Error.ValueTooWide, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    parsed.areas[0] = std.math.inf(f64);
    try std.testing.expectError(Error.NonFiniteValue, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    parsed.areas[0] = 1;

    const rows = @constCast(parsed.input.source_row.?);
    const saved = rows[1];
    rows[1] = 3; // the ANISOU record of atom 1
    try std.testing.expectError(Error.SourceRowMismatch, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    rows[1] = 0;
    try std.testing.expectError(Error.SourceRowMismatch, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    rows[1] = saved;

    try std.testing.expectError(Error.LengthMismatch, annotate(allocator, test_source, parsed.input, parsed.areas[1..], .{}));
    var no_rows = parsed.input;
    no_rows.source_row = null;
    try std.testing.expectError(Error.MissingSourceRows, annotate(allocator, test_source, no_rows, parsed.areas, .{}));
}

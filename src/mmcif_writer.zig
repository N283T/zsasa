//! Write per-atom SASA back into the mmCIF file it was calculated from.
//!
//! `writeAnnotated` copies the source text unchanged except for:
//!
//! - two columns appended to the first `_atom_site` loop,
//!   `_atom_site.zsasa_value` (area in Å²) and `_atom_site.zsasa_radius`
//!   (radius in Å). Rows that were not part of the calculation (hydrogens,
//!   HETATM, water, other models, alternate locations, unselected chains)
//!   get `?`. When the loop already has these columns (a file written by an
//!   earlier run), their values are replaced instead.
//! - `_zsasa_parameters` and `_zsasa_results` categories added at the end of
//!   the data block that holds `_atom_site`. Existing `_zsasa_*` categories
//!   in that block and the blocks before it are removed first, so running
//!   zsasa on its own output gives the same file.
//!
//! The atoms are matched to `_atom_site` rows by `AtomInput.source_row`,
//! which the mmCIF parser sets when `record_source_rows` is true.

const std = @import("std");
const cif = @import("cif_tokenizer.zig");
const types = @import("types.zig");
const analysis = @import("analysis.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;

pub const value_item = "zsasa_value";
pub const radius_item = "zsasa_radius";
const category_prefix = "_zsasa_";

pub const Error = error{
    /// `AtomInput.source_row` is not set
    MissingSourceRows,
    /// The atom arrays, the areas or the source rows differ in length
    LengthMismatch,
    /// The source rows are not strictly increasing or do not all exist in
    /// the `_atom_site` loop of the source
    SourceRowMismatch,
    /// An area or radius is NaN or infinite
    NonFiniteValue,
    /// The source has no `_atom_site` items
    NoAtomSiteLoop,
    /// The source has an unterminated string or text field, a tag without a
    /// value or a truncated `_atom_site` loop
    MalformedCif,
};

pub const Options = struct {
    /// `_zsasa_parameters.algorithm`, e.g. `shrake_rupley`
    algorithm: []const u8 = "unknown",
    /// `_zsasa_parameters.classifier`
    classifier: []const u8 = "unknown",
    probe_radius: f64 = 1.4,
    /// Item name of the resolution parameter (`n_points`, `n_slices`), or
    /// null to leave it out
    detail_item: ?[]const u8 = null,
    detail_count: u32 = 0,
    /// Polarity class of every atom, for the polar/apolar rows of
    /// `_zsasa_results` (see `analysis.calculateAtomPolarSummary`)
    atom_classes: ?[]const analysis.AtomClass = null,
};

const Spanned = cif.Tokenizer.Spanned;

/// Tokenizer with one token of lookahead
const Lexer = struct {
    tokenizer: cif.Tokenizer,
    pending: ?Spanned = null,

    fn next(self: *Lexer) Spanned {
        if (self.pending) |t| {
            self.pending = null;
            return t;
        }
        return self.tokenizer.nextSpanned();
    }

    fn peek(self: *Lexer) Spanned {
        if (self.pending == null) self.pending = self.tokenizer.nextSpanned();
        return self.pending.?;
    }
};

/// Copies the source to `out`, replacing spans as it goes
const Emitter = struct {
    source: []const u8,
    out: *std.Io.Writer,
    cursor: usize = 0,

    /// Copy the source up to `pos`
    fn copyTo(self: *Emitter, pos: usize) !void {
        try self.out.writeAll(self.source[self.cursor..pos]);
        self.cursor = pos;
    }

    /// Replace `source[start..end]` with nothing; also drop the rest of the
    /// line when it is blank, and a `#` line that follows
    fn remove(self: *Emitter, start: usize, end: usize) !void {
        try self.copyTo(start);
        self.cursor = skipBlankLine(self.source, end);
        const after_hash = blk: {
            var p = self.cursor;
            while (p < self.source.len and (self.source[p] == ' ' or self.source[p] == '\t')) p += 1;
            if (p < self.source.len and self.source[p] == '#') {
                const q = skipBlankLine(self.source, p + 1);
                // Only a line that is nothing but `#`
                if (q == self.source.len or self.source[q - 1] == '\n') break :blk q;
            }
            break :blk self.cursor;
        };
        self.cursor = after_hash;
    }
};

/// Skip spaces, tabs and one line break after `pos` when the rest of the
/// line is blank; otherwise return `pos`
fn skipBlankLine(source: []const u8, pos: usize) usize {
    var p = pos;
    while (p < source.len and (source[p] == ' ' or source[p] == '\t' or source[p] == '\r')) p += 1;
    if (p == source.len) return p;
    if (source[p] == '\n') return p + 1;
    return pos;
}

fn isZsasaCategoryTag(tag: []const u8) bool {
    return std.ascii.startsWithIgnoreCase(tag, category_prefix);
}

fn isAtomSiteTag(tag: []const u8) bool {
    return std.ascii.startsWithIgnoreCase(tag, "_atom_site.");
}

/// The token is a semicolon-delimited text field (it ends with a `;` in
/// column 1, so whatever follows it on the line is safer on a new line)
fn isTextField(source: []const u8, t: Spanned) bool {
    const value = switch (t.token) {
        .value => |v| v,
        else => return false,
    };
    return source[t.start] == ';' and @intFromPtr(value.ptr) == @intFromPtr(source.ptr) + t.start + 1;
}

/// Atom areas and radii as the text of the new columns
const AtomValues = struct {
    rows: []const u32,
    areas: []const f64,
    radii: []const f64,
    /// Next atom to match
    next: usize = 0,

    /// The atom read from `row`, if any; rows must be asked for in order
    fn atomForRow(self: *AtomValues, row: usize) ?usize {
        if (self.next < self.rows.len and self.rows[self.next] == row) {
            self.next += 1;
            return self.next - 1;
        }
        return null;
    }

    fn allMatched(self: AtomValues) bool {
        return self.next == self.rows.len;
    }
};

const NumberText = struct {
    buf: [64]u8 = undefined,
    len: usize = 0,

    fn slice(self: *const NumberText) []const u8 {
        return self.buf[0..self.len];
    }
};

fn formatNumber(value: ?f64) NumberText {
    var text = NumberText{};
    if (value) |v| {
        const s = std.fmt.bufPrint(&text.buf, "{d:.3}", .{v}) catch unreachable;
        text.len = s.len;
    } else {
        text.buf[0] = '?';
        text.len = 1;
    }
    return text;
}

/// Write `value` as a CIF value, quoted when it has to be
pub fn writeCifValue(out: *std.Io.Writer, value: []const u8) !void {
    if (!needsQuoting(value)) return out.writeAll(value);
    const has_single = std.mem.indexOfScalar(u8, value, '\'') != null;
    const has_double = std.mem.indexOfScalar(u8, value, '"') != null;
    const has_newline = std.mem.indexOfScalar(u8, value, '\n') != null;
    if (!has_newline and !has_single) return out.print("'{s}'", .{value});
    if (!has_newline and !has_double) return out.print("\"{s}\"", .{value});
    try out.print("\n;{s}\n;", .{value});
}

fn needsQuoting(value: []const u8) bool {
    if (value.len == 0) return true;
    switch (value[0]) {
        '_', '#', '$', '\'', '"', '[', ']', ';' => return true,
        else => {},
    }
    if (std.mem.eql(u8, value, ".") or std.mem.eql(u8, value, "?")) return true;
    for ([_][]const u8{ "data_", "loop_", "save_", "global_", "stop_" }) |word| {
        if (std.ascii.startsWithIgnoreCase(value, word)) return true;
    }
    for (value) |c| {
        if (c == ' ' or c == '\t' or c == '\n' or c == '\r') return true;
    }
    return false;
}

/// Write `source` to `out` with the areas and radii of `input` added to its
/// `_atom_site` items and the zsasa summary categories added to its block.
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
    for (atom_areas, input.r) |area, radius| {
        if (!std.math.isFinite(area) or !std.math.isFinite(radius)) return Error.NonFiniteValue;
    }

    var values = AtomValues{ .rows = rows, .areas = atom_areas, .radii = input.r };
    var lexer = Lexer{ .tokenizer = cif.Tokenizer.init(source) };
    var emitter = Emitter{ .source = source, .out = out };
    var atom_site_done = false;

    while (true) {
        const t = lexer.next();
        switch (t.token) {
            .eof => {
                if (!atom_site_done) return Error.NoAtomSiteLoop;
                try emitter.copyTo(t.start);
                try writeSummary(allocator, out, source, emitter.cursor, input, atom_areas, options);
                return;
            },
            .data_block => {
                if (atom_site_done) {
                    try emitter.copyTo(t.start);
                    try writeSummary(allocator, out, source, emitter.cursor, input, atom_areas, options);
                    try emitter.copyTo(source.len);
                    return;
                }
            },
            .loop => {
                const first = lexer.peek();
                const first_tag = switch (first.token) {
                    .tag => |tag| tag,
                    else => continue,
                };
                if (isZsasaCategoryTag(first_tag)) {
                    const end = try skipLoop(&lexer, t.end);
                    try emitter.remove(t.start, end);
                } else if (!atom_site_done and isAtomSiteTag(first_tag)) {
                    try annotateAtomSiteLoop(&lexer, &emitter, &values);
                    atom_site_done = true;
                }
            },
            .tag => |tag| {
                if (isZsasaCategoryTag(tag)) {
                    const value = lexer.next();
                    if (value.token != .value) return Error.MalformedCif;
                    try emitter.remove(t.start, value.end);
                } else if (!atom_site_done and isAtomSiteTag(tag)) {
                    try annotateAtomSitePairs(&lexer, &emitter, &values, t);
                    atom_site_done = true;
                }
            },
            .err => return Error.MalformedCif,
            else => {},
        }
    }
}

/// Consume the tags and values of a loop; return the end of its last token
fn skipLoop(lexer: *Lexer, loop_end: usize) !usize {
    var end = loop_end;
    while (lexer.peek().token == .tag) end = lexer.next().end;
    while (true) {
        const t = lexer.peek();
        switch (t.token) {
            .value => end = lexer.next().end,
            .err => return Error.MalformedCif,
            else => return end,
        }
    }
}

/// Column indices of the zsasa items of an `_atom_site` loop
const ZsasaColumns = struct {
    value: ?usize = null,
    radius: ?usize = null,

    fn note(self: *ZsasaColumns, tag: []const u8, index: usize) void {
        const item = tag["_atom_site.".len..];
        if (std.ascii.eqlIgnoreCase(item, value_item)) self.value = index;
        if (std.ascii.eqlIgnoreCase(item, radius_item)) self.radius = index;
    }
};

fn annotateAtomSiteLoop(lexer: *Lexer, emitter: *Emitter, values: *AtomValues) !void {
    var columns = ZsasaColumns{};
    var n_cols: usize = 0;
    var header_end: usize = 0;
    while (true) {
        const t = lexer.peek();
        const tag = switch (t.token) {
            .tag => |tag| tag,
            else => break,
        };
        if (!isAtomSiteTag(tag)) break;
        _ = lexer.next();
        columns.note(tag, n_cols);
        n_cols += 1;
        header_end = t.end;
    }

    try emitter.copyTo(header_end);
    if (columns.value == null) try emitter.out.writeAll("\n_atom_site." ++ value_item);
    if (columns.radius == null) try emitter.out.writeAll("\n_atom_site." ++ radius_item);

    var row: usize = 0;
    var col: usize = 0;
    var atom: ?usize = values.atomForRow(0);
    while (true) {
        const t = lexer.peek();
        switch (t.token) {
            .value => {},
            .err => return Error.MalformedCif,
            else => break,
        }
        _ = lexer.next();
        const area = formatNumber(if (atom) |i| values.areas[i] else null);
        const radius = formatNumber(if (atom) |i| values.radii[i] else null);
        if (columns.value == col) {
            try emitter.copyTo(t.start);
            try emitter.out.writeAll(area.slice());
            emitter.cursor = t.end;
        } else if (columns.radius == col) {
            try emitter.copyTo(t.start);
            try emitter.out.writeAll(radius.slice());
            emitter.cursor = t.end;
        }
        col += 1;
        if (col == n_cols) {
            if (columns.value == null or columns.radius == null) {
                try emitter.copyTo(t.end);
                const sep: []const u8 = if (isTextField(emitter.source, t)) "\n" else " ";
                if (columns.value == null) try emitter.out.print("{s}{s}", .{ sep, area.slice() });
                if (columns.radius == null) {
                    const sep2: []const u8 = if (columns.value == null) " " else sep;
                    try emitter.out.print("{s}{s}", .{ sep2, radius.slice() });
                }
            }
            col = 0;
            row += 1;
            atom = values.atomForRow(row);
        }
    }
    if (col != 0) return Error.MalformedCif;
    if (!values.allMatched()) return Error.SourceRowMismatch;
}

/// `_atom_site` items written as tag/value pairs: one atom. `first` is the
/// first tag, already consumed.
fn annotateAtomSitePairs(lexer: *Lexer, emitter: *Emitter, values: *AtomValues, first: Spanned) !void {
    const atom = values.atomForRow(0);
    const area = formatNumber(if (atom) |i| values.areas[i] else null);
    const radius = formatNumber(if (atom) |i| values.radii[i] else null);
    var columns = ZsasaColumns{};
    var index: usize = 0;
    var last_end: usize = first.end;
    var tag_token = first;
    while (true) {
        const tag = tag_token.token.tag;
        const value = lexer.next();
        if (value.token != .value) return Error.MalformedCif;
        columns.note(tag, index);
        if (columns.value == index or columns.radius == index) {
            try emitter.copyTo(value.start);
            try emitter.out.writeAll(if (columns.value == index) area.slice() else radius.slice());
            emitter.cursor = value.end;
        }
        last_end = value.end;
        index += 1;

        const t = lexer.peek();
        switch (t.token) {
            .tag => |next_tag| if (!isAtomSiteTag(next_tag)) break,
            else => break,
        }
        tag_token = lexer.next();
    }
    try emitter.copyTo(last_end);
    if (columns.value == null) try emitter.out.print("\n_atom_site." ++ value_item ++ " {s}", .{area.slice()});
    if (columns.radius == null) try emitter.out.print("\n_atom_site." ++ radius_item ++ " {s}", .{radius.slice()});
    if (!values.allMatched()) return Error.SourceRowMismatch;
}

/// Write the `_zsasa_parameters` and `_zsasa_results` categories. `pos` is
/// where in the source they go, to start them on a line of their own.
fn writeSummary(
    allocator: Allocator,
    out: *std.Io.Writer,
    source: []const u8,
    pos: usize,
    input: AtomInput,
    atom_areas: []const f64,
    options: Options,
) !void {
    if (pos > 0 and source[pos - 1] != '\n') try out.writeAll("\n");

    try out.writeAll("_zsasa_parameters.algorithm ");
    try writeCifValue(out, options.algorithm);
    try out.writeAll("\n_zsasa_parameters.classifier ");
    try writeCifValue(out, options.classifier);
    try out.print("\n_zsasa_parameters.probe_radius {d:.3}\n", .{options.probe_radius});
    if (options.detail_item) |item| {
        try out.print("_zsasa_parameters.{s} {d}\n", .{ item, options.detail_count });
    }
    try out.writeAll("#\nloop_\n_zsasa_results.type\n_zsasa_results.chain_id\n_zsasa_results.surface_area\n");

    var total: f64 = 0;
    for (atom_areas) |area| total += area;
    try out.print("total . {d:.3}\n", .{total});

    if (input.atom_name != null or input.element != null or options.atom_classes != null) {
        const polar = try analysis.calculateAtomPolarSummary(input, atom_areas, options.atom_classes);
        try out.print("polar . {d:.3}\napolar . {d:.3}\n", .{ polar.polar_sasa, polar.apolar_sasa });
    }

    if (input.chain_id) |chain_ids| {
        // Areas by chain, in the order the chains first appear
        var chains = std.StringArrayHashMapUnmanaged(f64).empty;
        defer chains.deinit(allocator);
        for (atom_areas, 0..) |area, i| {
            const chain = if (input.chain_id_full) |full| full[i] else chain_ids[i].slice();
            const entry = try chains.getOrPut(allocator, chain);
            if (!entry.found_existing) entry.value_ptr.* = 0;
            entry.value_ptr.* += area;
        }
        var it = chains.iterator();
        while (it.next()) |c| {
            try out.writeAll("chain ");
            try writeCifValue(out, c.key_ptr.*);
            try out.print(" {d:.3}\n", .{c.value_ptr.*});
        }
    }
    try out.writeAll("#\n");
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

const mmcif_parser = @import("mmcif_parser.zig");

const test_source =
    \\data_TEST
    \\#
    \\_entry.id TEST
    \\#
    \\loop_
    \\_atom_site.group_PDB
    \\_atom_site.id
    \\_atom_site.type_symbol
    \\_atom_site.label_atom_id
    \\_atom_site.label_comp_id
    \\_atom_site.label_asym_id
    \\_atom_site.Cartn_x
    \\_atom_site.Cartn_y
    \\_atom_site.Cartn_z
    \\_atom_site.pdbx_PDB_model_num
    \\ATOM   1 N N   ALA A 10.0 20.0 30.0 1
    \\ATOM   2 C CA  ALA A 11.0 21.0 31.0 1
    \\HETATM 3 O O   HOH B 12.0 22.0 32.0 1
    \\ATOM   4 C CB  ALA C 13.0 23.0 33.0 1
    \\#
    \\_exptl.method 'X-RAY DIFFRACTION'
    \\#
    \\
;

/// Parse `source` with source rows recorded and give every atom the area
/// `10 * (index + 1)`
const Parsed = struct {
    input: AtomInput,
    areas: []f64,

    fn init(source: []const u8, configure: anytype) !Parsed {
        var parser = mmcif_parser.MmcifParser.init(std.testing.allocator);
        parser.record_source_rows = true;
        parser.parse_inline_ccd = false;
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
        return mmcif_writer_annotate(source, self.input, self.areas);
    }
};

fn mmcif_writer_annotate(source: []const u8, input: AtomInput, areas: []const f64) ![]u8 {
    return annotate(std.testing.allocator, source, input, areas, .{
        .algorithm = "shrake_rupley",
        .classifier = "CCD",
        .detail_item = "n_points",
        .detail_count = 100,
    });
}

fn defaults(_: *mmcif_parser.MmcifParser) void {}

test "annotate adds columns, leaves unused rows unknown and adds the summary" {
    var parsed = try Parsed.init(test_source, defaults);
    defer parsed.deinit();
    const out = try parsed.annotate(test_source);
    defer std.testing.allocator.free(out);

    const expected =
        \\data_TEST
        \\#
        \\_entry.id TEST
        \\#
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.pdbx_PDB_model_num
        \\_atom_site.zsasa_value
        \\_atom_site.zsasa_radius
        \\ATOM   1 N N   ALA A 10.0 20.0 30.0 1 10.000 1.500
        \\ATOM   2 C CA  ALA A 11.0 21.0 31.0 1 20.000 1.500
        \\HETATM 3 O O   HOH B 12.0 22.0 32.0 1 ? ?
        \\ATOM   4 C CB  ALA C 13.0 23.0 33.0 1 30.000 1.500
        \\#
        \\_exptl.method 'X-RAY DIFFRACTION'
        \\#
        \\_zsasa_parameters.algorithm shrake_rupley
        \\_zsasa_parameters.classifier CCD
        \\_zsasa_parameters.probe_radius 1.400
        \\_zsasa_parameters.n_points 100
        \\#
        \\loop_
        \\_zsasa_results.type
        \\_zsasa_results.chain_id
        \\_zsasa_results.surface_area
        \\total . 60.000
        \\polar . 10.000
        \\apolar . 50.000
        \\chain A 30.000
        \\chain C 30.000
        \\#
        \\
    ;
    try std.testing.expectEqualStrings(expected, out);
}

test "annotate on its own output gives the same file" {
    var parsed = try Parsed.init(test_source, defaults);
    defer parsed.deinit();
    const first = try parsed.annotate(test_source);
    defer std.testing.allocator.free(first);

    // The annotated file reads back to the same atoms and rows
    var reparsed = try Parsed.init(first, defaults);
    defer reparsed.deinit();
    try std.testing.expectEqualSlices(u32, parsed.input.source_row.?, reparsed.input.source_row.?);

    const second = try reparsed.annotate(first);
    defer std.testing.allocator.free(second);
    try std.testing.expectEqualStrings(first, second);
}

test "annotate replaces existing values" {
    var parsed = try Parsed.init(test_source, defaults);
    defer parsed.deinit();
    const first = try parsed.annotate(test_source);
    defer std.testing.allocator.free(first);

    // Different areas, and HETATM now part of the calculation
    var reparsed = try Parsed.init(first, struct {
        fn f(p: *mmcif_parser.MmcifParser) void {
            p.atom_only = false;
        }
    }.f);
    defer reparsed.deinit();
    const second = try reparsed.annotate(first);
    defer std.testing.allocator.free(second);

    try std.testing.expect(std.mem.indexOf(u8, second, "HETATM 3 O O   HOH B 12.0 22.0 32.0 1 30.000 1.500\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, second, "ATOM   4 C CB  ALA C 13.0 23.0 33.0 1 40.000 1.500\n") != null);
    try std.testing.expectEqual(@as(usize, 1), std.mem.count(u8, second, "_atom_site.zsasa_value"));
    try std.testing.expectEqual(@as(usize, 1), std.mem.count(u8, second, "_zsasa_parameters.algorithm"));
    try std.testing.expectEqual(@as(usize, 1), std.mem.count(u8, second, "_zsasa_results.type"));
    try std.testing.expect(std.mem.indexOf(u8, second, "total . 100.000\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, second, "chain B 30.000\n") != null);
}

test "annotate puts the summary in the atom_site block and leaves later blocks alone" {
    const source = test_source ++
        \\data_SECOND
        \\_zsasa_parameters.algorithm other
        \\
    ;
    var parsed = try Parsed.init(source, defaults);
    defer parsed.deinit();
    const out = try parsed.annotate(source);
    defer std.testing.allocator.free(out);

    const second_block = std.mem.indexOf(u8, out, "data_SECOND").?;
    const summary = std.mem.indexOf(u8, out, "_zsasa_parameters.algorithm shrake_rupley").?;
    try std.testing.expect(summary < second_block);
    try std.testing.expect(std.mem.endsWith(u8, out, "data_SECOND\n_zsasa_parameters.algorithm other\n"));
}

test "annotate handles atom_site written as tag/value pairs" {
    const source =
        \\data_ONE
        \\_atom_site.group_PDB ATOM
        \\_atom_site.type_symbol C
        \\_atom_site.label_atom_id CA
        \\_atom_site.label_comp_id GLY
        \\_atom_site.label_asym_id A
        \\_atom_site.Cartn_x 1.0
        \\_atom_site.Cartn_y 2.0
        \\_atom_site.Cartn_z 3.0
        \\#
        \\
    ;
    var parsed = try Parsed.init(source, defaults);
    defer parsed.deinit();
    const out = try parsed.annotate(source);
    defer std.testing.allocator.free(out);

    try std.testing.expect(std.mem.indexOf(u8, out, "_atom_site.Cartn_z 3.0\n_atom_site.zsasa_value 10.000\n_atom_site.zsasa_radius 1.500\n#\n") != null);

    var reparsed = try Parsed.init(out, defaults);
    defer reparsed.deinit();
    const again = try reparsed.annotate(out);
    defer std.testing.allocator.free(again);
    try std.testing.expectEqualStrings(out, again);
}

test "annotate starts new values on a new line after a text field" {
    const source =
        \\data_T
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.note
        \\ATOM C CA GLY A 1.0 2.0 3.0
        \\;a note
        \\;
        \\ATOM C CA GLY A 4.0 5.0 6.0 'quoted note'
        \\
    ;
    var parsed = try Parsed.init(source, defaults);
    defer parsed.deinit();
    const out = try parsed.annotate(source);
    defer std.testing.allocator.free(out);

    try std.testing.expect(std.mem.indexOf(u8, out, ";a note\n;\n10.000 1.500\n") != null);
    try std.testing.expect(std.mem.indexOf(u8, out, "'quoted note' 20.000 1.500\n") != null);

    var reparsed = try Parsed.init(out, defaults);
    defer reparsed.deinit();
    try std.testing.expectEqual(@as(usize, 2), reparsed.input.atomCount());
}

test "annotate keeps the rows of atoms that other options left out" {
    const source =
        \\data_M
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_alt_id
        \\_atom_site.occupancy
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.pdbx_PDB_model_num
        \\ATOM C CA GLY A . 1.0 1.0 2.0 3.0 1
        \\ATOM C CB SER A A 0.6 4.0 5.0 6.0 1
        \\ATOM C CB SER A B 0.4 4.5 5.0 6.0 1
        \\ATOM H H  GLY A . 1.0 1.5 2.0 3.0 1
        \\ATOM C CA GLY B . 1.0 9.0 9.0 9.0 1
        \\ATOM C CA GLY A . 1.0 1.0 2.0 3.0 2
        \\
    ;
    var parsed = try Parsed.init(source, struct {
        fn f(p: *mmcif_parser.MmcifParser) void {
            p.first_model_only = true;
            p.chain_filter = &.{"A"};
        }
    }.f);
    defer parsed.deinit();
    try std.testing.expectEqualSlices(u32, &.{ 0, 1 }, parsed.input.source_row.?);

    const out = try parsed.annotate(source);
    defer std.testing.allocator.free(out);
    const rows =
        \\ATOM C CA GLY A . 1.0 1.0 2.0 3.0 1 10.000 1.500
        \\ATOM C CB SER A A 0.6 4.0 5.0 6.0 1 20.000 1.500
        \\ATOM C CB SER A B 0.4 4.5 5.0 6.0 1 ? ?
        \\ATOM H H  GLY A . 1.0 1.5 2.0 3.0 1 ? ?
        \\ATOM C CA GLY B . 1.0 9.0 9.0 9.0 1 ? ?
        \\ATOM C CA GLY A . 1.0 1.0 2.0 3.0 2 ? ?
        \\
    ;
    try std.testing.expect(std.mem.indexOf(u8, out, rows) != null);
}

test "source rows follow removed residues" {
    var parsed = try Parsed.init(test_source, struct {
        fn f(p: *mmcif_parser.MmcifParser) void {
            p.atom_only = false;
        }
    }.f);
    defer parsed.deinit();
    try std.testing.expectEqualSlices(u32, &.{ 0, 1, 2, 3 }, parsed.input.source_row.?);
    _ = try parsed.input.removeResidues(.{ .water = true });
    try std.testing.expectEqualSlices(u32, &.{ 0, 1, 3 }, parsed.input.source_row.?);
}

test "annotate rejects atoms that do not match the source" {
    var parsed = try Parsed.init(test_source, defaults);
    defer parsed.deinit();
    const allocator = std.testing.allocator;

    // A row past the end of the loop
    const rows = @constCast(parsed.input.source_row.?);
    const saved = rows[2];
    rows[2] = 9;
    try std.testing.expectError(Error.SourceRowMismatch, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    // Rows out of order
    rows[2] = 0;
    try std.testing.expectError(Error.SourceRowMismatch, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    rows[2] = saved;

    parsed.areas[0] = std.math.nan(f64);
    try std.testing.expectError(Error.NonFiniteValue, annotate(allocator, test_source, parsed.input, parsed.areas, .{}));
    parsed.areas[0] = 1;

    try std.testing.expectError(Error.LengthMismatch, annotate(allocator, test_source, parsed.input, parsed.areas[1..], .{}));
    try std.testing.expectError(Error.NoAtomSiteLoop, annotate(allocator, "data_X\n_entry.id X\n", parsed.input, parsed.areas, .{}));

    var no_rows = parsed.input;
    no_rows.source_row = null;
    try std.testing.expectError(Error.MissingSourceRows, annotate(allocator, test_source, no_rows, parsed.areas, .{}));
}

test "writeCifValue quotes values that need it" {
    const cases = [_]struct { []const u8, []const u8 }{
        .{ "CCD", "CCD" },
        .{ "", "''" },
        .{ "?", "'?'" },
        .{ ".", "'.'" },
        .{ "_x", "'_x'" },
        .{ "data_x", "'data_x'" },
        .{ "my classifier.toml", "'my classifier.toml'" },
        .{ "it's here", "\"it's here\"" },
        .{ "a'b \"c", "\n;a'b \"c\n;" },
        .{ "two\nlines", "\n;two\nlines\n;" },
    };
    for (cases) |case| {
        var aw = std.Io.Writer.Allocating.init(std.testing.allocator);
        defer aw.deinit();
        try writeCifValue(&aw.writer, case[0]);
        try std.testing.expectEqualStrings(case[1], aw.written());
    }
}

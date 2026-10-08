const std = @import("std");
const analysis = @import("analysis.zig");
const types = @import("types.zig");

const Allocator = std.mem.Allocator;
const SasaResult = types.SasaResult;
const AtomInput = types.AtomInput;

pub const TextOutputOptions = struct {
    input_name: []const u8 = "input",
    classifier_name: []const u8 = "unknown",
    algorithm_name: []const u8 = "unknown",
    probe_radius: f64 = 1.4,
    detail_count: u32 = 0,
    detail_label: []const u8 = "Detail",
};

const AreaBreakdown = struct {
    total: f64 = 0,
    side_chain: f64 = 0,
    main_chain: f64 = 0,
    apolar: f64 = 0,
    polar: f64 = 0,
};

/// One residue of the RSA table. `first_atom` is the first atom of the residue
/// in the input, from which the labels of the residue are read.
const ResidueArea = struct {
    first_atom: usize,
    area: AreaBreakdown,
};

const ChainArea = struct {
    /// Full chain ID, borrowed from the input
    chain: []const u8,
    area: AreaBreakdown,
};

/// Output format options
pub const OutputFormat = enum {
    json, // Pretty-printed JSON (default)
    compact, // Single-line JSON
    csv, // CSV format
    jsonl, // JSON Lines (one JSON object per line, for batch)
    freesasa, // FreeSASA-compatible human-readable text (single calc only)
    rsa, // FreeSASA/NACCESS-compatible residue RSA text (single calc only)
};

/// JSON structure for output
const JsonOutput = struct {
    total_area: f64,
    atom_areas: []const f64,
};

/// Convert SasaResult to compact JSON string (single line)
/// Caller must free the returned slice
pub fn sasaResultToJson(allocator: Allocator, result: SasaResult) ![]u8 {
    const output = JsonOutput{
        .total_area = result.total_area,
        .atom_areas = result.atom_areas,
    };

    return std.json.Stringify.valueAlloc(allocator, output, .{});
}

/// Convert SasaResult to pretty-printed JSON string
/// Caller must free the returned slice
pub fn sasaResultToJsonPretty(allocator: Allocator, result: SasaResult) ![]u8 {
    const output = JsonOutput{
        .total_area = result.total_area,
        .atom_areas = result.atom_areas,
    };

    return std.json.Stringify.valueAlloc(allocator, output, .{
        .whitespace = .indent_2,
    });
}

/// Convert SasaResult to CSV string (basic format)
/// Format: atom_index,area (with total at end)
/// Caller must free the returned slice
pub fn sasaResultToCsv(allocator: Allocator, result: SasaResult) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    // Header
    try writer.writeAll("atom_index,area\n");

    // Atom areas
    for (result.atom_areas, 0..) |area, i| {
        try writer.print("{d},{d:.6}\n", .{ i, area });
    }

    // Total
    try writer.print("total,{d:.6}\n", .{result.total_area});

    return aw.toOwnedSlice();
}

pub fn sasaResultToFreesasa(allocator: Allocator, result: SasaResult, options: TextOutputOptions) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    try writer.writeAll("## zsasa FreeSASA-compatible output ##\n\n");
    try writer.writeAll("PARAMETERS\n");
    try writer.print("algorithm    : {s}\n", .{options.algorithm_name});
    try writer.print("classifier   : {s}\n", .{options.classifier_name});
    try writer.print("probe-radius : {d:.2}\n", .{options.probe_radius});
    if (options.detail_count > 0) {
        try writer.print("{s:<13}: {d}\n", .{ options.detail_label, options.detail_count });
    }
    try writer.print("input        : {s}\n\n", .{options.input_name});
    try writer.writeAll("RESULTS (A^2)\n");
    try writer.print("Total   : {d:10.2}\n", .{result.total_area});

    return aw.toOwnedSlice();
}

fn isMainChainAtom(atom_name: []const u8) bool {
    return std.mem.eql(u8, atom_name, "N") or
        std.mem.eql(u8, atom_name, "CA") or
        std.mem.eql(u8, atom_name, "C") or
        std.mem.eql(u8, atom_name, "O") or
        std.mem.eql(u8, atom_name, "OXT");
}

fn isPolarAtom(atom_name: []const u8, element: ?u8) bool {
    if (element) |atomic_number| {
        return atomic_number == 7 or atomic_number == 8 or atomic_number == 15 or atomic_number == 16;
    }
    const trimmed = std.mem.trim(u8, atom_name, " ");
    if (trimmed.len == 0) return false;
    const c = std.ascii.toUpper(trimmed[0]);
    return c == 'N' or c == 'O' or c == 'P' or c == 'S';
}

fn addAtomArea(area: *AreaBreakdown, atom_area: f64, atom_name: ?[]const u8, element: ?u8) void {
    area.total += atom_area;
    if (atom_name) |name| {
        if (isMainChainAtom(name)) {
            area.main_chain += atom_area;
        } else {
            area.side_chain += atom_area;
        }
        if (isPolarAtom(name, element)) {
            area.polar += atom_area;
        } else {
            area.apolar += atom_area;
        }
    } else {
        area.side_chain += atom_area;
        area.apolar += atom_area;
    }
}

/// Sum atom areas per residue. The residues are those of
/// `analysis.ResidueIdentity`, the same ones as in the `--per-residue` table
/// and in the JSONL residue map.
fn collectResidueAreas(allocator: Allocator, input: AtomInput, atom_areas: []const f64) ![]ResidueArea {
    if (!input.hasResidueInfo()) return error.MissingResidueInfo;
    if (atom_areas.len != input.atomCount()) return error.LengthMismatch;

    const identity = try analysis.ResidueIdentity.init(input);
    const residues = try allocator.alloc(ResidueArea, identity.residueCount());
    errdefer allocator.free(residues);

    const atom_names = input.atom_name;
    const elements = input.element;

    var it = identity.residues();
    var residue_idx: usize = 0;
    while (it.next()) |range| : (residue_idx += 1) {
        var area = AreaBreakdown{};
        for (range.start..range.end) |i| {
            addAtomArea(
                &area,
                atom_areas[i],
                if (atom_names) |names| names[i].slice() else null,
                if (elements) |elem| elem[i] else null,
            );
        }
        residues[residue_idx] = .{ .first_atom = range.start, .area = area };
    }

    return residues;
}

/// Sum residue areas per chain. Chains are told apart by their full ID, and
/// a chain whose residues are not contiguous in the input still gets one entry.
fn collectChainAreas(allocator: Allocator, identity: analysis.ResidueIdentity, residues: []const ResidueArea) ![]ChainArea {
    var chains = std.ArrayListUnmanaged(ChainArea).empty;
    errdefer chains.deinit(allocator);

    for (residues) |residue| {
        const residue_chain = identity.chainLabel(residue.first_atom);
        var found_idx: ?usize = null;
        for (chains.items, 0..) |*chain, i| {
            if (std.mem.eql(u8, chain.chain, residue_chain)) {
                found_idx = i;
                break;
            }
        }

        const idx = found_idx orelse blk: {
            try chains.append(allocator, .{ .chain = residue_chain, .area = .{} });
            break :blk chains.items.len - 1;
        };
        chains.items[idx].area.total += residue.area.total;
        chains.items[idx].area.side_chain += residue.area.side_chain;
        chains.items[idx].area.main_chain += residue.area.main_chain;
        chains.items[idx].area.apolar += residue.area.apolar;
        chains.items[idx].area.polar += residue.area.polar;
    }

    return chains.toOwnedSlice(allocator);
}

// RSA text (`--format=rsa`)
//
// The fixed columns of a NACCESS `.rsa` file, as NACCESS and FreeSASA write
// them (FreeSASA prints `RES %s %c%s ` and five `%7.2f%6.1f` pairs, where the
// second `%s` is the residue number field of a PDB line, columns 23-27) and
// as fixed-column readers such as Biopython's `Bio.PDB.NACCESS` slice them.
// Columns are 1-based:
//
//   RES rows
//      1-3   `RES`
//      5-7   residue name, right-justified
//      9     chain ID
//     10-13  residue number, right-justified
//     14     insertion code
//     16-80  five pairs of an absolute (F7.2) and a relative (F6.1) value:
//            all atoms, side chain, main chain, non-polar, polar
//   CHAIN rows
//      1-5   `CHAIN`
//      6-8   number of the chain, right-justified
//     10     chain ID
//     12-21, 25-34, 38-47, 51-60, 64-73
//            absolute sums (F10.1) in the order of the RES rows
//   TOTAL row
//      1-5   `TOTAL`, then the sums in the columns of the CHAIN rows
//
// A row follows these columns whenever its labels and values fit. zsasa asks
// more of a value than its field width: the first column of the field must
// stay blank, so an absolute value is at most 999.99 and a relative value
// between -99.9 and 999.9. Biopython reads only the other columns of a field
// (`line[16:22]`, `line[23:28]`, ...), and the blank keeps neighboring values
// apart for readers that split a row at blanks. A label or value that does
// not fit is written in full all the same, with a blank between it and its
// neighbors, which moves the rest of the row to the right;
// `rsaResultNeedsLegacyWidthWarning` tells when that happens.

/// Widest absolute value of a RES row that keeps its fixed columns (`999.99`).
const rsa_abs_width = 6;
/// Widest relative value of a RES row that keeps its fixed columns (`999.9`).
const rsa_rel_width = 5;
/// Width of a sum in the CHAIN and TOTAL rows.
const rsa_sum_width = 10;

/// Labels of an RSA residue row.
const RsaResidueLabel = struct {
    residue_name: []const u8,
    chain: []const u8,
    /// Residue number in decimal, without the insertion code
    number: []const u8,
    insertion_code: []const u8,

    /// `number_buf` holds the residue number and must outlive the result.
    fn init(identity: analysis.ResidueIdentity, atom: usize, number_buf: *[16]u8) RsaResidueLabel {
        return .{
            .residue_name = identity.residue_names[atom].slice(),
            .chain = identity.chainLabel(atom),
            // An i32 has at most 11 characters
            .number = std.fmt.bufPrint(number_buf, "{d}", .{identity.residue_nums[atom]}) catch unreachable,
            .insertion_code = identity.insertion_codes[atom].slice(),
        };
    }

    fn fitsFixedColumns(self: RsaResidueLabel) bool {
        return self.residue_name.len <= 3 and
            self.chain.len <= 1 and
            self.number.len <= 4 and
            self.insertion_code.len <= 1;
    }

    fn write(self: RsaResidueLabel, writer: *std.Io.Writer) !void {
        if (self.fitsFixedColumns()) {
            // Chain ID, residue number and insertion code follow each other
            // without a blank (`A1000B`), as in a PDB line
            try writer.print("RES {s:>3} {s:1}{s:>4}{s:1} ", .{ self.residue_name, self.chain, self.number, self.insertion_code });
        } else {
            try writer.print("RES {s:>3} {s} {s:>4}{s:1} ", .{ self.residue_name, self.chain, self.number, self.insertion_code });
        }
    }
};

/// Relative all-atom SASA in percent of the maximum SASA of the residue type
/// (Tien et al. 2013), or null for a residue without a reference value.
fn rsaRelativeTotal(residue_name: []const u8, total: f64) ?f64 {
    const max_sasa = analysis.MaxSASA.get(residue_name) orelse return null;
    return if (max_sasa > 0) total * 100.0 / max_sasa else null;
}

/// Write an absolute and a relative value of a RES row. Each is right-justified
/// in its field (F7.2 and F6.1) behind at least one blank; `N/A` stands for a
/// relative value without a reference value, as in FreeSASA.
fn writeAbsRel(writer: *std.Io.Writer, abs: f64, rel: ?f64) !void {
    try writer.print(" {d:>" ++ std.fmt.comptimePrint("{d}", .{rsa_abs_width}) ++ ".2}", .{abs});
    if (rel) |value| {
        try writer.print(" {d:>" ++ std.fmt.comptimePrint("{d}", .{rsa_rel_width}) ++ ".1}", .{value});
    } else {
        try writer.writeAll("   N/A");
    }
}

fn formattedExceedsWidth(comptime fmt: []const u8, args: anytype, width: usize) bool {
    return std.fmt.count(fmt, args) > width;
}

fn absRelNeedsLegacyWidthWarning(abs: f64, rel: ?f64) bool {
    if (formattedExceedsWidth("{d:.2}", .{abs}, rsa_abs_width)) return true;
    if (rel) |value| {
        if (formattedExceedsWidth("{d:.1}", .{value}, rsa_rel_width)) return true;
    }
    return false;
}

fn areaBreakdownNeedsResidueWidthWarning(area: AreaBreakdown, rel_total: ?f64) bool {
    return absRelNeedsLegacyWidthWarning(area.total, rel_total) or
        absRelNeedsLegacyWidthWarning(area.side_chain, null) or
        absRelNeedsLegacyWidthWarning(area.main_chain, null) or
        absRelNeedsLegacyWidthWarning(area.apolar, null) or
        absRelNeedsLegacyWidthWarning(area.polar, null);
}

fn areaBreakdownNeedsSummaryWidthWarning(area: AreaBreakdown) bool {
    return formattedExceedsWidth("{d:.1}", .{area.total}, rsa_sum_width) or
        formattedExceedsWidth("{d:.1}", .{area.side_chain}, rsa_sum_width) or
        formattedExceedsWidth("{d:.1}", .{area.main_chain}, rsa_sum_width) or
        formattedExceedsWidth("{d:.1}", .{area.apolar}, rsa_sum_width) or
        formattedExceedsWidth("{d:.1}", .{area.polar}, rsa_sum_width);
}

/// Whether some row of the RSA text for `result` leaves the fixed columns
/// described above: a residue name of more than three characters, a chain ID
/// of more than one, a residue number of more than four, an insertion code of
/// more than one, more than 999 chains, or a value too wide for its field.
fn rsaResultNeedsLegacyWidthWarning(allocator: Allocator, result: SasaResult, input: AtomInput) !bool {
    const residues = try collectResidueAreas(allocator, input, result.atom_areas);
    defer allocator.free(residues);
    const identity = try analysis.ResidueIdentity.init(input);
    const chains = try collectChainAreas(allocator, identity, residues);
    defer allocator.free(chains);

    var total = AreaBreakdown{};
    for (residues) |residue| {
        var number_buf: [16]u8 = undefined;
        const label = RsaResidueLabel.init(identity, residue.first_atom, &number_buf);
        const rel_total = rsaRelativeTotal(label.residue_name, residue.area.total);

        if (!label.fitsFixedColumns()) return true;
        if (areaBreakdownNeedsResidueWidthWarning(residue.area, rel_total)) return true;

        total.total += residue.area.total;
        total.side_chain += residue.area.side_chain;
        total.main_chain += residue.area.main_chain;
        total.apolar += residue.area.apolar;
        total.polar += residue.area.polar;
    }

    for (chains, 0..) |chain, i| {
        if (i + 1 > 999 or chain.chain.len > 1) return true;
        if (areaBreakdownNeedsSummaryWidthWarning(chain.area)) return true;
    }

    return areaBreakdownNeedsSummaryWidthWarning(total);
}

pub fn sasaResultToRsa(allocator: Allocator, result: SasaResult, input: AtomInput, options: TextOutputOptions) ![]u8 {
    const residues = try collectResidueAreas(allocator, input, result.atom_areas);
    defer allocator.free(residues);
    const identity = try analysis.ResidueIdentity.init(input);
    const chains = try collectChainAreas(allocator, identity, residues);
    defer allocator.free(chains);

    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    try writer.writeAll("REM  zsasa FreeSASA/NACCESS-compatible RSA\n");
    try writer.print("REM  Absolute and relative SASAs for {s}\n", .{options.input_name});
    try writer.print("REM  Atomic radii: {s}\n", .{options.classifier_name});
    // The reference values do not depend on the classifier (analysis.MaxSASA)
    try writer.writeAll("REM  Reference values for relative SASA: Tien et al. 2013\n");
    try writer.print("REM  Algorithm: {s}\n", .{options.algorithm_name});
    try writer.print("REM  Probe-radius: {d:.2}\n", .{options.probe_radius});
    if (options.detail_count > 0) {
        try writer.print("REM  {s}: {d}\n", .{ options.detail_label, options.detail_count });
    }
    try writer.writeAll("REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar\n");
    try writer.writeAll("REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL\n");

    var total = AreaBreakdown{};
    for (residues) |residue| {
        var number_buf: [16]u8 = undefined;
        const label = RsaResidueLabel.init(identity, residue.first_atom, &number_buf);

        try label.write(writer);
        try writeAbsRel(writer, residue.area.total, rsaRelativeTotal(label.residue_name, residue.area.total));
        try writeAbsRel(writer, residue.area.side_chain, null);
        try writeAbsRel(writer, residue.area.main_chain, null);
        try writeAbsRel(writer, residue.area.apolar, null);
        try writeAbsRel(writer, residue.area.polar, null);
        try writer.writeAll("\n");

        total.total += residue.area.total;
        total.side_chain += residue.area.side_chain;
        total.main_chain += residue.area.main_chain;
        total.apolar += residue.area.apolar;
        total.polar += residue.area.polar;
    }

    try writer.writeAll("END  Absolute sums over single chains surface\n");
    for (chains, 0..) |chain, i| {
        try writer.print("CHAIN{d:3} {s:1} {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}\n", .{
            i + 1,
            chain.chain,
            chain.area.total,
            chain.area.side_chain,
            chain.area.main_chain,
            chain.area.apolar,
            chain.area.polar,
        });
    }

    try writer.writeAll("END  Absolute sums over all chains\n");
    try writer.print("TOTAL      {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}\n", .{
        total.total,
        total.side_chain,
        total.main_chain,
        total.apolar,
        total.polar,
    });

    return aw.toOwnedSlice();
}

/// Write one text field of a CSV row as RFC 4180 specifies: a field that
/// contains a comma, a double quote, CR or LF is enclosed in double quotes
/// with every double quote in it doubled, and any other field is written as
/// it is. A PDB chain ID can be `,` or `"`, and the residue name of an SDF
/// molecule is taken from its title.
fn writeCsvField(writer: *std.Io.Writer, field: []const u8) !void {
    if (std.mem.findAny(u8, field, ",\"\r\n") == null) return writer.writeAll(field);

    try writer.writeByte('"');
    for (field) |c| {
        if (c == '"') try writer.writeByte('"');
        try writer.writeByte(c);
    }
    try writer.writeByte('"');
}

/// Header of the rich CSV, the CSV written for input with residue information.
pub const rich_csv_header = "chain,residue,resnum,insertion_code,atom_name,x,y,z,radius,area";

/// Convert SasaResult to rich CSV string with structural information
/// Format: `rich_csv_header`, one row per atom, and a last row that holds
/// only the total area. `insertion_code` is empty for a residue without one.
/// A column that the input does not have at all is written as `-`.
/// Text fields are quoted where RFC 4180 requires it (`writeCsvField`).
/// Caller must free the returned slice
pub fn sasaResultToRichCsv(allocator: Allocator, input: AtomInput, atom_areas: []const f64) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    // Header
    try writer.writeAll(rich_csv_header ++ "\n");

    // Atom rows
    const n = input.atomCount();
    for (0..n) |i| {
        // Chain: the full ID where the parser kept one (mmCIF chain IDs can
        // be longer than the four characters of `chain_id`)
        if (input.chain_id_full) |chains| {
            try writeCsvField(writer, chains[i]);
        } else if (input.chain_id) |chains| {
            try writeCsvField(writer, chains[i].slice());
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Residue name
        if (input.residue) |residues| {
            try writeCsvField(writer, residues[i].slice());
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Residue number
        if (input.residue_num) |nums| {
            try writer.print("{d}", .{nums[i]});
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Insertion code (empty for a residue without one)
        if (input.insertion_code) |codes| {
            try writeCsvField(writer, codes[i].slice());
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Atom name
        if (input.atom_name) |names| {
            try writeCsvField(writer, names[i].slice());
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Coordinates and radius
        try writer.print("{d:.3},{d:.3},{d:.3},{d:.3},{d:.6}\n", .{
            input.x[i],
            input.y[i],
            input.z[i],
            input.r[i],
            atom_areas[i],
        });
    }

    // Total row: every column but the area is empty
    var total: f64 = 0;
    for (atom_areas) |a| total += a;
    try writer.print(",,,,,,,,,{d:.6}\n", .{total});

    return aw.toOwnedSlice();
}

/// Write SasaResult to file with specified format.
/// Builds the output string in memory and writes to file in a single syscall
/// to avoid per-atom write overhead (thousands of syscalls per file).
pub fn writeSasaResultWithFormat(
    allocator: Allocator,
    io: std.Io,
    result: SasaResult,
    path: []const u8,
    format: OutputFormat,
) !void {
    const output_str = switch (format) {
        .json => try sasaResultToJsonPretty(allocator, result),
        .compact => try sasaResultToJson(allocator, result),
        .csv => try sasaResultToCsv(allocator, result),
        .freesasa => try sasaResultToFreesasa(allocator, result, .{}),
        .rsa => return error.MissingResidueInfo,
        .jsonl => unreachable, // JSONL is handled at batch level, not per-file
    };
    defer allocator.free(output_str);

    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    try file.writeStreamingAll(io, output_str);
}

/// Write SasaResult to file with specified format, using rich CSV when input has structural info
pub fn writeSasaResultWithFormatAndInput(
    allocator: Allocator,
    io: std.Io,
    result: SasaResult,
    input: AtomInput,
    path: []const u8,
    format: OutputFormat,
) !void {
    return writeSasaResultWithFormatAndInputOptions(allocator, io, result, input, path, format, .{});
}

/// Write SasaResult to file with specified format, using caller-provided metadata for text formats.
pub fn writeSasaResultWithFormatAndInputOptions(
    allocator: Allocator,
    io: std.Io,
    result: SasaResult,
    input: AtomInput,
    path: []const u8,
    format: OutputFormat,
    options: TextOutputOptions,
) !void {
    if (format == .rsa and try rsaResultNeedsLegacyWidthWarning(allocator, result, input)) {
        std.debug.print(
            "Warning: --format=rsa output exceeds legacy NACCESS fixed-width columns; columns may be misaligned. Use --format=json for machine-readable output.\n",
            .{},
        );
    }

    const output_str = switch (format) {
        .json => try sasaResultToJsonPretty(allocator, result),
        .compact => try sasaResultToJson(allocator, result),
        .csv => if (input.hasResidueInfo())
            try sasaResultToRichCsv(allocator, input, result.atom_areas)
        else
            try sasaResultToCsv(allocator, result),
        .freesasa => try sasaResultToFreesasa(allocator, result, options),
        .rsa => try sasaResultToRsa(allocator, result, input, options),
        .jsonl => unreachable, // JSONL is handled at batch level, not per-file
    };
    defer allocator.free(output_str);

    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    try file.writeStreamingAll(io, output_str);
}

/// Write SasaResult to JSON file (default: compact for backward compatibility)
pub fn writeSasaResult(allocator: Allocator, io: std.Io, result: SasaResult, path: []const u8) !void {
    return writeSasaResultWithFormat(allocator, io, result, path, .compact);
}

pub const ResidueMap = struct {
    allocator: Allocator,
    residue_chain: []const types.FixedString4,
    residue_chain_full: ?[]const []const u8 = null,
    residue_name: []const types.FixedString5,
    residue_number: []const i32,
    residue_insertion_code: []const types.FixedString4,
    residue_atom_start: []const usize,
    residue_atom_count: []const usize,
    residue_sasa: []const f64,

    pub fn len(self: ResidueMap) usize {
        return self.residue_chain.len;
    }

    pub fn deinit(self: *ResidueMap) void {
        self.allocator.free(self.residue_chain);
        if (self.residue_chain_full) |chains| {
            for (chains) |chain| self.allocator.free(chain);
            self.allocator.free(chains);
        }
        self.allocator.free(self.residue_name);
        self.allocator.free(self.residue_number);
        self.allocator.free(self.residue_insertion_code);
        self.allocator.free(self.residue_atom_start);
        self.allocator.free(self.residue_atom_count);
        self.allocator.free(self.residue_sasa);
        self.* = undefined;
    }
};

pub const JsonlOptions = struct {
    decimals: ?u8 = null,
    include_atom_areas: bool = true,
    include_atom_identity: bool = false,
    include_total_area: bool = true,
};

fn roundJsonlFloat(value: f64, decimals: u8) f64 {
    if (!std.math.isFinite(value)) return value;
    const factor = std.math.pow(f64, 10.0, @floatFromInt(decimals));
    const rounded = @round(value * factor) / factor;
    return if (rounded == 0.0) 0.0 else rounded;
}

fn maybeRoundJsonlFloat(value: f64, options: JsonlOptions) f64 {
    return if (options.decimals) |decimals| roundJsonlFloat(value, decimals) else value;
}

fn maybeRoundJsonlFloatSlice(allocator: Allocator, values: []const f64, options: JsonlOptions) ![]const f64 {
    const decimals = options.decimals orelse return values;
    const rounded = try allocator.alloc(f64, values.len);
    for (values, 0..) |value, i| {
        rounded[i] = roundJsonlFloat(value, decimals);
    }
    return rounded;
}

/// Build the residue map of the JSONL output. The residues are those of
/// `analysis.ResidueIdentity`: one entry per run of consecutive atoms with
/// the same chain ID, residue number, insertion code and residue name.
pub fn buildResidueMap(allocator: Allocator, input: AtomInput, atom_areas: []const f64) !ResidueMap {
    if (atom_areas.len != input.atomCount()) return error.LengthMismatch;

    const identity = try analysis.ResidueIdentity.init(input);
    const chain_ids_full = identity.chain_ids_full;
    const residue_count = identity.residueCount();

    const residue_chain = try allocator.alloc(types.FixedString4, residue_count);
    errdefer allocator.free(residue_chain);
    const residue_chain_full: ?[][]const u8 = if (chain_ids_full != null)
        try allocator.alloc([]const u8, residue_count)
    else
        null;
    errdefer if (residue_chain_full) |chains| {
        allocator.free(chains);
    };
    const residue_name = try allocator.alloc(types.FixedString5, residue_count);
    errdefer allocator.free(residue_name);
    const residue_number = try allocator.alloc(i32, residue_count);
    errdefer allocator.free(residue_number);
    const residue_insertion_code = try allocator.alloc(types.FixedString4, residue_count);
    errdefer allocator.free(residue_insertion_code);
    const residue_atom_start = try allocator.alloc(usize, residue_count);
    errdefer allocator.free(residue_atom_start);
    const residue_atom_count = try allocator.alloc(usize, residue_count);
    errdefer allocator.free(residue_atom_count);
    const residue_sasa = try allocator.alloc(f64, residue_count);
    errdefer allocator.free(residue_sasa);

    var residue_idx: usize = 0;
    var residue_full_initialized: usize = 0;
    errdefer if (residue_chain_full) |chains| {
        for (chains[0..residue_full_initialized]) |chain| allocator.free(chain);
    };
    var it = identity.residues();
    while (it.next()) |range| : (residue_idx += 1) {
        const start = range.start;
        var sasa = atom_areas[start];
        for (atom_areas[start + 1 .. range.end]) |area| sasa += area;

        residue_chain[residue_idx] = identity.chain_ids[start];
        if (residue_chain_full) |chains| {
            const full = chain_ids_full.?[start];
            chains[residue_idx] = try allocator.dupe(u8, full);
            residue_full_initialized += 1;
        }
        residue_name[residue_idx] = identity.residue_names[start];
        residue_number[residue_idx] = identity.residue_nums[start];
        residue_insertion_code[residue_idx] = identity.insertion_codes[start];
        residue_atom_start[residue_idx] = start;
        residue_atom_count[residue_idx] = range.atomCount();
        residue_sasa[residue_idx] = sasa;
    }

    return .{
        .allocator = allocator,
        .residue_chain = residue_chain,
        .residue_chain_full = residue_chain_full,
        .residue_name = residue_name,
        .residue_number = residue_number,
        .residue_insertion_code = residue_insertion_code,
        .residue_atom_start = residue_atom_start,
        .residue_atom_count = residue_atom_count,
        .residue_sasa = residue_sasa,
    };
}

/// Serialize a single batch result as a JSONL line: {"filename":"...","total_area":...,"atom_areas":[...]}
pub fn fileResultToJsonlLine(allocator: Allocator, filename: []const u8, total_area: f64, atom_areas: []const f64) ![]u8 {
    return fileResultToJsonlLineOptions(allocator, filename, total_area, atom_areas, .{});
}

pub fn fileResultToJsonlLineOptions(allocator: Allocator, filename: []const u8, total_area: f64, atom_areas: []const f64, options: JsonlOptions) ![]u8 {
    const output_areas = try maybeRoundJsonlFloatSlice(allocator, atom_areas, options);
    defer if (options.decimals != null) allocator.free(output_areas);

    if (options.include_total_area and !options.include_atom_areas) {
        const JsonlEntry = struct {
            status: []const u8,
            filename: []const u8,
            total_area: f64,
        };
        return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
            .status = "ok",
            .filename = filename,
            .total_area = maybeRoundJsonlFloat(total_area, options),
        }, .{});
    }

    if (!options.include_total_area and options.include_atom_areas) {
        const JsonlEntry = struct {
            status: []const u8,
            filename: []const u8,
            atom_areas: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
            .status = "ok",
            .filename = filename,
            .atom_areas = output_areas,
        }, .{});
    }

    if (!options.include_total_area and !options.include_atom_areas) {
        const JsonlEntry = struct {
            status: []const u8,
            filename: []const u8,
        };
        return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
            .status = "ok",
            .filename = filename,
        }, .{});
    }

    const JsonlEntry = struct {
        status: []const u8,
        filename: []const u8,
        total_area: f64,
        atom_areas: []const f64,
    };

    const entry = JsonlEntry{
        .status = "ok",
        .filename = filename,
        .total_area = maybeRoundJsonlFloat(total_area, options),
        .atom_areas = output_areas,
    };

    return std.json.Stringify.valueAlloc(allocator, entry, .{});
}

pub fn fileErrorToJsonlLine(allocator: Allocator, filename: []const u8, error_msg: []const u8) ![]u8 {
    const JsonlEntry = struct {
        status: []const u8,
        filename: []const u8,
        @"error": []const u8,
    };

    return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
        .status = "err",
        .filename = filename,
        .@"error" = error_msg,
    }, .{});
}

pub fn fileResultWithResidueMapToJsonlLine(
    allocator: Allocator,
    filename: []const u8,
    total_area: f64,
    atom_areas: []const f64,
    residue_map: ResidueMap,
) ![]u8 {
    return fileResultWithResidueMapToJsonlLineOptions(allocator, filename, total_area, atom_areas, residue_map, .{});
}

pub fn fileResultWithResidueMapToJsonlLineOptions(
    allocator: Allocator,
    filename: []const u8,
    total_area: f64,
    atom_areas: []const f64,
    residue_map: ResidueMap,
    options: JsonlOptions,
) ![]u8 {
    const residue_chain = try allocator.alloc([]const u8, residue_map.len());
    defer allocator.free(residue_chain);
    const residue_name = try allocator.alloc([]const u8, residue_map.len());
    defer allocator.free(residue_name);
    const residue_insertion_code = try allocator.alloc([]const u8, residue_map.len());
    defer allocator.free(residue_insertion_code);
    const output_areas = try maybeRoundJsonlFloatSlice(allocator, atom_areas, options);
    defer if (options.decimals != null) allocator.free(output_areas);
    const output_residue_sasa = try maybeRoundJsonlFloatSlice(allocator, residue_map.residue_sasa, options);
    defer if (options.decimals != null) allocator.free(output_residue_sasa);

    for (0..residue_map.len()) |i| {
        residue_chain[i] = if (residue_map.residue_chain_full) |chains|
            chains[i]
        else
            residue_map.residue_chain[i].slice();
        residue_name[i] = residue_map.residue_name[i].slice();
        residue_insertion_code[i] = residue_map.residue_insertion_code[i].slice();
    }

    if (options.include_total_area and options.include_atom_areas) {
        const JsonlEntry = struct {
            status: []const u8,
            filename: []const u8,
            total_area: f64,
            atom_areas: []const f64,
            residue_chain: []const []const u8,
            residue_name: []const []const u8,
            residue_number: []const i32,
            residue_insertion_code: []const []const u8,
            residue_atom_start: []const usize,
            residue_atom_count: []const usize,
            residue_sasa: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
            .status = "ok",
            .filename = filename,
            .total_area = maybeRoundJsonlFloat(total_area, options),
            .atom_areas = output_areas,
            .residue_chain = residue_chain,
            .residue_name = residue_name,
            .residue_number = residue_map.residue_number,
            .residue_insertion_code = residue_insertion_code,
            .residue_atom_start = residue_map.residue_atom_start,
            .residue_atom_count = residue_map.residue_atom_count,
            .residue_sasa = output_residue_sasa,
        }, .{});
    }

    if (options.include_total_area and !options.include_atom_areas) {
        const JsonlEntry = struct {
            status: []const u8,
            filename: []const u8,
            total_area: f64,
            residue_chain: []const []const u8,
            residue_name: []const []const u8,
            residue_number: []const i32,
            residue_insertion_code: []const []const u8,
            residue_atom_start: []const usize,
            residue_atom_count: []const usize,
            residue_sasa: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
            .status = "ok",
            .filename = filename,
            .total_area = maybeRoundJsonlFloat(total_area, options),
            .residue_chain = residue_chain,
            .residue_name = residue_name,
            .residue_number = residue_map.residue_number,
            .residue_insertion_code = residue_insertion_code,
            .residue_atom_start = residue_map.residue_atom_start,
            .residue_atom_count = residue_map.residue_atom_count,
            .residue_sasa = output_residue_sasa,
        }, .{});
    }

    if (!options.include_total_area and options.include_atom_areas) {
        const JsonlEntry = struct {
            status: []const u8,
            filename: []const u8,
            atom_areas: []const f64,
            residue_chain: []const []const u8,
            residue_name: []const []const u8,
            residue_number: []const i32,
            residue_insertion_code: []const []const u8,
            residue_atom_start: []const usize,
            residue_atom_count: []const usize,
            residue_sasa: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
            .status = "ok",
            .filename = filename,
            .atom_areas = output_areas,
            .residue_chain = residue_chain,
            .residue_name = residue_name,
            .residue_number = residue_map.residue_number,
            .residue_insertion_code = residue_insertion_code,
            .residue_atom_start = residue_map.residue_atom_start,
            .residue_atom_count = residue_map.residue_atom_count,
            .residue_sasa = output_residue_sasa,
        }, .{});
    }

    const JsonlEntry = struct {
        status: []const u8,
        filename: []const u8,
        residue_chain: []const []const u8,
        residue_name: []const []const u8,
        residue_number: []const i32,
        residue_insertion_code: []const []const u8,
        residue_atom_start: []const usize,
        residue_atom_count: []const usize,
        residue_sasa: []const f64,
    };
    return std.json.Stringify.valueAlloc(allocator, JsonlEntry{
        .status = "ok",
        .filename = filename,
        .residue_chain = residue_chain,
        .residue_name = residue_name,
        .residue_number = residue_map.residue_number,
        .residue_insertion_code = residue_insertion_code,
        .residue_atom_start = residue_map.residue_atom_start,
        .residue_atom_count = residue_map.residue_atom_count,
        .residue_sasa = output_residue_sasa,
    }, .{});
}

pub const SelectionAtomIdentity = struct {
    source_atom_index: []const usize,
    atom_chain: []const []const u8,
    atom_residue_name: []const []const u8,
    atom_residue_number: []const i32,
    atom_insertion_code: []const []const u8,
    atom_name: []const []const u8,
    atom_element: []const []const u8,
};

pub const SelectionResultJsonl = struct {
    filename: []const u8,
    id: []const u8,
    chains: []const []const u8,
    total_area: f64,
    atom_areas: []const f64 = &.{},
    residue_map: ?ResidueMap = null,
    atom_identity: ?SelectionAtomIdentity = null,
};

pub const SelectionErrorJsonl = struct {
    filename: []const u8,
    id: []const u8,
    chains: []const []const u8,
    error_message: []const u8,
};

pub fn selectionErrorToJsonlLine(allocator: Allocator, row: SelectionErrorJsonl) ![]u8 {
    const Entry = struct {
        status: []const u8,
        filename: []const u8,
        id: []const u8,
        chains: []const []const u8,
        @"error": []const u8,
    };
    return std.json.Stringify.valueAlloc(allocator, Entry{
        .status = "err",
        .filename = row.filename,
        .id = row.id,
        .chains = row.chains,
        .@"error" = row.error_message,
    }, .{});
}

pub fn selectionResultToJsonlLineOptions(
    allocator: Allocator,
    row: SelectionResultJsonl,
    options: JsonlOptions,
) ![]u8 {
    if (row.atom_identity != null and !options.include_atom_identity) return error.UnexpectedAtomIdentity;
    const output_areas = try maybeRoundJsonlFloatSlice(allocator, row.atom_areas, options);
    defer if (options.decimals != null) allocator.free(output_areas);

    var residue_chain: [][]const u8 = &.{};
    var residue_name: [][]const u8 = &.{};
    var residue_insertion_code: [][]const u8 = &.{};
    if (row.residue_map) |map| {
        residue_chain = try allocator.alloc([]const u8, map.len());
        residue_name = try allocator.alloc([]const u8, map.len());
        residue_insertion_code = try allocator.alloc([]const u8, map.len());
        for (0..map.len()) |i| {
            residue_chain[i] = if (map.residue_chain_full) |chains| chains[i] else map.residue_chain[i].slice();
            residue_name[i] = map.residue_name[i].slice();
            residue_insertion_code[i] = map.residue_insertion_code[i].slice();
        }
    }
    defer {
        if (row.residue_map != null) {
            allocator.free(residue_chain);
            allocator.free(residue_name);
            allocator.free(residue_insertion_code);
        }
    }
    const output_residue_sasa = if (row.residue_map) |map|
        try maybeRoundJsonlFloatSlice(allocator, map.residue_sasa, options)
    else
        &.{};
    defer if (row.residue_map != null and options.decimals != null) allocator.free(output_residue_sasa);

    if (row.atom_identity) |identity| {
        if (!options.include_atom_areas) return error.AtomIdentityRequiresAtomAreas;
        if (row.residue_map) |map| {
            const Entry = struct {
                status: []const u8,
                filename: []const u8,
                id: []const u8,
                chains: []const []const u8,
                total_area: f64,
                atom_areas: []const f64,
                source_atom_index: []const usize,
                atom_chain: []const []const u8,
                atom_residue_name: []const []const u8,
                atom_residue_number: []const i32,
                atom_insertion_code: []const []const u8,
                atom_name: []const []const u8,
                atom_element: []const []const u8,
                residue_chain: []const []const u8,
                residue_name: []const []const u8,
                residue_number: []const i32,
                residue_insertion_code: []const []const u8,
                residue_atom_start: []const usize,
                residue_atom_count: []const usize,
                residue_sasa: []const f64,
            };
            return std.json.Stringify.valueAlloc(allocator, Entry{
                .status = "ok",
                .filename = row.filename,
                .id = row.id,
                .chains = row.chains,
                .total_area = maybeRoundJsonlFloat(row.total_area, options),
                .atom_areas = output_areas,
                .source_atom_index = identity.source_atom_index,
                .atom_chain = identity.atom_chain,
                .atom_residue_name = identity.atom_residue_name,
                .atom_residue_number = identity.atom_residue_number,
                .atom_insertion_code = identity.atom_insertion_code,
                .atom_name = identity.atom_name,
                .atom_element = identity.atom_element,
                .residue_chain = residue_chain,
                .residue_name = residue_name,
                .residue_number = map.residue_number,
                .residue_insertion_code = residue_insertion_code,
                .residue_atom_start = map.residue_atom_start,
                .residue_atom_count = map.residue_atom_count,
                .residue_sasa = output_residue_sasa,
            }, .{});
        }
        const Entry = struct {
            status: []const u8,
            filename: []const u8,
            id: []const u8,
            chains: []const []const u8,
            total_area: f64,
            atom_areas: []const f64,
            source_atom_index: []const usize,
            atom_chain: []const []const u8,
            atom_residue_name: []const []const u8,
            atom_residue_number: []const i32,
            atom_insertion_code: []const []const u8,
            atom_name: []const []const u8,
            atom_element: []const []const u8,
        };
        return std.json.Stringify.valueAlloc(allocator, Entry{
            .status = "ok",
            .filename = row.filename,
            .id = row.id,
            .chains = row.chains,
            .total_area = maybeRoundJsonlFloat(row.total_area, options),
            .atom_areas = output_areas,
            .source_atom_index = identity.source_atom_index,
            .atom_chain = identity.atom_chain,
            .atom_residue_name = identity.atom_residue_name,
            .atom_residue_number = identity.atom_residue_number,
            .atom_insertion_code = identity.atom_insertion_code,
            .atom_name = identity.atom_name,
            .atom_element = identity.atom_element,
        }, .{});
    }

    if (row.residue_map) |map| {
        if (options.include_atom_areas) {
            const Entry = struct {
                status: []const u8,
                filename: []const u8,
                id: []const u8,
                chains: []const []const u8,
                total_area: f64,
                atom_areas: []const f64,
                residue_chain: []const []const u8,
                residue_name: []const []const u8,
                residue_number: []const i32,
                residue_insertion_code: []const []const u8,
                residue_atom_start: []const usize,
                residue_atom_count: []const usize,
                residue_sasa: []const f64,
            };
            return std.json.Stringify.valueAlloc(allocator, Entry{
                .status = "ok",
                .filename = row.filename,
                .id = row.id,
                .chains = row.chains,
                .total_area = maybeRoundJsonlFloat(row.total_area, options),
                .atom_areas = output_areas,
                .residue_chain = residue_chain,
                .residue_name = residue_name,
                .residue_number = map.residue_number,
                .residue_insertion_code = residue_insertion_code,
                .residue_atom_start = map.residue_atom_start,
                .residue_atom_count = map.residue_atom_count,
                .residue_sasa = output_residue_sasa,
            }, .{});
        }
        const Entry = struct {
            status: []const u8,
            filename: []const u8,
            id: []const u8,
            chains: []const []const u8,
            total_area: f64,
            residue_chain: []const []const u8,
            residue_name: []const []const u8,
            residue_number: []const i32,
            residue_insertion_code: []const []const u8,
            residue_atom_start: []const usize,
            residue_atom_count: []const usize,
            residue_sasa: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, Entry{
            .status = "ok",
            .filename = row.filename,
            .id = row.id,
            .chains = row.chains,
            .total_area = maybeRoundJsonlFloat(row.total_area, options),
            .residue_chain = residue_chain,
            .residue_name = residue_name,
            .residue_number = map.residue_number,
            .residue_insertion_code = residue_insertion_code,
            .residue_atom_start = map.residue_atom_start,
            .residue_atom_count = map.residue_atom_count,
            .residue_sasa = output_residue_sasa,
        }, .{});
    }

    if (options.include_atom_areas) {
        const Entry = struct {
            status: []const u8,
            filename: []const u8,
            id: []const u8,
            chains: []const []const u8,
            total_area: f64,
            atom_areas: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, Entry{
            .status = "ok",
            .filename = row.filename,
            .id = row.id,
            .chains = row.chains,
            .total_area = maybeRoundJsonlFloat(row.total_area, options),
            .atom_areas = output_areas,
        }, .{});
    }

    const Entry = struct {
        status: []const u8,
        filename: []const u8,
        id: []const u8,
        chains: []const []const u8,
        total_area: f64,
    };
    return std.json.Stringify.valueAlloc(allocator, Entry{
        .status = "ok",
        .filename = row.filename,
        .id = row.id,
        .chains = row.chains,
        .total_area = maybeRoundJsonlFloat(row.total_area, options),
    }, .{});
}

pub const BsaAnalysisJsonl = struct {
    filename: []const u8,
    id: []const u8,
    name: []const u8,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
    sasa_partner_a: f64,
    sasa_partner_b: f64,
    sasa_complex: f64,
    delta_sasa_total: f64,
    bsa: f64,
    delta_sasa_level: []const u8,
    atom_output: bool = false,
    residue_partner: []const []const u8 = &.{},
    residue_chain: []const []const u8 = &.{},
    residue_name: []const []const u8 = &.{},
    residue_number: []const i32 = &.{},
    residue_insertion_code: []const []const u8 = &.{},
    residue_sasa_isolated: []const f64 = &.{},
    residue_sasa_complex: []const f64 = &.{},
    residue_delta_sasa: []const f64 = &.{},
    atom_index: []const usize = &.{},
    atom_partner: []const []const u8 = &.{},
    atom_chain: []const []const u8 = &.{},
    atom_residue_name: []const []const u8 = &.{},
    atom_residue_number: []const i32 = &.{},
    atom_insertion_code: []const []const u8 = &.{},
    atom_name: []const []const u8 = &.{},
    atom_element: []const []const u8 = &.{},
    atom_sasa_isolated: []const f64 = &.{},
    atom_sasa_complex: []const f64 = &.{},
    atom_delta_sasa: []const f64 = &.{},
};

pub const BsaAnalysisErrorJsonl = struct {
    filename: []const u8,
    id: []const u8,
    name: []const u8,
    error_message: []const u8,
};

pub fn bsaAnalysisErrorToJsonlLine(allocator: Allocator, row: BsaAnalysisErrorJsonl) ![]u8 {
    const Entry = struct {
        status: []const u8,
        filename: []const u8,
        id: []const u8,
        analysis: []const u8,
        name: []const u8,
        @"error": []const u8,
    };
    return std.json.Stringify.valueAlloc(allocator, Entry{
        .status = "err",
        .filename = row.filename,
        .id = row.id,
        .analysis = "bsa",
        .name = row.name,
        .@"error" = row.error_message,
    }, .{});
}

pub fn bsaAnalysisToJsonlLine(allocator: Allocator, row: BsaAnalysisJsonl) ![]u8 {
    return bsaAnalysisToJsonlLineOptions(allocator, row, .{});
}

pub fn bsaAnalysisToJsonlLineOptions(allocator: Allocator, row: BsaAnalysisJsonl, options: JsonlOptions) ![]u8 {
    const output_residue_sasa_isolated = try maybeRoundJsonlFloatSlice(allocator, row.residue_sasa_isolated, options);
    defer if (options.decimals != null) allocator.free(output_residue_sasa_isolated);
    const output_residue_sasa_complex = try maybeRoundJsonlFloatSlice(allocator, row.residue_sasa_complex, options);
    defer if (options.decimals != null) allocator.free(output_residue_sasa_complex);
    const output_residue_delta_sasa = try maybeRoundJsonlFloatSlice(allocator, row.residue_delta_sasa, options);
    defer if (options.decimals != null) allocator.free(output_residue_delta_sasa);
    const output_atom_sasa_isolated = try maybeRoundJsonlFloatSlice(allocator, row.atom_sasa_isolated, options);
    defer if (options.decimals != null) allocator.free(output_atom_sasa_isolated);
    const output_atom_sasa_complex = try maybeRoundJsonlFloatSlice(allocator, row.atom_sasa_complex, options);
    defer if (options.decimals != null) allocator.free(output_atom_sasa_complex);
    const output_atom_delta_sasa = try maybeRoundJsonlFloatSlice(allocator, row.atom_delta_sasa, options);
    defer if (options.decimals != null) allocator.free(output_atom_delta_sasa);

    if (std.mem.eql(u8, row.delta_sasa_level, "residue")) {
        if (row.atom_output) {
            const Entry = struct {
                status: []const u8,
                filename: []const u8,
                id: []const u8,
                analysis: []const u8,
                name: []const u8,
                partner_a: []const []const u8,
                partner_b: []const []const u8,
                sasa_partner_a: f64,
                sasa_partner_b: f64,
                sasa_complex: f64,
                delta_sasa_total: f64,
                bsa: f64,
                delta_sasa_level: []const u8,
                residue_partner: []const []const u8,
                residue_chain: []const []const u8,
                residue_name: []const []const u8,
                residue_number: []const i32,
                residue_insertion_code: []const []const u8,
                residue_sasa_isolated: []const f64,
                residue_sasa_complex: []const f64,
                residue_delta_sasa: []const f64,
                atom_index: []const usize,
                atom_partner: []const []const u8,
                atom_chain: []const []const u8,
                atom_residue_name: []const []const u8,
                atom_residue_number: []const i32,
                atom_insertion_code: []const []const u8,
                atom_name: []const []const u8,
                atom_element: []const []const u8,
                atom_sasa_isolated: []const f64,
                atom_sasa_complex: []const f64,
                atom_delta_sasa: []const f64,
            };
            return std.json.Stringify.valueAlloc(allocator, Entry{
                .status = "ok",
                .filename = row.filename,
                .id = row.id,
                .analysis = "bsa",
                .name = row.name,
                .partner_a = row.partner_a,
                .partner_b = row.partner_b,
                .sasa_partner_a = maybeRoundJsonlFloat(row.sasa_partner_a, options),
                .sasa_partner_b = maybeRoundJsonlFloat(row.sasa_partner_b, options),
                .sasa_complex = maybeRoundJsonlFloat(row.sasa_complex, options),
                .delta_sasa_total = maybeRoundJsonlFloat(row.delta_sasa_total, options),
                .bsa = maybeRoundJsonlFloat(row.bsa, options),
                .delta_sasa_level = row.delta_sasa_level,
                .residue_partner = row.residue_partner,
                .residue_chain = row.residue_chain,
                .residue_name = row.residue_name,
                .residue_number = row.residue_number,
                .residue_insertion_code = row.residue_insertion_code,
                .residue_sasa_isolated = output_residue_sasa_isolated,
                .residue_sasa_complex = output_residue_sasa_complex,
                .residue_delta_sasa = output_residue_delta_sasa,
                .atom_index = row.atom_index,
                .atom_partner = row.atom_partner,
                .atom_chain = row.atom_chain,
                .atom_residue_name = row.atom_residue_name,
                .atom_residue_number = row.atom_residue_number,
                .atom_insertion_code = row.atom_insertion_code,
                .atom_name = row.atom_name,
                .atom_element = row.atom_element,
                .atom_sasa_isolated = output_atom_sasa_isolated,
                .atom_sasa_complex = output_atom_sasa_complex,
                .atom_delta_sasa = output_atom_delta_sasa,
            }, .{});
        }
        const Entry = struct {
            status: []const u8,
            filename: []const u8,
            id: []const u8,
            analysis: []const u8,
            name: []const u8,
            partner_a: []const []const u8,
            partner_b: []const []const u8,
            sasa_partner_a: f64,
            sasa_partner_b: f64,
            sasa_complex: f64,
            delta_sasa_total: f64,
            bsa: f64,
            delta_sasa_level: []const u8,
            residue_partner: []const []const u8,
            residue_chain: []const []const u8,
            residue_name: []const []const u8,
            residue_number: []const i32,
            residue_insertion_code: []const []const u8,
            residue_sasa_isolated: []const f64,
            residue_sasa_complex: []const f64,
            residue_delta_sasa: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, Entry{
            .status = "ok",
            .filename = row.filename,
            .id = row.id,
            .analysis = "bsa",
            .name = row.name,
            .partner_a = row.partner_a,
            .partner_b = row.partner_b,
            .sasa_partner_a = maybeRoundJsonlFloat(row.sasa_partner_a, options),
            .sasa_partner_b = maybeRoundJsonlFloat(row.sasa_partner_b, options),
            .sasa_complex = maybeRoundJsonlFloat(row.sasa_complex, options),
            .delta_sasa_total = maybeRoundJsonlFloat(row.delta_sasa_total, options),
            .bsa = maybeRoundJsonlFloat(row.bsa, options),
            .delta_sasa_level = row.delta_sasa_level,
            .residue_partner = row.residue_partner,
            .residue_chain = row.residue_chain,
            .residue_name = row.residue_name,
            .residue_number = row.residue_number,
            .residue_insertion_code = row.residue_insertion_code,
            .residue_sasa_isolated = output_residue_sasa_isolated,
            .residue_sasa_complex = output_residue_sasa_complex,
            .residue_delta_sasa = output_residue_delta_sasa,
        }, .{});
    }

    const Entry = struct {
        status: []const u8,
        filename: []const u8,
        id: []const u8,
        analysis: []const u8,
        name: []const u8,
        partner_a: []const []const u8,
        partner_b: []const []const u8,
        sasa_partner_a: f64,
        sasa_partner_b: f64,
        sasa_complex: f64,
        delta_sasa_total: f64,
        bsa: f64,
        delta_sasa_level: []const u8,
    };
    return std.json.Stringify.valueAlloc(allocator, Entry{
        .status = "ok",
        .filename = row.filename,
        .id = row.id,
        .analysis = "bsa",
        .name = row.name,
        .partner_a = row.partner_a,
        .partner_b = row.partner_b,
        .sasa_partner_a = maybeRoundJsonlFloat(row.sasa_partner_a, options),
        .sasa_partner_b = maybeRoundJsonlFloat(row.sasa_partner_b, options),
        .sasa_complex = maybeRoundJsonlFloat(row.sasa_complex, options),
        .delta_sasa_total = maybeRoundJsonlFloat(row.delta_sasa_total, options),
        .bsa = maybeRoundJsonlFloat(row.bsa, options),
        .delta_sasa_level = row.delta_sasa_level,
    }, .{});
}

// Tests

/// One atom of a `TestStructure`.
const TestAtom = struct {
    chain: []const u8,
    residue: []const u8,
    number: i32,
    insertion: []const u8 = "",
    atom: []const u8 = "CA",
    area: f64 = 1.0,
};

/// Structure input with residue metadata, built from a list of atoms.
const TestStructure = struct {
    arena: std.heap.ArenaAllocator,
    input: AtomInput,
    areas: []f64,

    /// `full_chain_ids` also fills `chain_id_full`, as the mmCIF parser does.
    fn init(atoms: []const TestAtom, full_chain_ids: bool) !TestStructure {
        var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
        errdefer arena.deinit();
        const allocator = arena.allocator();

        const n = atoms.len;
        const x = try allocator.alloc(f64, n);
        const y = try allocator.alloc(f64, n);
        const z = try allocator.alloc(f64, n);
        const r = try allocator.alloc(f64, n);
        const chain_id = try allocator.alloc(types.FixedString4, n);
        const chain_id_full = try allocator.alloc([]const u8, n);
        const residue = try allocator.alloc(types.FixedString5, n);
        const residue_num = try allocator.alloc(i32, n);
        const insertion_code = try allocator.alloc(types.FixedString4, n);
        const atom_name = try allocator.alloc(types.FixedString4, n);
        const areas = try allocator.alloc(f64, n);
        for (atoms, 0..) |atom, i| {
            x[i] = @floatFromInt(i);
            y[i] = 0;
            z[i] = 0;
            r[i] = 1;
            chain_id[i] = types.FixedString4.fromSlice(atom.chain);
            chain_id_full[i] = atom.chain;
            residue[i] = types.FixedString5.fromSlice(atom.residue);
            residue_num[i] = atom.number;
            insertion_code[i] = types.FixedString4.fromSlice(atom.insertion);
            atom_name[i] = types.FixedString4.fromSlice(atom.atom);
            areas[i] = atom.area;
        }

        return .{
            .arena = arena,
            .input = .{
                .x = x,
                .y = y,
                .z = z,
                .r = r,
                .chain_id = chain_id,
                .chain_id_full = if (full_chain_ids) chain_id_full else null,
                .residue = residue,
                .residue_num = residue_num,
                .insertion_code = insertion_code,
                .atom_name = atom_name,
                .allocator = allocator,
            },
            .areas = areas,
        };
    }

    fn deinit(self: *TestStructure) void {
        self.arena.deinit();
    }

    fn result(self: *const TestStructure) SasaResult {
        var total: f64 = 0;
        for (self.areas) |area| total += area;
        return .{ .total_area = total, .atom_areas = self.areas, .allocator = self.arena.child_allocator };
    }
};

const ExpectedResidue = struct {
    chain: []const u8,
    residue: []const u8,
    number: i32,
    insertion: []const u8 = "",
    atom_count: usize,
    sasa: f64,
};

/// The `--per-residue` table, the RSA file and the JSONL residue map must
/// report the same residues, in the same order, with the same atoms and areas.
fn expectSameResiduesInAllOutputs(atoms: []const TestAtom, full_chain_ids: bool, expected: []const ExpectedResidue) !void {
    const allocator = std.testing.allocator;
    var structure = try TestStructure.init(atoms, full_chain_ids);
    defer structure.deinit();

    // --per-residue / --rsa table
    var table = try analysis.aggregateByResidue(allocator, structure.input, structure.areas);
    defer table.deinit();
    try std.testing.expectEqual(expected.len, table.residues.len);
    for (expected, table.residues) |want, got| {
        try std.testing.expectEqualStrings(want.chain, got.chainLabel());
        try std.testing.expectEqualStrings(want.residue, got.residue_name.slice());
        try std.testing.expectEqual(want.number, got.residue_num);
        try std.testing.expectEqualStrings(want.insertion, got.insertion_code.slice());
        try std.testing.expectEqual(want.atom_count, got.atom_count);
        try std.testing.expectApproxEqAbs(want.sasa, got.sasa, 1e-9);
    }

    // JSONL residue map
    var map = try buildResidueMap(allocator, structure.input, structure.areas);
    defer map.deinit();
    try std.testing.expectEqual(expected.len, map.len());
    for (expected, 0..) |want, i| {
        const chain = if (map.residue_chain_full) |full| full[i] else map.residue_chain[i].slice();
        try std.testing.expectEqualStrings(want.chain, chain);
        try std.testing.expectEqualStrings(want.residue, map.residue_name[i].slice());
        try std.testing.expectEqual(want.number, map.residue_number[i]);
        try std.testing.expectEqualStrings(want.insertion, map.residue_insertion_code[i].slice());
        try std.testing.expectEqual(want.atom_count, map.residue_atom_count[i]);
        try std.testing.expectApproxEqAbs(want.sasa, map.residue_sasa[i], 1e-9);
    }

    // RSA file: one RES row per residue, in the same order
    const rsa = try sasaResultToRsa(allocator, structure.result(), structure.input, .{});
    defer allocator.free(rsa);
    var lines = std.mem.splitScalar(u8, rsa, '\n');
    var row: usize = 0;
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "RES ")) continue;
        try std.testing.expect(row < expected.len);
        const want = expected[row];
        var fields = std.mem.tokenizeScalar(u8, line[4..], ' ');
        try std.testing.expectEqualStrings(want.residue, fields.next().?);
        var number_buf: [32]u8 = undefined;
        const number = try std.fmt.bufPrint(&number_buf, "{d}{s}", .{ want.number, want.insertion });
        if (want.chain.len == 1) {
            // Chain and residue number share a field when the number is four characters wide
            const label = fields.next().?;
            try std.testing.expectEqualStrings(want.chain, label[0..1]);
            try std.testing.expectEqualStrings(number, if (label.len > 1) label[1..] else fields.next().?);
        } else {
            try std.testing.expectEqualStrings(want.chain, fields.next().?);
            try std.testing.expectEqualStrings(number, fields.next().?);
        }
        try std.testing.expectApproxEqAbs(want.sasa, try std.fmt.parseFloat(f64, fields.next().?), 0.005);
        row += 1;
    }
    try std.testing.expectEqual(expected.len, row);
}

test "all residue outputs keep residues that differ only in the insertion code apart" {
    try expectSameResiduesInAllOutputs(&.{
        .{ .chain = "H", .residue = "GLY", .number = 10, .atom = "N", .area = 1 },
        .{ .chain = "H", .residue = "GLY", .number = 10, .atom = "CA", .area = 2 },
        .{ .chain = "H", .residue = "GLY", .number = 10, .insertion = "A", .atom = "N", .area = 4 },
        .{ .chain = "H", .residue = "GLY", .number = 10, .insertion = "B", .atom = "N", .area = 8 },
        .{ .chain = "H", .residue = "GLY", .number = 10, .insertion = "B", .atom = "CA", .area = 16 },
    }, false, &.{
        .{ .chain = "H", .residue = "GLY", .number = 10, .atom_count = 2, .sasa = 3 },
        .{ .chain = "H", .residue = "GLY", .number = 10, .insertion = "A", .atom_count = 1, .sasa = 4 },
        .{ .chain = "H", .residue = "GLY", .number = 10, .insertion = "B", .atom_count = 2, .sasa = 24 },
    });
}

test "all residue outputs keep residues with the same number and different names apart" {
    try expectSameResiduesInAllOutputs(&.{
        .{ .chain = "A", .residue = "GLY", .number = 10, .atom = "N", .area = 1 },
        .{ .chain = "A", .residue = "GLY", .number = 10, .atom = "CA", .area = 2 },
        .{ .chain = "A", .residue = "LYS", .number = 10, .atom = "N", .area = 4 },
    }, false, &.{
        .{ .chain = "A", .residue = "GLY", .number = 10, .atom_count = 2, .sasa = 3 },
        .{ .chain = "A", .residue = "LYS", .number = 10, .atom_count = 1, .sasa = 4 },
    });
}

test "all residue outputs report a non-contiguous residue once per run" {
    try expectSameResiduesInAllOutputs(&.{
        .{ .chain = "A", .residue = "ALA", .number = 1, .atom = "N", .area = 1 },
        .{ .chain = "B", .residue = "UNK", .number = 2, .atom = "C1", .area = 2 },
        .{ .chain = "A", .residue = "ALA", .number = 1, .atom = "CB", .area = 4 },
    }, false, &.{
        .{ .chain = "A", .residue = "ALA", .number = 1, .atom_count = 1, .sasa = 1 },
        .{ .chain = "B", .residue = "UNK", .number = 2, .atom_count = 1, .sasa = 2 },
        .{ .chain = "A", .residue = "ALA", .number = 1, .atom_count = 1, .sasa = 4 },
    });
}

test "all residue outputs report the residues of superimposed models once per model" {
    // Two models of MET 1 - GLY 2, as the parsers return them for the default --model
    try expectSameResiduesInAllOutputs(&.{
        .{ .chain = "A", .residue = "MET", .number = 1, .atom = "N", .area = 1 },
        .{ .chain = "A", .residue = "MET", .number = 1, .atom = "CA", .area = 2 },
        .{ .chain = "A", .residue = "GLY", .number = 2, .atom = "N", .area = 4 },
        .{ .chain = "A", .residue = "MET", .number = 1, .atom = "N", .area = 8 },
        .{ .chain = "A", .residue = "MET", .number = 1, .atom = "CA", .area = 16 },
        .{ .chain = "A", .residue = "GLY", .number = 2, .atom = "N", .area = 32 },
    }, false, &.{
        .{ .chain = "A", .residue = "MET", .number = 1, .atom_count = 2, .sasa = 3 },
        .{ .chain = "A", .residue = "GLY", .number = 2, .atom_count = 1, .sasa = 4 },
        .{ .chain = "A", .residue = "MET", .number = 1, .atom_count = 2, .sasa = 24 },
        .{ .chain = "A", .residue = "GLY", .number = 2, .atom_count = 1, .sasa = 32 },
    });
}

const long_chain_atoms = [_]TestAtom{
    .{ .chain = "AAAAA", .residue = "ALA", .number = 1, .atom = "N", .area = 1 },
    .{ .chain = "AAAAB", .residue = "ALA", .number = 1, .atom = "N", .area = 2 },
    .{ .chain = "AAAA", .residue = "ALA", .number = 1, .atom = "N", .area = 4 },
};

test "all residue outputs keep chains apart whose IDs share the first four characters" {
    try expectSameResiduesInAllOutputs(&long_chain_atoms, true, &.{
        .{ .chain = "AAAAA", .residue = "ALA", .number = 1, .atom_count = 1, .sasa = 1 },
        .{ .chain = "AAAAB", .residue = "ALA", .number = 1, .atom_count = 1, .sasa = 2 },
        .{ .chain = "AAAA", .residue = "ALA", .number = 1, .atom_count = 1, .sasa = 4 },
    });
}

test "sasaResultToRsa writes one CHAIN line per full chain ID" {
    const allocator = std.testing.allocator;
    var structure = try TestStructure.init(&long_chain_atoms, true);
    defer structure.deinit();

    const rsa = try sasaResultToRsa(allocator, structure.result(), structure.input, .{});
    defer allocator.free(rsa);

    var chain_lines: usize = 0;
    var lines = std.mem.splitScalar(u8, rsa, '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "CHAIN")) continue;
        var fields = std.mem.tokenizeScalar(u8, line["CHAIN".len..], ' ');
        try std.testing.expectEqual(chain_lines + 1, try std.fmt.parseInt(usize, fields.next().?, 10));
        try std.testing.expectEqualStrings(long_chain_atoms[chain_lines].chain, fields.next().?);
        try std.testing.expectApproxEqAbs(long_chain_atoms[chain_lines].area, try std.fmt.parseFloat(f64, fields.next().?), 0.05);
        chain_lines += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), chain_lines);
}

test "sasaResultToRichCsv writes full chain IDs" {
    const allocator = std.testing.allocator;
    var structure = try TestStructure.init(&long_chain_atoms, true);
    defer structure.deinit();

    const csv = try sasaResultToRichCsv(allocator, structure.input, structure.areas);
    defer allocator.free(csv);

    var lines = std.mem.splitScalar(u8, csv, '\n');
    _ = lines.next(); // header
    for (long_chain_atoms) |atom| {
        var fields = std.mem.splitScalar(u8, lines.next().?, ',');
        try std.testing.expectEqualStrings(atom.chain, fields.next().?);
    }
}

test "buildResidueMap groups consecutive atoms" {
    const allocator = std.testing.allocator;

    const x = [_]f64{ 0, 1, 2, 3, 4 };
    const y = [_]f64{ 0, 0, 0, 0, 0 };
    const z = [_]f64{ 0, 0, 0, 0, 0 };
    var r = [_]f64{ 1, 1, 1, 1, 1 };
    const chain = [_]types.FixedString4{
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("B"),
    };
    const residue = [_]types.FixedString5{
        types.FixedString5.fromSlice("MET"),
        types.FixedString5.fromSlice("MET"),
        types.FixedString5.fromSlice("GLY"),
        types.FixedString5.fromSlice("GLY"),
        types.FixedString5.fromSlice("ALA"),
    };
    const residue_num = [_]i32{ 1, 1, 2, 2, 7 };
    const insertion = [_]types.FixedString4{
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice(""),
    };
    const atom_areas = [_]f64{ 10.0, 2.5, 4.0, 6.0, 1.25 };

    const input = types.AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .allocator = allocator,
    };

    var map = try buildResidueMap(allocator, input, atom_areas[0..]);
    defer map.deinit();

    try std.testing.expectEqual(@as(usize, 3), map.len());
    try std.testing.expectEqualStrings("A", map.residue_chain[0].slice());
    try std.testing.expectEqualStrings("MET", map.residue_name[0].slice());
    try std.testing.expectEqual(@as(i32, 1), map.residue_number[0]);
    try std.testing.expectEqualStrings("", map.residue_insertion_code[0].slice());
    try std.testing.expectEqual(@as(usize, 0), map.residue_atom_start[0]);
    try std.testing.expectEqual(@as(usize, 2), map.residue_atom_count[0]);
    try std.testing.expectApproxEqAbs(@as(f64, 12.5), map.residue_sasa[0], 1e-9);

    try std.testing.expectEqualStrings("A", map.residue_chain[1].slice());
    try std.testing.expectEqualStrings("GLY", map.residue_name[1].slice());
    try std.testing.expectEqual(@as(i32, 2), map.residue_number[1]);
    try std.testing.expectEqualStrings("A", map.residue_insertion_code[1].slice());
    try std.testing.expectEqual(@as(usize, 2), map.residue_atom_start[1]);
    try std.testing.expectEqual(@as(usize, 2), map.residue_atom_count[1]);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), map.residue_sasa[1], 1e-9);
}

test "buildResidueMap keeps non-contiguous repeated residues as separate ranges" {
    const allocator = std.testing.allocator;

    const x = [_]f64{ 0, 1, 2 };
    const y = [_]f64{ 0, 0, 0 };
    const z = [_]f64{ 0, 0, 0 };
    var r = [_]f64{ 1, 1, 1 };
    const chain = [_]types.FixedString4{
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
    };
    const residue = [_]types.FixedString5{
        types.FixedString5.fromSlice("MET"),
        types.FixedString5.fromSlice("GLY"),
        types.FixedString5.fromSlice("MET"),
    };
    const residue_num = [_]i32{ 1, 2, 1 };
    const insertion = [_]types.FixedString4{
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
    };
    const atom_areas = [_]f64{ 1.0, 2.0, 3.0 };

    const input = types.AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .allocator = allocator,
    };

    var map = try buildResidueMap(allocator, input, atom_areas[0..]);
    defer map.deinit();

    try std.testing.expectEqual(@as(usize, 3), map.len());
    try std.testing.expectEqual(@as(usize, 0), map.residue_atom_start[0]);
    try std.testing.expectEqual(@as(usize, 1), map.residue_atom_count[0]);
    try std.testing.expectEqual(@as(usize, 2), map.residue_atom_start[2]);
    try std.testing.expectEqual(@as(usize, 1), map.residue_atom_count[2]);
    try std.testing.expectEqualStrings("MET", map.residue_name[0].slice());
    try std.testing.expectEqualStrings("MET", map.residue_name[2].slice());
}

test "buildResidueMap groups and outputs full chain IDs when present" {
    const allocator = std.testing.allocator;

    const x = [_]f64{ 0, 1 };
    const y = [_]f64{ 0, 0 };
    const z = [_]f64{ 0, 0 };
    var r = [_]f64{ 1, 1 };
    const chain = [_]types.FixedString4{
        types.FixedString4.fromSlice("ABCD"),
        types.FixedString4.fromSlice("ABCD"),
    };
    const chain_full = [_][]const u8{ "ABCD1", "ABCD2" };
    const residue = [_]types.FixedString5{
        types.FixedString5.fromSlice("ALA"),
        types.FixedString5.fromSlice("ALA"),
    };
    const residue_num = [_]i32{ 1, 1 };
    const insertion = [_]types.FixedString4{
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
    };
    const atom_areas = [_]f64{ 10.0, 20.0 };

    const input = types.AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .chain_id_full = chain_full[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .allocator = allocator,
    };

    var map = try buildResidueMap(allocator, input, atom_areas[0..]);
    defer map.deinit();

    try std.testing.expectEqual(@as(usize, 2), map.len());

    const json = try fileResultWithResidueMapToJsonlLine(allocator, "long.cif", 30.0, atom_areas[0..], map);
    defer allocator.free(json);

    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, json, .{});
    defer parsed.deinit();
    const chains = parsed.value.object.get("residue_chain").?.array;
    try std.testing.expectEqualStrings("ABCD1", chains.items[0].string);
    try std.testing.expectEqualStrings("ABCD2", chains.items[1].string);
}

test "fileResultWithResidueMapToJsonlLine serializes columnar residue arrays" {
    const allocator = std.testing.allocator;

    const atom_areas = [_]f64{ 10.0, 2.5, 1.25 };
    const residue_chain = [_]types.FixedString4{
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("B"),
    };
    const residue_name = [_]types.FixedString5{
        types.FixedString5.fromSlice("MET"),
        types.FixedString5.fromSlice("ALA"),
    };
    const residue_number = [_]i32{ 1, 7 };
    const residue_insertion_code = [_]types.FixedString4{
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
    };
    const residue_atom_start = [_]usize{ 0, 2 };
    const residue_atom_count = [_]usize{ 2, 1 };
    const residue_sasa = [_]f64{ 12.5, 1.25 };

    const map = ResidueMap{
        .allocator = allocator,
        .residue_chain = residue_chain[0..],
        .residue_name = residue_name[0..],
        .residue_number = residue_number[0..],
        .residue_insertion_code = residue_insertion_code[0..],
        .residue_atom_start = residue_atom_start[0..],
        .residue_atom_count = residue_atom_count[0..],
        .residue_sasa = residue_sasa[0..],
    };

    const line = try fileResultWithResidueMapToJsonlLine(allocator, "example.cif", 13.75, atom_areas[0..], map);
    defer allocator.free(line);

    try std.testing.expectEqualStrings(
        "{\"status\":\"ok\",\"filename\":\"example.cif\",\"total_area\":13.75,\"atom_areas\":[10,2.5,1.25],\"residue_chain\":[\"A\",\"B\"],\"residue_name\":[\"MET\",\"ALA\"],\"residue_number\":[1,7],\"residue_insertion_code\":[\"\",\"\"],\"residue_atom_start\":[0,2],\"residue_atom_count\":[2,1],\"residue_sasa\":[12.5,1.25]}",
        line,
    );
}

test "BSA analysis JSONL includes total and residue delta fields" {
    const allocator = std.testing.allocator;
    const partner_a = [_][]const u8{"A"};
    const partner_b = [_][]const u8{"B"};
    const residue_chain = [_][]const u8{ "A", "B" };
    const residue_name = [_][]const u8{ "GLY", "ALA" };
    const residue_number = [_]i32{ 1, 2 };
    const residue_insertion_code = [_][]const u8{ "", "" };
    const residue_partner = [_][]const u8{ "a", "b" };
    const residue_sasa_isolated = [_]f64{ 7.0, 9.0 };
    const residue_sasa_complex = [_]f64{ 4.0, 4.0 };
    const residue_delta_sasa = [_]f64{ 3.0, 5.0 };

    const line = try bsaAnalysisToJsonlLine(allocator, .{
        .filename = "tiny.pdb",
        .id = "interaction-001",
        .name = "interface_ab",
        .partner_a = partner_a[0..],
        .partner_b = partner_b[0..],
        .sasa_partner_a = 10.0,
        .sasa_partner_b = 20.0,
        .sasa_complex = 14.0,
        .delta_sasa_total = 16.0,
        .bsa = 8.0,
        .delta_sasa_level = "residue",
        .residue_partner = residue_partner[0..],
        .residue_chain = residue_chain[0..],
        .residue_name = residue_name[0..],
        .residue_number = residue_number[0..],
        .residue_insertion_code = residue_insertion_code[0..],
        .residue_sasa_isolated = residue_sasa_isolated[0..],
        .residue_sasa_complex = residue_sasa_complex[0..],
        .residue_delta_sasa = residue_delta_sasa[0..],
    });
    defer allocator.free(line);

    try std.testing.expect(std.mem.indexOf(u8, line, "\"analysis\":\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"status\":\"ok\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"id\":\"interaction-001\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"delta_sasa_total\":16") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"bsa\":8") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_delta_sasa\":[3,5]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_partner\":[\"a\",\"b\"]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_sasa_isolated\":[7,9]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_sasa_complex\":[4,4]") != null);
}

test "BSA analysis JSONL rounds floats when decimals option is set" {
    const allocator = std.testing.allocator;
    const partner_a = [_][]const u8{"A"};
    const partner_b = [_][]const u8{"B"};
    const residue_delta_sasa = [_]f64{ 3.14159, 5.55555 };

    const line = try bsaAnalysisToJsonlLineOptions(allocator, .{
        .filename = "tiny.pdb",
        .id = "tiny.pdb",
        .name = "interface_ab",
        .partner_a = partner_a[0..],
        .partner_b = partner_b[0..],
        .sasa_partner_a = 10.12345,
        .sasa_partner_b = 20.98765,
        .sasa_complex = 14.44444,
        .delta_sasa_total = 16.66666,
        .bsa = 8.33333,
        .delta_sasa_level = "residue",
        .residue_delta_sasa = residue_delta_sasa[0..],
    }, .{ .decimals = 2 });
    defer allocator.free(line);

    try std.testing.expect(std.mem.indexOf(u8, line, "\"sasa_partner_a\":10.12") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"sasa_partner_b\":20.99") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"sasa_complex\":14.44") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"delta_sasa_total\":16.67") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"bsa\":8.33") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_delta_sasa\":[3.14,5.56]") != null);
}

test "BSA analysis JSONL includes opt-in atom detail" {
    const allocator = std.testing.allocator;
    const partner_a = [_][]const u8{"A"};
    const partner_b = [_][]const u8{"B"};
    const atom_index = [_]usize{ 0, 1 };
    const atom_partner = [_][]const u8{ "a", "b" };
    const atom_chain = [_][]const u8{ "A", "B" };
    const atom_residue_name = [_][]const u8{ "GLY", "ALA" };
    const atom_residue_number = [_]i32{ 1, 2 };
    const atom_insertion_code = [_][]const u8{ "", "" };
    const atom_name = [_][]const u8{ "CA", "N" };
    const atom_element = [_][]const u8{ "C", "N" };
    const atom_sasa_isolated = [_]f64{ 7.0, 9.0 };
    const atom_sasa_complex = [_]f64{ 4.0, 4.0 };
    const atom_delta_sasa = [_]f64{ 3.0, 5.0 };

    const line = try bsaAnalysisToJsonlLine(allocator, .{
        .filename = "tiny.pdb",
        .id = "atoms",
        .name = "interface_ab",
        .partner_a = partner_a[0..],
        .partner_b = partner_b[0..],
        .sasa_partner_a = 10.0,
        .sasa_partner_b = 20.0,
        .sasa_complex = 14.0,
        .delta_sasa_total = 16.0,
        .bsa = 8.0,
        .delta_sasa_level = "residue",
        .atom_output = true,
        .atom_index = atom_index[0..],
        .atom_partner = atom_partner[0..],
        .atom_chain = atom_chain[0..],
        .atom_residue_name = atom_residue_name[0..],
        .atom_residue_number = atom_residue_number[0..],
        .atom_insertion_code = atom_insertion_code[0..],
        .atom_name = atom_name[0..],
        .atom_element = atom_element[0..],
        .atom_sasa_isolated = atom_sasa_isolated[0..],
        .atom_sasa_complex = atom_sasa_complex[0..],
        .atom_delta_sasa = atom_delta_sasa[0..],
    });
    defer allocator.free(line);

    try std.testing.expect(std.mem.indexOf(u8, line, "\"atom_partner\":[\"a\",\"b\"]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"atom_element\":[\"C\",\"N\"]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"atom_sasa_isolated\":[7,9]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"atom_sasa_complex\":[4,4]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"atom_delta_sasa\":[3,5]") != null);
}

test "BSA analysis error JSONL includes stable ID" {
    const line = try bsaAnalysisErrorToJsonlLine(std.testing.allocator, .{
        .filename = "bad.cif",
        .id = "interaction-bad",
        .name = "interfaces",
        .error_message = "read/parse failed: InvalidFormat",
    });
    defer std.testing.allocator.free(line);

    try std.testing.expectEqualStrings(
        "{\"status\":\"err\",\"filename\":\"bad.cif\",\"id\":\"interaction-bad\",\"analysis\":\"bsa\",\"name\":\"interfaces\",\"error\":\"read/parse failed: InvalidFormat\"}",
        line,
    );
}

test "selection JSONL includes stable ID chains and opt-in source atom identity" {
    const source_atom_index = [_]usize{ 0, 2 };
    const atom_chain = [_][]const u8{ "A", "C" };
    const atom_residue_name = [_][]const u8{ "GLY", "SER" };
    const atom_residue_number = [_]i32{ 1, 3 };
    const atom_insertion_code = [_][]const u8{ "", "" };
    const atom_name = [_][]const u8{ "N", "CA" };
    const atom_element = [_][]const u8{ "N", "C" };
    const chains = [_][]const u8{ "A", "C" };
    const areas = [_]f64{ 10.123, 20.456 };

    const line = try selectionResultToJsonlLineOptions(std.testing.allocator, .{
        .filename = "complex.cif",
        .id = "ac",
        .chains = chains[0..],
        .total_area = 30.579,
        .atom_areas = areas[0..],
        .atom_identity = .{
            .source_atom_index = source_atom_index[0..],
            .atom_chain = atom_chain[0..],
            .atom_residue_name = atom_residue_name[0..],
            .atom_residue_number = atom_residue_number[0..],
            .atom_insertion_code = atom_insertion_code[0..],
            .atom_name = atom_name[0..],
            .atom_element = atom_element[0..],
        },
    }, .{ .decimals = 2, .include_atom_identity = true });
    defer std.testing.allocator.free(line);

    const parsed = try std.json.parseFromSlice(std.json.Value, std.testing.allocator, line, .{});
    defer parsed.deinit();
    const object = parsed.value.object;
    try std.testing.expectEqualStrings("ok", object.get("status").?.string);
    try std.testing.expectEqualStrings("ac", object.get("id").?.string);
    try std.testing.expectEqual(@as(usize, 2), object.get("chains").?.array.items.len);
    try std.testing.expectEqual(@as(i64, 2), object.get("source_atom_index").?.array.items[1].integer);
    try std.testing.expectEqualStrings("N", object.get("atom_element").?.array.items[0].string);
}

test "selection error JSONL preserves requested ID and chains" {
    const chains = [_][]const u8{"Z"};
    const line = try selectionErrorToJsonlLine(std.testing.allocator, .{
        .filename = "complex.cif",
        .id = "missing",
        .chains = chains[0..],
        .error_message = "selected chain not found: Z",
    });
    defer std.testing.allocator.free(line);

    try std.testing.expectEqualStrings(
        "{\"status\":\"err\",\"filename\":\"complex.cif\",\"id\":\"missing\",\"chains\":[\"Z\"],\"error\":\"selected chain not found: Z\"}",
        line,
    );
}

test "sasaResultToJson basic" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 2);
    defer allocator.free(atom_areas);

    atom_areas[0] = 32.47;
    atom_areas[1] = 0.25;

    const result = SasaResult{
        .total_area = 18923.28,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const json = try sasaResultToJson(allocator, result);
    defer allocator.free(json);

    try std.testing.expectEqualStrings(
        "{\"total_area\":18923.28,\"atom_areas\":[32.47,0.25]}",
        json,
    );
}

test "sasaResultToJson empty atoms" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 0);
    defer allocator.free(atom_areas);

    const result = SasaResult{
        .total_area = 0.0,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const json = try sasaResultToJson(allocator, result);
    defer allocator.free(json);

    try std.testing.expectEqualStrings(
        "{\"total_area\":0,\"atom_areas\":[]}",
        json,
    );
}

test "sasaResultToJson single atom" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 1);
    defer allocator.free(atom_areas);

    atom_areas[0] = 123.45;

    const result = SasaResult{
        .total_area = 123.45,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const json = try sasaResultToJson(allocator, result);
    defer allocator.free(json);

    try std.testing.expectEqualStrings(
        "{\"total_area\":123.45,\"atom_areas\":[123.45]}",
        json,
    );
}

test "sasaResultToJsonPretty basic" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 2);
    defer allocator.free(atom_areas);

    atom_areas[0] = 10.5;
    atom_areas[1] = 20.3;

    const result = SasaResult{
        .total_area = 30.8,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const json = try sasaResultToJsonPretty(allocator, result);
    defer allocator.free(json);

    const expected =
        \\{
        \\  "total_area": 30.8,
        \\  "atom_areas": [
        \\    10.5,
        \\    20.3
        \\  ]
        \\}
    ;

    try std.testing.expectEqualStrings(expected, json);
}

test "sasaResultToCsv basic" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 3);
    defer allocator.free(atom_areas);

    atom_areas[0] = 10.5;
    atom_areas[1] = 20.3;
    atom_areas[2] = 5.0;

    const result = SasaResult{
        .total_area = 35.8,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const csv = try sasaResultToCsv(allocator, result);
    defer allocator.free(csv);

    const expected =
        \\atom_index,area
        \\0,10.500000
        \\1,20.300000
        \\2,5.000000
        \\total,35.800000
        \\
    ;

    try std.testing.expectEqualStrings(expected, csv);
}

test "sasaResultToFreesasa writes FreeSASA-compatible text summary" {
    const allocator = std.testing.allocator;

    var atom_areas = [_]f64{ 10.0, 20.0, 30.0 };
    const result = SasaResult{
        .total_area = 60.0,
        .atom_areas = atom_areas[0..],
        .allocator = allocator,
    };

    const output = try sasaResultToFreesasa(allocator, result, .{
        .input_name = "mini.pdb",
        .classifier_name = "naccess",
        .algorithm_name = "Shrake & Rupley",
        .probe_radius = 1.4,
        .detail_count = 100,
        .detail_label = "Points",
    });
    defer allocator.free(output);

    try std.testing.expectEqualStrings(
        \\## zsasa FreeSASA-compatible output ##
        \\
        \\PARAMETERS
        \\algorithm    : Shrake & Rupley
        \\classifier   : naccess
        \\probe-radius : 1.40
        \\Points       : 100
        \\input        : mini.pdb
        \\
        \\RESULTS (A^2)
        \\Total   :      60.00
        \\
    , output);
}

test "sasaResultToRsa writes residue, chain, and total rows" {
    const allocator = std.testing.allocator;

    const x = [_]f64{ 0, 1, 2 };
    const y = [_]f64{ 0, 0, 0 };
    const z = [_]f64{ 0, 0, 0 };
    var r = [_]f64{ 1, 1, 1 };
    const chain = [_]types.FixedString4{
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("B"),
    };
    const residue = [_]types.FixedString5{
        types.FixedString5.fromSlice("ALA"),
        types.FixedString5.fromSlice("ALA"),
        types.FixedString5.fromSlice("UNK"),
    };
    const residue_num = [_]i32{ 1, 1, 2 };
    const insertion = [_]types.FixedString4{
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice("A"),
    };
    const atom_names = [_]types.FixedString4{
        types.FixedString4.fromSlice("N"),
        types.FixedString4.fromSlice("CB"),
        types.FixedString4.fromSlice("C1"),
    };
    var atom_areas = [_]f64{ 10.0, 20.0, 30.0 };
    const result = SasaResult{
        .total_area = 60.0,
        .atom_areas = atom_areas[0..],
        .allocator = allocator,
    };
    const input = AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .atom_name = atom_names[0..],
        .allocator = allocator,
    };

    const output = try sasaResultToRsa(allocator, result, input, .{
        .input_name = "mini.pdb",
        .classifier_name = "naccess",
        .algorithm_name = "Shrake & Rupley",
        .probe_radius = 1.4,
        .detail_count = 100,
        .detail_label = "Test-points",
    });
    defer allocator.free(output);

    try std.testing.expectEqualStrings(
        \\REM  zsasa FreeSASA/NACCESS-compatible RSA
        \\REM  Absolute and relative SASAs for mini.pdb
        \\REM  Atomic radii: naccess
        \\REM  Reference values for relative SASA: Tien et al. 2013
        \\REM  Algorithm: Shrake & Rupley
        \\REM  Probe-radius: 1.40
        \\REM  Test-points: 100
        \\REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
        \\REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
        \\RES ALA A   1    30.00  23.3  20.00   N/A  10.00   N/A  20.00   N/A  10.00   N/A
        \\RES UNK B   2A   30.00   N/A  30.00   N/A   0.00   N/A  30.00   N/A   0.00   N/A
        \\END  Absolute sums over single chains surface
        \\CHAIN  1 A       30.0         20.0         10.0         20.0         10.0
        \\CHAIN  2 B       30.0         30.0          0.0         30.0          0.0
        \\END  Absolute sums over all chains
        \\TOTAL            60.0         50.0         10.0         50.0         10.0
        \\
    , output);
}

test "sasaResultToRsa writes one row per run of a non-contiguous residue" {
    const allocator = std.testing.allocator;

    const x = [_]f64{ 0, 1, 2 };
    const y = [_]f64{ 0, 0, 0 };
    const z = [_]f64{ 0, 0, 0 };
    var r = [_]f64{ 1, 1, 1 };
    const chain = [_]types.FixedString4{
        types.FixedString4.fromSlice("A"),
        types.FixedString4.fromSlice("B"),
        types.FixedString4.fromSlice("A"),
    };
    const residue = [_]types.FixedString5{
        types.FixedString5.fromSlice("ALA"),
        types.FixedString5.fromSlice("UNK"),
        types.FixedString5.fromSlice("ALA"),
    };
    const residue_num = [_]i32{ 1, 2, 1 };
    const insertion = [_]types.FixedString4{
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
        types.FixedString4.fromSlice(""),
    };
    const atom_names = [_]types.FixedString4{
        types.FixedString4.fromSlice("N"),
        types.FixedString4.fromSlice("C1"),
        types.FixedString4.fromSlice("CB"),
    };
    var atom_areas = [_]f64{ 10.0, 30.0, 20.0 };
    const result = SasaResult{
        .total_area = 60.0,
        .atom_areas = atom_areas[0..],
        .allocator = allocator,
    };
    const input = AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .atom_name = atom_names[0..],
        .allocator = allocator,
    };

    const output = try sasaResultToRsa(allocator, result, input, .{
        .input_name = "mini.pdb",
        .classifier_name = "naccess",
        .algorithm_name = "Shrake & Rupley",
        .probe_radius = 1.4,
        .detail_count = 100,
        .detail_label = "Test-points",
    });
    defer allocator.free(output);

    // The same rows as the JSONL residue map; see analysis.ResidueIdentity
    try std.testing.expect(std.mem.indexOf(u8, output, "RES ALA A   1    10.00   7.8") != null);
    try std.testing.expect(std.mem.indexOf(u8, output, "RES ALA A   1    20.00  15.5") != null);
    try std.testing.expectEqual(@as(?usize, null), std.mem.indexOf(u8, output, "RES ALA A   1    30.00"));
    // Chain totals still cover every residue of the chain
    try std.testing.expect(std.mem.indexOf(u8, output, "CHAIN  1 A       30.0") != null);
}

/// The `RES`, `CHAIN` and `TOTAL` rows of the RSA text for `atoms`.
/// Caller frees the list and the text it points into.
const RsaRows = struct {
    text: []u8,
    rows: std.ArrayListUnmanaged([]const u8),
    needs_warning: bool,

    fn init(atoms: []const TestAtom, full_chain_ids: bool) !RsaRows {
        const allocator = std.testing.allocator;
        var structure = try TestStructure.init(atoms, full_chain_ids);
        defer structure.deinit();

        const text = try sasaResultToRsa(allocator, structure.result(), structure.input, .{});
        errdefer allocator.free(text);
        var rows = std.ArrayListUnmanaged([]const u8).empty;
        errdefer rows.deinit(allocator);
        var lines = std.mem.splitScalar(u8, text, '\n');
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "RES ") or std.mem.startsWith(u8, line, "CHAIN") or std.mem.startsWith(u8, line, "TOTAL")) {
                try rows.append(allocator, line);
            }
        }
        return .{
            .text = text,
            .rows = rows,
            .needs_warning = try rsaResultNeedsLegacyWidthWarning(allocator, structure.result(), structure.input),
        };
    }

    fn deinit(self: *RsaRows) void {
        self.rows.deinit(std.testing.allocator);
        std.testing.allocator.free(self.text);
    }
};

test "sasaResultToRsa rows follow the NACCESS fixed columns" {
    var rsa = try RsaRows.init(&.{
        .{ .chain = "A", .residue = "MET", .number = 1, .atom = "N", .area = 54.39 },
        .{ .chain = "A", .residue = "GLY", .number = -5, .atom = "N", .area = 0 },
        .{ .chain = "A", .residue = "SER", .number = 10, .insertion = "A", .atom = "OG", .area = 7.5 },
        .{ .chain = "B", .residue = "THR", .number = 1000, .atom = "CB", .area = 999.99 },
        .{ .chain = "B", .residue = "THR", .number = 1000, .insertion = "B", .atom = "N", .area = 100.004 },
        .{ .chain = "B", .residue = "ALA", .number = -999, .atom = "CB", .area = 12.346 },
        // Two-character residue name and a blank chain ID
        .{ .chain = "", .residue = "DA", .number = 7, .atom = "P", .area = 1 },
    }, false);
    defer rsa.deinit();

    // 0-based columns: residue name [4:7], chain [8], residue number [9:13],
    // insertion code [13], then five pairs of F7.2 and F6.1 from column 15
    try std.testing.expectEqualStrings(
        \\REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
        \\REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
        \\RES MET A   1    54.39  24.3   0.00   N/A  54.39   N/A   0.00   N/A  54.39   N/A
        \\RES GLY A  -5     0.00   0.0   0.00   N/A   0.00   N/A   0.00   N/A   0.00   N/A
        \\RES SER A  10A    7.50   4.8   7.50   N/A   0.00   N/A   0.00   N/A   7.50   N/A
        \\RES THR B1000   999.99 581.4 999.99   N/A   0.00   N/A 999.99   N/A   0.00   N/A
        \\RES THR B1000B  100.00  58.1   0.00   N/A 100.00   N/A   0.00   N/A 100.00   N/A
        \\RES ALA B-999    12.35   9.6  12.35   N/A   0.00   N/A  12.35   N/A   0.00   N/A
        \\RES  DA     7     1.00   N/A   1.00   N/A   0.00   N/A   0.00   N/A   1.00   N/A
        \\END  Absolute sums over single chains surface
        \\CHAIN  1 A       61.9          7.5         54.4          0.0         61.9
        \\CHAIN  2 B     1112.3       1012.3        100.0       1012.3        100.0
        \\CHAIN  3          1.0          1.0          0.0          0.0          1.0
        \\END  Absolute sums over all chains
        \\TOTAL          1175.2       1020.8        154.4       1012.3        162.9
        \\
    , rsa.text[std.mem.indexOf(u8, rsa.text, "REM RES").?..]);
    try std.testing.expect(!rsa.needs_warning);

    const Expected = struct { name: []const u8, chain: u8, number: []const u8, insertion: u8 };
    const expected = [_]Expected{
        .{ .name = "MET", .chain = 'A', .number = "   1", .insertion = ' ' },
        .{ .name = "GLY", .chain = 'A', .number = "  -5", .insertion = ' ' },
        .{ .name = "SER", .chain = 'A', .number = "  10", .insertion = 'A' },
        .{ .name = "THR", .chain = 'B', .number = "1000", .insertion = ' ' },
        .{ .name = "THR", .chain = 'B', .number = "1000", .insertion = 'B' },
        .{ .name = "ALA", .chain = 'B', .number = "-999", .insertion = ' ' },
        .{ .name = " DA", .chain = ' ', .number = "   7", .insertion = ' ' },
    };
    for (expected, rsa.rows.items[0..expected.len]) |want, row| {
        try std.testing.expectEqual(@as(usize, 80), row.len);
        try std.testing.expectEqualStrings("RES ", row[0..4]);
        try std.testing.expectEqualStrings(want.name, row[4..7]);
        try std.testing.expectEqual(@as(u8, ' '), row[7]);
        try std.testing.expectEqual(want.chain, row[8]);
        try std.testing.expectEqualStrings(want.number, row[9..13]);
        try std.testing.expectEqual(want.insertion, row[13]);
        try std.testing.expectEqual(@as(u8, ' '), row[14]);
        // Every value field starts with a blank and holds a number or N/A
        for (0..5) |pair| {
            const abs = row[15 + 13 * pair ..][0..7];
            const rel = row[22 + 13 * pair ..][0..6];
            try std.testing.expectEqual(@as(u8, ' '), abs[0]);
            try std.testing.expectEqual(@as(u8, ' '), rel[0]);
            _ = try std.fmt.parseFloat(f64, std.mem.trim(u8, abs, " "));
            if (!std.mem.eql(u8, rel, "   N/A")) _ = try std.fmt.parseFloat(f64, std.mem.trim(u8, rel, " "));
        }
    }

    // CHAIN and TOTAL rows: sums in [11:21], [24:34], [37:47], [50:60], [63:73]
    for (rsa.rows.items[expected.len..], [_]u8{ 'A', 'B', ' ', 0 }) |row, chain| {
        try std.testing.expectEqual(@as(usize, 73), row.len);
        if (chain != 0) {
            try std.testing.expectEqualStrings("CHAIN", row[0..5]);
            try std.testing.expectEqual(@as(u8, ' '), row[8]);
            try std.testing.expectEqual(chain, row[9]);
        } else {
            try std.testing.expectEqualStrings("TOTAL      ", row[0..11]);
        }
        for (0..5) |i| {
            _ = try std.fmt.parseFloat(f64, std.mem.trim(u8, row[11 + 13 * i ..][0..10], " "));
            try std.testing.expectEqual(@as(u8, ' '), row[10 + 13 * i]);
        }
    }
}

test "sasaResultToRsa keeps labels and values whole when they do not fit the fixed columns" {
    var rsa = try RsaRows.init(&.{
        // Residue number of five characters
        .{ .chain = "A", .residue = "ALA", .number = 12345, .atom = "CB", .area = 1 },
        // Chain ID of more than one character, with a four-digit residue number
        .{ .chain = "AAAAA", .residue = "ALA", .number = 1000, .atom = "CB", .area = 2 },
        // Residue name of more than three characters
        .{ .chain = "B", .residue = "A1LXQ", .number = 1, .atom = "C1", .area = 4 },
        // Absolute values that would touch the value before them
        .{ .chain = "C", .residue = "UNK", .number = 1, .atom = "C1", .area = 2222.27 },
    }, true);
    defer rsa.deinit();

    try std.testing.expect(rsa.needs_warning);
    try std.testing.expectEqualStrings("RES ALA A 12345     1.00   0.8   1.00   N/A   0.00   N/A   1.00   N/A   0.00   N/A", rsa.rows.items[0]);
    try std.testing.expectEqualStrings("RES ALA AAAAA 1000     2.00   1.6   2.00   N/A   0.00   N/A   2.00   N/A   0.00   N/A", rsa.rows.items[1]);
    try std.testing.expectEqualStrings("RES A1LXQ B    1     4.00   N/A   4.00   N/A   0.00   N/A   4.00   N/A   0.00   N/A", rsa.rows.items[2]);
    // `N/A2222.27` would be read as one value by a reader that splits at blanks
    try std.testing.expectEqualStrings("RES UNK C   1   2222.27   N/A 2222.27   N/A   0.00   N/A 2222.27   N/A   0.00   N/A", rsa.rows.items[3]);
    try std.testing.expectEqualStrings("CHAIN  2 AAAAA        2.0          2.0          0.0          2.0          0.0", rsa.rows.items[5]);

    // Every row can be split at blanks into its labels and ten values
    for (rsa.rows.items[0..4]) |row| {
        var fields = std.mem.tokenizeScalar(u8, row, ' ');
        var count: usize = 0;
        while (fields.next()) |_| count += 1;
        try std.testing.expectEqual(@as(usize, 14), count);
    }
}

test "rsaResultNeedsLegacyWidthWarning is true exactly when a row leaves the fixed columns" {
    const Case = struct { atom: TestAtom, full_chain_ids: bool = false, warns: bool };
    const cases = [_]Case{
        .{ .atom = .{ .chain = "A", .residue = "UNK", .number = 9999, .insertion = "Z", .area = 999.994 }, .warns = false },
        .{ .atom = .{ .chain = "A", .residue = "UNK", .number = -999, .area = 0 }, .warns = false },
        .{ .atom = .{ .chain = "", .residue = "DA", .number = 1, .area = 1 }, .warns = false },
        // 999.995 is printed as 1000.00
        .{ .atom = .{ .chain = "A", .residue = "UNK", .number = 1, .area = 999.996 }, .warns = true },
        .{ .atom = .{ .chain = "A", .residue = "UNK", .number = 10000, .area = 1 }, .warns = true },
        .{ .atom = .{ .chain = "A", .residue = "UNK", .number = -1000, .area = 1 }, .warns = true },
        .{ .atom = .{ .chain = "AB", .residue = "UNK", .number = 1, .area = 1 }, .warns = true },
        .{ .atom = .{ .chain = "AAAAA", .residue = "UNK", .number = 1, .area = 1 }, .full_chain_ids = true, .warns = true },
        .{ .atom = .{ .chain = "A", .residue = "UNKX", .number = 1, .area = 1 }, .warns = true },
        .{ .atom = .{ .chain = "A", .residue = "UNK", .number = 1, .insertion = "AB", .area = 1 }, .warns = true },
        // The largest relative value of an absolute value that fits (GLY has
        // the smallest maximum SASA, 104)
        .{ .atom = .{ .chain = "A", .residue = "GLY", .number = 1, .area = 999.99 }, .warns = false },
    };
    for (cases) |case| {
        var rsa = try RsaRows.init(&.{case.atom}, case.full_chain_ids);
        defer rsa.deinit();
        errdefer std.debug.print("row: {s}\n", .{rsa.rows.items[0]});
        try std.testing.expectEqual(case.warns, rsa.needs_warning);
        // A row that needs no warning is exactly 80 columns wide
        if (!case.warns) try std.testing.expectEqual(@as(usize, 80), rsa.rows.items[0].len);
    }
}

test "rsaResultNeedsLegacyWidthWarning detects oversized RSA numeric columns" {
    const allocator = std.testing.allocator;

    const x = [_]f64{0};
    const y = [_]f64{0};
    const z = [_]f64{0};
    var r = [_]f64{1};
    const chain = [_]types.FixedString4{types.FixedString4.fromSlice("A")};
    const residue = [_]types.FixedString5{types.FixedString5.fromSlice("ALA")};
    const residue_num = [_]i32{1};
    const insertion = [_]types.FixedString4{types.FixedString4.fromSlice("")};
    const atom_names = [_]types.FixedString4{types.FixedString4.fromSlice("CB")};
    var atom_areas = [_]f64{10000.0};
    const result = SasaResult{
        .total_area = 10000.0,
        .atom_areas = atom_areas[0..],
        .allocator = allocator,
    };
    const input = AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .atom_name = atom_names[0..],
        .allocator = allocator,
    };

    try std.testing.expect(try rsaResultNeedsLegacyWidthWarning(allocator, result, input));
}

test "rsaResultNeedsLegacyWidthWarning accepts legacy-width RSA output" {
    const allocator = std.testing.allocator;

    const x = [_]f64{0};
    const y = [_]f64{0};
    const z = [_]f64{0};
    var r = [_]f64{1};
    const chain = [_]types.FixedString4{types.FixedString4.fromSlice("A")};
    const residue = [_]types.FixedString5{types.FixedString5.fromSlice("ALA")};
    const residue_num = [_]i32{1};
    const insertion = [_]types.FixedString4{types.FixedString4.fromSlice("")};
    const atom_names = [_]types.FixedString4{types.FixedString4.fromSlice("CB")};
    var atom_areas = [_]f64{30.0};
    const result = SasaResult{
        .total_area = 30.0,
        .atom_areas = atom_areas[0..],
        .allocator = allocator,
    };
    const input = AtomInput{
        .x = x[0..],
        .y = y[0..],
        .z = z[0..],
        .r = r[0..],
        .chain_id = chain[0..],
        .residue = residue[0..],
        .residue_num = residue_num[0..],
        .insertion_code = insertion[0..],
        .atom_name = atom_names[0..],
        .allocator = allocator,
    };

    try std.testing.expect(!try rsaResultNeedsLegacyWidthWarning(allocator, result, input));
}

test "sasaResultToCsv empty atoms" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 0);
    defer allocator.free(atom_areas);

    const result = SasaResult{
        .total_area = 0.0,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const csv = try sasaResultToCsv(allocator, result);
    defer allocator.free(csv);

    const expected =
        \\atom_index,area
        \\total,0.000000
        \\
    ;

    try std.testing.expectEqualStrings(expected, csv);
}

/// Absolute path of `name` inside the temporary directory. Caller frees.
fn tmpFilePath(allocator: Allocator, tmp: *std.testing.TmpDir, name: []const u8) ![]u8 {
    var buf: [std.fs.max_path_bytes]u8 = undefined;
    const len = try tmp.dir.realPath(std.testing.io, &buf);
    return std.fs.path.join(allocator, &.{ buf[0..len], name });
}

test "writeSasaResult creates file" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;

    const atom_areas = try allocator.alloc(f64, 2);
    defer allocator.free(atom_areas);

    atom_areas[0] = 10.5;
    atom_areas[1] = 20.3;

    const result = SasaResult{
        .total_area = 30.8,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const test_path = try tmpFilePath(allocator, &tmp_dir, "test_output.json");
    defer allocator.free(test_path);

    try writeSasaResult(allocator, io, result, test_path);

    // Read back and verify
    const file = try std.Io.Dir.cwd().openFile(io, test_path, .{});
    defer file.close(io);

    var read_buf: [4096]u8 = undefined;
    var r = file.reader(io, &read_buf);
    const content = try r.interface.allocRemaining(allocator, .unlimited);
    defer allocator.free(content);

    try std.testing.expectEqualStrings(
        "{\"total_area\":30.8,\"atom_areas\":[10.5,20.3]}",
        content,
    );
}

test "writeSasaResultWithFormat json" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;

    const atom_areas = try allocator.alloc(f64, 2);
    defer allocator.free(atom_areas);

    atom_areas[0] = 10.5;
    atom_areas[1] = 20.3;

    const result = SasaResult{
        .total_area = 30.8,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const test_path = try tmpFilePath(allocator, &tmp_dir, "test_format_json.json");
    defer allocator.free(test_path);

    try writeSasaResultWithFormat(allocator, io, result, test_path, .json);

    const file = try std.Io.Dir.cwd().openFile(io, test_path, .{});
    defer file.close(io);

    var read_buf: [4096]u8 = undefined;
    var r = file.reader(io, &read_buf);
    const content = try r.interface.allocRemaining(allocator, .unlimited);
    defer allocator.free(content);

    // Should be pretty-printed
    try std.testing.expect(std.mem.find(u8, content, "\n") != null);
}

test "writeSasaResultWithFormat csv" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;

    const atom_areas = try allocator.alloc(f64, 2);
    defer allocator.free(atom_areas);

    atom_areas[0] = 10.5;
    atom_areas[1] = 20.3;

    const result = SasaResult{
        .total_area = 30.8,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const test_path = try tmpFilePath(allocator, &tmp_dir, "test_format.csv");
    defer allocator.free(test_path);

    try writeSasaResultWithFormat(allocator, io, result, test_path, .csv);

    const file = try std.Io.Dir.cwd().openFile(io, test_path, .{});
    defer file.close(io);

    var read_buf: [4096]u8 = undefined;
    var r = file.reader(io, &read_buf);
    const content = try r.interface.allocRemaining(allocator, .unlimited);
    defer allocator.free(content);

    // Should start with header
    try std.testing.expect(std.mem.startsWith(u8, content, "atom_index,area\n"));
}

test "writeSasaResult overwrites existing file" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;

    const atom_areas = try allocator.alloc(f64, 1);
    defer allocator.free(atom_areas);

    atom_areas[0] = 50.0;

    const result = SasaResult{
        .total_area = 50.0,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    const test_path = try tmpFilePath(allocator, &tmp_dir, "test_overwrite.json");
    defer allocator.free(test_path);

    // Write first time
    try writeSasaResult(allocator, io, result, test_path);

    // Write second time (overwrite)
    atom_areas[0] = 99.9;
    const result2 = SasaResult{
        .total_area = 99.9,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
    try writeSasaResult(allocator, io, result2, test_path);

    // Verify overwrite
    const file = try std.Io.Dir.cwd().openFile(io, test_path, .{});
    defer file.close(io);

    var read_buf: [4096]u8 = undefined;
    var r = file.reader(io, &read_buf);
    const content = try r.interface.allocRemaining(allocator, .unlimited);
    defer allocator.free(content);

    try std.testing.expectEqualStrings(
        "{\"total_area\":99.9,\"atom_areas\":[99.9]}",
        content,
    );
}

test "sasaResultToRichCsv with full info" {
    const allocator = std.testing.allocator;

    // Create coordinate arrays
    const x = try allocator.alloc(f64, 2);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 2);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 2);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 2);
    defer allocator.free(r);

    x[0] = 1.0;
    x[1] = 2.0;
    y[0] = 3.0;
    y[1] = 4.0;
    z[0] = 5.0;
    z[1] = 6.0;
    r[0] = 1.5;
    r[1] = 1.7;

    // Create metadata arrays
    const chain_ids = try allocator.alloc(types.FixedString4, 2);
    defer allocator.free(chain_ids);
    chain_ids[0] = types.FixedString4.fromSlice("A");
    chain_ids[1] = types.FixedString4.fromSlice("A");

    const residues = try allocator.alloc(types.FixedString5, 2);
    defer allocator.free(residues);
    residues[0] = types.FixedString5.fromSlice("ALA");
    residues[1] = types.FixedString5.fromSlice("ALA");

    const atom_names = try allocator.alloc(types.FixedString4, 2);
    defer allocator.free(atom_names);
    atom_names[0] = types.FixedString4.fromSlice("N");
    atom_names[1] = types.FixedString4.fromSlice("CA");

    const residue_nums = try allocator.alloc(i32, 2);
    defer allocator.free(residue_nums);
    residue_nums[0] = 1;
    residue_nums[1] = 1;

    const insertion_codes = try allocator.alloc(types.FixedString4, 2);
    defer allocator.free(insertion_codes);
    insertion_codes[0] = types.FixedString4.fromSlice("");
    insertion_codes[1] = types.FixedString4.fromSlice("");

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .chain_id = chain_ids,
        .residue = residues,
        .atom_name = atom_names,
        .residue_num = residue_nums,
        .insertion_code = insertion_codes,
        .allocator = allocator,
    };

    const atom_areas = try allocator.alloc(f64, 2);
    defer allocator.free(atom_areas);
    atom_areas[0] = 10.5;
    atom_areas[1] = 20.3;

    const csv = try sasaResultToRichCsv(allocator, input, atom_areas);
    defer allocator.free(csv);

    // Header, one row per atom (the insertion code column is empty), total row
    try std.testing.expectEqualStrings(
        \\chain,residue,resnum,insertion_code,atom_name,x,y,z,radius,area
        \\A,ALA,1,,N,1.000,3.000,5.000,1.500,10.500000
        \\A,ALA,1,,CA,2.000,4.000,6.000,1.700,20.300000
        \\,,,,,,,,,30.800000
        \\
    , csv);
}

test "sasaResultToRichCsv writes the insertion code after the residue number" {
    const allocator = std.testing.allocator;
    var structure = try TestStructure.init(&.{
        .{ .chain = "H", .residue = "GLY", .number = 10, .atom = "N", .area = 1 },
        .{ .chain = "H", .residue = "SER", .number = 10, .insertion = "A", .atom = "N", .area = 2 },
        .{ .chain = "H", .residue = "THR", .number = 10, .insertion = "B", .atom = "OG1", .area = 4 },
        .{ .chain = "H", .residue = "ALA", .number = -3, .atom = "CB", .area = 8 },
    }, false);
    defer structure.deinit();

    const csv = try sasaResultToRichCsv(allocator, structure.input, structure.areas);
    defer allocator.free(csv);

    try std.testing.expectEqualStrings(
        \\chain,residue,resnum,insertion_code,atom_name,x,y,z,radius,area
        \\H,GLY,10,,N,0.000,0.000,0.000,1.000,1.000000
        \\H,SER,10,A,N,1.000,0.000,0.000,1.000,2.000000
        \\H,THR,10,B,OG1,2.000,0.000,0.000,1.000,4.000000
        \\H,ALA,-3,,CB,3.000,0.000,0.000,1.000,8.000000
        \\,,,,,,,,,15.000000
        \\
    , csv);

    // Every row has the ten fields of the header
    var lines = std.mem.tokenizeScalar(u8, csv, '\n');
    while (lines.next()) |line| {
        try std.testing.expectEqual(@as(usize, 9), std.mem.count(u8, line, ","));
    }
}

test "sasaResultToRichCsv quotes fields per RFC 4180" {
    const allocator = std.testing.allocator;
    var structure = try TestStructure.init(&.{
        // Nothing to quote: written exactly as without quoting support
        .{ .chain = "A", .residue = "ALA", .number = 1, .atom = "CA", .area = 1 },
        .{ .chain = " ", .residue = "A B", .number = 2, .atom = "C'", .area = 2 },
        // Comma, double quote, LF and CR
        .{ .chain = ",", .residue = "a,b", .number = 3, .insertion = ",", .atom = "C,1", .area = 4 },
        .{ .chain = "\"", .residue = "a\"b", .number = 4, .insertion = "\"", .atom = "\"\"", .area = 8 },
        .{ .chain = "A", .residue = "a\nb", .number = 5, .atom = "C\r1", .area = 16 },
    }, false);
    defer structure.deinit();

    const csv = try sasaResultToRichCsv(allocator, structure.input, structure.areas);
    defer allocator.free(csv);

    try std.testing.expectEqualStrings("chain,residue,resnum,insertion_code,atom_name,x,y,z,radius,area\n" ++
        "A,ALA,1,,CA,0.000,0.000,0.000,1.000,1.000000\n" ++
        " ,A B,2,,C',1.000,0.000,0.000,1.000,2.000000\n" ++
        "\",\",\"a,b\",3,\",\",\"C,1\",2.000,0.000,0.000,1.000,4.000000\n" ++
        "\"\"\"\",\"a\"\"b\",4,\"\"\"\",\"\"\"\"\"\",3.000,0.000,0.000,1.000,8.000000\n" ++
        "A,\"a\nb\",5,,\"C\r1\",4.000,0.000,0.000,1.000,16.000000\n" ++
        ",,,,,,,,,31.000000\n", csv);
}

test "sasaResultToRichCsv quotes a full chain ID" {
    const allocator = std.testing.allocator;
    var structure = try TestStructure.init(&.{
        .{ .chain = "A,\"long\"", .residue = "ALA", .number = 1, .atom = "CA", .area = 1 },
    }, true);
    defer structure.deinit();

    const csv = try sasaResultToRichCsv(allocator, structure.input, structure.areas);
    defer allocator.free(csv);

    try std.testing.expect(std.mem.find(u8, csv, "\n\"A,\"\"long\"\"\",ALA,1,,CA,") != null);
}

test "sasaResultToRichCsv without residue info uses dashes" {
    const allocator = std.testing.allocator;

    // Create coordinate arrays only (no metadata)
    const x = try allocator.alloc(f64, 1);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 1);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 1);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 1);
    defer allocator.free(r);

    x[0] = 1.0;
    y[0] = 2.0;
    z[0] = 3.0;
    r[0] = 1.5;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const atom_areas = try allocator.alloc(f64, 1);
    defer allocator.free(atom_areas);
    atom_areas[0] = 15.0;

    const csv = try sasaResultToRichCsv(allocator, input, atom_areas);
    defer allocator.free(csv);

    // Check that missing fields produce dashes
    try std.testing.expect(std.mem.find(u8, csv, "\n-,-,-,-,-,1.000,2.000,3.000,1.500,15.000000\n") != null);
}

test "fileResultToJsonlLine basic" {
    const allocator = std.testing.allocator;
    const areas = [_]f64{ 1.5, 2.0, 0.0 };
    const line = try fileResultToJsonlLine(allocator, "test.pdb", 6.789, &areas);
    defer allocator.free(line);

    // Parse back to verify valid JSON
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
    defer parsed.deinit();

    const obj = parsed.value.object;
    try std.testing.expectEqualStrings("ok", obj.get("status").?.string);
    try std.testing.expectEqualStrings("test.pdb", obj.get("filename").?.string);
}

test "fileResultToJsonlLine rounds floats when decimals option is set" {
    const allocator = std.testing.allocator;
    const areas = [_]f64{ 1.23456, 2.55555, 0.004 };
    const line = try fileResultToJsonlLineOptions(allocator, "test.pdb", 6.789, &areas, .{ .decimals = 2 });
    defer allocator.free(line);

    try std.testing.expectEqualStrings(
        "{\"status\":\"ok\",\"filename\":\"test.pdb\",\"total_area\":6.79,\"atom_areas\":[1.23,2.56,0]}",
        line,
    );
}

test "fileResultToJsonlLine can omit atom areas" {
    const allocator = std.testing.allocator;
    const areas = [_]f64{ 1.0, 2.0 };
    const line = try fileResultToJsonlLineOptions(allocator, "summary.pdb", 3.0, &areas, .{
        .include_atom_areas = false,
    });
    defer allocator.free(line);

    try std.testing.expectEqualStrings(
        "{\"status\":\"ok\",\"filename\":\"summary.pdb\",\"total_area\":3}",
        line,
    );
}

test "fileResultToJsonlLine can omit total area" {
    const allocator = std.testing.allocator;
    const areas = [_]f64{ 1.0, 2.0 };
    const line = try fileResultToJsonlLineOptions(allocator, "areas.pdb", 3.0, &areas, .{
        .include_total_area = false,
    });
    defer allocator.free(line);

    try std.testing.expectEqualStrings(
        "{\"status\":\"ok\",\"filename\":\"areas.pdb\",\"atom_areas\":[1,2]}",
        line,
    );
}

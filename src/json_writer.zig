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

const ResidueArea = struct {
    chain: types.FixedString4,
    residue: types.FixedString5,
    residue_number: i32,
    insertion_code: types.FixedString4,
    area: AreaBreakdown,
};

const ChainArea = struct {
    chain: types.FixedString4,
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

fn sameResidueKey(
    area: ResidueArea,
    chain: types.FixedString4,
    residue: types.FixedString5,
    residue_number: i32,
    insertion_code: types.FixedString4,
) bool {
    return area.residue_number == residue_number and
        std.mem.eql(u8, area.chain.slice(), chain.slice()) and
        std.mem.eql(u8, area.residue.slice(), residue.slice()) and
        std.mem.eql(u8, area.insertion_code.slice(), insertion_code.slice());
}

fn collectResidueAreas(allocator: Allocator, input: AtomInput, atom_areas: []const f64) ![]ResidueArea {
    if (!input.hasResidueInfo()) return error.MissingResidueInfo;
    if (atom_areas.len != input.atomCount()) return error.LengthMismatch;

    var residues = std.ArrayListUnmanaged(ResidueArea).empty;
    errdefer residues.deinit(allocator);

    const chains = input.chain_id.?;
    const residue_names = input.residue.?;
    const residue_nums = input.residue_num.?;
    const insertions = input.insertion_code.?;
    const atom_names = input.atom_name;
    const elements = input.element;

    for (0..input.atomCount()) |i| {
        const chain = chains[i];
        const residue = residue_names[i];
        const residue_num = residue_nums[i];
        const insertion = insertions[i];

        var found_idx: ?usize = null;
        for (residues.items, 0..) |residue_area, j| {
            if (sameResidueKey(residue_area, chain, residue, residue_num, insertion)) {
                found_idx = j;
                break;
            }
        }

        const idx = found_idx orelse blk: {
            try residues.append(allocator, .{
                .chain = chain,
                .residue = residue,
                .residue_number = residue_num,
                .insertion_code = insertion,
                .area = .{},
            });
            break :blk residues.items.len - 1;
        };

        addAtomArea(
            &residues.items[idx].area,
            atom_areas[i],
            if (atom_names) |names| names[i].slice() else null,
            if (elements) |elem| elem[i] else null,
        );
    }

    return residues.toOwnedSlice(allocator);
}

fn collectChainAreas(allocator: Allocator, residues: []const ResidueArea) ![]ChainArea {
    var chains = std.ArrayListUnmanaged(ChainArea).empty;
    errdefer chains.deinit(allocator);

    for (residues) |residue| {
        var found_idx: ?usize = null;
        for (chains.items, 0..) |*chain, i| {
            if (std.mem.eql(u8, chain.chain.slice(), residue.chain.slice())) {
                found_idx = i;
                break;
            }
        }

        const idx = found_idx orelse blk: {
            try chains.append(allocator, .{ .chain = residue.chain, .area = .{} });
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

fn writeAbsRel(writer: *std.Io.Writer, abs: f64, rel: ?f64) !void {
    try writer.print("{d:7.2}", .{abs});
    if (rel) |value| {
        try writer.print("{d:6.1}", .{value});
    } else {
        try writer.writeAll("   N/A");
    }
}

fn residueNumberString(buf: []u8, number: i32, insertion_code: types.FixedString4) []const u8 {
    if (insertion_code.len > 0) {
        return std.fmt.bufPrint(buf, "{d}{s}", .{ number, insertion_code.slice() }) catch "?";
    }
    return std.fmt.bufPrint(buf, "{d}", .{number}) catch "?";
}

pub fn sasaResultToRsa(allocator: Allocator, result: SasaResult, input: AtomInput, options: TextOutputOptions) ![]u8 {
    const residues = try collectResidueAreas(allocator, input, result.atom_areas);
    defer allocator.free(residues);
    const chains = try collectChainAreas(allocator, residues);
    defer allocator.free(chains);

    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    try writer.writeAll("REM  zsasa FreeSASA/NACCESS-compatible RSA\n");
    try writer.print("REM  Absolute and relative SASAs for {s}\n", .{options.input_name});
    try writer.print("REM  Atomic radii and reference values for relative SASA: {s}\n", .{options.classifier_name});
    try writer.print("REM  Algorithm: {s}\n", .{options.algorithm_name});
    try writer.print("REM  Probe-radius: {d:.2}\n", .{options.probe_radius});
    if (options.detail_count > 0) {
        try writer.print("REM  {s}: {d}\n", .{ options.detail_label, options.detail_count });
    }
    try writer.writeAll("REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar\n");
    try writer.writeAll("REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL\n");

    var total = AreaBreakdown{};
    for (residues) |residue| {
        var num_buf: [32]u8 = undefined;
        const num = residueNumberString(&num_buf, residue.residue_number, residue.insertion_code);
        const rel_total: ?f64 = if (analysis.MaxSASA.get(residue.residue.slice())) |max_sasa|
            if (max_sasa > 0) residue.area.total * 100.0 / max_sasa else null
        else
            null;

        try writer.print("RES {s:>3} {s:>3} {s:<4} ", .{ residue.residue.slice(), residue.chain.slice(), num });
        try writeAbsRel(writer, residue.area.total, rel_total);
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
        try writer.print("CHAIN{d:3} {s:>3} {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}\n", .{
            i + 1,
            chain.chain.slice(),
            chain.area.total,
            chain.area.side_chain,
            chain.area.main_chain,
            chain.area.apolar,
            chain.area.polar,
        });
    }

    try writer.writeAll("END  Absolute sums over all chains\n");
    try writer.print("TOTAL        {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}   {d:10.1}\n", .{
        total.total,
        total.side_chain,
        total.main_chain,
        total.apolar,
        total.polar,
    });

    return aw.toOwnedSlice();
}

/// Convert SasaResult to rich CSV string with structural information
/// Format: chain,residue,resnum,atom_name,x,y,z,radius,area
/// Caller must free the returned slice
pub fn sasaResultToRichCsv(allocator: Allocator, input: AtomInput, atom_areas: []const f64) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    errdefer aw.deinit();
    const writer = &aw.writer;

    // Header
    try writer.writeAll("chain,residue,resnum,atom_name,x,y,z,radius,area\n");

    // Atom rows
    const n = input.atomCount();
    for (0..n) |i| {
        // Chain
        if (input.chain_id) |chains| {
            try writer.writeAll(chains[i].slice());
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Residue name
        if (input.residue) |residues| {
            try writer.writeAll(residues[i].slice());
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

        // Atom name
        if (input.atom_name) |names| {
            try writer.writeAll(names[i].slice());
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

    // Total row
    var total: f64 = 0;
    for (atom_areas) |a| total += a;
    try writer.print(",,,,,,,,{d:.6}\n", .{total});

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
        self.allocator.free(self.residue_name);
        self.allocator.free(self.residue_number);
        self.allocator.free(self.residue_insertion_code);
        self.allocator.free(self.residue_atom_start);
        self.allocator.free(self.residue_atom_count);
        self.allocator.free(self.residue_sasa);
        self.* = undefined;
    }
};

fn sameResidue(
    chain_ids: []const types.FixedString4,
    residue_names: []const types.FixedString5,
    residue_nums: []const i32,
    insertion_codes: []const types.FixedString4,
    a: usize,
    b: usize,
) bool {
    return residue_nums[a] == residue_nums[b] and
        std.mem.eql(u8, chain_ids[a].slice(), chain_ids[b].slice()) and
        std.mem.eql(u8, residue_names[a].slice(), residue_names[b].slice()) and
        std.mem.eql(u8, insertion_codes[a].slice(), insertion_codes[b].slice());
}

pub fn buildResidueMap(allocator: Allocator, input: AtomInput, atom_areas: []const f64) !ResidueMap {
    const n = input.atomCount();
    if (atom_areas.len != n) return error.LengthMismatch;

    const chain_ids = input.chain_id orelse return error.MissingChainInfo;
    const residue_names = input.residue orelse return error.MissingResidueInfo;
    const residue_nums = input.residue_num orelse return error.MissingResidueNumInfo;
    const insertion_codes = input.insertion_code orelse return error.MissingInsertionCodeInfo;

    var residue_count: usize = 0;
    var i: usize = 0;
    while (i < n) {
        const start = i;
        residue_count += 1;
        i += 1;
        while (i < n and sameResidue(chain_ids, residue_names, residue_nums, insertion_codes, start, i)) : (i += 1) {}
    }

    const residue_chain = try allocator.alloc(types.FixedString4, residue_count);
    errdefer allocator.free(residue_chain);
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

    i = 0;
    var residue_idx: usize = 0;
    while (i < n) : (residue_idx += 1) {
        const start = i;
        var sasa = atom_areas[i];
        i += 1;
        while (i < n and sameResidue(chain_ids, residue_names, residue_nums, insertion_codes, start, i)) : (i += 1) {
            sasa += atom_areas[i];
        }

        residue_chain[residue_idx] = chain_ids[start];
        residue_name[residue_idx] = residue_names[start];
        residue_number[residue_idx] = residue_nums[start];
        residue_insertion_code[residue_idx] = insertion_codes[start];
        residue_atom_start[residue_idx] = start;
        residue_atom_count[residue_idx] = i - start;
        residue_sasa[residue_idx] = sasa;
    }

    return .{
        .allocator = allocator,
        .residue_chain = residue_chain,
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
    const JsonlEntry = struct {
        filename: []const u8,
        total_area: f64,
        atom_areas: []const f64,
    };

    const entry = JsonlEntry{
        .filename = filename,
        .total_area = total_area,
        .atom_areas = atom_areas,
    };

    return std.json.Stringify.valueAlloc(allocator, entry, .{});
}

pub fn fileResultWithResidueMapToJsonlLine(
    allocator: Allocator,
    filename: []const u8,
    total_area: f64,
    atom_areas: []const f64,
    residue_map: ResidueMap,
) ![]u8 {
    const residue_chain = try allocator.alloc([]const u8, residue_map.len());
    defer allocator.free(residue_chain);
    const residue_name = try allocator.alloc([]const u8, residue_map.len());
    defer allocator.free(residue_name);
    const residue_insertion_code = try allocator.alloc([]const u8, residue_map.len());
    defer allocator.free(residue_insertion_code);

    for (0..residue_map.len()) |i| {
        residue_chain[i] = residue_map.residue_chain[i].slice();
        residue_name[i] = residue_map.residue_name[i].slice();
        residue_insertion_code[i] = residue_map.residue_insertion_code[i].slice();
    }

    const JsonlEntry = struct {
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

    const entry = JsonlEntry{
        .filename = filename,
        .total_area = total_area,
        .atom_areas = atom_areas,
        .residue_chain = residue_chain,
        .residue_name = residue_name,
        .residue_number = residue_map.residue_number,
        .residue_insertion_code = residue_insertion_code,
        .residue_atom_start = residue_map.residue_atom_start,
        .residue_atom_count = residue_map.residue_atom_count,
        .residue_sasa = residue_map.residue_sasa,
    };

    return std.json.Stringify.valueAlloc(allocator, entry, .{});
}

pub const BsaAnalysisJsonl = struct {
    filename: []const u8,
    name: []const u8,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
    sasa_partner_a: f64,
    sasa_partner_b: f64,
    sasa_complex: f64,
    delta_sasa_total: f64,
    bsa: f64,
    delta_sasa_level: []const u8,
    residue_chain: []const []const u8 = &.{},
    residue_name: []const []const u8 = &.{},
    residue_number: []const i32 = &.{},
    residue_insertion_code: []const []const u8 = &.{},
    residue_delta_sasa: []const f64 = &.{},
};

pub fn bsaAnalysisToJsonlLine(allocator: Allocator, row: BsaAnalysisJsonl) ![]u8 {
    if (std.mem.eql(u8, row.delta_sasa_level, "residue")) {
        const Entry = struct {
            filename: []const u8,
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
            residue_chain: []const []const u8,
            residue_name: []const []const u8,
            residue_number: []const i32,
            residue_insertion_code: []const []const u8,
            residue_delta_sasa: []const f64,
        };
        return std.json.Stringify.valueAlloc(allocator, Entry{
            .filename = row.filename,
            .analysis = "bsa",
            .name = row.name,
            .partner_a = row.partner_a,
            .partner_b = row.partner_b,
            .sasa_partner_a = row.sasa_partner_a,
            .sasa_partner_b = row.sasa_partner_b,
            .sasa_complex = row.sasa_complex,
            .delta_sasa_total = row.delta_sasa_total,
            .bsa = row.bsa,
            .delta_sasa_level = row.delta_sasa_level,
            .residue_chain = row.residue_chain,
            .residue_name = row.residue_name,
            .residue_number = row.residue_number,
            .residue_insertion_code = row.residue_insertion_code,
            .residue_delta_sasa = row.residue_delta_sasa,
        }, .{});
    }

    const Entry = struct {
        filename: []const u8,
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
        .filename = row.filename,
        .analysis = "bsa",
        .name = row.name,
        .partner_a = row.partner_a,
        .partner_b = row.partner_b,
        .sasa_partner_a = row.sasa_partner_a,
        .sasa_partner_b = row.sasa_partner_b,
        .sasa_complex = row.sasa_complex,
        .delta_sasa_total = row.delta_sasa_total,
        .bsa = row.bsa,
        .delta_sasa_level = row.delta_sasa_level,
    }, .{});
}

// Tests
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
        "{\"filename\":\"example.cif\",\"total_area\":13.75,\"atom_areas\":[10,2.5,1.25],\"residue_chain\":[\"A\",\"B\"],\"residue_name\":[\"MET\",\"ALA\"],\"residue_number\":[1,7],\"residue_insertion_code\":[\"\",\"\"],\"residue_atom_start\":[0,2],\"residue_atom_count\":[2,1],\"residue_sasa\":[12.5,1.25]}",
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
    const residue_delta_sasa = [_]f64{ 3.0, 5.0 };

    const line = try bsaAnalysisToJsonlLine(allocator, .{
        .filename = "tiny.pdb",
        .name = "interface_ab",
        .partner_a = partner_a[0..],
        .partner_b = partner_b[0..],
        .sasa_partner_a = 10.0,
        .sasa_partner_b = 20.0,
        .sasa_complex = 14.0,
        .delta_sasa_total = 16.0,
        .bsa = 8.0,
        .delta_sasa_level = "residue",
        .residue_chain = residue_chain[0..],
        .residue_name = residue_name[0..],
        .residue_number = residue_number[0..],
        .residue_insertion_code = residue_insertion_code[0..],
        .residue_delta_sasa = residue_delta_sasa[0..],
    });
    defer allocator.free(line);

    try std.testing.expect(std.mem.indexOf(u8, line, "\"analysis\":\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"delta_sasa_total\":16") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"bsa\":8") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_delta_sasa\":[3,5]") != null);
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
        \\REM  Atomic radii and reference values for relative SASA: naccess
        \\REM  Algorithm: Shrake & Rupley
        \\REM  Probe-radius: 1.40
        \\REM  Test-points: 100
        \\REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
        \\REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
        \\RES ALA   A 1      30.00  23.3  20.00   N/A  10.00   N/A  20.00   N/A  10.00   N/A
        \\RES UNK   B 2A     30.00   N/A  30.00   N/A   0.00   N/A  30.00   N/A   0.00   N/A
        \\END  Absolute sums over single chains surface
        \\CHAIN  1   A       30.0         20.0         10.0         20.0         10.0
        \\CHAIN  2   B       30.0         30.0          0.0         30.0          0.0
        \\END  Absolute sums over all chains
        \\TOTAL              60.0         50.0         10.0         50.0         10.0
        \\
    , output);
}

test "sasaResultToRsa aggregates non-contiguous atoms from the same residue" {
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

    try std.testing.expect(std.mem.indexOf(u8, output, "RES ALA   A 1      30.00  23.3") != null);
    try std.testing.expectEqual(@as(?usize, null), std.mem.indexOf(u8, output, "RES ALA   A 1      10.00"));
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

    const test_path = "test_output.json";
    defer std.Io.Dir.cwd().deleteFile(io, test_path) catch {};

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

    const test_path = "test_format_json.json";
    defer std.Io.Dir.cwd().deleteFile(io, test_path) catch {};

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

    const test_path = "test_format.csv";
    defer std.Io.Dir.cwd().deleteFile(io, test_path) catch {};

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

    const test_path = "test_overwrite.json";
    defer std.Io.Dir.cwd().deleteFile(io, test_path) catch {};

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

    // Check header
    try std.testing.expect(std.mem.startsWith(u8, csv, "chain,residue,resnum,atom_name,x,y,z,radius,area\n"));

    // Check it contains expected data
    try std.testing.expect(std.mem.find(u8, csv, "A,ALA,1,N,") != null);
    try std.testing.expect(std.mem.find(u8, csv, "A,ALA,1,CA,") != null);

    // Check total row exists
    try std.testing.expect(std.mem.find(u8, csv, ",,,,,,,,30.800000\n") != null);
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
    try std.testing.expect(std.mem.find(u8, csv, "-,-,-,-,1.000,2.000,3.000,1.500,15.000000\n") != null);
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
    try std.testing.expectEqualStrings("test.pdb", obj.get("filename").?.string);
}

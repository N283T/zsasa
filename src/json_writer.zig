const std = @import("std");
const types = @import("types.zig");

const Allocator = std.mem.Allocator;
const SasaResult = types.SasaResult;
const AtomInput = types.AtomInput;

/// Output format options
pub const OutputFormat = enum {
    json, // Pretty-printed JSON (default)
    compact, // Single-line JSON
    csv, // CSV format
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
    var list = std.ArrayListUnmanaged(u8){};
    errdefer list.deinit(allocator);

    const writer = list.writer(allocator);

    // Header
    try writer.writeAll("atom_index,area\n");

    // Atom areas
    for (result.atom_areas, 0..) |area, i| {
        try writer.print("{d},{d:.6}\n", .{ i, area });
    }

    // Total
    try writer.print("total,{d:.6}\n", .{result.total_area});

    return list.toOwnedSlice(allocator);
}

/// Convert SasaResult to rich CSV string with structural information
/// Format: chain,residue,resnum,atom_name,x,y,z,radius,area
/// Caller must free the returned slice
pub fn sasaResultToRichCsv(allocator: Allocator, input: AtomInput, atom_areas: []const f64) ![]u8 {
    var list = std.ArrayListUnmanaged(u8){};
    errdefer list.deinit(allocator);

    const writer = list.writer(allocator);

    // Header
    try writer.writeAll("chain,residue,resnum,atom_name,x,y,z,radius,area\n");

    // Atom rows
    const n = input.atomCount();
    for (0..n) |i| {
        // Chain
        if (input.chain_id) |chains| {
            try writer.writeAll(chains[i]);
        } else {
            try writer.writeAll("-");
        }
        try writer.writeAll(",");

        // Residue name
        if (input.residue) |residues| {
            try writer.writeAll(residues[i]);
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
            try writer.writeAll(names[i]);
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

    return list.toOwnedSlice(allocator);
}

/// Write SasaResult to file with specified format
pub fn writeSasaResultWithFormat(
    allocator: Allocator,
    result: SasaResult,
    path: []const u8,
    format: OutputFormat,
) !void {
    const output_str = switch (format) {
        .json => try sasaResultToJsonPretty(allocator, result),
        .compact => try sasaResultToJson(allocator, result),
        .csv => try sasaResultToCsv(allocator, result),
    };
    defer allocator.free(output_str);

    const file = try std.fs.cwd().createFile(path, .{});
    defer file.close();

    try file.writeAll(output_str);
}

/// Write SasaResult to file with specified format, using rich CSV when input has structural info
pub fn writeSasaResultWithFormatAndInput(
    allocator: Allocator,
    result: SasaResult,
    input: AtomInput,
    path: []const u8,
    format: OutputFormat,
) !void {
    const output_str = switch (format) {
        .json => try sasaResultToJsonPretty(allocator, result),
        .compact => try sasaResultToJson(allocator, result),
        .csv => if (input.hasResidueInfo())
            try sasaResultToRichCsv(allocator, input, result.atom_areas)
        else
            try sasaResultToCsv(allocator, result),
    };
    defer allocator.free(output_str);

    const file = try std.fs.cwd().createFile(path, .{});
    defer file.close();

    try file.writeAll(output_str);
}

/// Write SasaResult to JSON file (default: compact for backward compatibility)
pub fn writeSasaResult(allocator: Allocator, result: SasaResult, path: []const u8) !void {
    return writeSasaResultWithFormat(allocator, result, path, .compact);
}

// Tests
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
    defer std.fs.cwd().deleteFile(test_path) catch {};

    try writeSasaResult(allocator, result, test_path);

    // Read back and verify
    const file = try std.fs.cwd().openFile(test_path, .{});
    defer file.close();

    const content = try file.readToEndAlloc(allocator, 1024);
    defer allocator.free(content);

    try std.testing.expectEqualStrings(
        "{\"total_area\":30.8,\"atom_areas\":[10.5,20.3]}",
        content,
    );
}

test "writeSasaResultWithFormat json" {
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

    const test_path = "test_format_json.json";
    defer std.fs.cwd().deleteFile(test_path) catch {};

    try writeSasaResultWithFormat(allocator, result, test_path, .json);

    const file = try std.fs.cwd().openFile(test_path, .{});
    defer file.close();

    const content = try file.readToEndAlloc(allocator, 1024);
    defer allocator.free(content);

    // Should be pretty-printed
    try std.testing.expect(std.mem.indexOf(u8, content, "\n") != null);
}

test "writeSasaResultWithFormat csv" {
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

    const test_path = "test_format.csv";
    defer std.fs.cwd().deleteFile(test_path) catch {};

    try writeSasaResultWithFormat(allocator, result, test_path, .csv);

    const file = try std.fs.cwd().openFile(test_path, .{});
    defer file.close();

    const content = try file.readToEndAlloc(allocator, 1024);
    defer allocator.free(content);

    // Should start with header
    try std.testing.expect(std.mem.startsWith(u8, content, "atom_index,area\n"));
}

test "writeSasaResult overwrites existing file" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 1);
    defer allocator.free(atom_areas);

    atom_areas[0] = 50.0;

    const result = SasaResult{
        .total_area = 50.0,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };

    const test_path = "test_overwrite.json";
    defer std.fs.cwd().deleteFile(test_path) catch {};

    // Write first time
    try writeSasaResult(allocator, result, test_path);

    // Write second time (overwrite)
    atom_areas[0] = 99.9;
    const result2 = SasaResult{
        .total_area = 99.9,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
    try writeSasaResult(allocator, result2, test_path);

    // Verify overwrite
    const file = try std.fs.cwd().openFile(test_path, .{});
    defer file.close();

    const content = try file.readToEndAlloc(allocator, 1024);
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
    const chain_ids = try allocator.alloc([]const u8, 2);
    defer allocator.free(chain_ids);
    chain_ids[0] = "A";
    chain_ids[1] = "A";

    const residues = try allocator.alloc([]const u8, 2);
    defer allocator.free(residues);
    residues[0] = "ALA";
    residues[1] = "ALA";

    const atom_names = try allocator.alloc([]const u8, 2);
    defer allocator.free(atom_names);
    atom_names[0] = "N";
    atom_names[1] = "CA";

    const residue_nums = try allocator.alloc(i32, 2);
    defer allocator.free(residue_nums);
    residue_nums[0] = 1;
    residue_nums[1] = 1;

    const insertion_codes = try allocator.alloc([]const u8, 2);
    defer allocator.free(insertion_codes);
    insertion_codes[0] = "";
    insertion_codes[1] = "";

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
    try std.testing.expect(std.mem.indexOf(u8, csv, "A,ALA,1,N,") != null);
    try std.testing.expect(std.mem.indexOf(u8, csv, "A,ALA,1,CA,") != null);

    // Check total row exists
    try std.testing.expect(std.mem.indexOf(u8, csv, ",,,,,,,,30.800000\n") != null);
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
    try std.testing.expect(std.mem.indexOf(u8, csv, "-,-,-,-,1.000,2.000,3.000,1.500,15.000000\n") != null);
}

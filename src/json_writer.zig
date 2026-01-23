const std = @import("std");
const types = @import("types.zig");

const Allocator = std.mem.Allocator;
const SasaResult = types.SasaResult;

/// JSON structure for output
const JsonOutput = struct {
    total_area: f64,
    atom_areas: []const f64,
};

/// Convert SasaResult to JSON string
/// Caller must free the returned slice
pub fn sasaResultToJson(allocator: Allocator, result: SasaResult) ![]u8 {
    const output = JsonOutput{
        .total_area = result.total_area,
        .atom_areas = result.atom_areas,
    };

    return std.json.Stringify.valueAlloc(allocator, output, .{});
}

/// Write SasaResult to JSON file
pub fn writeSasaResult(allocator: Allocator, result: SasaResult, path: []const u8) !void {
    const json_str = try sasaResultToJson(allocator, result);
    defer allocator.free(json_str);

    const file = try std.fs.cwd().createFile(path, .{});
    defer file.close();

    try file.writeAll(json_str);
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

const std = @import("std");
const types = @import("types.zig");
const AtomInput = types.AtomInput;
const Allocator = std.mem.Allocator;

/// JSON structure matching the input file format
const JsonInput = struct {
    x: []f64,
    y: []f64,
    z: []f64,
    r: []f64,
};

/// Parse atom input from JSON string
pub fn parseAtomInput(allocator: Allocator, json_str: []const u8) !AtomInput {
    const parsed = try std.json.parseFromSlice(
        JsonInput,
        allocator,
        json_str,
        .{},
    );
    defer parsed.deinit();

    const data = parsed.value;

    // Validate all arrays have same length
    const n = data.x.len;
    if (data.y.len != n or data.z.len != n or data.r.len != n) {
        return error.ArrayLengthMismatch;
    }

    if (n == 0) {
        return error.EmptyInput;
    }

    // Allocate and copy data
    const x = try allocator.alloc(f64, n);
    errdefer allocator.free(x);

    const y = try allocator.alloc(f64, n);
    errdefer allocator.free(y);

    const z = try allocator.alloc(f64, n);
    errdefer allocator.free(z);

    const r = try allocator.alloc(f64, n);
    errdefer allocator.free(r);

    @memcpy(x, data.x);
    @memcpy(y, data.y);
    @memcpy(z, data.z);
    @memcpy(r, data.r);

    return AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
}

/// Read atom input from JSON file
pub fn readAtomInputFromFile(allocator: Allocator, path: []const u8) !AtomInput {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();

    const max_size = 100 * 1024 * 1024; // 100 MB max
    const contents = try file.readToEndAlloc(allocator, max_size);
    defer allocator.free(contents);

    return try parseAtomInput(allocator, contents);
}

// Tests
test "parseAtomInput basic" {
    const allocator = std.testing.allocator;

    const json =
        \\{"x": [1.0, 2.0, 3.0], "y": [4.0, 5.0, 6.0], "z": [7.0, 8.0, 9.0], "r": [1.5, 1.6, 1.7]}
    ;

    var input = try parseAtomInput(allocator, json);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
    try std.testing.expectEqual(@as(f64, 1.0), input.x[0]);
    try std.testing.expectEqual(@as(f64, 2.0), input.x[1]);
    try std.testing.expectEqual(@as(f64, 3.0), input.x[2]);
    try std.testing.expectEqual(@as(f64, 4.0), input.y[0]);
    try std.testing.expectEqual(@as(f64, 1.5), input.r[0]);
    try std.testing.expectEqual(@as(f64, 1.7), input.r[2]);
}

test "parseAtomInput empty arrays" {
    const allocator = std.testing.allocator;

    const json =
        \\{"x": [], "y": [], "z": [], "r": []}
    ;

    const result = parseAtomInput(allocator, json);
    try std.testing.expectError(error.EmptyInput, result);
}

test "parseAtomInput mismatched lengths" {
    const allocator = std.testing.allocator;

    const json =
        \\{"x": [1.0, 2.0], "y": [4.0, 5.0], "z": [7.0, 8.0], "r": [1.5]}
    ;

    const result = parseAtomInput(allocator, json);
    try std.testing.expectError(error.ArrayLengthMismatch, result);
}

test "parseAtomInput missing field" {
    const allocator = std.testing.allocator;

    const json =
        \\{"x": [1.0, 2.0], "y": [4.0, 5.0], "z": [7.0, 8.0]}
    ;

    const result = parseAtomInput(allocator, json);
    try std.testing.expectError(error.MissingField, result);
}

test "parseAtomInput invalid JSON" {
    const allocator = std.testing.allocator;

    const json = "not valid json";

    const result = parseAtomInput(allocator, json);
    try std.testing.expect(std.meta.isError(result));
}

test "readAtomInputFromFile with real file" {
    const allocator = std.testing.allocator;

    const path = "examples/input_1a0q.json";

    // Check if file exists
    const file = std.fs.cwd().openFile(path, .{}) catch {
        // File doesn't exist, skip test
        return error.SkipZigTest;
    };
    file.close();

    var input = try readAtomInputFromFile(allocator, path);
    defer input.deinit();

    // Verify we got 3183 atoms as expected
    try std.testing.expectEqual(@as(usize, 3183), input.atomCount());

    // Verify some values match the file
    try std.testing.expectApproxEqAbs(@as(f64, 27.234), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 26.259), input.x[1], 0.001);

    // Verify all arrays are valid
    for (input.x, 0..) |_, i| {
        _ = input.x[i];
        _ = input.y[i];
        _ = input.z[i];
        _ = input.r[i];
    }
}

test "readAtomInputFromFile nonexistent file" {
    const allocator = std.testing.allocator;

    const result = readAtomInputFromFile(allocator, "nonexistent_file.json");
    try std.testing.expectError(error.FileNotFound, result);
}

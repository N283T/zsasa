const std = @import("std");
const types = @import("types.zig");
const AtomInput = types.AtomInput;
const Allocator = std.mem.Allocator;
const flate = std.compress.flate;

/// JSON structure matching the input file format
const JsonInput = struct {
    x: []f64,
    y: []f64,
    z: []f64,
    r: []f64,
    residue: ?[][]const u8 = null,
    atom_name: ?[][]const u8 = null,
    /// Element atomic numbers (e.g., 6=C, 7=N, 8=O, 20=Ca)
    element: ?[]u8 = null,
};

/// Validation error details
pub const ValidationError = struct {
    message: []const u8,
    index: ?usize = null, // Atom index if applicable
    value: ?f64 = null, // Invalid value if applicable
};

/// Validation result
pub const ValidationResult = struct {
    valid: bool,
    errors: []ValidationError,
    allocator: Allocator,

    pub fn deinit(self: *ValidationResult) void {
        self.allocator.free(self.errors);
    }
};

/// Maximum allowed radius in Angstroms
const MAX_RADIUS_ANGSTROMS: f64 = 100.0;

/// Validate a single value is finite (not NaN or Inf)
fn isFinite(value: f64) bool {
    return std.math.isFinite(value);
}

/// Validate radius value (must be positive and finite)
fn isValidRadius(value: f64) bool {
    return isFinite(value) and value > 0 and value <= MAX_RADIUS_ANGSTROMS;
}

/// Validate input data and collect all errors
pub fn validateInput(allocator: Allocator, input: AtomInput) !ValidationResult {
    var errors = std.ArrayListUnmanaged(ValidationError){};
    errdefer errors.deinit(allocator);

    const n = input.atomCount();

    // Check each atom
    for (0..n) |i| {
        // Check coordinates are finite
        if (!isFinite(input.x[i])) {
            try errors.append(allocator, ValidationError{
                .message = "Invalid x coordinate (NaN or Inf)",
                .index = i,
                .value = input.x[i],
            });
        }
        if (!isFinite(input.y[i])) {
            try errors.append(allocator, ValidationError{
                .message = "Invalid y coordinate (NaN or Inf)",
                .index = i,
                .value = input.y[i],
            });
        }
        if (!isFinite(input.z[i])) {
            try errors.append(allocator, ValidationError{
                .message = "Invalid z coordinate (NaN or Inf)",
                .index = i,
                .value = input.z[i],
            });
        }

        // Check radius is valid
        if (!isValidRadius(input.r[i])) {
            if (!isFinite(input.r[i])) {
                try errors.append(allocator, ValidationError{
                    .message = "Invalid radius (NaN or Inf)",
                    .index = i,
                    .value = input.r[i],
                });
            } else if (input.r[i] <= 0) {
                try errors.append(allocator, ValidationError{
                    .message = "Radius must be positive",
                    .index = i,
                    .value = input.r[i],
                });
            } else {
                try errors.append(allocator, ValidationError{
                    .message = "Radius too large (max 100 Angstroms)",
                    .index = i,
                    .value = input.r[i],
                });
            }
        }
    }

    return ValidationResult{
        .valid = errors.items.len == 0,
        .errors = try errors.toOwnedSlice(allocator),
        .allocator = allocator,
    };
}

/// Print validation errors to stderr
pub fn printValidationErrors(errors: []const ValidationError) void {
    std.debug.print("Input validation failed with {} error(s):\n", .{errors.len});
    for (errors, 0..) |err, i| {
        if (i >= 10) {
            std.debug.print("  ... and {} more errors\n", .{errors.len - 10});
            break;
        }
        if (err.index) |idx| {
            if (err.value) |val| {
                std.debug.print("  - Atom {}: {s} (value: {d})\n", .{ idx, err.message, val });
            } else {
                std.debug.print("  - Atom {}: {s}\n", .{ idx, err.message });
            }
        } else {
            std.debug.print("  - {s}\n", .{err.message});
        }
    }
}

/// Check for duplicate coordinates and print warning if found.
/// Returns the number of duplicate coordinate sets found.
/// This is a warning only - duplicate atoms can cause SASA calculation discrepancies
/// but are not treated as validation errors.
pub fn checkDuplicateCoordinates(allocator: Allocator, input: AtomInput) !usize {
    const n = input.atomCount();
    if (n < 2) return 0;

    // Use a hash map to detect duplicates
    // Key: packed coordinate bytes, Value: first occurrence index
    const CoordKey = struct {
        x_bits: u64,
        y_bits: u64,
        z_bits: u64,

        fn fromCoords(x: f64, y: f64, z: f64) @This() {
            return .{
                .x_bits = @bitCast(x),
                .y_bits = @bitCast(y),
                .z_bits = @bitCast(z),
            };
        }
    };

    var seen = std.AutoHashMap(CoordKey, usize).init(allocator);
    defer seen.deinit();

    var duplicate_count: usize = 0;

    for (0..n) |i| {
        const key = CoordKey.fromCoords(input.x[i], input.y[i], input.z[i]);
        const result = try seen.getOrPut(key);
        if (result.found_existing) {
            duplicate_count += 1;
        } else {
            result.value_ptr.* = i;
        }
    }

    if (duplicate_count > 0) {
        std.debug.print(
            "Warning: Found {} duplicate coordinate(s) in {} atoms. " ++
                "This may cause SASA calculation discrepancies with other tools.\n",
            .{ duplicate_count, n },
        );
    }

    return duplicate_count;
}

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

    // Validate optional arrays if present
    if (data.residue) |res| {
        if (res.len != n) return error.ArrayLengthMismatch;
    }
    if (data.atom_name) |names| {
        if (names.len != n) return error.ArrayLengthMismatch;
    }
    if (data.element) |elem| {
        if (elem.len != n) return error.ArrayLengthMismatch;
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

    // Copy optional residue names (using FixedString4)
    var residue: ?[]types.FixedString4 = null;
    if (data.residue) |res| {
        const res_copy = try allocator.alloc(types.FixedString4, n);
        for (res, 0..) |s, i| {
            res_copy[i] = types.FixedString4.fromSlice(s);
        }
        residue = res_copy;
    }
    errdefer if (residue) |res| allocator.free(res);

    // Copy optional atom names (using FixedString4)
    var atom_name: ?[]types.FixedString4 = null;
    if (data.atom_name) |names| {
        const names_copy = try allocator.alloc(types.FixedString4, n);
        for (names, 0..) |s, i| {
            names_copy[i] = types.FixedString4.fromSlice(s);
        }
        atom_name = names_copy;
    }

    // Copy optional element atomic numbers
    var element: ?[]u8 = null;
    if (data.element) |elem| {
        element = try allocator.dupe(u8, elem);
    }
    errdefer if (element) |elem| allocator.free(elem);

    return AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .residue = residue,
        .atom_name = atom_name,
        .element = element,
        .allocator = allocator,
    };
}

/// Read atom input from JSON file (handles both .json and .json.gz)
pub fn readAtomInputFromFile(allocator: Allocator, path: []const u8) !AtomInput {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();

    const max_size = 200 * 1024 * 1024; // 200 MB max

    // Check if file is gzip compressed
    if (std.mem.endsWith(u8, path, ".gz")) {
        // Read compressed file into memory first
        const compressed = try file.readToEndAlloc(allocator, max_size);
        defer allocator.free(compressed);

        // Decompress gzip using flate with fixed reader
        var input_reader: std.Io.Reader = .fixed(compressed);
        var decompress: flate.Decompress = .init(&input_reader, .gzip, &.{});

        // Read all decompressed data
        const contents = try decompress.reader.allocRemaining(allocator, std.Io.Limit.limited(max_size));
        defer allocator.free(contents);
        return try parseAtomInput(allocator, contents);
    } else {
        const contents = try file.readToEndAlloc(allocator, max_size);
        defer allocator.free(contents);
        return try parseAtomInput(allocator, contents);
    }
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
    try std.testing.expect(!input.hasClassificationInfo());
}

test "parseAtomInput with residue and atom_name" {
    const allocator = std.testing.allocator;

    const json =
        \\{"x": [1.0, 2.0], "y": [3.0, 4.0], "z": [5.0, 6.0], "r": [1.5, 1.6], "residue": ["ALA", "GLY"], "atom_name": ["CA", "N"]}
    ;

    var input = try parseAtomInput(allocator, json);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expect(input.hasClassificationInfo());

    const residue = input.residue.?;
    const atom_name = input.atom_name.?;

    try std.testing.expectEqualStrings("ALA", residue[0]);
    try std.testing.expectEqualStrings("GLY", residue[1]);
    try std.testing.expectEqualStrings("CA", atom_name[0]);
    try std.testing.expectEqualStrings("N", atom_name[1]);
}

test "parseAtomInput with element atomic numbers" {
    const allocator = std.testing.allocator;

    // 6=Carbon, 7=Nitrogen, 8=Oxygen
    const json =
        \\{"x": [1.0, 2.0, 3.0], "y": [4.0, 5.0, 6.0], "z": [7.0, 8.0, 9.0], "r": [1.5, 1.6, 1.7], "residue": ["ALA", "ALA", "ALA"], "atom_name": ["CA", "N", "O"], "element": [6, 7, 8]}
    ;

    var input = try parseAtomInput(allocator, json);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
    try std.testing.expect(input.hasClassificationInfo());
    try std.testing.expect(input.hasElementInfo());

    const element = input.element.?;
    try std.testing.expectEqual(@as(u8, 6), element[0]); // Carbon
    try std.testing.expectEqual(@as(u8, 7), element[1]); // Nitrogen
    try std.testing.expectEqual(@as(u8, 8), element[2]); // Oxygen
}

test "parseAtomInput with element distinguishes CA (Carbon) from Ca (Calcium)" {
    const allocator = std.testing.allocator;

    // First atom: CA = Carbon alpha (atomic number 6)
    // Second atom: CA = Calcium ion (atomic number 20)
    const json =
        \\{"x": [1.0, 2.0], "y": [3.0, 4.0], "z": [5.0, 6.0], "r": [1.7, 2.31], "residue": ["ALA", "CA"], "atom_name": ["CA", "CA"], "element": [6, 20]}
    ;

    var input = try parseAtomInput(allocator, json);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expect(input.hasElementInfo());

    const element = input.element.?;
    try std.testing.expectEqual(@as(u8, 6), element[0]); // Carbon (CA in amino acid)
    try std.testing.expectEqual(@as(u8, 20), element[1]); // Calcium (CA ion)
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

test "validateInput valid data" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 3);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 3);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 3);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 3);
    defer allocator.free(r);

    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    y[0] = 1.0;
    y[1] = 2.0;
    y[2] = 3.0;
    z[0] = 1.0;
    z[1] = 2.0;
    z[2] = 3.0;
    r[0] = 1.5;
    r[1] = 1.6;
    r[2] = 1.7;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    var result = try validateInput(allocator, input);
    defer result.deinit();

    try std.testing.expect(result.valid);
    try std.testing.expectEqual(@as(usize, 0), result.errors.len);
}

test "validateInput negative radius" {
    const allocator = std.testing.allocator;

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
    y[0] = 1.0;
    y[1] = 2.0;
    z[0] = 1.0;
    z[1] = 2.0;
    r[0] = 1.5;
    r[1] = -1.0; // Invalid

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    var result = try validateInput(allocator, input);
    defer result.deinit();

    try std.testing.expect(!result.valid);
    try std.testing.expectEqual(@as(usize, 1), result.errors.len);
    try std.testing.expectEqual(@as(usize, 1), result.errors[0].index.?);
}

test "validateInput NaN coordinate" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 2);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 2);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 2);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 2);
    defer allocator.free(r);

    x[0] = std.math.nan(f64); // Invalid
    x[1] = 2.0;
    y[0] = 1.0;
    y[1] = 2.0;
    z[0] = 1.0;
    z[1] = 2.0;
    r[0] = 1.5;
    r[1] = 1.6;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    var result = try validateInput(allocator, input);
    defer result.deinit();

    try std.testing.expect(!result.valid);
    try std.testing.expectEqual(@as(usize, 1), result.errors.len);
    try std.testing.expectEqual(@as(usize, 0), result.errors[0].index.?);
}

test "validateInput infinity radius" {
    const allocator = std.testing.allocator;

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
    y[0] = 1.0;
    y[1] = 2.0;
    z[0] = 1.0;
    z[1] = 2.0;
    r[0] = std.math.inf(f64); // Invalid
    r[1] = 1.6;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    var result = try validateInput(allocator, input);
    defer result.deinit();

    try std.testing.expect(!result.valid);
    try std.testing.expectEqual(@as(usize, 1), result.errors.len);
}

test "validateInput zero radius" {
    const allocator = std.testing.allocator;

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
    y[0] = 1.0;
    y[1] = 2.0;
    z[0] = 1.0;
    z[1] = 2.0;
    r[0] = 0.0; // Invalid
    r[1] = 1.6;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    var result = try validateInput(allocator, input);
    defer result.deinit();

    try std.testing.expect(!result.valid);
    try std.testing.expectEqual(@as(usize, 1), result.errors.len);
}

test "validateInput radius too large" {
    const allocator = std.testing.allocator;

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
    y[0] = 1.0;
    y[1] = 2.0;
    z[0] = 1.0;
    z[1] = 2.0;
    r[0] = 150.0; // Exceeds MAX_RADIUS_ANGSTROMS (100.0)
    r[1] = 1.6;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    var result = try validateInput(allocator, input);
    defer result.deinit();

    try std.testing.expect(!result.valid);
    try std.testing.expectEqual(@as(usize, 1), result.errors.len);
    try std.testing.expectEqual(@as(usize, 0), result.errors[0].index.?);
}

test "checkDuplicateCoordinates no duplicates" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 3);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 3);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 3);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 3);
    defer allocator.free(r);

    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 3.0;
    y[0] = 1.0;
    y[1] = 2.0;
    y[2] = 3.0;
    z[0] = 1.0;
    z[1] = 2.0;
    z[2] = 3.0;
    r[0] = 1.5;
    r[1] = 1.5;
    r[2] = 1.5;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const count = try checkDuplicateCoordinates(allocator, input);
    try std.testing.expectEqual(@as(usize, 0), count);
}

test "checkDuplicateCoordinates with duplicates" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 4);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 4);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 4);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 4);
    defer allocator.free(r);

    // Atoms 0 and 2 have identical coordinates
    // Atoms 1 and 3 have identical coordinates
    x[0] = 1.0;
    x[1] = 2.0;
    x[2] = 1.0; // duplicate of 0
    x[3] = 2.0; // duplicate of 1
    y[0] = 1.0;
    y[1] = 2.0;
    y[2] = 1.0;
    y[3] = 2.0;
    z[0] = 1.0;
    z[1] = 2.0;
    z[2] = 1.0;
    z[3] = 2.0;
    r[0] = 1.5;
    r[1] = 1.5;
    r[2] = 1.5;
    r[3] = 1.5;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const count = try checkDuplicateCoordinates(allocator, input);
    try std.testing.expectEqual(@as(usize, 2), count);
}

test "checkDuplicateCoordinates single atom" {
    const allocator = std.testing.allocator;

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

    const count = try checkDuplicateCoordinates(allocator, input);
    try std.testing.expectEqual(@as(usize, 0), count);
}

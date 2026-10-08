const std = @import("std");
const types = @import("types.zig");
const compressed = @import("compressed.zig");
const AtomInput = types.AtomInput;
const Allocator = std.mem.Allocator;

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
pub fn isValidRadius(value: f64) bool {
    return isFinite(value) and value > 0 and value <= MAX_RADIUS_ANGSTROMS;
}

/// Validate input data and collect all errors
pub fn validateInput(allocator: Allocator, input: AtomInput) !ValidationResult {
    var errors = std.ArrayListUnmanaged(ValidationError).empty;
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

const DuplicateCoordinateOptions = struct {
    log_warning: bool = true,
};

/// Check for duplicate coordinates and print warning if found.
/// Returns the number of duplicate coordinate sets found.
/// This is a warning only - duplicate atoms can cause SASA calculation discrepancies
/// but are not treated as validation errors.
pub fn checkDuplicateCoordinates(allocator: Allocator, input: AtomInput) !usize {
    return checkDuplicateCoordinatesOptions(allocator, input, .{});
}

fn checkDuplicateCoordinatesOptions(allocator: Allocator, input: AtomInput, options: DuplicateCoordinateOptions) !usize {
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

    var seen: std.AutoHashMapUnmanaged(CoordKey, usize) = .empty;
    defer seen.deinit(allocator);

    var duplicate_count: usize = 0;

    for (0..n) |i| {
        const key = CoordKey.fromCoords(input.x[i], input.y[i], input.z[i]);
        const result = try seen.getOrPut(allocator, key);
        if (result.found_existing) {
            duplicate_count += 1;
        } else {
            result.value_ptr.* = i;
        }
    }

    if (duplicate_count > 0 and options.log_warning) {
        std.debug.print(
            "Warning: Found {} duplicate coordinate(s) in {} atoms. " ++
                "This may cause SASA calculation discrepancies with other tools.\n",
            .{ duplicate_count, n },
        );
    }

    return duplicate_count;
}

/// Field values of the JSON input while it is being read. Owns every list
/// until `parseAtomInput` moves them into the result.
const JsonFields = struct {
    x: ?std.ArrayList(f64) = null,
    y: ?std.ArrayList(f64) = null,
    z: ?std.ArrayList(f64) = null,
    r: ?std.ArrayList(f64) = null,
    residue: ?std.ArrayList(types.FixedString5) = null,
    atom_name: ?std.ArrayList(types.FixedString4) = null,
    element: ?std.ArrayList(u8) = null,

    fn deinit(self: *JsonFields, allocator: Allocator) void {
        inline for (std.meta.fields(JsonFields)) |f| {
            if (@field(self, f.name)) |*list| list.deinit(allocator);
        }
    }
};

fn freeToken(allocator: Allocator, token: std.json.Token) void {
    switch (token) {
        .allocated_number, .allocated_string => |s| allocator.free(s),
        else => {},
    }
}

/// Read the next token of `scanner`; the caller releases it with `freeToken`.
fn nextToken(allocator: Allocator, scanner: *std.json.Scanner) !std.json.Token {
    return scanner.nextAlloc(allocator, .alloc_if_needed);
}

/// Read a JSON number. A string such as `"3"` is not a number.
fn readNumber(allocator: Allocator, scanner: *std.json.Scanner) !f64 {
    const token = try nextToken(allocator, scanner);
    defer freeToken(allocator, token);
    return switch (token) {
        .number, .allocated_number => |s| std.fmt.parseFloat(f64, s) catch error.ExpectedNumber,
        else => error.ExpectedNumber,
    };
}

/// Read an atomic number: a JSON number that is a whole number from 0 to 255.
/// Strings such as `"CN"` (whose bytes used to be taken as atomic numbers) are not accepted.
fn readAtomicNumber(allocator: Allocator, scanner: *std.json.Scanner) !u8 {
    const token = try nextToken(allocator, scanner);
    defer freeToken(allocator, token);
    const text = switch (token) {
        .number, .allocated_number => |s| s,
        else => return error.InvalidElement,
    };
    if (std.fmt.parseInt(u8, text, 10)) |n| return n else |_| {}
    const f = std.fmt.parseFloat(f64, text) catch return error.InvalidElement;
    if (!(f >= 0 and f <= 255) or f != @floor(f)) return error.InvalidElement;
    return @intFromFloat(f);
}

fn readString(allocator: Allocator, scanner: *std.json.Scanner, comptime T: type) !T {
    const token = try nextToken(allocator, scanner);
    defer freeToken(allocator, token);
    return switch (token) {
        .string, .allocated_string => |s| T.fromSlice(s),
        else => error.ExpectedString,
    };
}

/// Read a JSON array of items into a new list. `null` is accepted only when
/// `nullable` is set (optional fields), and yields `null`.
fn readArray(
    comptime T: type,
    comptime readItem: anytype,
    allocator: Allocator,
    scanner: *std.json.Scanner,
    comptime nullable: bool,
) !?std.ArrayList(T) {
    if (nullable and try scanner.peekNextTokenType() == .null) {
        _ = try scanner.next();
        return null;
    }
    const begin = try nextToken(allocator, scanner);
    defer freeToken(allocator, begin);
    if (begin != .array_begin) return error.ExpectedArray;

    var list: std.ArrayList(T) = .empty;
    errdefer list.deinit(allocator);
    while (try scanner.peekNextTokenType() != .array_end) {
        try list.append(allocator, try readItem(allocator, scanner));
    }
    _ = try scanner.next(); // array_end
    return list;
}

fn readResidue(allocator: Allocator, scanner: *std.json.Scanner) !types.FixedString5 {
    return readString(allocator, scanner, types.FixedString5);
}

fn readAtomName(allocator: Allocator, scanner: *std.json.Scanner) !types.FixedString4 {
    return readString(allocator, scanner, types.FixedString4);
}

/// Parse atom input from a JSON string.
///
/// `x`, `y`, `z` and `r` must be arrays of JSON numbers, `residue` and
/// `atom_name` arrays of strings and `element` an array of atomic numbers
/// (numbers from 0 to 255). Values of another JSON type are rejected with
/// `ExpectedNumber`, `ExpectedString`, `ExpectedArray` or `InvalidElement`
/// instead of being converted.
pub fn parseAtomInput(allocator: Allocator, json_str: []const u8) !AtomInput {
    var scanner = std.json.Scanner.initCompleteInput(allocator, json_str);
    defer scanner.deinit();

    var fields = JsonFields{};
    defer fields.deinit(allocator);

    const open = try nextToken(allocator, &scanner);
    defer freeToken(allocator, open);
    if (open != .object_begin) return error.UnexpectedToken;

    while (true) {
        const key_token = try nextToken(allocator, &scanner);
        defer freeToken(allocator, key_token);
        const key = switch (key_token) {
            .string, .allocated_string => |s| s,
            .object_end => break,
            else => return error.UnexpectedToken,
        };

        if (std.mem.eql(u8, key, "x")) {
            if (fields.x != null) return error.DuplicateField;
            fields.x = (try readArray(f64, readNumber, allocator, &scanner, false)).?;
        } else if (std.mem.eql(u8, key, "y")) {
            if (fields.y != null) return error.DuplicateField;
            fields.y = (try readArray(f64, readNumber, allocator, &scanner, false)).?;
        } else if (std.mem.eql(u8, key, "z")) {
            if (fields.z != null) return error.DuplicateField;
            fields.z = (try readArray(f64, readNumber, allocator, &scanner, false)).?;
        } else if (std.mem.eql(u8, key, "r")) {
            if (fields.r != null) return error.DuplicateField;
            fields.r = (try readArray(f64, readNumber, allocator, &scanner, false)).?;
        } else if (std.mem.eql(u8, key, "residue")) {
            if (fields.residue != null) return error.DuplicateField;
            fields.residue = try readArray(types.FixedString5, readResidue, allocator, &scanner, true);
        } else if (std.mem.eql(u8, key, "atom_name")) {
            if (fields.atom_name != null) return error.DuplicateField;
            fields.atom_name = try readArray(types.FixedString4, readAtomName, allocator, &scanner, true);
        } else if (std.mem.eql(u8, key, "element")) {
            if (fields.element != null) return error.DuplicateField;
            fields.element = try readArray(u8, readAtomicNumber, allocator, &scanner, true);
        } else {
            return error.UnknownField;
        }
    }

    const end = try nextToken(allocator, &scanner);
    defer freeToken(allocator, end);
    if (end != .end_of_document) return error.UnexpectedToken;

    const x = fields.x orelse return error.MissingField;
    const y = fields.y orelse return error.MissingField;
    const z = fields.z orelse return error.MissingField;
    const r = fields.r orelse return error.MissingField;

    // Validate all arrays have same length
    const n = x.items.len;
    if (y.items.len != n or z.items.len != n or r.items.len != n) {
        return error.ArrayLengthMismatch;
    }
    if (n == 0) {
        return error.EmptyInput;
    }
    if (fields.residue) |res| {
        if (res.items.len != n) return error.ArrayLengthMismatch;
    }
    if (fields.atom_name) |names| {
        if (names.items.len != n) return error.ArrayLengthMismatch;
    }
    if (fields.element) |elem| {
        if (elem.items.len != n) return error.ArrayLengthMismatch;
    }

    // Move the lists into the result. A list that has been taken is empty, so
    // `fields.deinit` stays correct if a later step fails.
    const x_out = try fields.x.?.toOwnedSlice(allocator);
    errdefer allocator.free(x_out);
    const y_out = try fields.y.?.toOwnedSlice(allocator);
    errdefer allocator.free(y_out);
    const z_out = try fields.z.?.toOwnedSlice(allocator);
    errdefer allocator.free(z_out);
    const r_out = try fields.r.?.toOwnedSlice(allocator);
    errdefer allocator.free(r_out);

    var residue: ?[]types.FixedString5 = null;
    errdefer if (residue) |s| allocator.free(s);
    if (fields.residue != null) residue = try fields.residue.?.toOwnedSlice(allocator);

    var atom_name: ?[]types.FixedString4 = null;
    errdefer if (atom_name) |s| allocator.free(s);
    if (fields.atom_name != null) atom_name = try fields.atom_name.?.toOwnedSlice(allocator);

    var element: ?[]u8 = null;
    errdefer if (element) |s| allocator.free(s);
    if (fields.element != null) element = try fields.element.?.toOwnedSlice(allocator);

    return AtomInput{
        .x = x_out,
        .y = y_out,
        .z = z_out,
        .r = r_out,
        .residue = residue,
        .atom_name = atom_name,
        .element = element,
        .allocator = allocator,
    };
}

/// Read atom input from JSON file (handles plain, .gz, and .zst files)
pub fn readAtomInputFromFile(allocator: Allocator, io: std.Io, path: []const u8) !AtomInput {
    const max_size = 200 * 1024 * 1024; // 200 MB max

    const contents = if (compressed.isCompressed(path))
        try compressed.read(allocator, path)
    else blk: {
        const file = try std.Io.Dir.cwd().openFile(io, path, .{});
        defer file.close(io);
        var read_buf: [65536]u8 = undefined;
        var r = file.reader(io, &read_buf);
        break :blk try r.interface.allocRemaining(allocator, .limited64(max_size));
    };
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

    try std.testing.expectEqualStrings("ALA", residue[0].slice());
    try std.testing.expectEqualStrings("GLY", residue[1].slice());
    try std.testing.expectEqualStrings("CA", atom_name[0].slice());
    try std.testing.expectEqualStrings("N", atom_name[1].slice());
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
    const io = std.testing.io;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(io, .{
        .sub_path = "atoms.json",
        .data =
        \\{"x": [27.234, 26.259, 25.5],
        \\ "y": [14.262, 13.883, 12.7],
        \\ "z": [5.595, 6.312, 7.0],
        \\ "r": [1.7, 1.55, 1.52]}
        ,
    });

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(io, &root_buf);
    const path = try std.fs.path.join(allocator, &.{ root_buf[0..root_len], "atoms.json" });
    defer allocator.free(path);

    var input = try readAtomInputFromFile(allocator, io, path);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
    try std.testing.expectEqualSlices(f64, &.{ 27.234, 26.259, 25.5 }, input.x);
    try std.testing.expectEqualSlices(f64, &.{ 14.262, 13.883, 12.7 }, input.y);
    try std.testing.expectEqualSlices(f64, &.{ 5.595, 6.312, 7.0 }, input.z);
    try std.testing.expectEqualSlices(f64, &.{ 1.7, 1.55, 1.52 }, input.r);
}

test "readAtomInputFromFile nonexistent file" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;

    const result = readAtomInputFromFile(allocator, io, "nonexistent_file.json");
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

    const count = try checkDuplicateCoordinatesOptions(allocator, input, .{ .log_warning = false });
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

fn expectJsonError(expected: anyerror, json: []const u8) !void {
    try std.testing.expectError(expected, parseAtomInput(std.testing.allocator, json));
}

test "parseAtomInput rejects strings where numbers are required" {
    // Coordinates and radii: a numeric string used to be converted.
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [0, "3"], "y": [0, 0], "z": [0, 0], "r": [1, 1]}
    );
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [0, 0], "y": ["1.5", 0], "z": [0, 0], "r": [1, 1]}
    );
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [0, 0], "y": [0, 0], "z": [0, "abc"], "r": [1, 1]}
    );
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [0, 0], "y": [0, 0], "z": [0, 0], "r": ["1.7", "1.7"]}
    );
    // Other JSON types are no numbers either.
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [true, 0], "y": [0, 0], "z": [0, 0], "r": [1, 1]}
    );
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [null, 0], "y": [0, 0], "z": [0, 0], "r": [1, 1]}
    );
    try expectJsonError(error.ExpectedNumber,
        \\{"x": [[0], 0], "y": [0, 0], "z": [0, 0], "r": [1, 1]}
    );
    // A field that is not an array.
    try expectJsonError(error.ExpectedArray,
        \\{"x": 1, "y": [0], "z": [0], "r": [1]}
    );
    try expectJsonError(error.ExpectedArray,
        \\{"x": "1,2", "y": [0], "z": [0], "r": [1]}
    );
    try expectJsonError(error.ExpectedArray,
        \\{"x": null, "y": [0], "z": [0], "r": [1]}
    );
}

test "parseAtomInput rejects element values that are not atomic numbers" {
    const head = "{\"x\": [0, 0], \"y\": [0, 0], \"z\": [0, 0], \"r\": [1, 1], \"element\": ";
    // A string was read as its bytes: "CN" became the atomic numbers 67 and 78.
    try expectJsonError(error.ExpectedArray, head ++ "\"CN\"}");
    try expectJsonError(error.InvalidElement, head ++ "[6, \"7\"]}");
    try expectJsonError(error.InvalidElement, head ++ "[\"C\", \"N\"]}");
    try expectJsonError(error.InvalidElement, head ++ "[6, 7.5]}");
    try expectJsonError(error.InvalidElement, head ++ "[6, -1]}");
    try expectJsonError(error.InvalidElement, head ++ "[6, 256]}");
    try expectJsonError(error.InvalidElement, head ++ "[6, true]}");
    try expectJsonError(error.InvalidElement, head ++ "[6, null]}");
    try expectJsonError(error.ArrayLengthMismatch, head ++ "[6]}");
}

test "parseAtomInput accepts the documented element forms" {
    const allocator = std.testing.allocator;
    const head = "{\"x\": [0, 1, 2], \"y\": [0, 0, 0], \"z\": [0, 0, 0], \"r\": [1, 1, 1]";

    var ints = try parseAtomInput(allocator, head ++ ", \"element\": [7, 6, 118]}");
    defer ints.deinit();
    try std.testing.expectEqualSlices(u8, &.{ 7, 6, 118 }, ints.element.?);

    // A whole number written with a fraction or exponent is still that number.
    var floats = try parseAtomInput(allocator, head ++ ", \"element\": [6.0, 7e0, 0]}");
    defer floats.deinit();
    try std.testing.expectEqualSlices(u8, &.{ 6, 7, 0 }, floats.element.?);

    // null and a missing field both mean "no element information".
    var null_element = try parseAtomInput(allocator, head ++ ", \"element\": null}");
    defer null_element.deinit();
    try std.testing.expect(null_element.element == null);
}

test "parseAtomInput rejects non-string residue and atom names" {
    const head = "{\"x\": [0], \"y\": [0], \"z\": [0], \"r\": [1], ";
    try expectJsonError(error.ExpectedString, head ++ "\"residue\": [1]}");
    try expectJsonError(error.ExpectedString, head ++ "\"atom_name\": [null]}");
    try expectJsonError(error.ExpectedArray, head ++ "\"residue\": \"ALA\"}");
}

test "parseAtomInput keeps rejecting unknown and repeated fields and trailing text" {
    try expectJsonError(error.UnknownField,
        \\{"x": [0], "y": [0], "z": [0], "r": [1], "radius": [1]}
    );
    try expectJsonError(error.DuplicateField,
        \\{"x": [0], "x": [1], "y": [0], "z": [0], "r": [1]}
    );
    try expectJsonError(error.MissingField,
        \\{"x": [0], "y": [0], "z": [0]}
    );
    try expectJsonError(error.UnexpectedToken,
        \\[0, 1]
    );
    try expectJsonError(error.SyntaxError,
        \\{"x": [0], "y": [0], "z": [0], "r": [1]} {}
    );
}

test "parseAtomInput reads escaped strings and numbers in any order" {
    const allocator = std.testing.allocator;
    var input = try parseAtomInput(allocator,
        \\{"atom_name": ["CA"], "residue": ["ALA"], "r": [1.5e0], "z": [-3], "y": [2], "x": [1E1]}
    );
    defer input.deinit();
    try std.testing.expectEqualStrings("CA", input.atom_name.?[0].slice());
    try std.testing.expectEqualStrings("ALA", input.residue.?[0].slice());
    try std.testing.expectEqual(@as(f64, 10.0), input.x[0]);
    try std.testing.expectEqual(@as(f64, -3.0), input.z[0]);
    try std.testing.expectEqual(@as(f64, 1.5), input.r[0]);
}

test "parseAtomInput releases everything on allocation failure and on errors" {
    const Check = struct {
        fn ok(allocator: Allocator) !void {
            var input = try parseAtomInput(allocator,
                \\{"x": [1, 2], "y": [3, 4], "z": [5, 6], "r": [1.5, 1.6],
                \\ "residue": ["ALA", "GLY"], "atom_name": ["CA", "N"], "element": [6, 7]}
            );
            input.deinit();
        }
        fn bad(allocator: Allocator) !void {
            var input = parseAtomInput(allocator,
                \\{"x": [1, 2], "y": [3, 4], "z": [5, 6], "r": [1.5, 1.6],
                \\ "residue": ["ALA", "GLY"], "atom_name": ["CA", "N"], "element": [6, "7"]}
            ) catch |err| switch (err) {
                error.InvalidElement => return,
                else => |e| return e,
            };
            input.deinit();
            return error.TestUnexpectedResult;
        }
    };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Check.ok, .{});
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Check.bad, .{});
}

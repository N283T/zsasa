//! TOML-based classifier configuration parser.
//!
//! Converts parsed TOML documents into `classifier.Classifier` instances,
//! using the same data structures and error types as the existing FreeSASA
//! configuration parser.
//!
//! ## TOML Format
//!
//! ```toml
//! name = "NACCESS"
//!
//! [types]
//! C_ALI = { radius = 1.87, class = "apolar" }
//! O     = { radius = 1.40, class = "polar" }
//!
//! [[atoms]]
//! residue = "ANY"
//! atom = "CA"
//! type = "C_ALI"
//!
//! [[atoms]]
//! residue = "ALA"
//! atom = "CB"
//! type = "C_ALI"
//! ```
//!
//! ## Usage
//!
//! ```zig
//! const toml_classifier_parser = @import("toml_classifier_parser.zig");
//!
//! var cls = try toml_classifier_parser.parseConfig(allocator, toml_content);
//! defer cls.deinit();
//!
//! const radius = cls.getRadius("ALA", "CA");
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const Classifier = classifier.Classifier;
const AtomClass = classifier.AtomClass;
const toml_parser = @import("toml_parser.zig");
const Value = toml_parser.Value;
const classifier_parser = @import("classifier_parser.zig");

pub const ParseError = classifier_parser.ParseError;
pub const Error = ParseError || Allocator.Error || toml_parser.Error;

/// Intermediate type definition from the [types] table.
const TypeDef = struct {
    radius: f64,
    class: AtomClass,
};

/// Parse a TOML classifier configuration string into a Classifier.
///
/// The TOML document must contain:
/// - An optional top-level `name` string (defaults to "custom")
/// - A `[types]` table mapping type names to `{ radius, class }` inline tables
/// - Zero or more `[[atoms]]` array-of-tables entries with `residue`, `atom`,
///   and `type` string fields
///
/// Returns a Classifier that must be freed with `deinit()`.
pub fn parseConfig(allocator: Allocator, content: []const u8) Error!Classifier {
    var doc = try toml_parser.parse(allocator, content);
    defer doc.deinit();

    const name = doc.getString("name") orelse "custom";

    // Parse [types] section into temporary map
    var types = std.StringHashMap(TypeDef).init(allocator);
    defer types.deinit();

    if (doc.getTable("types")) |types_table| {
        for (types_table.entries) |entry| {
            switch (entry.value) {
                .inline_table => |fields| {
                    const radius = getFloat(fields, "radius") orelse return error.InvalidRadius;
                    const class_str = getString(fields, "class") orelse return error.InvalidClass;
                    const class = parseClass(class_str) orelse return error.InvalidClass;
                    if (types.contains(entry.key)) return error.DuplicateType;
                    try types.put(entry.key, .{ .radius = radius, .class = class });
                },
                else => return error.InvalidTypeDefinition,
            }
        }
    }

    var result = try Classifier.init(allocator, name);
    errdefer result.deinit();

    // Iterate raw array_tables and filter by name. Document.getArrayTables()
    // currently returns all entries regardless of name, so we filter manually.
    for (doc.array_tables) |at| {
        if (!std.mem.eql(u8, at.name, "atoms")) continue;

        const residue = getString(at.entries, "residue") orelse return error.InvalidAtomDefinition;
        const atom_name = getString(at.entries, "atom") orelse return error.InvalidAtomDefinition;
        const type_name = getString(at.entries, "type") orelse return error.InvalidAtomDefinition;

        const type_def = types.get(type_name) orelse return error.UndefinedType;
        try result.addAtom(residue, atom_name, type_def.radius, type_def.class);
    }

    return result;
}

/// Parse a class string into an AtomClass.
fn parseClass(s: []const u8) ?AtomClass {
    if (std.mem.eql(u8, s, "polar")) return .polar;
    if (std.mem.eql(u8, s, "apolar")) return .apolar;
    return null;
}

/// Look up a float value by key in an entry slice.
/// Also accepts integer values, converting them to float.
fn getFloat(entries: []const Value.Entry, key: []const u8) ?f64 {
    for (entries) |entry| {
        if (std.mem.eql(u8, entry.key, key)) {
            switch (entry.value) {
                .float => |f| return f,
                .integer => |i| return @floatFromInt(i),
                else => return null,
            }
        }
    }
    return null;
}

/// Look up a string value by key in an entry slice.
fn getString(entries: []const Value.Entry, key: []const u8) ?[]const u8 {
    for (entries) |entry| {
        if (std.mem.eql(u8, entry.key, key)) {
            switch (entry.value) {
                .string => |s| return s,
                else => return null,
            }
        }
    }
    return null;
}

// =============================================================================
// Tests
// =============================================================================

test "parseConfig minimal TOML" {
    const config =
        \\name = "test"
        \\
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ALA"
        \\atom = "CA"
        \\type = "C_ALI"
    ;
    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();
    try std.testing.expectEqualStrings("test", result.name);
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(AtomClass.apolar, result.getClass("ALA", "CA"));
}

test "parseConfig TOML with ANY fallback" {
    const config =
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\C_CAR = { radius = 1.76, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "CA"
        \\type = "C_ALI"
        \\
        \\[[atoms]]
        \\residue = "CYS"
        \\atom = "CA"
        \\type = "C_CAR"
    ;
    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();
    try std.testing.expectEqual(@as(?f64, 1.76), result.getRadius("CYS", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
}

test "parseConfig TOML default name" {
    const config =
        \\[types]
        \\C = { radius = 1.70, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "C"
        \\type = "C"
    ;
    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();
    try std.testing.expectEqualStrings("custom", result.name);
}

test "parseConfig TOML error: undefined type" {
    const config =
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "CA"
        \\type = "C_UNDEFINED"
    ;
    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.UndefinedType, result);
}

test "parseConfig TOML error: invalid class" {
    const config =
        \\[types]
        \\C = { radius = 1.87, class = "neither" }
    ;
    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.InvalidClass, result);
}

test "parseConfig TOML error: missing atom field" {
    const config =
        \\[types]
        \\C = { radius = 1.70, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ALA"
        \\atom = "CA"
    ;
    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.InvalidAtomDefinition, result);
}

test "parseConfig TOML error: duplicate type" {
    const config =
        \\[types]
        \\C = { radius = 1.70, class = "apolar" }
        \\C = { radius = 2.00, class = "polar" }
    ;
    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.DuplicateType, result);
}

test "parseConfig TOML error: invalid radius (missing)" {
    const config =
        \\[types]
        \\C = { class = "apolar" }
    ;
    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.InvalidRadius, result);
}

test "parseConfig TOML error: invalid type definition (not inline table)" {
    const config =
        \\[types]
        \\C = "not a table"
    ;
    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.InvalidTypeDefinition, result);
}

test "parseConfig TOML multiple types and atoms" {
    const config =
        \\name = "NACCESS"
        \\
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\C_CAR = { radius = 1.76, class = "apolar" }
        \\N_AMD = { radius = 1.65, class = "polar" }
        \\O = { radius = 1.40, class = "polar" }
        \\S = { radius = 1.85, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "C"
        \\type = "C_CAR"
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "O"
        \\type = "O"
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "CA"
        \\type = "C_ALI"
        \\
        \\[[atoms]]
        \\residue = "CYS"
        \\atom = "SG"
        \\type = "S"
        \\
        \\[[atoms]]
        \\residue = "ARG"
        \\atom = "NE"
        \\type = "N_AMD"
    ;
    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();
    try std.testing.expectEqualStrings("NACCESS", result.name);
    try std.testing.expectEqual(@as(?f64, 1.76), result.getRadius("ALA", "C"));
    try std.testing.expectEqual(@as(?f64, 1.40), result.getRadius("ALA", "O"));
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.85), result.getRadius("CYS", "SG"));
    try std.testing.expectEqual(@as(?f64, 1.65), result.getRadius("ARG", "NE"));
    try std.testing.expectEqual(AtomClass.apolar, result.getClass("ALA", "C"));
    try std.testing.expectEqual(AtomClass.polar, result.getClass("ALA", "O"));
}

test "TOML and FreeSASA configs produce identical classifier" {
    const classifier_parser_mod = @import("classifier_parser.zig");

    const freesasa_config =
        \\name: roundtrip
        \\
        \\types:
        \\C_ALI 1.87 apolar
        \\C_CAR 1.76 apolar
        \\N 1.65 polar
        \\O 1.40 polar
        \\S 1.85 apolar
        \\
        \\atoms:
        \\ANY C   C_CAR
        \\ANY O   O
        \\ANY CA  C_ALI
        \\ANY N   N
        \\ALA CB  C_ALI
        \\CYS SG  S
    ;

    const toml_config =
        \\name = "roundtrip"
        \\
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\C_CAR = { radius = 1.76, class = "apolar" }
        \\N = { radius = 1.65, class = "polar" }
        \\O = { radius = 1.40, class = "polar" }
        \\S = { radius = 1.85, class = "apolar" }
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "C"
        \\type = "C_CAR"
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "O"
        \\type = "O"
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "CA"
        \\type = "C_ALI"
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "N"
        \\type = "N"
        \\
        \\[[atoms]]
        \\residue = "ALA"
        \\atom = "CB"
        \\type = "C_ALI"
        \\
        \\[[atoms]]
        \\residue = "CYS"
        \\atom = "SG"
        \\type = "S"
    ;

    const allocator = std.testing.allocator;
    var fs_result = try classifier_parser_mod.parseConfig(allocator, freesasa_config);
    defer fs_result.deinit();
    var toml_result = try parseConfig(allocator, toml_config);
    defer toml_result.deinit();

    // Names must match
    try std.testing.expectEqualStrings(fs_result.name, toml_result.name);

    // Test same atoms produce same radii and classes
    const test_cases = [_]struct { []const u8, []const u8 }{
        .{ "ALA", "C" },
        .{ "ALA", "O" },
        .{ "ALA", "CA" },
        .{ "ALA", "N" },
        .{ "ALA", "CB" },
        .{ "CYS", "SG" },
        .{ "GLY", "C" },
        .{ "GLY", "CA" },
    };

    for (test_cases) |tc| {
        const residue = tc[0];
        const atom = tc[1];
        try std.testing.expectEqual(
            fs_result.getRadius(residue, atom),
            toml_result.getRadius(residue, atom),
        );
        try std.testing.expectEqual(
            fs_result.getClass(residue, atom),
            toml_result.getClass(residue, atom),
        );
    }
}

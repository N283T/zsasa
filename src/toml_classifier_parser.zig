//! TOML-based classifier configuration parser.
//!
//! Converts parsed TOML documents into `classifier.Classifier` instances,
//! independent of the legacy FreeSASA-style custom classifier parser.
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
const json_parser = @import("json_parser.zig");
const Value = toml_parser.Value;

/// Where a parse error was found: the 1-based line, or 0 when unknown.
pub const Diagnostic = toml_parser.Diagnostic;

pub const ParseError = error{
    /// Type definition is missing required fields or is not an inline table.
    InvalidTypeDefinition,
    /// Atom definition is missing residue, atom name, or type reference.
    InvalidAtomDefinition,
    /// Radius value is missing or not a number.
    InvalidRadius,
    /// Radius is not finite, not positive or larger than 100 Angstroms
    /// (the range the JSON input accepts).
    RadiusOutOfRange,
    /// Class value is not "polar" or "apolar".
    InvalidClass,
    /// Referenced type not defined in the [types] section.
    UndefinedType,
    /// Duplicate type definition.
    DuplicateType,
    /// A key appears twice in the same table, array-of-tables entry or inline table.
    DuplicateKey,
    /// A `[table]` header appears twice.
    DuplicateTable,
};
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
    return parseConfigDiag(allocator, content, null);
}

/// Like `parseConfig`, but stores the line of the offending input in `diag`
/// when the content is rejected.
pub fn parseConfigDiag(allocator: Allocator, content: []const u8, diag: ?*Diagnostic) Error!Classifier {
    var local: Diagnostic = .{};
    const d = diag orelse &local;
    d.* = .{};

    var doc = try toml_parser.parseDiag(allocator, content, d);
    defer doc.deinit();

    try rejectDuplicates(doc, d);

    const name = doc.getString("name") orelse "custom";

    // Parse [types] section into temporary map
    var types: std.StringHashMapUnmanaged(TypeDef) = .empty;
    defer types.deinit(allocator);

    if (doc.getTable("types")) |types_table| {
        for (types_table.entries) |entry| {
            switch (entry.value) {
                .inline_table => |fields| {
                    d.line = entry.line;
                    const radius = getFloat(fields, "radius") orelse return error.InvalidRadius;
                    // Same range as JSON input radii; also rules out values
                    // that overflow areas to infinity.
                    if (!json_parser.isValidRadius(radius)) return error.RadiusOutOfRange;
                    const class_str = getString(fields, "class") orelse return error.InvalidClass;
                    const class = parseClass(class_str) orelse return error.InvalidClass;
                    if (types.contains(entry.key)) return error.DuplicateType;
                    try types.put(allocator, entry.key, .{ .radius = radius, .class = class });
                },
                else => {
                    d.line = entry.line;
                    return error.InvalidTypeDefinition;
                },
            }
        }
    }

    var result = try Classifier.init(allocator, name);
    errdefer result.deinit();

    // Iterate raw array_tables and filter by name. Document.getArrayTables()
    // currently returns all entries regardless of name, so we filter manually.
    for (doc.array_tables) |at| {
        if (!std.mem.eql(u8, at.name, "atoms")) continue;
        d.line = at.line;

        const residue = getString(at.entries, "residue") orelse return error.InvalidAtomDefinition;
        const atom_name = getString(at.entries, "atom") orelse return error.InvalidAtomDefinition;
        const type_name = getString(at.entries, "type") orelse return error.InvalidAtomDefinition;

        const type_def = types.get(type_name) orelse return error.UndefinedType;
        try result.addAtom(residue, atom_name, type_def.radius, type_def.class);
    }

    return result;
}

/// Reject keys that appear twice in one table, array-of-tables entry or
/// inline table, and `[table]` headers that appear twice. The first value
/// used to win silently.
fn rejectDuplicates(doc: toml_parser.Document, diag: *Diagnostic) ParseError!void {
    for (doc.tables, 0..) |table, i| {
        for (doc.tables[0..i]) |earlier| {
            if (std.mem.eql(u8, earlier.name, table.name)) {
                diag.line = table.line;
                return error.DuplicateTable;
            }
        }
        // A type name defined twice in [types] is reported as DuplicateType.
        try rejectDuplicateKeys(table.entries, diag, std.mem.eql(u8, table.name, "types"));
    }
    for (doc.array_tables) |at| try rejectDuplicateKeys(at.entries, diag, false);
}

fn rejectDuplicateKeys(entries: []const Value.Entry, diag: *Diagnostic, skip_entry_keys: bool) ParseError!void {
    for (entries, 0..) |entry, i| {
        if (!skip_entry_keys) {
            for (entries[0..i]) |earlier| {
                if (std.mem.eql(u8, earlier.key, entry.key)) {
                    diag.line = entry.line;
                    return error.DuplicateKey;
                }
            }
        }
        switch (entry.value) {
            .inline_table => |fields| {
                for (fields, 0..) |field, j| {
                    for (fields[0..j]) |earlier_field| {
                        if (std.mem.eql(u8, earlier_field.key, field.key)) {
                            diag.line = entry.line;
                            return error.DuplicateKey;
                        }
                    }
                }
            },
            else => {},
        }
    }
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

// Line 1 is `name`, line 3 `[types]`, line 4 the C_ALI type, line 6 `[[atoms]]`.
const test_head = "name = \"t\"\n\n[types]\nC_ALI = { radius = 1.87, class = \"apolar\" }\n";
const test_atoms = "\n[[atoms]]\nresidue = \"ANY\"\natom = \"CA\"\ntype = \"C_ALI\"\n";

fn typeLine(comptime definition: []const u8) []const u8 {
    return "name = \"t\"\n\n[types]\n" ++ definition ++ "\n" ++ test_atoms;
}

fn expectConfigError(expected: anyerror, expected_line: usize, content: []const u8) !void {
    var diag: Diagnostic = .{};
    try std.testing.expectError(expected, parseConfigDiag(std.testing.allocator, content, &diag));
    try std.testing.expectEqual(expected_line, diag.line);
}

test "parseConfig rejects radii that are not finite or outside 0 to 100" {
    // 1e300 overflowed the areas to infinity, which was written to the JSON output as `inf`.
    const bad = [_][]const u8{ "1e300", "inf", "-inf", "nan", "-1.0", "0", "0.0", "100.5", "-0.5" };
    inline for (bad) |radius| {
        try expectConfigError(
            error.RadiusOutOfRange,
            4,
            comptime typeLine("C_ALI = { radius = " ++ radius ++ ", class = \"apolar\" }"),
        );
    }

    // The bounds of the accepted range, and an integer radius.
    const ok = [_][]const u8{ "100.0", "0.5", "2", "1e1" };
    inline for (ok) |radius| {
        var cls = try parseConfig(
            std.testing.allocator,
            comptime typeLine("C_ALI = { radius = " ++ radius ++ ", class = \"apolar\" }"),
        );
        cls.deinit();
    }
}

test "parseConfig keeps reporting a missing or non-numeric radius as InvalidRadius" {
    try expectConfigError(error.InvalidRadius, 4, typeLine("C_ALI = { class = \"apolar\" }"));
    try expectConfigError(error.InvalidRadius, 4, typeLine("C_ALI = { radius = \"1.8\", class = \"apolar\" }"));
}

test "parseConfig rejects text after a value" {
    // After a string
    try expectConfigError(error.UnexpectedCharacter, 1, "name = \"t\" junk\n\n[types]\nC_ALI = { radius = 1.87, class = \"apolar\" }\n" ++ test_atoms);
    // After an inline table
    try expectConfigError(error.UnexpectedCharacter, 4, typeLine("C_ALI = { radius = 1.87, class = \"apolar\" } extra"));
    // After a string inside an inline table
    try expectConfigError(error.UnexpectedCharacter, 4, typeLine("C_ALI = { radius = 1.87, class = \"apolar\" junk }"));
    // In an [[atoms]] entry
    try expectConfigError(
        error.UnexpectedCharacter,
        7,
        test_head ++ "\n[[atoms]]\nresidue = \"ANY\" x\natom = \"CA\"\ntype = \"C_ALI\"\n",
    );
    // A comment after a value is not text after a value.
    var commented = try parseConfig(
        std.testing.allocator,
        "name = \"t\" # note\n\n[types]\nC_ALI = { radius = 1.87, class = \"apolar\" } # note\n" ++ test_atoms ++ "# end\n",
    );
    commented.deinit();
}

test "parseConfig rejects duplicate keys and sections with their line" {
    try expectConfigError(error.DuplicateKey, 2, "name = \"a\"\nname = \"b\"\n\n[types]\nC_ALI = { radius = 1.87, class = \"apolar\" }\n" ++ test_atoms);
    try expectConfigError(error.DuplicateTable, 6, test_head ++ "\n[types]\nO = { radius = 1.4, class = \"polar\" }\n" ++ test_atoms);
    try expectConfigError(error.DuplicateKey, 4, typeLine("C_ALI = { radius = 1.87, radius = 1.0, class = \"apolar\" }"));
    try expectConfigError(
        error.DuplicateKey,
        9,
        test_head ++ "\n[[atoms]]\nresidue = \"ANY\"\natom = \"CA\"\natom = \"CB\"\ntype = \"C_ALI\"\n",
    );
    // A type name that is defined twice
    try expectConfigError(error.DuplicateType, 5, test_head ++ "C_ALI = { radius = 1.5, class = \"polar\" }\n" ++ test_atoms);
    // Repeated [[atoms]] entries are the normal way to list atoms.
    var cls = try parseConfig(std.testing.allocator, test_head ++ test_atoms ++ test_atoms);
    cls.deinit();
}

test "parseConfig syntax errors carry the line" {
    try expectConfigError(error.UnexpectedCharacter, 3, "name = \"t\"\n\nnot a pair\n");
    try expectConfigError(error.UnterminatedString, 2, "name = \"t\"\nname2 = \"open\n");
}

test "parseConfig releases memory when it rejects a document" {
    const Check = struct {
        fn ok(allocator: Allocator) !void {
            var cls = try parseConfig(allocator, test_head ++ test_atoms ++ test_atoms);
            cls.deinit();
        }
        fn rejected(allocator: Allocator) !void {
            var cls = parseConfig(allocator, test_head ++ "O = { radius = 1e300, class = \"polar\" }\n" ++ test_atoms) catch |err| switch (err) {
                error.RadiusOutOfRange => return,
                else => |e| return e,
            };
            cls.deinit();
            return error.TestUnexpectedResult;
        }
        fn duplicate(allocator: Allocator) !void {
            var cls = parseConfig(allocator, test_head ++ test_atoms ++ "\n[types]\nO = { radius = 1.4, class = \"polar\" }\n") catch |err| switch (err) {
                error.DuplicateTable => return,
                else => |e| return e,
            };
            cls.deinit();
            return error.TestUnexpectedResult;
        }
    };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Check.ok, .{});
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Check.rejected, .{});
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Check.duplicate, .{});
}

//! FreeSASA-compatible configuration file parser.
//!
//! This module parses classifier configuration files in the FreeSASA format,
//! allowing users to define custom atom types and radii.
//!
//! ## File Format
//!
//! ```text
//! # Comment
//! name: NACCESS
//!
//! types:
//! C_ALI 1.87 apolar
//! C_CAR 1.76 apolar
//! O     1.40 polar
//!
//! atoms:
//! ANY C   C_CAR
//! ANY CA  C_ALI
//! ALA CB  C_ALI
//! ```
//!
//! - Lines starting with `#` are comments
//! - `name:` section defines the classifier name (optional)
//! - `types:` section defines atom types with radius and polarity class
//! - `atoms:` section maps (residue, atom_name) pairs to defined types
//!
//! ## Usage
//!
//! ```zig
//! const parser = @import("classifier_parser.zig");
//!
//! // Parse from string
//! var classifier = try parser.parseConfig(allocator, config_content);
//! defer classifier.deinit();
//!
//! // Parse from file
//! var classifier = try parser.parseConfigFile(allocator, "naccess.config");
//! defer classifier.deinit();
//!
//! // Use the classifier
//! const radius = classifier.getRadius("ALA", "CA");
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const Classifier = classifier.Classifier;
const AtomClass = classifier.AtomClass;
const AtomProperties = classifier.AtomProperties;
const toml_classifier_parser = @import("toml_classifier_parser.zig");
const toml_parser = @import("toml_parser.zig");

/// Parse error with line information
pub const ParseError = error{
    /// Unknown section header (not "name:", "types:", or "atoms:")
    UnknownSection,
    /// Type definition is missing radius or class
    InvalidTypeDefinition,
    /// Atom definition is missing residue, atom name, or type reference
    InvalidAtomDefinition,
    /// Radius value is not a valid float
    InvalidRadius,
    /// Class value is not "polar" or "apolar"
    InvalidClass,
    /// Referenced type not defined in types section
    UndefinedType,
    /// Duplicate type definition
    DuplicateType,
    /// atoms section appeared before types section
    AtomsBeforeTypes,
    /// File read error
    FileReadError,
    /// Out of memory
    OutOfMemory,
};

/// All possible errors from parsing (includes TOML parser errors for .toml files)
pub const Error = ParseError || Allocator.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError || toml_parser.TomlError;

/// Result of parsing with line number context for error reporting
pub const ParseResult = struct {
    classifier: Classifier,
    /// Line number of last successfully parsed line (for debugging)
    lines_parsed: usize,
};

/// Atom type definition (intermediate storage during parsing)
const TypeDef = struct {
    radius: f64,
    class: AtomClass,
};

/// Current section being parsed
const Section = enum {
    none,
    types,
    atoms,
};

/// Parse a FreeSASA configuration file from a string.
///
/// The configuration file format:
/// - Lines starting with `#` are comments
/// - `name: NAME` defines the classifier name (default: "custom")
/// - `types:` section header, followed by type definitions
/// - `atoms:` section header, followed by atom mappings
///
/// Returns a Classifier that must be freed with `deinit()`.
pub fn parseConfig(allocator: Allocator, content: []const u8) Error!Classifier {
    const result = try parseConfigWithLineInfo(allocator, content);
    return result.classifier;
}

/// Parse configuration with line number tracking for better error context.
pub fn parseConfigWithLineInfo(allocator: Allocator, content: []const u8) Error!ParseResult {
    // Temporary storage for type definitions
    var types = std.StringHashMap(TypeDef).init(allocator);
    defer types.deinit();

    // We need to track allocated keys to free them
    var type_keys = std.ArrayListUnmanaged([]u8).empty;
    defer {
        for (type_keys.items) |key| {
            allocator.free(key);
        }
        type_keys.deinit(allocator);
    }

    var name: []const u8 = "custom";
    var section: Section = .none;
    var line_num: usize = 0;
    var types_seen = false;

    // First pass: collect types and name
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |raw_line| {
        line_num += 1;

        // Strip comments and whitespace
        const line = stripComment(raw_line);
        const trimmed = std.mem.trim(u8, line, " \t\r");

        if (trimmed.len == 0) continue;

        // Check for section headers
        if (std.mem.startsWith(u8, trimmed, "name:")) {
            const value = std.mem.trim(u8, trimmed[5..], " \t");
            if (value.len > 0) {
                name = value;
            }
            continue;
        }

        if (std.mem.eql(u8, trimmed, "types:")) {
            section = .types;
            types_seen = true;
            continue;
        }

        if (std.mem.eql(u8, trimmed, "atoms:")) {
            section = .atoms;
            continue;
        }

        // Parse based on current section
        switch (section) {
            .none => {
                // Unknown line outside section - ignore or could be error
                continue;
            },
            .types => {
                // Parse type definition: TYPE_NAME RADIUS CLASS
                var parts = std.mem.tokenizeAny(u8, trimmed, " \t");

                const type_name = parts.next() orelse continue;
                const radius_str = parts.next() orelse return error.InvalidTypeDefinition;
                const class_str = parts.next() orelse return error.InvalidTypeDefinition;

                const radius = std.fmt.parseFloat(f64, radius_str) catch {
                    return error.InvalidRadius;
                };

                const class = parseClass(class_str) orelse return error.InvalidClass;

                // Check for duplicate
                if (types.contains(type_name)) {
                    return error.DuplicateType;
                }

                // Store type (need to dupe the key since it's from content slice)
                const key = try allocator.dupe(u8, type_name);
                try type_keys.append(allocator, key);
                try types.put(key, TypeDef{ .radius = radius, .class = class });
            },
            .atoms => {
                // Skip atoms in first pass
                continue;
            },
        }
    }

    // Create classifier
    var result = try Classifier.init(allocator, name);
    errdefer result.deinit();

    // Second pass: process atoms section
    section = .none;
    lines = std.mem.splitScalar(u8, content, '\n');
    line_num = 0;

    while (lines.next()) |raw_line| {
        line_num += 1;

        const line = stripComment(raw_line);
        const trimmed = std.mem.trim(u8, line, " \t\r");

        if (trimmed.len == 0) continue;

        if (std.mem.eql(u8, trimmed, "types:")) {
            section = .types;
            continue;
        }

        if (std.mem.eql(u8, trimmed, "atoms:")) {
            if (!types_seen) {
                return error.AtomsBeforeTypes;
            }
            section = .atoms;
            continue;
        }

        if (section == .atoms) {
            // Parse atom definition: RESIDUE ATOM_NAME TYPE_NAME
            var parts = std.mem.tokenizeAny(u8, trimmed, " \t");

            const residue = parts.next() orelse continue;
            const atom_name = parts.next() orelse return error.InvalidAtomDefinition;
            const type_name = parts.next() orelse return error.InvalidAtomDefinition;

            // Look up type
            const type_def = types.get(type_name) orelse return error.UndefinedType;

            // Add to classifier
            try result.addAtom(residue, atom_name, type_def.radius, type_def.class);
        }
    }

    return ParseResult{
        .classifier = result,
        .lines_parsed = line_num,
    };
}

/// Parse a classifier configuration file from disk.
///
/// Format is auto-detected by file extension:
/// - `.toml` files are parsed as TOML format
/// - All other extensions use the FreeSASA config format
///
/// The file is read entirely into memory before parsing.
pub fn parseConfigFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Classifier {
    // Auto-detect format by extension
    if (std.mem.endsWith(u8, path, ".toml")) {
        return parseTomlConfigFile(allocator, io, path);
    }

    // Default: FreeSASA format
    return parseFreeSasaConfigFile(allocator, io, path);
}

/// Parse a TOML classifier configuration file from disk.
fn parseTomlConfigFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Classifier {
    const content = try readFileContent(allocator, io, path);
    defer allocator.free(content);
    return toml_classifier_parser.parseConfig(allocator, content);
}

/// Parse a FreeSASA classifier configuration file from disk.
fn parseFreeSasaConfigFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Classifier {
    const content = try readFileContent(allocator, io, path);
    defer allocator.free(content);
    return parseConfig(allocator, content);
}

/// Read entire file into a heap-allocated buffer.
fn readFileContent(allocator: Allocator, io: std.Io, path: []const u8) Error![]u8 {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch {
        return error.FileReadError;
    };
    defer file.close(io);

    var read_buf: [65536]u8 = undefined;
    var r = file.reader(io, &read_buf);
    return r.interface.allocRemaining(allocator, .unlimited) catch {
        return error.FileReadError;
    };
}

/// Strip comment from a line (everything after # is removed)
fn stripComment(line: []const u8) []const u8 {
    if (std.mem.findScalar(u8, line, '#')) |idx| {
        return line[0..idx];
    }
    return line;
}

/// Parse class string to AtomClass
fn parseClass(s: []const u8) ?AtomClass {
    if (std.mem.eql(u8, s, "polar")) return .polar;
    if (std.mem.eql(u8, s, "apolar")) return .apolar;
    return null;
}

// =============================================================================
// Tests
// =============================================================================

test "parseConfig minimal" {
    const config =
        \\name: test
        \\types:
        \\C_ALI 1.87 apolar
        \\atoms:
        \\ALA CA C_ALI
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    try std.testing.expectEqualStrings("test", result.name);
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(AtomClass.apolar, result.getClass("ALA", "CA"));
}

test "parseConfig with comments" {
    const config =
        \\# This is a comment
        \\name: test
        \\
        \\types:
        \\C_ALI 1.87 apolar  # inline comment
        \\# Another comment
        \\O 1.40 polar
        \\
        \\atoms:
        \\ANY CA C_ALI
        \\ANY O O
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.40), result.getRadius("GLY", "O"));
}

test "parseConfig ANY fallback" {
    const config =
        \\types:
        \\C_ALI 1.87 apolar
        \\C_CAR 1.76 apolar
        \\atoms:
        \\ANY CA C_ALI
        \\CYS CA C_CAR
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    // CYS has specific entry
    try std.testing.expectEqual(@as(?f64, 1.76), result.getRadius("CYS", "CA"));
    // Other residues use ANY fallback
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("GLY", "CA"));
}

test "parseConfig default name" {
    const config =
        \\types:
        \\C 1.70 apolar
        \\atoms:
        \\ANY C C
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    try std.testing.expectEqualStrings("custom", result.name);
}

test "parseConfig error: undefined type" {
    const config =
        \\types:
        \\C_ALI 1.87 apolar
        \\atoms:
        \\ANY CA C_UNDEFINED
    ;

    const allocator = std.testing.allocator;
    const result = parseConfig(allocator, config);
    try std.testing.expectError(error.UndefinedType, result);
}

test "parseConfig error: invalid radius" {
    const config =
        \\types:
        \\C_ALI notanumber apolar
    ;

    const allocator = std.testing.allocator;
    const result = parseConfig(allocator, config);
    try std.testing.expectError(error.InvalidRadius, result);
}

test "parseConfig error: invalid class" {
    const config =
        \\types:
        \\C_ALI 1.87 neither
    ;

    const allocator = std.testing.allocator;
    const result = parseConfig(allocator, config);
    try std.testing.expectError(error.InvalidClass, result);
}

test "parseConfig error: duplicate type" {
    const config =
        \\types:
        \\C_ALI 1.87 apolar
        \\C_ALI 2.00 polar
    ;

    const allocator = std.testing.allocator;
    const result = parseConfig(allocator, config);
    try std.testing.expectError(error.DuplicateType, result);
}

test "parseConfig error: invalid type definition (missing radius)" {
    const config =
        \\types:
        \\C_ALI
    ;

    const allocator = std.testing.allocator;
    const result = parseConfig(allocator, config);
    try std.testing.expectError(error.InvalidTypeDefinition, result);
}

test "parseConfig error: invalid atom definition (missing type)" {
    const config =
        \\types:
        \\C_ALI 1.87 apolar
        \\atoms:
        \\ANY CA
    ;

    const allocator = std.testing.allocator;
    const result = parseConfig(allocator, config);
    try std.testing.expectError(error.InvalidAtomDefinition, result);
}

test "parseConfig multiple types and atoms" {
    const config =
        \\name: NACCESS
        \\
        \\types:
        \\C_ALI 1.87 apolar
        \\C_CAR 1.76 apolar
        \\N_AMD 1.65 polar
        \\O 1.40 polar
        \\S 1.85 apolar
        \\
        \\atoms:
        \\ANY C   C_CAR
        \\ANY O   O
        \\ANY CA  C_ALI
        \\ANY N   N_AMD
        \\ANY CB  C_ALI
        \\ALA CB  C_ALI
        \\CYS SG  S
        \\ARG NE  N_AMD
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    try std.testing.expectEqualStrings("NACCESS", result.name);

    // Check various lookups
    try std.testing.expectEqual(@as(?f64, 1.76), result.getRadius("ALA", "C"));
    try std.testing.expectEqual(@as(?f64, 1.40), result.getRadius("ALA", "O"));
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.65), result.getRadius("ALA", "N"));
    try std.testing.expectEqual(@as(?f64, 1.85), result.getRadius("CYS", "SG"));
    try std.testing.expectEqual(@as(?f64, 1.65), result.getRadius("ARG", "NE"));

    // Check classes
    try std.testing.expectEqual(AtomClass.apolar, result.getClass("ALA", "C"));
    try std.testing.expectEqual(AtomClass.polar, result.getClass("ALA", "O"));
    try std.testing.expectEqual(AtomClass.polar, result.getClass("ALA", "N"));
    try std.testing.expectEqual(AtomClass.apolar, result.getClass("CYS", "SG"));
}

test "parseConfig with extra whitespace" {
    // Test with extra spaces (tabs are tokenized the same way via tokenizeAny)
    const config =
        \\types:
        \\  C_ALI   1.87   apolar
        \\  O     1.40    polar
        \\atoms:
        \\  ANY   CA   C_ALI
        \\  ANY   O   O
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.40), result.getRadius("ALA", "O"));
}

test "parseConfig empty sections" {
    const config =
        \\types:
        \\atoms:
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfig(allocator, config);
    defer result.deinit();

    try std.testing.expectEqual(@as(usize, 0), result.count());
}

test "stripComment" {
    try std.testing.expectEqualStrings("hello", stripComment("hello"));
    try std.testing.expectEqualStrings("hello ", stripComment("hello # world"));
    try std.testing.expectEqualStrings("", stripComment("# comment only"));
    try std.testing.expectEqualStrings("test ", stripComment("test # inline # multiple"));
}

test "parseClass" {
    try std.testing.expectEqual(AtomClass.polar, parseClass("polar").?);
    try std.testing.expectEqual(AtomClass.apolar, parseClass("apolar").?);
    try std.testing.expectEqual(@as(?AtomClass, null), parseClass("unknown"));
}

test "parseConfigWithLineInfo" {
    const config =
        \\name: test
        \\types:
        \\C 1.70 apolar
        \\atoms:
        \\ANY C C
    ;

    const allocator = std.testing.allocator;
    var result = try parseConfigWithLineInfo(allocator, config);
    defer result.classifier.deinit();

    try std.testing.expectEqual(@as(usize, 5), result.lines_parsed);
}

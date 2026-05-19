//! TOML-only custom classifier configuration parser.
//!
//! This module keeps the historical custom classifier parser entry points while
//! requiring custom classifier files and inline config content to use the TOML
//! format. Legacy FreeSASA-style custom classifier configs are rejected with a
//! dedicated error so callers can print migration guidance.
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
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const Classifier = classifier.Classifier;
const AtomClass = classifier.AtomClass;
const toml_classifier_parser = @import("toml_classifier_parser.zig");

pub const ParseError = toml_classifier_parser.ParseError;
pub const Error = toml_classifier_parser.Error || error{
    UnsupportedLegacyFormat,
    UnsupportedConfigExtension,
    FileReadError,
};

/// Parse a TOML classifier configuration string.
///
/// Legacy FreeSASA-style configs whose first meaningful line is `name:`,
/// `types:`, or `atoms:` are rejected with `error.UnsupportedLegacyFormat`.
pub fn parseConfig(allocator: Allocator, content: []const u8) Error!Classifier {
    if (isLegacyFreeSasaConfig(content)) {
        return error.UnsupportedLegacyFormat;
    }
    return toml_classifier_parser.parseConfig(allocator, content);
}

/// Parse a TOML classifier configuration file from disk.
///
/// Custom classifier files are TOML-only; paths not ending in `.toml` are
/// rejected before any file I/O is attempted.
pub fn parseConfigFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Classifier {
    if (!isTomlPath(path)) {
        return error.UnsupportedConfigExtension;
    }

    const content = try readFileContent(allocator, io, path);
    defer allocator.free(content);
    return parseConfig(allocator, content);
}

fn isTomlPath(path: []const u8) bool {
    return std.mem.endsWith(u8, path, ".toml");
}

fn isLegacyFreeSasaConfig(content: []const u8) bool {
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |raw_line| {
        const line = stripComment(raw_line);
        const trimmed = std.mem.trim(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        if (std.mem.startsWith(u8, trimmed, "name:") or
            std.mem.eql(u8, trimmed, "types:") or
            std.mem.eql(u8, trimmed, "atoms:"))
        {
            return true;
        }
    }
    return false;
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

/// Strip comment from a line (everything after # is removed).
fn stripComment(line: []const u8) []const u8 {
    if (std.mem.findScalar(u8, line, '#')) |idx| {
        return line[0..idx];
    }
    return line;
}

// =============================================================================
// Tests
// =============================================================================

test "custom classifier parser accepts TOML content" {
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

test "custom classifier parser rejects legacy FreeSASA content" {
    const config =
        \\# legacy FreeSASA-style config
        \\name: test
        \\types:
        \\C_ALI 1.87 apolar
        \\atoms:
        \\ALA CA C_ALI
    ;

    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.UnsupportedLegacyFormat, result);
}

test "custom classifier parser rejects legacy marker after other content" {
    const config =
        \\legacy-header
        \\types:
        \\C_ALI 1.87 apolar
        \\atoms:
        \\ALA CA C_ALI
    ;

    const result = parseConfig(std.testing.allocator, config);
    try std.testing.expectError(error.UnsupportedLegacyFormat, result);
}

test "custom classifier file extension must be TOML" {
    try std.testing.expect(isTomlPath("my_radii.toml"));
    try std.testing.expect(!isTomlPath("my_radii.config"));
    try std.testing.expect(!isTomlPath("my_radii.conf"));

    const result_config = parseConfigFile(std.testing.allocator, std.testing.io, "my_radii.config");
    try std.testing.expectError(error.UnsupportedConfigExtension, result_config);

    const result_conf = parseConfigFile(std.testing.allocator, std.testing.io, "my_radii.conf");
    try std.testing.expectError(error.UnsupportedConfigExtension, result_conf);
}

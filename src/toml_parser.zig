//! Minimal TOML subset parser.
//!
//! Supports the subset of TOML needed for classifier configuration files:
//! - Top-level key-value pairs
//! - `[table]` sections
//! - `[[array_of_tables]]` sections
//! - Inline tables `{ key = value, ... }`
//! - `"string"` values (with `\"` and `\\` escape handling)
//! - Float and integer values
//! - `#` comments (line-end and full-line)
//!
//! All string slices point into the original input buffer; only inline table
//! entry slices are heap-allocated.
//!
//! ## Usage
//!
//! ```zig
//! const toml = @import("toml_parser.zig");
//! var doc = try toml.parse(allocator, content);
//! defer doc.deinit();
//!
//! const name = doc.getString("name");
//! const types = doc.getTable("types");
//! const atoms = doc.getArrayTables("atoms");
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const TomlError = error{
    UnterminatedString,
    InvalidNumber,
    InvalidInlineTable,
    ExpectedEquals,
    ExpectedValue,
    UnexpectedCharacter,
};

pub const Error = TomlError || Allocator.Error;

/// A TOML value.
pub const Value = union(enum) {
    string: []const u8,
    float: f64,
    integer: i64,
    inline_table: []const Entry,

    pub const Entry = struct {
        key: []const u8,
        value: Value,
    };
};

/// A parsed TOML document section.
pub const Table = struct {
    name: []const u8,
    entries: []const Value.Entry,
};

/// Result of parsing a TOML document.
pub const Document = struct {
    /// Top-level tables: root (name=""), [types], etc.
    tables: []Table,
    /// Array-of-tables entries: [[atoms]] produces multiple entries.
    array_tables: []ArrayTable,
    allocator: Allocator,

    pub const ArrayTable = struct {
        name: []const u8,
        entries: []const Value.Entry,
    };

    pub fn deinit(self: *Document) void {
        for (self.tables) |table| {
            for (table.entries) |entry| {
                freeValue(self.allocator, entry.value);
            }
            self.allocator.free(table.entries);
        }
        self.allocator.free(self.tables);
        for (self.array_tables) |at| {
            for (at.entries) |entry| {
                freeValue(self.allocator, entry.value);
            }
            self.allocator.free(at.entries);
        }
        self.allocator.free(self.array_tables);
    }

    /// Look up a top-level (root table) string value by key.
    pub fn getString(self: *const Document, key: []const u8) ?[]const u8 {
        for (self.tables) |table| {
            if (table.name.len == 0) {
                for (table.entries) |entry| {
                    if (std.mem.eql(u8, entry.key, key)) {
                        switch (entry.value) {
                            .string => |s| return s,
                            else => return null,
                        }
                    }
                }
            }
        }
        return null;
    }

    /// Look up a named table.
    pub fn getTable(self: *const Document, name: []const u8) ?Table {
        for (self.tables) |table| {
            if (std.mem.eql(u8, table.name, name)) return table;
        }
        return null;
    }

    /// Return all parsed array-of-tables entries.
    /// Currently returns all entries regardless of name, which works
    /// correctly when only one [[array_table_name]] appears in the document.
    pub fn getArrayTables(self: *const Document, _: []const u8) []const ArrayTable {
        return self.array_tables;
    }
};

/// Recursively free heap-allocated memory owned by a `Value`.
fn freeValue(allocator: Allocator, value: Value) void {
    switch (value) {
        .inline_table => |entries| {
            for (entries) |entry| {
                freeValue(allocator, entry.value);
            }
            allocator.free(entries);
        },
        .string, .float, .integer => {},
    }
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

/// Section kind currently being parsed.
const SectionKind = enum {
    table,
    array_of_tables,
};

/// Parse a TOML document from `content`.
///
/// All returned string slices point into `content`; only inline-table entry
/// arrays are heap-allocated. Call `Document.deinit()` to free all owned
/// memory.
pub fn parse(allocator: Allocator, content: []const u8) Error!Document {
    var tables = std.ArrayListUnmanaged(Table){};
    errdefer {
        for (tables.items) |table| {
            for (table.entries) |entry| {
                freeValue(allocator, entry.value);
            }
            allocator.free(table.entries);
        }
        tables.deinit(allocator);
    }

    var array_tables = std.ArrayListUnmanaged(Document.ArrayTable){};
    errdefer {
        for (array_tables.items) |at| {
            for (at.entries) |entry| {
                freeValue(allocator, entry.value);
            }
            allocator.free(at.entries);
        }
        array_tables.deinit(allocator);
    }

    // Current section state
    var current_name: []const u8 = "";
    var current_kind: SectionKind = .table;
    var current_entries = std.ArrayListUnmanaged(Value.Entry){};
    errdefer {
        for (current_entries.items) |entry| {
            freeValue(allocator, entry.value);
        }
        current_entries.deinit(allocator);
    }

    var line_iter = std.mem.splitScalar(u8, content, '\n');
    while (line_iter.next()) |raw_line| {
        const line = stripComment(raw_line);
        const trimmed = std.mem.trim(u8, line, " \t\r");
        if (trimmed.len == 0) continue;

        // [[array_of_tables]]
        if (std.mem.startsWith(u8, trimmed, "[[")) {
            const end = std.mem.indexOf(u8, trimmed[2..], "]]") orelse
                return error.UnexpectedCharacter;
            // Flush current section
            try flushSection(
                allocator,
                &tables,
                &array_tables,
                current_name,
                current_kind,
                &current_entries,
            );
            current_name = std.mem.trim(u8, trimmed[2 .. 2 + end], " \t");
            current_kind = .array_of_tables;
            continue;
        }

        // [table]
        if (trimmed[0] == '[') {
            if (std.mem.indexOfScalar(u8, trimmed, ']')) |end| {
                // Flush current section
                try flushSection(
                    allocator,
                    &tables,
                    &array_tables,
                    current_name,
                    current_kind,
                    &current_entries,
                );
                current_name = std.mem.trim(u8, trimmed[1..end], " \t");
                current_kind = .table;
                continue;
            }
        }

        // key = value
        if (std.mem.indexOf(u8, trimmed, "=")) |eq_pos| {
            // Split on the first '=' to separate key and value.
            const key_raw = trimmed[0..eq_pos];
            const key = std.mem.trim(u8, key_raw, " \t");
            if (key.len == 0) return error.UnexpectedCharacter;

            const value_raw = trimmed[eq_pos + 1 ..];
            const value_trimmed = std.mem.trim(u8, value_raw, " \t");
            if (value_trimmed.len == 0) return error.ExpectedValue;

            const value = try parseValue(allocator, value_trimmed);
            try current_entries.append(allocator, .{ .key = key, .value = value });
            continue;
        }

        return error.UnexpectedCharacter;
    }

    // Flush the final section
    try flushSection(
        allocator,
        &tables,
        &array_tables,
        current_name,
        current_kind,
        &current_entries,
    );

    const owned_tables = try tables.toOwnedSlice(allocator);
    errdefer {
        for (owned_tables) |table| {
            for (table.entries) |entry| freeValue(allocator, entry.value);
            allocator.free(table.entries);
        }
        allocator.free(owned_tables);
    }
    const owned_array_tables = try array_tables.toOwnedSlice(allocator);

    return Document{
        .tables = owned_tables,
        .array_tables = owned_array_tables,
        .allocator = allocator,
    };
}

/// Flush accumulated entries into the appropriate list and reset for the
/// next section.
fn flushSection(
    allocator: Allocator,
    tables: *std.ArrayListUnmanaged(Table),
    array_tables: *std.ArrayListUnmanaged(Document.ArrayTable),
    name: []const u8,
    kind: SectionKind,
    entries: *std.ArrayListUnmanaged(Value.Entry),
) Allocator.Error!void {
    // Even if entries is empty, create the section so getTable("") works
    // for root when there are no top-level keys (edge case). However, to
    // avoid empty root tables cluttering the output, only emit a section
    // if there are entries or it is the first time for root.
    const is_empty_root_first_flush = (name.len == 0 and tables.items.len == 0);
    if (entries.items.len == 0 and !is_empty_root_first_flush) {
        // No entries and not the implicit root -- skip.
        return;
    }

    const owned_entries = try entries.toOwnedSlice(allocator);
    errdefer allocator.free(owned_entries);
    switch (kind) {
        .table => try tables.append(allocator, .{
            .name = name,
            .entries = owned_entries,
        }),
        .array_of_tables => try array_tables.append(allocator, .{
            .name = name,
            .entries = owned_entries,
        }),
    }
}

// ---------------------------------------------------------------------------
// Value parsing
// ---------------------------------------------------------------------------

/// Parse a single TOML value from the beginning of `raw`.
fn parseValue(allocator: Allocator, raw: []const u8) Error!Value {
    if (raw.len == 0) return error.ExpectedValue;

    // String
    if (raw[0] == '"') return parseString(raw);

    // Inline table
    if (raw[0] == '{') return parseInlineTable(allocator, raw);

    // Number (float or integer)
    return parseNumber(raw);
}

/// Parse a `"string"` value.  Handles `\"` and `\\` escapes.
/// Returns a slice into the original input (no allocation).
fn parseString(raw: []const u8) Error!Value {
    std.debug.assert(raw[0] == '"');
    var i: usize = 1;
    while (i < raw.len) {
        if (raw[i] == '\\') {
            // Skip the escaped character; bounds-check before advancing
            if (i + 1 >= raw.len) return error.UnterminatedString;
            i += 2;
            continue;
        }
        if (raw[i] == '"') {
            // Return the content between the quotes (excluding quotes).
            // Note: this returns the raw content including escape sequences.
            // For our use case (identifiers and simple strings), this is fine.
            return Value{ .string = raw[1..i] };
        }
        i += 1;
    }
    return error.UnterminatedString;
}

/// Parse an inline table `{ key = value, key = value }`.
fn parseInlineTable(allocator: Allocator, raw: []const u8) Error!Value {
    std.debug.assert(raw[0] == '{');

    // Find the closing '}'
    const close = std.mem.indexOfScalar(u8, raw, '}') orelse
        return error.InvalidInlineTable;

    const inner = std.mem.trim(u8, raw[1..close], " \t");
    if (inner.len == 0) {
        // Empty inline table
        const empty = try allocator.alloc(Value.Entry, 0);
        return Value{ .inline_table = empty };
    }

    var entries = std.ArrayListUnmanaged(Value.Entry){};
    errdefer {
        for (entries.items) |entry| {
            freeValue(allocator, entry.value);
        }
        entries.deinit(allocator);
    }

    // Split on ',' but we must be careful about strings containing commas.
    // For our use case, inline tables contain simple key = value pairs with
    // no nested inline tables or comma-in-string values, so a simple split
    // on ',' works.
    var pair_iter = std.mem.splitScalar(u8, inner, ',');
    while (pair_iter.next()) |pair_raw| {
        const pair = std.mem.trim(u8, pair_raw, " \t");
        if (pair.len == 0) continue;

        const eq = std.mem.indexOfScalar(u8, pair, '=') orelse
            return error.ExpectedEquals;

        const key = std.mem.trim(u8, pair[0..eq], " \t");
        const val_raw = std.mem.trim(u8, pair[eq + 1 ..], " \t");
        if (val_raw.len == 0) return error.ExpectedValue;

        const value = try parseValue(allocator, val_raw);
        try entries.append(allocator, .{ .key = key, .value = value });
    }

    return Value{ .inline_table = try entries.toOwnedSlice(allocator) };
}

/// Parse a number value (float or integer).
fn parseNumber(raw: []const u8) Error!Value {
    // Try float first (catches "1.87", "3.14159", etc.)
    if (std.fmt.parseFloat(f64, raw)) |f| {
        // Distinguish floats from integers: if the raw text contains a '.'
        // or 'e'/'E', treat it as float.
        if (std.mem.indexOfScalar(u8, raw, '.') != null or
            std.mem.indexOfScalar(u8, raw, 'e') != null or
            std.mem.indexOfScalar(u8, raw, 'E') != null)
        {
            return Value{ .float = f };
        }
        // Otherwise it parsed as float but is actually an integer
        if (std.fmt.parseInt(i64, raw, 10)) |i| {
            return Value{ .integer = i };
        } else |_| {
            // Fallback: return as float
            return Value{ .float = f };
        }
    } else |_| {}

    // Try integer
    if (std.fmt.parseInt(i64, raw, 10)) |i| {
        return Value{ .integer = i };
    } else |_| {}

    return error.InvalidNumber;
}

/// Strip a `#` comment from the end of a line, being careful not to strip
/// inside a quoted string.
fn stripComment(line: []const u8) []const u8 {
    var in_string = false;
    var i: usize = 0;
    while (i < line.len) {
        if (line[i] == '\\' and in_string) {
            // Skip escaped character inside string; bounds-check
            if (i + 1 >= line.len) break;
            i += 2;
            continue;
        }
        if (line[i] == '"') {
            in_string = !in_string;
        } else if (line[i] == '#' and !in_string) {
            return line[0..i];
        }
        i += 1;
    }
    return line;
}

// ===========================================================================
// Tests
// ===========================================================================

test "parse top-level string" {
    const input = "name = \"my-classifier\"";
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    try std.testing.expectEqualStrings("my-classifier", doc.getString("name").?);
}

test "parse table with inline tables" {
    const input =
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\O = { radius = 1.40, class = "polar" }
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    const types_table = doc.getTable("types").?;
    try std.testing.expectEqual(@as(usize, 2), types_table.entries.len);
}

test "parse array of tables" {
    const input =
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "CA"
        \\type = "C_ALI"
        \\
        \\[[atoms]]
        \\residue = "ALA"
        \\atom = "CB"
        \\type = "C_ALI"
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    try std.testing.expectEqual(@as(usize, 2), doc.array_tables.len);
}

test "parse comments and blank lines" {
    const input =
        \\# This is a comment
        \\name = "test"
        \\
        \\# Another comment
        \\[types]
        \\C = { radius = 1.70, class = "apolar" }
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    try std.testing.expectEqualStrings("test", doc.getString("name").?);
}

test "parse float and integer values" {
    const input =
        \\pi = 3.14159
        \\count = 42
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    const root = doc.getTable("").?;
    try std.testing.expectEqual(@as(usize, 2), root.entries.len);
}

test "parse error: unterminated string" {
    const input = "name = \"unterminated";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.UnterminatedString, result);
}

test "parse inline table with string and float" {
    const input =
        \\[types]
        \\C = { radius = 1.70, class = "apolar" }
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    const types = doc.getTable("types").?;
    try std.testing.expectEqual(@as(usize, 1), types.entries.len);
    switch (types.entries[0].value) {
        .inline_table => |entries| {
            try std.testing.expectEqual(@as(usize, 2), entries.len);
        },
        else => return error.TestUnexpectedResult,
    }
}

test "parse full classifier config" {
    const input =
        \\name = "test-classifier"
        \\
        \\[types]
        \\C_ALI = { radius = 1.87, class = "apolar" }
        \\O = { radius = 1.40, class = "polar" }
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "CA"
        \\type = "C_ALI"
        \\
        \\[[atoms]]
        \\residue = "ANY"
        \\atom = "O"
        \\type = "O"
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();

    try std.testing.expectEqualStrings("test-classifier", doc.getString("name").?);
    try std.testing.expect(doc.getTable("types") != null);
    try std.testing.expectEqual(@as(usize, 2), doc.array_tables.len);
}

test "parse error: invalid inline table" {
    const input = "x = { key = 1";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.InvalidInlineTable, result);
}

test "parse empty document" {
    const input = "# just a comment\n";
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    // Empty doc should have an empty root table
    try std.testing.expectEqual(@as(usize, 1), doc.tables.len);
}

test "parse error: expected value (empty rhs)" {
    const input = "key = ";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.ExpectedValue, result);
}

test "parse error: expected equals in inline table" {
    const input = "[types]\nC = { radius 1.87 }";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.ExpectedEquals, result);
}

test "parse error: unexpected character (bare word line)" {
    const input = "not a valid line";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.UnexpectedCharacter, result);
}

test "parse error: invalid number" {
    const input = "key = abc_not_number";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.InvalidNumber, result);
}

test "parse error: malformed [[ without ]]" {
    const input = "[[atoms\nresidue = \"ALA\"";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.UnexpectedCharacter, result);
}

test "parse error: trailing backslash in string" {
    const input = "key = \"trailing\\";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.UnterminatedString, result);
}

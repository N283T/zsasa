# TOML Classifier Config Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Allow users to define custom atom classifiers in TOML format via `--config=file.toml`, auto-detected by `.toml` extension.

**Architecture:** A minimal TOML subset parser (`toml_parser.zig`) handles generic TOML tokenization. A TOML classifier parser (`toml_classifier_parser.zig`) maps TOML data to the existing `classifier.Classifier` struct. The existing `parseConfigFile` function dispatches to the correct parser based on file extension.

**Tech Stack:** Zig 0.15.x, no external dependencies. TOML subset parser handles only what classifier configs need.

---

### Task 1: Minimal TOML Subset Parser

Build a generic TOML parser that handles: top-level key-value pairs, `[table]` sections, inline tables `{ k = v }`, `[[array_of_tables]]`, `"strings"`, floats, integers, and `#` comments.

**Files:**
- Create: `src/toml_parser.zig`

**Step 1: Create the TOML parser with data structures and tests**

Create `src/toml_parser.zig` with these types and functions:

```zig
const std = @import("std");
const Allocator = std.mem.Allocator;

pub const TomlError = error{
    UnterminatedString,
    InvalidNumber,
    InvalidInlineTable,
    ExpectedEquals,
    ExpectedValue,
    UnexpectedCharacter,
    OutOfMemory,
};

pub const Error = TomlError || Allocator.Error;

/// A TOML value
pub const Value = union(enum) {
    string: []const u8,     // "hello" -> slice into input
    float: f64,             // 1.87
    integer: i64,           // 42
    inline_table: []Entry,  // { radius = 1.87, class = "apolar" }

    pub const Entry = struct {
        key: []const u8,
        value: Value,
    };
};

/// A parsed TOML document section
pub const Table = struct {
    name: []const u8,       // "" for root, "types" for [types]
    entries: []const Value.Entry,
};

/// Result of parsing a TOML document
pub const Document = struct {
    /// Top-level tables: root (name=""), [types], etc.
    tables: []Table,
    /// Array-of-tables entries: [[atoms]] produces multiple entries
    array_tables: []ArrayTable,
    allocator: Allocator,

    pub const ArrayTable = struct {
        name: []const u8,       // "atoms"
        entries: []const Value.Entry,  // { residue = "ALA", atom = "CA", type = "C_ALI" }
    };

    pub fn deinit(self: *Document) void {
        // Free all allocated slices
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

    /// Look up a top-level string key (root table)
    pub fn getString(self: *const Document, key: []const u8) ?[]const u8 {
        for (self.tables) |table| {
            if (table.name.len == 0) { // root table
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

    /// Look up a named table
    pub fn getTable(self: *const Document, name: []const u8) ?Table {
        for (self.tables) |table| {
            if (std.mem.eql(u8, table.name, name)) return table;
        }
        return null;
    }

    /// Get all array-of-tables entries for a given name
    pub fn getArrayTables(self: *const Document, name: []const u8) []const ArrayTable {
        // Count matching entries
        var count: usize = 0;
        for (self.array_tables) |at| {
            if (std.mem.eql(u8, at.name, name)) count += 1;
        }
        // For simplicity, return the full slice - caller filters by name
        // (all array_tables in this parser are for a single name typically)
        _ = count;
        return self.array_tables;
    }
};

fn freeValue(allocator: Allocator, value: Value) void {
    switch (value) {
        .inline_table => |entries| {
            for (entries) |entry| {
                freeValue(allocator, entry.value);
            }
            allocator.free(entries);
        },
        else => {},
    }
}
```

The main `parse` function should:
1. Split input by lines
2. Skip empty lines and `#` comments
3. Detect `[table]` headers, `[[array_of_tables]]` headers, and `key = value` pairs
4. For values: parse `"string"`, numbers (float/integer), inline tables `{ ... }`
5. Group entries into tables and array tables
6. Return a `Document` that owns all allocated memory

Implementation approach for `pub fn parse(allocator: Allocator, content: []const u8) Error!Document`:

```zig
pub fn parse(allocator: Allocator, content: []const u8) Error!Document {
    var tables = std.ArrayList(Table).init(allocator);
    defer tables.deinit();
    var array_tables = std.ArrayList(Document.ArrayTable).init(allocator);
    defer array_tables.deinit();

    // Current table being built
    var current_name: []const u8 = ""; // root table
    var current_entries = std.ArrayList(Value.Entry).init(allocator);
    defer current_entries.deinit();

    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |raw_line| {
        const line = std.mem.trim(u8, stripComment(raw_line), " \t\r");
        if (line.len == 0) continue;

        // Check for [[array_of_tables]]
        if (std.mem.startsWith(u8, line, "[[") and std.mem.endsWith(u8, line, "]]")) {
            // Flush current table
            try flushTable(&tables, current_name, &current_entries, allocator);
            // Start new array table entry
            const aot_name = std.mem.trim(u8, line[2 .. line.len - 2], " \t");
            // Parse entries until next header
            // ... collect entries into current_entries, then flush as ArrayTable
            // Actually: set current context to array-of-tables mode
            // For simplicity: flush, then start collecting for this AoT entry
            try flushCurrentAsArrayTable(&array_tables, aot_name, &current_entries, allocator);
            current_name = aot_name; // mark we're in AoT mode
            // Actually need a flag to know if we're in AoT or table mode
            continue;
        }

        // Check for [table]
        if (line[0] == '[' and line[line.len - 1] == ']') {
            try flushTable(&tables, current_name, &current_entries, allocator);
            current_name = std.mem.trim(u8, line[1 .. line.len - 1], " \t");
            continue;
        }

        // Parse key = value
        if (std.mem.indexOfScalar(u8, line, '=')) |eq_idx| {
            const key = std.mem.trim(u8, line[0..eq_idx], " \t");
            const val_str = std.mem.trim(u8, line[eq_idx + 1 ..], " \t");
            const value = try parseValue(allocator, val_str);
            try current_entries.append(.{ .key = key, .value = value });
        }
    }
    // Flush last section
    // ... (flush to tables or array_tables depending on context)

    return Document{
        .tables = try tables.toOwnedSlice(),
        .array_tables = try array_tables.toOwnedSlice(),
        .allocator = allocator,
    };
}
```

The value parser should handle:

```zig
fn parseValue(allocator: Allocator, s: []const u8) Error!Value {
    if (s.len == 0) return error.ExpectedValue;

    // String: "..."
    if (s[0] == '"') {
        // Find closing quote (handle \" escapes)
        // Return .{ .string = content_between_quotes }
    }

    // Inline table: { ... }
    if (s[0] == '{') {
        // Parse comma-separated key = value pairs until }
        // Return .{ .inline_table = entries }
    }

    // Number: try float first, then integer
    if (std.fmt.parseFloat(f64, s)) |f| {
        return .{ .float = f };
    } else |_| {}
    if (std.fmt.parseInt(i64, s, 10)) |i| {
        return .{ .integer = i };
    } else |_| {}

    return error.ExpectedValue;
}
```

Key implementation details:
- String parsing: find matching `"`, handle `\"` and `\\` escapes within
- Inline table: tokenize on `,` between `{` and `}`, recursively parse each `key = value`
- Numbers: try `std.fmt.parseFloat` first (catches `1.87`), then `std.fmt.parseInt`
- All string slices point into the original `content` buffer (no allocation for strings)
- Only inline table `[]Entry` slices are heap-allocated
- Track whether current section is `[table]` or `[[array_of_tables]]` with an enum flag

**Tests to include in the file:**

```zig
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
    // Verify both values parsed
    const root = doc.getTable("").?;
    try std.testing.expectEqual(@as(usize, 2), root.entries.len);
}

test "parse error: unterminated string" {
    const input = "name = \"unterminated";
    const result = parse(std.testing.allocator, input);
    try std.testing.expectError(error.UnterminatedString, result);
}

test "parse inline table with string escape" {
    const input =
        \\[types]
        \\C = { radius = 1.70, class = "apolar" }
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    const types = doc.getTable("types").?;
    // C entry should have inline table with radius and class
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
```

**Step 2: Run tests**

Run: `zig build test 2>&1 | tail -20`
Expected: All existing tests still pass, new TOML parser tests pass.

**Step 3: Commit**

```bash
git add src/toml_parser.zig
git commit -m "feat: add minimal TOML subset parser (#158)"
```

---

### Task 2: TOML Classifier Parser

Convert parsed TOML documents into `classifier.Classifier` instances using the same data structures and error types as the existing FreeSASA parser.

**Files:**
- Create: `src/toml_classifier_parser.zig`

**Step 1: Create the TOML classifier parser with tests**

Create `src/toml_classifier_parser.zig`:

```zig
const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const Classifier = classifier.Classifier;
const AtomClass = classifier.AtomClass;
const toml_parser = @import("toml_parser.zig");
const Value = toml_parser.Value;

// Re-use error types from classifier_parser for consistency
const classifier_parser = @import("classifier_parser.zig");
pub const ParseError = classifier_parser.ParseError;
pub const Error = ParseError || Allocator.Error || toml_parser.Error;

/// Parse a TOML classifier config from string content.
///
/// Expected TOML schema:
/// ```toml
/// name = "classifier-name"
///
/// [types]
/// TYPE_NAME = { radius = 1.87, class = "apolar" }
///
/// [[atoms]]
/// residue = "ALA"
/// atom = "CA"
/// type = "TYPE_NAME"
/// ```
///
/// Returns a Classifier that must be freed with `deinit()`.
pub fn parseConfig(allocator: Allocator, content: []const u8) Error!Classifier {
    var doc = try toml_parser.parse(allocator, content);
    defer doc.deinit();

    // Get classifier name (default: "custom")
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

    // Create classifier
    var result = try Classifier.init(allocator, name);
    errdefer result.deinit();

    // Parse [[atoms]] entries
    for (doc.array_tables) |at| {
        if (!std.mem.eql(u8, at.name, "atoms")) continue;

        const residue = getStringFromEntries(at.entries, "residue") orelse return error.InvalidAtomDefinition;
        const atom = getStringFromEntries(at.entries, "atom") orelse return error.InvalidAtomDefinition;
        const type_name = getStringFromEntries(at.entries, "type") orelse return error.InvalidAtomDefinition;

        const type_def = types.get(type_name) orelse return error.UndefinedType;
        try result.addAtom(residue, atom, type_def.radius, type_def.class);
    }

    return result;
}

const TypeDef = struct {
    radius: f64,
    class: AtomClass,
};

fn parseClass(s: []const u8) ?AtomClass {
    if (std.mem.eql(u8, s, "polar")) return .polar;
    if (std.mem.eql(u8, s, "apolar")) return .apolar;
    return null;
}

/// Get a float value from inline table entries
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

/// Get a string value from inline table entries
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

/// Get a string value from a flat entries list (for [[atoms]] sections)
fn getStringFromEntries(entries: []const Value.Entry, key: []const u8) ?[]const u8 {
    return getString(entries, key);
}
```

Tests to include:

```zig
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

    // CYS has specific entry
    try std.testing.expectEqual(@as(?f64, 1.76), result.getRadius("CYS", "CA"));
    // Other residues use ANY fallback
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
    // Missing "type" field

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
```

**Step 2: Run tests**

Run: `zig build test 2>&1 | tail -20`
Expected: All tests pass.

**Step 3: Commit**

```bash
git add src/toml_classifier_parser.zig
git commit -m "feat: add TOML classifier parser (#158)"
```

---

### Task 3: Wire Up Format Detection

Modify `parseConfigFile` in `classifier_parser.zig` to detect `.toml` extension and dispatch to the TOML parser.

**Files:**
- Modify: `src/classifier_parser.zig:256-279` (the `parseConfigFile` function)

**Step 1: Add the TOML dispatch**

At the top of `classifier_parser.zig`, add the import:

```zig
const toml_classifier_parser = @import("toml_classifier_parser.zig");
```

Then modify `parseConfigFile` to detect `.toml`:

```zig
/// Parse a classifier config file from disk.
///
/// Format is auto-detected by file extension:
/// - `.toml` files are parsed as TOML
/// - All other extensions use the FreeSASA config format
pub fn parseConfigFile(allocator: Allocator, path: []const u8) Error!Classifier {
    // Auto-detect format by extension
    if (std.mem.endsWith(u8, path, ".toml")) {
        return parseTomlConfigFile(allocator, path);
    }

    // Default: FreeSASA format
    // ... (existing code stays the same)
}

/// Parse a TOML classifier config file from disk.
fn parseTomlConfigFile(allocator: Allocator, path: []const u8) Error!Classifier {
    const file = std.fs.cwd().openFile(path, .{}) catch {
        return error.FileReadError;
    };
    defer file.close();

    const stat = file.stat() catch {
        return error.FileReadError;
    };

    const content = allocator.alloc(u8, stat.size) catch {
        return error.OutOfMemory;
    };
    defer allocator.free(content);

    const bytes_read = file.readAll(content) catch {
        return error.FileReadError;
    };

    return toml_classifier_parser.parseConfig(allocator, content[0..bytes_read]);
}
```

Note: The `Error` type union in `classifier_parser.zig` may need to include `toml_parser.TomlError`. Check if the error sets are compatible. If not, widen the `Error` type:

```zig
pub const Error = ParseError || Allocator.Error || std.fs.File.OpenError || std.fs.File.ReadError || toml_parser.TomlError;
```

**Step 2: Run tests**

Run: `zig build test 2>&1 | tail -20`
Expected: All tests pass. No behavioral change for existing `.config` files.

**Step 3: Commit**

```bash
git add src/classifier_parser.zig
git commit -m "feat: auto-detect TOML classifier by .toml extension (#158)"
```

---

### Task 4: Round-Trip Equivalence Test

Verify that equivalent FreeSASA and TOML configs produce identical `Classifier` output.

**Files:**
- Modify: `src/toml_classifier_parser.zig` (add test at bottom)

**Step 1: Add the round-trip test**

Add this test to `src/toml_classifier_parser.zig`:

```zig
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
    const test_cases = [_][2][]const u8{
        .{ "ALA", "C" },  .{ "ALA", "O" },  .{ "ALA", "CA" },
        .{ "ALA", "N" },  .{ "ALA", "CB" }, .{ "CYS", "SG" },
        .{ "GLY", "C" },  .{ "GLY", "CA" }, // ANY fallback
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
```

**Step 2: Run tests**

Run: `zig build test 2>&1 | tail -20`
Expected: Round-trip test passes.

**Step 3: Commit**

```bash
git add src/toml_classifier_parser.zig
git commit -m "test: add round-trip equivalence test for TOML vs FreeSASA (#158)"
```

---

### Task 5: Register New Modules and Update Documentation

Add the new modules to `root.zig`, update `docs/cli.md` with TOML format documentation, and add a CHANGELOG entry.

**Files:**
- Modify: `src/root.zig` - Add `toml_parser` and `toml_classifier_parser` exports
- Modify: `docs/cli.md:372-385` - Add TOML config format documentation
- Modify: `CHANGELOG.md:10-18` - Add entry

**Step 1: Update root.zig**

Add after the existing `pub const stream_writer` line:

```zig
pub const toml_parser = @import("toml_parser.zig");
pub const toml_classifier_parser = @import("toml_classifier_parser.zig");
```

And in the module doc comment, add:

```zig
//! - `toml_parser` - Minimal TOML subset parser
//! - `toml_classifier_parser` - TOML-to-Classifier converter
```

And in the test block:

```zig
_ = toml_parser;
_ = toml_classifier_parser;
```

**Step 2: Update docs/cli.md**

Replace the "Custom Config" section (around line 372-385) with:

```markdown
### Custom Config (`--config=FILE`)

Load a custom classifier from a config file. The format is auto-detected by extension:

- `.toml` files use [TOML format](#toml-format)
- All other extensions use [FreeSASA format](#freesasa-format)

#### TOML Format

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }
C_CAR = { radius = 1.76, class = "apolar" }
N     = { radius = 1.65, class = "polar" }
O     = { radius = 1.40, class = "polar" }
S     = { radius = 1.85, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"

[[atoms]]
residue = "ALA"
atom = "CB"
type = "C_ALI"
```

- `name` - Classifier name (optional, default: "custom")
- `[types]` - Define atom types with radius (Å) and class (`"polar"` or `"apolar"`)
- `[[atoms]]` - Map (residue, atom) pairs to defined types. Use `"ANY"` for fallback entries.

#### FreeSASA Format

```
name: my-classifier

types:
C_ALI 1.87 apolar
C_CAR 1.76 apolar
O     1.40 polar

atoms:
ANY CA  C_ALI
ALA CB  C_ALI
```
```

**Step 3: Update CHANGELOG.md**

Add under `### Added` in the `[Unreleased]` section:

```markdown
- TOML format support for custom classifier configs (`--config=file.toml`): human-friendly alternative to FreeSASA format with auto-detection by file extension (#158)
```

**Step 4: Run tests**

Run: `zig build test 2>&1 | tail -20`
Expected: All tests pass.

**Step 5: Commit**

```bash
git add src/root.zig docs/cli.md CHANGELOG.md
git commit -m "docs: add TOML classifier documentation and changelog (#158)"
```

---

## Summary

| Task | What | New/Modified |
|------|------|--------------|
| 1 | Minimal TOML subset parser | Create `src/toml_parser.zig` |
| 2 | TOML classifier parser | Create `src/toml_classifier_parser.zig` |
| 3 | Format detection dispatch | Modify `src/classifier_parser.zig` |
| 4 | Round-trip equivalence test | Modify `src/toml_classifier_parser.zig` |
| 5 | Module registration + docs | Modify `src/root.zig`, `docs/cli.md`, `CHANGELOG.md` |

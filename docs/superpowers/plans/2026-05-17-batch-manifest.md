# Batch TOML Manifest Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add TOML manifest support to `zsasa batch` so a single command can run named A/B/AB chain-selection batch jobs with job-separated output.

**Architecture:** Extend the existing minimal TOML parser only enough for this manifest shape, add a focused `src/batch_manifest.zig` parser/resolver, then wire resolved manifest jobs into the existing `src/batch.zig` execution path. Existing non-manifest batch behavior remains the default path.

**Tech Stack:** Zig 0.16, existing `std.Io` file APIs, existing `src/toml_parser.zig`, existing `src/batch.zig`, `zig build test`, `zig build run` smoke tests.

---

## File Structure

- Modify `/Users/nagaet/zsasa/src/toml_parser.zig`
  - Add boolean values and arrays of strings so TOML manifest files can express `use_bitmask = true` and `chains = ["A", "B"]`.
- Create `/Users/nagaet/zsasa/src/batch_manifest.zig`
  - Own all manifest parsing, validation, path-safe job name checks, and conversion of TOML values into manifest data structures.
- Modify `/Users/nagaet/zsasa/src/batch.zig`
  - Add `--manifest`, non-manifest `--chain`, `--auth-chain`, explicit CLI presence flags, chain-filter propagation, and manifest job execution.
- Modify `/Users/nagaet/zsasa/src/root.zig`
  - Export/import `batch_manifest` for tests and compilation coverage.
- Modify `/Users/nagaet/zsasa/website/docs/cli/commands.md`
  - Document `--manifest`, TOML example, precedence rules, and output layout.
- Optionally modify `/Users/nagaet/zsasa/CHANGELOG.md`
  - Add an unreleased entry if this repository expects every user-facing CLI feature to be listed immediately.

---

### Task 1: Extend TOML Parser for Manifest Values

**Files:**
- Modify: `/Users/nagaet/zsasa/src/toml_parser.zig`

- [ ] **Step 1: Add failing tests for booleans and string arrays**

Append these tests near the existing TOML parser tests in `/Users/nagaet/zsasa/src/toml_parser.zig`:

```zig
test "parse boolean values" {
    const input =
        \\use_bitmask = true
        \\quiet = false
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    const root = doc.getTable("").?;
    try std.testing.expectEqual(@as(usize, 2), root.entries.len);
    switch (root.entries[0].value) {
        .boolean => |b| try std.testing.expectEqual(true, b),
        else => return error.TestUnexpectedResult,
    }
    switch (root.entries[1].value) {
        .boolean => |b| try std.testing.expectEqual(false, b),
        else => return error.TestUnexpectedResult,
    }
}

test "parse string array values" {
    const input =
        \\chains = ["A", "B"]
        \\sdf = ["lig1.sdf", "lig2.sdf"]
    ;
    var doc = try parse(std.testing.allocator, input);
    defer doc.deinit();
    const root = doc.getTable("").?;
    try std.testing.expectEqual(@as(usize, 2), root.entries.len);
    switch (root.entries[0].value) {
        .string_array => |items| {
            try std.testing.expectEqual(@as(usize, 2), items.len);
            try std.testing.expectEqualStrings("A", items[0]);
            try std.testing.expectEqualStrings("B", items[1]);
        },
        else => return error.TestUnexpectedResult,
    }
    switch (root.entries[1].value) {
        .string_array => |items| {
            try std.testing.expectEqual(@as(usize, 2), items.len);
            try std.testing.expectEqualStrings("lig1.sdf", items[0]);
            try std.testing.expectEqualStrings("lig2.sdf", items[1]);
        },
        else => return error.TestUnexpectedResult,
    }
}
```

- [ ] **Step 2: Run the focused test and verify it fails**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task1-fail.log
```

Expected: FAIL with compile errors referencing missing `boolean` and `string_array` union fields.

- [ ] **Step 3: Add boolean and string-array TOML values**

In `/Users/nagaet/zsasa/src/toml_parser.zig`, update `Value`:

```zig
pub const Value = union(enum) {
    string: []const u8,
    float: f64,
    integer: i64,
    boolean: bool,
    string_array: []const []const u8,
    inline_table: []const Entry,

    pub const Entry = struct {
        key: []const u8,
        value: Value,
    };
};
```

Update `freeValue`:

```zig
fn freeValue(allocator: Allocator, value: Value) void {
    switch (value) {
        .inline_table => |entries| {
            for (entries) |entry| {
                freeValue(allocator, entry.value);
            }
            allocator.free(entries);
        },
        .string_array => |items| allocator.free(items),
        .string, .float, .integer, .boolean => {},
    }
}
```

Update `parseValue` before number parsing:

```zig
fn parseValue(allocator: Allocator, raw: []const u8) Error!Value {
    if (raw.len == 0) return error.ExpectedValue;

    if (raw[0] == '"') return parseString(raw);
    if (raw[0] == '{') return parseInlineTable(allocator, raw);
    if (raw[0] == '[') return parseStringArray(allocator, raw);
    if (std.mem.eql(u8, raw, "true")) return Value{ .boolean = true };
    if (std.mem.eql(u8, raw, "false")) return Value{ .boolean = false };

    return parseNumber(raw);
}
```

Add this helper near `parseInlineTable`:

```zig
fn parseStringArray(allocator: Allocator, raw: []const u8) Error!Value {
    std.debug.assert(raw[0] == '[');
    const close = std.mem.findScalar(u8, raw, ']') orelse return error.UnexpectedCharacter;
    const trailing = std.mem.trim(u8, raw[close + 1 ..], " \t");
    if (trailing.len != 0) return error.UnexpectedCharacter;

    const inner = std.mem.trim(u8, raw[1..close], " \t");
    var items = std.ArrayListUnmanaged([]const u8).empty;
    errdefer items.deinit(allocator);

    if (inner.len != 0) {
        var iter = std.mem.splitScalar(u8, inner, ',');
        while (iter.next()) |part_raw| {
            const part = std.mem.trim(u8, part_raw, " \t");
            if (part.len == 0) return error.ExpectedValue;
            const parsed = try parseString(part);
            switch (parsed) {
                .string => |s| try items.append(allocator, s),
                else => unreachable,
            }
        }
    }

    return Value{ .string_array = try items.toOwnedSlice(allocator) };
}
```

- [ ] **Step 4: Run the TOML tests and verify they pass**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task1-pass.log
```

Expected: PASS for the full Zig test suite.

- [ ] **Step 5: Commit Task 1**

Run:

```bash
git add /Users/nagaet/zsasa/src/toml_parser.zig
git commit -m "feat: parse manifest TOML values"
```

---

### Task 2: Add Batch Manifest Parser

**Files:**
- Create: `/Users/nagaet/zsasa/src/batch_manifest.zig`
- Modify: `/Users/nagaet/zsasa/src/root.zig`

- [ ] **Step 1: Create failing parser tests and public types**

Create `/Users/nagaet/zsasa/src/batch_manifest.zig` with this initial content:

```zig
const std = @import("std");
const Allocator = std.mem.Allocator;
const toml_parser = @import("toml_parser.zig");

pub const ManifestError = error{
    UnsupportedVersion,
    MissingJobName,
    DuplicateJobName,
    UnsafeJobName,
    InvalidFieldType,
    NoJobs,
};

pub const Error = ManifestError || Allocator.Error || toml_parser.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError;

pub const Globals = struct {
    input_dir: ?[]const u8 = null,
    output_dir: ?[]const u8 = null,
    algorithm: ?[]const u8 = null,
    classifier: ?[]const u8 = null,
    ccd: ?[]const u8 = null,
    sdf: ?[]const []const u8 = null,
    threads: ?usize = null,
    probe_radius: ?f64 = null,
    n_points: ?u32 = null,
    n_slices: ?u32 = null,
    precision: ?[]const u8 = null,
    format: ?[]const u8 = null,
    include_hydrogens: ?bool = null,
    include_hetatm: ?bool = null,
    use_bitmask: ?bool = null,
    timing: ?bool = null,
    quiet: ?bool = null,
    auth_chain: ?bool = null,
};

pub const Job = struct {
    name: []const u8,
    chains: ?[]const []const u8 = null,
    auth_chain: ?bool = null,
};

pub const Manifest = struct {
    allocator: Allocator,
    content: []const u8,
    globals: Globals,
    jobs: []Job,

    pub fn deinit(self: *Manifest) void {
        if (self.globals.sdf) |items| self.allocator.free(items);
        for (self.jobs) |job| {
            if (job.chains) |chains| self.allocator.free(chains);
        }
        self.allocator.free(self.jobs);
        self.allocator.free(self.content);
        self.* = undefined;
    }
};

pub fn parse(allocator: Allocator, content: []const u8) Error!Manifest {
    _ = allocator;
    _ = content;
    return error.NoJobs;
}

pub fn parseFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Manifest {
    const f = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer f.close(io);
    var read_buf: [65536]u8 = undefined;
    var reader = f.reader(io, &read_buf);
    const content = try reader.interface.allocRemaining(allocator, .unlimited);
    errdefer allocator.free(content);
    return parseOwned(allocator, content);
}

fn parseOwned(allocator: Allocator, content: []const u8) Error!Manifest {
    _ = allocator;
    _ = content;
    return error.NoJobs;
}

test "parse A B AB manifest" {
    const input =
        \\version = 1
        \\input_dir = "structures"
        \\output_dir = "results"
        \\format = "jsonl"
        \\use_bitmask = true
        \\n_points = 128
        \\classifier = "ccd"
        \\
        \\[[jobs]]
        \\name = "chain_A"
        \\chains = ["A"]
        \\
        \\[[jobs]]
        \\name = "chain_B"
        \\chains = ["B"]
        \\
        \\[[jobs]]
        \\name = "complex_AB"
        \\chains = ["A", "B"]
    ;
    var manifest = try parse(std.testing.allocator, input);
    defer manifest.deinit();
    try std.testing.expectEqualStrings("structures", manifest.globals.input_dir.?);
    try std.testing.expectEqualStrings("results", manifest.globals.output_dir.?);
    try std.testing.expectEqualStrings("jsonl", manifest.globals.format.?);
    try std.testing.expectEqual(true, manifest.globals.use_bitmask.?);
    try std.testing.expectEqual(@as(u32, 128), manifest.globals.n_points.?);
    try std.testing.expectEqual(@as(usize, 3), manifest.jobs.len);
    try std.testing.expectEqualStrings("chain_A", manifest.jobs[0].name);
    try std.testing.expectEqualStrings("A", manifest.jobs[0].chains.?[0]);
    try std.testing.expectEqualStrings("complex_AB", manifest.jobs[2].name);
    try std.testing.expectEqual(@as(usize, 2), manifest.jobs[2].chains.?.len);
}

test "reject duplicate and unsafe job names" {
    const duplicate =
        \\version = 1
        \\[[jobs]]
        \\name = "dup"
        \\[[jobs]]
        \\name = "dup"
    ;
    try std.testing.expectError(error.DuplicateJobName, parse(std.testing.allocator, duplicate));

    const unsafe =
        \\version = 1
        \\[[jobs]]
        \\name = "../escape"
    ;
    try std.testing.expectError(error.UnsafeJobName, parse(std.testing.allocator, unsafe));
}
```

- [ ] **Step 2: Export the module from root**

Add this to `/Users/nagaet/zsasa/src/root.zig` near the other public imports:

```zig
pub const batch_manifest = @import("batch_manifest.zig");
```

Add this to the root test block near other `_ = ...;` lines:

```zig
    _ = batch_manifest;
```

- [ ] **Step 3: Run tests and verify parser tests fail**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task2-fail.log
```

Expected: FAIL in `batch_manifest` tests with `error.NoJobs`.

- [ ] **Step 4: Implement manifest parsing helpers**

Replace the initial `parse` and `parseOwned` bodies in `/Users/nagaet/zsasa/src/batch_manifest.zig` with these helpers and implementation:

```zig
pub fn parse(allocator: Allocator, content: []const u8) Error!Manifest {
    const owned = try allocator.dupe(u8, content);
    errdefer allocator.free(owned);
    return parseOwned(allocator, owned);
}

fn parseOwned(allocator: Allocator, content: []const u8) Error!Manifest {
    var doc = try toml_parser.parse(allocator, content);
    defer doc.deinit();

    const root = doc.getTable("") orelse return error.UnsupportedVersion;
    const version = getInteger(root.entries, "version") orelse return error.UnsupportedVersion;
    if (version != 1) return error.UnsupportedVersion;

    var globals = Globals{};
    globals.input_dir = getString(root.entries, "input_dir");
    globals.output_dir = getString(root.entries, "output_dir");
    globals.algorithm = getString(root.entries, "algorithm");
    globals.classifier = getString(root.entries, "classifier");
    globals.ccd = getString(root.entries, "ccd");
    globals.threads = if (getInteger(root.entries, "threads")) |v| @intCast(v) else null;
    globals.probe_radius = getFloat(root.entries, "probe_radius");
    globals.n_points = if (getInteger(root.entries, "n_points")) |v| @intCast(v) else null;
    globals.n_slices = if (getInteger(root.entries, "n_slices")) |v| @intCast(v) else null;
    globals.precision = getString(root.entries, "precision");
    globals.format = getString(root.entries, "format");
    globals.include_hydrogens = getBool(root.entries, "include_hydrogens");
    globals.include_hetatm = getBool(root.entries, "include_hetatm");
    globals.use_bitmask = getBool(root.entries, "use_bitmask");
    globals.timing = getBool(root.entries, "timing");
    globals.quiet = getBool(root.entries, "quiet");
    globals.auth_chain = getBool(root.entries, "auth_chain");
    if (getStringArray(root.entries, "sdf")) |items| {
        globals.sdf = try allocator.dupe([]const u8, items);
    } else if (getString(root.entries, "sdf")) |one| {
        const items = try allocator.alloc([]const u8, 1);
        items[0] = one;
        globals.sdf = items;
    }
    errdefer if (globals.sdf) |items| allocator.free(items);

    var jobs = std.ArrayListUnmanaged(Job).empty;
    errdefer {
        for (jobs.items) |job| if (job.chains) |chains| allocator.free(chains);
        jobs.deinit(allocator);
    }

    for (doc.array_tables) |at| {
        if (!std.mem.eql(u8, at.name, "jobs")) continue;
        const name = getString(at.entries, "name") orelse return error.MissingJobName;
        if (!isSafeJobName(name)) return error.UnsafeJobName;
        for (jobs.items) |existing| {
            if (std.mem.eql(u8, existing.name, name)) return error.DuplicateJobName;
        }
        var chains_copy: ?[]const []const u8 = null;
        if (getStringArray(at.entries, "chains")) |chains| {
            chains_copy = try allocator.dupe([]const u8, chains);
        }
        errdefer if (chains_copy) |chains| allocator.free(chains);
        try jobs.append(allocator, .{
            .name = name,
            .chains = chains_copy,
            .auth_chain = getBool(at.entries, "auth_chain"),
        });
    }

    if (jobs.items.len == 0) return error.NoJobs;

    return .{
        .allocator = allocator,
        .content = content,
        .globals = globals,
        .jobs = try jobs.toOwnedSlice(allocator),
    };
}

fn getString(entries: []const toml_parser.Value.Entry, key: []const u8) ?[]const u8 {
    for (entries) |entry| if (std.mem.eql(u8, entry.key, key)) {
        return switch (entry.value) { .string => |s| s, else => null };
    };
    return null;
}

fn getStringArray(entries: []const toml_parser.Value.Entry, key: []const u8) ?[]const []const u8 {
    for (entries) |entry| if (std.mem.eql(u8, entry.key, key)) {
        return switch (entry.value) { .string_array => |items| items, else => null };
    };
    return null;
}

fn getBool(entries: []const toml_parser.Value.Entry, key: []const u8) ?bool {
    for (entries) |entry| if (std.mem.eql(u8, entry.key, key)) {
        return switch (entry.value) { .boolean => |b| b, else => null };
    };
    return null;
}

fn getInteger(entries: []const toml_parser.Value.Entry, key: []const u8) ?i64 {
    for (entries) |entry| if (std.mem.eql(u8, entry.key, key)) {
        return switch (entry.value) { .integer => |i| i, else => null };
    };
    return null;
}

fn getFloat(entries: []const toml_parser.Value.Entry, key: []const u8) ?f64 {
    for (entries) |entry| if (std.mem.eql(u8, entry.key, key)) {
        return switch (entry.value) {
            .float => |f| f,
            .integer => |i| @floatFromInt(i),
            else => null,
        };
    };
    return null;
}

fn isSafeJobName(name: []const u8) bool {
    if (name.len == 0) return false;
    if (std.mem.indexOfScalar(u8, name, '/') != null) return false;
    if (std.mem.indexOfScalar(u8, name, '\\') != null) return false;
    var parts = std.mem.splitScalar(u8, name, '.');
    _ = parts;
    return std.mem.indexOf(u8, name, "..") == null;
}
```

- [ ] **Step 5: Run parser tests and verify they pass**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task2-pass.log
```

Expected: PASS for the full Zig test suite.

- [ ] **Step 6: Commit Task 2**

Run:

```bash
git add /Users/nagaet/zsasa/src/batch_manifest.zig /Users/nagaet/zsasa/src/root.zig
git commit -m "feat: parse batch manifests"
```

---

### Task 3: Add Batch CLI State for Manifest and Chain Filtering

**Files:**
- Modify: `/Users/nagaet/zsasa/src/batch.zig`

- [ ] **Step 1: Add failing argument parser tests**

Append these tests near the existing `BatchArgs` tests in `/Users/nagaet/zsasa/src/batch.zig`:

```zig
test "BatchArgs --manifest" {
    const args = [_][]const u8{ "zsasa", "batch", "--manifest", "bsa.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("bsa.toml", parsed.manifest_path.?);
}

test "BatchArgs --chain=A" {
    const args = [_][]const u8{ "zsasa", "batch", "--chain=A", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("A", parsed.chain_filter.?);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
}

test "BatchArgs --auth-chain" {
    const args = [_][]const u8{ "zsasa", "batch", "--auth-chain", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_auth_chain);
}

test "BatchArgs explicit option flags" {
    const args = [_][]const u8{
        "zsasa", "batch", "--threads=8", "--n-points=128", "--format=jsonl", "--use-bitmask", "input_dir/",
    };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.threads_explicit);
    try std.testing.expectEqual(true, parsed.n_points_explicit);
    try std.testing.expectEqual(true, parsed.format_explicit);
    try std.testing.expectEqual(true, parsed.use_bitmask_explicit);
}
```

- [ ] **Step 2: Run tests and verify they fail**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task3-fail.log
```

Expected: FAIL with missing `BatchArgs` fields such as `manifest_path` and `chain_filter`.

- [ ] **Step 3: Extend imports and `BatchArgs` fields**

At the top of `/Users/nagaet/zsasa/src/batch.zig`, add:

```zig
const batch_manifest = @import("batch_manifest.zig");
```

Extend `BatchArgs` with these fields:

```zig
    manifest_path: ?[]const u8 = null,
    chain_filter: ?[]const u8 = null,
    use_auth_chain: bool = false,
    threads_explicit: bool = false,
    probe_radius_explicit: bool = false,
    n_points_explicit: bool = false,
    n_slices_explicit: bool = false,
    algorithm_explicit: bool = false,
    precision_explicit: bool = false,
    format_explicit: bool = false,
    classifier_explicit: bool = false,
    include_hydrogens_explicit: bool = false,
    include_hetatm_explicit: bool = false,
    use_bitmask_explicit: bool = false,
    ccd_explicit: bool = false,
    sdf_explicit: bool = false,
    quiet_explicit: bool = false,
    timing_explicit: bool = false,
```

- [ ] **Step 4: Update `parseArgs` for manifest, chain, auth-chain, and presence flags**

In each existing option branch in `parseArgs`, set the matching explicit flag when the user provided that option. For example:

```zig
if (std.mem.startsWith(u8, arg, "--threads=")) {
    const value = arg["--threads=".len..];
    result.n_threads = std.fmt.parseInt(usize, value, 10) catch {
        std.debug.print("Error: Invalid thread count: {s}\n", .{value});
        std.process.exit(1);
    };
    result.threads_explicit = true;
} else if (std.mem.eql(u8, arg, "--threads")) {
    i += 1;
    if (i >= args.len) {
        std.debug.print("Error: Missing value for --threads\n", .{});
        std.process.exit(1);
    }
    result.n_threads = std.fmt.parseInt(usize, args[i], 10) catch {
        std.debug.print("Error: Invalid thread count: {s}\n", .{args[i]});
        std.process.exit(1);
    };
    result.threads_explicit = true;
}
```

Add these new option branches before the unknown-option branch:

```zig
        // --manifest=PATH or --manifest PATH
        else if (std.mem.startsWith(u8, arg, "--manifest=")) {
            result.manifest_path = arg["--manifest=".len..];
        } else if (std.mem.eql(u8, arg, "--manifest")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --manifest\n", .{});
                std.process.exit(1);
            }
            result.manifest_path = args[i];
        }
        // --chain=ID or --chain ID
        else if (std.mem.startsWith(u8, arg, "--chain=")) {
            result.chain_filter = arg["--chain=".len..];
        } else if (std.mem.eql(u8, arg, "--chain")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --chain\n", .{});
                std.process.exit(1);
            }
            result.chain_filter = args[i];
        }
        // --auth-chain
        else if (std.mem.eql(u8, arg, "--auth-chain")) {
            result.use_auth_chain = true;
        }
```

Set the matching explicit flags in existing branches:

```zig
result.probe_radius_explicit = true;
result.n_points_explicit = true;
result.n_slices_explicit = true;
result.format_explicit = true;
result.algorithm_explicit = true;
result.classifier_explicit = true;
result.precision_explicit = true;
result.include_hydrogens_explicit = true;
result.include_hetatm_explicit = true;
result.use_bitmask_explicit = true;
result.ccd_explicit = true;
result.sdf_explicit = true;
result.quiet_explicit = true;
result.timing_explicit = true;
```

- [ ] **Step 5: Run parser tests and verify they pass**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task3-pass.log
```

Expected: PASS for the full Zig test suite.

- [ ] **Step 6: Commit Task 3**

Run:

```bash
git add /Users/nagaet/zsasa/src/batch.zig
git commit -m "feat: parse batch manifest options"
```

---

### Task 4: Add Chain Filtering to Existing Batch Execution

**Files:**
- Modify: `/Users/nagaet/zsasa/src/batch.zig`

- [ ] **Step 1: Add a failing non-manifest chain filtering test**

Add this test near the `BatchArgs --chain=A` test in `/Users/nagaet/zsasa/src/batch.zig`:

```zig
test "parseBatchChainFilter splits comma-separated chains" {
    const chains = try parseBatchChainFilter(std.testing.allocator, "A, B,AB");
    defer std.testing.allocator.free(chains);
    try std.testing.expectEqual(@as(usize, 3), chains.len);
    try std.testing.expectEqualStrings("A", chains[0]);
    try std.testing.expectEqualStrings("B", chains[1]);
    try std.testing.expectEqualStrings("AB", chains[2]);
}
```

- [ ] **Step 2: Run tests and verify they fail**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task4-fail.log
```

Expected: FAIL with missing `parseBatchChainFilter`.

- [ ] **Step 3: Add `BatchConfig` chain fields and parser helper**

Extend `BatchConfig` in `/Users/nagaet/zsasa/src/batch.zig`:

```zig
    chain_filter: ?[]const []const u8 = null,
    use_auth_chain: bool = false,
```

Add this helper near the parse helper functions:

```zig
fn parseBatchChainFilter(allocator: Allocator, filter_str: []const u8) ![]const []const u8 {
    var chains = std.ArrayListUnmanaged([]const u8).empty;
    errdefer chains.deinit(allocator);

    var iter = std.mem.splitScalar(u8, filter_str, ',');
    while (iter.next()) |chain| {
        const trimmed = std.mem.trim(u8, chain, " ");
        if (trimmed.len > 0) {
            try chains.append(allocator, trimmed);
        }
    }

    return chains.toOwnedSlice(allocator);
}
```

Update `readInputFile` parser setup:

```zig
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            parser.chain_filter = config.chain_filter;
            parser.use_auth_chain = config.use_auth_chain;
            break :blk parser.parseFile(io, path);
        },
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            parser.chain_filter = config.chain_filter;
            break :blk parser.parseFile(io, path);
        },
```

Update the non-manifest `config` construction in `run` to parse the CLI chain filter:

```zig
    var chain_filter_slice: ?[]const []const u8 = null;
    if (args.chain_filter) |filter_str| {
        chain_filter_slice = try parseBatchChainFilter(allocator, filter_str);
    }
    defer if (chain_filter_slice) |s| allocator.free(s);
```

Set these `BatchConfig` fields:

```zig
        .chain_filter = chain_filter_slice,
        .use_auth_chain = args.use_auth_chain,
```

- [ ] **Step 4: Run tests and verify they pass**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task4-pass.log
```

Expected: PASS for the full Zig test suite.

- [ ] **Step 5: Smoke check current issue symptom is fixed for non-manifest mode**

Run:

```bash
zig build run -- batch --chain=A --help 2>&1 | tee /tmp/zsasa-task4-help.log
```

Expected: the command prints batch help or exits through the help path, not `Unknown option: --chain=A`.

- [ ] **Step 6: Commit Task 4**

Run:

```bash
git add /Users/nagaet/zsasa/src/batch.zig
git commit -m "feat: filter batch inputs by chain"
```

---

### Task 5: Resolve and Run Manifest Jobs

**Files:**
- Modify: `/Users/nagaet/zsasa/src/batch.zig`

- [ ] **Step 1: Add failing manifest conflict and output tests**

Append these tests near the batch tests in `/Users/nagaet/zsasa/src/batch.zig`:

```zig
test "manifestJsonlOutputPath uses job file under output dir" {
    const path = try manifestJsonlOutputPath(std.testing.allocator, "results", "chain_A");
    defer std.testing.allocator.free(path);
    try std.testing.expectEqualStrings("results/chain_A.jsonl", path);
}

test "manifestPerFileOutputDir uses job directory under output dir" {
    const path = try manifestPerFileOutputDir(std.testing.allocator, "results", "complex_AB");
    defer std.testing.allocator.free(path);
    try std.testing.expectEqualStrings("results/complex_AB", path);
}
```

- [ ] **Step 2: Run tests and verify they fail**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task5-fail.log
```

Expected: FAIL with missing output helper functions.

- [ ] **Step 3: Add output path helpers**

Add these helpers near `getOutputExtension` in `/Users/nagaet/zsasa/src/batch.zig`:

```zig
fn manifestJsonlOutputPath(allocator: Allocator, output_dir: []const u8, job_name: []const u8) ![]const u8 {
    const filename = try std.fmt.allocPrint(allocator, "{s}.jsonl", .{job_name});
    defer allocator.free(filename);
    return std.fs.path.join(allocator, &.{ output_dir, filename });
}

fn manifestPerFileOutputDir(allocator: Allocator, output_dir: []const u8, job_name: []const u8) ![]const u8 {
    return std.fs.path.join(allocator, &.{ output_dir, job_name });
}
```

- [ ] **Step 4: Add manifest execution functions**

Add these functions above `run` in `/Users/nagaet/zsasa/src/batch.zig`:

```zig
fn applyManifestGlobals(config: *BatchConfig, args: BatchArgs, globals: batch_manifest.Globals) void {
    if (!args.threads_explicit) if (globals.threads) |v| config.n_threads = v;
    if (!args.algorithm_explicit) if (globals.algorithm) |v| config.algorithm = parseAlgorithm(v);
    if (!args.n_points_explicit) if (globals.n_points) |v| config.n_points = v;
    if (!args.n_slices_explicit) if (globals.n_slices) |v| config.n_slices = v;
    if (!args.probe_radius_explicit) if (globals.probe_radius) |v| config.probe_radius = v;
    if (!args.precision_explicit) if (globals.precision) |v| config.precision = parsePrecision(v);
    if (!args.format_explicit) if (globals.format) |v| config.output_format = parseOutputFormat(v);
    if (!args.timing_explicit) if (globals.timing) |v| config.show_timing = v;
    if (!args.quiet_explicit) if (globals.quiet) |v| {
        config.quiet = v;
        config.show_progress = !v;
    };
    if (!args.classifier_explicit) if (globals.classifier) |v| config.classifier_type = parseClassifierType(v);
    if (!args.include_hydrogens_explicit) if (globals.include_hydrogens) |v| config.include_hydrogens = v;
    if (!args.include_hetatm_explicit) if (globals.include_hetatm) |v| config.include_hetatm = v;
    if (!args.use_bitmask_explicit) if (globals.use_bitmask) |v| config.use_bitmask = v;
    if (globals.auth_chain) |v| config.use_auth_chain = v;
}

fn applyCliOverrides(config: *BatchConfig, args: BatchArgs) void {
    if (args.threads_explicit) config.n_threads = args.n_threads;
    if (args.algorithm_explicit) config.algorithm = args.algorithm;
    if (args.n_points_explicit) config.n_points = args.n_points;
    if (args.n_slices_explicit) config.n_slices = args.n_slices;
    if (args.probe_radius_explicit) config.probe_radius = args.probe_radius;
    if (args.precision_explicit) config.precision = args.precision;
    if (args.format_explicit) config.output_format = args.output_format;
    if (args.timing_explicit) config.show_timing = args.show_timing;
    if (args.quiet_explicit) {
        config.quiet = args.quiet;
        config.show_progress = args.show_progress;
    }
    if (args.classifier_explicit) config.classifier_type = args.classifier_type;
    if (args.include_hydrogens_explicit) config.include_hydrogens = args.include_hydrogens;
    if (args.include_hetatm_explicit) config.include_hetatm = args.include_hetatm;
    if (args.use_bitmask_explicit) config.use_bitmask = args.use_bitmask;
}

fn runManifest(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (args.chain_filter != null) {
        std.debug.print("Error: --manifest cannot be combined with --chain; use [[jobs]].chains in the manifest\n", .{});
        return error.InvalidArgument;
    }

    const path = args.manifest_path.?;
    var manifest = batch_manifest.parseFile(allocator, io, path) catch |err| {
        std.debug.print("Error reading manifest '{s}': {s}\n", .{ path, @errorName(err) });
        return err;
    };
    defer manifest.deinit();

    const input_dir = args.input_path orelse manifest.globals.input_dir orelse {
        std.debug.print("Error: Missing input directory; set input_dir in manifest or pass <input_dir>\n", .{});
        return error.MissingArgument;
    };
    const output_dir = args.output_path orelse manifest.globals.output_dir;
    if (manifest.jobs.len > 1 and output_dir == null) {
        std.debug.print("Error: Manifest with multiple jobs requires output_dir or positional [output_dir]\n", .{});
        return error.MissingArgument;
    }

    var aggregate_success: usize = 0;
    var aggregate_failed: usize = 0;

    for (manifest.jobs) |job| {
        var config = BatchConfig{};
        applyManifestGlobals(&config, args, manifest.globals);
        applyCliOverrides(&config, args);
        config.include_hetatm = config.include_hetatm or (config.classifier_type == .ccd);
        config.store_atom_areas = (config.output_format == .jsonl);
        if (job.auth_chain) |v| config.use_auth_chain = v;
        config.chain_filter = job.chains;

        const job_output_dir = if (output_dir) |out| blk: {
            if (config.output_format == .jsonl) {
                break :blk out;
            }
            break :blk try manifestPerFileOutputDir(allocator, out, job.name);
        } else null;
        defer if (output_dir != null and config.output_format != .jsonl) if (job_output_dir) |p| allocator.free(p);

        const jsonl_output_path = if (output_dir) |out| blk: {
            if (config.output_format == .jsonl) break :blk try manifestJsonlOutputPath(allocator, out, job.name);
            break :blk null;
        } else null;
        defer if (jsonl_output_path) |p| allocator.free(p);

        if (!config.quiet) {
            std.debug.print("Manifest job: {s}\n", .{job.name});
        }

        var result = try runBatch(allocator, io, input_dir, job_output_dir, config, jsonl_output_path);
        defer result.deinit();
        aggregate_success += result.successful;
        aggregate_failed += result.failed;
        if (!config.quiet) result.printSummary(config.show_timing);
    }

    std.debug.print("Manifest complete: {d} successful, {d} failed\n", .{ aggregate_success, aggregate_failed });
}
```

If `runBatch` has a different signature in the current file, use the existing public entry point that accepts `(allocator, io, input_dir, output_dir, config, jsonl_output_path)`; do not duplicate sequential/parallel dispatch logic.

- [ ] **Step 5: Route `run` to manifest mode**

At the top of `pub fn run(allocator: Allocator, io: std.Io, args: BatchArgs) !void`, add:

```zig
    if (args.manifest_path != null) {
        return runManifest(allocator, io, args);
    }
```

Keep the rest of `run` as the non-manifest path.

- [ ] **Step 6: Run tests and verify manifest execution compiles**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task5-pass.log
```

Expected: PASS. The existing dispatcher is `runBatch(allocator, io, input_dir, output_dir, config, jsonl_output_path)`, so `runManifest` should call that function exactly once per manifest job with the resolved `BatchConfig`.

- [ ] **Step 7: Commit Task 5**

Run:

```bash
git add /Users/nagaet/zsasa/src/batch.zig
git commit -m "feat: run batch manifests"
```

---

### Task 6: Add CLI Help and Website Docs

**Files:**
- Modify: `/Users/nagaet/zsasa/src/batch.zig`
- Modify: `/Users/nagaet/zsasa/website/docs/cli/commands.md`

- [ ] **Step 1: Update batch help text**

In `printHelp` in `/Users/nagaet/zsasa/src/batch.zig`, add option lines:

```text
        \    --manifest=PATH     TOML manifest with one or more named batch jobs
        \    --chain=ID          Filter by chain ID for non-manifest batch (e.g. A or A,B)
        \    --auth-chain        Use auth_asym_id instead of label_asym_id for mmCIF chain matching
```

Add examples:

```text
        \    {s} batch --manifest bsa.toml
        \    {s} batch structures/ results/ --manifest bsa.toml
        \    {s} batch structures/ results_A/ --chain=A --format=jsonl
```

Update the format string argument list to include `program_name` for each new `{s}`.

- [ ] **Step 2: Update CLI docs with manifest example**

In `/Users/nagaet/zsasa/website/docs/cli/commands.md`, add this subsection under `batch`:

```markdown
#### Batch TOML manifest

Use `--manifest` to run several named batch jobs over the same input directory. CLI positional paths and explicit options override manifest values.

```bash
zsasa batch --manifest bsa.toml
zsasa batch structures/ results/ --manifest bsa.toml --threads=8
```

```toml
version = 1
input_dir = "structures"
output_dir = "results"
format = "jsonl"
use_bitmask = true
n_points = 128
classifier = "ccd"

[[jobs]]
name = "chain_A"
chains = ["A"]

[[jobs]]
name = "chain_B"
chains = ["B"]

[[jobs]]
name = "complex_AB"
chains = ["A", "B"]
```

For `format = "jsonl"`, each job writes one file such as `results/chain_A.jsonl`. For `json`, `compact`, and `csv`, each job writes a directory such as `results/chain_A/`.

Precedence is: built-in defaults < manifest globals < job settings < explicit CLI options. `--chain` is for non-manifest single-job batch mode; manifest jobs should use `[[jobs]].chains`.
```

- [ ] **Step 3: Run docs/build-neutral checks**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-task6-test.log
```

Expected: PASS for the full Zig test suite.

- [ ] **Step 4: Commit Task 6**

Run:

```bash
git add /Users/nagaet/zsasa/src/batch.zig /Users/nagaet/zsasa/website/docs/cli/commands.md
git commit -m "docs: document batch manifests"
```

---

### Task 7: End-to-End Manifest Smoke Test

**Files:**
- Create temporary files under `/tmp/zsasa-manifest-smoke/` only
- No source changes expected unless the smoke test finds a bug

- [ ] **Step 1: Create a small multi-chain PDB fixture**

Run:

```bash
rm -rf /tmp/zsasa-manifest-smoke
mkdir -p /tmp/zsasa-manifest-smoke/structures
cat > /tmp/zsasa-manifest-smoke/structures/two_chains.pdb <<'PDB'
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  N   ALA B   1       8.000   0.000   0.000  1.00 20.00           N
ATOM      5  CA  ALA B   1       9.500   0.000   0.000  1.00 20.00           C
ATOM      6  C   ALA B   1      10.000   1.400   0.000  1.00 20.00           C
TER
END
PDB
cat > /tmp/zsasa-manifest-smoke/bsa.toml <<'TOML'
version = 1
input_dir = "/tmp/zsasa-manifest-smoke/structures"
output_dir = "/tmp/zsasa-manifest-smoke/results"
format = "jsonl"
n_points = 32
classifier = "naccess"
quiet = true

[[jobs]]
name = "chain_A"
chains = ["A"]

[[jobs]]
name = "chain_B"
chains = ["B"]

[[jobs]]
name = "complex_AB"
chains = ["A", "B"]
TOML
```

- [ ] **Step 2: Run the manifest**

Run:

```bash
zig build run -- batch --manifest /tmp/zsasa-manifest-smoke/bsa.toml 2>&1 | tee /tmp/zsasa-manifest-smoke/run.log
```

Expected: exit code 0 and output includes `Manifest complete`.

- [ ] **Step 3: Verify three JSONL outputs exist**

Run:

```bash
ls -1 /tmp/zsasa-manifest-smoke/results
wc -l /tmp/zsasa-manifest-smoke/results/*.jsonl
```

Expected:

```text
chain_A.jsonl
chain_B.jsonl
complex_AB.jsonl
```

Each JSONL file should have one line.

- [ ] **Step 4: Verify atom counts reflect chain selection**

Run:

```bash
python3 - <<'PY'
import json
from pathlib import Path
base = Path('/tmp/zsasa-manifest-smoke/results')
counts = {}
for name in ['chain_A', 'chain_B', 'complex_AB']:
    row = json.loads((base / f'{name}.jsonl').read_text().strip())
    counts[name] = len(row['atom_areas'])
print(counts)
assert counts['chain_A'] == 3, counts
assert counts['chain_B'] == 3, counts
assert counts['complex_AB'] == 6, counts
PY
```

Expected: prints `{'chain_A': 3, 'chain_B': 3, 'complex_AB': 6}` and exits 0.

- [ ] **Step 5: Verify CLI positional path override**

Run:

```bash
cat > /tmp/zsasa-manifest-smoke/no_paths.toml <<'TOML'
version = 1
format = "jsonl"
n_points = 32
classifier = "naccess"
quiet = true

[[jobs]]
name = "chain_A"
chains = ["A"]
TOML
zig build run -- batch /tmp/zsasa-manifest-smoke/structures /tmp/zsasa-manifest-smoke/override-results --manifest /tmp/zsasa-manifest-smoke/no_paths.toml
ls -1 /tmp/zsasa-manifest-smoke/override-results
```

Expected: `chain_A.jsonl` exists under `/tmp/zsasa-manifest-smoke/override-results`.

- [ ] **Step 6: Commit source changes made during the smoke test**

When the smoke test leads to a source change, run:

```bash
git add /Users/nagaet/zsasa/src /Users/nagaet/zsasa/website/docs/cli/commands.md
git commit -m "fix: stabilize batch manifest execution"
```

When no code changes were needed, do not create an empty commit.

---

### Task 8: Final Verification and Review Prep

**Files:**
- No source changes expected unless verification finds a bug

- [ ] **Step 1: Run full Zig tests**

Run:

```bash
zig build test 2>&1 | tee /tmp/zsasa-final-zig-test.log
```

Expected: PASS.

- [ ] **Step 2: Run manifest smoke test again from Task 7**

Run:

```bash
zig build run -- batch --manifest /tmp/zsasa-manifest-smoke/bsa.toml 2>&1 | tee /tmp/zsasa-final-manifest-smoke.log
python3 - <<'PY'
import json
from pathlib import Path
base = Path('/tmp/zsasa-manifest-smoke/results')
counts = {}
for name in ['chain_A', 'chain_B', 'complex_AB']:
    row = json.loads((base / f'{name}.jsonl').read_text().strip())
    counts[name] = len(row['atom_areas'])
print(counts)
assert counts == {'chain_A': 3, 'chain_B': 3, 'complex_AB': 6}, counts
PY
```

Expected: both commands exit 0.

- [ ] **Step 3: Inspect git history and status**

Run:

```bash
git log --oneline --decorate -8
git status --short
```

Expected: recent commits include the design commit and task commits. `git status --short` should show no modified tracked files. It may still show the pre-existing untracked `/Users/nagaet/zsasa/zig-pkg/` directory; do not add it unless the user explicitly asks.

- [ ] **Step 4: Summarize implementation for review**

Prepare a concise review summary containing:

```markdown
## Summary
- Added TOML batch manifests with named jobs and A/B/AB chain selections.
- Added batch chain filtering for non-manifest mode.
- Documented manifest precedence and output layout.

## Verification
- `zig build test`
- `zig build run -- batch --manifest /tmp/zsasa-manifest-smoke/bsa.toml`
- JSONL atom count smoke: chain_A=3, chain_B=3, complex_AB=6
```

Do not claim verification passed unless the exact commands above passed in the current session.

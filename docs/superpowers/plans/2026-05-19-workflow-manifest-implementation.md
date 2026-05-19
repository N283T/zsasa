# Workflow Manifest Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Promote batch manifests into reusable workflow files for `calc` and `batch`, and make custom classifier configuration TOML-only.

**Architecture:** Add a parser-only `workflow_manifest` module that normalizes sectioned workflow files and legacy flat batch manifests into one data model. `calc.zig` and `batch.zig` resolve that parsed workflow into their existing argument/config structures, preserving current direct CLI behavior and explicit CLI override precedence. `classifier_parser.zig` becomes a TOML-only compatibility facade around `toml_classifier_parser.zig`.

**Tech Stack:** Zig 0.15-style std APIs, existing `toml_parser.zig`, existing `toml_classifier_parser.zig`, existing CLI modules `calc.zig` and `batch.zig`, Docusaurus docs under `website/docs`.

---

## Scope Check

The approved spec covers one coherent subsystem: workflow-file parsing plus command consumers and classifier format cleanup. The work is large enough to split by files but each task leaves the repository buildable and testable.

## File Structure

- Create `src/workflow_manifest.zig`: TOML workflow parser, normalized data structs, validation, legacy flat batch manifest mapping, parser tests.
- Modify `src/root.zig`: export `workflow_manifest` and remove or stop exporting `batch_manifest` after consumers move.
- Modify `src/classifier_parser.zig`: reduce to TOML-only file loading facade; remove FreeSASA-format parser implementation and tests.
- Modify `src/calc.zig`: add `--workflow`, explicit-option flags, workflow resolution, TOML-only config wording.
- Modify `src/batch.zig`: import `workflow_manifest`, add `--workflow` alias, resolve new workflow sections, support referenced custom classifier for batch, preserve `--manifest` alias.
- Modify docs: `README.md`, `CHANGELOG.md`, `website/docs/cli/commands.md`, `website/docs/cli/input.md`, `website/docs/guide/classifiers.mdx`.

---

### Task 1: Add normalized workflow manifest parser

**Files:**
- Create: `src/workflow_manifest.zig`
- Modify: `src/root.zig`
- Reference: `src/batch_manifest.zig`
- Test: Zig unit tests in `src/workflow_manifest.zig`

- [ ] **Step 1: Write failing parser tests**

Create `src/workflow_manifest.zig` with only imports, empty public structs if needed for compilation, and these tests at the bottom. The tests should fail because `parse()` and fields are not implemented yet.

```zig
const std = @import("std");
const Allocator = std.mem.Allocator;
const toml_parser = @import("toml_parser.zig");

pub const WorkflowError = error{
    UnsupportedVersion,
    InvalidKind,
    MissingJobName,
    DuplicateJobName,
    UnsafeJobName,
    InvalidFieldType,
    InvalidClassifierConfig,
    UnknownField,
    NoJobs,
};

pub const Error = WorkflowError || Allocator.Error || toml_parser.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError;

pub const Input = struct { path: ?[]const u8 = null, dir: ?[]const u8 = null, chain: ?[]const u8 = null, model: ?u32 = null, mol: ?[]const u8 = null };
pub const Output = struct { path: ?[]const u8 = null, dir: ?[]const u8 = null, format: ?[]const u8 = null };
pub const Calculation = struct { algorithm: ?[]const u8 = null, threads: ?usize = null, probe_radius: ?f64 = null, n_points: ?u32 = null, n_slices: ?u32 = null, precision: ?[]const u8 = null, include_hydrogens: ?bool = null, include_hetatm: ?bool = null, use_bitmask: ?bool = null, timing: ?bool = null, quiet: ?bool = null, auth_chain: ?bool = null, residue_map: ?bool = null, per_residue: ?bool = null, rsa: ?bool = null, polar: ?bool = null, validate_only: ?bool = null };
pub const ClassifierConfig = struct { type: ?[]const u8 = null, config: ?[]const u8 = null, ccd: ?[]const u8 = null, sdf: ?[]const []const u8 = null };
pub const Job = struct { name: []const u8, chains: ?[]const []const u8 = null, auth_chain: ?bool = null };
pub const Workflow = struct {
    allocator: Allocator,
    content: []const u8,
    input: Input = .{},
    output: Output = .{},
    calculation: Calculation = .{},
    classifier: ClassifierConfig = .{},
    jobs: []Job = &.{},
    is_legacy_batch_manifest: bool = false,
    pub fn deinit(self: *Workflow) void { _ = self; }
};

pub fn parse(allocator: Allocator, content: []const u8) Error!Workflow { _ = allocator; _ = content; return error.UnsupportedVersion; }

test "parse sectioned calc workflow" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\path = "structure.cif"
        \\chain = "A,B"
        \\model = 1
        \\mol = "LIG"
        \\
        \\[output]
        \\path = "result.json"
        \\format = "json"
        \\
        \\[calculation]
        \\algorithm = "sr"
        \\probe_radius = 1.4
        \\n_points = 128
        \\n_slices = 20
        \\precision = "f64"
        \\use_bitmask = true
        \\include_hydrogens = false
        \\include_hetatm = true
        \\auth_chain = true
        \\per_residue = true
        \\rsa = true
        \\polar = true
        \\validate_only = false
        \\
        \\[classifier]
        \\type = "custom"
        \\config = "my_radii.toml"
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqual(false, workflow.is_legacy_batch_manifest);
    try std.testing.expectEqualStrings("structure.cif", workflow.input.path.?);
    try std.testing.expectEqualStrings("A,B", workflow.input.chain.?);
    try std.testing.expectEqual(@as(u32, 1), workflow.input.model.?);
    try std.testing.expectEqualStrings("LIG", workflow.input.mol.?);
    try std.testing.expectEqualStrings("result.json", workflow.output.path.?);
    try std.testing.expectEqualStrings("json", workflow.output.format.?);
    try std.testing.expectEqualStrings("sr", workflow.calculation.algorithm.?);
    try std.testing.expectEqual(@as(f64, 1.4), workflow.calculation.probe_radius.?);
    try std.testing.expectEqual(@as(u32, 128), workflow.calculation.n_points.?);
    try std.testing.expectEqual(@as(u32, 20), workflow.calculation.n_slices.?);
    try std.testing.expectEqualStrings("f64", workflow.calculation.precision.?);
    try std.testing.expectEqual(true, workflow.calculation.use_bitmask.?);
    try std.testing.expectEqual(false, workflow.calculation.include_hydrogens.?);
    try std.testing.expectEqual(true, workflow.calculation.include_hetatm.?);
    try std.testing.expectEqual(true, workflow.calculation.auth_chain.?);
    try std.testing.expectEqual(true, workflow.calculation.per_residue.?);
    try std.testing.expectEqual(true, workflow.calculation.rsa.?);
    try std.testing.expectEqual(true, workflow.calculation.polar.?);
    try std.testing.expectEqual(false, workflow.calculation.validate_only.?);
    try std.testing.expectEqualStrings("custom", workflow.classifier.type.?);
    try std.testing.expectEqualStrings("my_radii.toml", workflow.classifier.config.?);
    try std.testing.expectEqual(@as(usize, 0), workflow.jobs.len);
}

test "parse sectioned batch workflow with jobs" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "structures"
        \\
        \\[output]
        \\dir = "results"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\threads = 8
        \\n_points = 128
        \\residue_map = true
        \\
        \\[classifier]
        \\type = "ccd"
        \\ccd = "components.zsdc"
        \\sdf = ["ligand.sdf", "cofactor.sdf"]
        \\
        \\[[jobs]]
        \\name = "chain_A"
        \\chains = ["A"]
        \\
        \\[[jobs]]
        \\name = "complex_AB"
        \\chains = ["A", "B"]
        \\auth_chain = true
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqualStrings("structures", workflow.input.dir.?);
    try std.testing.expectEqualStrings("results", workflow.output.dir.?);
    try std.testing.expectEqualStrings("jsonl", workflow.output.format.?);
    try std.testing.expectEqual(@as(usize, 8), workflow.calculation.threads.?);
    try std.testing.expectEqual(true, workflow.calculation.residue_map.?);
    try std.testing.expectEqualStrings("ccd", workflow.classifier.type.?);
    try std.testing.expectEqualStrings("components.zsdc", workflow.classifier.ccd.?);
    try std.testing.expectEqual(@as(usize, 2), workflow.classifier.sdf.?.len);
    try std.testing.expectEqualStrings("ligand.sdf", workflow.classifier.sdf.?[0]);
    try std.testing.expectEqual(@as(usize, 2), workflow.jobs.len);
    try std.testing.expectEqualStrings("chain_A", workflow.jobs[0].name);
    try std.testing.expectEqualStrings("A", workflow.jobs[0].chains.?[0]);
    try std.testing.expectEqualStrings("complex_AB", workflow.jobs[1].name);
    try std.testing.expectEqualStrings("B", workflow.jobs[1].chains.?[1]);
    try std.testing.expectEqual(true, workflow.jobs[1].auth_chain.?);
}

test "parse legacy flat batch manifest" {
    const input =
        \\version = 1
        \\input_dir = "structures"
        \\output_dir = "results"
        \\format = "jsonl"
        \\classifier = "ccd"
        \\ccd = "components.zsdc"
        \\sdf = "ligand.sdf"
        \\threads = 4
        \\n_points = 128
        \\use_bitmask = true
        \\
        \\[[jobs]]
        \\name = "chain_A"
        \\chains = ["A"]
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqual(true, workflow.is_legacy_batch_manifest);
    try std.testing.expectEqualStrings("structures", workflow.input.dir.?);
    try std.testing.expectEqualStrings("results", workflow.output.dir.?);
    try std.testing.expectEqualStrings("jsonl", workflow.output.format.?);
    try std.testing.expectEqualStrings("ccd", workflow.classifier.type.?);
    try std.testing.expectEqualStrings("components.zsdc", workflow.classifier.ccd.?);
    try std.testing.expectEqual(@as(usize, 1), workflow.classifier.sdf.?.len);
    try std.testing.expectEqualStrings("ligand.sdf", workflow.classifier.sdf.?[0]);
    try std.testing.expectEqual(@as(usize, 4), workflow.calculation.threads.?);
    try std.testing.expectEqual(@as(u32, 128), workflow.calculation.n_points.?);
    try std.testing.expectEqual(true, workflow.calculation.use_bitmask.?);
    try std.testing.expectEqualStrings("chain_A", workflow.jobs[0].name);
}

test "reject invalid classifier combinations" {
    const custom_with_ccd =
        \\version = 1
        \\[classifier]
        \\type = "custom"
        \\config = "my_radii.toml"
        \\ccd = "components.zsdc"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, custom_with_ccd));

    const custom_without_config =
        \\version = 1
        \\[classifier]
        \\type = "custom"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, custom_without_config));

    const builtin_with_config =
        \\version = 1
        \\[classifier]
        \\type = "ccd"
        \\config = "my_radii.toml"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, builtin_with_config));
}
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
zig build test
```

Expected: FAIL with errors from `src/workflow_manifest.zig` because `parse()` returns `UnsupportedVersion` and `Workflow.deinit()` does not free owned memory.

- [ ] **Step 3: Implement parser ownership and section parsing**

Replace the placeholder `Workflow.deinit()` and `parse()` with a real implementation. Port the helper style from `src/batch_manifest.zig`, but normalize into the new structs. The important implementation shape is:

```zig
pub const Workflow = struct {
    allocator: Allocator,
    content: []const u8,
    input: Input = .{},
    output: Output = .{},
    calculation: Calculation = .{},
    classifier: ClassifierConfig = .{},
    jobs: []Job = &.{},
    is_legacy_batch_manifest: bool = false,

    pub fn deinit(self: *Workflow) void {
        if (self.classifier.sdf) |items| self.allocator.free(items);
        for (self.jobs) |job| {
            if (job.chains) |chains| self.allocator.free(chains);
        }
        self.allocator.free(self.jobs);
        self.allocator.free(self.content);
        self.* = undefined;
    }
};

pub fn parse(allocator: Allocator, content: []const u8) Error!Workflow {
    const owned_content = try allocator.dupe(u8, content);
    errdefer allocator.free(owned_content);
    return parseOwned(allocator, owned_content);
}

fn parseOwned(allocator: Allocator, owned_content: []const u8) Error!Workflow {
    var doc = try toml_parser.parse(allocator, owned_content);
    defer doc.deinit();

    const root = doc.getTable("") orelse return error.UnsupportedVersion;
    try validateVersion(root);
    try validateKind(root);

    const is_legacy = isLegacyBatchManifest(root);
    var workflow = Workflow{
        .allocator = allocator,
        .content = owned_content,
        .is_legacy_batch_manifest = is_legacy,
    };
    errdefer workflow.deinit();

    if (is_legacy) {
        try parseLegacyRootIntoWorkflow(allocator, root, &workflow);
    } else {
        if (doc.getTable("input")) |table| workflow.input = try parseInput(table);
        if (doc.getTable("output")) |table| workflow.output = try parseOutput(table);
        if (doc.getTable("calculation")) |table| workflow.calculation = try parseCalculation(table);
        if (doc.getTable("classifier")) |table| workflow.classifier = try parseClassifier(allocator, table);
    }

    workflow.jobs = try parseJobs(allocator, doc.array_tables);
    try validateClassifier(workflow.classifier);
    return workflow;
}
```

Add the helper functions with these exact responsibilities:

```zig
fn validateVersion(root: toml_parser.Table) WorkflowError!void;
fn validateKind(root: toml_parser.Table) WorkflowError!void;
fn isLegacyBatchManifest(root: toml_parser.Table) bool;
fn parseLegacyRootIntoWorkflow(allocator: Allocator, root: toml_parser.Table, workflow: *Workflow) Error!void;
fn parseInput(table: toml_parser.Table) WorkflowError!Input;
fn parseOutput(table: toml_parser.Table) WorkflowError!Output;
fn parseCalculation(table: toml_parser.Table) WorkflowError!Calculation;
fn parseClassifier(allocator: Allocator, table: toml_parser.Table) Error!ClassifierConfig;
fn validateClassifier(config: ClassifierConfig) WorkflowError!void;
fn parseJobs(allocator: Allocator, array_tables: []const toml_parser.Document.ArrayTable) Error![]Job;
```

Keep the value helpers from `src/batch_manifest.zig` and expand them only as needed:

```zig
fn findValue(entries: []const toml_parser.Value.Entry, key: []const u8) ?toml_parser.Value;
fn optionalString(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?[]const u8;
fn optionalBool(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?bool;
fn optionalFloat(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?f64;
fn optionalUsize(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?usize;
fn optionalU32(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?u32;
fn optionalStringArray(allocator: Allocator, entries: []const toml_parser.Value.Entry, key: []const u8) Error!?[]const []const u8;
fn optionalStringOrStringArray(allocator: Allocator, entries: []const toml_parser.Value.Entry, key: []const u8) Error!?[]const []const u8;
fn isSafeJobName(name: []const u8) bool;
```

- [ ] **Step 4: Export the new module**

Modify `src/root.zig`:

```zig
pub const workflow_manifest = @import("workflow_manifest.zig");
```

In the root test block, add:

```zig
_ = workflow_manifest;
```

Keep `batch_manifest` exported until `src/batch.zig` is migrated in Task 4; removing it now would break the build.

- [ ] **Step 5: Run tests to verify pass**

Run:

```bash
zig build test
```

Expected: PASS for the new `workflow_manifest` tests, with no regressions in existing tests.

- [ ] **Step 6: Commit**

```bash
git add src/workflow_manifest.zig src/root.zig
git commit -m "feat: add workflow manifest parser"
```

---

### Task 2: Make custom classifier config TOML-only

**Files:**
- Modify: `src/classifier_parser.zig`
- Reference: `src/toml_classifier_parser.zig`
- Test: Zig unit tests in `src/classifier_parser.zig`

- [ ] **Step 1: Replace FreeSASA parser tests with TOML-only tests**

In `src/classifier_parser.zig`, remove tests named `parseConfig minimal`, `parseConfig with comments`, `parseConfig ANY fallback`, and other FreeSASA-format parser tests. Add these TOML-only facade tests:

```zig
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

    var result = try parseConfig(std.testing.allocator, config);
    defer result.deinit();

    try std.testing.expectEqualStrings("test", result.name);
    try std.testing.expectEqual(@as(?f64, 1.87), result.getRadius("ALA", "CA"));
    try std.testing.expectEqual(AtomClass.apolar, result.getClass("ALA", "CA"));
}

test "custom classifier parser rejects legacy FreeSASA content" {
    const config =
        \\name: legacy
        \\types:
        \\C_ALI 1.87 apolar
        \\atoms:
        \\ALA CA C_ALI
    ;

    try std.testing.expectError(error.UnsupportedLegacyFormat, parseConfig(std.testing.allocator, config));
}

test "custom classifier file extension must be TOML" {
    try std.testing.expectEqual(true, isTomlPath("my_radii.toml"));
    try std.testing.expectEqual(false, isTomlPath("my_radii.config"));
    try std.testing.expectEqual(false, isTomlPath("my_radii.conf"));
}
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
zig build test
```

Expected: FAIL because `UnsupportedLegacyFormat` and `isTomlPath()` do not exist yet, and `parseConfig()` still parses FreeSASA content.

- [ ] **Step 3: Reduce `classifier_parser.zig` to a TOML-only facade**

Replace the FreeSASA parser implementation with a small wrapper around `toml_classifier_parser.parseConfig`. Keep the public module name to minimize changes in `calc.zig` and future call sites.

```zig
//! TOML-only custom classifier configuration facade.
//!
//! Custom classifier configs use the TOML format parsed by
//! `toml_classifier_parser.zig`. The legacy FreeSASA-style text format was
//! removed; use TOML `[types]` and `[[atoms]]` entries instead.

const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const Classifier = classifier.Classifier;
const AtomClass = classifier.AtomClass;
const toml_classifier_parser = @import("toml_classifier_parser.zig");
const toml_parser = @import("toml_parser.zig");

pub const ParseError = error{
    UnsupportedLegacyFormat,
    UnsupportedConfigExtension,
    FileReadError,
};

pub const Error = ParseError || Allocator.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError || toml_parser.TomlError || @import("classifier_parser.zig").ParseError;
```

Do not use the self-import line shown above literally. Instead, use the parse errors already exported by `toml_classifier_parser`:

```zig
pub const Error = ParseError || toml_classifier_parser.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError;
```

Implement:

```zig
pub fn parseConfig(allocator: Allocator, content: []const u8) Error!Classifier {
    if (looksLikeLegacyFreeSasa(content)) return error.UnsupportedLegacyFormat;
    return toml_classifier_parser.parseConfig(allocator, content);
}

pub fn parseConfigFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Classifier {
    if (!isTomlPath(path)) return error.UnsupportedConfigExtension;
    const content = try readFileContent(allocator, io, path);
    defer allocator.free(content);
    return parseConfig(allocator, content);
}

fn isTomlPath(path: []const u8) bool {
    return std.mem.endsWith(u8, path, ".toml");
}

fn looksLikeLegacyFreeSasa(content: []const u8) bool {
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |raw_line| {
        const line = std.mem.trim(u8, raw_line, " \t\r");
        if (line.len == 0 or std.mem.startsWith(u8, line, "#")) continue;
        if (std.mem.eql(u8, line, "types:") or std.mem.eql(u8, line, "atoms:")) return true;
        if (std.mem.startsWith(u8, line, "name:")) return true;
    }
    return false;
}

fn readFileContent(allocator: Allocator, io: std.Io, path: []const u8) Error![]u8 {
    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch return error.FileReadError;
    defer file.close(io);
    var read_buf: [65536]u8 = undefined;
    var r = file.reader(io, &read_buf);
    return r.interface.allocRemaining(allocator, .unlimited) catch return error.FileReadError;
}
```

Ensure `toml_classifier_parser.zig` no longer imports `classifier_parser.zig` for `ParseError`; otherwise this facade creates an import cycle. Move the shared parse error enum into `toml_classifier_parser.zig`:

```zig
pub const ParseError = error{
    UnknownSection,
    InvalidTypeDefinition,
    InvalidAtomDefinition,
    InvalidRadius,
    InvalidClass,
    UndefinedType,
    DuplicateType,
    AtomsBeforeTypes,
    FileReadError,
    OutOfMemory,
};
```

Then remove this line from `toml_classifier_parser.zig`:

```zig
const classifier_parser = @import("classifier_parser.zig");
```

- [ ] **Step 4: Update calc error messages for TOML-only config**

In `src/calc.zig`, change the catch block for custom config loading to show migration guidance for the two new errors:

```zig
var custom_classifier = classifier_parser.parseConfigFile(allocator, io, config_path) catch |err| {
    switch (err) {
        error.UnsupportedConfigExtension => std.debug.print(
            "Error loading config file '{s}': custom classifier configs are TOML-only; rename or convert the file to .toml\n",
            .{config_path},
        ),
        error.UnsupportedLegacyFormat => std.debug.print(
            "Error loading config file '{s}': FreeSASA-style custom classifier configs are no longer supported; convert to TOML [types] and [[atoms]]\n",
            .{config_path},
        ),
        else => std.debug.print("Error loading config file '{s}': {s}\n", .{ config_path, @errorName(err) }),
    }
    std.process.exit(1);
};
```

- [ ] **Step 5: Run tests to verify pass**

Run:

```bash
zig build test
```

Expected: PASS. Existing `toml_classifier_parser` tests still pass. FreeSASA parser tests no longer exist.

- [ ] **Step 6: Commit**

```bash
git add src/classifier_parser.zig src/toml_classifier_parser.zig src/calc.zig
git commit -m "feat: require TOML custom classifiers"
```

---

### Task 3: Add calc workflow parsing and resolution

**Files:**
- Modify: `src/calc.zig`
- Reference: `src/workflow_manifest.zig`
- Test: Zig unit tests in `src/calc.zig`

- [ ] **Step 1: Add failing CLI parse and resolver tests**

At the bottom of `src/calc.zig`, add:

```zig
test "CalcArgs --workflow=FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "--workflow=sasa.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("sasa.toml", parsed.workflow_path.?);
}

test "CalcArgs --workflow FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "--workflow", "sasa.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("sasa.toml", parsed.workflow_path.?);
}

test "calc workflow applies fields when CLI did not override" {
    var args = CalcArgs{ .workflow_path = "sasa.toml" };
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .input = .{ .path = "structure.cif", .chain = "A", .model = 1, .mol = "LIG" },
        .output = .{ .path = "result.json", .format = "compact" },
        .calculation = .{ .algorithm = "lr", .n_slices = 44, .probe_radius = 1.5, .precision = "f32", .per_residue = true, .rsa = true, .polar = true, .auth_chain = true, .include_hetatm = true, .quiet = true },
        .classifier = .{ .type = "custom", .config = "my_radii.toml" },
        .jobs = &.{},
    };

    try applyWorkflowToCalcArgs(&args, workflow);

    try std.testing.expectEqualStrings("structure.cif", args.input_path.?);
    try std.testing.expectEqualStrings("result.json", args.output_path);
    try std.testing.expectEqual(OutputFormat.compact, args.output_format);
    try std.testing.expectEqual(Algorithm.lr, args.algorithm);
    try std.testing.expectEqual(@as(u32, 44), args.n_slices);
    try std.testing.expectEqual(@as(f64, 1.5), args.probe_radius);
    try std.testing.expectEqual(Precision.f32, args.precision);
    try std.testing.expectEqualStrings("A", args.chain_filter.?);
    try std.testing.expectEqual(@as(u32, 1), args.model_num.?);
    try std.testing.expectEqualStrings("LIG", args.mol_selector.?);
    try std.testing.expectEqual(true, args.per_residue);
    try std.testing.expectEqual(true, args.rsa);
    try std.testing.expectEqual(true, args.polar);
    try std.testing.expectEqual(true, args.use_auth_chain);
    try std.testing.expectEqual(true, args.include_hetatm);
    try std.testing.expectEqual(true, args.quiet);
    try std.testing.expectEqualStrings("my_radii.toml", args.config_path.?);
}

test "calc CLI explicit classifier overrides workflow classifier" {
    var args = CalcArgs{ .workflow_path = "sasa.toml", .classifier_type = .naccess, .classifier_explicit = true };
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{ .type = "custom", .config = "my_radii.toml" },
        .jobs = &.{},
    };

    try applyWorkflowToCalcArgs(&args, workflow);
    try std.testing.expectEqual(ClassifierType.naccess, args.classifier_type.?);
    try std.testing.expect(args.config_path == null);
}
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
zig build test
```

Expected: FAIL because `workflow_path`, explicit flags, `workflow_manifest` import, and `applyWorkflowToCalcArgs()` do not exist.

- [ ] **Step 3: Import workflow module and add explicit flags**

At the top of `src/calc.zig`, add:

```zig
const workflow_manifest = @import("workflow_manifest.zig");
```

Extend `CalcArgs` with:

```zig
workflow_path: ?[]const u8 = null,
threads_explicit: bool = false,
probe_radius_explicit: bool = false,
n_points_explicit: bool = false,
n_slices_explicit: bool = false,
algorithm_explicit: bool = false,
precision_explicit: bool = false,
format_explicit: bool = false,
classifier_explicit: bool = false,
config_explicit: bool = false,
chain_explicit: bool = false,
model_explicit: bool = false,
auth_chain_explicit: bool = false,
include_hydrogens_explicit: bool = false,
include_hetatm_explicit: bool = false,
per_residue_explicit: bool = false,
rsa_explicit: bool = false,
polar_explicit: bool = false,
use_bitmask_explicit: bool = false,
ccd_explicit: bool = false,
sdf_explicit: bool = false,
mol_explicit: bool = false,
quiet_explicit: bool = false,
validate_explicit: bool = false,
timing_explicit: bool = false,
```

When parsing each CLI option, set the matching explicit flag before assigning the value. For example:

```zig
else if (std.mem.startsWith(u8, arg, "--algorithm=")) {
    result.algorithm_explicit = true;
    const value = arg["--algorithm=".len..];
    result.algorithm = parseAlgorithm(value);
}
```

Apply this pattern to every workflow-overridable CLI option listed in the new fields.

- [ ] **Step 4: Parse `--workflow`**

Add to `parseArgs()` before unknown-option handling:

```zig
else if (std.mem.startsWith(u8, arg, "--workflow=")) {
    result.workflow_path = arg["--workflow=".len..];
} else if (std.mem.eql(u8, arg, "--workflow")) {
    i += 1;
    if (i >= args.len) {
        std.debug.print("Error: Missing value for --workflow\n", .{});
        std.process.exit(1);
    }
    result.workflow_path = args[i];
}
```

- [ ] **Step 5: Implement calc workflow resolver**

Add helper functions near the parse helpers:

```zig
fn applyWorkflowToCalcArgs(args: *CalcArgs, workflow: workflow_manifest.Workflow) !void {
    if (!args.*.chain_explicit) if (workflow.input.chain) |v| args.*.chain_filter = v;
    if (!args.*.model_explicit) if (workflow.input.model) |v| args.*.model_num = v;
    if (!args.*.mol_explicit) if (workflow.input.mol) |v| args.*.mol_selector = v;
    if (args.*.input_path == null) args.*.input_path = workflow.input.path;
    if (!args.*.output_path_explicit) if (workflow.output.path) |v| args.*.output_path = v;
    if (!args.*.format_explicit) if (workflow.output.format) |v| args.*.output_format = parseOutputFormat(v);

    try applyWorkflowCalculationToCalcArgs(args, workflow.calculation);
    try applyWorkflowClassifierToCalcArgs(args, workflow.classifier);
}

fn applyWorkflowCalculationToCalcArgs(args: *CalcArgs, c: workflow_manifest.Calculation) !void {
    if (!args.*.threads_explicit) if (c.threads) |v| args.*.n_threads = v;
    if (!args.*.algorithm_explicit) if (c.algorithm) |v| args.*.algorithm = parseAlgorithm(v);
    if (!args.*.probe_radius_explicit) if (c.probe_radius) |v| args.*.probe_radius = validateWorkflowProbeRadius(v) catch |err| {
        std.debug.print("Error: workflow probe_radius must be finite and between 0 and 10 Angstroms: {d}\n", .{v});
        return err;
    };
    if (!args.*.n_points_explicit) if (c.n_points) |v| args.*.n_points = validateWorkflowNPoints(v) catch |err| {
        std.debug.print("Error: workflow n_points must be between 1 and 10000: {d}\n", .{v});
        return err;
    };
    if (!args.*.n_slices_explicit) if (c.n_slices) |v| args.*.n_slices = validateWorkflowNSlices(v) catch |err| {
        std.debug.print("Error: workflow n_slices must be between 1 and 1000: {d}\n", .{v});
        return err;
    };
    if (!args.*.precision_explicit) if (c.precision) |v| args.*.precision = parsePrecision(v);
    if (!args.*.use_bitmask_explicit) if (c.use_bitmask) |v| args.*.use_bitmask = v;
    if (!args.*.include_hydrogens_explicit) if (c.include_hydrogens) |v| args.*.include_hydrogens = v;
    if (!args.*.include_hetatm_explicit) if (c.include_hetatm) |v| args.*.include_hetatm = v;
    if (!args.*.auth_chain_explicit) if (c.auth_chain) |v| args.*.use_auth_chain = v;
    if (!args.*.per_residue_explicit) if (c.per_residue) |v| args.*.per_residue = v;
    if (!args.*.rsa_explicit) if (c.rsa) |v| args.*.rsa = v;
    if (!args.*.polar_explicit) if (c.polar) |v| args.*.polar = v;
    if (!args.*.validate_explicit) if (c.validate_only) |v| args.*.validate_only = v;
    if (!args.*.quiet_explicit) if (c.quiet) |v| args.*.quiet = v;
    if (!args.*.timing_explicit) if (c.timing) |v| args.*.show_timing = v;
}

fn applyWorkflowClassifierToCalcArgs(args: *CalcArgs, c: workflow_manifest.ClassifierConfig) !void {
    const classifier_type = c.type orelse return;
    if (std.mem.eql(u8, classifier_type, "custom")) {
        if (!args.*.classifier_explicit and !args.*.config_explicit) args.*.config_path = c.config;
        return;
    }
    if (!args.*.classifier_explicit and !args.*.config_explicit) args.*.classifier_type = parseClassifierType(classifier_type);
    if (!args.*.ccd_explicit) args.*.ccd_path = c.ccd;
    if (!args.*.sdf_explicit) if (c.sdf) |items| {
        args.*.sdf_paths = .{};
        for (items) |item| args.*.sdf_paths.append(item) catch return error.InvalidArgument;
    };
}
```

Add reusable validators mirroring the current parse helpers:

```zig
fn validateWorkflowProbeRadius(radius: f64) !f64;
fn validateWorkflowNPoints(n: u32) !u32;
fn validateWorkflowNSlices(n: u32) !u32;
```

- [ ] **Step 6: Load workflow before required argument validation**

At the start of `run()` before checking `input_path`, convert `args` to mutable and apply the workflow:

```zig
pub fn run(allocator: std.mem.Allocator, io: std.Io, args: CalcArgs) !void {
    var effective_args = args;
    if (effective_args.workflow_path) |workflow_path| {
        var workflow = workflow_manifest.parseFile(allocator, io, workflow_path) catch |err| {
            std.debug.print("Error reading workflow '{s}': {s}\n", .{ workflow_path, @errorName(err) });
            return err;
        };
        defer workflow.deinit();
        try applyWorkflowToCalcArgs(&effective_args, workflow);
    }

    const input_path = effective_args.input_path orelse {
        std.debug.print("Error: Missing input file\n", .{});
        std.debug.print("Usage: zsasa calc [OPTIONS] <input> [output.json]\n", .{});
        return error.MissingArgument;
    };
```

Then use `effective_args` throughout `run()` instead of `args`. Avoid creating a second `var effective_args = args;` later in the function; merge the existing CCD/HETATM adjustment into the same mutable variable.

- [ ] **Step 7: Run tests to verify pass**

Run:

```bash
zig build test
```

Expected: PASS, including new calc workflow tests.

- [ ] **Step 8: Commit**

```bash
git add src/calc.zig
git commit -m "feat: support calc workflow files"
```

---

### Task 4: Migrate batch to workflow files and add custom classifier support

**Files:**
- Modify: `src/batch.zig`
- Modify: `src/root.zig`
- Delete: `src/batch_manifest.zig` only after imports no longer need it
- Test: Zig unit tests in `src/batch.zig` and `src/workflow_manifest.zig`

- [ ] **Step 1: Add failing batch CLI and resolver tests**

In `src/batch.zig`, add or update tests:

```zig
test "BatchArgs --workflow" {
    const args = [_][]const u8{ "zsasa", "batch", "--workflow", "bsa.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("bsa.toml", parsed.workflow_path.?);
}

test "BatchArgs --manifest remains alias" {
    const args = [_][]const u8{ "zsasa", "batch", "--manifest", "bsa.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("bsa.toml", parsed.workflow_path.?);
}

test "workflow residue_map applies when CLI does not override format" {
    var config = BatchConfig{};
    const args = BatchArgs{};
    const calculation = workflow_manifest.Calculation{ .residue_map = true };
    const output = workflow_manifest.Output{ .format = "jsonl" };

    try applyWorkflowToBatchConfig(&config, args, calculation, output, .{});

    try std.testing.expectEqual(OutputFormat.jsonl, config.output_format);
    try std.testing.expectEqual(true, config.residue_map);
}

test "workflow custom classifier config path resolves for batch" {
    var config = BatchConfig{};
    const args = BatchArgs{};
    const classifier_config = workflow_manifest.ClassifierConfig{ .type = "custom", .config = "my_radii.toml" };

    try applyWorkflowClassifierToBatchConfig(&config, args, classifier_config);

    try std.testing.expect(config.classifier_type == null);
    try std.testing.expectEqualStrings("my_radii.toml", config.custom_classifier_path.?);
}
```

- [ ] **Step 2: Run tests to verify failure**

Run:

```bash
zig build test
```

Expected: FAIL because `workflow_path`, `workflow_manifest`, `applyWorkflowToBatchConfig()`, `applyWorkflowClassifierToBatchConfig()`, and `custom_classifier_path` do not exist.

- [ ] **Step 3: Switch batch import and argument field**

At the top of `src/batch.zig`, replace:

```zig
const batch_manifest = @import("batch_manifest.zig");
```

with:

```zig
const workflow_manifest = @import("workflow_manifest.zig");
const classifier_parser = @import("classifier_parser.zig");
```

In `BatchArgs`, replace `manifest_path` with:

```zig
workflow_path: ?[]const u8 = null,
```

Parse both spellings into `workflow_path`:

```zig
else if (std.mem.startsWith(u8, arg, "--workflow=")) {
    result.workflow_path = arg["--workflow=".len..];
} else if (std.mem.eql(u8, arg, "--workflow")) {
    i += 1;
    if (i >= args.len) {
        std.debug.print("Error: Missing value for --workflow\n", .{});
        std.process.exit(1);
    }
    result.workflow_path = args[i];
} else if (std.mem.startsWith(u8, arg, "--manifest=")) {
    result.workflow_path = arg["--manifest=".len..];
} else if (std.mem.eql(u8, arg, "--manifest")) {
    i += 1;
    if (i >= args.len) {
        std.debug.print("Error: Missing value for --manifest\n", .{});
        std.process.exit(1);
    }
    result.workflow_path = args[i];
}
```

- [ ] **Step 4: Extend BatchConfig for custom classifier**

Add to `BatchConfig`:

```zig
custom_classifier: ?*const classifier.Classifier = null,
custom_classifier_path: ?[]const u8 = null,
```

Add a batch custom-classifier application helper near `applyBuiltinClassifier()`:

```zig
fn applyCustomClassifier(input: *AtomInput, custom_classifier: *const classifier.Classifier) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    for (0..n) |i| {
        if (custom_classifier.getRadius(residues[i].slice(), atom_names[i].slice())) |r| {
            new_radii[i] = r;
        } else if (input.element) |elements| {
            new_radii[i] = classifier.guessRadiusFromAtomicNumber(elements[i]) orelse input.r[i];
        } else {
            new_radii[i] = classifier.guessRadiusFromAtomName(atom_names[i].slice()) orelse input.r[i];
        }
    }

    input.allocator.free(input.r);
    input.r = new_radii;
}
```

Update both classifier application sites in `processOneFile()` and `processOneSdfMoleculeInner()`:

```zig
if (config.custom_classifier) |custom| {
    if (input.hasClassificationInfo()) {
        applyCustomClassifier(&input, custom) catch |err| {
            result.status = .err;
            result.error_msg = std.fmt.allocPrint(result_allocator, "custom classifier failed: {s}", .{@errorName(err)}) catch null;
            return result;
        };
    }
} else if (config.classifier_type) |ct| {
    // existing built-in classifier branch
}
```

Use `res` instead of `result` in `processOneSdfMoleculeInner()` because that function works with a local result copy.

- [ ] **Step 5: Implement workflow-to-batch config helpers**

Replace `applyManifestGlobals()` with helpers that consume workflow sections:

```zig
fn applyWorkflowToBatchConfig(
    config: *BatchConfig,
    args: BatchArgs,
    calculation: workflow_manifest.Calculation,
    output: workflow_manifest.Output,
    classifier_config: workflow_manifest.ClassifierConfig,
) !void {
    if (!args.threads_explicit) if (calculation.threads) |v| config.n_threads = v;
    if (!args.algorithm_explicit) if (calculation.algorithm) |v| config.algorithm = parseAlgorithm(v);
    if (!args.n_points_explicit) if (calculation.n_points) |v| config.n_points = validateManifestNPoints(v) catch |err| {
        std.debug.print("Error: workflow n_points must be between 1 and 10000: {d}\n", .{v});
        return err;
    };
    if (!args.n_slices_explicit) if (calculation.n_slices) |v| config.n_slices = validateManifestNSlices(v) catch |err| {
        std.debug.print("Error: workflow n_slices must be between 1 and 1000: {d}\n", .{v});
        return err;
    };
    if (!args.probe_radius_explicit) if (calculation.probe_radius) |v| config.probe_radius = validateManifestProbeRadius(v) catch |err| {
        std.debug.print("Error: workflow probe_radius must be finite and between 0 and 10 Angstroms: {d}\n", .{v});
        return err;
    };
    if (!args.precision_explicit) if (calculation.precision) |v| config.precision = parsePrecision(v);
    if (!args.format_explicit) if (output.format) |v| config.output_format = parseOutputFormat(v);
    if (!args.timing_explicit) if (calculation.timing) |v| config.show_timing = v;
    if (!args.quiet_explicit) if (calculation.quiet) |v| {
        config.quiet = v;
        config.show_progress = !v;
    };
    if (!args.include_hydrogens_explicit) if (calculation.include_hydrogens) |v| config.include_hydrogens = v;
    if (!args.include_hetatm_explicit) if (calculation.include_hetatm) |v| config.include_hetatm = v;
    if (!args.use_bitmask_explicit) if (calculation.use_bitmask) |v| config.use_bitmask = v;
    if (calculation.auth_chain) |v| config.use_auth_chain = v;
    if (calculation.residue_map) |v| config.residue_map = v;
    try applyWorkflowClassifierToBatchConfig(config, args, classifier_config);
}

fn applyWorkflowClassifierToBatchConfig(config: *BatchConfig, args: BatchArgs, c: workflow_manifest.ClassifierConfig) !void {
    const classifier_name = c.type orelse return;
    if (std.mem.eql(u8, classifier_name, "custom")) {
        if (!args.classifier_explicit) {
            config.classifier_type = null;
            config.custom_classifier_path = c.config;
        }
        return;
    }
    if (!args.classifier_explicit) config.classifier_type = parseClassifierType(classifier_name);
}
```

Keep `applyCliOverrides()` and adapt it so explicit `--classifier` clears custom classifier state:

```zig
if (args.classifier_explicit) {
    config.classifier_type = args.classifier_type;
    config.custom_classifier = null;
    config.custom_classifier_path = null;
}
```

- [ ] **Step 6: Replace manifest run path with workflow run path**

Rename `runManifest()` to `runWorkflow()` and load `workflow_manifest.Workflow`:

```zig
fn parseWorkflowFile(allocator: Allocator, io: std.Io, path: []const u8) !workflow_manifest.Workflow {
    return workflow_manifest.parseFile(allocator, io, path);
}

fn runWorkflow(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (args.chain_filter != null) {
        std.debug.print("Error: --workflow/--manifest cannot be combined with --chain; use [[jobs]].chains in the workflow\n", .{});
        return error.InvalidArgument;
    }

    const workflow_path = args.workflow_path.?;
    var workflow = parseWorkflowFile(allocator, io, workflow_path) catch |err| {
        std.debug.print("Error reading workflow '{s}': {s}\n", .{ workflow_path, @errorName(err) });
        return err;
    };
    defer workflow.deinit();

    const input_dir = args.input_path orelse workflow.input.dir orelse {
        std.debug.print("Error: Missing input directory (provide positional input_dir or [input].dir in workflow)\n", .{});
        return error.MissingArgument;
    };
    const output_dir = args.output_path orelse workflow.output.dir;

    if (workflow.jobs.len > 1 and output_dir == null) {
        std.debug.print("Error: workflow with multiple jobs requires an output directory\n", .{});
        return error.InvalidArgument;
    }

    // Load CCD, SDF, and custom classifier before iterating jobs.
}
```

In `run()`, dispatch on `workflow_path`:

```zig
if (args.workflow_path != null) return runWorkflow(allocator, io, args);
```

- [ ] **Step 7: Load workflow classifier resources**

In `runWorkflow()`, compute paths from workflow unless CLI explicitly overrides:

```zig
const load_quiet = if (args.quiet_explicit) args.quiet else (workflow.calculation.quiet orelse args.quiet);
const ccd_path = if (args.ccd_explicit) args.ccd_path else workflow.classifier.ccd;
const workflow_sdf_paths: []const []const u8 = workflow.classifier.sdf orelse &.{};
const sdf_paths = if (args.sdf_explicit) args.sdf_paths.constSlice() else workflow_sdf_paths;
```

Load custom classifier once if configured:

```zig
var custom_classifier: ?classifier.Classifier = null;
if (workflow.classifier.type) |name| {
    if (std.mem.eql(u8, name, "custom") and !args.classifier_explicit) {
        const config_path = workflow.classifier.config.?;
        custom_classifier = classifier_parser.parseConfigFile(allocator, io, config_path) catch |err| {
            switch (err) {
                error.UnsupportedConfigExtension => std.debug.print("Error loading config file '{s}': custom classifier configs are TOML-only; rename or convert the file to .toml\n", .{config_path}),
                error.UnsupportedLegacyFormat => std.debug.print("Error loading config file '{s}': FreeSASA-style custom classifier configs are no longer supported; convert to TOML [types] and [[atoms]]\n", .{config_path}),
                else => std.debug.print("Error loading config file '{s}': {s}\n", .{ config_path, @errorName(err) }),
            }
            return err;
        };
    }
}
defer if (custom_classifier) |*c| c.deinit();
```

When building each job config:

```zig
config.custom_classifier = if (custom_classifier != null) &custom_classifier.? else null;
```

- [ ] **Step 8: Remove old batch_manifest module from root**

After `src/batch.zig` compiles without `batch_manifest`, modify `src/root.zig`:

```zig
// Remove this line:
pub const batch_manifest = @import("batch_manifest.zig");

// Remove this test-block line:
_ = batch_manifest;
```

Delete the old file only if no imports remain:

```bash
rg "batch_manifest" src
```

Expected output: no matches. Then:

```bash
git rm src/batch_manifest.zig
```

- [ ] **Step 9: Run tests to verify pass**

Run:

```bash
zig build test
```

Expected: PASS. Existing batch `--manifest` tests should pass after updating expected field names from `manifest_path` to `workflow_path`.

- [ ] **Step 10: Commit**

```bash
git add src/batch.zig src/root.zig src/workflow_manifest.zig
git add -u src/batch_manifest.zig
git commit -m "feat: support batch workflow files"
```

---

### Task 5: Update command help and user docs

**Files:**
- Modify: `src/calc.zig`
- Modify: `src/batch.zig`
- Modify: `website/docs/cli/commands.md`
- Modify: `website/docs/cli/input.md`
- Modify: `website/docs/guide/classifiers.mdx`
- Modify: `README.md`
- Modify: `CHANGELOG.md`

- [ ] **Step 1: Update CLI help text**

In `src/calc.zig` help, add `--workflow` and change `--config` wording:

```text
--workflow=PATH    TOML workflow file for input, output, calculation, and classifier settings
--config=FILE      Custom classifier TOML file
```

Add example:

```text
{s} calc --workflow sasa.toml
{s} calc --config=custom.toml input.cif output.json
```

In `src/batch.zig` help, prefer workflow and document manifest alias:

```text
USAGE:
    {s} batch --workflow <workflow.toml>
    {s} batch [OPTIONS] <input_dir> [output_dir]

--workflow=PATH    TOML workflow file with one or more named batch jobs
--manifest=PATH    Compatibility alias for --workflow
```

Add examples:

```text
{s} batch --workflow bsa.toml
{s} batch structures/ results/ --workflow bsa.toml
{s} batch --manifest legacy-bsa.toml
```

- [ ] **Step 2: Update workflow docs in `website/docs/cli/commands.md`**

Replace the batch-only manifest section with a workflow section containing these examples:

```toml
version = 1
kind = "workflow"

[input]
path = "structure.cif"

[output]
path = "result.json"
format = "json"

[calculation]
algorithm = "sr"
n_points = 128
use_bitmask = true

[classifier]
type = "ccd"
```

and batch:

```toml
version = 1
kind = "workflow"

[input]
dir = "structures"

[output]
dir = "results"
format = "jsonl"

[calculation]
residue_map = true
n_points = 128

[classifier]
type = "custom"
config = "my_radii.toml"

[[jobs]]
name = "chain_A"
chains = ["A"]

[[jobs]]
name = "complex_AB"
chains = ["A", "B"]
```

State precedence exactly:

```text
Precedence is: built-in defaults < workflow settings < job settings < explicit CLI options.
```

- [ ] **Step 3: Update classifier docs for TOML-only configs**

In `website/docs/guide/classifiers.mdx` and `website/docs/cli/input.md`, remove FreeSASA-format examples. Keep this TOML example:

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }
C_CAR = { radius = 1.76, class = "apolar" }
N     = { radius = 1.65, class = "polar" }
O     = { radius = 1.40, class = "polar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"

[[atoms]]
residue = "ALA"
atom = "CB"
type = "C_ALI"
```

Add a migration subsection:

```markdown
### Migrating old FreeSASA-style configs

Old custom classifier files used this shape:

```text
name: my-classifier

types:
C_ALI 1.87 apolar

atoms:
ANY CA C_ALI
```

Convert each type line to `[types]` inline-table entries and each atom mapping to a `[[atoms]]` table:

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"
```
```

- [ ] **Step 4: Update README and changelog**

In `README.md`, mention workflow files near CLI examples:

```markdown
# Reproducible workflow file
zsasa calc --workflow sasa.toml
zsasa batch --workflow bsa.toml
```

In `CHANGELOG.md`, under the unreleased or next-version section, add:

```markdown
- Added TOML workflow files via `calc --workflow` and `batch --workflow`; `batch --manifest` remains as a compatibility alias.
- Breaking: custom classifier configs are now TOML-only. Legacy FreeSASA-style custom classifier files are no longer supported; convert them to `[types]` and `[[atoms]]` TOML format.
```

- [ ] **Step 5: Run documentation text checks**

Run:

```bash
rg -n "FreeSASA format|FreeSASA-style|--manifest|--workflow|custom\.config|my_classifier\.conf" README.md CHANGELOG.md website/docs src
```

Expected:

- `FreeSASA-style` appears only in migration/breaking-change explanations and error/help text.
- `--manifest` appears only as a compatibility alias.
- `custom.config` and `my_classifier.conf` no longer appear as recommended examples.

- [ ] **Step 6: Run tests**

Run:

```bash
zig build test
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/calc.zig src/batch.zig README.md CHANGELOG.md website/docs/cli/commands.md website/docs/cli/input.md website/docs/guide/classifiers.mdx
git commit -m "docs: document workflow files and TOML classifiers"
```

---

### Task 6: Add smoke fixtures and end-to-end verification

**Files:**
- Create: `test_data/workflow-calc.toml`
- Create: `test_data/workflow-batch.toml`
- Create: `test_data/custom-classifier.toml`
- Create: `test_data/legacy-batch-manifest.toml`
- Modify: no source files unless smoke commands reveal bugs

- [ ] **Step 1: Create smoke fixtures**

Create `test_data/custom-classifier.toml`:

```toml
name = "smoke-custom"

[types]
C_ALI = { radius = 1.87, class = "apolar" }
O     = { radius = 1.40, class = "polar" }
N     = { radius = 1.65, class = "polar" }
S     = { radius = 1.85, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"

[[atoms]]
residue = "ANY"
atom = "O"
type = "O"

[[atoms]]
residue = "ANY"
atom = "N"
type = "N"

[[atoms]]
residue = "ANY"
atom = "SG"
type = "S"
```

Create `test_data/workflow-calc.toml`:

```toml
version = 1
kind = "workflow"

[input]
path = "examples/1crn.pdb"

[output]
path = "zig-out/workflow-smoke/calc.json"
format = "json"

[calculation]
n_points = 32
quiet = true

[classifier]
type = "custom"
config = "test_data/custom-classifier.toml"
```

Create `test_data/workflow-batch.toml`:

```toml
version = 1
kind = "workflow"

[input]
dir = "examples"

[output]
dir = "zig-out/workflow-smoke/batch"
format = "jsonl"

[calculation]
n_points = 32
quiet = true
residue_map = true

[classifier]
type = "ccd"

[[jobs]]
name = "all"
```

Create `test_data/legacy-batch-manifest.toml`:

```toml
version = 1
input_dir = "examples"
output_dir = "zig-out/workflow-smoke/legacy"
format = "jsonl"
classifier = "ccd"
n_points = 32
quiet = true

[[jobs]]
name = "all"
```

- [ ] **Step 2: Run full tests**

Run:

```bash
zig build test
```

Expected: PASS.

- [ ] **Step 3: Run calc workflow smoke**

Run:

```bash
zig build run -- calc --workflow test_data/workflow-calc.toml
```

Expected: command exits 0 and creates `zig-out/workflow-smoke/calc.json`.

Verify output:

```bash
test -s zig-out/workflow-smoke/calc.json
```

Expected: exits 0.

- [ ] **Step 4: Run batch workflow smoke**

Run:

```bash
zig build run -- batch --workflow test_data/workflow-batch.toml
```

Expected: command exits 0 and creates `zig-out/workflow-smoke/batch/all.jsonl`.

Verify output:

```bash
test -s zig-out/workflow-smoke/batch/all.jsonl
```

Expected: exits 0.

- [ ] **Step 5: Run legacy manifest compatibility smoke**

Run:

```bash
zig build run -- batch --manifest test_data/legacy-batch-manifest.toml
```

Expected: command exits 0 and creates `zig-out/workflow-smoke/legacy/all.jsonl`.

Verify output:

```bash
test -s zig-out/workflow-smoke/legacy/all.jsonl
```

Expected: exits 0.

- [ ] **Step 6: Verify no old parser references remain**

Run:

```bash
rg -n "batch_manifest|parseFreeSasaConfigFile|types:|atoms:|custom\.config|my_classifier\.conf" src README.md website/docs docs/superpowers/specs/2026-05-19-workflow-manifest-design.md
```

Expected:

- No `batch_manifest` references in `src`.
- No `parseFreeSasaConfigFile` references.
- `types:` and `atoms:` appear only in migration examples that explain the removed legacy format.
- Old `.config` / `.conf` examples are absent from recommended usage.

- [ ] **Step 7: Commit fixtures**

```bash
git add test_data/workflow-calc.toml test_data/workflow-batch.toml test_data/custom-classifier.toml test_data/legacy-batch-manifest.toml
git commit -m "test: add workflow smoke fixtures"
```

---

### Task 7: Final verification and review cleanup

**Files:**
- Modify only files that verification identifies as broken.

- [ ] **Step 1: Check working tree**

Run:

```bash
git status --short
```

Expected: only known pre-existing untracked files may remain, such as `traj_sasa.csv` and `zig-pkg/`. No source/doc changes should be unstaged.

- [ ] **Step 2: Run full test suite**

Run:

```bash
zig build test
```

Expected: PASS.

- [ ] **Step 3: Run workflow smoke commands again**

Run:

```bash
zig build run -- calc --workflow test_data/workflow-calc.toml
zig build run -- batch --workflow test_data/workflow-batch.toml
zig build run -- batch --manifest test_data/legacy-batch-manifest.toml
```

Expected: all three commands exit 0.

- [ ] **Step 4: Inspect output files**

Run:

```bash
test -s zig-out/workflow-smoke/calc.json
test -s zig-out/workflow-smoke/batch/all.jsonl
test -s zig-out/workflow-smoke/legacy/all.jsonl
```

Expected: all three commands exit 0.

- [ ] **Step 5: Review diff**

Run:

```bash
git diff --stat main...HEAD
git diff main...HEAD -- src/workflow_manifest.zig src/calc.zig src/batch.zig src/classifier_parser.zig src/toml_classifier_parser.zig src/root.zig
```

Expected: diff shows only workflow parsing, CLI resolution, TOML-only classifier cleanup, docs, and fixtures. No unrelated formatting or generated-file changes.

- [ ] **Step 6: Prepare final implementation summary**

Write a concise summary for the user with:

```text
Implemented:
- calc/batch --workflow support
- batch --manifest compatibility alias
- sectioned workflow schema and legacy flat batch mapping
- TOML-only custom classifier configs
- docs and migration guidance

Verification:
- zig build test
- zig build run -- calc --workflow test_data/workflow-calc.toml
- zig build run -- batch --workflow test_data/workflow-batch.toml
- zig build run -- batch --manifest test_data/legacy-batch-manifest.toml
```

Do not claim completion until every verification command above has actually passed.

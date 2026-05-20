# Workflow File-First Batch Execution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `zsasa batch --workflow` process each input file once and reuse the parsed/classified structure across `[[jobs]]` chain selections.

**Architecture:** Keep the workflow TOML schema and public outputs unchanged. Add a workflow-specific file-first path in `src/batch.zig`: scan once, prepare shared LUTs once, parse/classify each non-SDF file once, derive job-specific `AtomInput` copies, then run the existing SASA dispatcher for each job. Preserve the existing job-first workflow path as a fallback for SDF inputs and as a reference during implementation.

**Tech Stack:** Zig 0.16, existing `src/batch.zig` batch runner, existing parsers/classifiers, existing `json_writer` output functions, `zig build test`.

---

## File Structure

- Modify: `src/batch.zig`
  - Add `copySelectedAtomInput` and small helpers for chain-filtered `AtomInput` copies.
  - Extract the calculation/output part of `processOneFile` into a reusable helper that accepts an already prepared `AtomInput`.
  - Rename the current workflow loop to `runWorkflowJobFirst` and add a new `runWorkflowFileFirst` path for non-SDF workflow directories.
  - Add focused unit tests in the existing test section.
- Modify: `website/docs/guide/workflows.md`
  - Add one short note that workflow batch jobs reuse parsed input internally and still compute each chain selection independently.
- Optionally modify: `website/docs/guide/batch.md`
  - Add one sentence pointing dimer/complex named jobs to workflow files if the workflows page note feels insufficient.

## Task 1: Add chain-selected AtomInput copy helper

**Files:**
- Modify: `src/batch.zig`

- [ ] **Step 1: Add failing tests for chain selection**

Add these tests near existing `src/batch.zig` tests. Use owned arrays so `AtomInput.deinit()` can free them.

```zig
test "copySelectedAtomInput filters one chain" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "A", "B", "A" });
    defer input.deinit();

    var selected = try copySelectedAtomInput(allocator, input, &.{"A"});
    defer selected.deinit();

    try std.testing.expectEqual(@as(usize, 2), selected.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), selected.x[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), selected.x[1], 1e-12);
    try std.testing.expectEqualStrings("A", selected.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("A", selected.chain_id.?[1].slice());
}

test "copySelectedAtomInput filters multiple chains in source order" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "A", "B", "C", "A" });
    defer input.deinit();

    var selected = try copySelectedAtomInput(allocator, input, &.{ "C", "A" });
    defer selected.deinit();

    try std.testing.expectEqual(@as(usize, 3), selected.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), selected.x[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), selected.x[1], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), selected.x[2], 1e-12);
}

test "copySelectedAtomInput null chains duplicates all atoms" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "A", "B" });
    defer input.deinit();

    var selected = try copySelectedAtomInput(allocator, input, null);
    defer selected.deinit();

    try std.testing.expectEqual(@as(usize, 2), selected.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), selected.x[0], 1e-12);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), selected.x[1], 1e-12);
}

test "copySelectedAtomInput no matching chains returns empty input" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "A", "B" });
    defer input.deinit();

    var selected = try copySelectedAtomInput(allocator, input, &.{"Z"});
    defer selected.deinit();

    try std.testing.expectEqual(@as(usize, 0), selected.atomCount());
    try std.testing.expect(selected.chain_id != null);
    try std.testing.expectEqual(@as(usize, 0), selected.chain_id.?.len);
}
```

Also add this test fixture helper above those tests:

```zig
fn makeTestAtomInput(allocator: Allocator, chains: []const []const u8) !AtomInput {
    const n = chains.len;
    const x = try allocator.alloc(f64, n);
    errdefer allocator.free(x);
    const y = try allocator.alloc(f64, n);
    errdefer allocator.free(y);
    const z = try allocator.alloc(f64, n);
    errdefer allocator.free(z);
    const r = try allocator.alloc(f64, n);
    errdefer allocator.free(r);
    const residue = try allocator.alloc(types.FixedString5, n);
    errdefer allocator.free(residue);
    const atom_name = try allocator.alloc(types.FixedString4, n);
    errdefer allocator.free(atom_name);
    const element = try allocator.alloc(u8, n);
    errdefer allocator.free(element);
    const chain_id = try allocator.alloc(types.FixedString4, n);
    errdefer allocator.free(chain_id);
    const residue_num = try allocator.alloc(i32, n);
    errdefer allocator.free(residue_num);
    const insertion_code = try allocator.alloc(types.FixedString4, n);
    errdefer allocator.free(insertion_code);

    for (chains, 0..) |chain, i| {
        x[i] = @floatFromInt(i);
        y[i] = @floatFromInt(i + 10);
        z[i] = @floatFromInt(i + 20);
        r[i] = 1.5;
        residue[i] = types.FixedString5.fromSlice("GLY");
        atom_name[i] = types.FixedString4.fromSlice("CA");
        element[i] = 6;
        chain_id[i] = types.FixedString4.fromSlice(chain);
        residue_num[i] = @intCast(i + 1);
        insertion_code[i] = types.FixedString4.fromSlice("");
    }

    return AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .residue = residue,
        .atom_name = atom_name,
        .element = element,
        .chain_id = chain_id,
        .residue_num = residue_num,
        .insertion_code = insertion_code,
        .allocator = allocator,
    };
}
```

- [ ] **Step 2: Run the new tests and verify they fail**

Run:

```bash
zig build test --summary all 2>&1 | rg "copySelectedAtomInput|makeTestAtomInput|error"
```

Expected: compile failure mentioning `copySelectedAtomInput` is undeclared.

- [ ] **Step 3: Implement the helper**

Add these helpers in `src/batch.zig` near `readInputFile` or before `processOneFile`:

```zig
fn chainMatchesFilter(chain: types.FixedString4, chains: []const []const u8) bool {
    for (chains) |target| {
        if (chain.eqlSlice(target)) return true;
    }
    return false;
}

fn atomSelectedByChains(input: AtomInput, index: usize, chains: ?[]const []const u8) bool {
    const filter = chains orelse return true;
    if (filter.len == 0) return true;
    const chain_ids = input.chain_id orelse return true;
    return chainMatchesFilter(chain_ids[index], filter);
}

fn countSelectedAtoms(input: AtomInput, chains: ?[]const []const u8) usize {
    var count: usize = 0;
    for (0..input.atomCount()) |i| {
        if (atomSelectedByChains(input, i, chains)) count += 1;
    }
    return count;
}

fn copySelectedAtomInput(allocator: Allocator, input: AtomInput, chains: ?[]const []const u8) !AtomInput {
    const selected_count = countSelectedAtoms(input, chains);

    const x = try allocator.alloc(f64, selected_count);
    errdefer allocator.free(x);
    const y = try allocator.alloc(f64, selected_count);
    errdefer allocator.free(y);
    const z = try allocator.alloc(f64, selected_count);
    errdefer allocator.free(z);
    const r = try allocator.alloc(f64, selected_count);
    errdefer allocator.free(r);

    const residue = if (input.residue != null) try allocator.alloc(types.FixedString5, selected_count) else null;
    errdefer if (residue) |v| allocator.free(v);
    const atom_name = if (input.atom_name != null) try allocator.alloc(types.FixedString4, selected_count) else null;
    errdefer if (atom_name) |v| allocator.free(v);
    const element = if (input.element != null) try allocator.alloc(u8, selected_count) else null;
    errdefer if (element) |v| allocator.free(v);
    const chain_id = if (input.chain_id != null) try allocator.alloc(types.FixedString4, selected_count) else null;
    errdefer if (chain_id) |v| allocator.free(v);
    const residue_num = if (input.residue_num != null) try allocator.alloc(i32, selected_count) else null;
    errdefer if (residue_num) |v| allocator.free(v);
    const insertion_code = if (input.insertion_code != null) try allocator.alloc(types.FixedString4, selected_count) else null;
    errdefer if (insertion_code) |v| allocator.free(v);

    var out_i: usize = 0;
    for (0..input.atomCount()) |i| {
        if (!atomSelectedByChains(input, i, chains)) continue;
        x[out_i] = input.x[i];
        y[out_i] = input.y[i];
        z[out_i] = input.z[i];
        r[out_i] = input.r[i];
        if (residue) |v| v[out_i] = input.residue.?[i];
        if (atom_name) |v| v[out_i] = input.atom_name.?[i];
        if (element) |v| v[out_i] = input.element.?[i];
        if (chain_id) |v| v[out_i] = input.chain_id.?[i];
        if (residue_num) |v| v[out_i] = input.residue_num.?[i];
        if (insertion_code) |v| v[out_i] = input.insertion_code.?[i];
        out_i += 1;
    }

    return AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .residue = residue,
        .atom_name = atom_name,
        .element = element,
        .chain_id = chain_id,
        .residue_num = residue_num,
        .insertion_code = insertion_code,
        .allocator = allocator,
    };
}
```

The `input.chain_id == null` path intentionally selects all atoms. That preserves existing JSON behavior because `readInputFile` has never applied `BatchConfig.chain_filter` to JSON inputs.

- [ ] **Step 4: Run tests for the helper**

Run:

```bash
zig build test --summary all 2>&1 | rg "copySelectedAtomInput|Build Summary|error"
```

Expected: helper tests pass and no `copySelectedAtomInput` error remains.

- [ ] **Step 5: Format and commit Task 1**

Run:

```bash
zig fmt src/batch.zig
zig build test --summary all
```

Expected: all tests pass.

Commit:

```bash
git add src/batch.zig
git commit -m "test: add workflow chain selection helper"
```

## Task 2: Extract reusable SASA calculation/output helper

**Files:**
- Modify: `src/batch.zig`

- [ ] **Step 1: Add a regression test for calculating a prepared input**

Add this test near the new helper tests:

```zig
test "calculatePreparedInputResult computes total for selected input" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "A", "A" });
    defer input.deinit();

    var result = calculatePreparedInputResult(
        f64,
        allocator,
        std.testing.io,
        allocator,
        input,
        null,
        "prepared.pdb",
        BatchConfig{ .n_points = 16, .quiet = true },
        1,
        null,
        null,
        null,
    );
    defer if (result.atom_areas) |areas| allocator.free(areas);
    defer if (result.error_msg) |msg| allocator.free(msg);

    try std.testing.expectEqual(.ok, result.status);
    try std.testing.expectEqual(@as(usize, 2), result.n_atoms);
    try std.testing.expect(result.total_sasa > 0.0);
}
```

- [ ] **Step 2: Run the test and verify it fails**

Run:

```bash
zig build test --summary all 2>&1 | rg "calculatePreparedInputResult|error"
```

Expected: compile failure mentioning `calculatePreparedInputResult` is undeclared.

- [ ] **Step 3: Extract f64/f32 calculation body into helper**

Add this generic helper above `processOneFile`:

```zig
fn calculatePreparedInputResult(
    comptime T: type,
    arena: Allocator,
    io: std.Io,
    result_allocator: Allocator,
    input: AtomInput,
    output_dir: ?[]const u8,
    filename: []const u8,
    config: BatchConfig,
    n_threads: usize,
    lut: ?*const bitmask_lut.BitmaskLutGen(T),
    coarse_lut: ?*const bitmask_lut.BitmaskLutGen(T),
    fine_lut: ?*const bitmask_lut.BitmaskLutGen(T),
) FileResult {
    var result = FileResult{
        .filename = filename,
        .n_atoms = input.atomCount(),
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .ok,
    };

    var sasa_timer = std.Io.Timestamp.now(io, .awake);
    const probe_radius: T = if (T == f64) config.probe_radius else @as(T, @floatCast(config.probe_radius));
    var sasa_result = calculateSasaDispatch(
        T,
        arena,
        input,
        config,
        probe_radius,
        n_threads,
        lut,
        coarse_lut,
        fine_lut,
    ) catch |err| {
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "SASA calculation failed: {s}", .{@errorName(err)}) catch null;
        return result;
    };
    defer sasa_result.deinit();
    result.sasa_time_ns = @intCast(sasa_timer.untilNow(io, .awake).nanoseconds);
    result.total_sasa = if (T == f64) sasa_result.total_area else @as(f64, @floatCast(sasa_result.total_area));

    if (config.store_atom_areas) {
        if (T == f64) {
            result.atom_areas = arena.dupe(f64, sasa_result.atom_areas) catch {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                return result;
            };
        } else {
            const areas_f64 = arena.alloc(f64, sasa_result.atom_areas.len) catch {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                return result;
            };
            for (sasa_result.atom_areas, 0..) |v, j| areas_f64[j] = @floatCast(v);
            result.atom_areas = areas_f64;
        }
    }

    if (config.residue_map) {
        const areas = result.atom_areas orelse blk: {
            if (T == f64) break :blk sasa_result.atom_areas;
            const areas_f64 = arena.alloc(f64, sasa_result.atom_areas.len) catch {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                return result;
            };
            for (sasa_result.atom_areas, 0..) |v, j| areas_f64[j] = @floatCast(v);
            break :blk areas_f64;
        };
        if (!attachResidueMap(arena, result_allocator, input, areas, &result)) return result;
    }

    if (output_dir) |out_dir| {
        if (!config.store_atom_areas) {
            writeSasaOutput(T, arena, io, &sasa_result, out_dir, filename, config.output_format) catch |err| {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "output write failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
        }
    }

    return result;
}
```

Then replace the duplicated f64 and f32 blocks in `processOneFile` with calls to `calculatePreparedInputResult`. Preserve the parse and classifier code above it.

- [ ] **Step 4: Run tests**

Run:

```bash
zig fmt src/batch.zig
zig build test --summary all
```

Expected: all tests pass.

- [ ] **Step 5: Commit Task 2**

```bash
git add src/batch.zig
git commit -m "refactor: reuse prepared batch SASA calculation"
```

## Task 3: Add workflow file-first runner for non-SDF inputs

**Files:**
- Modify: `src/batch.zig`

- [ ] **Step 1: Preserve old workflow runner under a new name**

Rename the current `runWorkflow` function to `runWorkflowJobFirst` without changing its body. Then add a new small dispatcher named `runWorkflow`:

```zig
fn runWorkflow(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (try workflowInputContainsSdf(allocator, io, args)) {
        return runWorkflowJobFirst(allocator, io, args);
    }
    return runWorkflowFileFirst(allocator, io, args);
}
```

Add the new helper declaration with a temporary body that delegates to job-first:

```zig
fn runWorkflowFileFirst(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    return runWorkflowJobFirst(allocator, io, args);
}
```

- [ ] **Step 2: Add SDF fallback detector**

Add:

```zig
fn workflowInputContainsSdf(allocator: Allocator, io: std.Io, args: BatchArgs) !bool {
    const workflow_path = args.workflow_path orelse return false;
    var workflow = try parseWorkflowFile(allocator, io, workflow_path);
    defer workflow.deinit();
    const input_dir = args.input_path orelse workflow.input.dir orelse return false;
    const files = try scanDirectory(allocator, io, input_dir);
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }
    for (files) |filename| {
        if (format_detect.detectInputFormat(filename) == .sdf) return true;
    }
    return false;
}
```

This helper can parse the workflow once before the real runner parses it again. That overhead only happens once and keeps the fallback simple.

- [ ] **Step 3: Run tests to ensure rename is safe**

Run:

```bash
zig fmt src/batch.zig
zig build test --summary all
```

Expected: all tests pass.

- [ ] **Step 4: Implement workflow result accumulation types**

Add these small structs near `BatchResult`:

```zig
const WorkflowJobState = struct {
    name: []const u8,
    config: BatchConfig,
    output_dir: ?[]const u8 = null,
    jsonl_output_path: ?[]const u8 = null,
    successful: usize = 0,
    failed: usize = 0,
    total_sasa_time_ns: u64 = 0,

    fn deinit(self: *WorkflowJobState, allocator: Allocator) void {
        if (self.output_dir) |path| allocator.free(path);
        if (self.jsonl_output_path) |path| allocator.free(path);
    }
};
```

Add a JSONL write helper that writes one result to a path by creating/appending through the existing JSONL serializer. Use truncation only when initializing each job output path, not for every result.

```zig
fn truncateJsonlOutput(io: std.Io, path: []const u8) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    file.close(io);
}

fn appendJsonlResultToFile(io: std.Io, file: std.Io.File, allocator: Allocator, result: *FileResult) !void {
    var write_buf: [64 * 1024]u8 = undefined;
    var writer = std.Io.File.Writer.initStreaming(file, io, &write_buf);
    const line = try fileResultToJsonlLine(allocator, result);
    defer allocator.free(line);
    try writer.interface.writeAll(line);
    try writer.interface.writeByte('\n');
    try writer.interface.flush();
}
```

If Zig's `std.Io.File.Writer.initStreaming` cannot be re-created safely for the same file handle, replace this helper with a `WorkflowJsonlWriter` struct that owns `file`, `buffer`, and `writer` for the full workflow run.

- [ ] **Step 5: Implement shared workflow setup in `runWorkflowFileFirst`**

Replace the temporary body with setup equivalent to `runWorkflowJobFirst`:

```zig
fn runWorkflowFileFirst(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (args.chain_filter != null) {
        std.debug.print("Error: --workflow cannot be combined with --chain; --manifest is a compatibility alias for --workflow; use [[jobs]].chains in the workflow\n", .{});
        return error.InvalidArgument;
    }

    const workflow_path = args.workflow_path.?;
    var workflow = parseWorkflowFile(allocator, io, workflow_path) catch |err| {
        std.debug.print("Error reading workflow file '{s}': {s}\n", .{ workflow_path, @errorName(err) });
        return err;
    };
    defer workflow.deinit();

    if (workflow.jobs.len == 0) {
        std.debug.print("Error: batch workflow requires at least one [[jobs]] entry\n", .{});
        return error.NoJobs;
    }

    const input_dir = args.input_path orelse workflow.input.dir orelse {
        std.debug.print("Error: Missing input directory (provide positional input_dir or [input].dir in workflow)\n", .{});
        return error.MissingArgument;
    };
    const output_dir = args.output_path orelse workflow.output.dir;
    if (workflow.jobs.len > 1 and output_dir == null) {
        std.debug.print("Error: workflow with multiple jobs requires an output directory\n", .{});
        return error.InvalidArgument;
    }

    // Copy resource loading logic from runWorkflowJobFirst exactly:
    // load_quiet, resource_config, ext_ccd, sdf_ccd, custom_classifier.
    // Then build WorkflowJobState entries with per-job config and output paths.
}
```

When building `WorkflowJobState` entries, use the same rules as `runWorkflowJobFirst`:

```zig
var states = try allocator.alloc(WorkflowJobState, workflow.jobs.len);
errdefer allocator.free(states);
for (workflow.jobs, 0..) |job, i| {
    var config = BatchConfig{};
    try applyWorkflowToBatchConfig(&config, args, workflow.calculation, workflow.output, workflow.classifier);
    applyCliOverrides(&config, args);
    config.include_hetatm = config.include_hetatm or (config.classifier_type == .ccd);
    config.store_atom_areas = (config.output_format == .jsonl);
    config.external_ccd = if (ext_ccd != null) &ext_ccd.? else null;
    config.sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null;
    if (custom_classifier) |*c| config.custom_classifier = c;
    applyWorkflowJobOverrides(&config, args, job);
    try validateResidueMapFormat(config.output_format, config.residue_map);

    var job_output_dir: ?[]const u8 = null;
    var jsonl_output_path: ?[]const u8 = null;
    if (config.output_format == .jsonl) {
        if (output_dir) |out| {
            try std.Io.Dir.cwd().createDirPath(io, out);
            jsonl_output_path = try workflowJsonlOutputPath(allocator, out, job.name);
            try truncateJsonlOutput(io, jsonl_output_path.?);
        }
    } else if (output_dir) |out| {
        job_output_dir = try workflowPerFileOutputDir(allocator, out, job.name);
        try std.Io.Dir.cwd().createDirPath(io, job_output_dir.?);
    }

    states[i] = .{
        .name = job.name,
        .config = config,
        .output_dir = job_output_dir,
        .jsonl_output_path = jsonl_output_path,
    };
}
defer {
    for (states) |*state| state.deinit(allocator);
    allocator.free(states);
}
```

- [ ] **Step 6: Add the file-first processing loop**

Inside `runWorkflowFileFirst`, after state setup:

```zig
const files = try scanDirectory(allocator, io, input_dir);
defer {
    for (files) |f| allocator.free(f);
    allocator.free(files);
}

var luts = try BatchLuts.init(allocator, resource_config);
defer luts.deinit();

var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
defer arena.deinit();

for (files) |filename| {
    const input_path = try std.fs.path.join(arena.allocator(), &.{ input_dir, filename });

    var source_config = resource_config;
    source_config.chain_filter = null;
    var source_input = readInputFile(arena.allocator(), io, input_path, source_config) catch |err| {
        for (states) |*state| {
            state.failed += 1;
            std.debug.print("Error running workflow job '{s}' on '{s}': read/parse failed: {s}\n", .{ state.name, filename, @errorName(err) });
        }
        _ = arena.reset(.retain_capacity);
        continue;
    };
    defer source_input.deinit();

    if (source_config.custom_classifier) |custom_classifier| {
        if (source_input.hasClassificationInfo()) try applyCustomClassifier(&source_input, custom_classifier, source_config.quiet);
    } else if (source_config.classifier_type) |ct| {
        const format = format_detect.detectInputFormat(input_path);
        if (format != .json and source_input.hasClassificationInfo()) try applyBuiltinClassifier(&source_input, ct, source_config.sdf_ccd, source_config.external_ccd);
    }

    for (workflow.jobs, 0..) |job, job_index| {
        var state = &states[job_index];
        if (!state.config.quiet) std.debug.print("Workflow job: {s}\n", .{state.name});

        const selected_chains: ?[]const []const u8 = if (format_detect.detectInputFormat(input_path) == .json) null else job.chains;
        var selected_input = copySelectedAtomInput(arena.allocator(), source_input, selected_chains) catch |err| {
            state.failed += 1;
            std.debug.print("Error running workflow job '{s}' on '{s}': selection failed: {s}\n", .{ state.name, filename, @errorName(err) });
            continue;
        };
        defer selected_input.deinit();

        var result = switch (state.config.precision) {
            .f64 => calculatePreparedInputResult(f64, arena.allocator(), io, allocator, selected_input, state.output_dir, filename, state.config, 1, luts.f64Ptr(), luts.coarseF64Ptr(), luts.fineF64Ptr()),
            .f32 => calculatePreparedInputResult(f32, arena.allocator(), io, allocator, selected_input, state.output_dir, filename, state.config, 1, luts.f32Ptr(), luts.coarseF32Ptr(), luts.fineF32Ptr()),
        };

        if (result.status == .ok) {
            state.successful += 1;
            state.total_sasa_time_ns += result.sasa_time_ns;
        } else {
            state.failed += 1;
        }

        if (state.config.output_format == .jsonl) {
            if (state.jsonl_output_path) |path| {
                // Prefer a persistent file handle if implemented in Step 4.
                const file = try std.Io.Dir.cwd().openFile(io, path, .{});
                defer file.close(io);
                try appendJsonlResultToFile(io, file, arena.allocator(), &result);
            } else {
                var stdout_writer = std.Io.File.Writer.initStreaming(std.Io.File.stdout(), io, &[_]u8{});
                writeJsonlResult(&stdout_writer, arena.allocator(), &result);
            }
        }

        result.atom_areas = null;
        result.residue_map = null;
    }

    _ = arena.reset(.retain_capacity);
}

var successful: usize = 0;
var failed: usize = 0;
for (states) |state| {
    successful += state.successful;
    failed += state.failed;
}
std.debug.print("Workflow complete: {d} successful, {d} failed\n", .{ successful, failed });
```

Adjust the JSONL file handling to the exact Zig file open flags available in this repository. The required behavior is: create/truncate each job JSONL file once, then append one line per successful or failed file result in deterministic scan order.

- [ ] **Step 7: Run tests and fix compile errors**

Run:

```bash
zig fmt src/batch.zig
zig build test --summary all
```

Expected: all tests pass. If the JSONL helper requires persistent writer state, implement `WorkflowJsonlWriter` before proceeding.

- [ ] **Step 8: Commit Task 3**

```bash
git add src/batch.zig
git commit -m "feat: process workflow batch files once"
```

## Task 4: Add output compatibility tests for workflow jobs

**Files:**
- Modify: `src/batch.zig`
- Create test fixtures under `test_data/` only if existing examples are insufficient.

- [ ] **Step 1: Add a focused test for job state output path generation**

Add this unit test if no direct file-first integration test is practical:

```zig
test "workflow job state keeps existing output layout" {
    const allocator = std.testing.allocator;

    const jsonl_path = try workflowJsonlOutputPath(allocator, "results", "chain_a");
    defer allocator.free(jsonl_path);
    try std.testing.expectEqualStrings("results/chain_a.jsonl", jsonl_path);

    const per_file_dir = try workflowPerFileOutputDir(allocator, "results", "complex_ab");
    defer allocator.free(per_file_dir);
    try std.testing.expectEqualStrings("results/complex_ab", per_file_dir);
}
```

- [ ] **Step 2: Add a file-first vs job-first smoke script command to manual verification**

Use existing workflow fixture first:

```bash
rm -rf zig-out/workflow-smoke
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa batch --workflow test_data/workflow-batch.toml
ls zig-out/workflow-smoke/batch
```

Expected: command succeeds and creates `zig-out/workflow-smoke/batch/all.jsonl` for the existing fixture.

- [ ] **Step 3: Add a temporary local dimer workflow for manual comparison**

Create `/tmp/zsasa-workflow-ab.toml` during verification, not committed:

```toml
version = 1
kind = "workflow"

[input]
dir = "examples"

[output]
dir = "zig-out/workflow-ab"
format = "jsonl"

[calculation]
n_points = 32
quiet = true
residue_map = true

[classifier]
type = "ccd"

[[jobs]]
name = "chain_a"
chains = ["A"]

[[jobs]]
name = "chain_b"
chains = ["B"]

[[jobs]]
name = "complex_ab"
chains = ["A", "B"]
```

Run:

```bash
rm -rf zig-out/workflow-ab
./zig-out/bin/zsasa batch --workflow /tmp/zsasa-workflow-ab.toml
find zig-out/workflow-ab -type f -maxdepth 1 -print | sort
```

Expected: `chain_a.jsonl`, `chain_b.jsonl`, and `complex_ab.jsonl` exist. If the examples directory has no B chain, `chain_b` and `complex_ab` should still behave consistently with the previous chain-filtered parser behavior for the fixture.

- [ ] **Step 4: Run full test suite**

```bash
zig fmt --check src/
zig build test --summary all
```

Expected: both commands pass.

- [ ] **Step 5: Commit Task 4**

```bash
git add src/batch.zig
git commit -m "test: cover workflow file-first outputs"
```

## Task 5: Document workflow execution reuse

**Files:**
- Modify: `website/docs/guide/workflows.md`

- [ ] **Step 1: Add documentation note**

In `website/docs/guide/workflows.md`, after the batch workflow example, add:

```markdown
Workflow batch runs reuse parsed structures across jobs internally. For named chain analyses such as chain A, chain B, and complex AB, list only the jobs you want; zsasa will parse each input structure once and then compute each requested chain selection independently.
```

- [ ] **Step 2: Run docs-adjacent check**

Run:

```bash
rg -n "reuse parsed structures|chain A, chain B" website/docs/guide/workflows.md
```

Expected: the new note appears exactly once.

- [ ] **Step 3: Commit docs**

```bash
git add website/docs/guide/workflows.md
git commit -m "docs: note workflow job reuse"
```

## Task 6: Final verification and cleanup

**Files:**
- No new files unless previous tasks require test fixtures.

- [ ] **Step 1: Run required verification**

```bash
zig fmt --check src/
zig build test --summary all
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa batch --workflow test_data/workflow-batch.toml
```

Expected: all commands pass. The workflow command prints `Workflow complete:` and returns exit code 0.

- [ ] **Step 2: Inspect git state**

```bash
git status --short --branch
git log --oneline --decorate -6
```

Expected: branch is `perf/workflow-file-first`; working tree is clean after commits.

- [ ] **Step 3: Prepare final summary**

Summarize:

- Behavior unchanged for workflow schema and outputs.
- `batch --workflow` now uses file-first processing for non-SDF workflow inputs.
- SDF workflows use the preserved job-first fallback.
- Verification commands and results.

## Self-Review

- Spec coverage: The plan keeps `[[jobs]]`, preserves output layout, shares scan/LUT/parse/classifier work, keeps SASA per job, avoids `[[complexes]]`, and includes tests/docs.
- Placeholder scan: No placeholder markers or unspecified edge-case steps remain.
- Type consistency: Helper names used consistently: `copySelectedAtomInput`, `calculatePreparedInputResult`, `runWorkflowJobFirst`, `runWorkflowFileFirst`, `WorkflowJobState`.

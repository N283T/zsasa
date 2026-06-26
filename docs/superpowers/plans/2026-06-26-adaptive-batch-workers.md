# Adaptive Batch Workers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in `zsasa batch --adaptive-workers` mode that calibrates file-worker count on a single machine and uses the smallest near-best worker count for the remaining batch run.

**Architecture:** Keep batch mode file-level parallelism. Refactor the existing parallel runner so a fixed worker count can process a contiguous work-item range, then use that primitive for calibration windows and the final main range. Add small pure helper functions for candidate generation and worker selection so scheduler policy is unit-testable without large datasets.

**Tech Stack:** Zig 0.16, existing `src/batch.zig` batch runner, Docusaurus Markdown docs in `website/docs/`.

---

## File Structure

- Modify `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`
  - Add CLI parsing for `--adaptive-workers`.
  - Add `BatchConfig.adaptive_workers` and `BatchArgs.adaptive_workers` fields.
  - Add pure scheduler helpers and tests.
  - Refactor parallel work execution into a reusable range helper.
  - Add adaptive calibration path inside `runBatchParallel`.
- Modify `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/batch.md`
  - Document experimental adaptive worker mode for I/O-bound large directory runs.
- Modify `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/cli/commands.md`
  - Add `--adaptive-workers` to batch options.

No C API or Python API changes are planned in the first implementation. The C API will keep default `adaptive_workers = false` through the new `BatchConfig` default field.

---

### Task 1: Add scheduler policy helpers and unit tests

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`

- [ ] **Step 1: Write failing tests for adaptive worker candidates and selection**

Add these tests near the existing batch tests in `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`:

```zig
test "adaptive worker candidates include powers of two and max" {
    const allocator = std.testing.allocator;

    const one = try adaptiveWorkerCandidates(allocator, 1);
    defer allocator.free(one);
    try std.testing.expectEqualSlices(usize, &.{1}, one);

    const two = try adaptiveWorkerCandidates(allocator, 2);
    defer allocator.free(two);
    try std.testing.expectEqualSlices(usize, &.{ 1, 2 }, two);

    const ten = try adaptiveWorkerCandidates(allocator, 10);
    defer allocator.free(ten);
    try std.testing.expectEqualSlices(usize, &.{ 1, 2, 4, 8, 10 }, ten);
}

test "adaptive worker selection chooses smallest near best throughput" {
    const samples = [_]AdaptiveWorkerSample{
        .{ .workers = 1, .items = 100, .elapsed_ns = 100_000_000 },
        .{ .workers = 2, .items = 100, .elapsed_ns = 60_000_000 },
        .{ .workers = 4, .items = 100, .elapsed_ns = 40_000_000 },
        .{ .workers = 8, .items = 100, .elapsed_ns = 39_000_000 },
    };
    try std.testing.expectEqual(@as(usize, 4), chooseAdaptiveWorkerCount(&samples));
}

test "adaptive calibration plan falls back for small datasets" {
    const candidates_len: usize = 5;
    try std.testing.expectEqual(@as(?usize, null), adaptiveCalibrationWindowSize(100, candidates_len));
    try std.testing.expectEqual(@as(?usize, 64), adaptiveCalibrationWindowSize(640, candidates_len));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
zig build test --summary all 2>&1 | rg "adaptive worker|error:|FAIL"
```

Expected: compile/test failure because `adaptiveWorkerCandidates`, `AdaptiveWorkerSample`, `chooseAdaptiveWorkerCount`, and `adaptiveCalibrationWindowSize` do not exist yet.

- [ ] **Step 3: Add minimal scheduler helper implementation**

Add this code after `BatchLuts` in `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`:

```zig
const adaptive_worker_min_items: usize = 64;
const adaptive_worker_max_items: usize = 256;
const adaptive_worker_near_best_ratio: f64 = 0.95;

const AdaptiveWorkerSample = struct {
    workers: usize,
    items: usize,
    elapsed_ns: u64,

    fn throughput(self: AdaptiveWorkerSample) f64 {
        if (self.items == 0 or self.elapsed_ns == 0) return 0.0;
        const seconds = @as(f64, @floatFromInt(self.elapsed_ns)) / 1_000_000_000.0;
        return @as(f64, @floatFromInt(self.items)) / seconds;
    }
};

fn adaptiveWorkerCandidates(allocator: Allocator, max_workers: usize) ![]usize {
    const capped_max = @max(max_workers, 1);
    var candidates = std.ArrayListUnmanaged(usize).empty;
    errdefer candidates.deinit(allocator);

    var workers: usize = 1;
    while (workers < capped_max) : (workers *= 2) {
        try candidates.append(allocator, workers);
    }
    if (candidates.items.len == 0 or candidates.items[candidates.items.len - 1] != capped_max) {
        try candidates.append(allocator, capped_max);
    }

    return candidates.toOwnedSlice(allocator);
}

fn adaptiveCalibrationWindowSize(total_items: usize, candidate_count: usize) ?usize {
    if (candidate_count == 0) return null;
    const required = candidate_count * adaptive_worker_min_items + adaptive_worker_min_items;
    if (total_items < required) return null;

    const per_candidate = total_items / (candidate_count + 1);
    return @min(adaptive_worker_max_items, @max(adaptive_worker_min_items, per_candidate));
}

fn chooseAdaptiveWorkerCount(samples: []const AdaptiveWorkerSample) usize {
    if (samples.len == 0) return 1;

    var best_throughput: f64 = 0.0;
    for (samples) |sample| {
        best_throughput = @max(best_throughput, sample.throughput());
    }

    const threshold = best_throughput * adaptive_worker_near_best_ratio;
    for (samples) |sample| {
        if (sample.throughput() >= threshold) return sample.workers;
    }

    return samples[samples.len - 1].workers;
}
```

- [ ] **Step 4: Run tests to verify helpers pass**

Run:

```bash
zig build test --summary all 2>&1 | rg "adaptive worker|Build Summary|error:|FAIL"
```

Expected: no adaptive helper failures. If the full command prints no `error:` or `FAIL` lines and `Build Summary` reports success, continue.

- [ ] **Step 5: Commit scheduler helpers**

Run:

```bash
git add /Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig
git commit -m "test: cover adaptive batch worker policy"
```

---

### Task 2: Add CLI and config plumbing for `--adaptive-workers`

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`

- [ ] **Step 1: Write failing parse test**

Add this test next to the other `BatchArgs` parser tests in `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`:

```zig
test "BatchArgs --adaptive-workers" {
    const args = [_][]const u8{ "zsasa", "batch", "--adaptive-workers", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.adaptive_workers);
    try std.testing.expectEqual(true, parsed.adaptive_workers_explicit);
}
```

- [ ] **Step 2: Run parser test to verify it fails**

Run:

```bash
zig build test --summary all 2>&1 | rg "adaptive-workers|adaptive_workers|error:|FAIL"
```

Expected: compile/test failure because the fields and parser branch do not exist.

- [ ] **Step 3: Add config and args fields**

In `BatchConfig`, add this field after `adaptive_high`:

```zig
    adaptive_workers: bool = false, // Experimental worker-count calibration for I/O-bound batch runs
```

In `BatchArgs`, add this field after `adaptive_high`:

```zig
    adaptive_workers: bool = false,
```

In `BatchArgs`, add this explicit flag after `adaptive_high_explicit`:

```zig
    adaptive_workers_explicit: bool = false,
```

- [ ] **Step 4: Parse `--adaptive-workers`**

In `parseArgs`, add this branch after the `--adaptive-sr` branch and before `--coarse-points` handling:

```zig
        } else if (std.mem.eql(u8, arg, "--adaptive-workers")) {
            result.adaptive_workers = true;
            result.adaptive_workers_explicit = true;
```

- [ ] **Step 5: Add workflow guard and config assignment**

In `run`, extend the workflow guard so workflow mode rejects the new option:

```zig
    if (args.workflow_path != null) {
        if (args.adaptive_sr_explicit or args.coarse_points_explicit or args.fine_points_explicit or args.adaptive_low_explicit or args.adaptive_high_explicit) {
            std.debug.print("Error: adaptive SR options are not supported with --workflow yet\n", .{});
            return error.InvalidArgument;
        }
        if (args.adaptive_workers_explicit) {
            std.debug.print("Error: --adaptive-workers is not supported with --workflow yet\n", .{});
            return error.InvalidArgument;
        }
        return runWorkflow(allocator, io, args);
    }
```

In the `BatchConfig` literal in `run`, add:

```zig
        .adaptive_workers = args.adaptive_workers,
```

- [ ] **Step 6: Update CLI help text**

In `printHelp`, add this option near `--threads`:

```zig
        \\    --adaptive-workers  Experimental: calibrate file-worker count for I/O-bound large batch runs
```

- [ ] **Step 7: Run parser tests**

Run:

```bash
zig build test --summary all 2>&1 | rg "BatchArgs --adaptive-workers|adaptive-workers|Build Summary|error:|FAIL"
```

Expected: parser test passes and no compile errors.

- [ ] **Step 8: Commit CLI plumbing**

Run:

```bash
git add /Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig
git commit -m "feat: add adaptive workers batch flag"
```

---

### Task 3: Refactor parallel execution into a reusable range helper

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`

- [ ] **Step 1: Add a named build-work-items result type**

Add this type immediately before `buildWorkItems`:

```zig
const BuildWorkItemsResult = struct {
    items: []WorkItem,
    sdf_sources: std.StringHashMapUnmanaged([]const u8),
};
```

Change the `buildWorkItems` signature from:

```zig
) !struct { items: []WorkItem, sdf_sources: std.StringHashMapUnmanaged([]const u8) } {
```

into:

```zig
) !BuildWorkItemsResult {
```

- [ ] **Step 2: Run baseline JSONL parallel test before refactor**

Run:

```bash
zig build test --summary all 2>&1 | rg "runBatchParallel writes parseable JSONL|Build Summary|error:|FAIL"
```

Expected: existing `runBatchParallel writes parseable JSONL with multiple threads` test passes before the refactor.

- [ ] **Step 3: Change `ParallelContext.processed_count` to a pointer**

Replace this field in `ParallelContext`:

```zig
    processed_count: std.atomic.Value(usize),
```

with:

```zig
    processed_count: *std.atomic.Value(usize),
```

No worker code changes are needed because Zig pointer field access still supports `ctx.processed_count.fetchAdd(...)` and `ctx.processed_count.load(...)`.

- [ ] **Step 4: Add range execution helper**

Add this function before `runBatchParallel`:

```zig
fn runWorkItemRangeParallel(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    work_items: []const WorkItem,
    file_results: []FileResult,
    build_result: *const BuildWorkItemsResult,
    luts: *const BatchLuts,
    jsonl_stream_ptr: ?*JsonlStreamWriter,
    n_workers: usize,
    processed_count: *std.atomic.Value(usize),
    expected_processed: usize,
    progress_node: ?*std.Progress.Node,
) !void {
    if (work_items.len == 0) return;

    var ctx = ParallelContext{
        .work_items = work_items,
        .input_dir = input_dir,
        .output_dir = output_dir,
        .config = config,
        .results = file_results,
        .result_allocator = allocator,
        .next_item = std.atomic.Value(usize).init(0),
        .processed_count = processed_count,
        .lut_f64 = luts.f64Ptr(),
        .lut_f32 = luts.f32Ptr(),
        .coarse_lut_f64 = luts.coarseF64Ptr(),
        .fine_lut_f64 = luts.fineF64Ptr(),
        .coarse_lut_f32 = luts.coarseF32Ptr(),
        .fine_lut_f32 = luts.fineF32Ptr(),
        .jsonl_stream = jsonl_stream_ptr,
        .io = io,
        .sdf_sources = build_result.sdf_sources,
    };

    const actual_workers = @min(@max(n_workers, 1), work_items.len);
    const threads = try allocator.alloc(std.Thread, actual_workers);
    defer allocator.free(threads);

    for (threads) |*thread| {
        thread.* = try std.Thread.spawn(.{}, parallelWorker, .{&ctx});
    }

    if (progress_node) |node| {
        while (processed_count.load(.acquire) < expected_processed) {
            node.setCompletedItems(processed_count.load(.acquire));
            std.Io.sleep(io, .fromMilliseconds(50), .awake) catch {};
        }
        node.setCompletedItems(expected_processed);
    }

    for (threads) |thread| {
        thread.join();
    }
}
```

- [ ] **Step 5: Replace duplicated thread spawning in `runBatchParallel`**

Inside `runBatchParallel`, keep scanning, work-item building, result allocation, thread-count selection, LUT creation, and JSONL setup. Replace the existing `ParallelContext` creation, thread spawning, progress loop, and join block with:

```zig
    var processed_count = std.atomic.Value(usize).init(0);

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing items", .estimated_total_items = work_items.len })
    else
        .none;
    defer progress_root.end();

    try runWorkItemRangeParallel(
        allocator,
        io,
        input_dir,
        output_dir,
        config,
        work_items,
        file_results,
        &build_result,
        &luts,
        jsonl_stream_ptr,
        actual_threads,
        &processed_count,
        work_items.len,
        if (shouldShowProgress(config)) &progress_root else null,
    );

    if (shouldShowProgress(config)) {
        progress_root.setCompletedItems(work_items.len);
    }
```

- [ ] **Step 6: Run regression tests**

Run:

```bash
zig build test --summary all 2>&1 | rg "runBatchParallel writes parseable JSONL|Build Summary|error:|FAIL"
```

Expected: existing JSONL parallel test still passes and no compile errors.

- [ ] **Step 7: Commit refactor**

Run:

```bash
git add /Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig
git commit -m "refactor: isolate batch work-item range execution"
```

---

### Task 4: Implement adaptive calibration inside `runBatchParallel`

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`

- [ ] **Step 1: Write integration test for adaptive JSONL batch**

Add this test after `runBatchParallel writes parseable JSONL with multiple threads`:

```zig
test "runBatchParallel adaptive workers writes parseable JSONL" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_path = try std.fs.path.join(allocator, &.{ root_path, "adaptive-results.jsonl" });
    defer allocator.free(output_path);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const pdb_data =
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n" ++
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C\n" ++
        "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00 20.00           C\n" ++
        "END\n";

    var name_buf: [32]u8 = undefined;
    for (0..10) |i| {
        const filename = try std.fmt.bufPrint(&name_buf, "tiny-{d}.pdb", .{i});
        const path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = path, .data = pdb_data });
    }

    var result = try runBatchParallel(allocator, std.testing.io, input_dir, null, .{
        .n_threads = 4,
        .adaptive_workers = true,
        .algorithm = .sr,
        .n_points = 8,
        .quiet = true,
        .show_progress = false,
        .output_format = .jsonl,
        .store_atom_areas = true,
        .classifier_type = .naccess,
    }, output_path);
    defer result.deinit();

    try std.testing.expectEqual(@as(usize, 10), result.total_files);
    try std.testing.expectEqual(@as(usize, 10), result.successful);
    try std.testing.expectEqual(@as(usize, 0), result.failed);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(64 * 1024));
    defer allocator.free(content);
    try std.testing.expectEqual(@as(usize, 10), std.mem.count(u8, content, "\n"));
}
```

This small test exercises the adaptive flag and the small-dataset fallback path.

- [ ] **Step 2: Run integration test to verify it fails before implementation**

Run:

```bash
zig build test --summary all 2>&1 | rg "adaptive workers writes parseable|adaptive_workers|error:|FAIL"
```

Expected: failure if Task 2 was not completed, or pass-through without adaptive behavior if the flag exists but no adaptive branch exists. If it passes because fallback behavior already routes normally, continue and rely on the next adaptive branch code plus helper tests for policy coverage.

- [ ] **Step 3: Add adaptive phase runner**

Add this function before `runBatchParallel`:

```zig
fn runAdaptiveWorkItemsParallel(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    work_items: []const WorkItem,
    file_results: []FileResult,
    build_result: *const BuildWorkItemsResult,
    luts: *const BatchLuts,
    jsonl_stream_ptr: ?*JsonlStreamWriter,
    max_workers: usize,
    progress_node: ?*std.Progress.Node,
) !void {
    const candidates = try adaptiveWorkerCandidates(allocator, max_workers);
    defer allocator.free(candidates);

    const maybe_window = adaptiveCalibrationWindowSize(work_items.len, candidates.len);
    if (maybe_window == null) {
        if (!config.quiet) {
            std.debug.print("Adaptive workers: dataset too small for calibration; using {d} workers\n", .{max_workers});
        }
        var processed_count = std.atomic.Value(usize).init(0);
        return runWorkItemRangeParallel(allocator, io, input_dir, output_dir, config, work_items, file_results, build_result, luts, jsonl_stream_ptr, max_workers, &processed_count, work_items.len, progress_node);
    }
    const window_size = maybe_window.?;

    var processed_count = std.atomic.Value(usize).init(0);
    var offset: usize = 0;
    var samples = try allocator.alloc(AdaptiveWorkerSample, candidates.len);
    defer allocator.free(samples);

    for (candidates, 0..) |workers, i| {
        const end = @min(offset + window_size, work_items.len);
        var timer = std.Io.Timestamp.now(io, .awake);
        try runWorkItemRangeParallel(allocator, io, input_dir, output_dir, config, work_items[offset..end], file_results[offset..end], build_result, luts, jsonl_stream_ptr, workers, &processed_count, end, progress_node);
        samples[i] = .{
            .workers = workers,
            .items = end - offset,
            .elapsed_ns = @intCast(timer.untilNow(io, .awake).nanoseconds),
        };
        offset = end;
    }

    const selected_workers = chooseAdaptiveWorkerCount(samples);
    if (!config.quiet) {
        std.debug.print("Adaptive workers: selected {d} of {d} workers\n", .{ selected_workers, max_workers });
        for (samples) |sample| {
            std.debug.print("  {d} workers: {d:.1} files/sec over {d} items\n", .{ sample.workers, sample.throughput(), sample.items });
        }
    }

    if (offset < work_items.len) {
        try runWorkItemRangeParallel(allocator, io, input_dir, output_dir, config, work_items[offset..], file_results[offset..], build_result, luts, jsonl_stream_ptr, selected_workers, &processed_count, work_items.len, progress_node);
    }
}
```

- [ ] **Step 4: Wire adaptive mode into `runBatchParallel`**

In `runBatchParallel`, after JSONL stream setup and before aggregate results, route execution like this:

```zig
    if (config.adaptive_workers) {
        try runAdaptiveWorkItemsParallel(
            allocator,
            io,
            input_dir,
            output_dir,
            config,
            work_items,
            file_results,
            &build_result,
            &luts,
            jsonl_stream_ptr,
            actual_threads,
            if (shouldShowProgress(config)) &progress_root else null,
        );
    } else {
        var processed_count = std.atomic.Value(usize).init(0);
        try runWorkItemRangeParallel(
            allocator,
            io,
            input_dir,
            output_dir,
            config,
            work_items,
            file_results,
            &build_result,
            &luts,
            jsonl_stream_ptr,
            actual_threads,
            &processed_count,
            work_items.len,
            if (shouldShowProgress(config)) &progress_root else null,
        );
    }
```

Keep the existing aggregate loop after this block unchanged.

- [ ] **Step 5: Run integration and full Zig tests**

Run:

```bash
zig build test --summary all
```

Expected: all Zig tests pass.

- [ ] **Step 6: Commit adaptive runner**

Run:

```bash
git add /Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig
git commit -m "feat: calibrate adaptive batch workers"
```

---

### Task 5: Update user documentation

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/batch.md`
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/cli/commands.md`

- [ ] **Step 1: Update batch guide**

In `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/batch.md`, add this section after “JSONL for Large Runs” and before “Experimental Adaptive Bitmask SR”:

```markdown
## Experimental Adaptive Workers

For very large single-machine runs, storage can become the bottleneck before CPU. `--adaptive-workers` runs a short calibration pass and then uses the smallest file-worker count that is close to the best observed throughput:

```bash
zsasa batch structures/ results.jsonl \
  --format=jsonl \
  --threads=10 \
  --adaptive-workers
```

`--threads` remains the maximum worker count. This mode is intended for I/O-bound large directory runs on laptops and workstations; it does not distribute work across machines and does not change the output schema. If the input set is too small to calibrate meaningfully, `zsasa` falls back to the requested worker count.
```

- [ ] **Step 2: Update CLI options table**

In `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/cli/commands.md`, add this row near `--threads=N` in the batch/general options table:

```markdown
| `--adaptive-workers` | Experimental batch mode: calibrate file-worker count for I/O-bound large directory runs | off |
```

- [ ] **Step 3: Run docs diff check**

Run:

```bash
git diff --check -- /Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/batch.md /Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/cli/commands.md
```

Expected: no whitespace errors.

- [ ] **Step 4: Commit documentation**

Run:

```bash
git add /Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/batch.md /Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/cli/commands.md
git commit -m "docs: document adaptive batch workers"
```

---

### Task 6: Final verification and performance smoke

**Files:**
- Verify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`
- Verify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/batch.md`
- Verify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/cli/commands.md`

- [ ] **Step 1: Format Zig files**

Run:

```bash
zig fmt src/batch.zig
```

Expected: command exits 0.

- [ ] **Step 2: Run focused Zig checks**

Run:

```bash
zig fmt --check src/
zig build test
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa batch --help | rg "adaptive-workers"
```

Expected: all commands exit 0 and help output includes `--adaptive-workers`.

- [ ] **Step 3: Run small functional CLI check**

Run:

```bash
rm -rf /tmp/zsasa-adaptive-check
mkdir -p /tmp/zsasa-adaptive-check/input
cp examples/1ubq.pdb /tmp/zsasa-adaptive-check/input/1ubq-a.pdb
cp examples/1crn.pdb /tmp/zsasa-adaptive-check/input/1crn-a.pdb
./zig-out/bin/zsasa batch /tmp/zsasa-adaptive-check/input /tmp/zsasa-adaptive-check/results.jsonl \
  --format=jsonl \
  --threads=4 \
  --adaptive-workers \
  --classifier=naccess \
  --quiet \
  --timing
python3 - <<'PY'
from pathlib import Path
import json
path = Path('/tmp/zsasa-adaptive-check/results.jsonl')
rows = [json.loads(line) for line in path.read_text().splitlines() if line]
assert len(rows) == 2, len(rows)
assert all('filename' in row for row in rows)
assert all('total_area' in row for row in rows)
print('jsonl rows', len(rows))
PY
```

Expected: CLI exits 0 and Python prints `jsonl rows 2`.

- [ ] **Step 4: Optional symlink sample performance smoke**

Only run this if time permits and local data is available:

```bash
rm -rf /tmp/zsasa-swissprot-1k
mkdir -p /tmp/zsasa-swissprot-1k
find /Users/nagaet/pdb/afdb/swissprot_pdb_v6 -maxdepth 1 -type f -name '*.pdb' | head -1000 | while read -r f; do ln -s "$f" "/tmp/zsasa-swissprot-1k/$(basename "$f")"; done
./zig-out/bin/zsasa batch /tmp/zsasa-swissprot-1k /tmp/zsasa-swissprot-baseline.jsonl --format=jsonl --threads=10 --use-bitmask --precision=f32 --n-points=128 --classifier=naccess --quiet --timing
./zig-out/bin/zsasa batch /tmp/zsasa-swissprot-1k /tmp/zsasa-swissprot-adaptive.jsonl --format=jsonl --threads=10 --adaptive-workers --use-bitmask --precision=f32 --n-points=128 --classifier=naccess --quiet --timing
```

Expected: both commands exit 0. Record `BATCH_TOTAL_TIME_MS`, `BATCH_SUCCESS`, and whether adaptive selected a worker count when not run with `--quiet`.

- [ ] **Step 5: Final commit if verification caused formatting changes**

If `zig fmt` changed files, run:

```bash
git add /Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig
git commit -m "style: format adaptive batch workers"
```

If there are no changes, do not create an empty commit.

# Batch Chunked I/O Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add experimental chunk-range scheduling and chunked JSONL writing so SwissProt-scale batch runs can compare current behavior, simple chunking, and batched writer I/O.

**Architecture:** Keep the current per-file parallel worker as the default path. Add chunked worker context/functions that claim contiguous ranges with an atomic `next_index`, optionally buffer JSONL lines per chunk, and write them through the existing JSONL stream mutex once per chunk.

**Tech Stack:** Zig 0.16, existing `src/batch.zig` batch pipeline, existing JSONL writer helpers, Zig unit tests, CLI smoke tests.

---

### Task 1: CLI/config plumbing for experimental chunk options

**Files:**
- Modify: `src/batch.zig`
- Test: `src/batch.zig`

- [ ] **Step 1: Add failing parser/config tests**

Add tests near existing `BatchArgs` tests:

```zig
test "BatchArgs --chunk-size" {
    const args = [_][]const u8{ "zsasa", "batch", "--chunk-size=256", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(usize, 256), parsed.chunk_size.?);
    try std.testing.expectEqual(true, parsed.chunk_size_explicit);
}

test "BatchArgs --chunked-jsonl" {
    const args = [_][]const u8{ "zsasa", "batch", "--chunked-jsonl", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.chunked_jsonl);
    try std.testing.expectEqual(true, parsed.chunked_jsonl_explicit);
}
```

- [ ] **Step 2: Verify RED**

Run: `zig build test`

Expected: compile failure because `chunk_size`, `chunk_size_explicit`, `chunked_jsonl`, and `chunked_jsonl_explicit` do not exist.

- [ ] **Step 3: Add fields and parser cases**

Add to `BatchConfig`:

```zig
chunk_size: ?usize = null,
chunked_jsonl: bool = false,
```

Add to `BatchArgs`:

```zig
chunk_size: ?usize = null,
chunked_jsonl: bool = false,
chunk_size_explicit: bool = false,
chunked_jsonl_explicit: bool = false,
```

Add parse support:

```zig
} else if (std.mem.startsWith(u8, arg, "--chunk-size=")) {
    result.chunk_size_explicit = true;
    result.chunk_size = std.fmt.parseInt(usize, arg["--chunk-size=".len..], 10) catch blk: {
        std.debug.print("Error: Invalid --chunk-size value: {s}\n", .{arg["--chunk-size=".len..]});
        break :blk null;
    };
} else if (std.mem.eql(u8, arg, "--chunked-jsonl")) {
    result.chunked_jsonl = true;
    result.chunked_jsonl_explicit = true;
```

Copy parsed values into `BatchConfig` in `run()`.

- [ ] **Step 4: Verify GREEN**

Run: `zig build test`

Expected: tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/batch.zig
git commit -m "feat: add experimental batch chunk flags"
```

### Task 2: Validate experimental option combinations

**Files:**
- Modify: `src/batch.zig`
- Test: `src/batch.zig`

- [ ] **Step 1: Add failing validation tests**

Add tests for a helper named `validateChunkOptions`:

```zig
test "validateChunkOptions rejects zero chunk size" {
    try std.testing.expectError(error.InvalidArgument, validateChunkOptions(.{ .chunk_size = 0 }));
}

test "validateChunkOptions rejects chunked jsonl without chunk size" {
    try std.testing.expectError(error.InvalidArgument, validateChunkOptions(.{ .chunked_jsonl = true, .output_format = .jsonl }));
}

test "validateChunkOptions rejects chunked jsonl for non-jsonl output" {
    try std.testing.expectError(error.InvalidArgument, validateChunkOptions(.{ .chunk_size = 256, .chunked_jsonl = true, .output_format = .json }));
}

test "validateChunkOptions accepts chunked jsonl with jsonl and chunk size" {
    try validateChunkOptions(.{ .chunk_size = 256, .chunked_jsonl = true, .output_format = .jsonl });
}
```

- [ ] **Step 2: Verify RED**

Run: `zig build test`

Expected: compile failure because `validateChunkOptions` does not exist.

- [ ] **Step 3: Add helper and call it**

Add:

```zig
fn validateChunkOptions(config: BatchConfig) !void {
    if (config.chunk_size) |size| {
        if (size == 0) return error.InvalidArgument;
    }
    if (config.chunked_jsonl) {
        if (config.chunk_size == null) return error.InvalidArgument;
        if (config.output_format != .jsonl) return error.InvalidArgument;
    }
}
```

Call it from `runBatchSequential` and `runBatchParallel` after `validateBatchOutputFormat(config.output_format)`.

In CLI `run()`, print useful errors before running batch:

```zig
if (config.chunk_size != null and config.chunk_size.? == 0) {
    std.debug.print("Error: --chunk-size must be greater than zero\n", .{});
    return error.InvalidArgument;
}
if (config.chunked_jsonl and config.chunk_size == null) {
    std.debug.print("Error: --chunked-jsonl requires --chunk-size=N\n", .{});
    return error.InvalidArgument;
}
if (config.chunked_jsonl and config.output_format != .jsonl) {
    std.debug.print("Error: --chunked-jsonl requires --format=jsonl\n", .{});
    return error.InvalidArgument;
}
```

- [ ] **Step 4: Verify GREEN**

Run: `zig build test`

Expected: tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/batch.zig
git commit -m "feat: validate batch chunk options"
```

### Task 3: Add chunked range scheduler with existing per-file JSONL writes

**Files:**
- Modify: `src/batch.zig`
- Test: `src/batch.zig`

- [ ] **Step 1: Add failing JSONL parseability test for simple chunk mode**

Add a test modeled on existing parallel JSONL tests:

```zig
test "runBatchParallel chunk size writes parseable JSONL" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_path = try std.fs.path.join(allocator, &.{ root_path, "chunk-results.jsonl" });
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
        .chunk_size = 3,
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

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(64 * 1024));
    defer allocator.free(content);
    try std.testing.expectEqual(@as(usize, 10), std.mem.count(u8, content, "\n"));
}
```

- [ ] **Step 2: Verify RED**

Run: `zig build test`

Expected: test fails because `chunk_size` is ignored or not yet supported by a chunk scheduler observation. If it passes because output behavior is unchanged, add an internal helper test for `claimChunkRange` that fails until implemented.

- [ ] **Step 3: Implement chunk range scheduler**

Add a `ChunkParallelContext` similar to `ParallelContext` but with `chunk_size`, `next_index`, and `chunked_jsonl` fields. Add `claimChunkRange`:

```zig
const ChunkRange = struct { start: usize, end: usize };

fn claimChunkRange(next_index: *std.atomic.Value(usize), total: usize, chunk_size: usize) ?ChunkRange {
    const start = next_index.fetchAdd(chunk_size, .monotonic);
    if (start >= total) return null;
    return .{ .start = start, .end = @min(start + chunk_size, total) };
}
```

Implement `chunkParallelWorker` that loops over claimed ranges and processes each item using the same logic as `parallelWorker`.

In `runBatchParallel`, use the chunked scheduler when `config.chunk_size != null`; otherwise keep the existing path.

- [ ] **Step 4: Verify GREEN**

Run: `zig build test`

Expected: tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/batch.zig
git commit -m "feat: add chunked batch range scheduler"
```

### Task 4: Add chunked JSONL writer path

**Files:**
- Modify: `src/batch.zig`
- Test: `src/batch.zig`

- [ ] **Step 1: Add failing parseability test for chunked writer mode**

Add a test like Task 3 but with `.chunked_jsonl = true`, `.chunk_size = 4`, and 12 input files. It should parse each JSONL line with `std.json.parseFromSlice` and assert `filename`, `total_area`, and `atom_areas` exist.

- [ ] **Step 2: Verify RED**

Run: `zig build test`

Expected: failure because `chunked_jsonl` is not implemented or validation prevents it from working.

- [ ] **Step 3: Implement chunked JSONL buffering**

Extend `JsonlStreamWriter` with:

```zig
pub fn writeBytes(self: *JsonlStreamWriter, bytes: []const u8) void {
    if (bytes.len == 0) return;
    self.mutex.lockUncancelable(self.io);
    defer self.mutex.unlock(self.io);

    var buf: [64 * 1024]u8 = undefined;
    var w = std.Io.File.Writer.initStreaming(self.file, self.io, &buf);
    w.interface.writeAll(bytes) catch |err| {
        logWarning("JSONL chunk write failed: {s}", .{@errorName(err)});
        self.write_failed = true;
        return;
    };
    w.interface.flush() catch |err| {
        logWarning("JSONL chunk flush failed: {s}", .{@errorName(err)});
        self.write_failed = true;
    };
}
```

In `chunkParallelWorker`, when `chunked_jsonl` is true, append each successful line and newline to a per-worker `std.ArrayListUnmanaged(u8)` and call `stream.writeBytes(buffer.items)` once after finishing each claimed range. Reset the buffer retaining capacity after each chunk.

- [ ] **Step 4: Verify GREEN**

Run: `zig build test`

Expected: tests pass.

- [ ] **Step 5: Commit**

```bash
git add src/batch.zig
git commit -m "feat: buffer batch jsonl writes by chunk"
```

### Task 5: Help/docs and smoke benchmark harness

**Files:**
- Modify: `src/batch.zig`
- Modify: `website/docs/cli/commands.md`
- Modify: `website/docs/guide/batch.md`

- [ ] **Step 1: Update help text**

Add options to `printHelp`:

```zig
\\    --chunk-size=N      Experimental: process parallel batch work in N-item chunks
\\    --chunked-jsonl     Experimental: buffer JSONL writes by chunk; requires --chunk-size and --format=jsonl
```

- [ ] **Step 2: Update docs concisely**

Document that the options are experimental and intended for large-directory benchmarking. Mention that simple chunking and chunked JSONL can be compared separately.

- [ ] **Step 3: Run verification**

Run:

```bash
zig fmt --check src/
zig build test
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa batch --help 2>&1 | rg 'chunk-size|chunked-jsonl'
```

- [ ] **Step 4: Run small smoke comparisons**

Create 12 small PDB inputs under `/tmp/zsasa-chunk-smoke/input` and run:

```bash
./zig-out/bin/zsasa batch /tmp/zsasa-chunk-smoke/input /tmp/zsasa-chunk-smoke/base.jsonl --format=jsonl --threads=4 --classifier=naccess --quiet --timing
./zig-out/bin/zsasa batch /tmp/zsasa-chunk-smoke/input /tmp/zsasa-chunk-smoke/simple.jsonl --format=jsonl --threads=4 --chunk-size=4 --classifier=naccess --quiet --timing
./zig-out/bin/zsasa batch /tmp/zsasa-chunk-smoke/input /tmp/zsasa-chunk-smoke/chunked.jsonl --format=jsonl --threads=4 --chunk-size=4 --chunked-jsonl --classifier=naccess --quiet --timing
```

Parse all three JSONL files and assert each has 12 rows and the same keys.

- [ ] **Step 5: Commit**

```bash
git add src/batch.zig website/docs/cli/commands.md website/docs/guide/batch.md
git commit -m "docs: document experimental batch chunk options"
```

### Task 6: SwissProt benchmark comparison

**Files:**
- No source changes unless benchmark reveals a correctness issue.

- [ ] **Step 1: Build ReleaseFast**

Run: `zig build -Doptimize=ReleaseFast`

- [ ] **Step 2: Run SwissProt full comparison**

Run the three benchmark variants against `/Users/nagaet/pdb/afdb/swissprot_pdb_v6` with output under `/tmp/zsasa-swissprot-chunk-bench`:

```bash
./zig-out/bin/zsasa batch /Users/nagaet/pdb/afdb/swissprot_pdb_v6 /tmp/zsasa-swissprot-chunk-bench/base.jsonl --format=jsonl --threads=10 --precision=f32 --n-points=100 --use-bitmask --timing --quiet
./zig-out/bin/zsasa batch /Users/nagaet/pdb/afdb/swissprot_pdb_v6 /tmp/zsasa-swissprot-chunk-bench/simple.jsonl --format=jsonl --threads=10 --chunk-size=256 --precision=f32 --n-points=100 --use-bitmask --timing --quiet
./zig-out/bin/zsasa batch /Users/nagaet/pdb/afdb/swissprot_pdb_v6 /tmp/zsasa-swissprot-chunk-bench/chunked.jsonl --format=jsonl --threads=10 --chunk-size=256 --chunked-jsonl --precision=f32 --n-points=100 --use-bitmask --timing --quiet
```

Use `/usr/bin/time -l` around each command and delete the 15 GB JSONL output after collecting each log.

- [ ] **Step 3: Summarize results**

Report wall time, files/s, `BATCH_TOTAL_TIME_MS`, and max RSS. If neither simple chunk nor chunked JSONL improves SwissProt wall time, recommend not merging the feature.

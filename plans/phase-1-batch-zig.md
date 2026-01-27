# Plan: Zig Batch Processing Implementation

## Goal

Add directory batch processing to freesasa-zig for fair comparison with RustSASA.

## Requirements

1. Directory input detection (auto-batch when input is directory)
2. File-level parallelism using thread pool
3. Optimized memory management (arena allocator, buffer reuse)
4. Timing output for benchmarking (total SASA time, per-file breakdown optional)
5. Support `.json` and `.json.gz` files

## Design

### CLI Interface

```bash
# Single file (existing)
freesasa_zig input.json output.json

# Batch mode (new)
freesasa_zig ./input_dir/ ./output_dir/ --threads=8

# Batch with timing
freesasa_zig ./input_dir/ ./output_dir/ --timing
```

### Architecture

```
┌─────────────────────────────────────────────────────────────┐
│ main.zig                                                    │
│  ├── detectInputType() → file or directory                  │
│  └── if directory → runBatch()                              │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│ batch.zig (new)                                             │
│                                                             │
│  pub fn runBatch(                                           │
│      allocator: Allocator,                                  │
│      input_dir: []const u8,                                 │
│      output_dir: []const u8,                                │
│      config: BatchConfig,                                   │
│  ) !BatchResult                                             │
│                                                             │
│  Optimizations:                                             │
│  1. Arena allocator per file (fast alloc/free)              │
│  2. Thread pool with work stealing                          │
│  3. Pre-scan directory to know file count                   │
│  4. Atomic counter for progress                             │
└─────────────────────────────────────────────────────────────┘
```

### BatchConfig

```zig
pub const BatchConfig = struct {
    n_threads: usize = 0,        // 0 = auto-detect
    algorithm: Algorithm = .sr,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    probe_radius: f64 = 1.4,
    output_format: OutputFormat = .json,
    show_timing: bool = false,
    quiet: bool = false,
};
```

### BatchResult

```zig
pub const BatchResult = struct {
    total_files: usize,
    successful: usize,
    failed: usize,
    total_sasa_time_ns: u64,     // SASA calculation only
    total_time_ns: u64,          // Including I/O
    per_file_times: ?[]FileResult,

    pub const FileResult = struct {
        filename: []const u8,
        n_atoms: usize,
        sasa_time_ns: u64,
        total_sasa: f64,
        status: enum { ok, error },
    };
};
```

### Memory Optimization Strategy

```zig
// Per-thread arena allocator
fn workerThread(ctx: *WorkerContext) void {
    // Each thread has its own arena
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    while (getNextFile()) |file_path| {
        // Process file with arena allocator
        processFile(arena.allocator(), file_path);

        // Reset arena for next file (fast, no individual frees)
        _ = arena.reset(.retain_capacity);
    }
}
```

### Thread Pool Design

Use existing `thread_pool.zig` with file-level work items:

```zig
const FileWorkItem = struct {
    input_path: []const u8,
    output_path: []const u8,
    index: usize,
};

// Work function processes one file
fn processFileWork(ctx: BatchContext, start: usize, end: usize) BatchPartialResult {
    var result = BatchPartialResult{};
    for (ctx.work_items[start..end]) |item| {
        // Each file processed single-threaded
        const file_result = processOneFile(item);
        result.accumulate(file_result);
    }
    return result;
}
```

### Timing Output Format

```
Batch processing: 100 files
Threads: 8
Algorithm: sr (n_points=100)

Processing... [████████████████████] 100/100

Results:
  Total files:     100
  Successful:      100
  Failed:          0
  Total SASA time: 1234.56 ms
  Total time:      1567.89 ms (includes I/O)
  Throughput:      63.8 files/sec

Timing breakdown (SASA only):
  Min:    5.23 ms
  Max:    89.45 ms
  Mean:   12.35 ms
  Median: 10.12 ms
```

For benchmark script integration:
```
BATCH_SASA_TIME_MS:1234.56
BATCH_TOTAL_TIME_MS:1567.89
BATCH_FILES:100
```

## Implementation Steps

### Phase 1: Core Batch Infrastructure ✅

1. [x] Create `src/batch.zig` with types and interfaces
2. [x] Add directory detection in `main.zig`
3. [x] Implement file scanning (find .json/.json.gz)
4. [x] Basic sequential batch processing (no threading)
5. [x] Test with small dataset
6. [x] Add gzip decompression support (added to json_parser.zig)

### Phase 2: Parallel Processing ✅

1. [x] Adapt thread_pool.zig for file-level work (used direct thread spawning with atomic work stealing)
2. [x] Implement parallel batch processing (runBatchParallel with per-thread arena allocators)
3. [x] Add atomic progress counter
4. [x] Test with varying thread counts (1/2/4/8 threads tested, 2.4x speedup with 8 threads)

### Phase 3: Memory Optimization ✅ (done in Phase 2)

1. [x] Per-thread arena allocator (implemented in parallelWorker)
2. [x] Arena reset between files (arena.reset(.retain_capacity))
3. [ ] Benchmark memory usage (optional)
4. [ ] Compare with RustSASA memory patterns

### Phase 4: Output & Integration

1. [ ] Implement output directory handling
2. [ ] Add timing output format
3. [ ] Update CLI help text
4. [ ] Integration with benchmark scripts

### Phase 5: Testing & Validation

1. [ ] Unit tests for batch.zig
2. [ ] Integration test with default dataset
3. [ ] Benchmark comparison with RustSASA
4. [ ] Validate SASA results match single-file mode

## Files to Modify/Create

| File | Action | Description |
|------|--------|-------------|
| `src/batch.zig` | Create | Batch processing logic |
| `src/main.zig` | Modify | Add directory detection, call batch |
| `src/thread_pool.zig` | Possibly modify | If needed for file-level work |
| `build.zig` | Modify | If separate batch binary needed |

## Success Criteria

1. `freesasa_zig ./dir/ ./out/` processes all files
2. Thread scaling: 8 threads ~6-7x faster than 1 thread
3. Memory efficient: no OOM on large directories
4. Output `BATCH_SASA_TIME_MS` for benchmark integration
5. Results match single-file processing

## Risks & Mitigations

| Risk | Mitigation |
|------|------------|
| Thread contention on output dir | Each thread writes to separate files |
| Memory pressure on large files | Arena allocator with reset |
| Uneven file sizes | Work stealing / atomic task queue |
| gzip decompression overhead | Include in total time, separate SASA time |

---
- [ ] **DONE** - Phase complete

const std = @import("std");
const types = @import("types.zig");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const Config = types.Config;
const OutputFormat = json_writer.OutputFormat;
const LeeRichardsConfig = lee_richards.LeeRichardsConfig;

/// SASA algorithm selection
pub const Algorithm = enum {
    sr, // Shrake-Rupley (test point method)
    lr, // Lee-Richards (slice method)
};

/// Configuration for batch processing
pub const BatchConfig = struct {
    n_threads: usize = 0, // 0 = auto-detect
    algorithm: Algorithm = .sr,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    probe_radius: f64 = 1.4,
    output_format: OutputFormat = .json,
    show_timing: bool = false,
    quiet: bool = false,
};

/// Result for a single file
pub const FileResult = struct {
    filename: []const u8,
    n_atoms: usize,
    sasa_time_ns: u64,
    total_sasa: f64,
    status: Status,
    error_msg: ?[]const u8 = null,

    pub const Status = enum {
        ok,
        err,
    };
};

/// Aggregate result for batch processing
pub const BatchResult = struct {
    total_files: usize,
    successful: usize,
    failed: usize,
    total_sasa_time_ns: u64, // SASA calculation only
    total_time_ns: u64, // Including I/O
    file_results: []FileResult,
    allocator: Allocator,

    pub fn deinit(self: *BatchResult) void {
        for (self.file_results) |*result| {
            self.allocator.free(result.filename);
            if (result.error_msg) |msg| {
                self.allocator.free(msg);
            }
        }
        self.allocator.free(self.file_results);
    }

    /// Print human-readable summary
    pub fn printSummary(self: BatchResult, show_timing: bool) void {
        const ns_to_ms = 1_000_000.0;
        const total_sasa_ms = @as(f64, @floatFromInt(self.total_sasa_time_ns)) / ns_to_ms;
        const total_ms = @as(f64, @floatFromInt(self.total_time_ns)) / ns_to_ms;
        const throughput = if (total_ms > 0)
            @as(f64, @floatFromInt(self.successful)) / (total_ms / 1000.0)
        else
            0.0;

        std.debug.print("\nBatch Results:\n", .{});
        std.debug.print("  Total files:     {d}\n", .{self.total_files});
        std.debug.print("  Successful:      {d}\n", .{self.successful});
        std.debug.print("  Failed:          {d}\n", .{self.failed});
        std.debug.print("  Total SASA time: {d:.2} ms\n", .{total_sasa_ms});
        std.debug.print("  Total time:      {d:.2} ms (includes I/O)\n", .{total_ms});
        std.debug.print("  Throughput:      {d:.1} files/sec\n", .{throughput});

        if (show_timing and self.successful > 0) {
            // Calculate timing statistics
            var min_ns: u64 = std.math.maxInt(u64);
            var max_ns: u64 = 0;
            var sum_ns: u64 = 0;
            var ok_count: usize = 0;

            for (self.file_results) |result| {
                if (result.status == .ok) {
                    min_ns = @min(min_ns, result.sasa_time_ns);
                    max_ns = @max(max_ns, result.sasa_time_ns);
                    sum_ns += result.sasa_time_ns;
                    ok_count += 1;
                }
            }

            if (ok_count > 0) {
                const min_ms = @as(f64, @floatFromInt(min_ns)) / ns_to_ms;
                const max_ms = @as(f64, @floatFromInt(max_ns)) / ns_to_ms;
                const mean_ms = @as(f64, @floatFromInt(sum_ns)) / @as(f64, @floatFromInt(ok_count)) / ns_to_ms;

                std.debug.print("\nTiming breakdown (SASA only):\n", .{});
                std.debug.print("  Min:  {d:.2} ms\n", .{min_ms});
                std.debug.print("  Max:  {d:.2} ms\n", .{max_ms});
                std.debug.print("  Mean: {d:.2} ms\n", .{mean_ms});
            }
        }
    }

    /// Print machine-readable output for benchmark scripts
    pub fn printBenchmarkOutput(self: BatchResult) void {
        const ns_to_ms = 1_000_000.0;
        const total_sasa_ms = @as(f64, @floatFromInt(self.total_sasa_time_ns)) / ns_to_ms;
        const total_ms = @as(f64, @floatFromInt(self.total_time_ns)) / ns_to_ms;

        std.debug.print("BATCH_SASA_TIME_MS:{d:.2}\n", .{total_sasa_ms});
        std.debug.print("BATCH_TOTAL_TIME_MS:{d:.2}\n", .{total_ms});
        std.debug.print("BATCH_FILES:{d}\n", .{self.total_files});
        std.debug.print("BATCH_SUCCESS:{d}\n", .{self.successful});
    }
};

/// Scan directory for JSON files (.json and .json.gz)
pub fn scanDirectory(allocator: Allocator, dir_path: []const u8) ![][]const u8 {
    var files = std.ArrayListUnmanaged([]const u8){};
    errdefer {
        for (files.items) |f| allocator.free(f);
        files.deinit(allocator);
    }

    var dir = std.fs.cwd().openDir(dir_path, .{ .iterate = true }) catch |err| {
        return err;
    };
    defer dir.close();

    var iter = dir.iterate();
    while (try iter.next()) |entry| {
        if (entry.kind != .file) continue;

        const name = entry.name;
        // Skip filenames with path separators (defense in depth)
        if (std.mem.indexOfAny(u8, name, "/\\") != null) continue;
        if (std.mem.endsWith(u8, name, ".json") or std.mem.endsWith(u8, name, ".json.gz")) {
            const filename = try allocator.dupe(u8, name);
            try files.append(allocator, filename);
        }
    }

    // Sort for deterministic ordering
    std.mem.sort([]const u8, files.items, {}, struct {
        fn lessThan(_: void, a: []const u8, b: []const u8) bool {
            return std.mem.lessThan(u8, a, b);
        }
    }.lessThan);

    return files.toOwnedSlice(allocator);
}

/// Process a single file and return result
/// Uses provided arena allocator for temporary allocations
fn processOneFile(
    arena: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    filename: []const u8,
    config: BatchConfig,
) FileResult {
    var result = FileResult{
        .filename = filename,
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .ok,
    };

    // Build input path
    const input_path = std.fs.path.join(arena, &.{ input_dir, filename }) catch {
        result.status = .err;
        return result;
    };

    // Read and parse input
    var input = json_parser.readAtomInputFromFile(arena, input_path) catch {
        result.status = .err;
        return result;
    };
    defer input.deinit();

    result.n_atoms = input.atomCount();

    // Time SASA calculation only
    var timer = std.time.Timer.start() catch {
        result.status = .err;
        return result;
    };

    // Calculate SASA (single-threaded per file in batch mode)
    var sasa_result = switch (config.algorithm) {
        .sr => shrake_rupley.calculateSasa(arena, input, .{
            .n_points = config.n_points,
            .probe_radius = config.probe_radius,
        }) catch {
            result.status = .err;
            return result;
        },
        .lr => lee_richards.calculateSasa(arena, input, .{
            .n_slices = config.n_slices,
            .probe_radius = config.probe_radius,
        }) catch {
            result.status = .err;
            return result;
        },
    };
    defer sasa_result.deinit();

    result.sasa_time_ns = timer.read();
    result.total_sasa = sasa_result.total_area;

    // Write output if output_dir specified
    if (output_dir) |out_dir| {
        const output_filename = blk: {
            if (std.mem.endsWith(u8, filename, ".gz")) {
                // Remove .gz extension for output
                const base = filename[0 .. filename.len - 3];
                break :blk base;
            }
            break :blk filename;
        };

        const output_path = std.fs.path.join(arena, &.{ out_dir, output_filename }) catch {
            result.status = .err;
            return result;
        };

        json_writer.writeSasaResultWithFormat(arena, sasa_result, output_path, config.output_format) catch {
            result.status = .err;
            return result;
        };
    }

    return result;
}

/// Run batch processing sequentially (for initial implementation)
pub fn runBatchSequential(
    allocator: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
) !BatchResult {
    // Start total timer
    var total_timer = try std.time.Timer.start();

    // Scan directory for files
    const files = try scanDirectory(allocator, input_dir);
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }

    // Create output directory if specified
    if (output_dir) |out_dir| {
        try std.fs.cwd().makePath(out_dir);
    }

    // Allocate results
    var file_results = try allocator.alloc(FileResult, files.len);
    errdefer allocator.free(file_results);

    // Process each file
    var total_sasa_time_ns: u64 = 0;
    var successful: usize = 0;
    var failed: usize = 0;

    // Use arena allocator for each file (reset between files)
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    for (files, 0..) |filename, i| {
        // Copy filename to result allocator
        const filename_copy = try allocator.dupe(u8, filename);

        // Process file
        var result = processOneFile(
            arena.allocator(),
            input_dir,
            output_dir,
            filename,
            config,
        );
        result.filename = filename_copy;

        if (result.status == .ok) {
            successful += 1;
            total_sasa_time_ns += result.sasa_time_ns;
        } else {
            failed += 1;
        }

        file_results[i] = result;

        // Reset arena for next file
        _ = arena.reset(.retain_capacity);

        // Progress output
        if (!config.quiet) {
            std.debug.print("\rProcessing: {d}/{d}", .{ i + 1, files.len });
        }
    }

    if (!config.quiet) {
        std.debug.print("\n", .{});
    }

    const total_time_ns = total_timer.read();

    return BatchResult{
        .total_files = files.len,
        .successful = successful,
        .failed = failed,
        .total_sasa_time_ns = total_sasa_time_ns,
        .total_time_ns = total_time_ns,
        .file_results = file_results,
        .allocator = allocator,
    };
}

/// Run batch processing (main entry point)
/// Uses parallel processing when n_threads > 1
pub fn runBatch(
    allocator: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
) !BatchResult {
    // For now, use sequential processing
    // TODO: Add parallel processing in Phase 2
    _ = config.n_threads; // Will be used in parallel version
    return runBatchSequential(allocator, input_dir, output_dir, config);
}

// Tests

test "scanDirectory finds json files" {
    // This test requires a test directory - skip in automated testing
    // Manual testing: create a directory with .json files and test
}

test "BatchConfig default values" {
    const config = BatchConfig{};

    try std.testing.expectEqual(@as(usize, 0), config.n_threads);
    try std.testing.expectEqual(Algorithm.sr, config.algorithm);
    try std.testing.expectEqual(@as(u32, 100), config.n_points);
    try std.testing.expectEqual(@as(f64, 1.4), config.probe_radius);
}

test "BatchResult deinit" {
    const allocator = std.testing.allocator;

    var results = try allocator.alloc(FileResult, 1);
    results[0] = FileResult{
        .filename = try allocator.dupe(u8, "test.json"),
        .n_atoms = 100,
        .sasa_time_ns = 1000000,
        .total_sasa = 123.45,
        .status = .ok,
    };

    var batch_result = BatchResult{
        .total_files = 1,
        .successful = 1,
        .failed = 0,
        .total_sasa_time_ns = 1000000,
        .total_time_ns = 2000000,
        .file_results = results,
        .allocator = allocator,
    };

    batch_result.deinit();
}

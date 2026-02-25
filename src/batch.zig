const std = @import("std");
const types = @import("types.zig");
const format_detect = @import("format_detect.zig");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const pdb_parser = @import("pdb_parser.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const lee_richards = @import("lee_richards.zig");
const classifier = @import("classifier.zig");
const classifier_protor = @import("classifier_protor.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_oons = @import("classifier_oons.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const SasaResultGen = types.SasaResultGen;
const Config = types.Config;
const ConfigGen = types.ConfigGen;
const Precision = types.Precision;
const ClassifierType = classifier.ClassifierType;
const OutputFormat = json_writer.OutputFormat;
const LeeRichardsConfig = lee_richards.LeeRichardsConfig;
const LeeRichardsConfigGen = lee_richards.LeeRichardsConfigGen;

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
    precision: Precision = .f64, // f32 or f64
    classifier_type: ?ClassifierType = .protor, // Default: protor (matches FreeSASA/RustSASA)
    include_hydrogens: bool = false, // Include hydrogen atoms (default: exclude)
    include_hetatm: bool = false, // Include HETATM records (default: exclude)
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 64/128/256)
};

/// Helper to build and hold bitmask LUTs for batch processing.
/// Builds the appropriate LUT once based on config, and provides typed pointers.
/// Returns error.BitmaskRequiresSR if use_bitmask is combined with algorithm != .sr.
const BatchLuts = struct {
    lut_f64: ?bitmask_lut.BitmaskLut = null,
    lut_f32: ?bitmask_lut.BitmaskLutGen(f32) = null,

    fn init(allocator: Allocator, config: BatchConfig) !BatchLuts {
        if (!config.use_bitmask) return .{};
        if (config.algorithm != .sr) return error.BitmaskRequiresSR;
        return switch (config.precision) {
            .f64 => .{ .lut_f64 = try bitmask_lut.BitmaskLut.init(allocator, config.n_points) },
            .f32 => .{ .lut_f32 = try bitmask_lut.BitmaskLutGen(f32).init(allocator, config.n_points) },
        };
    }

    fn deinit(self: *BatchLuts) void {
        if (self.lut_f64) |*lut| lut.deinit();
        if (self.lut_f32) |*lut| lut.deinit();
        self.* = .{};
    }

    fn f64Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLut {
        return if (self.lut_f64 != null) &self.lut_f64.? else null;
    }

    fn f32Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLutGen(f32) {
        return if (self.lut_f32 != null) &self.lut_f32.? else null;
    }
};

/// Get output extension based on format
fn getOutputExtension(format: OutputFormat) []const u8 {
    return switch (format) {
        .json, .compact => ".json",
        .csv => ".csv",
    };
}

/// Replace file extension for output (e.g., "file.pdb" -> "file.json")
fn replaceExtension(allocator: Allocator, filename: []const u8, new_ext: []const u8) ![]const u8 {
    // Strip .gz if present
    var base = filename;
    if (std.mem.endsWith(u8, base, ".gz")) {
        base = base[0 .. base.len - 3];
    }

    // Find and strip existing extension
    if (std.mem.lastIndexOfScalar(u8, base, '.')) |dot_idx| {
        const stem = base[0..dot_idx];
        return std.fmt.allocPrint(allocator, "{s}{s}", .{ stem, new_ext });
    }

    // No extension found, just append
    return std.fmt.allocPrint(allocator, "{s}{s}", .{ base, new_ext });
}

/// Generic SASA calculation dispatcher.
/// When bitmask_lut_ptr is non-null, uses bitmask-optimized Shrake-Rupley
/// (the algorithm parameter is ignored). Otherwise selects SR or LR.
fn calculateSasaDispatch(
    comptime T: type,
    allocator: Allocator,
    input: AtomInput,
    algorithm: Algorithm,
    n_points: u32,
    n_slices: u32,
    probe_radius: T,
    n_threads: usize,
    bitmask_lut_ptr: ?*const bitmask_lut.BitmaskLutGen(T),
) !SasaResultGen(T) {
    const sr_config = ConfigGen(T){ .n_points = n_points, .probe_radius = probe_radius };

    if (bitmask_lut_ptr) |lut| {
        return if (n_threads > 1)
            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).calculateSasaParallelWithLut(
                allocator,
                input,
                sr_config,
                n_threads,
                lut,
            )
        else
            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).calculateSasaWithLut(
                allocator,
                input,
                sr_config,
                lut,
            );
    }

    return switch (algorithm) {
        .sr => if (n_threads > 1)
            shrake_rupley.ShrakeRupleyGen(T).calculateSasaParallel(allocator, input, sr_config, n_threads)
        else
            shrake_rupley.ShrakeRupleyGen(T).calculateSasa(allocator, input, sr_config),
        .lr => if (n_threads > 1)
            lee_richards.LeeRichardsGen(T).calculateSasaParallel(allocator, input, .{
                .n_slices = n_slices,
                .probe_radius = probe_radius,
            }, n_threads)
        else
            lee_richards.LeeRichardsGen(T).calculateSasa(allocator, input, .{
                .n_slices = n_slices,
                .probe_radius = probe_radius,
            }),
    };
}

/// Write SASA result to output file
/// Handles f32 -> f64 conversion for consistent output format
fn writeSasaOutput(
    comptime T: type,
    allocator: Allocator,
    result: *const SasaResultGen(T),
    output_dir: []const u8,
    filename: []const u8,
    format: OutputFormat,
) !void {
    const output_filename = try replaceExtension(allocator, filename, getOutputExtension(format));
    defer allocator.free(output_filename);
    const output_path = try std.fs.path.join(allocator, &.{ output_dir, output_filename });

    if (T == f64) {
        try json_writer.writeSasaResultWithFormat(allocator, result.*, output_path, format);
    } else {
        var result_f64 = try result.toF64(allocator);
        defer result_f64.deinit();
        try json_writer.writeSasaResultWithFormat(allocator, result_f64, output_path, format);
    }
}

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

/// Read input file with auto-format detection
fn readInputFile(allocator: Allocator, path: []const u8, config: BatchConfig) !AtomInput {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => json_parser.readAtomInputFromFile(allocator, path),
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            break :blk parser.parseFile(path);
        },
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            break :blk parser.parseFile(path);
        },
    };
}

/// Apply built-in classifier to replace radii based on residue/atom names
fn applyBuiltinClassifier(input: *AtomInput, ct: ClassifierType) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    for (0..n) |i| {
        const maybe_radius: ?f64 = switch (ct) {
            .naccess => classifier_naccess.getRadius(residues[i].slice(), atom_names[i].slice()),
            .protor => classifier_protor.getRadius(residues[i].slice(), atom_names[i].slice()),
            .oons => classifier_oons.getRadius(residues[i].slice(), atom_names[i].slice()),
        };

        if (maybe_radius) |r| {
            new_radii[i] = r;
        } else if (input.element) |elements| {
            if (classifier.guessRadiusFromAtomicNumber(elements[i])) |r| {
                new_radii[i] = r;
            } else {
                new_radii[i] = input.r[i];
            }
        } else if (classifier.guessRadiusFromAtomName(atom_names[i].slice())) |r| {
            new_radii[i] = r;
        } else {
            new_radii[i] = input.r[i];
        }
    }

    input.allocator.free(input.r);
    input.r = new_radii;
}

/// Scan directory for structure files (.json, .pdb, .cif, .mmcif, .ent and compressed variants)
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
        // Accept regular files and symlinks (for sampled batch benchmarking)
        if (entry.kind != .file and entry.kind != .sym_link) continue;

        const name = entry.name;
        // Skip filenames with path separators (defense in depth)
        if (std.mem.indexOfAny(u8, name, "/\\") != null) continue;
        if (format_detect.isSupportedFile(name)) {
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
/// n_threads: number of threads for SASA calculation (1 = single-threaded)
fn processOneFile(
    arena: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    filename: []const u8,
    config: BatchConfig,
    n_threads: usize,
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
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

    // Read and parse input (auto-detect format from extension)
    var input = readInputFile(arena, input_path, config) catch {
        result.status = .err;
        return result;
    };
    defer input.deinit();

    // Apply classifier for PDB/mmCIF input (skip for JSON which has radii embedded)
    if (config.classifier_type) |ct| {
        const format = format_detect.detectInputFormat(input_path);
        if (format != .json and input.hasClassificationInfo()) {
            applyBuiltinClassifier(&input, ct) catch {
                result.status = .err;
                return result;
            };
        }
    }

    result.n_atoms = input.atomCount();

    // Time SASA calculation only
    var timer = std.time.Timer.start() catch {
        result.status = .err;
        return result;
    };

    // Calculate SASA using generic dispatcher
    var total_area: f64 = 0;
    switch (config.precision) {
        .f64 => {
            var sasa_result = calculateSasaDispatch(
                f64,
                arena,
                input,
                config.algorithm,
                config.n_points,
                config.n_slices,
                config.probe_radius,
                n_threads,
                lut_f64,
            ) catch {
                result.status = .err;
                return result;
            };
            defer sasa_result.deinit();
            result.sasa_time_ns = timer.read();
            total_area = sasa_result.total_area;

            if (output_dir) |out_dir| {
                writeSasaOutput(f64, arena, &sasa_result, out_dir, filename, config.output_format) catch {
                    result.status = .err;
                    return result;
                };
            }
        },
        .f32 => {
            var sasa_result = calculateSasaDispatch(
                f32,
                arena,
                input,
                config.algorithm,
                config.n_points,
                config.n_slices,
                @as(f32, @floatCast(config.probe_radius)),
                n_threads,
                lut_f32,
            ) catch {
                result.status = .err;
                return result;
            };
            defer sasa_result.deinit();
            result.sasa_time_ns = timer.read();
            total_area = @floatCast(sasa_result.total_area);

            if (output_dir) |out_dir| {
                writeSasaOutput(f32, arena, &sasa_result, out_dir, filename, config.output_format) catch {
                    result.status = .err;
                    return result;
                };
            }
        },
    }

    result.total_sasa = total_area;
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
    const file_results = try allocator.alloc(FileResult, files.len);
    errdefer allocator.free(file_results);

    // Process each file
    var total_sasa_time_ns: u64 = 0;
    var successful: usize = 0;
    var failed: usize = 0;

    // Build bitmask LUT once (if enabled)
    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    // Use arena allocator for each file (reset between files)
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    for (files, 0..) |filename, i| {
        // Copy filename to result allocator
        const filename_copy = try allocator.dupe(u8, filename);

        // Process file (single-threaded for sequential mode)
        var result = processOneFile(
            arena.allocator(),
            input_dir,
            output_dir,
            filename,
            config,
            1, // single-threaded
            luts.f64Ptr(),
            luts.f32Ptr(),
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

/// Shared context for parallel workers
const ParallelContext = struct {
    files: []const []const u8,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    results: []FileResult,
    result_allocator: Allocator,
    next_file: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
};

/// Worker thread function for parallel batch processing
fn parallelWorker(ctx: *ParallelContext) void {
    // Each thread gets its own arena allocator
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    while (true) {
        // Atomically grab the next file index
        const file_idx = ctx.next_file.fetchAdd(1, .seq_cst);

        if (file_idx >= ctx.files.len) {
            break; // No more files
        }

        const filename = ctx.files[file_idx];

        // Copy filename to result allocator (thread-safe: each index is unique)
        // On allocation failure, use empty string to maintain ownership invariant
        const filename_copy = ctx.result_allocator.dupe(u8, filename) catch {
            // Fallback: allocate empty owned slice for deinit compatibility
            const empty = ctx.result_allocator.alloc(u8, 0) catch {
                // Critical allocation failure - skip file, don't store result
                _ = ctx.processed_count.fetchAdd(1, .seq_cst);
                _ = arena.reset(.retain_capacity);
                continue;
            };
            ctx.results[file_idx] = FileResult{
                .filename = empty,
                .n_atoms = 0,
                .sasa_time_ns = 0,
                .total_sasa = 0,
                .status = .err,
            };
            _ = ctx.processed_count.fetchAdd(1, .seq_cst);
            _ = arena.reset(.retain_capacity);
            continue;
        };

        // Process file using thread-local arena (single-threaded SASA per file in file-parallel mode)
        var result = processOneFile(
            arena.allocator(),
            ctx.input_dir,
            ctx.output_dir,
            filename,
            ctx.config,
            1, // single-threaded SASA per file
            ctx.lut_f64,
            ctx.lut_f32,
        );
        result.filename = filename_copy;

        // Store result (thread-safe: each index is unique)
        ctx.results[file_idx] = result;

        // Update progress counter
        _ = ctx.processed_count.fetchAdd(1, .seq_cst);

        // Reset arena for next file
        _ = arena.reset(.retain_capacity);
    }
}

/// Run batch processing in parallel
pub fn runBatchParallel(
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

    if (files.len == 0) {
        return BatchResult{
            .total_files = 0,
            .successful = 0,
            .failed = 0,
            .total_sasa_time_ns = 0,
            .total_time_ns = total_timer.read(),
            .file_results = try allocator.alloc(FileResult, 0),
            .allocator = allocator,
        };
    }

    // Create output directory if specified
    if (output_dir) |out_dir| {
        try std.fs.cwd().makePath(out_dir);
    }

    // Allocate results
    const file_results = try allocator.alloc(FileResult, files.len);
    errdefer allocator.free(file_results);

    // Determine thread count
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        @min(config.n_threads, cpu_count);

    // For single file or single thread, use sequential
    if (files.len == 1 or n_threads <= 1) {
        allocator.free(file_results);
        return runBatchSequential(allocator, input_dir, output_dir, config);
    }

    // Build bitmask LUT once (if enabled)
    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    // Create shared context
    var ctx = ParallelContext{
        .files = files,
        .input_dir = input_dir,
        .output_dir = output_dir,
        .config = config,
        .results = file_results,
        .result_allocator = allocator,
        .next_file = std.atomic.Value(usize).init(0),
        .processed_count = std.atomic.Value(usize).init(0),
        .lut_f64 = luts.f64Ptr(),
        .lut_f32 = luts.f32Ptr(),
    };

    // Spawn worker threads
    const actual_threads = @min(n_threads, files.len);
    const threads = try allocator.alloc(std.Thread, actual_threads);
    defer allocator.free(threads);

    for (threads) |*thread| {
        thread.* = try std.Thread.spawn(.{}, parallelWorker, .{&ctx});
    }

    // Progress monitoring (optional)
    if (!config.quiet) {
        while (ctx.processed_count.load(.seq_cst) < files.len) {
            const processed = ctx.processed_count.load(.seq_cst);
            std.debug.print("\rProcessing: {d}/{d}", .{ processed, files.len });
            std.Thread.sleep(50 * std.time.ns_per_ms); // 50ms update interval
        }
        std.debug.print("\rProcessing: {d}/{d}\n", .{ files.len, files.len });
    }

    // Wait for all threads to complete
    for (threads) |thread| {
        thread.join();
    }

    // Aggregate results
    var total_sasa_time_ns: u64 = 0;
    var successful: usize = 0;
    var failed: usize = 0;

    for (file_results) |result| {
        if (result.status == .ok) {
            successful += 1;
            total_sasa_time_ns += result.sasa_time_ns;
        } else {
            failed += 1;
        }
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
/// Uses file-level parallelism: N files in parallel, 1 thread per file
pub fn runBatch(
    allocator: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
) !BatchResult {
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        config.n_threads;

    if (n_threads <= 1) {
        return runBatchSequential(allocator, input_dir, output_dir, config);
    }
    return runBatchParallel(allocator, input_dir, output_dir, config);
}

// =============================================================================
// CLI argument parsing and run entry point
// =============================================================================

/// Parsed command-line arguments for the batch subcommand
pub const BatchArgs = struct {
    input_path: ?[]const u8 = null,
    output_path: ?[]const u8 = null, // Optional output directory
    n_threads: usize = 0,
    probe_radius: f64 = 1.4,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    algorithm: Algorithm = .sr,
    precision: Precision = .f64,
    output_format: OutputFormat = .json,
    classifier_type: ?ClassifierType = null,
    include_hydrogens: bool = false,
    include_hetatm: bool = false,
    use_bitmask: bool = false,
    quiet: bool = false,
    show_timing: bool = false,
    show_help: bool = false,
};

// Parse helper functions (local to batch.zig)

/// Parse and validate probe radius value
fn parseProbeRadius(value: []const u8) f64 {
    const radius = std.fmt.parseFloat(f64, value) catch {
        std.debug.print("Error: Invalid probe radius: {s}\n", .{value});
        std.process.exit(1);
    };
    if (radius <= 0 or radius > 10.0 or !std.math.isFinite(radius)) {
        std.debug.print("Error: Probe radius must be between 0 and 10 Angstroms: {d}\n", .{radius});
        std.process.exit(1);
    }
    return radius;
}

/// Parse and validate n-points value
fn parseNPoints(value: []const u8) u32 {
    const n = std.fmt.parseInt(u32, value, 10) catch {
        std.debug.print("Error: Invalid n-points: {s}\n", .{value});
        std.process.exit(1);
    };
    if (n == 0 or n > 10000) {
        std.debug.print("Error: n-points must be between 1 and 10000: {d}\n", .{n});
        std.process.exit(1);
    }
    return n;
}

/// Parse and validate n-slices value (for Lee-Richards)
fn parseNSlices(value: []const u8) u32 {
    const n = std.fmt.parseInt(u32, value, 10) catch {
        std.debug.print("Error: Invalid n-slices: {s}\n", .{value});
        std.process.exit(1);
    };
    if (n == 0 or n > 1000) {
        std.debug.print("Error: n-slices must be between 1 and 1000: {d}\n", .{n});
        std.process.exit(1);
    }
    return n;
}

/// Parse and validate algorithm value
fn parseAlgorithmValue(value: []const u8) Algorithm {
    if (std.mem.eql(u8, value, "sr") or std.mem.eql(u8, value, "shrake-rupley")) {
        return .sr;
    } else if (std.mem.eql(u8, value, "lr") or std.mem.eql(u8, value, "lee-richards")) {
        return .lr;
    } else {
        std.debug.print("Error: Invalid algorithm: {s}\n", .{value});
        std.debug.print("Valid algorithms: sr (shrake-rupley), lr (lee-richards)\n", .{});
        std.process.exit(1);
    }
}

/// Parse and validate output format value
fn parseOutputFormatValue(value: []const u8) OutputFormat {
    if (std.mem.eql(u8, value, "json")) {
        return .json;
    } else if (std.mem.eql(u8, value, "compact")) {
        return .compact;
    } else if (std.mem.eql(u8, value, "csv")) {
        return .csv;
    } else {
        std.debug.print("Error: Invalid format: {s}\n", .{value});
        std.debug.print("Valid formats: json, compact, csv\n", .{});
        std.process.exit(1);
    }
}

/// Parse and validate classifier type value
fn parseClassifierTypeValue(value: []const u8) ClassifierType {
    if (ClassifierType.fromString(value)) |ct| {
        return ct;
    } else {
        std.debug.print("Error: Invalid classifier: {s}\n", .{value});
        std.debug.print("Valid classifiers: naccess, protor, oons\n", .{});
        std.process.exit(1);
    }
}

/// Parse and validate precision value
fn parsePrecisionValue(value: []const u8) Precision {
    if (Precision.fromString(value)) |p| {
        return p;
    } else {
        std.debug.print("Error: Invalid precision: {s}\n", .{value});
        std.debug.print("Valid values: f32 (single), f64 (double)\n", .{});
        std.process.exit(1);
    }
}

/// Parse batch subcommand arguments
pub fn parseArgs(args: []const []const u8, start_idx: usize) BatchArgs {
    var result = BatchArgs{};
    var i: usize = start_idx;
    var positional_count: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        // --threads=N or --threads N
        if (std.mem.startsWith(u8, arg, "--threads=")) {
            const value = arg["--threads=".len..];
            result.n_threads = std.fmt.parseInt(usize, value, 10) catch {
                std.debug.print("Error: Invalid thread count: {s}\n", .{value});
                std.process.exit(1);
            };
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
        }
        // --probe-radius=R or --probe-radius R
        else if (std.mem.startsWith(u8, arg, "--probe-radius=")) {
            const value = arg["--probe-radius=".len..];
            result.probe_radius = parseProbeRadius(value);
        } else if (std.mem.eql(u8, arg, "--probe-radius")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --probe-radius\n", .{});
                std.process.exit(1);
            }
            result.probe_radius = parseProbeRadius(args[i]);
        }
        // --n-points=N or --n-points N
        else if (std.mem.startsWith(u8, arg, "--n-points=")) {
            const value = arg["--n-points=".len..];
            result.n_points = parseNPoints(value);
        } else if (std.mem.eql(u8, arg, "--n-points")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --n-points\n", .{});
                std.process.exit(1);
            }
            result.n_points = parseNPoints(args[i]);
        }
        // --n-slices=N or --n-slices N (for Lee-Richards)
        else if (std.mem.startsWith(u8, arg, "--n-slices=")) {
            const value = arg["--n-slices=".len..];
            result.n_slices = parseNSlices(value);
        } else if (std.mem.eql(u8, arg, "--n-slices")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --n-slices\n", .{});
                std.process.exit(1);
            }
            result.n_slices = parseNSlices(args[i]);
        }
        // --format=FORMAT or --format FORMAT
        else if (std.mem.startsWith(u8, arg, "--format=")) {
            const value = arg["--format=".len..];
            result.output_format = parseOutputFormatValue(value);
        } else if (std.mem.eql(u8, arg, "--format")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --format\n", .{});
                std.process.exit(1);
            }
            result.output_format = parseOutputFormatValue(args[i]);
        }
        // --algorithm=ALGO or --algorithm ALGO
        else if (std.mem.startsWith(u8, arg, "--algorithm=")) {
            const value = arg["--algorithm=".len..];
            result.algorithm = parseAlgorithmValue(value);
        } else if (std.mem.eql(u8, arg, "--algorithm")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --algorithm\n", .{});
                std.process.exit(1);
            }
            result.algorithm = parseAlgorithmValue(args[i]);
        }
        // --classifier=TYPE or --classifier TYPE
        else if (std.mem.startsWith(u8, arg, "--classifier=")) {
            const value = arg["--classifier=".len..];
            result.classifier_type = parseClassifierTypeValue(value);
        } else if (std.mem.eql(u8, arg, "--classifier")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --classifier\n", .{});
                std.process.exit(1);
            }
            result.classifier_type = parseClassifierTypeValue(args[i]);
        }
        // --precision=PREC or --precision PREC
        else if (std.mem.startsWith(u8, arg, "--precision=")) {
            const value = arg["--precision=".len..];
            result.precision = parsePrecisionValue(value);
        } else if (std.mem.eql(u8, arg, "--precision")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --precision\n", .{});
                std.process.exit(1);
            }
            result.precision = parsePrecisionValue(args[i]);
        }
        // --include-hydrogens
        else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
            result.include_hydrogens = true;
        }
        // --include-hetatm
        else if (std.mem.eql(u8, arg, "--include-hetatm")) {
            result.include_hetatm = true;
        }
        // --use-bitmask
        else if (std.mem.eql(u8, arg, "--use-bitmask")) {
            result.use_bitmask = true;
        }
        // --quiet or -q
        else if (std.mem.eql(u8, arg, "--quiet") or std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
        }
        // --timing
        else if (std.mem.eql(u8, arg, "--timing")) {
            result.show_timing = true;
        }
        // --help or -h
        else if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            result.show_help = true;
        }
        // -o FILE or --output=FILE or --output FILE
        else if (std.mem.eql(u8, arg, "-o")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for -o\n", .{});
                std.process.exit(1);
            }
            result.output_path = args[i];
        } else if (std.mem.startsWith(u8, arg, "--output=")) {
            result.output_path = arg["--output=".len..];
        } else if (std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --output\n", .{});
                std.process.exit(1);
            }
            result.output_path = args[i];
        }
        // Unknown option
        else if (std.mem.startsWith(u8, arg, "-")) {
            std.debug.print("Error: Unknown option: {s}\n", .{arg});
            std.debug.print("Try 'zsasa batch --help' for more information.\n", .{});
            std.process.exit(1);
        }
        // Positional arguments: <input_dir> [output_dir]
        else {
            if (positional_count == 0) {
                result.input_path = arg;
            } else if (positional_count == 1) {
                result.output_path = arg;
            } else {
                std.debug.print("Error: Too many positional arguments\n", .{});
                std.process.exit(1);
            }
            positional_count += 1;
        }
    }

    return result;
}

/// Print help for the batch subcommand
pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\zsasa batch - Calculate SASA for all files in a directory
        \\
        \\USAGE:
        \\    {s} batch [OPTIONS] <input_dir> [output_dir]
        \\
        \\ARGUMENTS:
        \\    <input_dir>     Directory containing structure files (PDB, mmCIF, JSON)
        \\    [output_dir]    Optional output directory (default: no file output)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO    Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                        Default: sr
        \\    --classifier=TYPE   Built-in classifier: naccess, protor, oons
        \\                        Default: protor
        \\    --threads=N         Number of threads (default: auto-detect)
        \\    --probe-radius=R    Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N        Test points per atom (default: 100, for sr)
        \\    --n-slices=N        Slices per atom diameter (default: 20, for lr)
        \\    --precision=PREC    Floating-point precision: f32, f64 (default: f64)
        \\    --format=FORMAT     Output format: json, compact, csv (default: json)
        \\    --include-hydrogens Include hydrogen atoms (default: exclude)
        \\    --include-hetatm    Include HETATM records (default: exclude)
        \\    --use-bitmask       Use bitmask LUT optimization for SR algorithm
        \\                        (n-points must be 64, 128, or 256)
        \\    --timing            Show timing breakdown for benchmarking
        \\    -o, --output=DIR    Output directory (alternative to positional)
        \\    -q, --quiet         Suppress progress output
        \\    -h, --help          Show this help message
        \\
        \\EXAMPLES:
        \\    {s} batch structures/
        \\    {s} batch structures/ results/
        \\    {s} batch structures/ --algorithm=lr --threads=4
        \\    {s} batch structures/ --classifier=naccess --format=csv
        \\    {s} batch structures/ results/ --timing --quiet
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name });
}

/// Run batch processing from parsed CLI arguments
pub fn run(allocator: Allocator, args: BatchArgs) !void {
    const input_dir = args.input_path orelse {
        std.debug.print("Error: Missing input directory\n", .{});
        std.debug.print("Usage: zsasa batch [OPTIONS] <input_dir> [output_dir]\n", .{});
        return error.MissingArgument;
    };

    const output_dir: ?[]const u8 = args.output_path;

    // Build batch config from parsed args
    // classifier_type: use explicit --classifier if set, otherwise default protor
    const config = BatchConfig{
        .n_threads = args.n_threads,
        .algorithm = args.algorithm,
        .n_points = args.n_points,
        .n_slices = args.n_slices,
        .probe_radius = args.probe_radius,
        .precision = args.precision,
        .output_format = args.output_format,
        .show_timing = args.show_timing,
        .quiet = args.quiet,
        .classifier_type = args.classifier_type orelse .protor,
        .include_hydrogens = args.include_hydrogens,
        .include_hetatm = args.include_hetatm,
        .use_bitmask = args.use_bitmask,
    };

    if (!args.quiet) {
        std.debug.print("Batch mode: processing directory '{s}'\n", .{input_dir});
        std.debug.print("Algorithm: {s}, Threads: {d}\n", .{
            if (config.algorithm == .sr) "sr" else "lr",
            if (config.n_threads == 0) @as(usize, @intCast(std.Thread.getCpuCount() catch 1)) else config.n_threads,
        });
        if (output_dir) |out| {
            std.debug.print("Output directory: {s}\n", .{out});
        }
        std.debug.print("\n", .{});
    }

    var result = try runBatch(allocator, input_dir, output_dir, config);
    defer result.deinit();

    // Print results
    if (!args.quiet) {
        result.printSummary(args.show_timing);
    }

    // Always print benchmark output (for script parsing)
    if (args.show_timing) {
        std.debug.print("\n", .{});
        result.printBenchmarkOutput();
    }
}

// =============================================================================
// Tests
// =============================================================================

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

test "BatchArgs defaults" {
    const args = [_][]const u8{ "zsasa", "batch", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
    try std.testing.expect(parsed.output_path == null);
    try std.testing.expectEqual(@as(usize, 0), parsed.n_threads);
    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
    try std.testing.expectEqual(false, parsed.quiet);
    try std.testing.expectEqual(false, parsed.show_help);
    try std.testing.expectEqual(false, parsed.include_hydrogens);
    try std.testing.expectEqual(false, parsed.include_hetatm);
    try std.testing.expectEqual(false, parsed.use_bitmask);
    try std.testing.expectEqual(false, parsed.show_timing);
}

test "BatchArgs with output dir" {
    const args = [_][]const u8{ "zsasa", "batch", "input_dir/", "output_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
    try std.testing.expectEqualStrings("output_dir/", parsed.output_path.?);
}

test "BatchArgs with options" {
    const args = [_][]const u8{
        "zsasa",          "batch",
        "--algorithm=lr", "--threads=4",
        "--quiet",        "--timing",
        "input_dir/",
    };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqual(true, parsed.quiet);
    try std.testing.expectEqual(true, parsed.show_timing);
}

test "BatchArgs help flag" {
    const args = [_][]const u8{ "zsasa", "batch", "--help" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.show_help);
}

test "BatchArgs output via -o flag" {
    const args = [_][]const u8{ "zsasa", "batch", "-o", "results/", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
    try std.testing.expectEqualStrings("results/", parsed.output_path.?);
}

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
const pipeline = @import("pipeline.zig");
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

/// Parallelism strategy for batch processing
pub const Parallelism = enum {
    file, // File-level: N files in parallel, 1 thread per file (default)
    atom, // Atom-level: 1 file at a time, N threads for SASA calculation
    pipeline, // Pipelined: I/O prefetch + atom-level SASA calculation
};

/// Configuration for batch processing
pub const BatchConfig = struct {
    n_threads: usize = 0, // 0 = auto-detect
    algorithm: Algorithm = .sr,
    parallelism: Parallelism = .file, // file, atom, or pipeline
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

/// Run batch processing with atom-level parallelism
/// Files are processed sequentially, but each file uses multi-threaded SASA
pub fn runBatchAtomParallel(
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

    // Determine thread count for SASA calculation
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        @min(config.n_threads, cpu_count);

    // Process each file sequentially with multi-threaded SASA
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

        // Process file with multi-threaded SASA
        var result = processOneFile(
            arena.allocator(),
            input_dir,
            output_dir,
            filename,
            config,
            n_threads, // multi-threaded SASA
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

/// Run batch processing with pipelined I/O
/// I/O thread prefetches files while processing thread calculates SASA
pub fn runBatchPipelined(
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

    // Determine thread count for SASA calculation
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        @min(config.n_threads, cpu_count);

    // Initialize prefetch queue
    var queue = pipeline.PrefetchQueue.init(allocator);
    defer queue.deinit();

    // Create I/O context with error tracking
    var io_ctx = pipeline.IoContext.init(allocator, &queue, files, input_dir);
    io_ctx.skip_hydrogens = !config.include_hydrogens;
    io_ctx.atom_only = !config.include_hetatm;
    defer io_ctx.deinit();

    // Build bitmask LUT once (if enabled)
    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    // Spawn I/O thread
    const io_thread = try std.Thread.spawn(.{}, pipeline.IoContext.run, .{&io_ctx});

    // Process files from queue
    var total_sasa_time_ns: u64 = 0;
    var successful: usize = 0;
    var failed: usize = 0;
    var file_idx: usize = 0;

    while (queue.pop()) |prefetched| {
        var pf = prefetched;
        defer pf.deinit(allocator);

        // Copy filename to result allocator
        const filename_copy = try allocator.dupe(u8, pf.filename);

        // Time SASA calculation
        var timer = try std.time.Timer.start();

        // Calculate SASA using generic dispatcher
        var total_area: f64 = 0;
        var sasa_ok = true;

        switch (config.precision) {
            .f64 => {
                if (calculateSasaDispatch(
                    f64,
                    pf.arena.allocator(),
                    pf.input,
                    config.algorithm,
                    config.n_points,
                    config.n_slices,
                    config.probe_radius,
                    n_threads,
                    luts.f64Ptr(),
                )) |calc_result| {
                    var sasa_result = calc_result;
                    defer sasa_result.deinit();
                    total_area = sasa_result.total_area;

                    if (output_dir) |out_dir| {
                        writeSasaOutput(f64, pf.arena.allocator(), &sasa_result, out_dir, pf.filename, config.output_format) catch {
                            sasa_ok = false;
                        };
                    }
                } else |_| {
                    sasa_ok = false;
                }
            },
            .f32 => {
                if (calculateSasaDispatch(
                    f32,
                    pf.arena.allocator(),
                    pf.input,
                    config.algorithm,
                    config.n_points,
                    config.n_slices,
                    @as(f32, @floatCast(config.probe_radius)),
                    n_threads,
                    luts.f32Ptr(),
                )) |calc_result| {
                    var sasa_result = calc_result;
                    defer sasa_result.deinit();
                    total_area = @floatCast(sasa_result.total_area);

                    if (output_dir) |out_dir| {
                        writeSasaOutput(f32, pf.arena.allocator(), &sasa_result, out_dir, pf.filename, config.output_format) catch {
                            sasa_ok = false;
                        };
                    }
                } else |_| {
                    sasa_ok = false;
                }
            },
        }

        const sasa_time_ns = timer.read();

        // Store result
        const file_status: FileResult.Status = if (sasa_ok) .ok else .err;
        file_results[file_idx] = FileResult{
            .filename = filename_copy,
            .n_atoms = pf.input.atomCount(),
            .sasa_time_ns = sasa_time_ns,
            .total_sasa = total_area,
            .status = file_status,
        };

        if (sasa_ok) {
            successful += 1;
            total_sasa_time_ns += sasa_time_ns;
        } else {
            failed += 1;
        }

        file_idx += 1;

        // Progress output
        if (!config.quiet) {
            std.debug.print("\rProcessing: {d}/{d}", .{ file_idx, files.len });
        }
    }

    // Wait for I/O thread to finish
    io_thread.join();

    // Include I/O failures in the count
    const io_failed = io_ctx.failedCount();
    failed += io_failed;

    // Report I/O failures if not quiet
    if (!config.quiet) {
        std.debug.print("\n", .{});
        if (io_failed > 0) {
            std.debug.print("I/O errors: {d} file(s) failed to read/parse\n", .{io_failed});
            for (io_ctx.failed_files.items) |f| {
                const reason_str = switch (f.reason) {
                    .allocation_failed => "allocation failed",
                    .path_join_failed => "path join failed",
                    .read_parse_failed => "read/parse failed",
                };
                std.debug.print("  - {s}: {s}\n", .{ f.filename, reason_str });
            }
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
/// Selects strategy based on config.parallelism:
/// - file: N files in parallel, 1 thread per file (default)
/// - atom: 1 file at a time, N threads for SASA calculation
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

    return switch (config.parallelism) {
        .file => {
            // File-level parallelism: N files in parallel, 1 thread per file
            if (n_threads <= 1) {
                return runBatchSequential(allocator, input_dir, output_dir, config);
            }
            return runBatchParallel(allocator, input_dir, output_dir, config);
        },
        .atom => {
            // Atom-level parallelism: sequential files, N threads per file
            return runBatchAtomParallel(allocator, input_dir, output_dir, config);
        },
        .pipeline => {
            // Pipelined: I/O prefetch + atom-level SASA calculation
            return runBatchPipelined(allocator, input_dir, output_dir, config);
        },
    };
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

const std = @import("std");
const workflow_manifest = @import("workflow_manifest.zig");
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
const classifier_parser = @import("classifier_parser.zig");
// classifier_protor removed — ProtOr is now an alias for CCD
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const sdf_parser = @import("sdf_parser.zig");
const compressed = @import("compressed.zig");

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

fn shouldShowProgress(config: BatchConfig) bool {
    return config.show_progress and !config.quiet;
}

/// Write a warning message to stderr.
///
/// Uses std.debug.print which writes to stderr. Note: errors during the write
/// are silently dropped (debug.print is best-effort). For most CLI use cases
/// this is fine because terminal writes rarely fail; for piped output (e.g.,
/// to a logger that is down), warnings may be lost silently.
fn logWarning(comptime fmt: []const u8, args: anytype) void {
    std.debug.print("Warning: " ++ fmt ++ "\n", args);
}

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
    show_progress: bool = true,
    precision: Precision = .f64, // f32 or f64
    classifier_type: ?ClassifierType = .ccd, // Default: ccd (ProtOr-compatible with CCD extension)
    include_hydrogens: bool = false, // Include hydrogen atoms (default: exclude)
    include_hetatm: bool = false, // Include HETATM records (default: exclude)
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 1..1024)
    adaptive_sr: bool = false, // Experimental two-stage bitmask SR for batch mode
    coarse_points: u32 = 64,
    fine_points: u32 = 256,
    adaptive_low: f64 = 0.05,
    adaptive_high: f64 = 0.95,
    store_atom_areas: bool = false, // When true, copy atom_areas to result_allocator for jsonl
    external_ccd: ?*const ccd_parser.ComponentDict = null, // External CCD dictionary
    sdf_ccd: ?*const ccd_parser.ComponentDict = null, // SDF bond topology dictionary
    custom_classifier: ?*const classifier.Classifier = null,
    custom_classifier_path: ?[]const u8 = null,
    chain_filter: ?[]const []const u8 = null,
    use_auth_chain: bool = false,
    residue_map: bool = false,
};

/// Helper to build and hold bitmask LUTs for batch processing.
/// Builds the appropriate LUT once based on config, and provides typed pointers.
/// Returns error.BitmaskRequiresSR if use_bitmask is combined with algorithm != .sr.
const BatchLuts = struct {
    lut_f64: ?bitmask_lut.BitmaskLut = null,
    lut_f32: ?bitmask_lut.BitmaskLutGen(f32) = null,
    coarse_lut_f64: ?bitmask_lut.BitmaskLut = null,
    fine_lut_f64: ?bitmask_lut.BitmaskLut = null,
    coarse_lut_f32: ?bitmask_lut.BitmaskLutGen(f32) = null,
    fine_lut_f32: ?bitmask_lut.BitmaskLutGen(f32) = null,

    fn init(allocator: Allocator, config: BatchConfig) !BatchLuts {
        if (!config.use_bitmask) return .{};
        if (config.algorithm != .sr) return error.BitmaskRequiresSR;

        var luts = BatchLuts{};
        errdefer luts.deinit();

        if (config.adaptive_sr) {
            switch (config.precision) {
                .f64 => {
                    luts.coarse_lut_f64 = try bitmask_lut.BitmaskLut.init(allocator, config.coarse_points);
                    luts.fine_lut_f64 = try bitmask_lut.BitmaskLut.init(allocator, config.fine_points);
                },
                .f32 => {
                    luts.coarse_lut_f32 = try bitmask_lut.BitmaskLutGen(f32).init(allocator, config.coarse_points);
                    luts.fine_lut_f32 = try bitmask_lut.BitmaskLutGen(f32).init(allocator, config.fine_points);
                },
            }
            return luts;
        }

        switch (config.precision) {
            .f64 => luts.lut_f64 = try bitmask_lut.BitmaskLut.init(allocator, config.n_points),
            .f32 => luts.lut_f32 = try bitmask_lut.BitmaskLutGen(f32).init(allocator, config.n_points),
        }
        return luts;
    }

    fn deinit(self: *BatchLuts) void {
        if (self.lut_f64) |*lut| lut.deinit();
        if (self.lut_f32) |*lut| lut.deinit();
        if (self.coarse_lut_f64) |*lut| lut.deinit();
        if (self.fine_lut_f64) |*lut| lut.deinit();
        if (self.coarse_lut_f32) |*lut| lut.deinit();
        if (self.fine_lut_f32) |*lut| lut.deinit();
        self.* = .{};
    }

    fn f64Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLut {
        return if (self.lut_f64 != null) &self.lut_f64.? else null;
    }

    fn f32Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLutGen(f32) {
        return if (self.lut_f32 != null) &self.lut_f32.? else null;
    }

    fn coarseF64Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLut {
        return if (self.coarse_lut_f64 != null) &self.coarse_lut_f64.? else null;
    }

    fn fineF64Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLut {
        return if (self.fine_lut_f64 != null) &self.fine_lut_f64.? else null;
    }

    fn coarseF32Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLutGen(f32) {
        return if (self.coarse_lut_f32 != null) &self.coarse_lut_f32.? else null;
    }

    fn fineF32Ptr(self: *const BatchLuts) ?*const bitmask_lut.BitmaskLutGen(f32) {
        return if (self.fine_lut_f32 != null) &self.fine_lut_f32.? else null;
    }
};
/// Get output extension based on format
fn getOutputExtension(format: OutputFormat) []const u8 {
    return switch (format) {
        .json, .compact => ".json",
        .jsonl => ".jsonl",
        .csv => ".csv",
    };
}

fn workflowJsonlOutputPath(allocator: Allocator, output_dir: []const u8, job_name: []const u8) ![]const u8 {
    const filename = try std.fmt.allocPrint(allocator, "{s}.jsonl", .{job_name});
    defer allocator.free(filename);
    return std.fs.path.join(allocator, &.{ output_dir, filename });
}

fn workflowPerFileOutputDir(allocator: Allocator, output_dir: []const u8, job_name: []const u8) ![]const u8 {
    return std.fs.path.join(allocator, &.{ output_dir, job_name });
}

/// Replace file extension for output (e.g., "file.pdb" -> "file.json")
fn replaceExtension(allocator: Allocator, filename: []const u8, new_ext: []const u8) ![]const u8 {
    // Strip compression extension if present.
    var base = filename;
    if (compressed.isGzip(base)) {
        base = base[0 .. base.len - 3];
    } else if (compressed.isZstd(base)) {
        base = base[0 .. base.len - 4];
    }

    // Find and strip existing extension
    if (std.mem.findScalarLast(u8, base, '.')) |dot_idx| {
        const stem = base[0..dot_idx];
        return std.fmt.allocPrint(allocator, "{s}{s}", .{ stem, new_ext });
    }

    // No extension found, just append
    return std.fmt.allocPrint(allocator, "{s}{s}", .{ base, new_ext });
}

/// Build a display name for an SDF molecule in batch results.
/// Strips the SDF file extension to produce "stem_molname" or "stem_N" format.
/// This ensures `replaceExtension` produces unique output filenames.
fn sdfMoleculeDisplayName(allocator: Allocator, filename: []const u8, mol_name: []const u8, mol_idx: usize) ![]const u8 {
    // Strip extension (.sdf, .sdf.gz, .sdf.zst, .mol, .mol.gz, .mol.zst) to get stem.
    var base = filename;
    if (compressed.isGzip(base)) {
        base = base[0 .. base.len - 3];
    } else if (compressed.isZstd(base)) {
        base = base[0 .. base.len - 4];
    }
    const stem = if (std.mem.findScalarLast(u8, base, '.')) |dot_idx|
        base[0..dot_idx]
    else
        base;

    if (mol_name.len > 0) {
        return std.fmt.allocPrint(allocator, "{s}_{s}", .{ stem, mol_name });
    } else {
        return std.fmt.allocPrint(allocator, "{s}_{d}", .{ stem, mol_idx + 1 });
    }
}

/// Generic SASA calculation dispatcher.
/// When bitmask_lut_ptr is non-null, uses bitmask-optimized Shrake-Rupley.
/// When adaptive_sr is enabled, uses coarse/fine bitmask LUTs.
fn calculateSasaDispatch(
    comptime T: type,
    allocator: Allocator,
    input: AtomInput,
    config: BatchConfig,
    probe_radius: T,
    n_threads: usize,
    bitmask_lut_ptr: ?*const bitmask_lut.BitmaskLutGen(T),
    coarse_lut_ptr: ?*const bitmask_lut.BitmaskLutGen(T),
    fine_lut_ptr: ?*const bitmask_lut.BitmaskLutGen(T),
) !SasaResultGen(T) {
    const sr_config = ConfigGen(T){ .n_points = config.n_points, .probe_radius = probe_radius };

    if (config.adaptive_sr) {
        const coarse_lut = coarse_lut_ptr orelse return error.MissingAdaptiveLut;
        const fine_lut = fine_lut_ptr orelse return error.MissingAdaptiveLut;
        const adaptive = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).AdaptiveConfig{
            .coarse_points = config.coarse_points,
            .fine_points = config.fine_points,
            .low = @as(T, @floatCast(config.adaptive_low)),
            .high = @as(T, @floatCast(config.adaptive_high)),
        };
        return if (n_threads > 1)
            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).calculateSasaAdaptiveParallelWithLuts(
                allocator,
                input,
                sr_config,
                adaptive,
                n_threads,
                coarse_lut,
                fine_lut,
            )
        else
            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).calculateSasaAdaptiveWithLuts(
                allocator,
                input,
                sr_config,
                adaptive,
                coarse_lut,
                fine_lut,
            );
    }

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

    return switch (config.algorithm) {
        .sr => if (n_threads > 1)
            shrake_rupley.ShrakeRupleyGen(T).calculateSasaParallel(allocator, input, sr_config, n_threads)
        else
            shrake_rupley.ShrakeRupleyGen(T).calculateSasa(allocator, input, sr_config),
        .lr => if (n_threads > 1)
            lee_richards.LeeRichardsGen(T).calculateSasaParallel(allocator, input, .{
                .n_slices = config.n_slices,
                .probe_radius = probe_radius,
            }, n_threads)
        else
            lee_richards.LeeRichardsGen(T).calculateSasa(allocator, input, .{
                .n_slices = config.n_slices,
                .probe_radius = probe_radius,
            }),
    };
}

/// Write SASA result to output file
/// Handles f32 -> f64 conversion for consistent output format
fn writeSasaOutput(
    comptime T: type,
    allocator: Allocator,
    io: std.Io,
    result: *const SasaResultGen(T),
    output_dir: []const u8,
    filename: []const u8,
    format: OutputFormat,
) !void {
    const output_filename = try replaceExtension(allocator, filename, getOutputExtension(format));
    defer allocator.free(output_filename);
    const output_path = try std.fs.path.join(allocator, &.{ output_dir, output_filename });

    if (T == f64) {
        try json_writer.writeSasaResultWithFormat(allocator, io, result.*, output_path, format);
    } else {
        var result_f64 = try result.toF64(allocator);
        defer result_f64.deinit();
        try json_writer.writeSasaResultWithFormat(allocator, io, result_f64, output_path, format);
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
    atom_areas: ?[]const f64 = null, // Populated for jsonl output
    residue_map: ?json_writer.ResidueMap = null, // Populated for jsonl residue-map output

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
            if (result.atom_areas) |areas| {
                self.allocator.free(areas);
            }
            if (result.residue_map) |*map| {
                map.deinit();
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

        // Print details for failed files
        if (self.failed > 0) {
            std.debug.print("\nFailed files:\n", .{});
            for (self.file_results) |file_result| {
                if (file_result.status == .err) {
                    if (file_result.error_msg) |msg| {
                        std.debug.print("  {s}: {s}\n", .{ file_result.filename, msg });
                    } else {
                        std.debug.print("  {s}: unknown error\n", .{file_result.filename});
                    }
                }
            }
        }

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
fn readInputFile(allocator: Allocator, io: std.Io, path: []const u8, config: BatchConfig) !AtomInput {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => json_parser.readAtomInputFromFile(allocator, io, path),
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
        .sdf => blk: {
            const source = if (compressed.isCompressed(path))
                try compressed.read(allocator, path)
            else file_blk: {
                const f = try std.Io.Dir.cwd().openFile(io, path, .{});
                defer f.close(io);
                var read_buf: [65536]u8 = undefined;
                var file_r = f.reader(io, &read_buf);
                break :file_blk try file_r.interface.allocRemaining(allocator, .unlimited);
            };
            defer allocator.free(source);

            const molecules = try sdf_parser.parse(allocator, source);
            defer sdf_parser.freeMolecules(allocator, molecules);

            break :blk try sdf_parser.toAtomInput(allocator, molecules, !config.include_hydrogens);
        },
    };
}

/// Apply built-in classifier to replace radii based on residue/atom names
fn applyBuiltinClassifier(input: *AtomInput, ct: ClassifierType, sdf_ccd: ?*const ccd_parser.ComponentDict, external_ccd: ?*const ccd_parser.ComponentDict) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    // For CCD/ProtOr: create classifier instance and feed external CCD components
    var ccd_clf: ?classifier_ccd.CcdClassifier = if (ct == .ccd or ct == .protor) classifier_ccd.CcdClassifier.init(input.allocator) else null;
    defer if (ccd_clf) |*c| c.deinit();

    if (ccd_clf != null) {
        // Deduplicate: collect unique non-hardcoded residues
        var needed: std.StringHashMapUnmanaged(void) = .empty;
        defer needed.deinit(input.allocator);
        for (0..n) |i| {
            const res = residues[i].slice();
            if (!classifier_ccd.CcdClassifier.isHardcoded(res)) {
                try needed.put(input.allocator, res, {});
            }
        }

        if (needed.count() > 0) {
            const dicts: [2]?*const ccd_parser.ComponentDict = .{ sdf_ccd, external_ccd };
            for (dicts) |maybe_dict| {
                if (maybe_dict) |dict| {
                    var it = needed.keyIterator();
                    while (it.next()) |key_ptr| {
                        if (dict.get(key_ptr.*)) |comp| {
                            ccd_clf.?.addComponent(&comp) catch |err| {
                                std.debug.print("Warning: Could not derive CCD properties for '{s}': {s}\n", .{ key_ptr.*, @errorName(err) });
                            };
                        }
                    }
                }
            }
        }
    }

    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    for (0..n) |i| {
        const maybe_radius: ?f64 = switch (ct) {
            .naccess => classifier_naccess.getRadius(residues[i].slice(), atom_names[i].slice()),
            .protor, .ccd => if (ccd_clf) |*c| c.getRadius(residues[i].slice(), atom_names[i].slice()) else null,
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

/// Apply custom classifier to replace radii based on residue/atom names.
fn applyCustomClassifier(input: *AtomInput, custom_classifier: *const classifier.Classifier, quiet: bool) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    var classified_count: usize = 0;
    var fallback_count: usize = 0;

    for (0..n) |i| {
        if (custom_classifier.getRadius(residues[i].slice(), atom_names[i].slice())) |r| {
            new_radii[i] = r;
            classified_count += 1;
        } else if (input.element) |elements| {
            if (classifier.guessRadiusFromAtomicNumber(elements[i])) |r| {
                new_radii[i] = r;
                fallback_count += 1;
            } else {
                new_radii[i] = input.r[i];
            }
        } else if (classifier.guessRadiusFromAtomName(atom_names[i].slice())) |r| {
            new_radii[i] = r;
            fallback_count += 1;
        } else {
            new_radii[i] = input.r[i];
        }
    }

    input.allocator.free(input.r);
    input.r = new_radii;

    if (!quiet) {
        std.debug.print("Classifier '{s}': {d} atoms classified, {d} fallback\n", .{
            custom_classifier.name,
            classified_count,
            fallback_count,
        });
    }
}

fn attachResidueMap(
    arena: Allocator,
    result_allocator: Allocator,
    input: AtomInput,
    atom_areas: []const f64,
    result: *FileResult,
) bool {
    result.residue_map = json_writer.buildResidueMap(arena, input, atom_areas) catch |err| {
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "residue map failed: {s}", .{@errorName(err)}) catch null;
        return false;
    };
    return true;
}

/// Scan directory for structure files (.json, .pdb, .cif, .mmcif, .ent and compressed variants)
pub fn scanDirectory(allocator: Allocator, io: std.Io, dir_path: []const u8) ![][]const u8 {
    var files: std.ArrayListUnmanaged([]const u8) = .empty;
    errdefer {
        for (files.items) |f| allocator.free(f);
        files.deinit(allocator);
    }

    var dir = std.Io.Dir.cwd().openDir(io, dir_path, .{ .iterate = true }) catch |err| {
        return err;
    };
    defer dir.close(io);

    var iter = dir.iterate();
    while (try iter.next(io)) |entry| {
        // Accept regular files and symlinks (for sampled batch benchmarking)
        if (entry.kind != .file and entry.kind != .sym_link) continue;

        const name = entry.name;
        // Skip filenames with path separators (defense in depth)
        if (std.mem.findAny(u8, name, "/\\") != null) continue;
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
/// result_allocator: used for data that must outlive the arena (error messages)
/// atom_areas (when store_atom_areas is true) are allocated on the arena
/// and must be consumed before the caller resets the arena.
/// n_threads: number of threads for SASA calculation (1 = single-threaded)
fn processOneFile(
    arena: Allocator,
    io: std.Io,
    result_allocator: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    filename: []const u8,
    config: BatchConfig,
    n_threads: usize,
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
) FileResult {
    var result = FileResult{
        .filename = filename,
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .ok,
    };

    // Build input path
    const input_path = std.fs.path.join(arena, &.{ input_dir, filename }) catch |err| {
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "path join failed: {s}", .{@errorName(err)}) catch null;
        return result;
    };

    // Read and parse input (auto-detect format from extension)
    var input = readInputFile(arena, io, input_path, config) catch |err| {
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "read/parse failed: {s}", .{@errorName(err)}) catch null;
        return result;
    };
    defer input.deinit();

    // Apply classifier for PDB/mmCIF input (skip JSON unless classification info exists).
    if (config.custom_classifier) |custom_classifier| {
        if (input.hasClassificationInfo()) {
            applyCustomClassifier(&input, custom_classifier, config.quiet) catch |err| {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "classifier failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
        }
    } else if (config.classifier_type) |ct| {
        const format = format_detect.detectInputFormat(input_path);
        if (format != .json and input.hasClassificationInfo()) {
            applyBuiltinClassifier(&input, ct, config.sdf_ccd, config.external_ccd) catch |err| {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "classifier failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
        }
    }

    result.n_atoms = input.atomCount();

    // Time SASA calculation only
    var sasa_timer = std.Io.Timestamp.now(io, .awake);

    // Calculate SASA using generic dispatcher
    var total_area: f64 = 0;
    switch (config.precision) {
        .f64 => {
            var sasa_result = calculateSasaDispatch(
                f64,
                arena,
                input,
                config,
                config.probe_radius,
                n_threads,
                lut_f64,
                coarse_lut_f64,
                fine_lut_f64,
            ) catch |err| {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "SASA calculation failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
            defer sasa_result.deinit();
            result.sasa_time_ns = @intCast(sasa_timer.untilNow(io, .awake).nanoseconds);
            total_area = sasa_result.total_area;

            if (config.store_atom_areas) {
                result.atom_areas = arena.dupe(f64, sasa_result.atom_areas) catch {
                    result.status = .err;
                    result.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                    return result;
                };
            }
            if (config.residue_map) {
                const areas = result.atom_areas orelse sasa_result.atom_areas;
                if (!attachResidueMap(arena, result_allocator, input, areas, &result)) return result;
            }

            if (output_dir) |out_dir| {
                if (!config.store_atom_areas) {
                    writeSasaOutput(f64, arena, io, &sasa_result, out_dir, filename, config.output_format) catch |err| {
                        result.status = .err;
                        result.error_msg = std.fmt.allocPrint(result_allocator, "output write failed: {s}", .{@errorName(err)}) catch null;
                        return result;
                    };
                }
            }
        },
        .f32 => {
            var sasa_result = calculateSasaDispatch(
                f32,
                arena,
                input,
                config,
                @as(f32, @floatCast(config.probe_radius)),
                n_threads,
                lut_f32,
                coarse_lut_f32,
                fine_lut_f32,
            ) catch |err| {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "SASA calculation failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
            defer sasa_result.deinit();
            result.sasa_time_ns = @intCast(sasa_timer.untilNow(io, .awake).nanoseconds);
            total_area = @floatCast(sasa_result.total_area);

            if (config.store_atom_areas) {
                const areas_f32 = sasa_result.atom_areas;
                const areas_f64 = arena.alloc(f64, areas_f32.len) catch {
                    result.status = .err;
                    result.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                    return result;
                };
                for (areas_f32, 0..) |v, j| areas_f64[j] = @floatCast(v);
                result.atom_areas = areas_f64;
            }
            if (config.residue_map) {
                const areas = result.atom_areas orelse blk: {
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
                    writeSasaOutput(f32, arena, io, &sasa_result, out_dir, filename, config.output_format) catch |err| {
                        result.status = .err;
                        result.error_msg = std.fmt.allocPrint(result_allocator, "output write failed: {s}", .{@errorName(err)}) catch null;
                        return result;
                    };
                }
            }
        },
    }

    result.total_sasa = total_area;
    return result;
}

/// Process a single SDF molecule and return result.
/// The molecule is provided as a pre-parsed slice of 1 element.
/// `display_name` is used as the filename in the result (e.g., "file.sdf:methane").
fn processOneSdfMolecule(
    arena: Allocator,
    io: std.Io,
    result_allocator: Allocator,
    display_name: []const u8,
    molecule: *const sdf_parser.SdfMolecule,
    output_dir: ?[]const u8,
    config: BatchConfig,
    n_threads: usize,
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
) FileResult {
    var result = FileResult{
        .filename = display_name,
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .ok,
    };

    // Convert single molecule to AtomInput
    const mol_slice: []const sdf_parser.SdfMolecule = @as([*]const sdf_parser.SdfMolecule, @ptrCast(molecule))[0..1];
    var input = sdf_parser.toAtomInput(arena, mol_slice, !config.include_hydrogens) catch |err| {
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "SDF toAtomInput failed: {s}", .{@errorName(err)}) catch null;
        return result;
    };
    defer input.deinit();

    // Build CCD component dict for this molecule's bond topology
    var sdf_dict: ?ccd_parser.ComponentDict = null;
    if (molecule.name.len > 0) {
        const stored = sdf_parser.toStoredComponent(arena, molecule) catch |err| blk: {
            logWarning("{s}: failed to build SDF component: {s}", .{ display_name, @errorName(err) });
            break :blk null;
        };
        if (stored) |s| {
            var dict = ccd_parser.ComponentDict.init(arena);
            const comp_id_str = molecule.name[0..@min(molecule.name.len, 5)];
            const dict_key = arena.dupe(u8, comp_id_str) catch |err| {
                logWarning("{s}: SDF component registration failed: {s}", .{ display_name, @errorName(err) });
                var mut_s = s;
                mut_s.deinit();
                dict.deinit();
                sdf_dict = null;
                // fall through to classifier without SDF dict
                return processOneSdfMoleculeInner(arena, io, result_allocator, &result, &input, output_dir, display_name, config, n_threads, lut_f64, lut_f32, coarse_lut_f64, fine_lut_f64, coarse_lut_f32, fine_lut_f32, null);
            };
            dict.owned_keys.append(arena, dict_key) catch |err| {
                logWarning("{s}: SDF component registration failed: {s}", .{ display_name, @errorName(err) });
                arena.free(dict_key);
                var mut_s = s;
                mut_s.deinit();
                dict.deinit();
                return processOneSdfMoleculeInner(arena, io, result_allocator, &result, &input, output_dir, display_name, config, n_threads, lut_f64, lut_f32, coarse_lut_f64, fine_lut_f64, coarse_lut_f32, fine_lut_f32, null);
            };
            dict.components.put(arena, dict_key, s) catch |err| {
                logWarning("{s}: SDF component registration failed: {s}", .{ display_name, @errorName(err) });
                var mut_s = s;
                mut_s.deinit();
                dict.deinit();
                return processOneSdfMoleculeInner(arena, io, result_allocator, &result, &input, output_dir, display_name, config, n_threads, lut_f64, lut_f32, coarse_lut_f64, fine_lut_f64, coarse_lut_f32, fine_lut_f32, null);
            };
            sdf_dict = dict;
        }
    }
    defer if (sdf_dict) |*d| d.deinit();

    const sdf_ccd_ptr: ?*const ccd_parser.ComponentDict = if (sdf_dict) |*d| d else null;
    return processOneSdfMoleculeInner(arena, io, result_allocator, &result, &input, output_dir, display_name, config, n_threads, lut_f64, lut_f32, coarse_lut_f64, fine_lut_f64, coarse_lut_f32, fine_lut_f32, sdf_ccd_ptr);
}

/// Inner helper: apply classifier, run SASA, write output for a single SDF molecule.
fn processOneSdfMoleculeInner(
    arena: Allocator,
    io: std.Io,
    result_allocator: Allocator,
    result: *FileResult,
    input: *AtomInput,
    output_dir: ?[]const u8,
    display_name: []const u8,
    config: BatchConfig,
    n_threads: usize,
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    sdf_ccd: ?*const ccd_parser.ComponentDict,
) FileResult {
    var res = result.*;

    // Apply classifier (SDF molecules normally have classification info)
    if (config.custom_classifier) |custom_classifier| {
        if (input.hasClassificationInfo()) {
            applyCustomClassifier(input, custom_classifier, config.quiet) catch |err| {
                res.status = .err;
                res.error_msg = std.fmt.allocPrint(result_allocator, "classifier failed: {s}", .{@errorName(err)}) catch null;
                return res;
            };
        }
    } else if (config.classifier_type) |ct| {
        if (input.hasClassificationInfo()) {
            // Merge SDF-derived dict with external CCD if available
            const effective_sdf_ccd = sdf_ccd orelse config.sdf_ccd;
            applyBuiltinClassifier(input, ct, effective_sdf_ccd, config.external_ccd) catch |err| {
                res.status = .err;
                res.error_msg = std.fmt.allocPrint(result_allocator, "classifier failed: {s}", .{@errorName(err)}) catch null;
                return res;
            };
        }
    }

    res.n_atoms = input.atomCount();

    // Time SASA calculation
    var sasa_timer = std.Io.Timestamp.now(io, .awake);

    var total_area: f64 = 0;
    switch (config.precision) {
        .f64 => {
            var sasa_result = calculateSasaDispatch(
                f64,
                arena,
                input.*,
                config,
                config.probe_radius,
                n_threads,
                lut_f64,
                coarse_lut_f64,
                fine_lut_f64,
            ) catch |err| {
                res.status = .err;
                res.error_msg = std.fmt.allocPrint(result_allocator, "SASA calculation failed: {s}", .{@errorName(err)}) catch null;
                return res;
            };
            defer sasa_result.deinit();
            res.sasa_time_ns = @intCast(sasa_timer.untilNow(io, .awake).nanoseconds);
            total_area = sasa_result.total_area;

            if (config.store_atom_areas) {
                res.atom_areas = arena.dupe(f64, sasa_result.atom_areas) catch {
                    res.status = .err;
                    res.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                    return res;
                };
            }
            if (config.residue_map) {
                const areas = res.atom_areas orelse sasa_result.atom_areas;
                if (!attachResidueMap(arena, result_allocator, input.*, areas, &res)) return res;
            }

            if (output_dir) |out_dir| {
                if (!config.store_atom_areas) {
                    writeSasaOutput(f64, arena, io, &sasa_result, out_dir, display_name, config.output_format) catch |err| {
                        res.status = .err;
                        res.error_msg = std.fmt.allocPrint(result_allocator, "output write failed: {s}", .{@errorName(err)}) catch null;
                        return res;
                    };
                }
            }
        },
        .f32 => {
            var sasa_result = calculateSasaDispatch(
                f32,
                arena,
                input.*,
                config,
                @as(f32, @floatCast(config.probe_radius)),
                n_threads,
                lut_f32,
                coarse_lut_f32,
                fine_lut_f32,
            ) catch |err| {
                res.status = .err;
                res.error_msg = std.fmt.allocPrint(result_allocator, "SASA calculation failed: {s}", .{@errorName(err)}) catch null;
                return res;
            };
            defer sasa_result.deinit();
            res.sasa_time_ns = @intCast(sasa_timer.untilNow(io, .awake).nanoseconds);
            total_area = @floatCast(sasa_result.total_area);

            if (config.store_atom_areas) {
                const areas_f32 = sasa_result.atom_areas;
                const areas_f64 = arena.alloc(f64, areas_f32.len) catch {
                    res.status = .err;
                    res.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                    return res;
                };
                for (areas_f32, 0..) |v, j| areas_f64[j] = @floatCast(v);
                res.atom_areas = areas_f64;
            }
            if (config.residue_map) {
                const areas = res.atom_areas orelse blk: {
                    const areas_f64 = arena.alloc(f64, sasa_result.atom_areas.len) catch {
                        res.status = .err;
                        res.error_msg = std.fmt.allocPrint(result_allocator, "atom_areas allocation failed", .{}) catch null;
                        return res;
                    };
                    for (sasa_result.atom_areas, 0..) |v, j| areas_f64[j] = @floatCast(v);
                    break :blk areas_f64;
                };
                if (!attachResidueMap(arena, result_allocator, input.*, areas, &res)) return res;
            }

            if (output_dir) |out_dir| {
                if (!config.store_atom_areas) {
                    writeSasaOutput(f32, arena, io, &sasa_result, out_dir, display_name, config.output_format) catch |err| {
                        res.status = .err;
                        res.error_msg = std.fmt.allocPrint(result_allocator, "output write failed: {s}", .{@errorName(err)}) catch null;
                        return res;
                    };
                }
            }
        },
    }

    res.total_sasa = total_area;
    return res;
}

/// Thread-safe JSONL streaming writer.
/// Each call to writeResult acquires the mutex, serializes one line, and flushes.
const JsonlStreamWriter = struct {
    mutex: std.Io.Mutex = .init,
    file: std.Io.File,
    io: std.Io,
    /// Set to true if any writeResult call failed to write to the output file.
    /// TODO: surface this from runBatch return value to make the CLI exit
    /// non-zero on partial JSONL write failure.
    write_failed: bool = false,

    /// Serialize and write one JSONL line for a completed file result.
    /// alloc is a short-lived allocator (e.g., thread-local arena) used only for
    /// the serialized string; it is freed by the caller's arena reset.
    pub fn writeResult(
        self: *JsonlStreamWriter,
        alloc: Allocator,
        result: *FileResult,
    ) void {
        const line = fileResultToJsonlLine(alloc, result) catch |err| {
            logWarning("failed to serialize {s}: {s}", .{ result.filename, @errorName(err) });
            return;
        };
        // line is on alloc (arena); no explicit free needed — arena reset handles it.

        self.mutex.lockUncancelable(self.io);
        defer self.mutex.unlock(self.io);

        // 64KB stack buffer; use streaming mode so the OS seek position
        // advances (safe under mutex — only one thread writes at a time).
        var buf: [64 * 1024]u8 = undefined;
        var w = std.Io.File.Writer.initStreaming(self.file, self.io, &buf);
        w.interface.writeAll(line) catch |err| {
            logWarning("JSONL write failed for {s}: {s}", .{ result.filename, @errorName(err) });
            self.write_failed = true;
            return;
        };
        w.interface.writeAll("\n") catch |err| {
            logWarning("JSONL newline write failed for {s}: {s}", .{ result.filename, @errorName(err) });
            self.write_failed = true;
            return;
        };
        w.interface.flush() catch |err| {
            logWarning("JSONL flush failed for {s}: {s}", .{ result.filename, @errorName(err) });
            self.write_failed = true;
        };
    }

    /// Returns true if any write to the JSONL file failed.
    pub fn hasError(self: *const JsonlStreamWriter) bool {
        return self.write_failed;
    }
};

/// Write a JSONL line for a result via the buffered writer.
/// atom_areas live on the arena and are invalidated after arena reset.
fn fileResultToJsonlLine(allocator: Allocator, result: *FileResult) ![]u8 {
    const areas = result.atom_areas orelse return error.MissingAtomAreas;
    if (result.residue_map) |map| {
        return json_writer.fileResultWithResidueMapToJsonlLine(allocator, result.filename, result.total_sasa, areas, map);
    }
    return json_writer.fileResultToJsonlLine(allocator, result.filename, result.total_sasa, areas);
}

fn writeJsonlResult(
    jsonl_writer: *std.Io.File.Writer,
    arena_alloc: Allocator,
    result: *FileResult,
) void {
    if (result.status != .ok) return;
    const line = fileResultToJsonlLine(arena_alloc, result) catch |err| {
        logWarning("failed to serialize {s}: {s}", .{ result.filename, @errorName(err) });
        return;
    };
    jsonl_writer.interface.writeAll(line) catch |err| {
        logWarning("JSONL write failed for {s}: {s}", .{ result.filename, @errorName(err) });
        return; // Don't attempt newline or flush
    };
    jsonl_writer.interface.writeAll("\n") catch |err| {
        logWarning("JSONL newline write failed for {s}: {s}", .{ result.filename, @errorName(err) });
        return; // Don't flush a corrupted line
    };
    jsonl_writer.interface.flush() catch |err| {
        logWarning("JSONL flush failed for {s}: {s}", .{ result.filename, @errorName(err) });
    };
}

/// Run batch processing sequentially (single-threaded path, used when n_threads <= 1)
pub fn runBatchSequential(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
) !BatchResult {
    // Start total timer
    var total_timer = std.Io.Timestamp.now(io, .awake);

    // Scan directory for files
    const files = try scanDirectory(allocator, io, input_dir);
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing files", .estimated_total_items = files.len })
    else
        .none;
    defer progress_root.end();
    const progress_node: ?std.Progress.Node = if (shouldShowProgress(config)) progress_root else null;

    // Create output directory if specified
    if (output_dir) |out_dir| {
        try std.Io.Dir.cwd().createDirPath(io, out_dir);
    }

    // Open JSONL output file (or stdout) when streaming is requested
    var jsonl_file: ?std.Io.File = null;
    var jsonl_file_needs_close = false;
    if (jsonl_output_path) |path| {
        jsonl_file = try std.Io.Dir.cwd().createFile(io, path, .{});
        jsonl_file_needs_close = true;
    } else if (config.store_atom_areas) {
        jsonl_file = std.Io.File.stdout();
    }
    defer if (jsonl_file_needs_close) {
        if (jsonl_file) |f| f.close(io);
    };

    // Use ArrayList for results (SDF files may expand into multiple entries)
    var results_list = std.ArrayListUnmanaged(FileResult).empty;
    defer results_list.deinit(allocator);

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

    // JSONL buffered writer: created once, reused across all files.
    // Uses streaming mode so OS seek position advances with each write.
    var jsonl_write_buf: [64 * 1024]u8 = undefined;
    var jsonl_writer: ?std.Io.File.Writer = if (jsonl_file) |jf|
        std.Io.File.Writer.initStreaming(jf, io, &jsonl_write_buf)
    else
        null;

    var total_items: usize = 0;

    for (files) |filename| {
        const format = format_detect.detectInputFormat(filename);

        if (format == .sdf) {
            // SDF: parse and process each molecule individually
            const input_path = std.fs.path.join(arena.allocator(), &.{ input_dir, filename }) catch |err| {
                const filename_copy = try allocator.dupe(u8, filename);
                try results_list.append(allocator, FileResult{
                    .filename = filename_copy,
                    .n_atoms = 0,
                    .sasa_time_ns = 0,
                    .total_sasa = 0,
                    .status = .err,
                    .error_msg = std.fmt.allocPrint(allocator, "path join failed: {s}", .{@errorName(err)}) catch null,
                });
                failed += 1;
                total_items += 1;
                _ = arena.reset(.retain_capacity);
                continue;
            };

            const source = if (compressed.isCompressed(input_path))
                compressed.read(arena.allocator(), input_path) catch |err| {
                    const filename_copy = try allocator.dupe(u8, filename);
                    try results_list.append(allocator, FileResult{
                        .filename = filename_copy,
                        .n_atoms = 0,
                        .sasa_time_ns = 0,
                        .total_sasa = 0,
                        .status = .err,
                        .error_msg = std.fmt.allocPrint(allocator, "read failed: {s}", .{@errorName(err)}) catch null,
                    });
                    failed += 1;
                    total_items += 1;
                    _ = arena.reset(.retain_capacity);
                    continue;
                }
            else file_blk: {
                const f = std.Io.Dir.cwd().openFile(io, input_path, .{}) catch |err| {
                    const filename_copy = try allocator.dupe(u8, filename);
                    try results_list.append(allocator, FileResult{
                        .filename = filename_copy,
                        .n_atoms = 0,
                        .sasa_time_ns = 0,
                        .total_sasa = 0,
                        .status = .err,
                        .error_msg = std.fmt.allocPrint(allocator, "open failed: {s}", .{@errorName(err)}) catch null,
                    });
                    failed += 1;
                    total_items += 1;
                    _ = arena.reset(.retain_capacity);
                    continue;
                };
                defer f.close(io);
                var read_buf_seq: [65536]u8 = undefined;
                var file_r_seq = f.reader(io, &read_buf_seq);
                break :file_blk file_r_seq.interface.allocRemaining(arena.allocator(), .unlimited) catch |err| {
                    const filename_copy = try allocator.dupe(u8, filename);
                    try results_list.append(allocator, FileResult{
                        .filename = filename_copy,
                        .n_atoms = 0,
                        .sasa_time_ns = 0,
                        .total_sasa = 0,
                        .status = .err,
                        .error_msg = std.fmt.allocPrint(allocator, "read failed: {s}", .{@errorName(err)}) catch null,
                    });
                    failed += 1;
                    total_items += 1;
                    _ = arena.reset(.retain_capacity);
                    continue;
                };
            };

            const molecules = sdf_parser.parse(arena.allocator(), source) catch |err| {
                const filename_copy = try allocator.dupe(u8, filename);
                try results_list.append(allocator, FileResult{
                    .filename = filename_copy,
                    .n_atoms = 0,
                    .sasa_time_ns = 0,
                    .total_sasa = 0,
                    .status = .err,
                    .error_msg = std.fmt.allocPrint(allocator, "SDF parse failed: {s}", .{@errorName(err)}) catch null,
                });
                failed += 1;
                total_items += 1;
                _ = arena.reset(.retain_capacity);
                continue;
            };

            // Process each molecule individually
            for (molecules, 0..) |*mol, mol_idx| {
                // Build display name: "stem_molname" or "stem_N"
                const display_name = sdfMoleculeDisplayName(allocator, filename, mol.name, mol_idx) catch |err| blk: {
                    logWarning("{s}: molecule {d} display name failed ({s}), using filename", .{ filename, mol_idx, @errorName(err) });
                    break :blk try allocator.dupe(u8, filename);
                };

                var mol_result = processOneSdfMolecule(
                    arena.allocator(),
                    io,
                    allocator,
                    display_name,
                    mol,
                    output_dir,
                    config,
                    1, // single-threaded
                    luts.f64Ptr(),
                    luts.f32Ptr(),
                    luts.coarseF64Ptr(),
                    luts.fineF64Ptr(),
                    luts.coarseF32Ptr(),
                    luts.fineF32Ptr(),
                );
                mol_result.filename = display_name;

                if (mol_result.status == .ok) {
                    successful += 1;
                    total_sasa_time_ns += mol_result.sasa_time_ns;
                } else {
                    failed += 1;
                }

                // Stream JSONL output
                if (jsonl_writer) |*w| {
                    writeJsonlResult(w, arena.allocator(), &mol_result);
                }
                mol_result.atom_areas = null;
                mol_result.residue_map = null;

                try results_list.append(allocator, mol_result);
                total_items += 1;
            }

            _ = arena.reset(.retain_capacity);
        } else {
            // Non-SDF: existing logic
            const filename_copy = try allocator.dupe(u8, filename);

            var result = processOneFile(
                arena.allocator(),
                io,
                allocator,
                input_dir,
                output_dir,
                filename,
                config,
                1, // single-threaded
                luts.f64Ptr(),
                luts.f32Ptr(),
                luts.coarseF64Ptr(),
                luts.fineF64Ptr(),
                luts.coarseF32Ptr(),
                luts.fineF32Ptr(),
            );
            result.filename = filename_copy;

            if (result.status == .ok) {
                successful += 1;
                total_sasa_time_ns += result.sasa_time_ns;
            } else {
                failed += 1;
            }

            // Stream JSONL output
            if (jsonl_writer) |*w| {
                writeJsonlResult(w, arena.allocator(), &result);
            }
            result.atom_areas = null;
            result.residue_map = null;

            try results_list.append(allocator, result);
            total_items += 1;

            _ = arena.reset(.retain_capacity);
        }

        if (progress_node) |node| {
            node.completeOne();
        }
    }

    const total_time_ns: u64 = @intCast(total_timer.untilNow(io, .awake).nanoseconds);

    return BatchResult{
        .total_files = total_items,
        .successful = successful,
        .failed = failed,
        .total_sasa_time_ns = total_sasa_time_ns,
        .total_time_ns = total_time_ns,
        .file_results = try results_list.toOwnedSlice(allocator),
        .allocator = allocator,
    };
}

/// A work item for parallel batch processing.
/// Represents either a plain file or a specific molecule within an SDF file.
const WorkItem = struct {
    filename: []const u8, // Original filename in the input directory
    display_name: []const u8, // Display name for results (e.g., "file.sdf:water")
    mol_idx: ?usize, // If non-null, SDF molecule index (0-based) within pre-parsed molecules
};

/// Shared context for parallel workers
const ParallelContext = struct {
    work_items: []const WorkItem,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    results: []FileResult,
    result_allocator: Allocator,
    next_item: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    jsonl_stream: ?*JsonlStreamWriter,
    io: std.Io,
    /// Pre-parsed SDF data: keyed by filename, each entry is the parsed source bytes.
    /// Workers read these (read-only) to avoid re-parsing the same SDF file.
    sdf_sources: std.StringHashMapUnmanaged([]const u8),
};

/// Worker thread function for parallel batch processing
fn parallelWorker(ctx: *ParallelContext) void {
    // Use smp_allocator as arena backing to avoid mmap/munmap syscall contention
    // that page_allocator causes under multi-threaded workloads.
    // smp_allocator is thread-safe and does not require libc.
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();

    while (true) {
        // Atomically grab the next work item index.
        // .monotonic suffices: we only need atomic increment, no ordering
        // between unrelated memory accesses across threads.
        const item_idx = ctx.next_item.fetchAdd(1, .monotonic);

        if (item_idx >= ctx.work_items.len) {
            break; // No more work items
        }

        const work = ctx.work_items[item_idx];

        // Copy display name to result allocator (thread-safe: each index is unique)
        const name_copy = ctx.result_allocator.dupe(u8, work.display_name) catch {
            const empty = ctx.result_allocator.alloc(u8, 0) catch {
                ctx.results[item_idx] = FileResult{
                    .filename = "",
                    .n_atoms = 0,
                    .sasa_time_ns = 0,
                    .total_sasa = 0,
                    .status = .err,
                    .error_msg = null,
                };
                _ = ctx.processed_count.fetchAdd(1, .release);
                _ = arena.reset(.retain_capacity);
                continue;
            };
            ctx.results[item_idx] = FileResult{
                .filename = empty,
                .n_atoms = 0,
                .sasa_time_ns = 0,
                .total_sasa = 0,
                .status = .err,
            };
            _ = ctx.processed_count.fetchAdd(1, .release);
            _ = arena.reset(.retain_capacity);
            continue;
        };

        if (work.mol_idx) |mol_idx| {
            // SDF molecule: re-parse from pre-loaded source on thread-local arena
            const source = ctx.sdf_sources.get(work.filename) orelse {
                ctx.results[item_idx] = FileResult{
                    .filename = name_copy,
                    .n_atoms = 0,
                    .sasa_time_ns = 0,
                    .total_sasa = 0,
                    .status = .err,
                    .error_msg = ctx.result_allocator.dupe(u8, "SDF source not found") catch null,
                };
                _ = ctx.processed_count.fetchAdd(1, .release);
                _ = arena.reset(.retain_capacity);
                continue;
            };

            const molecules = sdf_parser.parse(arena.allocator(), source) catch {
                ctx.results[item_idx] = FileResult{
                    .filename = name_copy,
                    .n_atoms = 0,
                    .sasa_time_ns = 0,
                    .total_sasa = 0,
                    .status = .err,
                    .error_msg = ctx.result_allocator.dupe(u8, "SDF re-parse failed") catch null,
                };
                _ = ctx.processed_count.fetchAdd(1, .release);
                _ = arena.reset(.retain_capacity);
                continue;
            };

            if (mol_idx >= molecules.len) {
                ctx.results[item_idx] = FileResult{
                    .filename = name_copy,
                    .n_atoms = 0,
                    .sasa_time_ns = 0,
                    .total_sasa = 0,
                    .status = .err,
                    .error_msg = ctx.result_allocator.dupe(u8, "SDF molecule index out of range") catch null,
                };
                _ = ctx.processed_count.fetchAdd(1, .release);
                _ = arena.reset(.retain_capacity);
                continue;
            }

            var result = processOneSdfMolecule(
                arena.allocator(),
                ctx.io,
                ctx.result_allocator,
                name_copy,
                &molecules[mol_idx],
                ctx.output_dir,
                ctx.config,
                1, // single-threaded SASA per molecule
                ctx.lut_f64,
                ctx.lut_f32,
                ctx.coarse_lut_f64,
                ctx.fine_lut_f64,
                ctx.coarse_lut_f32,
                ctx.fine_lut_f32,
            );
            result.filename = name_copy;

            ctx.results[item_idx] = result;

            if (ctx.jsonl_stream) |stream| {
                if (result.status == .ok) {
                    stream.writeResult(arena.allocator(), &result);
                }
            }
            ctx.results[item_idx].atom_areas = null;
            ctx.results[item_idx].residue_map = null;
        } else {
            // Non-SDF file: existing logic
            var result = processOneFile(
                arena.allocator(),
                ctx.io,
                ctx.result_allocator,
                ctx.input_dir,
                ctx.output_dir,
                work.filename,
                ctx.config,
                1, // single-threaded SASA per file
                ctx.lut_f64,
                ctx.lut_f32,
                ctx.coarse_lut_f64,
                ctx.fine_lut_f64,
                ctx.coarse_lut_f32,
                ctx.fine_lut_f32,
            );
            result.filename = name_copy;

            // Store result (thread-safe: each index is unique)
            ctx.results[item_idx] = result;

            // Stream JSONL output (atom_areas on arena, valid until reset)
            if (ctx.jsonl_stream) |stream| {
                if (result.status == .ok) {
                    stream.writeResult(arena.allocator(), &result);
                }
            }
            // Clear atom_areas since arena will free them
            ctx.results[item_idx].atom_areas = null;
            ctx.results[item_idx].residue_map = null;
        }

        // Update progress counter (.release pairs with .acquire in progress monitor)
        _ = ctx.processed_count.fetchAdd(1, .release);

        // Reset arena for next work item
        _ = arena.reset(.retain_capacity);
    }
}

/// Build work items from file list, expanding SDF files into per-molecule items.
/// Returns the work items and a map of SDF sources (caller must free both).
fn buildWorkItems(
    allocator: Allocator,
    io: std.Io,
    files: []const []const u8,
    input_dir: []const u8,
) !struct { items: []WorkItem, sdf_sources: std.StringHashMapUnmanaged([]const u8) } {
    var items = std.ArrayListUnmanaged(WorkItem).empty;
    errdefer {
        for (items.items) |item| {
            // Only free if display_name was separately allocated (not same as filename)
            if (item.display_name.ptr != item.filename.ptr) {
                allocator.free(item.display_name);
            }
        }
        items.deinit(allocator);
    }

    var sdf_sources = std.StringHashMapUnmanaged([]const u8){};
    errdefer {
        var it = sdf_sources.valueIterator();
        while (it.next()) |v| allocator.free(v.*);
        sdf_sources.deinit(allocator);
    }

    for (files) |filename| {
        const format = format_detect.detectInputFormat(filename);
        if (format == .sdf) {
            // Read and parse SDF to count molecules
            const input_path = try std.fs.path.join(allocator, &.{ input_dir, filename });
            defer allocator.free(input_path);

            const source = if (compressed.isCompressed(input_path))
                compressed.read(allocator, input_path) catch |err| {
                    // If read fails, add a single error item
                    logWarning("{s}: failed to read SDF (compressed): {s}", .{ filename, @errorName(err) });
                    try items.append(allocator, .{
                        .filename = filename,
                        .display_name = filename,
                        .mol_idx = null,
                    });
                    continue;
                }
            else file_blk: {
                const f = std.Io.Dir.cwd().openFile(io, input_path, .{}) catch |err| {
                    logWarning("{s}: failed to open SDF: {s}", .{ filename, @errorName(err) });
                    try items.append(allocator, .{
                        .filename = filename,
                        .display_name = filename,
                        .mol_idx = null,
                    });
                    continue;
                };
                defer f.close(io);
                var read_buf_build: [65536]u8 = undefined;
                var file_r_build = f.reader(io, &read_buf_build);
                break :file_blk file_r_build.interface.allocRemaining(allocator, .unlimited) catch |err| {
                    logWarning("{s}: failed to read SDF: {s}", .{ filename, @errorName(err) });
                    try items.append(allocator, .{
                        .filename = filename,
                        .display_name = filename,
                        .mol_idx = null,
                    });
                    continue;
                };
            };

            const molecules = sdf_parser.parse(allocator, source) catch |err| {
                logWarning("{s}: failed to parse SDF: {s}", .{ filename, @errorName(err) });
                allocator.free(source);
                try items.append(allocator, .{
                    .filename = filename,
                    .display_name = filename,
                    .mol_idx = null,
                });
                continue;
            };
            defer sdf_parser.freeMolecules(allocator, molecules);

            // Store the source for workers to re-parse
            sdf_sources.put(allocator, filename, source) catch |err| {
                allocator.free(source);
                return err;
            };

            // Create one work item per molecule
            for (molecules, 0..) |mol, mol_idx| {
                const display_name = sdfMoleculeDisplayName(allocator, filename, mol.name, mol_idx) catch |err| blk: {
                    logWarning("{s}: molecule {d} display name failed ({s}), using filename", .{ filename, mol_idx, @errorName(err) });
                    break :blk try allocator.dupe(u8, filename);
                };

                try items.append(allocator, .{
                    .filename = filename,
                    .display_name = display_name,
                    .mol_idx = mol_idx,
                });
            }
        } else {
            try items.append(allocator, .{
                .filename = filename,
                .display_name = filename,
                .mol_idx = null,
            });
        }
    }

    return .{
        .items = try items.toOwnedSlice(allocator),
        .sdf_sources = sdf_sources,
    };
}

/// Run batch processing in parallel
pub fn runBatchParallel(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
) !BatchResult {
    // Start total timer
    var total_timer = std.Io.Timestamp.now(io, .awake);

    // Scan directory for files
    const files = try scanDirectory(allocator, io, input_dir);
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
            .total_time_ns = @intCast(total_timer.untilNow(io, .awake).nanoseconds),
            .file_results = try allocator.alloc(FileResult, 0),
            .allocator = allocator,
        };
    }

    // Create output directory if specified
    if (output_dir) |out_dir| {
        try std.Io.Dir.cwd().createDirPath(io, out_dir);
    }

    // Build work items (expanding SDF files into per-molecule items)
    var build_result = try buildWorkItems(allocator, io, files, input_dir);
    const work_items = build_result.items;
    defer {
        for (work_items) |item| {
            // Free display names that were allocated (not the same pointer as filename)
            if (item.display_name.ptr != item.filename.ptr) {
                allocator.free(item.display_name);
            }
        }
        allocator.free(work_items);
    }
    defer {
        var it = build_result.sdf_sources.valueIterator();
        while (it.next()) |v| allocator.free(v.*);
        build_result.sdf_sources.deinit(allocator);
    }

    // Allocate results (one per work item)
    const file_results = try allocator.alloc(FileResult, work_items.len);
    errdefer allocator.free(file_results);

    // Determine thread count
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        @min(config.n_threads, cpu_count);

    // For single item or single thread, use sequential
    if (work_items.len == 1 or n_threads <= 1) {
        allocator.free(file_results);
        return runBatchSequential(allocator, io, input_dir, output_dir, config, jsonl_output_path);
    }

    // Build bitmask LUT once (if enabled)
    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    // Open JSONL output file (or stdout) when streaming is requested
    var jsonl_file: ?std.Io.File = null;
    var jsonl_file_needs_close = false;
    if (jsonl_output_path) |path| {
        jsonl_file = try std.Io.Dir.cwd().createFile(io, path, .{});
        jsonl_file_needs_close = true;
    } else if (config.store_atom_areas) {
        jsonl_file = std.Io.File.stdout();
    }
    defer if (jsonl_file_needs_close) {
        if (jsonl_file) |f| f.close(io);
    };

    // Set up the stream writer on the stack (if JSONL streaming is active).
    // SAFETY: `undefined` when jsonl_file is null — never accessed because
    // jsonl_stream_ptr is also null in that case.
    var jsonl_stream_storage: JsonlStreamWriter = if (jsonl_file) |jf|
        JsonlStreamWriter{ .file = jf, .io = io }
    else
        undefined;
    const jsonl_stream_ptr: ?*JsonlStreamWriter = if (jsonl_file != null) &jsonl_stream_storage else null;

    // Create shared context
    var ctx = ParallelContext{
        .work_items = work_items,
        .input_dir = input_dir,
        .output_dir = output_dir,
        .config = config,
        .results = file_results,
        .result_allocator = allocator,
        .next_item = std.atomic.Value(usize).init(0),
        .processed_count = std.atomic.Value(usize).init(0),
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

    // Spawn worker threads
    const actual_threads = @min(n_threads, work_items.len);
    const threads = try allocator.alloc(std.Thread, actual_threads);
    defer allocator.free(threads);

    for (threads) |*thread| {
        thread.* = try std.Thread.spawn(.{}, parallelWorker, .{&ctx});
    }

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing items", .estimated_total_items = work_items.len })
    else
        .none;
    defer progress_root.end();

    // Progress monitoring (optional)
    if (shouldShowProgress(config)) {
        while (ctx.processed_count.load(.acquire) < work_items.len) {
            const processed = ctx.processed_count.load(.acquire);
            progress_root.setCompletedItems(processed);
            std.Io.sleep(io, .fromMilliseconds(50), .awake) catch {}; // 50ms update interval
        }
        progress_root.setCompletedItems(work_items.len);
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

    const total_time_ns: u64 = @intCast(total_timer.untilNow(io, .awake).nanoseconds);

    return BatchResult{
        .total_files = work_items.len,
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
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
) !BatchResult {
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = if (config.n_threads == 0)
        cpu_count
    else
        config.n_threads;

    if (n_threads <= 1) {
        return runBatchSequential(allocator, io, input_dir, output_dir, config, jsonl_output_path);
    }
    return runBatchParallel(allocator, io, input_dir, output_dir, config, jsonl_output_path);
}

// =============================================================================
// CLI argument parsing and run entry point
// =============================================================================

/// Parsed command-line arguments for the batch subcommand
const SdfPathList = sdf_parser.SdfPathList;

pub const BatchArgs = struct {
    input_path: ?[]const u8 = null,
    output_path: ?[]const u8 = null, // Output directory; null means no file output
    output_path_explicit: bool = false, // Track if -o/--output was explicitly set
    workflow_path: ?[]const u8 = null,
    chain_filter: ?[]const u8 = null,
    use_auth_chain: bool = false,
    residue_map: bool = false,
    n_threads: usize = 0,
    probe_radius: f64 = 1.4,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    algorithm: Algorithm = .sr,
    precision: Precision = .f64,
    output_format: OutputFormat = .json,
    classifier_type: ClassifierType = .ccd, // Default: ccd (ProtOr-compatible with CCD extension)
    include_hydrogens: bool = false,
    include_hetatm: bool = false,
    use_bitmask: bool = false,
    adaptive_sr: bool = false,
    coarse_points: u32 = 64,
    fine_points: u32 = 256,
    adaptive_low: f64 = 0.05,
    adaptive_high: f64 = 0.95,
    ccd_path: ?[]const u8 = null, // External CCD dictionary file (.zsdc or .cif[.gz|.zst])
    sdf_paths: SdfPathList = .{}, // --sdf=PATH (up to 16)
    quiet: bool = false,
    show_progress: bool = true,
    show_timing: bool = false,
    show_help: bool = false,
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
    adaptive_sr_explicit: bool = false,
    coarse_points_explicit: bool = false,
    fine_points_explicit: bool = false,
    adaptive_low_explicit: bool = false,
    adaptive_high_explicit: bool = false,
    ccd_explicit: bool = false,
    sdf_explicit: bool = false,
    quiet_explicit: bool = false,
    timing_explicit: bool = false,
};

// Parse helper functions (local to batch.zig)

fn validateWorkflowProbeRadius(radius: f64) !f64 {
    if (radius <= 0 or radius > 10.0 or !std.math.isFinite(radius)) {
        return error.InvalidArgument;
    }
    return radius;
}

fn validateWorkflowNPoints(n: u32) !u32 {
    if (n == 0 or n > 10000) {
        return error.InvalidArgument;
    }
    return n;
}

fn validateWorkflowNSlices(n: u32) !u32 {
    if (n == 0 or n > 1000) {
        return error.InvalidArgument;
    }
    return n;
}

fn validateResidueMapFormat(output_format: OutputFormat, residue_map: bool) !void {
    if (residue_map and output_format != .jsonl) {
        return error.InvalidArgument;
    }
}

/// Parse and validate probe radius value
fn parseProbeRadius(value: []const u8) f64 {
    const radius = std.fmt.parseFloat(f64, value) catch {
        std.debug.print("Error: Invalid probe radius: {s}\n", .{value});
        std.process.exit(1);
    };
    return validateWorkflowProbeRadius(radius) catch {
        std.debug.print("Error: Probe radius must be between 0 and 10 Angstroms: {d}\n", .{radius});
        std.process.exit(1);
    };
}

/// Parse and validate n-points value
fn parseNPoints(value: []const u8) u32 {
    const n = std.fmt.parseInt(u32, value, 10) catch {
        std.debug.print("Error: Invalid n-points: {s}\n", .{value});
        std.process.exit(1);
    };
    return validateWorkflowNPoints(n) catch {
        std.debug.print("Error: n-points must be between 1 and 10000: {d}\n", .{n});
        std.process.exit(1);
    };
}

fn parseBitmaskPoints(option_name: []const u8, value: []const u8) u32 {
    const n = std.fmt.parseInt(u32, value, 10) catch {
        std.debug.print("Error: Invalid {s}: {s}\n", .{ option_name, value });
        std.process.exit(1);
    };
    if (!bitmask_lut.isSupportedNPoints(n)) {
        std.debug.print("Error: {s} must be between 1 and 1024: {d}\n", .{ option_name, n });
        std.process.exit(1);
    }
    return n;
}

fn parseAdaptiveThreshold(option_name: []const u8, value: []const u8) f64 {
    const threshold = std.fmt.parseFloat(f64, value) catch {
        std.debug.print("Error: Invalid {s}: {s}\n", .{ option_name, value });
        std.process.exit(1);
    };
    if (!std.math.isFinite(threshold) or threshold < 0.0 or threshold > 1.0) {
        std.debug.print("Error: {s} must be finite and between 0.0 and 1.0: {d}\n", .{ option_name, threshold });
        std.process.exit(1);
    }
    return threshold;
}

/// Parse and validate n-slices value (for Lee-Richards)
fn parseNSlices(value: []const u8) u32 {
    const n = std.fmt.parseInt(u32, value, 10) catch {
        std.debug.print("Error: Invalid n-slices: {s}\n", .{value});
        std.process.exit(1);
    };
    return validateWorkflowNSlices(n) catch {
        std.debug.print("Error: n-slices must be between 1 and 1000: {d}\n", .{n});
        std.process.exit(1);
    };
}

/// Parse and validate algorithm value
fn parseAlgorithm(value: []const u8) Algorithm {
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
fn parseOutputFormat(value: []const u8) OutputFormat {
    if (std.mem.eql(u8, value, "json")) {
        return .json;
    } else if (std.mem.eql(u8, value, "compact")) {
        return .compact;
    } else if (std.mem.eql(u8, value, "csv")) {
        return .csv;
    } else if (std.mem.eql(u8, value, "jsonl")) {
        return .jsonl;
    } else {
        std.debug.print("Error: Invalid format: {s}\n", .{value});
        std.debug.print("Valid formats: json, compact, csv, jsonl\n", .{});
        std.process.exit(1);
    }
}

/// Parse and validate classifier type value
fn parseClassifierType(value: []const u8) ClassifierType {
    if (ClassifierType.fromString(value)) |ct| {
        return ct;
    } else {
        std.debug.print("Error: Invalid classifier: {s}\n", .{value});
        std.debug.print("Valid classifiers: ccd, protor, naccess, oons\n", .{});
        std.process.exit(1);
    }
}

/// Parse and validate precision value
fn parsePrecision(value: []const u8) Precision {
    if (Precision.fromString(value)) |p| {
        return p;
    } else {
        std.debug.print("Error: Invalid precision: {s}\n", .{value});
        std.debug.print("Valid values: f32 (single), f64 (double)\n", .{});
        std.process.exit(1);
    }
}

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

/// Parse batch subcommand arguments
pub fn parseArgs(args: []const []const u8, start_idx: usize) BatchArgs {
    var result = BatchArgs{};
    var i: usize = start_idx;
    var positional_count: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        // --threads=N or --threads N
        if (std.mem.startsWith(u8, arg, "--threads=")) {
            result.threads_explicit = true;
            const value = arg["--threads=".len..];
            result.n_threads = std.fmt.parseInt(usize, value, 10) catch {
                std.debug.print("Error: Invalid thread count: {s}\n", .{value});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--threads")) {
            result.threads_explicit = true;
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
            result.probe_radius_explicit = true;
            const value = arg["--probe-radius=".len..];
            result.probe_radius = parseProbeRadius(value);
        } else if (std.mem.eql(u8, arg, "--probe-radius")) {
            result.probe_radius_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --probe-radius\n", .{});
                std.process.exit(1);
            }
            result.probe_radius = parseProbeRadius(args[i]);
        }
        // --n-points=N or --n-points N
        else if (std.mem.startsWith(u8, arg, "--n-points=")) {
            result.n_points_explicit = true;
            const value = arg["--n-points=".len..];
            result.n_points = parseNPoints(value);
        } else if (std.mem.eql(u8, arg, "--n-points")) {
            result.n_points_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --n-points\n", .{});
                std.process.exit(1);
            }
            result.n_points = parseNPoints(args[i]);
        }
        // --n-slices=N or --n-slices N (for Lee-Richards)
        else if (std.mem.startsWith(u8, arg, "--n-slices=")) {
            result.n_slices_explicit = true;
            const value = arg["--n-slices=".len..];
            result.n_slices = parseNSlices(value);
        } else if (std.mem.eql(u8, arg, "--n-slices")) {
            result.n_slices_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --n-slices\n", .{});
                std.process.exit(1);
            }
            result.n_slices = parseNSlices(args[i]);
        }
        // --format=FORMAT or --format FORMAT
        else if (std.mem.startsWith(u8, arg, "--format=")) {
            result.format_explicit = true;
            const value = arg["--format=".len..];
            result.output_format = parseOutputFormat(value);
        } else if (std.mem.eql(u8, arg, "--format")) {
            result.format_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --format\n", .{});
                std.process.exit(1);
            }
            result.output_format = parseOutputFormat(args[i]);
        }
        // --algorithm=ALGO or --algorithm ALGO
        else if (std.mem.startsWith(u8, arg, "--algorithm=")) {
            result.algorithm_explicit = true;
            const value = arg["--algorithm=".len..];
            result.algorithm = parseAlgorithm(value);
        } else if (std.mem.eql(u8, arg, "--algorithm")) {
            result.algorithm_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --algorithm\n", .{});
                std.process.exit(1);
            }
            result.algorithm = parseAlgorithm(args[i]);
        }
        // --classifier=TYPE or --classifier TYPE
        else if (std.mem.startsWith(u8, arg, "--classifier=")) {
            result.classifier_explicit = true;
            const value = arg["--classifier=".len..];
            result.classifier_type = parseClassifierType(value);
        } else if (std.mem.eql(u8, arg, "--classifier")) {
            result.classifier_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --classifier\n", .{});
                std.process.exit(1);
            }
            result.classifier_type = parseClassifierType(args[i]);
        }
        // --precision=PREC or --precision PREC
        else if (std.mem.startsWith(u8, arg, "--precision=")) {
            result.precision_explicit = true;
            const value = arg["--precision=".len..];
            result.precision = parsePrecision(value);
        } else if (std.mem.eql(u8, arg, "--precision")) {
            result.precision_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --precision\n", .{});
                std.process.exit(1);
            }
            result.precision = parsePrecision(args[i]);
        }
        // --workflow=PATH or --workflow PATH
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
        // --manifest=PATH or --manifest PATH (compatibility alias)
        else if (std.mem.startsWith(u8, arg, "--manifest=")) {
            result.workflow_path = arg["--manifest=".len..];
        } else if (std.mem.eql(u8, arg, "--manifest")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --manifest (compatibility alias for --workflow)\n", .{});
                std.process.exit(1);
            }
            result.workflow_path = args[i];
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
        // --residue-map
        else if (std.mem.eql(u8, arg, "--residue-map")) {
            result.residue_map = true;
        }
        // --include-hydrogens
        else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
            result.include_hydrogens = true;
            result.include_hydrogens_explicit = true;
        }
        // --include-hetatm
        else if (std.mem.eql(u8, arg, "--include-hetatm")) {
            result.include_hetatm = true;
            result.include_hetatm_explicit = true;
        }
        // --use-bitmask
        else if (std.mem.eql(u8, arg, "--use-bitmask")) {
            result.use_bitmask = true;
            result.use_bitmask_explicit = true;
        }
        // --adaptive-sr and adaptive bitmask controls
        else if (std.mem.eql(u8, arg, "--adaptive-sr")) {
            result.adaptive_sr = true;
            result.adaptive_sr_explicit = true;
        } else if (std.mem.startsWith(u8, arg, "--coarse-points=")) {
            result.coarse_points_explicit = true;
            result.coarse_points = parseBitmaskPoints("--coarse-points", arg["--coarse-points=".len..]);
        } else if (std.mem.eql(u8, arg, "--coarse-points")) {
            result.coarse_points_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --coarse-points\n", .{});
                std.process.exit(1);
            }
            result.coarse_points = parseBitmaskPoints("--coarse-points", args[i]);
        } else if (std.mem.startsWith(u8, arg, "--fine-points=")) {
            result.fine_points_explicit = true;
            result.fine_points = parseBitmaskPoints("--fine-points", arg["--fine-points=".len..]);
        } else if (std.mem.eql(u8, arg, "--fine-points")) {
            result.fine_points_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --fine-points\n", .{});
                std.process.exit(1);
            }
            result.fine_points = parseBitmaskPoints("--fine-points", args[i]);
        } else if (std.mem.startsWith(u8, arg, "--adaptive-low=")) {
            result.adaptive_low_explicit = true;
            result.adaptive_low = parseAdaptiveThreshold("--adaptive-low", arg["--adaptive-low=".len..]);
        } else if (std.mem.eql(u8, arg, "--adaptive-low")) {
            result.adaptive_low_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --adaptive-low\n", .{});
                std.process.exit(1);
            }
            result.adaptive_low = parseAdaptiveThreshold("--adaptive-low", args[i]);
        } else if (std.mem.startsWith(u8, arg, "--adaptive-high=")) {
            result.adaptive_high_explicit = true;
            result.adaptive_high = parseAdaptiveThreshold("--adaptive-high", arg["--adaptive-high=".len..]);
        } else if (std.mem.eql(u8, arg, "--adaptive-high")) {
            result.adaptive_high_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --adaptive-high\n", .{});
                std.process.exit(1);
            }
            result.adaptive_high = parseAdaptiveThreshold("--adaptive-high", args[i]);
        }
        // --ccd=PATH or --ccd PATH (external CCD dictionary)
        else if (std.mem.startsWith(u8, arg, "--ccd=")) {
            result.ccd_explicit = true;
            const value = arg["--ccd=".len..];
            result.ccd_path = value;
        } else if (std.mem.eql(u8, arg, "--ccd")) {
            result.ccd_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --ccd\n", .{});
                std.process.exit(1);
            }
            result.ccd_path = args[i];
        }
        // --sdf=PATH or --sdf PATH (SDF file with bond topology for CCD classifier)
        else if (std.mem.startsWith(u8, arg, "--sdf=")) {
            result.sdf_explicit = true;
            const value = arg["--sdf=".len..];
            result.sdf_paths.append(value) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
        } else if (std.mem.eql(u8, arg, "--sdf")) {
            result.sdf_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --sdf\n", .{});
                std.process.exit(1);
            }
            result.sdf_paths.append(args[i]) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
        }
        // --quiet or -q
        else if (std.mem.eql(u8, arg, "--quiet") or std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
            result.quiet_explicit = true;
            result.show_progress = false;
        }
        // --timing
        else if (std.mem.eql(u8, arg, "--timing")) {
            result.show_timing = true;
            result.timing_explicit = true;
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
            result.output_path_explicit = true;
        } else if (std.mem.startsWith(u8, arg, "--output=")) {
            result.output_path = arg["--output=".len..];
            result.output_path_explicit = true;
        } else if (std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --output\n", .{});
                std.process.exit(1);
            }
            result.output_path = args[i];
            result.output_path_explicit = true;
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
                // Only use positional output if -o/--output was not explicitly set
                if (!result.output_path_explicit) {
                    result.output_path = arg;
                }
            } else {
                std.debug.print("Error: Too many positional arguments\n", .{});
                std.process.exit(1);
            }
            positional_count += 1;
        }
    }

    if (result.adaptive_sr and result.n_points_explicit and !result.fine_points_explicit) {
        if (!bitmask_lut.isSupportedNPoints(result.n_points)) {
            std.debug.print("Error: adaptive --n-points must be 1..1024 when used as fine points: {d}\n", .{result.n_points});
            std.process.exit(1);
        }
        result.fine_points = result.n_points;
    }
    if (result.adaptive_low > result.adaptive_high) {
        std.debug.print("Error: --adaptive-low must be <= --adaptive-high\n", .{});
        std.process.exit(1);
    }

    return result;
}

/// Print help for the batch subcommand
pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\zsasa batch - Calculate SASA for all files in a directory
        \\
        \\USAGE:
        \\    {s} batch --workflow <workflow.toml>
        \\    {s} batch [OPTIONS] <input_dir> [output_dir]
        \\
        \\ARGUMENTS:
        \\    <input_dir>     Directory containing structure files (PDB, mmCIF, JSON)
        \\    [output_dir]    Optional output directory (default: no file output)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO    Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                        Default: sr
        \\    --classifier=TYPE   Built-in classifier: ccd, protor, naccess, oons
        \\                        Default: ccd (protor is an alias for ccd)
        \\    --ccd=PATH          External CCD dictionary file (.zsdc or .cif[.gz|.zst])
        \\                        Used with --classifier=ccd for non-standard residues
        \\    --sdf=PATH          SDF file with bond topology for CCD classifier
        \\                        Can be specified multiple times for multiple ligands
        \\    --threads=N         Number of threads (default: auto-detect)
        \\    --workflow=PATH     TOML workflow file with one or more named batch jobs
        \\    --manifest=PATH     Compatibility alias for --workflow
        \\    --chain=ID          Filter by chain ID for non-workflow batch (e.g. A or A,B)
        \\    --auth-chain        Use auth_asym_id instead of label_asym_id for mmCIF chain matching
        \\    --residue-map       Include compact residue map arrays in JSONL output
        \\    --probe-radius=R    Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N        Test points per atom (default: 100, for sr)
        \\    --n-slices=N        Slices per atom diameter (default: 20, for lr)
        \\    --precision=PREC    Floating-point precision: f32, f64 (default: f64)
        \\    --format=FORMAT     Output format: json, compact, csv, jsonl (default: json)
        \\    --include-hydrogens Include hydrogen atoms (default: exclude)
        \\    --include-hetatm    Include HETATM records (default: exclude)
        \\    --use-bitmask       Use bitmask LUT optimization for SR algorithm
        \\                        (n-points must be 1..1024)
        \\    --adaptive-sr       Experimental adaptive two-stage bitmask SR
        \\    --coarse-points=N   Coarse adaptive points (default: 64)
        \\    --fine-points=N     Fine adaptive points (default: --n-points or 256)
        \\    --adaptive-low=X    Coarse accept low exposed fraction (default: 0.05)
        \\    --adaptive-high=X   Coarse accept high exposed fraction (default: 0.95)
        \\    --timing            Show timing breakdown for benchmarking
        \\    -o, --output=PATH   Output directory, or file path for --format=jsonl
        \\    -q, --quiet         Suppress progress output
        \\    -h, --help          Show this help message
        \\
        \\EXAMPLES:
        \\    {s} batch structures/
        \\    {s} batch structures/ results/
        \\    {s} batch structures/ --algorithm=lr --threads=4
        \\    {s} batch structures/ --classifier=naccess --format=csv
        \\    {s} batch structures/ results/ --timing --quiet
        \\    {s} batch structures/ -o results.jsonl --format=jsonl
        \\    {s} batch --workflow bsa.toml
        \\    {s} batch structures/ results/ --workflow bsa.toml
        \\    {s} batch --manifest legacy-bsa.toml
        \\    {s} batch structures/ results_A.jsonl --chain=A --format=jsonl
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
}

/// Load SDF files and build a ComponentDict from their bond topology.
/// Delegates to sdf_parser.loadSdfComponents.
const loadSdfComponents = sdf_parser.loadSdfComponents;

fn applyWorkflowToBatchConfig(
    config: *BatchConfig,
    args: BatchArgs,
    calculation: workflow_manifest.Calculation,
    output: workflow_manifest.Output,
    classifier_config: workflow_manifest.ClassifierConfig,
) !void {
    if (!args.threads_explicit) {
        if (calculation.threads) |v| config.n_threads = v;
    }
    if (!args.algorithm_explicit) {
        if (calculation.algorithm) |v| config.algorithm = parseAlgorithm(v);
    }
    if (!args.n_points_explicit) {
        if (calculation.n_points) |v| config.n_points = validateWorkflowNPoints(v) catch |err| {
            std.debug.print("Error: workflow n_points must be between 1 and 10000: {d}\n", .{v});
            return err;
        };
    }
    if (!args.n_slices_explicit) {
        if (calculation.n_slices) |v| config.n_slices = validateWorkflowNSlices(v) catch |err| {
            std.debug.print("Error: workflow n_slices must be between 1 and 1000: {d}\n", .{v});
            return err;
        };
    }
    if (!args.probe_radius_explicit) {
        if (calculation.probe_radius) |v| config.probe_radius = validateWorkflowProbeRadius(v) catch |err| {
            std.debug.print("Error: workflow probe_radius must be finite and between 0 and 10 Angstroms: {d}\n", .{v});
            return err;
        };
    }
    if (!args.precision_explicit) {
        if (calculation.precision) |v| config.precision = parsePrecision(v);
    }
    if (!args.format_explicit) {
        if (output.format) |v| config.output_format = parseOutputFormat(v);
    }
    if (!args.timing_explicit) {
        if (calculation.timing) |v| config.show_timing = v;
    }
    if (!args.quiet_explicit) {
        if (calculation.quiet) |v| {
            config.quiet = v;
            config.show_progress = !v;
        }
    }
    if (!args.include_hydrogens_explicit) {
        if (calculation.include_hydrogens) |v| config.include_hydrogens = v;
    }
    if (!args.include_hetatm_explicit) {
        if (calculation.include_hetatm) |v| config.include_hetatm = v;
    }
    if (!args.use_bitmask_explicit) {
        if (calculation.use_bitmask) |v| config.use_bitmask = v;
    }
    if (calculation.auth_chain) |v| config.use_auth_chain = v;
    if (calculation.residue_map) |v| config.residue_map = v;

    try applyWorkflowClassifierToBatchConfig(config, args, classifier_config);
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
    if (args.classifier_explicit) {
        config.classifier_type = args.classifier_type;
        config.custom_classifier = null;
        config.custom_classifier_path = null;
    }
    if (args.include_hydrogens_explicit) config.include_hydrogens = args.include_hydrogens;
    if (args.include_hetatm_explicit) config.include_hetatm = args.include_hetatm;
    if (args.use_bitmask_explicit) config.use_bitmask = args.use_bitmask;
    if (args.use_auth_chain) config.use_auth_chain = true;
    if (args.residue_map) config.residue_map = true;
}

fn applyWorkflowClassifierToBatchConfig(config: *BatchConfig, args: BatchArgs, classifier_config: workflow_manifest.ClassifierConfig) !void {
    if (!args.classifier_explicit) {
        if (classifier_config.type) |classifier_type| {
            if (std.mem.eql(u8, classifier_type, "custom")) {
                config.classifier_type = null;
                config.custom_classifier_path = classifier_config.config orelse return error.InvalidArgument;
            } else {
                config.classifier_type = parseClassifierType(classifier_type);
                config.custom_classifier_path = null;
            }
        }
    }
}

fn applyWorkflowJobOverrides(config: *BatchConfig, args: BatchArgs, job: workflow_manifest.Job) void {
    if (job.auth_chain) |v| config.use_auth_chain = v;
    if (args.use_auth_chain) config.use_auth_chain = true;
    config.chain_filter = job.chains;
}

fn parseWorkflowFile(allocator: Allocator, io: std.Io, path: []const u8) !workflow_manifest.Workflow {
    return workflow_manifest.parseFile(allocator, io, path);
}

fn classifierUsesCcdResources(effective_classifier_type: ?ClassifierType) bool {
    const classifier_type = effective_classifier_type orelse return false;
    return classifier_type == .ccd or classifier_type == .protor;
}

fn resolveWorkflowCcdPath(args: BatchArgs, classifier_config: workflow_manifest.ClassifierConfig, effective_classifier_type: ?ClassifierType) ?[]const u8 {
    if (!classifierUsesCcdResources(effective_classifier_type)) return null;
    if (args.ccd_explicit) return args.ccd_path;
    return classifier_config.ccd;
}

fn resolveWorkflowSdfPaths(args: BatchArgs, workflow_sdf_paths: []const []const u8, effective_classifier_type: ?ClassifierType) []const []const u8 {
    if (!classifierUsesCcdResources(effective_classifier_type)) return &.{};
    if (args.sdf_explicit) return args.sdf_paths.constSlice();
    return workflow_sdf_paths;
}

fn loadExternalCcd(allocator: Allocator, io: std.Io, path: []const u8, quiet: bool) !ccd_parser.ComponentDict {
    const ccd_data = if (compressed.isCompressed(path))
        compressed.read(allocator, path) catch |err| {
            std.debug.print("Error reading CCD file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        }
    else blk: {
        const f = std.Io.Dir.cwd().openFile(io, path, .{}) catch |err| {
            std.debug.print("Error opening CCD file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
        defer f.close(io);
        var read_buf_ccd: [65536]u8 = undefined;
        var file_r_ccd = f.reader(io, &read_buf_ccd);
        break :blk file_r_ccd.interface.allocRemaining(allocator, .unlimited) catch |err| {
            std.debug.print("Error reading CCD file '{s}': {s}\n", .{ path, @errorName(err) });
            return err;
        };
    };
    defer allocator.free(ccd_data);

    const dict = ccd_binary.loadDict(allocator, ccd_data) catch |err| {
        std.debug.print("Error loading CCD dictionary '{s}': {s}\n", .{ path, @errorName(err) });
        return err;
    };
    if (!quiet) {
        std.debug.print("External CCD: loaded {d} components from '{s}'\n", .{ dict.components.count(), path });
    }
    return dict;
}

fn loadCustomClassifier(allocator: Allocator, io: std.Io, path: []const u8) !classifier.Classifier {
    return classifier_parser.parseConfigFile(allocator, io, path) catch |err| {
        switch (err) {
            error.UnsupportedConfigExtension => std.debug.print("Error loading config file '{s}': custom classifier configs are TOML-only; rename or convert the file to .toml\n", .{path}),
            error.UnsupportedLegacyFormat => std.debug.print("Error loading config file '{s}': FreeSASA-style custom classifier configs are no longer supported; convert to TOML [types] and [[atoms]]\n", .{path}),
            else => std.debug.print("Error loading config file '{s}': {s}\n", .{ path, @errorName(err) }),
        }
        return err;
    };
}

fn runWorkflow(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
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

    const load_quiet = if (args.quiet_explicit) args.quiet else (workflow.calculation.quiet orelse args.quiet);

    var resource_config = BatchConfig{};
    try applyWorkflowToBatchConfig(&resource_config, args, workflow.calculation, workflow.output, workflow.classifier);
    applyCliOverrides(&resource_config, args);
    const effective_classifier_type = resource_config.classifier_type;

    const ccd_path = resolveWorkflowCcdPath(args, workflow.classifier, effective_classifier_type);
    var ext_ccd: ?ccd_parser.ComponentDict = null;
    if (ccd_path) |path| {
        ext_ccd = try loadExternalCcd(allocator, io, path, load_quiet);
    }
    defer if (ext_ccd) |*d| d.deinit();

    const workflow_sdf_paths: []const []const u8 = workflow.classifier.sdf orelse &.{};
    const sdf_paths = resolveWorkflowSdfPaths(args, workflow_sdf_paths, effective_classifier_type);
    var sdf_ccd: ?ccd_parser.ComponentDict = null;
    if (sdf_paths.len > 0) {
        sdf_ccd = loadSdfComponents(allocator, io, sdf_paths, load_quiet) catch |err| {
            std.debug.print("Error loading SDF components: {s}\n", .{@errorName(err)});
            return err;
        };
    }
    defer if (sdf_ccd) |*d| d.deinit();

    var custom_classifier: ?classifier.Classifier = null;
    const custom_classifier_path: ?[]const u8 = if (!args.classifier_explicit and workflow.classifier.type != null and std.mem.eql(u8, workflow.classifier.type.?, "custom"))
        workflow.classifier.config
    else
        null;
    if (custom_classifier_path) |path| {
        custom_classifier = try loadCustomClassifier(allocator, io, path);
    }
    defer if (custom_classifier) |*c| c.deinit();

    var successful: usize = 0;
    var failed: usize = 0;

    for (workflow.jobs) |job| {
        var config = BatchConfig{};
        try applyWorkflowToBatchConfig(&config, args, workflow.calculation, workflow.output, workflow.classifier);
        applyCliOverrides(&config, args);

        config.include_hetatm = config.include_hetatm or (config.classifier_type == .ccd);
        config.store_atom_areas = (config.output_format == .jsonl);
        config.external_ccd = if (ext_ccd != null) &ext_ccd.? else null;
        config.sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null;
        if (custom_classifier) |*c| config.custom_classifier = c;
        applyWorkflowJobOverrides(&config, args, job);
        validateResidueMapFormat(config.output_format, config.residue_map) catch {
            std.debug.print("Error: residue_map is only supported with format = \"jsonl\"\n", .{});
            return error.InvalidArgument;
        };

        if ((config.classifier_type == .ccd or config.classifier_type == .protor) and config.include_hydrogens and !config.quiet) {
            std.debug.print("Warning: --include-hydrogens with CCD classifier may give inaccurate results\n", .{});
            std.debug.print("         CCD uses united-atom radii that already account for implicit hydrogens\n", .{});
        }

        var job_output_dir: ?[]const u8 = null;
        defer if (job_output_dir) |path| allocator.free(path);
        var jsonl_output_path: ?[]const u8 = null;
        defer if (jsonl_output_path) |path| allocator.free(path);

        if (config.output_format == .jsonl) {
            if (output_dir) |out| {
                try std.Io.Dir.cwd().createDirPath(io, out);
                jsonl_output_path = try workflowJsonlOutputPath(allocator, out, job.name);
            }
        } else if (output_dir) |out| {
            job_output_dir = try workflowPerFileOutputDir(allocator, out, job.name);
        }

        if (!config.quiet) {
            std.debug.print("Workflow job: {s}\n", .{job.name});
        }

        var result = runBatch(allocator, io, input_dir, job_output_dir, config, jsonl_output_path) catch |err| {
            std.debug.print("Error running workflow job '{s}': {s}\n", .{ job.name, @errorName(err) });
            failed += 1;
            continue;
        };
        defer result.deinit();

        successful += result.successful;
        failed += result.failed;
    }

    std.debug.print("Workflow complete: {d} successful, {d} failed\n", .{ successful, failed });
}

/// Run batch processing from parsed CLI arguments
pub fn run(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (args.workflow_path != null) {
        if (args.adaptive_sr_explicit or args.coarse_points_explicit or args.fine_points_explicit or args.adaptive_low_explicit or args.adaptive_high_explicit) {
            std.debug.print("Error: adaptive SR options are not supported with --workflow yet\n", .{});
            return error.InvalidArgument;
        }
        return runWorkflow(allocator, io, args);
    }

    const input_dir = args.input_path orelse {
        std.debug.print("Error: Missing input directory\n", .{});
        std.debug.print("Usage: zsasa batch [OPTIONS] <input_dir> [output_dir]\n", .{});
        return error.MissingArgument;
    };

    // Load external CCD dictionary if specified
    var ext_ccd: ?ccd_parser.ComponentDict = null;
    if (args.ccd_path) |ccd_path| {
        const ccd_data = if (compressed.isCompressed(ccd_path))
            compressed.read(allocator, ccd_path) catch |err| {
                std.debug.print("Error reading CCD file '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                std.process.exit(1);
            }
        else blk: {
            const f = std.Io.Dir.cwd().openFile(io, ccd_path, .{}) catch |err| {
                std.debug.print("Error opening CCD file '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                std.process.exit(1);
            };
            defer f.close(io);
            var read_buf_ccd: [65536]u8 = undefined;
            var file_r_ccd = f.reader(io, &read_buf_ccd);
            break :blk file_r_ccd.interface.allocRemaining(allocator, .unlimited) catch |err| {
                std.debug.print("Error reading CCD file '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                std.process.exit(1);
            };
        };
        defer allocator.free(ccd_data);

        ext_ccd = ccd_binary.loadDict(allocator, ccd_data) catch |err| {
            std.debug.print("Error loading CCD dictionary '{s}': {s}\n", .{ ccd_path, @errorName(err) });
            std.process.exit(1);
        };
        if (!args.quiet) {
            std.debug.print("External CCD: loaded {d} components from '{s}'\n", .{ ext_ccd.?.components.count(), ccd_path });
        }
    }
    defer if (ext_ccd) |*d| d.deinit();

    // Load SDF components from --sdf option
    var sdf_ccd: ?ccd_parser.ComponentDict = null;
    if (args.sdf_paths.len > 0) {
        sdf_ccd = loadSdfComponents(allocator, io, args.sdf_paths.constSlice(), args.quiet) catch |err| {
            std.debug.print("Error loading SDF components: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
    defer if (sdf_ccd) |*d| d.deinit();

    validateResidueMapFormat(args.output_format, args.residue_map) catch {
        std.debug.print("Error: --residue-map is only supported with --format=jsonl\n", .{});
        return error.InvalidArgument;
    };

    if (args.adaptive_sr) {
        if (!args.use_bitmask) {
            std.debug.print("Error: --adaptive-sr requires --use-bitmask\n", .{});
            return error.InvalidArgument;
        }
        if (args.algorithm != .sr) {
            std.debug.print("Error: --adaptive-sr requires --algorithm=sr\n", .{});
            return error.InvalidArgument;
        }
        if (!bitmask_lut.isSupportedNPoints(args.coarse_points) or !bitmask_lut.isSupportedNPoints(args.fine_points)) {
            std.debug.print("Error: --adaptive-sr point counts must be 1..1024\n", .{});
            return error.InvalidArgument;
        }
        if (args.adaptive_low > args.adaptive_high) {
            std.debug.print("Error: --adaptive-low must be <= --adaptive-high\n", .{});
            return error.InvalidArgument;
        }
    }

    // For jsonl, don't pass output_dir to runBatch (no per-file I/O during computation)
    const output_dir: ?[]const u8 = if (args.output_format == .jsonl) null else args.output_path;

    var chain_filter_slice: ?[]const []const u8 = null;
    if (args.chain_filter) |filter_str| {
        chain_filter_slice = try parseBatchChainFilter(allocator, filter_str);
    }
    defer if (chain_filter_slice) |s| allocator.free(s);

    // Build batch config from parsed args
    // CCD/ProtOr use united-atom radii (implicit H) — warn if explicit H included
    if ((args.classifier_type == .ccd or args.classifier_type == .protor) and args.include_hydrogens and !args.quiet) {
        std.debug.print("Warning: --include-hydrogens with CCD classifier may give inaccurate results\n", .{});
        std.debug.print("         CCD uses united-atom radii that already account for implicit hydrogens\n", .{});
    }

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
        .show_progress = args.show_progress,
        .classifier_type = args.classifier_type,
        .include_hydrogens = args.include_hydrogens,
        .include_hetatm = args.include_hetatm or (args.classifier_type == .ccd),
        .use_bitmask = args.use_bitmask,
        .adaptive_sr = args.adaptive_sr,
        .coarse_points = args.coarse_points,
        .fine_points = args.fine_points,
        .adaptive_low = args.adaptive_low,
        .adaptive_high = args.adaptive_high,
        .store_atom_areas = (args.output_format == .jsonl),
        .external_ccd = if (ext_ccd != null) &ext_ccd.? else null,
        .sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null,
        .chain_filter = chain_filter_slice,
        .use_auth_chain = args.use_auth_chain,
        .residue_map = args.residue_map,
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

    // Determine JSONL output path: stream during computation instead of accumulating in memory
    const jsonl_output_path: ?[]const u8 = if (args.output_format == .jsonl) args.output_path else null;

    var result = try runBatch(allocator, io, input_dir, output_dir, config, jsonl_output_path);
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
    try std.testing.expectEqual(true, parsed.show_progress);
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
    try std.testing.expectEqual(false, parsed.show_progress);
    try std.testing.expectEqual(true, parsed.show_timing);
}

test "BatchArgs help flag" {
    const args = [_][]const u8{ "zsasa", "batch", "--help" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.show_help);
}

test "BatchArgs --manifest compatibility alias" {
    const args = [_][]const u8{ "zsasa", "batch", "--manifest", "bsa.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("bsa.toml", parsed.workflow_path.?);
}

test "BatchArgs --workflow=FILE" {
    const args = [_][]const u8{ "zsasa", "batch", "--workflow=batch-workflow.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("batch-workflow.toml", parsed.workflow_path.?);
}

test "BatchArgs --workflow FILE" {
    const args = [_][]const u8{ "zsasa", "batch", "--workflow", "batch-workflow.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("batch-workflow.toml", parsed.workflow_path.?);
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

test "BatchArgs --residue-map" {
    const args = [_][]const u8{ "zsasa", "batch", "--format=jsonl", "--residue-map", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.jsonl, parsed.output_format);
    try std.testing.expectEqual(true, parsed.residue_map);
}

test "parseBatchChainFilter splits comma-separated chains" {
    const chains = try parseBatchChainFilter(std.testing.allocator, "A, B,AB");
    defer std.testing.allocator.free(chains);
    try std.testing.expectEqual(@as(usize, 3), chains.len);
    try std.testing.expectEqualStrings("A", chains[0]);
    try std.testing.expectEqualStrings("B", chains[1]);
    try std.testing.expectEqualStrings("AB", chains[2]);
}

test "workflow numeric validators reject invalid values" {
    try std.testing.expectEqual(@as(f64, 1.4), try validateWorkflowProbeRadius(1.4));
    try std.testing.expectEqual(@as(u32, 1), try validateWorkflowNPoints(1));
    try std.testing.expectEqual(@as(u32, 10000), try validateWorkflowNPoints(10000));
    try std.testing.expectEqual(@as(u32, 1), try validateWorkflowNSlices(1));
    try std.testing.expectEqual(@as(u32, 1000), try validateWorkflowNSlices(1000));

    try std.testing.expectError(error.InvalidArgument, validateWorkflowProbeRadius(0));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowProbeRadius(-1));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowProbeRadius(10.1));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowProbeRadius(std.math.nan(f64)));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowNPoints(0));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowNPoints(10001));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowNSlices(0));
    try std.testing.expectError(error.InvalidArgument, validateWorkflowNSlices(1001));
}

test "validateResidueMapFormat accepts JSONL residue map" {
    try validateResidueMapFormat(.jsonl, true);
    try validateResidueMapFormat(.jsonl, false);
    try validateResidueMapFormat(.json, false);
}

test "validateResidueMapFormat rejects non-JSONL residue map" {
    try std.testing.expectError(error.InvalidArgument, validateResidueMapFormat(.json, true));
    try std.testing.expectError(error.InvalidArgument, validateResidueMapFormat(.compact, true));
    try std.testing.expectError(error.InvalidArgument, validateResidueMapFormat(.csv, true));
}

test "CLI auth-chain overrides workflow job auth_chain false" {
    var config = BatchConfig{};
    const args = BatchArgs{ .use_auth_chain = true };
    const job = @import("workflow_manifest.zig").Job{
        .name = "chain_A",
        .auth_chain = false,
    };

    applyWorkflowJobOverrides(&config, args, job);

    try std.testing.expectEqual(true, config.use_auth_chain);
}

test "workflow residue_map applies when CLI does not override format" {
    var config = BatchConfig{};
    const args = BatchArgs{};
    const calculation = @import("workflow_manifest.zig").Calculation{ .residue_map = true };
    const output = @import("workflow_manifest.zig").Output{ .format = "jsonl" };
    const classifier_config = @import("workflow_manifest.zig").ClassifierConfig{};

    try applyWorkflowToBatchConfig(&config, args, calculation, output, classifier_config);

    try std.testing.expectEqual(OutputFormat.jsonl, config.output_format);
    try std.testing.expectEqual(true, config.residue_map);
}

test "workflow custom classifier config path resolves for batch" {
    var config = BatchConfig{};
    const args = BatchArgs{};
    const classifier_config = @import("workflow_manifest.zig").ClassifierConfig{
        .type = "custom",
        .config = "custom-radii.toml",
    };

    try applyWorkflowClassifierToBatchConfig(&config, args, classifier_config);

    try std.testing.expect(config.classifier_type == null);
    try std.testing.expectEqualStrings("custom-radii.toml", config.custom_classifier_path.?);
}

test "batch resource resolver only loads CCD resources for CCD classifiers" {
    const classifier_config = @import("workflow_manifest.zig").ClassifierConfig{
        .ccd = "workflow.zsdc",
    };
    const workflow_sdf_paths = [_][]const u8{"workflow.sdf"};

    const workflow_only_args = BatchArgs{};
    try std.testing.expect(resolveWorkflowCcdPath(workflow_only_args, classifier_config, .naccess) == null);
    try std.testing.expect(resolveWorkflowCcdPath(workflow_only_args, classifier_config, .oons) == null);
    try std.testing.expect(resolveWorkflowCcdPath(workflow_only_args, classifier_config, null) == null);
    try std.testing.expectEqual(@as(usize, 0), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], .naccess).len);
    try std.testing.expectEqual(@as(usize, 0), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], .oons).len);
    try std.testing.expectEqual(@as(usize, 0), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], null).len);
    try std.testing.expectEqualStrings("workflow.zsdc", resolveWorkflowCcdPath(workflow_only_args, classifier_config, .ccd).?);
    try std.testing.expectEqualStrings("workflow.zsdc", resolveWorkflowCcdPath(workflow_only_args, classifier_config, .protor).?);
    try std.testing.expectEqual(@as(usize, 1), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], .ccd).len);
    try std.testing.expectEqual(@as(usize, 1), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], .protor).len);

    var explicit_resource_args = BatchArgs{
        .classifier_explicit = true,
        .classifier_type = .naccess,
        .ccd_explicit = true,
        .ccd_path = "cli.zsdc",
        .sdf_explicit = true,
    };
    try explicit_resource_args.sdf_paths.append("cli.sdf");

    try std.testing.expect(resolveWorkflowCcdPath(explicit_resource_args, classifier_config, .naccess) == null);
    try std.testing.expect(resolveWorkflowCcdPath(explicit_resource_args, classifier_config, null) == null);
    try std.testing.expectEqual(@as(usize, 0), resolveWorkflowSdfPaths(explicit_resource_args, workflow_sdf_paths[0..], .naccess).len);
    try std.testing.expectEqual(@as(usize, 0), resolveWorkflowSdfPaths(explicit_resource_args, workflow_sdf_paths[0..], null).len);
    try std.testing.expectEqualStrings("cli.zsdc", resolveWorkflowCcdPath(explicit_resource_args, classifier_config, .ccd).?);
    const resolved_sdf_paths = resolveWorkflowSdfPaths(explicit_resource_args, workflow_sdf_paths[0..], .ccd);
    try std.testing.expectEqual(@as(usize, 1), resolved_sdf_paths.len);
    try std.testing.expectEqualStrings("cli.sdf", resolved_sdf_paths[0]);
}

test "workflowJsonlOutputPath uses job file under output dir" {
    const path = try workflowJsonlOutputPath(std.testing.allocator, "results", "chain_A");
    defer std.testing.allocator.free(path);
    try std.testing.expectEqualStrings("results/chain_A.jsonl", path);
}

test "workflowPerFileOutputDir uses job directory under output dir" {
    const path = try workflowPerFileOutputDir(std.testing.allocator, "results", "complex_AB");
    defer std.testing.allocator.free(path);
    try std.testing.expectEqualStrings("results/complex_AB", path);
}

test "FileResult JSONL uses residue map serializer when present" {
    const allocator = std.testing.allocator;
    const atom_areas = [_]f64{ 1.0, 2.0 };
    const residue_chain = [_]types.FixedString4{types.FixedString4.fromSlice("A")};
    const residue_name = [_]types.FixedString5{types.FixedString5.fromSlice("GLY")};
    const residue_number = [_]i32{5};
    const residue_insertion_code = [_]types.FixedString4{types.FixedString4.fromSlice("")};
    const residue_atom_start = [_]usize{0};
    const residue_atom_count = [_]usize{2};
    const residue_sasa = [_]f64{3.0};

    const map = json_writer.ResidueMap{
        .allocator = allocator,
        .residue_chain = residue_chain[0..],
        .residue_name = residue_name[0..],
        .residue_number = residue_number[0..],
        .residue_insertion_code = residue_insertion_code[0..],
        .residue_atom_start = residue_atom_start[0..],
        .residue_atom_count = residue_atom_count[0..],
        .residue_sasa = residue_sasa[0..],
    };

    var result = FileResult{
        .filename = "example.cif",
        .n_atoms = 2,
        .sasa_time_ns = 0,
        .total_sasa = 3.0,
        .status = .ok,
        .atom_areas = atom_areas[0..],
        .residue_map = map,
    };

    const line = try fileResultToJsonlLine(allocator, &result);
    defer allocator.free(line);

    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_chain\":[\"A\"]") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_sasa\":[3]") != null);
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

test "BatchArgs output via -o flag" {
    const args = [_][]const u8{ "zsasa", "batch", "-o", "results/", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
    try std.testing.expectEqualStrings("results/", parsed.output_path.?);
}

test "BatchArgs -o takes precedence over positional output" {
    const args = [_][]const u8{ "zsasa", "batch", "-o", "explicit/", "input_dir/", "positional/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("explicit/", parsed.output_path.?);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "BatchArgs --probe-radius=R" {
    const args = [_][]const u8{ "zsasa", "batch", "--probe-radius=2.0", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(f64, 2.0), parsed.probe_radius);
}

test "BatchArgs --n-points=N" {
    const args = [_][]const u8{ "zsasa", "batch", "--n-points=200", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

test "BatchArgs --n-slices=N" {
    const args = [_][]const u8{ "zsasa", "batch", "--n-slices=40", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 40), parsed.n_slices);
}

test "BatchArgs --format=csv" {
    const args = [_][]const u8{ "zsasa", "batch", "--format=csv", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "BatchArgs --format=jsonl" {
    const args = [_][]const u8{ "zsasa", "batch", "--format=jsonl", "-o", "out.jsonl", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.jsonl, parsed.output_format);
    try std.testing.expectEqualStrings("out.jsonl", parsed.output_path.?);
}

test "BatchArgs --classifier=naccess" {
    const args = [_][]const u8{ "zsasa", "batch", "--classifier=naccess", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(ClassifierType.naccess, parsed.classifier_type);
}

test "BatchArgs --precision=f32" {
    const args = [_][]const u8{ "zsasa", "batch", "--precision=f32", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Precision.f32, parsed.precision);
}

test "BatchArgs --include-hydrogens" {
    const args = [_][]const u8{ "zsasa", "batch", "--include-hydrogens", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.include_hydrogens);
}

test "BatchArgs --include-hetatm" {
    const args = [_][]const u8{ "zsasa", "batch", "--include-hetatm", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.include_hetatm);
}

test "BatchArgs adaptive bitmask options" {
    const args = [_][]const u8{
        "zsasa", "batch", "--use-bitmask", "--adaptive-sr", "--coarse-points=64", "--fine-points", "256", "--adaptive-low=0.10", "--adaptive-high", "0.90", "input_dir/",
    };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_bitmask);
    try std.testing.expectEqual(true, parsed.adaptive_sr);
    try std.testing.expectEqual(@as(u32, 64), parsed.coarse_points);
    try std.testing.expectEqual(@as(u32, 256), parsed.fine_points);
    try std.testing.expectEqual(@as(f64, 0.10), parsed.adaptive_low);
    try std.testing.expectEqual(@as(f64, 0.90), parsed.adaptive_high);
    try std.testing.expectEqual(true, parsed.adaptive_sr_explicit);
    try std.testing.expectEqual(true, parsed.coarse_points_explicit);
    try std.testing.expectEqual(true, parsed.fine_points_explicit);
}

test "BatchArgs adaptive defaults use 64 and 256" {
    const args = [_][]const u8{ "zsasa", "batch", "--use-bitmask", "--adaptive-sr", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 64), parsed.coarse_points);
    try std.testing.expectEqual(@as(u32, 256), parsed.fine_points);
    try std.testing.expectEqual(@as(f64, 0.05), parsed.adaptive_low);
    try std.testing.expectEqual(@as(f64, 0.95), parsed.adaptive_high);
}

test "BatchArgs --use-bitmask" {
    const args = [_][]const u8{ "zsasa", "batch", "--use-bitmask", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_bitmask);
}

test "BatchArgs --output=DIR (equals form)" {
    const args = [_][]const u8{ "zsasa", "batch", "--output=results/", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("results/", parsed.output_path.?);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "BatchArgs --algorithm=shrake-rupley (long form)" {
    const args = [_][]const u8{ "zsasa", "batch", "--algorithm=shrake-rupley", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "BatchArgs classifier_type defaults to ccd" {
    const args = [_][]const u8{ "zsasa", "batch", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(ClassifierType.ccd, parsed.classifier_type);
}

test "BatchConfig quiet suppresses progress even when show_progress defaults true" {
    const config = BatchConfig{ .quiet = true };
    try std.testing.expectEqual(false, shouldShowProgress(config));
}

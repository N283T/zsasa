const std = @import("std");
const builtin = @import("builtin");
const workflow_manifest = @import("workflow_manifest.zig");
const chain_map = @import("chain_map.zig");
const types = @import("types.zig");
const format_detect = @import("format_detect.zig");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const bcif_parser = @import("bcif_parser.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const af_model_parser = @import("af_model_parser.zig");
const pdb_parser = @import("pdb_parser.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const lee_richards = @import("lee_richards.zig");
const classifier = @import("classifier.zig");
const classifier_parser = @import("classifier_parser.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const sdf_parser = @import("sdf_parser.zig");
const compressed = @import("compressed.zig");
const input_io = @import("input_io.zig");
const element_module = @import("element.zig");
const test_support = @import("test_support.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const SasaResultGen = types.SasaResultGen;
const Config = types.Config;
const ConfigGen = types.ConfigGen;
const Precision = types.Precision;
const ClassifierType = classifier.ClassifierType;
const OutputFormat = json_writer.OutputFormat;
pub const InputIoMode = input_io.InputIoMode;
const LeeRichardsConfig = lee_richards.LeeRichardsConfig;
const LeeRichardsConfigGen = lee_richards.LeeRichardsConfigGen;
const TrigMode = lee_richards.TrigMode;

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

const JsonlMetadataMode = enum {
    none,
    sidecar,
};

/// Configuration for batch processing
pub const BatchConfig = struct {
    n_threads: usize = 0, // 0 = auto-detect
    algorithm: Algorithm = .sr,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    lr_trig: TrigMode = .exact, // Lee-Richards arc angles: exact or approximate
    probe_radius: f64 = 1.4,
    output_format: OutputFormat = .json,
    show_timing: bool = false,
    profile_stages: bool = false,
    input_io: InputIoMode = .auto,
    quiet: bool = false,
    show_progress: bool = true,
    precision: Precision = .f64, // f32 or f64
    classifier_type: ?ClassifierType = .ccd, // Default: ccd (ProtOr-compatible with CCD extension)
    include_hydrogens: bool = false, // Include hydrogen atoms (default: exclude)
    include_hetatm: bool = false, // Include HETATM records (default: exclude)
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 1..1024)
    bitmask_correction: bool = false, // Experimental exposed-fraction correction for bitmask SR
    bitmask_correction_coeff: f64 = shrake_rupley_bitmask.default_bitmask_correction_coeff,
    adaptive_sr: bool = false, // Experimental two-stage bitmask SR for batch mode
    coarse_points: u32 = 64,
    fine_points: u32 = 256,
    adaptive_low: f64 = 0.10,
    adaptive_high: f64 = 0.90,
    store_atom_areas: bool = false, // When true, keep atom_areas in each result (for JSONL rows)
    external_ccd: ?*const ccd_parser.ComponentDict = null, // External CCD dictionary
    sdf_ccd: ?*const ccd_parser.ComponentDict = null, // SDF bond topology dictionary
    custom_classifier: ?*const classifier.Classifier = null,
    custom_classifier_path: ?[]const u8 = null,
    chain_filter: ?[]const []const u8 = null,
    chain_map: ?*const chain_map.ChainMap = null,
    use_auth_chain: bool = false,
    alt_loc_mode: mmcif_parser.AltLocMode = .auto,
    alt_loc_id: u8 = 'A',
    af_model_fast: bool = false,
    residue_map: bool = false,
    jsonl_decimals: ?u8 = null,
    jsonl_include_atom_areas: bool = true,
    jsonl_include_atom_identity: bool = false,
    jsonl_include_total_area: bool = true,
    jsonl_metadata: JsonlMetadataMode = .none,
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
        .freesasa, .rsa => unreachable,
    };
}

fn workflowJsonlOutputPath(allocator: Allocator, output_dir: []const u8, job_name: []const u8) ![]const u8 {
    const filename = try std.fmt.allocPrint(allocator, "{s}.jsonl", .{job_name});
    defer allocator.free(filename);
    return std.fs.path.join(allocator, &.{ output_dir, filename });
}

fn workflowAnalysisJsonlOutputPath(allocator: Allocator, output_dir: []const u8, analysis_name: []const u8) ![]const u8 {
    const filename = try std.fmt.allocPrint(allocator, "{s}.jsonl", .{analysis_name});
    defer allocator.free(filename);
    return std.fs.path.join(allocator, &.{ output_dir, filename });
}

fn workflowJsonlMetadataPath(allocator: Allocator, output_dir: []const u8, name: []const u8) ![]const u8 {
    const filename = try std.fmt.allocPrint(allocator, "{s}.meta.json", .{name});
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

/// Stem of an SDF/MOL input file name: the name without its compression and
/// format extensions ("lig.v2.sdf.gz" -> "lig.v2").
fn sdfFileStem(filename: []const u8) []const u8 {
    var base = filename;
    if (compressed.isGzip(base)) {
        base = base[0 .. base.len - 3];
    } else if (compressed.isZstd(base)) {
        base = base[0 .. base.len - 4];
    }
    return if (std.mem.findScalarLast(u8, base, '.')) |dot_idx|
        base[0..dot_idx]
    else
        base;
}

/// Whether `c` cannot be used in a file name component. Path separators and
/// control characters are replaced on every platform so that output names do
/// not depend on where the batch runs; the other characters Windows reserves
/// are only replaced there, where such a name cannot be created anyway.
fn isUnsafeFileNameByte(c: u8, os_tag: std.Target.Os.Tag) bool {
    return switch (c) {
        '/', '\\' => true,
        '<', '>', ':', '"', '|', '?', '*' => os_tag == .windows,
        else => std.ascii.isControl(c),
    };
}

/// Per-file output name of an SDF molecule: its display name with every byte
/// that cannot be part of a file name replaced by '_', followed by `ext`.
///
/// The extension is appended, never substituted for a "previous extension",
/// so dots in the file stem or the molecule title are kept. The result holds
/// no path separator and always ends in `ext`, so it is never "", "." or "..":
/// joined to the output directory it cannot name anything outside of it.
fn sdfMoleculeOutputName(allocator: Allocator, display_name: []const u8, ext: []const u8) ![]u8 {
    const name = try std.mem.concat(allocator, u8, &.{ display_name, ext });
    for (name[0..display_name.len]) |*c| {
        if (isUnsafeFileNameByte(c.*, builtin.os.tag)) c.* = '_';
    }
    return name;
}

/// What a per-file output is named after.
const OutputNameSource = enum {
    /// An input file name; its extension is replaced ("1ubq.cif.gz" -> "1ubq.json").
    input_file,
    /// The display name of an SDF molecule; see `sdfMoleculeOutputName`.
    sdf_molecule,
};

/// Name of the per-file output written for `name`. The collision check and
/// the writer both go through this function, so they agree on every name.
fn perFileOutputName(allocator: Allocator, source: OutputNameSource, name: []const u8, ext: []const u8) ![]const u8 {
    return switch (source) {
        .input_file => replaceExtension(allocator, name, ext),
        .sdf_molecule => sdfMoleculeOutputName(allocator, name, ext),
    };
}

/// Build the display names of the molecules of one SDF/MOL file, in file
/// order. A display name is the `filename` of the molecule in batch results
/// and JSONL rows, and its per-file output is named after it.
///
/// - A molecule is named "stem_title", or "stem_N" (N = 1-based position in
///   the file) when its title is blank.
/// - When that name is shared by several molecules of the file, each of them
///   gets its position appended: "stem_title_N".
/// - While the result is still the name of another molecule of the file,
///   "_N" is appended again: titles "x", "x", "x_2" give "stem_x_1",
///   "stem_x_2_2" and "stem_x_2".
///
/// Names are compared as the output names they produce and without regard to
/// ASCII case, so the per-file outputs of one file never overwrite each other,
/// also on a case-insensitive filesystem. A molecule whose name is unique in
/// this sense always keeps the plain "stem_title" or "stem_N". Titles are kept
/// verbatim: display names are not file names, see `sdfMoleculeOutputName`.
///
/// Caller owns the returned slice and every name in it.
fn sdfMoleculeDisplayNames(
    allocator: Allocator,
    filename: []const u8,
    molecules: []const sdf_parser.SdfMolecule,
) ![][]const u8 {
    var scratch_arena = std.heap.ArenaAllocator.init(allocator);
    defer scratch_arena.deinit();
    const scratch = scratch_arena.allocator();

    const stem = sdfFileStem(filename);

    // Plain names, and how many molecules share each of them.
    const plain = try scratch.alloc([]const u8, molecules.len);
    const plain_keys = try scratch.alloc([]const u8, molecules.len);
    var key_counts = std.StringHashMapUnmanaged(usize).empty;
    for (molecules, 0..) |mol, i| {
        plain[i] = if (mol.name.len > 0)
            try std.fmt.allocPrint(scratch, "{s}_{s}", .{ stem, mol.name })
        else
            try std.fmt.allocPrint(scratch, "{s}_{d}", .{ stem, i + 1 });
        plain_keys[i] = try sdfMoleculeNameKey(scratch, plain[i]);
        const entry = try key_counts.getOrPutValue(scratch, plain_keys[i], 0);
        entry.value_ptr.* += 1;
    }

    const names = try allocator.alloc([]const u8, molecules.len);
    var n_names: usize = 0;
    errdefer {
        for (names[0..n_names]) |name| allocator.free(name);
        allocator.free(names);
    }

    for (plain, plain_keys, 0..) |plain_name, plain_key, i| {
        var name = plain_name;
        if (key_counts.get(plain_key).? > 1) {
            // A name ending in this molecule's position cannot equal that of
            // another shared name, only the plain name of a unique molecule.
            while (true) {
                name = try std.fmt.allocPrint(scratch, "{s}_{d}", .{ name, i + 1 });
                const count = key_counts.get(try sdfMoleculeNameKey(scratch, name)) orelse break;
                if (count != 1) break;
            }
        }
        names[i] = try allocator.dupe(u8, name);
        n_names += 1;
    }
    return names;
}

/// Key under which two SDF molecule display names produce the same per-file
/// output: the output name without extension, in ASCII lower case.
fn sdfMoleculeNameKey(allocator: Allocator, display_name: []const u8) ![]const u8 {
    const key = try sdfMoleculeOutputName(allocator, display_name, "");
    for (key) |*c| c.* = std.ascii.toLower(c.*);
    return key;
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
        const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(T){
            .enabled = config.bitmask_correction,
            .coeff = @floatCast(config.bitmask_correction_coeff),
        };
        return if (n_threads > 1)
            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).calculateSasaParallelWithLutAndCorrection(
                allocator,
                input,
                sr_config,
                n_threads,
                lut,
                correction,
            )
        else
            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(T).calculateSasaWithLutAndCorrection(
                allocator,
                input,
                sr_config,
                lut,
                correction,
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
                .trig = config.lr_trig,
            }, n_threads)
        else
            lee_richards.LeeRichardsGen(T).calculateSasa(allocator, input, .{
                .n_slices = config.n_slices,
                .probe_radius = probe_radius,
                .trig = config.lr_trig,
            }),
    };
}

/// Write SASA result to output file
/// Handles f32 -> f64 conversion for consistent output format
/// `name` is an input file name or an SDF molecule display name, as told by
/// `name_source`; the output file name is derived by `perFileOutputName`.
fn writeSasaOutput(
    comptime T: type,
    allocator: Allocator,
    io: std.Io,
    result: *const SasaResultGen(T),
    output_dir: []const u8,
    name_source: OutputNameSource,
    name: []const u8,
    format: OutputFormat,
) !void {
    try validateBatchOutputFormat(format);

    const output_filename = try perFileOutputName(allocator, name_source, name, getOutputExtension(format));
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
    read_parse_time_ns: u64 = 0,
    classifier_time_ns: u64 = 0,

    pub const Status = enum {
        ok,
        err,
    };
};

/// Most failed inputs a failure report lists by name; the rest are counted.
const max_listed_failures = 20;

/// One failed input of a batch run or workflow job: what failed, and why.
const Failure = struct {
    name: []const u8,
    reason: []const u8,

    fn lessThan(_: void, a: Failure, b: Failure) bool {
        return switch (std.mem.order(u8, a.name, b.name)) {
            .lt => true,
            .gt => false,
            .eq => std.mem.lessThan(u8, a.reason, b.reason),
        };
    }
};

/// Where the JSONL rows of a run go, if it writes any. They hold every
/// failure, so a report that lists only some of them points there.
pub const JsonlDestination = union(enum) {
    none,
    stdout,
    file: []const u8,

    fn of(config: BatchConfig, jsonl_output_path: ?[]const u8) JsonlDestination {
        if (!batchWritesJsonl(config)) return .none;
        return if (jsonl_output_path) |path| .{ .file = path } else .stdout;
    }
};

/// The failed inputs of one batch run, or of one job of a workflow.
///
/// This is what a run says about the inputs it could not process. It is
/// printed to stderr at the end of the run whether or not `--quiet` is set:
/// quiet mode suppresses progress and the summary, not errors.
const FailureReport = struct {
    /// Workflow job the inputs belong to; null for a run without jobs.
    job: ?[]const u8 = null,
    /// What is counted: one input file or SDF molecule, one chain selection
    /// of a selection map, or one interface of a BSA analysis.
    unit: Unit = .input,
    total: usize,
    failed: usize,
    /// The failures in the order they are listed. At most
    /// `max_listed_failures` of them are; the slice may also hold fewer than
    /// `failed` entries when a failure could not be recorded.
    failures: []const Failure,
    jsonl: JsonlDestination = .none,

    const Unit = enum {
        input,
        selection,
        interface,

        fn noun(self: Unit, count: usize) []const u8 {
            return switch (self) {
                .input => if (count == 1) "input" else "inputs",
                .selection => if (count == 1) "selection" else "selections",
                .interface => if (count == 1) "interface" else "interfaces",
            };
        }
    };

    /// Write the report; nothing when no input failed.
    ///
    ///     2 of 40 inputs failed:
    ///       bad.pdb: read/parse failed: NoAtomsFound
    ///       empty.cif: read/parse failed: NoAtomSiteLoop
    ///
    /// A workflow job starts with "Job 'name': ". When more inputs failed
    /// than are listed, a last line counts the rest and, for JSONL output,
    /// says where all of them are.
    fn write(self: FailureReport, w: *std.Io.Writer) std.Io.Writer.Error!void {
        if (self.failed == 0) return;

        if (self.job) |job| try w.print("Job '{s}': ", .{job});
        try w.print("{d} of {d} {s} failed:\n", .{ self.failed, self.total, self.unit.noun(self.total) });

        const listed = @min(self.failures.len, max_listed_failures, self.failed);
        for (self.failures[0..listed]) |failure| {
            try w.print("  {s}: {s}\n", .{ failure.name, failure.reason });
        }
        if (self.failed > listed) {
            try w.print("  ... and {d} more", .{self.failed - listed});
            switch (self.jsonl) {
                .none => {},
                .stdout => try w.writeAll(" (every failure is a \"status\":\"err\" row in the JSONL output)"),
                .file => |path| try w.print(" (every failure is a \"status\":\"err\" row in {s})", .{path}),
            }
            try w.writeByte('\n');
        }
    }

    /// The report as text. Caller frees the result.
    fn format(self: FailureReport, allocator: Allocator) ![]u8 {
        var aw = std.Io.Writer.Allocating.init(allocator);
        defer aw.deinit();
        try self.write(&aw.writer);
        return aw.toOwnedSlice();
    }

    /// Print the report to stderr.
    fn print(self: FailureReport, allocator: Allocator) void {
        if (self.failed == 0) return;
        const text = self.format(allocator) catch {
            // Out of memory: the counts at least.
            std.debug.print("{d} of {d} {s} failed\n", .{ self.failed, self.total, self.unit.noun(self.total) });
            return;
        };
        defer allocator.free(text);
        std.debug.print("{s}", .{text});
    }

    /// A copy that owns its strings and holds only the failures that are
    /// listed, for a report that outlives the results it was made from.
    fn dupe(self: FailureReport, arena: Allocator) !FailureReport {
        const listed = @min(self.failures.len, max_listed_failures);
        const failures = try arena.alloc(Failure, listed);
        for (self.failures[0..listed], failures) |failure, *copy| {
            copy.* = .{
                .name = try arena.dupe(u8, failure.name),
                .reason = try arena.dupe(u8, failure.reason),
            };
        }
        var copy = self;
        copy.failures = failures;
        if (self.job) |job| copy.job = try arena.dupe(u8, job);
        switch (self.jsonl) {
            .file => |path| copy.jsonl = .{ .file = try arena.dupe(u8, path) },
            .none, .stdout => {},
        }
        return copy;
    }
};

/// The failures of a run whose rows are written by the workers as they go,
/// without a result per input to read them from afterwards. Safe to use from
/// several threads; `allocator` must be thread-safe.
const FailureLog = struct {
    allocator: Allocator,
    mutex: std.Io.Mutex = .init,
    entries: std.ArrayListUnmanaged(Failure) = .empty,

    fn deinit(self: *FailureLog) void {
        for (self.entries.items) |failure| {
            self.allocator.free(failure.name);
            self.allocator.free(failure.reason);
        }
        self.entries.deinit(self.allocator);
    }

    /// Record that `filename` failed. `id` tells the rows of one file apart
    /// (a selection or an interface); it is left out when it is the file name
    /// itself. Best effort: a failure that cannot be stored is still counted
    /// by the runner and still has its JSONL row.
    fn record(self: *FailureLog, io: std.Io, filename: []const u8, id: ?[]const u8, reason: []const u8) void {
        const row_id: ?[]const u8 = if (id) |value| (if (std.mem.eql(u8, value, filename)) null else value) else null;
        const name = (if (row_id) |value|
            std.fmt.allocPrint(self.allocator, "{s} [{s}]", .{ filename, value })
        else
            self.allocator.dupe(u8, filename)) catch return;
        const owned_reason = self.allocator.dupe(u8, reason) catch {
            self.allocator.free(name);
            return;
        };

        self.mutex.lockUncancelable(io);
        defer self.mutex.unlock(io);
        self.entries.append(self.allocator, .{ .name = name, .reason = owned_reason }) catch {
            self.allocator.free(name);
            self.allocator.free(owned_reason);
        };
    }

    /// The recorded failures ordered by name, so that a report does not
    /// depend on the order in which the workers finished. Call after the
    /// workers have joined.
    fn sorted(self: *FailureLog) []const Failure {
        std.mem.sort(Failure, self.entries.items, {}, Failure.lessThan);
        return self.entries.items;
    }
};

/// Aggregate result for batch processing
pub const BatchResult = struct {
    total_files: usize,
    successful: usize,
    failed: usize,
    total_sasa_time_ns: u64, // SASA calculation only
    total_time_ns: u64, // Including I/O
    scan_time_ns: u64 = 0,
    build_items_time_ns: u64 = 0,
    process_time_ns: u64 = 0,
    read_parse_time_ns: u64 = 0,
    classifier_time_ns: u64 = 0,
    jsonl_write_time_ns: u64 = 0,
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

    /// The failure report of this run: the failed results in input order.
    /// `buffer` holds the listed failures; they refer to the results.
    fn failureReport(
        self: BatchResult,
        buffer: *[max_listed_failures]Failure,
        job: ?[]const u8,
        jsonl: JsonlDestination,
    ) FailureReport {
        var listed: usize = 0;
        for (self.file_results) |file_result| {
            if (listed == buffer.len) break;
            if (file_result.status != .err) continue;
            buffer[listed] = .{ .name = file_result.filename, .reason = file_result.error_msg orelse "unknown error" };
            listed += 1;
        }
        return .{
            .job = job,
            .total = self.total_files,
            .failed = self.failed,
            .failures = buffer[0..listed],
            .jsonl = jsonl,
        };
    }

    /// Print the failed inputs to stderr; nothing when every input succeeded.
    /// Not part of the summary: it is printed in quiet mode too.
    pub fn printFailures(self: BatchResult, jsonl: JsonlDestination) void {
        var buffer: [max_listed_failures]Failure = undefined;
        self.failureReport(&buffer, null, jsonl).print(self.allocator);
    }

    /// Print human-readable summary, including the failed inputs
    pub fn printSummary(self: BatchResult, show_timing: bool, jsonl: JsonlDestination) void {
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

        // The failed inputs, as quiet mode reports them
        if (self.failed > 0) {
            std.debug.print("\n", .{});
            self.printFailures(jsonl);
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
        const scan_ms = @as(f64, @floatFromInt(self.scan_time_ns)) / ns_to_ms;
        const build_ms = @as(f64, @floatFromInt(self.build_items_time_ns)) / ns_to_ms;
        const process_ms = @as(f64, @floatFromInt(self.process_time_ns)) / ns_to_ms;

        std.debug.print("BATCH_SASA_TIME_MS:{d:.2}\n", .{total_sasa_ms});
        std.debug.print("BATCH_TOTAL_TIME_MS:{d:.2}\n", .{total_ms});
        std.debug.print("BATCH_SCAN_TIME_MS:{d:.2}\n", .{scan_ms});
        std.debug.print("BATCH_BUILD_ITEMS_TIME_MS:{d:.2}\n", .{build_ms});
        std.debug.print("BATCH_PROCESS_TIME_MS:{d:.2}\n", .{process_ms});
        std.debug.print("BATCH_FILES:{d}\n", .{self.total_files});
        std.debug.print("BATCH_SUCCESS:{d}\n", .{self.successful});
    }

    pub fn printStageProfile(self: BatchResult) void {
        const ns_to_ms = 1_000_000.0;
        const read_parse_ms = @as(f64, @floatFromInt(self.read_parse_time_ns)) / ns_to_ms;
        const classifier_ms = @as(f64, @floatFromInt(self.classifier_time_ns)) / ns_to_ms;
        const jsonl_write_ms = @as(f64, @floatFromInt(self.jsonl_write_time_ns)) / ns_to_ms;
        std.debug.print("BATCH_READ_PARSE_TIME_MS:{d:.2}\n", .{read_parse_ms});
        std.debug.print("BATCH_CLASSIFIER_TIME_MS:{d:.2}\n", .{classifier_ms});
        std.debug.print("BATCH_JSONL_WRITE_TIME_MS:{d:.2}\n", .{jsonl_write_ms});
    }
};

const WorkflowJobState = struct {
    name: []const u8,
    config: BatchConfig,
    output_dir: ?[]const u8 = null,
    jsonl_output_path: ?[]const u8 = null,
    successful: usize = 0,
    failed: usize = 0,
    total_sasa_time_ns: u64 = 0,
    /// The failed inputs of the job, for the report at the end of the run.
    failures: FailureLog,

    fn deinit(self: *WorkflowJobState, allocator: Allocator) void {
        if (self.output_dir) |path| allocator.free(path);
        if (self.jsonl_output_path) |path| allocator.free(path);
        self.failures.deinit();
    }

    /// Call when every input has been processed and the counters are final.
    fn failureReport(self: *WorkflowJobState) FailureReport {
        return .{
            .job = self.name,
            .total = self.successful + self.failed,
            .failed = self.failed,
            .failures = self.failures.sorted(),
            .jsonl = JsonlDestination.of(self.config, self.jsonl_output_path),
        };
    }
};

/// The last lines of a file-first workflow: the totals, then the failed
/// inputs of each job.
fn printWorkflowJobStates(allocator: Allocator, states: []WorkflowJobState) void {
    var successful: usize = 0;
    var failed: usize = 0;
    for (states) |state| {
        successful += state.successful;
        failed += state.failed;
    }
    std.debug.print("Workflow complete: {d} successful, {d} failed\n", .{ successful, failed });
    for (states) |*state| state.failureReport().print(allocator);
}

const WorkflowJobCounter = struct {
    successful: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    failed: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    total_sasa_time_ns: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),
};

const WorkflowJobRuntime = struct {
    state: *WorkflowJobState,
    jsonl_stream: ?JsonlStreamWriter = null,
    jsonl_buffer: [64 * 1024]u8 = undefined,
    jsonl_file: ?std.Io.File = null,
    jsonl_file_needs_close: bool = false,
    counter: WorkflowJobCounter = .{},

    fn close(self: *WorkflowJobRuntime, io: std.Io) void {
        if (self.jsonl_stream) |*stream| {
            stream.flush() catch |err| {
                logWarning("workflow JSONL flush failed for {s}: {s}", .{ self.state.name, @errorName(err) });
            };
        }
        if (self.jsonl_file_needs_close) {
            if (self.jsonl_file) |file| file.close(io);
        }
    }
};

fn chainMatchesFilter(chain: types.FixedString4, chains: []const []const u8) bool {
    for (chains) |target| {
        if (chain.eqlSlice(target)) return true;
    }
    return false;
}

fn chainFullMatchesFilter(chain: []const u8, chains: []const []const u8) bool {
    for (chains) |target| {
        if (std.mem.eql(u8, chain, target)) return true;
    }
    return false;
}

fn atomSelectedByChains(input: AtomInput, index: usize, chains: ?[]const []const u8) bool {
    const filter = chains orelse return true;
    if (filter.len == 0) return true;
    if (input.chain_id_full) |chain_ids_full| {
        return chainFullMatchesFilter(chain_ids_full[index], filter);
    }
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
    const chain_id_full = if (input.chain_id_full != null) try allocator.alloc([]const u8, selected_count) else null;
    var chain_id_full_copied: usize = 0;
    errdefer if (chain_id_full) |v| {
        for (v[0..chain_id_full_copied]) |chain| allocator.free(chain);
        allocator.free(v);
    };
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
        if (chain_id_full) |v| {
            v[out_i] = try allocator.dupe(u8, input.chain_id_full.?[i]);
            chain_id_full_copied += 1;
        }
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
        .chain_id_full = chain_id_full,
        .residue_num = residue_num,
        .insertion_code = insertion_code,
        .allocator = allocator,
    };
}

/// Read input file with auto-format detection
const ParsedInput = struct {
    input: AtomInput,
    inline_ccd: ?ccd_parser.ComponentDict = null,

    fn deinit(self: *ParsedInput) void {
        self.input.deinit();
        if (self.inline_ccd) |*dict| {
            dict.deinit();
            self.inline_ccd = null;
        }
    }

    fn inlineCcdPtr(self: *const ParsedInput) ?*const ccd_parser.ComponentDict {
        if (self.inline_ccd != null) return &self.inline_ccd.?;
        return null;
    }
};

fn readInputFile(allocator: Allocator, io: std.Io, path: []const u8, config: BatchConfig) !ParsedInput {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => .{ .input = try json_parser.readAtomInputFromFile(allocator, io, path) },
        .bcif => blk: {
            var parser = bcif_parser.BcifParser.init(allocator);
            errdefer parser.deinitCcd();
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            parser.chain_filter = config.chain_filter;
            parser.use_auth_chain = config.use_auth_chain;
            parser.alt_loc_mode = config.alt_loc_mode;
            parser.alt_loc_id = config.alt_loc_id;
            parser.parse_inline_ccd = classifierUsesCcdResources(config.classifier_type);
            const input = try parser.parseFileWithInputIo(io, path, config.input_io);
            break :blk .{ .input = input, .inline_ccd = parser.takeInlineCcd() };
        },
        .mmcif => blk: {
            if (shouldTryAfModelFastParser(config)) {
                const input = af_model_parser.parseFileWithOptions(allocator, io, path, .{
                    .io_mode = config.input_io.resolve(.read),
                }) catch |err| if (shouldFallbackAfModelFastError(err)) null else return err;
                if (input) |fast_input| {
                    break :blk .{ .input = fast_input };
                }
            }
            var parser = mmcif_parser.MmcifParser.init(allocator);
            errdefer parser.deinitCcd();
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            parser.chain_filter = config.chain_filter;
            parser.use_auth_chain = config.use_auth_chain;
            parser.alt_loc_mode = config.alt_loc_mode;
            parser.alt_loc_id = config.alt_loc_id;
            parser.parse_inline_ccd = classifierUsesCcdResources(config.classifier_type);
            const input = try parser.parseFileWithInputIo(io, path, config.input_io);
            break :blk .{ .input = input, .inline_ccd = parser.takeInlineCcd() };
        },
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.skip_hydrogens = !config.include_hydrogens;
            parser.atom_only = !config.include_hetatm;
            parser.chain_filter = config.chain_filter;
            parser.alt_loc_mode = config.alt_loc_mode;
            parser.alt_loc_id = config.alt_loc_id;
            break :blk .{ .input = try parser.parseFileWithInputIo(io, path, config.input_io) };
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

            break :blk .{ .input = try sdf_parser.toAtomInput(allocator, molecules, !config.include_hydrogens) };
        },
    };
}

/// Whether the AlphaFold-model fast parser may read an mmCIF input. It only
/// reads the leading ATOM rows by label chain, so every option that asks for
/// anything else goes to the generic parser: HETATM rows that follow the
/// ATOM rows would otherwise be dropped although they were requested.
fn shouldTryAfModelFastParser(config: BatchConfig) bool {
    return config.af_model_fast and
        config.chain_filter == null and
        !config.use_auth_chain and
        !config.include_hydrogens and
        !config.include_hetatm and
        config.alt_loc_mode == .auto;
}

fn shouldFallbackAfModelFastError(err: anyerror) bool {
    return err == error.UnsupportedLayout;
}

/// Apply built-in classifier to replace radii based on residue/atom names
fn applyBuiltinClassifier(
    input: *AtomInput,
    ct: ClassifierType,
    sdf_ccd: ?*const ccd_parser.ComponentDict,
    inline_ccd: ?*const ccd_parser.ComponentDict,
    external_ccd: ?*const ccd_parser.ComponentDict,
) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    // CCD and ProtOr share the static ProtOr-compatible table. Only CCD may
    // extend it with runtime component topology.
    var ccd_clf: ?classifier_ccd.CcdClassifier = if (ct == .ccd or ct == .protor) classifier_ccd.CcdClassifier.init(input.allocator) else null;
    defer if (ccd_clf) |*c| c.deinit();

    if (ct == .ccd and ccd_clf != null) {
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
            const dicts: [3]?*const ccd_parser.ComponentDict = .{ sdf_ccd, inline_ccd, external_ccd };
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

        // Not in the tables: element-based radius, atom name-based without an
        // element, else the input radius
        new_radii[i] = maybe_radius orelse classifier.guessFallbackRadius(
            ct,
            if (input.element) |elements| elements[i] else null,
            residues[i].slice(),
            atom_names[i].slice(),
        ) orelse input.r[i];
    }

    input.allocator.free(input.r);
    input.r = new_radii;
}

/// Apply the CCD classifier to one SDF/MOL molecule: radii from the
/// molecule's own bond topology (`sdf_parser.applyTopologyRadii`).
fn applySdfTopologyClassifier(input: *AtomInput, own_component: *const ccd_parser.StoredComponent) !void {
    const view = own_component.view();
    _ = try sdf_parser.applyTopologyRadii(input, &view);
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
        } else if (classifier.guessRadiusFromResidueAtom(residues[i].slice(), atom_names[i].slice())) |r| {
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

fn jsonlOptions(config: BatchConfig) json_writer.JsonlOptions {
    return .{
        .decimals = config.jsonl_decimals,
        .include_atom_areas = config.jsonl_include_atom_areas,
        .include_atom_identity = config.jsonl_include_atom_identity,
        .include_total_area = config.jsonl_include_total_area,
    };
}

/// Whether the atom areas of each result are kept for its JSONL row.
fn batchShouldStoreAtomAreas(config: BatchConfig) bool {
    return config.output_format == .jsonl and config.jsonl_include_atom_areas;
}

/// Whether a run writes one JSONL row per input instead of one output file
/// per input. This is a property of the output format alone: JSONL without
/// atom areas (`[output.jsonl] atom_areas = false`) still writes rows, to the
/// JSONL file or to stdout, and never writes per-file outputs.
fn batchWritesJsonl(config: BatchConfig) bool {
    return config.output_format == .jsonl;
}

fn classifierTypeName(classifier_type: ?ClassifierType) []const u8 {
    const ct = classifier_type orelse return "none";
    return switch (ct) {
        .naccess => "naccess",
        .protor => "protor",
        .oons => "oons",
        .ccd => "ccd",
    };
}

fn precisionName(precision: Precision) []const u8 {
    return switch (precision) {
        .f32 => "f32",
        .f64 => "f64",
    };
}

fn algorithmName(algorithm: Algorithm) []const u8 {
    return switch (algorithm) {
        .sr => "sr",
        .lr => "lr",
    };
}

fn writeWorkflowJsonlMetadata(
    allocator: Allocator,
    io: std.Io,
    path: []const u8,
    name: []const u8,
    config: BatchConfig,
) !void {
    if (config.jsonl_metadata != .sidecar) return;

    const JsonlMeta = struct {
        atom_areas: bool,
        atom_identity: bool,
        total_area: bool,
        decimals: ?u8,
    };
    const CalculationMeta = struct {
        algorithm: []const u8,
        n_points: u32,
        n_slices: u32,
        probe_radius: f64,
        precision: []const u8,
        classifier: []const u8,
    };
    const Metadata = struct {
        name: []const u8,
        job: []const u8,
        format: []const u8,
        metadata: []const u8,
        jsonl: JsonlMeta,
        calculation: CalculationMeta,
    };

    const meta = Metadata{
        .name = name,
        .job = name,
        .format = "jsonl",
        .metadata = "sidecar",
        .jsonl = .{
            .atom_areas = config.jsonl_include_atom_areas,
            .atom_identity = config.jsonl_include_atom_identity,
            .total_area = config.jsonl_include_total_area,
            .decimals = config.jsonl_decimals,
        },
        .calculation = .{
            .algorithm = algorithmName(config.algorithm),
            .n_points = config.n_points,
            .n_slices = config.n_slices,
            .probe_radius = config.probe_radius,
            .precision = precisionName(config.precision),
            .classifier = classifierTypeName(config.classifier_type),
        },
    };

    const text = try json_writer.stringifyFinite(allocator, meta, .{ .whitespace = .indent_2 });
    defer allocator.free(text);
    try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = path, .data = text });
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

/// Calculate SASA and optionally write output for an already prepared input.
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
    result.total_sasa = sasa_result.total_area;

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
        if (!batchWritesJsonl(config)) {
            writeSasaOutput(T, arena, io, &sasa_result, out_dir, .input_file, filename, config.output_format) catch |err| {
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "output write failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
        }
    }

    return result;
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

    var file_config = config;
    if (config.chain_map) |map| {
        const selection = map.get(filename) orelse {
            result.status = .err;
            result.error_msg = std.fmt.allocPrint(result_allocator, "chain map entry not found", .{}) catch null;
            return result;
        };
        file_config.chain_filter = selection.chains;
        file_config.use_auth_chain = selection.asym_id_type == .auth;
    }

    // Build input path
    const input_path = std.fs.path.join(arena, &.{ input_dir, filename }) catch |err| {
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "path join failed: {s}", .{@errorName(err)}) catch null;
        return result;
    };

    // Read and parse input (auto-detect format from extension)
    var read_parse_timer: std.Io.Timestamp = undefined;
    if (file_config.profile_stages) read_parse_timer = std.Io.Timestamp.now(io, .awake);
    var parsed = readInputFile(arena, io, input_path, file_config) catch |err| {
        if (file_config.profile_stages) result.read_parse_time_ns = @intCast(read_parse_timer.untilNow(io, .awake).nanoseconds);
        result.status = .err;
        result.error_msg = std.fmt.allocPrint(result_allocator, "read/parse failed: {s}", .{@errorName(err)}) catch null;
        return result;
    };
    if (file_config.profile_stages) result.read_parse_time_ns = @intCast(read_parse_timer.untilNow(io, .awake).nanoseconds);
    defer parsed.deinit();

    // Apply classifier for PDB/mmCIF input (skip JSON unless classification info exists).
    var classifier_timer: std.Io.Timestamp = undefined;
    if (file_config.profile_stages) classifier_timer = std.Io.Timestamp.now(io, .awake);
    if (file_config.custom_classifier) |custom_classifier| {
        if (parsed.input.hasClassificationInfo()) {
            applyCustomClassifier(&parsed.input, custom_classifier, file_config.quiet) catch |err| {
                if (file_config.profile_stages) result.classifier_time_ns = @intCast(classifier_timer.untilNow(io, .awake).nanoseconds);
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "classifier failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
        }
    } else if (file_config.classifier_type) |ct| {
        const format = format_detect.detectInputFormat(input_path);
        if (format != .json and parsed.input.hasClassificationInfo()) {
            applyBuiltinClassifier(&parsed.input, ct, file_config.sdf_ccd, parsed.inlineCcdPtr(), file_config.external_ccd) catch |err| {
                if (file_config.profile_stages) result.classifier_time_ns = @intCast(classifier_timer.untilNow(io, .awake).nanoseconds);
                result.status = .err;
                result.error_msg = std.fmt.allocPrint(result_allocator, "classifier failed: {s}", .{@errorName(err)}) catch null;
                return result;
            };
        }
    }
    if (file_config.profile_stages) result.classifier_time_ns = @intCast(classifier_timer.untilNow(io, .awake).nanoseconds);

    var calc_result = switch (file_config.precision) {
        .f64 => calculatePreparedInputResult(
            f64,
            arena,
            io,
            result_allocator,
            parsed.input,
            output_dir,
            filename,
            file_config,
            n_threads,
            lut_f64,
            coarse_lut_f64,
            fine_lut_f64,
        ),
        .f32 => calculatePreparedInputResult(
            f32,
            arena,
            io,
            result_allocator,
            parsed.input,
            output_dir,
            filename,
            file_config,
            n_threads,
            lut_f32,
            coarse_lut_f32,
            fine_lut_f32,
        ),
    };
    calc_result.read_parse_time_ns = result.read_parse_time_ns;
    calc_result.classifier_time_ns = result.classifier_time_ns;
    return calc_result;
}

/// Process a single SDF molecule and return result.
/// `display_name` is used as the filename in the result and names the
/// per-file output (see `sdfMoleculeDisplayNames`).
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

    // The CCD classifier takes the radii of the molecule from its own bond
    // topology, whatever its title is (it may be blank)
    var sdf_component: ?ccd_parser.StoredComponent = null;
    if (classifierUsesCcdResources(config.classifier_type)) {
        sdf_component = sdf_parser.toStoredComponent(arena, molecule) catch |err| blk: {
            logWarning("{s}: failed to build SDF component: {s}", .{ display_name, @errorName(err) });
            break :blk null;
        };
    }
    defer if (sdf_component) |*c| c.deinit();

    const sdf_component_ptr: ?*const ccd_parser.StoredComponent = if (sdf_component) |*c| c else null;
    return processOneSdfMoleculeInner(arena, io, result_allocator, &result, &input, output_dir, display_name, config, n_threads, lut_f64, lut_f32, coarse_lut_f64, fine_lut_f64, coarse_lut_f32, fine_lut_f32, sdf_component_ptr);
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
    sdf_component: ?*const ccd_parser.StoredComponent,
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
            // The molecule is its own component definition: its atoms are
            // matched to its own bond topology, not looked up by residue
            // name, so --sdf and --ccd have nothing to add.
            const classified: anyerror!void = if (sdf_component) |own_component|
                applySdfTopologyClassifier(input, own_component)
            else
                applyBuiltinClassifier(input, ct, config.sdf_ccd, null, config.external_ccd);
            classified catch |err| {
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
                if (!batchWritesJsonl(config)) {
                    writeSasaOutput(f64, arena, io, &sasa_result, out_dir, .sdf_molecule, display_name, config.output_format) catch |err| {
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
            total_area = sasa_result.total_area;

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
                if (!batchWritesJsonl(config)) {
                    writeSasaOutput(f32, arena, io, &sasa_result, out_dir, .sdf_molecule, display_name, config.output_format) catch |err| {
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
/// Each call to writeResult acquires the mutex and serializes one line into a
/// persistent buffered file writer. The owner must call flush after workers join.
const JsonlStreamWriter = struct {
    const buffer_size = 64 * 1024;
    const large_line_threshold = buffer_size / 2;

    mutex: std.Io.Mutex = .init,
    file: std.Io.File,
    io: std.Io,
    options: json_writer.JsonlOptions = .{},
    writer: std.Io.File.Writer,
    /// Set to true if any JSONL serialization or write failed.
    write_failed: std.atomic.Value(bool) = std.atomic.Value(bool).init(false),

    pub fn init(
        file: std.Io.File,
        io: std.Io,
        options: json_writer.JsonlOptions,
        buffer: *[buffer_size]u8,
    ) JsonlStreamWriter {
        return .{
            .file = file,
            .io = io,
            .options = options,
            .writer = std.Io.File.Writer.initStreaming(file, io, buffer),
        };
    }

    fn recordFailure(self: *JsonlStreamWriter) void {
        self.write_failed.store(true, .release);
    }

    fn writeLine(self: *JsonlStreamWriter, filename: []const u8, line: []const u8) !void {
        self.mutex.lockUncancelable(self.io);
        defer self.mutex.unlock(self.io);

        // Human AFDB PDB JSONL lines are often close to the 64 KiB buffer size.
        // Keeping such lines in a persistent shared buffer causes extra buffer
        // churn; use the previous per-line streaming writer behavior for large
        // lines while retaining persistent buffering for smaller rounded JSONL.
        if (line.len > large_line_threshold) {
            self.writer.interface.flush() catch |err| {
                logWarning("JSONL flush failed before large write for {s}: {s}", .{ filename, @errorName(err) });
                self.recordFailure();
                return err;
            };
            var local_buf: [buffer_size]u8 = undefined;
            var local_writer = std.Io.File.Writer.initStreaming(self.file, self.io, &local_buf);
            local_writer.interface.writeAll(line) catch |err| {
                logWarning("JSONL large write failed for {s}: {s}", .{ filename, @errorName(err) });
                self.recordFailure();
                return err;
            };
            local_writer.interface.writeAll("\n") catch |err| {
                logWarning("JSONL large newline write failed for {s}: {s}", .{ filename, @errorName(err) });
                self.recordFailure();
                return err;
            };
            local_writer.interface.flush() catch |err| {
                logWarning("JSONL large flush failed for {s}: {s}", .{ filename, @errorName(err) });
                self.recordFailure();
                return err;
            };
            return;
        }

        self.writer.interface.writeAll(line) catch |err| {
            logWarning("JSONL write failed for {s}: {s}", .{ filename, @errorName(err) });
            self.recordFailure();
            return err;
        };
        self.writer.interface.writeAll("\n") catch |err| {
            logWarning("JSONL newline write failed for {s}: {s}", .{ filename, @errorName(err) });
            self.recordFailure();
            return err;
        };
    }

    /// Serialize and write one JSONL line for a completed file result.
    /// alloc is a short-lived allocator (e.g., thread-local arena) used only for
    /// the serialized string; it is freed by the caller's arena reset.
    pub fn writeResult(
        self: *JsonlStreamWriter,
        alloc: Allocator,
        result: *FileResult,
    ) void {
        const line = fileResultToJsonlLineOptions(alloc, result, self.options) catch |err| {
            logWarning("failed to serialize {s}: {s}", .{ result.filename, @errorName(err) });
            self.recordFailure();
            return;
        };
        // line is on alloc (arena); no explicit free needed — arena reset handles it.

        self.writeLine(result.filename, line) catch {
            self.recordFailure();
        };
    }

    pub fn flush(self: *JsonlStreamWriter) !void {
        self.mutex.lockUncancelable(self.io);
        defer self.mutex.unlock(self.io);
        try self.writer.interface.flush();
    }

    /// Returns true if any write to the JSONL file failed.
    pub fn hasError(self: *const JsonlStreamWriter) bool {
        return self.write_failed.load(.acquire);
    }
};

/// Write a JSONL line for a result via the buffered writer.
/// atom_areas live on the arena and are invalidated after arena reset.
fn fileResultToJsonlLine(allocator: Allocator, result: *FileResult) ![]u8 {
    return fileResultToJsonlLineOptions(allocator, result, .{});
}

fn fileResultToJsonlLineOptions(allocator: Allocator, result: *FileResult, options: json_writer.JsonlOptions) ![]u8 {
    if (result.status == .err) {
        return json_writer.fileErrorToJsonlLine(allocator, result.filename, result.error_msg orelse "unknown error");
    }
    if (result.residue_map) |map| {
        const areas = result.atom_areas orelse if (options.include_atom_areas) return error.MissingAtomAreas else &.{};
        return json_writer.fileResultWithResidueMapToJsonlLineOptions(allocator, result.filename, result.total_sasa, areas, map, options);
    }
    const areas = result.atom_areas orelse if (options.include_atom_areas) return error.MissingAtomAreas else &.{};
    return json_writer.fileResultToJsonlLineOptions(allocator, result.filename, result.total_sasa, areas, options);
}

fn writeJsonlResult(
    jsonl_writer: *std.Io.File.Writer,
    arena_alloc: Allocator,
    result: *FileResult,
    options: json_writer.JsonlOptions,
) !void {
    const line = fileResultToJsonlLineOptions(arena_alloc, result, options) catch |err| {
        logWarning("failed to serialize {s}: {s}", .{ result.filename, @errorName(err) });
        return error.JsonlWriteFailed;
    };
    try jsonl_writer.interface.writeAll(line);
    try jsonl_writer.interface.writeAll("\n");
    try jsonl_writer.interface.flush();
}

fn truncateJsonlOutput(io: std.Io, path: []const u8) !void {
    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
    file.close(io);
}

/// Opens an existing JSONL file for `appendJsonlResultToFile`. The file is
/// opened for reading as well as writing: finding its end queries the file
/// size, which Windows refuses on a write-only handle (`AccessDenied`).
fn openJsonlForAppend(io: std.Io, path: []const u8) !std.Io.File {
    return std.Io.Dir.cwd().openFile(io, path, .{ .mode = .read_write });
}

/// Appends one row to `file`, which must come from `openJsonlForAppend`.
fn appendJsonlResultToFile(io: std.Io, file: std.Io.File, allocator: Allocator, result: *FileResult, options: json_writer.JsonlOptions) !void {
    const line = try fileResultToJsonlLineOptions(allocator, result, options);
    defer allocator.free(line);

    const file_len = try file.length(io);
    var write_buf: [64 * 1024]u8 = undefined;
    var writer = std.Io.File.Writer.init(file, io, &write_buf);
    try writer.seekTo(file_len);
    try writer.interface.writeAll(line);
    try writer.interface.writeByte('\n');
    try writer.interface.flush();
}

fn writeBsaAnalysisJsonl(
    stream: *JsonlStreamWriter,
    allocator: Allocator,
    row: json_writer.BsaAnalysisJsonl,
    options: json_writer.JsonlOptions,
) !void {
    const line = try json_writer.bsaAnalysisToJsonlLineOptions(allocator, row, options);
    defer allocator.free(line);
    try stream.writeLine(row.filename, line);
}

fn writeBsaAnalysisErrorJsonl(
    stream: *JsonlStreamWriter,
    allocator: Allocator,
    row: json_writer.BsaAnalysisErrorJsonl,
) !void {
    const line = try json_writer.bsaAnalysisErrorToJsonlLine(allocator, row);
    defer allocator.free(line);
    try stream.writeLine(row.filename, line);
}

fn analysisName(analysis: workflow_manifest.Analysis) []const u8 {
    return analysis.name orelse "bsa";
}

fn analysisLevel(analysis: workflow_manifest.Analysis) []const u8 {
    return analysis.level orelse "total";
}

fn appendChainGroups(allocator: Allocator, a: []const []const u8, b: []const []const u8) ![]const []const u8 {
    const out = try allocator.alloc([]const u8, a.len + b.len);
    @memcpy(out[0..a.len], a);
    @memcpy(out[a.len..], b);
    return out;
}

const BsaResidueDeltaArrays = struct {
    residue_partner: []const []const u8 = &.{},
    residue_chain: []const []const u8 = &.{},
    residue_name: []const []const u8 = &.{},
    residue_number: []const i32 = &.{},
    residue_insertion_code: []const []const u8 = &.{},
    residue_sasa_isolated: []const f64 = &.{},
    residue_sasa_complex: []const f64 = &.{},
    residue_delta_sasa: []const f64 = &.{},
};

fn buildBsaResidueDeltaArrays(
    allocator: Allocator,
    input: AtomInput,
    atom_sasa_isolated: []const f64,
    atom_sasa_complex: []const f64,
    partner_a: []const []const u8,
) !BsaResidueDeltaArrays {
    var isolated_map = try json_writer.buildResidueMap(allocator, input, atom_sasa_isolated);
    defer isolated_map.deinit();
    var complex_map = try json_writer.buildResidueMap(allocator, input, atom_sasa_complex);
    defer complex_map.deinit();
    if (isolated_map.len() != complex_map.len()) return error.InvalidResidueMap;

    const residue_partner = try allocator.alloc([]const u8, isolated_map.len());
    const residue_chain = try allocator.alloc([]const u8, isolated_map.len());
    const residue_name = try allocator.alloc([]const u8, isolated_map.len());
    const residue_insertion_code = try allocator.alloc([]const u8, isolated_map.len());
    const residue_number = try allocator.dupe(i32, isolated_map.residue_number);
    const residue_sasa_isolated = try allocator.dupe(f64, isolated_map.residue_sasa);
    const residue_sasa_complex = try allocator.dupe(f64, complex_map.residue_sasa);
    const residue_delta_sasa = try allocator.alloc(f64, isolated_map.len());

    for (0..isolated_map.len()) |i| {
        residue_partner[i] = if (atomSelectedByChains(input, isolated_map.residue_atom_start[i], partner_a)) "a" else "b";
        // The full chain ID when the input has one, as in the atom arrays:
        // the fixed-width chain holds at most four characters.
        residue_chain[i] = try allocator.dupe(u8, if (isolated_map.residue_chain_full) |chains|
            chains[i]
        else
            isolated_map.residue_chain[i].slice());
        residue_name[i] = try allocator.dupe(u8, isolated_map.residue_name[i].slice());
        residue_insertion_code[i] = try allocator.dupe(u8, isolated_map.residue_insertion_code[i].slice());
        residue_delta_sasa[i] = residue_sasa_isolated[i] - residue_sasa_complex[i];
    }

    return .{
        .residue_partner = residue_partner,
        .residue_chain = residue_chain,
        .residue_name = residue_name,
        .residue_number = residue_number,
        .residue_insertion_code = residue_insertion_code,
        .residue_sasa_isolated = residue_sasa_isolated,
        .residue_sasa_complex = residue_sasa_complex,
        .residue_delta_sasa = residue_delta_sasa,
    };
}

const BsaAtomArrays = struct {
    atom_index: []const usize = &.{},
    atom_partner: []const []const u8 = &.{},
    atom_chain: []const []const u8 = &.{},
    atom_residue_name: []const []const u8 = &.{},
    atom_residue_number: []const i32 = &.{},
    atom_insertion_code: []const []const u8 = &.{},
    atom_name: []const []const u8 = &.{},
    atom_element: []const []const u8 = &.{},
};

fn buildBsaAtomArrays(allocator: Allocator, input: AtomInput, partner_a: []const []const u8) !BsaAtomArrays {
    if (!input.hasResidueInfo() or input.atom_name == null or input.element == null) return error.MissingAtomMetadata;

    const atom_count = input.atomCount();
    const atom_index = try allocator.alloc(usize, atom_count);
    const atom_partner = try allocator.alloc([]const u8, atom_count);
    const atom_chain = try allocator.alloc([]const u8, atom_count);
    const atom_residue_name = try allocator.alloc([]const u8, atom_count);
    const atom_residue_number = try allocator.dupe(i32, input.residue_num.?);
    const atom_insertion_code = try allocator.alloc([]const u8, atom_count);
    const atom_name = try allocator.alloc([]const u8, atom_count);
    const atom_element = try allocator.alloc([]const u8, atom_count);

    for (0..atom_count) |i| {
        atom_index[i] = i;
        atom_partner[i] = if (atomSelectedByChains(input, i, partner_a)) "a" else "b";
        atom_chain[i] = try allocator.dupe(u8, if (input.chain_id_full) |chains| chains[i] else input.chain_id.?[i].slice());
        atom_residue_name[i] = try allocator.dupe(u8, input.residue.?[i].slice());
        atom_insertion_code[i] = try allocator.dupe(u8, input.insertion_code.?[i].slice());
        atom_name[i] = try allocator.dupe(u8, input.atom_name.?[i].slice());
        atom_element[i] = element_module.fromAtomicNumber(input.element.?[i]).symbol();
    }

    return .{
        .atom_index = atom_index,
        .atom_partner = atom_partner,
        .atom_chain = atom_chain,
        .atom_residue_name = atom_residue_name,
        .atom_residue_number = atom_residue_number,
        .atom_insertion_code = atom_insertion_code,
        .atom_name = atom_name,
        .atom_element = atom_element,
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
    try validateBatchOutputFormat(config.output_format);

    var prepared = try prepareBatch(allocator, io, input_dir, output_dir, config);
    defer prepared.deinit(allocator);

    return runPreparedSequential(allocator, io, input_dir, output_dir, config, jsonl_output_path, &prepared);
}

/// Process the work items of a prepared batch one after another.
fn runPreparedSequential(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
    prepared: *const PreparedBatch,
) !BatchResult {
    const work_items = prepared.work.items.items;

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing files", .estimated_total_items = work_items.len })
    else
        .none;
    defer progress_root.end();
    const progress_node: ?std.Progress.Node = if (shouldShowProgress(config)) progress_root else null;

    // Create output directory if specified
    if (output_dir) |out_dir| {
        try std.Io.Dir.cwd().createDirPath(io, out_dir);
    }

    // Sequential mode writes JSONL as each file completes.
    var jsonl_file: ?std.Io.File = null;
    var jsonl_file_needs_close = false;
    if (jsonl_output_path) |path| {
        jsonl_file = try std.Io.Dir.cwd().createFile(io, path, .{});
        jsonl_file_needs_close = true;
    } else if (batchWritesJsonl(config)) {
        jsonl_file = std.Io.File.stdout();
    }
    defer if (jsonl_file_needs_close) {
        if (jsonl_file) |f| f.close(io);
    };

    // One result per work item (SDF files expand into one item per molecule)
    var results_list = std.ArrayListUnmanaged(FileResult).empty;
    defer results_list.deinit(allocator);
    try results_list.ensureTotalCapacity(allocator, work_items.len);

    // Process each item
    var total_sasa_time_ns: u64 = 0;
    var total_read_parse_time_ns: u64 = 0;
    var total_classifier_time_ns: u64 = 0;
    var total_jsonl_write_time_ns: u64 = 0;
    var successful: usize = 0;
    var failed: usize = 0;

    // Build bitmask LUT once (if enabled)
    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    // Use arena allocator for each item (reset between items)
    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    // JSONL buffered writer: created once, reused across all files.
    // Uses streaming mode so OS seek position advances with each write.
    var jsonl_write_buf: [64 * 1024]u8 = undefined;
    var jsonl_writer: ?std.Io.File.Writer = if (jsonl_file) |jf|
        std.Io.File.Writer.initStreaming(jf, io, &jsonl_write_buf)
    else
        null;

    var process_timer = std.Io.Timestamp.now(io, .awake);
    for (work_items) |item| {
        const name_copy = try allocator.dupe(u8, item.display_name);

        var result = processWorkItem(
            arena.allocator(),
            io,
            allocator,
            input_dir,
            output_dir,
            item,
            name_copy,
            config,
            luts.f64Ptr(),
            luts.f32Ptr(),
            luts.coarseF64Ptr(),
            luts.fineF64Ptr(),
            luts.coarseF32Ptr(),
            luts.fineF32Ptr(),
        );

        if (result.status == .ok) {
            successful += 1;
            total_sasa_time_ns += result.sasa_time_ns;
            total_read_parse_time_ns += result.read_parse_time_ns;
            total_classifier_time_ns += result.classifier_time_ns;
        } else {
            failed += 1;
        }

        // Stream JSONL output
        if (jsonl_writer) |*w| {
            var write_timer: std.Io.Timestamp = undefined;
            if (config.profile_stages) write_timer = std.Io.Timestamp.now(io, .awake);
            try writeJsonlResult(w, arena.allocator(), &result, jsonlOptions(config));
            if (config.profile_stages) total_jsonl_write_time_ns += @intCast(write_timer.untilNow(io, .awake).nanoseconds);
        }
        result.atom_areas = null;
        result.residue_map = null;

        results_list.appendAssumeCapacity(result);

        _ = arena.reset(.retain_capacity);

        if (progress_node) |node| {
            node.completeOne();
        }
    }

    const process_time_ns: u64 = @intCast(process_timer.untilNow(io, .awake).nanoseconds);
    const total_time_ns: u64 = @intCast(prepared.total_timer.untilNow(io, .awake).nanoseconds);

    return BatchResult{
        .total_files = work_items.len,
        .successful = successful,
        .failed = failed,
        .total_sasa_time_ns = total_sasa_time_ns,
        .total_time_ns = total_time_ns,
        .scan_time_ns = prepared.scan_time_ns,
        .build_items_time_ns = prepared.build_items_time_ns,
        .process_time_ns = process_time_ns,
        .read_parse_time_ns = total_read_parse_time_ns,
        .classifier_time_ns = total_classifier_time_ns,
        .jsonl_write_time_ns = total_jsonl_write_time_ns,
        .file_results = try results_list.toOwnedSlice(allocator),
        .allocator = allocator,
    };
}

/// A work item for batch processing.
/// Represents either a plain file or a specific molecule within an SDF file.
const WorkItem = struct {
    filename: []const u8, // Original filename in the input directory
    display_name: []const u8, // Name in results: the filename, or "stem_molname" for an SDF molecule
    /// Pre-parsed molecule of an SDF input (owned by `BatchWork`), null for other inputs.
    molecule: ?*const sdf_parser.SdfMolecule = null,
    /// 0-based position of `molecule` in its file.
    mol_idx: usize = 0,
    /// Set when an SDF input could not be read or parsed. The item then only
    /// reports this error, under the file name.
    load_error: ?anyerror = null,

    fn outputNameSource(self: WorkItem) OutputNameSource {
        return if (self.molecule != null) .sdf_molecule else .input_file;
    }
};

/// The work items of a batch run together with the parsed SDF molecules they
/// point to. SDF inputs are read and parsed once, when the items are built:
/// molecule names are needed to check the output names before anything is
/// calculated, and the runners then process those same molecules.
const BatchWork = struct {
    items: std.ArrayListUnmanaged(WorkItem) = .empty,
    sdf_molecules: std.ArrayListUnmanaged([]const sdf_parser.SdfMolecule) = .empty,

    fn deinit(self: *BatchWork, allocator: Allocator) void {
        for (self.items.items) |item| {
            // The display name of a plain file is the file name itself.
            if (item.display_name.ptr != item.filename.ptr) allocator.free(item.display_name);
        }
        self.items.deinit(allocator);
        for (self.sdf_molecules.items) |molecules| sdf_parser.freeMolecules(allocator, molecules);
        self.sdf_molecules.deinit(allocator);
    }
};

fn resolveBatchThreadCount(config_threads: usize, cpu_count: usize) usize {
    const auto_threads = @max(cpu_count, 1);
    return if (config_threads == 0) auto_threads else config_threads;
}

fn joinSpawnedThreads(threads: []std.Thread, spawned_count: usize) void {
    const count = @min(spawned_count, threads.len);
    for (threads[0..count]) |thread| {
        thread.join();
    }
}

/// Shared context for parallel workers
const ParallelContext = struct {
    /// Read-only for workers, including the SDF molecules the items point to.
    work_items: []const WorkItem,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    results: []FileResult,
    result_allocator: Allocator,
    next_item: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    jsonl_write_time_ns: std.atomic.Value(u64),
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    jsonl_stream: ?*JsonlStreamWriter,
    io: std.Io,
};

/// Process one work item and return its result under `result_name`.
/// Allocations follow `processOneFile`; SASA runs single-threaded.
fn processWorkItem(
    arena: Allocator,
    io: std.Io,
    result_allocator: Allocator,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    item: WorkItem,
    result_name: []const u8,
    config: BatchConfig,
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
) FileResult {
    var result = if (item.load_error) |err|
        FileResult{
            .filename = result_name,
            .n_atoms = 0,
            .sasa_time_ns = 0,
            .total_sasa = 0,
            .status = .err,
            .error_msg = std.fmt.allocPrint(result_allocator, "read/parse failed: {s}", .{@errorName(err)}) catch null,
        }
    else if (item.molecule) |molecule|
        processOneSdfMolecule(
            arena,
            io,
            result_allocator,
            item.display_name,
            molecule,
            output_dir,
            config,
            1, // single-threaded SASA per molecule
            lut_f64,
            lut_f32,
            coarse_lut_f64,
            fine_lut_f64,
            coarse_lut_f32,
            fine_lut_f32,
        )
    else
        processOneFile(
            arena,
            io,
            result_allocator,
            input_dir,
            output_dir,
            item.filename,
            config,
            1, // single-threaded SASA per file
            lut_f64,
            lut_f32,
            coarse_lut_f64,
            fine_lut_f64,
            coarse_lut_f32,
            fine_lut_f32,
        );
    result.filename = result_name;
    return result;
}

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

        var result = processWorkItem(
            arena.allocator(),
            ctx.io,
            ctx.result_allocator,
            ctx.input_dir,
            ctx.output_dir,
            work,
            name_copy,
            ctx.config,
            ctx.lut_f64,
            ctx.lut_f32,
            ctx.coarse_lut_f64,
            ctx.fine_lut_f64,
            ctx.coarse_lut_f32,
            ctx.fine_lut_f32,
        );

        // Store result (thread-safe: each index is unique)
        ctx.results[item_idx] = result;

        // Stream JSONL output (atom_areas on arena, valid until reset)
        if (ctx.jsonl_stream) |stream| {
            var write_timer: std.Io.Timestamp = undefined;
            if (ctx.config.profile_stages) write_timer = std.Io.Timestamp.now(ctx.io, .awake);
            stream.writeResult(arena.allocator(), &result);
            if (ctx.config.profile_stages) {
                const elapsed: u64 = @intCast(write_timer.untilNow(ctx.io, .awake).nanoseconds);
                _ = ctx.jsonl_write_time_ns.fetchAdd(elapsed, .monotonic);
            }
        }
        // Clear arena-owned payloads after streaming.
        ctx.results[item_idx].atom_areas = null;
        ctx.results[item_idx].residue_map = null;

        // Update progress counter (.release pairs with .acquire in progress monitor)
        _ = ctx.processed_count.fetchAdd(1, .release);

        // Reset arena for next work item
        _ = arena.reset(.retain_capacity);
    }
}

/// Build work items from file list, expanding SDF files into per-molecule items.
/// Caller must release the result with `BatchWork.deinit`.
fn buildWorkItems(
    allocator: Allocator,
    io: std.Io,
    files: []const []const u8,
    input_dir: []const u8,
) !BatchWork {
    var work = BatchWork{};
    errdefer work.deinit(allocator);

    for (files) |filename| {
        if (format_detect.detectInputFormat(filename) != .sdf) {
            try work.items.append(allocator, .{ .filename = filename, .display_name = filename });
            continue;
        }

        // Read and parse the SDF file. A file that cannot be loaded becomes a
        // single item that reports the error.
        const input_path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(input_path);

        const source = if (compressed.isCompressed(input_path))
            compressed.read(allocator, input_path) catch |err| {
                logWarning("{s}: failed to read SDF (compressed): {s}", .{ filename, @errorName(err) });
                try work.items.append(allocator, .{ .filename = filename, .display_name = filename, .load_error = err });
                continue;
            }
        else file_blk: {
            const f = std.Io.Dir.cwd().openFile(io, input_path, .{}) catch |err| {
                logWarning("{s}: failed to open SDF: {s}", .{ filename, @errorName(err) });
                try work.items.append(allocator, .{ .filename = filename, .display_name = filename, .load_error = err });
                continue;
            };
            defer f.close(io);
            var read_buf_build: [65536]u8 = undefined;
            var file_r_build = f.reader(io, &read_buf_build);
            break :file_blk file_r_build.interface.allocRemaining(allocator, .unlimited) catch |err| {
                logWarning("{s}: failed to read SDF: {s}", .{ filename, @errorName(err) });
                try work.items.append(allocator, .{ .filename = filename, .display_name = filename, .load_error = err });
                continue;
            };
        };
        // Parsed molecules do not refer to the source text.
        defer allocator.free(source);

        const molecules = sdf_parser.parse(allocator, source) catch |err| {
            logWarning("{s}: failed to parse SDF: {s}", .{ filename, @errorName(err) });
            try work.items.append(allocator, .{ .filename = filename, .display_name = filename, .load_error = err });
            continue;
        };
        work.sdf_molecules.append(allocator, molecules) catch |err| {
            sdf_parser.freeMolecules(allocator, molecules);
            return err;
        };

        // Create one work item per molecule
        try work.items.ensureUnusedCapacity(allocator, molecules.len);
        const display_names = try sdfMoleculeDisplayNames(allocator, filename, molecules);
        defer allocator.free(display_names);
        for (molecules, display_names, 0..) |*mol, display_name, mol_idx| {
            work.items.appendAssumeCapacity(.{
                .filename = filename,
                .display_name = display_name,
                .molecule = mol,
                .mol_idx = mol_idx,
            });
        }
    }

    return work;
}

/// The inputs of a batch run after scanning and validation.
const PreparedBatch = struct {
    files: [][]const u8,
    work: BatchWork,
    total_timer: std.Io.Timestamp,
    scan_time_ns: u64,
    build_items_time_ns: u64,

    fn deinit(self: *PreparedBatch, allocator: Allocator) void {
        self.work.deinit(allocator);
        freeScannedFiles(allocator, self.files);
    }
};

/// Scan `input_dir`, expand the inputs into work items and validate them.
///
/// Both runners start here. Nothing is calculated and the output directory is
/// not created yet, so a batch that is rejected leaves nothing behind.
fn prepareBatch(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
) !PreparedBatch {
    // Start total timer
    const total_timer = std.Io.Timestamp.now(io, .awake);

    // Scan directory for files
    const scan_timer = std.Io.Timestamp.now(io, .awake);
    const files = try scanDirectory(allocator, io, input_dir);
    errdefer freeScannedFiles(allocator, files);
    const scan_time_ns: u64 = @intCast(scan_timer.untilNow(io, .awake).nanoseconds);
    try validateChainMapInputFormats(files, config);

    // Build work items (expanding SDF files into per-molecule items)
    const build_timer = std.Io.Timestamp.now(io, .awake);
    var work = try buildWorkItems(allocator, io, files, input_dir);
    errdefer work.deinit(allocator);
    const build_items_time_ns: u64 = @intCast(build_timer.untilNow(io, .awake).nanoseconds);
    try validateUniqueOutputNames(allocator, work.items.items, output_dir, config);

    return .{
        .files = files,
        .work = work,
        .total_timer = total_timer,
        .scan_time_ns = scan_time_ns,
        .build_items_time_ns = build_items_time_ns,
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
    try validateBatchOutputFormat(config.output_format);

    var prepared = try prepareBatch(allocator, io, input_dir, output_dir, config);
    defer prepared.deinit(allocator);

    return runPrepared(allocator, io, input_dir, output_dir, config, jsonl_output_path, &prepared);
}

/// Process the work items of a prepared batch: in parallel, N items at a time
/// with one SASA thread each, or one after another when there is nothing to
/// run in parallel.
fn runPrepared(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
    prepared: *const PreparedBatch,
) !BatchResult {
    const work_items = prepared.work.items.items;

    // Determine thread count
    const cpu_count = std.Thread.getCpuCount() catch 1;
    const n_threads = resolveBatchThreadCount(config.n_threads, cpu_count);
    const actual_threads = @min(n_threads, work_items.len);

    // Nothing to run in parallel for at most one item or a single thread. An
    // empty directory goes this way too, so that it leaves the same output
    // behind whatever the thread count: an empty JSONL file, not a stale one.
    if (work_items.len <= 1 or n_threads <= 1) {
        return runPreparedSequential(allocator, io, input_dir, output_dir, config, jsonl_output_path, prepared);
    }

    // Create output directory if specified
    if (output_dir) |out_dir| {
        try std.Io.Dir.cwd().createDirPath(io, out_dir);
    }

    // Allocate results (one per work item)
    const file_results = try allocator.alloc(FileResult, work_items.len);
    errdefer allocator.free(file_results);

    // Build bitmask LUT once (if enabled)
    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    // Open JSONL output file (or stdout) when requested.
    var jsonl_file: ?std.Io.File = null;
    var jsonl_file_needs_close = false;
    if (jsonl_output_path) |path| {
        jsonl_file = try std.Io.Dir.cwd().createFile(io, path, .{});
        jsonl_file_needs_close = true;
    } else if (batchWritesJsonl(config)) {
        jsonl_file = std.Io.File.stdout();
    }
    defer if (jsonl_file_needs_close) {
        if (jsonl_file) |f| f.close(io);
    };

    // Set up the stream writer on the stack (if JSONL streaming is active).
    // SAFETY: `undefined` when jsonl_file is null — never accessed because
    // jsonl_stream_ptr is also null in that case.
    var jsonl_stream_buffer: [64 * 1024]u8 = undefined;
    var jsonl_stream_storage: JsonlStreamWriter = if (jsonl_file) |jf|
        JsonlStreamWriter.init(jf, io, jsonlOptions(config), &jsonl_stream_buffer)
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
        .jsonl_write_time_ns = std.atomic.Value(u64).init(0),
        .lut_f64 = luts.f64Ptr(),
        .lut_f32 = luts.f32Ptr(),
        .coarse_lut_f64 = luts.coarseF64Ptr(),
        .fine_lut_f64 = luts.fineF64Ptr(),
        .coarse_lut_f32 = luts.coarseF32Ptr(),
        .fine_lut_f32 = luts.fineF32Ptr(),
        .jsonl_stream = jsonl_stream_ptr,
        .io = io,
    };

    // Spawn worker threads
    const threads = try allocator.alloc(std.Thread, actual_threads);
    defer allocator.free(threads);

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing items", .estimated_total_items = work_items.len })
    else
        .none;
    defer progress_root.end();

    var process_timer = std.Io.Timestamp.now(io, .awake);
    {
        // Scoped so a later error cannot join the same threads again.
        var spawned_count: usize = 0;
        errdefer joinSpawnedThreads(threads, spawned_count);
        for (threads) |*thread| {
            thread.* = try std.Thread.spawn(.{}, parallelWorker, .{&ctx});
            spawned_count += 1;
        }

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
        joinSpawnedThreads(threads, spawned_count);
    }
    const process_time_ns: u64 = @intCast(process_timer.untilNow(io, .awake).nanoseconds);

    // Aggregate results
    var total_sasa_time_ns: u64 = 0;
    var total_read_parse_time_ns: u64 = 0;
    var total_classifier_time_ns: u64 = 0;
    var successful: usize = 0;
    var failed: usize = 0;

    for (file_results) |result| {
        if (result.status == .ok) {
            successful += 1;
            total_sasa_time_ns += result.sasa_time_ns;
            total_read_parse_time_ns += result.read_parse_time_ns;
            total_classifier_time_ns += result.classifier_time_ns;
        } else {
            failed += 1;
        }
    }

    if (jsonl_stream_ptr) |stream| {
        try stream.flush();
        if (stream.hasError()) return error.JsonlWriteFailed;
    }

    const total_time_ns: u64 = @intCast(prepared.total_timer.untilNow(io, .awake).nanoseconds);

    return BatchResult{
        .total_files = work_items.len,
        .successful = successful,
        .failed = failed,
        .total_sasa_time_ns = total_sasa_time_ns,
        .total_time_ns = total_time_ns,
        .scan_time_ns = prepared.scan_time_ns,
        .build_items_time_ns = prepared.build_items_time_ns,
        .process_time_ns = process_time_ns,
        .read_parse_time_ns = total_read_parse_time_ns,
        .classifier_time_ns = total_classifier_time_ns,
        .jsonl_write_time_ns = ctx.jsonl_write_time_ns.load(.monotonic),
        .file_results = file_results,
        .allocator = allocator,
    };
}

/// Run batch processing (main entry point)
/// Uses file-level parallelism: N files in parallel, 1 thread per file; with
/// one thread or at most one input the files are processed one after another.
pub fn runBatch(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
) !BatchResult {
    return runBatchParallel(allocator, io, input_dir, output_dir, config, jsonl_output_path);
}

/// The step of `runBatchReportingStage` that was running when it failed.
pub const BatchStage = enum {
    /// Reading the input directory and validating the inputs.
    scan_inputs,
    /// Creating the output directory.
    create_output_dir,
    /// Processing the files.
    process,
};

/// `runBatch`, taking the steps one at a time so that a caller can tell
/// which one failed: the same filesystem error means a bad input directory
/// while scanning and a bad output directory while creating it.
///
/// `stage` is set before each step; after an error it names the step that
/// returned it. As in the workflow runner, the output directory is created
/// only after the inputs are accepted.
pub fn runBatchReportingStage(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
    stage: *BatchStage,
) !BatchResult {
    stage.* = .scan_inputs;
    try validateBatchOutputFormat(config.output_format);

    var prepared = try prepareBatch(allocator, io, input_dir, output_dir, config);
    defer prepared.deinit(allocator);

    if (output_dir) |dir| {
        stage.* = .create_output_dir;
        try std.Io.Dir.cwd().createDirPath(io, dir);
    }

    stage.* = .process;
    return runPrepared(allocator, io, input_dir, output_dir, config, null, &prepared);
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
    alt_loc_mode: mmcif_parser.AltLocMode = .auto,
    alt_loc_id: u8 = 'A',
    af_model_fast: bool = false,
    residue_map: bool = false,
    jsonl_decimals: ?u8 = null,
    n_threads: usize = 0,
    probe_radius: f64 = 1.4,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    lr_trig: TrigMode = .exact,
    algorithm: Algorithm = .sr,
    precision: Precision = .f64,
    output_format: OutputFormat = .json,
    classifier_type: ClassifierType = .ccd, // Default: ccd (ProtOr-compatible with CCD extension)
    include_hydrogens: bool = false,
    include_hetatm: bool = false,
    use_bitmask: bool = false,
    bitmask_correction: bool = false,
    bitmask_correction_coeff: f64 = shrake_rupley_bitmask.default_bitmask_correction_coeff,
    adaptive_sr: bool = false,
    coarse_points: u32 = 64,
    fine_points: u32 = 256,
    adaptive_low: f64 = 0.10,
    adaptive_high: f64 = 0.90,
    ccd_path: ?[]const u8 = null, // External CCD dictionary file (.zsdc or .cif[.gz|.zst])
    sdf_paths: SdfPathList = .{}, // --sdf=PATH (up to 16)
    quiet: bool = false,
    show_progress: bool = true,
    show_timing: bool = false,
    profile_stages: bool = false,
    input_io: InputIoMode = .auto,
    show_help: bool = false,
    threads_explicit: bool = false,
    probe_radius_explicit: bool = false,
    n_points_explicit: bool = false,
    n_slices_explicit: bool = false,
    lr_trig_explicit: bool = false,
    algorithm_explicit: bool = false,
    precision_explicit: bool = false,
    format_explicit: bool = false,
    classifier_explicit: bool = false,
    include_hydrogens_explicit: bool = false,
    include_hetatm_explicit: bool = false,
    alt_loc_explicit: bool = false,
    af_model_fast_explicit: bool = false,
    use_bitmask_explicit: bool = false,
    bitmask_correction_explicit: bool = false,
    bitmask_correction_coeff_explicit: bool = false,
    adaptive_sr_explicit: bool = false,
    coarse_points_explicit: bool = false,
    fine_points_explicit: bool = false,
    adaptive_low_explicit: bool = false,
    adaptive_high_explicit: bool = false,
    ccd_explicit: bool = false,
    sdf_explicit: bool = false,
    quiet_explicit: bool = false,
    timing_explicit: bool = false,
    profile_stages_explicit: bool = false,
    input_io_explicit: bool = false,
    jsonl_decimals_explicit: bool = false,
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

fn validateBatchOutputFormat(output_format: OutputFormat) !void {
    switch (output_format) {
        .json, .compact, .csv, .jsonl => {},
        .freesasa, .rsa => return error.InvalidArgument,
    }
}

fn validateChainMapInputFormats(files: []const []const u8, config: BatchConfig) !void {
    return validateMappedChainInputFormats(files, config.chain_map != null);
}

fn validateMappedChainInputFormats(files: []const []const u8, enabled: bool) !void {
    if (!enabled) return;
    for (files) |filename| {
        switch (format_detect.detectInputFormat(filename)) {
            .pdb, .mmcif, .bcif => {},
            .json, .sdf => return error.UnsupportedChainMapInputFormat,
        }
    }
}

/// A per-file output together with the input it would be written for.
const OutputNameClaim = struct {
    output_name: []const u8,
    filename: []const u8,
    /// 1-based position of the molecule in `filename` for an SDF molecule
    /// output, null for the single output of any other input.
    molecule: ?usize = null,

    fn lessThan(_: void, a: OutputNameClaim, b: OutputNameClaim) bool {
        const orders = [_]std.math.Order{
            std.ascii.orderIgnoreCase(a.output_name, b.output_name),
            std.mem.order(u8, a.output_name, b.output_name),
            std.mem.order(u8, a.filename, b.filename),
        };
        for (orders) |order| {
            if (order != .eq) return order == .lt;
        }
        return (a.molecule orelse 0) < (b.molecule orelse 0);
    }

    /// Output names are compared without regard to ASCII case: names that
    /// differ only in case are one file on a case-insensitive filesystem (the
    /// macOS and Windows defaults). On a case-sensitive filesystem this
    /// rejects inputs that could be written side by side, deliberately, so
    /// that a batch gives the same answer everywhere.
    fn sharesOutputWith(a: OutputNameClaim, b: OutputNameClaim) bool {
        return std.ascii.eqlIgnoreCase(a.output_name, b.output_name);
    }
};

/// Find the work items whose per-file output name is shared with another item.
///
/// The names compared are the ones the runners write (`perFileOutputName`):
/// `1crn.pdb` and `1crn.cif.gz` both map to `1crn.json`, and an SDF molecule
/// output such as `lig_1.json` can equal the output of `lig_1.pdb` or of a
/// molecule in another SDF file. Items that fail before writing an output
/// (no chain map entry, unreadable SDF input) claim nothing.
///
/// Returns the colliding claims sorted so that inputs sharing an output name
/// are adjacent. All returned memory is allocated from `arena`.
fn findOutputNameCollisions(
    arena: Allocator,
    items: []const WorkItem,
    config: BatchConfig,
) ![]const OutputNameClaim {
    const ext = getOutputExtension(config.output_format);

    var claims = try std.ArrayListUnmanaged(OutputNameClaim).initCapacity(arena, items.len);
    for (items) |item| {
        if (item.load_error != null) continue;
        // Files without a chain map entry fail before any output is written.
        if (config.chain_map) |map| {
            if (map.get(item.filename) == null) continue;
        }

        claims.appendAssumeCapacity(.{
            .output_name = try perFileOutputName(arena, item.outputNameSource(), item.display_name, ext),
            .filename = item.filename,
            .molecule = if (item.molecule != null) item.mol_idx + 1 else null,
        });
    }

    std.mem.sort(OutputNameClaim, claims.items, {}, OutputNameClaim.lessThan);

    var collisions = std.ArrayListUnmanaged(OutputNameClaim).empty;
    for (claims.items, 0..) |claim, i| {
        const shares_prev = i > 0 and claim.sharesOutputWith(claims.items[i - 1]);
        const shares_next = i + 1 < claims.items.len and claim.sharesOutputWith(claims.items[i + 1]);
        if (shares_prev or shares_next) try collisions.append(arena, claim);
    }
    return collisions.items;
}

/// Reject work items whose per-file outputs would overwrite each other.
///
/// Only applies when one output file is written per item; JSONL output
/// keeps one record per item and cannot collide.
fn validateUniqueOutputNames(
    allocator: Allocator,
    items: []const WorkItem,
    output_dir: ?[]const u8,
    config: BatchConfig,
) !void {
    if (output_dir == null or batchWritesJsonl(config)) return;

    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();

    const collisions = try findOutputNameCollisions(arena.allocator(), items, config);
    if (collisions.len == 0) return;

    std.debug.print("{s}", .{try formatOutputNameCollisions(arena.allocator(), collisions)});
    return error.OutputNameCollision;
}

/// Describe `collisions` (as returned by `findOutputNameCollisions`) for the
/// user: one line per shared output name, listing the inputs that claim it.
fn formatOutputNameCollisions(allocator: Allocator, collisions: []const OutputNameClaim) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    defer aw.deinit();
    const w = &aw.writer;

    var n_names: usize = 0;
    for (collisions, 0..) |claim, i| {
        if (i == 0 or !claim.sharesOutputWith(collisions[i - 1])) n_names += 1;
    }
    try w.print("Error: {d} output name{s} shared by more than one input:", .{
        n_names,
        if (n_names == 1) " is" else "s are",
    });

    var group_start: usize = 0;
    while (group_start < collisions.len) {
        var group_end = group_start + 1;
        while (group_end < collisions.len and collisions[group_end].sharesOutputWith(collisions[group_start])) {
            group_end += 1;
        }
        const group = collisions[group_start..group_end];
        group_start = group_end;

        // Claims are sorted, so equal spellings of the name are adjacent.
        var n_spellings: usize = 1;
        for (group[1..], group[0 .. group.len - 1]) |claim, prev| {
            if (!std.mem.eql(u8, claim.output_name, prev.output_name)) n_spellings += 1;
        }

        try w.print("\n  {s} <- ", .{group[0].output_name});
        for (group, 0..) |claim, i| {
            if (i > 0) try w.writeAll(", ");
            try w.writeAll(claim.filename);
            if (claim.molecule) |n| try w.print(" (molecule {d})", .{n});
        }
        if (n_spellings == group.len) {
            try w.writeAll(" (the output names differ only in case)");
        } else if (n_spellings > 1) {
            try w.writeAll(" (some of the output names differ only in case)");
        }
    }
    try w.writeAll(
        "\nSplit these inputs into separate directories, or use JSONL output (--format=jsonl) to keep one record per input.\n",
    );
    return aw.toOwnedSlice();
}

/// `validateUniqueOutputNames` for a runner that writes one output per input
/// file whatever its format (the file-first workflow runner, which does not
/// expand SDF files into molecules).
fn validateUniqueFileOutputNames(
    allocator: Allocator,
    files: []const []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
) !void {
    if (output_dir == null or batchWritesJsonl(config)) return;

    const items = try plainWorkItems(allocator, files);
    defer allocator.free(items);
    return validateUniqueOutputNames(allocator, items, output_dir, config);
}

/// One work item per file, without expanding SDF files. Caller frees the slice.
fn plainWorkItems(allocator: Allocator, files: []const []const u8) ![]WorkItem {
    const items = try allocator.alloc(WorkItem, files.len);
    for (files, items) |filename, *item| {
        item.* = .{ .filename = filename, .display_name = filename };
    }
    return items;
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

fn parseBitmaskCorrectionCoeff(value: []const u8) f64 {
    const coeff = std.fmt.parseFloat(f64, value) catch {
        std.debug.print("Error: Invalid bitmask correction coefficient: {s}\n", .{value});
        std.process.exit(1);
    };
    if (!std.math.isFinite(coeff) or coeff < 0.0) {
        std.debug.print("Error: Bitmask correction coefficient must be finite and non-negative: {d}\n", .{coeff});
        std.process.exit(1);
    }
    return coeff;
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

/// Parse and validate lr-trig value (for Lee-Richards)
fn parseLrTrig(value: []const u8) TrigMode {
    return TrigMode.fromString(value) orelse {
        std.debug.print("Error: Invalid lr-trig: {s}\n", .{value});
        std.debug.print("Valid values: exact, fast\n", .{});
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

fn parseAltLocSetting(value: []const u8) mmcif_parser.AltLocSetting {
    return mmcif_parser.parseAltLocSetting(value) orelse {
        std.debug.print("Error: Invalid altloc mode: {s}\n", .{value});
        std.debug.print("Valid altloc modes: auto, none, all, highest-occupancy, or a single altLoc ID like A\n", .{});
        std.process.exit(1);
    };
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

fn parseInputIoMode(value: []const u8) InputIoMode {
    return InputIoMode.parse(value) orelse {
        std.debug.print("Error: Invalid input I/O mode: {s}\n", .{value});
        std.debug.print("Valid input I/O modes: auto, mmap, read\n", .{});
        std.process.exit(1);
    };
}

fn parseJsonlDecimals(value: []const u8) u8 {
    const decimals = std.fmt.parseInt(u8, value, 10) catch {
        std.debug.print("Error: Invalid jsonl decimals: {s}\n", .{value});
        std.process.exit(1);
    };
    if (decimals > 15) {
        std.debug.print("Error: --jsonl-decimals must be between 0 and 15: {d}\n", .{decimals});
        std.process.exit(1);
    }
    return decimals;
}

/// Split a `--chain` value such as "A" or "A,B" into chain IDs. A value
/// without any chain ID ("", "," or " , ") is `error.EmptyChainFilter`: it
/// would select nothing and fail every input.
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
    if (chains.items.len == 0) return error.EmptyChainFilter;

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
        // --lr-trig=MODE or --lr-trig MODE (for Lee-Richards)
        else if (std.mem.startsWith(u8, arg, "--lr-trig=")) {
            result.lr_trig_explicit = true;
            const value = arg["--lr-trig=".len..];
            result.lr_trig = parseLrTrig(value);
        } else if (std.mem.eql(u8, arg, "--lr-trig")) {
            result.lr_trig_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --lr-trig\n", .{});
                std.process.exit(1);
            }
            result.lr_trig = parseLrTrig(args[i]);
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
        // --jsonl-decimals=N or --jsonl-decimals N
        else if (std.mem.startsWith(u8, arg, "--jsonl-decimals=")) {
            result.jsonl_decimals_explicit = true;
            const value = arg["--jsonl-decimals=".len..];
            result.jsonl_decimals = parseJsonlDecimals(value);
        } else if (std.mem.eql(u8, arg, "--jsonl-decimals")) {
            result.jsonl_decimals_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --jsonl-decimals\n", .{});
                std.process.exit(1);
            }
            result.jsonl_decimals = parseJsonlDecimals(args[i]);
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
        // --input-io=MODE or --input-io MODE
        else if (std.mem.startsWith(u8, arg, "--input-io=")) {
            result.input_io_explicit = true;
            const value = arg["--input-io=".len..];
            result.input_io = parseInputIoMode(value);
        } else if (std.mem.eql(u8, arg, "--input-io")) {
            result.input_io_explicit = true;
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --input-io\n", .{});
                std.process.exit(1);
            }
            result.input_io = parseInputIoMode(args[i]);
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
        // --af-model-fast
        else if (std.mem.eql(u8, arg, "--af-model-fast")) {
            result.af_model_fast = true;
            result.af_model_fast_explicit = true;
        }
        // --altloc=MODE or --altloc MODE
        else if (std.mem.startsWith(u8, arg, "--altloc=")) {
            const setting = parseAltLocSetting(arg["--altloc=".len..]);
            result.alt_loc_mode = setting.mode;
            result.alt_loc_id = setting.id;
            result.alt_loc_explicit = true;
        } else if (std.mem.eql(u8, arg, "--altloc")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --altloc\n", .{});
                std.process.exit(1);
            }
            const setting = parseAltLocSetting(args[i]);
            result.alt_loc_mode = setting.mode;
            result.alt_loc_id = setting.id;
            result.alt_loc_explicit = true;
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
        // --bitmask-correction
        else if (std.mem.eql(u8, arg, "--bitmask-correction")) {
            result.bitmask_correction = true;
            result.bitmask_correction_explicit = true;
        }
        // --bitmask-correction-coeff=VALUE or --bitmask-correction-coeff VALUE
        else if (std.mem.startsWith(u8, arg, "--bitmask-correction-coeff=")) {
            result.bitmask_correction_coeff = parseBitmaskCorrectionCoeff(arg["--bitmask-correction-coeff=".len..]);
            result.bitmask_correction = true;
            result.bitmask_correction_explicit = true;
            result.bitmask_correction_coeff_explicit = true;
        } else if (std.mem.eql(u8, arg, "--bitmask-correction-coeff")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --bitmask-correction-coeff\n", .{});
                std.process.exit(1);
            }
            result.bitmask_correction_coeff = parseBitmaskCorrectionCoeff(args[i]);
            result.bitmask_correction = true;
            result.bitmask_correction_explicit = true;
            result.bitmask_correction_coeff_explicit = true;
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
        // --profile-stages
        else if (std.mem.eql(u8, arg, "--profile-stages")) {
            result.profile_stages = true;
            result.profile_stages_explicit = true;
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
        \\    <input_dir>     Directory containing structure files
        \\                    Supported: .json, .cif, .mmcif, .bcif, .pdb, .ent,
        \\                    .sdf, .mol, each also as .gz or .zst
        \\    [output_dir]    Optional output directory (default: no file output)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO    Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                        Default: sr
        \\    --classifier=TYPE   Built-in classifier: ccd, protor, naccess, oons
        \\                        Default: ccd (protor uses static ProtOr radii only)
        \\    --ccd=PATH          External CCD dictionary file (.zsdc or .cif[.gz|.zst])
        \\                        Used with --classifier=ccd for non-standard residues
        \\    --sdf=PATH          SDF file with bond topology for CCD classifier
        \\                        Can be specified multiple times for multiple ligands
        \\    --threads=N         Number of batch workers (default: auto-detect)
        \\                        Explicit values may exceed CPU count for I/O-bound inputs
        \\    --workflow=PATH     TOML workflow file with one or more named batch jobs
        \\    --manifest=PATH     Compatibility alias for --workflow
        \\    --chain=ID          Filter by chain ID for non-workflow batch (e.g. A or A,B)
        \\    --auth-chain        Use auth_asym_id instead of label_asym_id for mmCIF chain matching,
        \\                        and auth_seq_id instead of label_seq_id for residue numbers
        \\    --af-model-fast     Use an experimental AlphaFold-model mmCIF fast parser
        \\                        with multi-chain metadata and safe generic fallback
        \\    --altloc=MODE       Alternate-location handling (PDB/mmCIF/BCIF): auto, none,
        \\                        all, highest-occupancy, or a single ID like A
        \\                        (default: auto)
        \\    --residue-map       Include compact residue map arrays in JSONL output
        \\    --probe-radius=R    Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N        Test points per atom (default: 100, for sr)
        \\    --n-slices=N        Slices per atom diameter (default: 20, for lr)
        \\    --lr-trig=MODE      Arc angles for lr: exact (acos/atan2, default) or
        \\                        fast (polynomial approximation, the results of
        \\                        zsasa 0.9.1 and earlier; totals come out a few
        \\                        tenths of a percent too high)
        \\    --precision=PREC    Floating-point precision: f32, f64 (default: f64)
        \\    --input-io=MODE     File input strategy where supported: auto, mmap, read
        \\                        Default: auto (AF fast uses read; others keep defaults)
        \\    --format=FORMAT     Output format: json, compact, csv, jsonl (default: json)
        \\    --jsonl-decimals=N  Round JSONL floating-point values to N decimals (0..15)
        \\    --include-hydrogens Include hydrogen atoms (default: exclude)
        \\    --include-hetatm    Include HETATM records (default: exclude)
        \\    --use-bitmask       Use bitmask LUT optimization for SR algorithm
        \\                        (n-points must be 1..1024)
        \\    --bitmask-correction
        \\                        Experimental correction for bitmask quantization bias
        \\                        Requires --use-bitmask
        \\    --bitmask-correction-coeff=V
        \\                        Override correction coefficient (default: 0.020)
        \\    --adaptive-sr       Experimental adaptive two-stage bitmask SR
        \\    --coarse-points=N   Coarse adaptive points (default: 64)
        \\    --fine-points=N     Fine adaptive points (default: --n-points or 256)
        \\    --adaptive-low=X    Coarse accept low exposed fraction (default: 0.10)
        \\    --adaptive-high=X   Coarse accept high exposed fraction (default: 0.90)
        \\    --timing            Show timing breakdown for benchmarking
        \\    --profile-stages    Include read/parse, classifier, and JSONL stage timings
        \\    -o, --output=PATH   Output directory, or file path for --format=jsonl
        \\    -q, --quiet         Suppress progress output and the summary; inputs
        \\                        that fail are still reported on stderr
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
    if (!args.lr_trig_explicit) {
        if (calculation.lr_trig) |v| config.lr_trig = TrigMode.fromString(v) orelse {
            std.debug.print("Error: workflow lr_trig must be \"exact\" or \"fast\": {s}\n", .{v});
            return error.InvalidArgument;
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
    if (!args.jsonl_decimals_explicit) {
        if (output.jsonl.decimals) |v| config.jsonl_decimals = v;
    }
    if (output.jsonl.atom_areas) |v| config.jsonl_include_atom_areas = v;
    if (output.jsonl.atom_identity) |v| config.jsonl_include_atom_identity = v;
    if (output.jsonl.total_area) |v| config.jsonl_include_total_area = v;
    if (output.jsonl.metadata) |v| config.jsonl_metadata = parseWorkflowJsonlMetadata(v) catch |err| {
        std.debug.print("Error: output.jsonl.metadata must be \"none\" or \"sidecar\": {s}\n", .{v});
        return err;
    };
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
    if (!args.alt_loc_explicit) {
        if (calculation.altloc) |v| {
            config.alt_loc_mode = v.mode;
            config.alt_loc_id = v.id;
        }
    }
    if (args.af_model_fast_explicit) config.af_model_fast = args.af_model_fast;
    if (calculation.residue_map) |v| config.residue_map = v;

    try applyWorkflowClassifierToBatchConfig(config, args, classifier_config);
}

fn parseWorkflowJsonlMetadata(value: []const u8) !JsonlMetadataMode {
    if (std.mem.eql(u8, value, "none")) return .none;
    if (std.mem.eql(u8, value, "sidecar")) return .sidecar;
    return error.InvalidArgument;
}

fn applyCliOverrides(config: *BatchConfig, args: BatchArgs) void {
    if (args.threads_explicit) config.n_threads = args.n_threads;
    if (args.algorithm_explicit) config.algorithm = args.algorithm;
    if (args.n_points_explicit) config.n_points = args.n_points;
    if (args.n_slices_explicit) config.n_slices = args.n_slices;
    if (args.lr_trig_explicit) config.lr_trig = args.lr_trig;
    if (args.probe_radius_explicit) config.probe_radius = args.probe_radius;
    if (args.precision_explicit) config.precision = args.precision;
    if (args.format_explicit) config.output_format = args.output_format;
    if (args.timing_explicit) config.show_timing = args.show_timing;
    if (args.profile_stages_explicit) config.profile_stages = args.profile_stages;
    if (args.input_io_explicit) config.input_io = args.input_io;
    if (args.af_model_fast_explicit) config.af_model_fast = args.af_model_fast;
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
    if (args.bitmask_correction_explicit) config.bitmask_correction = args.bitmask_correction;
    if (args.bitmask_correction_coeff_explicit) config.bitmask_correction_coeff = args.bitmask_correction_coeff;
    if (args.use_auth_chain) config.use_auth_chain = true;
    if (args.alt_loc_explicit) {
        config.alt_loc_mode = args.alt_loc_mode;
        config.alt_loc_id = args.alt_loc_id;
    }
    if (args.residue_map) config.residue_map = true;
    if (args.jsonl_decimals_explicit) config.jsonl_decimals = args.jsonl_decimals;
}

fn validateBitmaskCorrectionConfig(config: BatchConfig) !void {
    if (config.bitmask_correction and !config.use_bitmask) {
        std.debug.print("Error: --bitmask-correction requires --use-bitmask\n", .{});
        return error.InvalidArgument;
    }
    if (config.bitmask_correction and config.adaptive_sr) {
        std.debug.print("Error: --bitmask-correction is not supported with --adaptive-sr\n", .{});
        return error.InvalidArgument;
    }
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

/// The keys and options that a batch workflow does not read: those of the
/// manifest (see `workflow_manifest.checkKeys`) and the command-line options
/// that only have an effect outside a workflow.
fn batchWorkflowFindings(workflow: workflow_manifest.Workflow, args: BatchArgs) workflow_manifest.Findings {
    const is_analysis = workflow.analysis != null;
    var findings = workflow_manifest.checkKeys(workflow, if (is_analysis) .batch_analysis else .batch_jobs);

    // An [analysis] workflow reports the SASA time; a workflow with jobs
    // reports through the per-job summary only.
    if (!is_analysis) {
        if (args.timing_explicit and args.show_timing) {
            findings.add(.warning, "--timing has no effect on a workflow with [[jobs]] " ++
                "(only an [analysis] workflow reports timing): remove the option");
        } else if (workflow.calculation.timing orelse false) {
            findings.add(.warning, "[calculation] timing = true has no effect on a workflow with [[jobs]] " ++
                "(only an [analysis] workflow reports timing): remove the key");
        }
    }
    if (args.profile_stages_explicit and args.profile_stages) {
        findings.add(.warning, "--profile-stages has no effect on a workflow " ++
            "(stage timings are reported only by batch without --workflow): remove the option");
    }
    if (is_analysis) {
        if (args.residue_map) {
            findings.add(.err, "--residue-map is not supported by an [analysis] workflow " ++
                "(its rows carry no residue map): remove the option");
        }
        if (args.format_explicit and args.output_format != .jsonl) {
            findings.add(.err, "--format other than jsonl is not supported by an [analysis] workflow " ++
                "(it always writes JSONL): remove the option");
        }
    }
    return findings;
}

fn parseWorkflowFile(allocator: Allocator, io: std.Io, path: []const u8) !workflow_manifest.Workflow {
    return workflow_manifest.parseFile(allocator, io, path);
}

fn printWorkflowReadError(path: []const u8, err: anyerror) void {
    if (builtin.is_test) return;
    std.debug.print("Error reading workflow file '{s}': {s}\n", .{ path, @errorName(err) });
    if (workflow_manifest.errorHint(err)) |hint| std.debug.print("  {s}\n", .{hint});
}

fn classifierUsesCcdResources(effective_classifier_type: ?ClassifierType) bool {
    const classifier_type = effective_classifier_type orelse return false;
    return classifier_type == .ccd;
}

fn batchArgsUseCcdResources(args: BatchArgs) bool {
    return classifierUsesCcdResources(args.classifier_type);
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
    var diag: classifier_parser.Diagnostic = .{};
    return classifier_parser.parseConfigFileDiag(allocator, io, path, &diag) catch |err| {
        switch (err) {
            error.UnsupportedConfigExtension => std.debug.print("Error loading config file '{s}': custom classifier configs are TOML-only; rename or convert the file to .toml\n", .{path}),
            error.UnsupportedLegacyFormat => std.debug.print("Error loading config file '{s}': FreeSASA-style custom classifier configs are no longer supported; convert to TOML [types] and [[atoms]]\n", .{path}),
            else => if (diag.line != 0)
                std.debug.print("Error loading config file '{s}' (line {d}): {s}\n", .{ path, diag.line, @errorName(err) })
            else
                std.debug.print("Error loading config file '{s}': {s}\n", .{ path, @errorName(err) }),
        }
        return err;
    };
}

fn runWorkflow(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (args.chain_filter != null) return runWorkflowJobFirst(allocator, io, args);

    const workflow_path = args.workflow_path orelse return runWorkflowFileFirst(allocator, io, args, null);
    var workflow = parseWorkflowFile(allocator, io, workflow_path) catch {
        return runWorkflowFileFirst(allocator, io, args, null);
    };
    defer workflow.deinit();

    if (workflow.analysis != null) return runWorkflowBsaAnalysis(allocator, io, args, workflow);

    if (workflow.jobs.len == 0) return runWorkflowFileFirst(allocator, io, args, null);
    if (workflowHasChainMap(workflow)) return runWorkflowJobFirst(allocator, io, args);
    const input_dir = args.input_path orelse workflow.input.dir orelse return runWorkflowFileFirst(allocator, io, args, null);
    const output_dir = args.output_path orelse workflow.output.dir;
    if (workflow.jobs.len > 1 and output_dir == null) return runWorkflowFileFirst(allocator, io, args, null);

    var resource_config = BatchConfig{};
    try applyWorkflowToBatchConfig(&resource_config, args, workflow.calculation, workflow.output, workflow.classifier);
    applyCliOverrides(&resource_config, args);
    try validateBitmaskCorrectionConfig(resource_config);
    if (workflowRequiresJobFirstForAuthChain(args, workflow, resource_config)) {
        return runWorkflowJobFirst(allocator, io, args);
    }

    const files = try scanDirectory(allocator, io, input_dir);
    defer freeScannedFiles(allocator, files);

    if (workflowFilesContainSdf(files)) {
        return runWorkflowJobFirst(allocator, io, args);
    }
    if (workflowRequiresJobFirstForLongChainFormats(files, workflow)) {
        return runWorkflowJobFirst(allocator, io, args);
    }
    return runWorkflowFileFirst(allocator, io, args, files);
}

fn workflowHasChainMap(workflow: workflow_manifest.Workflow) bool {
    for (workflow.jobs) |job| {
        if (job.chain_map != null) return true;
    }
    return false;
}

fn freeScannedFiles(allocator: Allocator, files: []const []const u8) void {
    for (files) |f| allocator.free(f);
    allocator.free(files);
}

fn workflowFilesContainSdf(files: []const []const u8) bool {
    for (files) |filename| {
        if (format_detect.detectInputFormat(filename) == .sdf) return true;
    }
    return false;
}

fn workflowRequiresJobFirstForLongChainFormats(files: []const []const u8, workflow: workflow_manifest.Workflow) bool {
    _ = files;
    _ = workflow;
    return false;
}

fn workflowRequiresJobFirstForAuthChain(args: BatchArgs, workflow: workflow_manifest.Workflow, shared_config: BatchConfig) bool {
    for (workflow.jobs) |job| {
        var job_use_auth_chain = shared_config.use_auth_chain;
        if (job.auth_chain) |v| job_use_auth_chain = v;
        if (args.use_auth_chain) job_use_auth_chain = true;
        if (job_use_auth_chain != shared_config.use_auth_chain) return true;
    }
    return false;
}

const BsaInterfaceSelection = struct {
    id: []const u8,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
};

const BsaInterfaceStats = struct {
    successful: bool,
    sasa_time_ns: u64 = 0,
    /// The error row of a failed interface, for the failure report. The
    /// strings live as long as the allocator the row was written with.
    failure: ?struct { filename: []const u8, id: []const u8, reason: []const u8 } = null,
};

const BsaWorkflowCounter = struct {
    successful: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    failed: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    total_sasa_time_ns: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),
};

fn bsaPartnersOverlap(partner_a: []const []const u8, partner_b: []const []const u8) bool {
    for (partner_a) |a| {
        for (partner_b) |b| {
            if (std.mem.eql(u8, a, b)) return true;
        }
    }
    return false;
}

fn bsaInputContainsChain(input: AtomInput, target: []const u8) bool {
    for (0..input.atomCount()) |i| {
        if (input.chain_id_full) |chains| {
            if (std.mem.eql(u8, chains[i], target)) return true;
        } else if (input.chain_id) |chains| {
            if (chains[i].eqlSlice(target)) return true;
        }
    }
    return false;
}

fn bsaInputContainsAllChains(input: AtomInput, chains: []const []const u8) bool {
    for (chains) |chain| {
        if (!bsaInputContainsChain(input, chain)) return false;
    }
    return true;
}

fn writeBsaInterfaceError(
    jsonl_stream: *JsonlStreamWriter,
    allocator: Allocator,
    filename: []const u8,
    id: []const u8,
    name: []const u8,
    error_message: []const u8,
) !BsaInterfaceStats {
    try writeBsaAnalysisErrorJsonl(jsonl_stream, allocator, .{
        .filename = filename,
        .id = id,
        .name = name,
        .error_message = error_message,
    });
    return .{
        .successful = false,
        .failure = .{ .filename = filename, .id = id, .reason = error_message },
    };
}

fn processBsaInterface(
    allocator: Allocator,
    io: std.Io,
    jsonl_stream: *JsonlStreamWriter,
    filename: []const u8,
    name: []const u8,
    level: []const u8,
    atom_output: bool,
    selection: BsaInterfaceSelection,
    source_input: AtomInput,
    config: BatchConfig,
    sasa_threads: usize,
    luts: *const BatchLuts,
) !BsaInterfaceStats {
    if (bsaPartnersOverlap(selection.partner_a, selection.partner_b)) {
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, "partner chain groups overlap");
    }
    if (!bsaInputContainsAllChains(source_input, selection.partner_a)) {
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, "partner A chain not found");
    }
    if (!bsaInputContainsAllChains(source_input, selection.partner_b)) {
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, "partner B chain not found");
    }

    const complex_chains = try appendChainGroups(allocator, selection.partner_a, selection.partner_b);
    const partner_a_input = try copySelectedAtomInput(allocator, source_input, selection.partner_a);
    const partner_b_input = try copySelectedAtomInput(allocator, source_input, selection.partner_b);
    const complex_input = try copySelectedAtomInput(allocator, source_input, complex_chains);

    const partner_a_result = switch (config.precision) {
        .f64 => calculatePreparedInputResult(f64, allocator, io, allocator, partner_a_input, null, filename, config, sasa_threads, luts.f64Ptr(), luts.coarseF64Ptr(), luts.fineF64Ptr()),
        .f32 => calculatePreparedInputResult(f32, allocator, io, allocator, partner_a_input, null, filename, config, sasa_threads, luts.f32Ptr(), luts.coarseF32Ptr(), luts.fineF32Ptr()),
    };
    const partner_b_result = switch (config.precision) {
        .f64 => calculatePreparedInputResult(f64, allocator, io, allocator, partner_b_input, null, filename, config, sasa_threads, luts.f64Ptr(), luts.coarseF64Ptr(), luts.fineF64Ptr()),
        .f32 => calculatePreparedInputResult(f32, allocator, io, allocator, partner_b_input, null, filename, config, sasa_threads, luts.f32Ptr(), luts.coarseF32Ptr(), luts.fineF32Ptr()),
    };
    const complex_result = switch (config.precision) {
        .f64 => calculatePreparedInputResult(f64, allocator, io, allocator, complex_input, null, filename, config, sasa_threads, luts.f64Ptr(), luts.coarseF64Ptr(), luts.fineF64Ptr()),
        .f32 => calculatePreparedInputResult(f32, allocator, io, allocator, complex_input, null, filename, config, sasa_threads, luts.f32Ptr(), luts.coarseF32Ptr(), luts.fineF32Ptr()),
    };

    if (partner_a_result.status != .ok or partner_b_result.status != .ok or complex_result.status != .ok) {
        const detail = partner_a_result.error_msg orelse partner_b_result.error_msg orelse complex_result.error_msg orelse "unknown error";
        const message = try std.fmt.allocPrint(allocator, "SASA calculation failed: {s}", .{detail});
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, message);
    }

    const partner_a_areas = partner_a_result.atom_areas orelse
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, "partner A atom areas unavailable");
    const partner_b_areas = partner_b_result.atom_areas orelse
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, "partner B atom areas unavailable");
    const complex_areas = complex_result.atom_areas orelse
        return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, "complex atom areas unavailable");

    const atom_sasa_isolated = try allocator.alloc(f64, complex_input.atomCount());
    const atom_delta_sasa = try allocator.alloc(f64, complex_input.atomCount());
    var a_index: usize = 0;
    var b_index: usize = 0;
    for (0..complex_input.atomCount()) |i| {
        atom_sasa_isolated[i] = if (atomSelectedByChains(complex_input, i, selection.partner_a)) blk: {
            defer a_index += 1;
            break :blk partner_a_areas[a_index];
        } else blk: {
            defer b_index += 1;
            break :blk partner_b_areas[b_index];
        };
        atom_delta_sasa[i] = atom_sasa_isolated[i] - complex_areas[i];
    }

    var residue_arrays = BsaResidueDeltaArrays{};
    if (std.mem.eql(u8, level, "residue")) {
        residue_arrays = buildBsaResidueDeltaArrays(
            allocator,
            complex_input,
            atom_sasa_isolated,
            complex_areas,
            selection.partner_a,
        ) catch |err| {
            const message = try std.fmt.allocPrint(allocator, "residue detail failed: {s}", .{@errorName(err)});
            return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, message);
        };
    }

    var atom_arrays = BsaAtomArrays{};
    if (atom_output) {
        atom_arrays = buildBsaAtomArrays(allocator, complex_input, selection.partner_a) catch |err| {
            const message = try std.fmt.allocPrint(allocator, "atom detail failed: {s}", .{@errorName(err)});
            return writeBsaInterfaceError(jsonl_stream, allocator, filename, selection.id, name, message);
        };
    }

    const delta_total = partner_a_result.total_sasa + partner_b_result.total_sasa - complex_result.total_sasa;
    try writeBsaAnalysisJsonl(jsonl_stream, allocator, .{
        .filename = filename,
        .id = selection.id,
        .name = name,
        .partner_a = selection.partner_a,
        .partner_b = selection.partner_b,
        .sasa_partner_a = partner_a_result.total_sasa,
        .sasa_partner_b = partner_b_result.total_sasa,
        .sasa_complex = complex_result.total_sasa,
        .delta_sasa_total = delta_total,
        .bsa = delta_total / 2.0,
        .delta_sasa_level = level,
        .atom_output = atom_output,
        .residue_partner = residue_arrays.residue_partner,
        .residue_chain = residue_arrays.residue_chain,
        .residue_name = residue_arrays.residue_name,
        .residue_number = residue_arrays.residue_number,
        .residue_insertion_code = residue_arrays.residue_insertion_code,
        .residue_sasa_isolated = residue_arrays.residue_sasa_isolated,
        .residue_sasa_complex = residue_arrays.residue_sasa_complex,
        .residue_delta_sasa = residue_arrays.residue_delta_sasa,
        .atom_index = atom_arrays.atom_index,
        .atom_partner = atom_arrays.atom_partner,
        .atom_chain = atom_arrays.atom_chain,
        .atom_residue_name = atom_arrays.atom_residue_name,
        .atom_residue_number = atom_arrays.atom_residue_number,
        .atom_insertion_code = atom_arrays.atom_insertion_code,
        .atom_name = atom_arrays.atom_name,
        .atom_element = atom_arrays.atom_element,
        .atom_sasa_isolated = if (atom_output) atom_sasa_isolated else &.{},
        .atom_sasa_complex = if (atom_output) complex_areas else &.{},
        .atom_delta_sasa = if (atom_output) atom_delta_sasa else &.{},
    }, jsonlOptions(config));

    return .{
        .successful = true,
        .sasa_time_ns = partner_a_result.sasa_time_ns + partner_b_result.sasa_time_ns + complex_result.sasa_time_ns,
    };
}

fn effectiveBsaSasaThreads(config: BatchConfig, file_threads: usize) usize {
    return if (file_threads > 1) 1 else effectiveWorkflowSasaThreads(config);
}

const BsaParallelContext = struct {
    files: []const []const u8,
    input_dir: []const u8,
    interface_map: ?*const chain_map.InterfaceMap,
    fixed_partner_a: ?[]const []const u8,
    fixed_partner_b: ?[]const []const u8,
    name: []const u8,
    level: []const u8,
    atom_output: bool,
    config: BatchConfig,
    sasa_threads: usize,
    luts: *const BatchLuts,
    jsonl_stream: *JsonlStreamWriter,
    next_file: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    progress_node: std.Progress.Node,
    counter: BsaWorkflowCounter = .{},
    /// Failed interfaces for the report at the end of the run.
    failures: *FailureLog,
    worker_failed: std.atomic.Value(bool) = std.atomic.Value(bool).init(false),
    io: std.Io,
};

fn markBsaFileCompleted(ctx: *BsaParallelContext) void {
    _ = ctx.processed_count.fetchAdd(1, .release);
    ctx.progress_node.completeOne();
}

fn recordBsaInterfaceStats(ctx: *BsaParallelContext, stats: BsaInterfaceStats) void {
    if (stats.successful) {
        _ = ctx.counter.successful.fetchAdd(1, .monotonic);
        _ = ctx.counter.total_sasa_time_ns.fetchAdd(stats.sasa_time_ns, .monotonic);
    } else {
        _ = ctx.counter.failed.fetchAdd(1, .monotonic);
        if (stats.failure) |failure| ctx.failures.record(ctx.io, failure.filename, failure.id, failure.reason);
    }
}

fn writeBsaSourceErrorRows(
    ctx: *BsaParallelContext,
    allocator: Allocator,
    filename: []const u8,
    map_indices: ?[]const usize,
    message: []const u8,
) !void {
    if (ctx.interface_map) |map| {
        if (map_indices) |indices| {
            for (indices) |entry_index| {
                const entry = map.entries[entry_index];
                const stats = try writeBsaInterfaceError(
                    ctx.jsonl_stream,
                    allocator,
                    filename,
                    entry.id orelse filename,
                    ctx.name,
                    message,
                );
                recordBsaInterfaceStats(ctx, stats);
            }
            return;
        }
    }

    const stats = try writeBsaInterfaceError(
        ctx.jsonl_stream,
        allocator,
        filename,
        filename,
        ctx.name,
        message,
    );
    recordBsaInterfaceStats(ctx, stats);
}

fn processBsaFile(ctx: *BsaParallelContext, allocator: Allocator, filename: []const u8) !void {
    const map_indices: ?[]const usize = if (ctx.interface_map) |map| map.getIndices(filename) else null;
    if (ctx.interface_map != null and map_indices == null) {
        return writeBsaSourceErrorRows(ctx, allocator, filename, null, "interface chain map entry not found");
    }

    var source_config = ctx.config;
    source_config.chain_filter = null;
    if (ctx.interface_map) |map| {
        const indices = map_indices.?;
        const first_type = map.entries[indices[0]].asym_id_type;
        for (indices[1..]) |entry_index| {
            if (map.entries[entry_index].asym_id_type != first_type) {
                return writeBsaSourceErrorRows(
                    ctx,
                    allocator,
                    filename,
                    indices,
                    "interfaces for one file must use the same asym_id_type",
                );
            }
        }
        source_config.use_auth_chain = first_type == .auth;
    }

    const input_path = try std.fs.path.join(allocator, &.{ ctx.input_dir, filename });
    var source_parsed = readInputFile(allocator, ctx.io, input_path, source_config) catch |err| {
        const message = try std.fmt.allocPrint(allocator, "read/parse failed: {s}", .{@errorName(err)});
        return writeBsaSourceErrorRows(ctx, allocator, filename, map_indices, message);
    };
    defer source_parsed.deinit();

    workflowClassifySourceInput(&source_parsed.input, source_parsed.inlineCcdPtr(), input_path, source_config) catch |err| {
        const message = try std.fmt.allocPrint(allocator, "classifier failed: {s}", .{@errorName(err)});
        return writeBsaSourceErrorRows(ctx, allocator, filename, map_indices, message);
    };

    if (ctx.interface_map) |map| {
        for (map_indices.?) |entry_index| {
            const entry = map.entries[entry_index];
            var interface_arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
            defer interface_arena.deinit();
            const stats = try processBsaInterface(
                interface_arena.allocator(),
                ctx.io,
                ctx.jsonl_stream,
                filename,
                ctx.name,
                ctx.level,
                ctx.atom_output,
                .{
                    .id = entry.id orelse filename,
                    .partner_a = entry.partner_a,
                    .partner_b = entry.partner_b,
                },
                source_parsed.input,
                ctx.config,
                ctx.sasa_threads,
                ctx.luts,
            );
            recordBsaInterfaceStats(ctx, stats);
        }
    } else {
        var interface_arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
        defer interface_arena.deinit();
        const stats = try processBsaInterface(
            interface_arena.allocator(),
            ctx.io,
            ctx.jsonl_stream,
            filename,
            ctx.name,
            ctx.level,
            ctx.atom_output,
            .{
                .id = filename,
                .partner_a = ctx.fixed_partner_a.?,
                .partner_b = ctx.fixed_partner_b.?,
            },
            source_parsed.input,
            ctx.config,
            ctx.sasa_threads,
            ctx.luts,
        );
        recordBsaInterfaceStats(ctx, stats);
    }
}

fn bsaParallelWorker(ctx: *BsaParallelContext) void {
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();

    while (true) {
        const file_index = ctx.next_file.fetchAdd(1, .monotonic);
        if (file_index >= ctx.files.len) break;
        const filename = ctx.files[file_index];
        defer markBsaFileCompleted(ctx);

        processBsaFile(ctx, arena.allocator(), filename) catch |err| {
            ctx.worker_failed.store(true, .release);
            std.debug.print("Error running BSA analysis on '{s}': {s}\n", .{ filename, @errorName(err) });
        };
        _ = arena.reset(.retain_capacity);
    }
}

fn runWorkflowBsaAnalysis(
    allocator: Allocator,
    io: std.Io,
    args: BatchArgs,
    workflow: workflow_manifest.Workflow,
) !void {
    const analysis = workflow.analysis orelse return error.InvalidArgument;
    try batchWorkflowFindings(workflow, args).report();
    const fixed_partner_a = analysis.partner_a;
    const fixed_partner_b = analysis.partner_b;
    const level = analysisLevel(analysis);
    const atom_output = analysis.atom_output orelse false;
    const name = analysisName(analysis);

    if (workflow.output.format) |format| {
        if (!std.mem.eql(u8, format, "jsonl")) {
            std.debug.print("Error: BSA analysis workflow output format must be \"jsonl\"\n", .{});
            return error.InvalidArgument;
        }
    }

    const input_dir = args.input_path orelse workflow.input.dir orelse {
        std.debug.print("Error: Missing input directory (provide positional input_dir or [input].dir in workflow)\n", .{});
        return error.MissingArgument;
    };
    const output_dir = args.output_path orelse workflow.output.dir orelse {
        std.debug.print("Error: BSA analysis workflow requires an output directory\n", .{});
        return error.InvalidArgument;
    };

    try std.Io.Dir.cwd().createDirPath(io, output_dir);
    const jsonl_output_path = try workflowAnalysisJsonlOutputPath(allocator, output_dir, name);
    defer allocator.free(jsonl_output_path);
    try truncateJsonlOutput(io, jsonl_output_path);
    const jsonl_file = try std.Io.Dir.cwd().openFile(io, jsonl_output_path, .{ .mode = .write_only });
    defer jsonl_file.close(io);
    var jsonl_buffer: [64 * 1024]u8 = undefined;
    var jsonl_stream = JsonlStreamWriter.init(jsonl_file, io, .{}, &jsonl_buffer);
    errdefer jsonl_stream.flush() catch {};

    const load_quiet = if (args.quiet_explicit) args.quiet else (workflow.calculation.quiet orelse args.quiet);

    var config = BatchConfig{};
    try applyWorkflowToBatchConfig(&config, args, workflow.calculation, workflow.output, workflow.classifier);
    applyCliOverrides(&config, args);
    try validateBitmaskCorrectionConfig(config);
    config.store_atom_areas = true;
    config.residue_map = false;
    config.output_format = .jsonl;
    config.chain_filter = null;
    const effective_classifier_type = config.classifier_type;

    if (config.jsonl_metadata == .sidecar) {
        const metadata_path = try workflowJsonlMetadataPath(allocator, output_dir, name);
        defer allocator.free(metadata_path);
        try writeWorkflowJsonlMetadata(allocator, io, metadata_path, name, config);
    }

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

    config.external_ccd = if (ext_ccd != null) &ext_ccd.? else null;
    config.sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null;
    if (custom_classifier) |*c| config.custom_classifier = c;

    var interface_map: ?chain_map.InterfaceMap = null;
    if (analysis.chain_map) |path| {
        interface_map = chain_map.loadInterfaceFile(allocator, io, path) catch |err| {
            std.debug.print("Error loading interface chain map '{s}' for BSA analysis '{s}': {s}\n", .{ path, name, @errorName(err) });
            return err;
        };
    }
    defer if (interface_map) |*map| map.deinit();

    if ((config.classifier_type == .ccd or config.classifier_type == .protor) and config.include_hydrogens and !config.quiet) {
        std.debug.print("Warning: --include-hydrogens with CCD classifier may give inaccurate results\n", .{});
        std.debug.print("         CCD uses united-atom radii that already account for implicit hydrogens\n", .{});
    }

    const files = try scanDirectory(allocator, io, input_dir);
    defer freeScannedFiles(allocator, files);
    try validateMappedChainInputFormats(files, interface_map != null);

    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    var discovered_files = std.StringHashMapUnmanaged(void){};
    defer discovered_files.deinit(allocator);
    for (files) |filename| try discovered_files.put(allocator, filename, {});

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing files", .estimated_total_items = files.len })
    else
        .none;
    defer progress_root.end();

    var failures = FailureLog{ .allocator = allocator };
    defer failures.deinit();

    const file_threads = effectiveWorkflowFileThreads(config, files.len);
    var ctx = BsaParallelContext{
        .files = files,
        .input_dir = input_dir,
        .interface_map = if (interface_map) |*map| map else null,
        .fixed_partner_a = fixed_partner_a,
        .fixed_partner_b = fixed_partner_b,
        .name = name,
        .level = level,
        .atom_output = atom_output,
        .config = config,
        .sasa_threads = effectiveBsaSasaThreads(config, file_threads),
        .luts = &luts,
        .jsonl_stream = &jsonl_stream,
        .next_file = std.atomic.Value(usize).init(0),
        .processed_count = std.atomic.Value(usize).init(0),
        .progress_node = progress_root,
        .failures = &failures,
        .io = io,
    };

    if (file_threads > 1 and files.len > 1) {
        const threads = try allocator.alloc(std.Thread, file_threads);
        defer allocator.free(threads);
        var spawned_count: usize = 0;
        errdefer joinSpawnedThreads(threads, spawned_count);
        for (threads) |*thread| {
            thread.* = try std.Thread.spawn(.{}, bsaParallelWorker, .{&ctx});
            spawned_count += 1;
        }
        joinSpawnedThreads(threads, spawned_count);
    } else {
        bsaParallelWorker(&ctx);
    }

    if (ctx.processed_count.load(.acquire) != files.len) return error.BsaWorkerFailed;

    if (interface_map) |*map| {
        for (map.entries) |entry| {
            if (discovered_files.contains(entry.filename)) continue;
            var error_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
            defer error_arena.deinit();
            const stats = try writeBsaInterfaceError(
                &jsonl_stream,
                error_arena.allocator(),
                entry.filename,
                entry.id orelse entry.filename,
                name,
                "input structure not found",
            );
            recordBsaInterfaceStats(&ctx, stats);
        }
    }

    try jsonl_stream.flush();
    if (jsonl_stream.hasError()) return error.JsonlWriteFailed;
    if (ctx.worker_failed.load(.acquire)) return error.BsaWorkerFailed;

    const successful = ctx.counter.successful.load(.monotonic);
    const failed = ctx.counter.failed.load(.monotonic);
    const total_sasa_time_ns = ctx.counter.total_sasa_time_ns.load(.monotonic);
    if (config.show_timing) {
        const total_sasa_ms = @as(f64, @floatFromInt(total_sasa_time_ns)) / 1_000_000.0;
        std.debug.print("BSA analysis SASA time: {d:.2} ms\n", .{total_sasa_ms});
    }
    std.debug.print("Workflow complete: {d} successful, {d} failed\n", .{ successful, failed });
    (FailureReport{
        .unit = .interface,
        .total = successful + failed,
        .failed = failed,
        .failures = failures.sorted(),
        .jsonl = .{ .file = jsonl_output_path },
    }).print(allocator);
}

const SelectionMapBatchStats = struct {
    successful: usize,
    failed: usize,
    total_sasa_time_ns: u64,
    read_parse_count: usize,
    classifier_count: usize,
    calculation_count: usize,
    file_threads: usize,
    sasa_threads: usize,
};

const SelectionMapCounter = struct {
    successful: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    failed: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    total_sasa_time_ns: std.atomic.Value(u64) = std.atomic.Value(u64).init(0),
    read_parse_count: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    classifier_count: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    calculation_count: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
};

const SelectionGroup = struct {
    representative_entry_index: usize,
    result: FileResult,
    atom_identity: ?json_writer.SelectionAtomIdentity = null,
};

fn selectionMapEstimatedFileCost(map: *const chain_map.ChainMap, filename: []const u8) usize {
    return if (map.getIndices(filename)) |indices| indices.len else 0;
}

fn sortSelectionMapFilesByEstimatedCost(files: [][]const u8, map: *const chain_map.ChainMap) void {
    std.mem.sort([]const u8, files, map, struct {
        fn lessThan(chain_map_value: *const chain_map.ChainMap, a: []const u8, b: []const u8) bool {
            const a_cost = selectionMapEstimatedFileCost(chain_map_value, a);
            const b_cost = selectionMapEstimatedFileCost(chain_map_value, b);
            if (a_cost != b_cost) return a_cost > b_cost;
            return std.mem.lessThan(u8, a, b);
        }
    }.lessThan);
}

fn scheduleSelectionMapFiles(files: [][]const u8, map: *const chain_map.ChainMap) void {
    if (map.hasMultipleSelections()) sortSelectionMapFilesByEstimatedCost(files, map);
}

fn claimSelectionMapFile(files: []const []const u8, next_file: *std.atomic.Value(usize)) ?[]const u8 {
    const file_index = next_file.fetchAdd(1, .monotonic);
    if (file_index >= files.len) return null;
    return files[file_index];
}

const SelectionMapContext = struct {
    files: []const []const u8,
    input_dir: []const u8,
    map: *const chain_map.ChainMap,
    config: BatchConfig,
    sasa_threads: usize,
    luts: *const BatchLuts,
    jsonl_stream: *JsonlStreamWriter,
    /// Failed selections for the report at the end of the run, if wanted.
    failures: ?*FailureLog = null,
    next_file: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    progress_node: std.Progress.Node,
    counter: SelectionMapCounter = .{},
    worker_failed: std.atomic.Value(bool) = std.atomic.Value(bool).init(false),
    io: std.Io,
};

fn chainsContain(chains: []const []const u8, target: []const u8) bool {
    for (chains) |chain| {
        if (std.mem.eql(u8, chain, target)) return true;
    }
    return false;
}

fn uniqueChainCount(chains: []const []const u8) usize {
    var count: usize = 0;
    for (chains, 0..) |chain, i| {
        var seen = false;
        for (chains[0..i]) |previous| {
            if (std.mem.eql(u8, chain, previous)) {
                seen = true;
                break;
            }
        }
        if (!seen) count += 1;
    }
    return count;
}

fn sameCanonicalChainSet(a: []const []const u8, b: []const []const u8) bool {
    if (uniqueChainCount(a) != uniqueChainCount(b)) return false;
    for (a) |chain| {
        if (!chainsContain(b, chain)) return false;
    }
    return true;
}

fn missingSelectedChain(input: AtomInput, chains: []const []const u8) ?[]const u8 {
    for (chains) |chain| {
        if (!bsaInputContainsChain(input, chain)) return chain;
    }
    return null;
}

fn selectedSourceAtomIndices(allocator: Allocator, source_input: AtomInput, chains: []const []const u8) ![]const usize {
    const count = countSelectedAtoms(source_input, chains);
    const indices = try allocator.alloc(usize, count);
    var out_index: usize = 0;
    for (0..source_input.atomCount()) |source_index| {
        if (!atomSelectedByChains(source_input, source_index, chains)) continue;
        indices[out_index] = source_index;
        out_index += 1;
    }
    return indices;
}

fn buildSelectionAtomIdentity(
    allocator: Allocator,
    selected_input: AtomInput,
    source_atom_index: []const usize,
) !json_writer.SelectionAtomIdentity {
    if (!selected_input.hasResidueInfo() or selected_input.atom_name == null or selected_input.element == null) {
        return error.MissingAtomMetadata;
    }
    if (source_atom_index.len != selected_input.atomCount()) return error.LengthMismatch;

    const atom_count = selected_input.atomCount();
    const atom_chain = try allocator.alloc([]const u8, atom_count);
    const atom_residue_name = try allocator.alloc([]const u8, atom_count);
    const atom_residue_number = try allocator.dupe(i32, selected_input.residue_num.?);
    const atom_insertion_code = try allocator.alloc([]const u8, atom_count);
    const atom_name = try allocator.alloc([]const u8, atom_count);
    const atom_element = try allocator.alloc([]const u8, atom_count);
    for (0..atom_count) |i| {
        atom_chain[i] = if (selected_input.chain_id_full) |chains|
            chains[i]
        else
            selected_input.chain_id.?[i].slice();
        atom_residue_name[i] = selected_input.residue.?[i].slice();
        atom_insertion_code[i] = selected_input.insertion_code.?[i].slice();
        atom_name[i] = selected_input.atom_name.?[i].slice();
        atom_element[i] = element_module.fromAtomicNumber(selected_input.element.?[i]).symbol();
    }
    return .{
        .source_atom_index = source_atom_index,
        .atom_chain = atom_chain,
        .atom_residue_name = atom_residue_name,
        .atom_residue_number = atom_residue_number,
        .atom_insertion_code = atom_insertion_code,
        .atom_name = atom_name,
        .atom_element = atom_element,
    };
}

fn selectionErrorResult(allocator: Allocator, filename: []const u8, comptime format: []const u8, args: anytype) FileResult {
    return .{
        .filename = filename,
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .err,
        .error_msg = std.fmt.allocPrint(allocator, format, args) catch null,
    };
}

fn writeSelectionError(
    stream: *JsonlStreamWriter,
    allocator: Allocator,
    filename: []const u8,
    id: []const u8,
    chains: []const []const u8,
    message: []const u8,
) !void {
    const line = try json_writer.selectionErrorToJsonlLine(allocator, .{
        .filename = filename,
        .id = id,
        .chains = chains,
        .error_message = message,
    });
    defer allocator.free(line);
    try stream.writeLine(filename, line);
}

fn writeSelectionSuccess(
    stream: *JsonlStreamWriter,
    allocator: Allocator,
    entry: chain_map.Entry,
    result: FileResult,
    identity: ?json_writer.SelectionAtomIdentity,
) !void {
    const line = try json_writer.selectionResultToJsonlLineOptions(allocator, .{
        .filename = entry.filename,
        .id = entry.id orelse entry.filename,
        .chains = entry.chains,
        .total_area = result.total_sasa,
        .atom_areas = result.atom_areas orelse &.{},
        .residue_map = result.residue_map,
        .atom_identity = identity,
    }, stream.options);
    defer allocator.free(line);
    try stream.writeLine(entry.filename, line);
}

/// Write the error row of one failed selection, count it and keep it for the
/// failure report.
fn failSelection(
    ctx: *SelectionMapContext,
    allocator: Allocator,
    filename: []const u8,
    id: []const u8,
    chains: []const []const u8,
    message: []const u8,
) !void {
    try writeSelectionError(ctx.jsonl_stream, allocator, filename, id, chains, message);
    _ = ctx.counter.failed.fetchAdd(1, .monotonic);
    if (ctx.failures) |log| log.record(ctx.io, filename, id, message);
}

fn writeSelectionSourceErrorRows(
    ctx: *SelectionMapContext,
    allocator: Allocator,
    filename: []const u8,
    map_indices: ?[]const usize,
    message: []const u8,
) !void {
    if (map_indices) |indices| {
        for (indices) |entry_index| {
            const entry = ctx.map.entries[entry_index];
            try failSelection(ctx, allocator, filename, entry.id orelse filename, entry.chains, message);
        }
        return;
    }
    try failSelection(ctx, allocator, filename, filename, &.{}, message);
}

fn calculateSelectionGroup(
    ctx: *SelectionMapContext,
    allocator: Allocator,
    filename: []const u8,
    source_input: AtomInput,
    representative_entry_index: usize,
) SelectionGroup {
    const entry = ctx.map.entries[representative_entry_index];
    if (missingSelectedChain(source_input, entry.chains)) |missing_chain| {
        return .{
            .representative_entry_index = representative_entry_index,
            .result = selectionErrorResult(allocator, filename, "selected chain not found: {s}", .{missing_chain}),
        };
    }

    const source_indices: []const usize = if (ctx.config.jsonl_include_atom_identity)
        selectedSourceAtomIndices(allocator, source_input, entry.chains) catch |err| {
            return .{
                .representative_entry_index = representative_entry_index,
                .result = selectionErrorResult(allocator, filename, "selection source-index mapping failed: {s}", .{@errorName(err)}),
            };
        }
    else
        &.{};
    const selected_input = copySelectedAtomInput(allocator, source_input, entry.chains) catch |err| {
        return .{
            .representative_entry_index = representative_entry_index,
            .result = selectionErrorResult(allocator, filename, "selection failed: {s}", .{@errorName(err)}),
        };
    };

    _ = ctx.counter.calculation_count.fetchAdd(1, .monotonic);
    const result = switch (ctx.config.precision) {
        .f64 => calculatePreparedInputResult(
            f64,
            allocator,
            ctx.io,
            allocator,
            selected_input,
            null,
            filename,
            ctx.config,
            ctx.sasa_threads,
            ctx.luts.f64Ptr(),
            ctx.luts.coarseF64Ptr(),
            ctx.luts.fineF64Ptr(),
        ),
        .f32 => calculatePreparedInputResult(
            f32,
            allocator,
            ctx.io,
            allocator,
            selected_input,
            null,
            filename,
            ctx.config,
            ctx.sasa_threads,
            ctx.luts.f32Ptr(),
            ctx.luts.coarseF32Ptr(),
            ctx.luts.fineF32Ptr(),
        ),
    };
    if (result.status != .ok) {
        return .{
            .representative_entry_index = representative_entry_index,
            .result = result,
        };
    }

    const identity: ?json_writer.SelectionAtomIdentity = if (ctx.config.jsonl_include_atom_identity)
        buildSelectionAtomIdentity(allocator, selected_input, source_indices) catch |err| {
            return .{
                .representative_entry_index = representative_entry_index,
                .result = selectionErrorResult(allocator, filename, "atom identity failed: {s}", .{@errorName(err)}),
            };
        }
    else
        null;

    return .{
        .representative_entry_index = representative_entry_index,
        .result = result,
        .atom_identity = identity,
    };
}

fn processSelectionMapFile(ctx: *SelectionMapContext, allocator: Allocator, filename: []const u8) !void {
    const map_indices = ctx.map.getIndices(filename) orelse
        return writeSelectionSourceErrorRows(ctx, allocator, filename, null, "selection chain map entry not found");

    const first_type = ctx.map.entries[map_indices[0]].asym_id_type;
    for (map_indices[1..]) |entry_index| {
        if (ctx.map.entries[entry_index].asym_id_type != first_type) {
            return writeSelectionSourceErrorRows(
                ctx,
                allocator,
                filename,
                map_indices,
                "selections for one file must use the same asym_id_type",
            );
        }
    }

    var source_config = ctx.config;
    source_config.chain_filter = null;
    source_config.use_auth_chain = first_type == .auth;
    const input_path = try std.fs.path.join(allocator, &.{ ctx.input_dir, filename });

    _ = ctx.counter.read_parse_count.fetchAdd(1, .monotonic);
    var source_parsed = readInputFile(allocator, ctx.io, input_path, source_config) catch |err| {
        const message = try std.fmt.allocPrint(allocator, "read/parse failed: {s}", .{@errorName(err)});
        return writeSelectionSourceErrorRows(ctx, allocator, filename, map_indices, message);
    };
    defer source_parsed.deinit();

    _ = ctx.counter.classifier_count.fetchAdd(1, .monotonic);
    workflowClassifySourceInput(&source_parsed.input, source_parsed.inlineCcdPtr(), input_path, source_config) catch |err| {
        const message = try std.fmt.allocPrint(allocator, "classifier failed: {s}", .{@errorName(err)});
        return writeSelectionSourceErrorRows(ctx, allocator, filename, map_indices, message);
    };

    var groups = std.ArrayListUnmanaged(SelectionGroup).empty;
    for (map_indices) |entry_index| {
        const entry = ctx.map.entries[entry_index];
        var group_index: ?usize = null;
        for (groups.items, 0..) |group, i| {
            const representative = ctx.map.entries[group.representative_entry_index];
            if (sameCanonicalChainSet(entry.chains, representative.chains)) {
                group_index = i;
                break;
            }
        }
        if (group_index == null) {
            try groups.append(allocator, calculateSelectionGroup(ctx, allocator, filename, source_parsed.input, entry_index));
        }
    }

    for (map_indices) |entry_index| {
        const entry = ctx.map.entries[entry_index];
        const group = for (groups.items) |candidate| {
            const representative = ctx.map.entries[candidate.representative_entry_index];
            if (sameCanonicalChainSet(entry.chains, representative.chains)) break candidate;
        } else unreachable;

        if (group.result.status == .ok) {
            try writeSelectionSuccess(ctx.jsonl_stream, allocator, entry, group.result, group.atom_identity);
            _ = ctx.counter.successful.fetchAdd(1, .monotonic);
        } else {
            try failSelection(
                ctx,
                allocator,
                filename,
                entry.id orelse filename,
                entry.chains,
                group.result.error_msg orelse "unknown error",
            );
        }
    }

    for (groups.items) |group| {
        if (group.result.status == .ok) {
            _ = ctx.counter.total_sasa_time_ns.fetchAdd(group.result.sasa_time_ns, .monotonic);
        }
    }
}

fn selectionMapWorker(ctx: *SelectionMapContext) void {
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();

    while (true) {
        const filename = claimSelectionMapFile(ctx.files, &ctx.next_file) orelse break;
        defer {
            _ = ctx.processed_count.fetchAdd(1, .release);
            ctx.progress_node.completeOne();
        }

        processSelectionMapFile(ctx, arena.allocator(), filename) catch |err| {
            ctx.worker_failed.store(true, .release);
            std.debug.print("Error running selection map on '{s}': {s}\n", .{ filename, @errorName(err) });
        };
        _ = arena.reset(.retain_capacity);
    }
}

fn runSelectionMapBatch(
    allocator: Allocator,
    io: std.Io,
    input_dir: []const u8,
    config: BatchConfig,
    jsonl_output_path: ?[]const u8,
    map: *const chain_map.ChainMap,
    failures: ?*FailureLog,
) !SelectionMapBatchStats {
    if (config.output_format != .jsonl) return error.MultiSelectionMapRequiresJsonl;
    if (!config.jsonl_include_total_area) return error.SelectionMapRequiresTotalArea;
    if (config.jsonl_include_atom_identity and !config.jsonl_include_atom_areas) {
        return error.AtomIdentityRequiresAtomAreas;
    }

    const files = try scanDirectory(allocator, io, input_dir);
    defer freeScannedFiles(allocator, files);
    try validateMappedChainInputFormats(files, true);
    scheduleSelectionMapFiles(files, map);

    var luts = try BatchLuts.init(allocator, config);
    defer luts.deinit();

    const jsonl_file = if (jsonl_output_path) |path|
        try std.Io.Dir.cwd().createFile(io, path, .{})
    else
        std.Io.File.stdout();
    defer if (jsonl_output_path != null) jsonl_file.close(io);
    var jsonl_buffer: [64 * 1024]u8 = undefined;
    var jsonl_stream = JsonlStreamWriter.init(jsonl_file, io, jsonlOptions(config), &jsonl_buffer);
    errdefer jsonl_stream.flush() catch {};

    var discovered_files = std.StringHashMapUnmanaged(void){};
    defer discovered_files.deinit(allocator);
    for (files) |filename| try discovered_files.put(allocator, filename, {});

    var progress_root: std.Progress.Node = if (shouldShowProgress(config))
        std.Progress.start(io, .{ .root_name = "Processing files", .estimated_total_items = files.len })
    else
        .none;
    defer progress_root.end();

    const file_threads = effectiveWorkflowFileThreads(config, files.len);
    const sasa_threads = if (file_threads > 1) 1 else effectiveWorkflowSasaThreads(config);
    var ctx = SelectionMapContext{
        .files = files,
        .input_dir = input_dir,
        .map = map,
        .config = config,
        .sasa_threads = sasa_threads,
        .luts = &luts,
        .jsonl_stream = &jsonl_stream,
        .failures = failures,
        .next_file = std.atomic.Value(usize).init(0),
        .processed_count = std.atomic.Value(usize).init(0),
        .progress_node = progress_root,
        .io = io,
    };

    if (file_threads > 1 and files.len > 1) {
        const threads = try allocator.alloc(std.Thread, file_threads);
        defer allocator.free(threads);
        var spawned_count: usize = 0;
        errdefer joinSpawnedThreads(threads, spawned_count);
        for (threads) |*thread| {
            thread.* = try std.Thread.spawn(.{}, selectionMapWorker, .{&ctx});
            spawned_count += 1;
        }
        joinSpawnedThreads(threads, spawned_count);
    } else {
        selectionMapWorker(&ctx);
    }

    if (ctx.processed_count.load(.acquire) != files.len) return error.SelectionMapWorkerFailed;

    for (map.entries) |entry| {
        if (discovered_files.contains(entry.filename)) continue;
        var error_arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
        defer error_arena.deinit();
        try failSelection(
            &ctx,
            error_arena.allocator(),
            entry.filename,
            entry.id orelse entry.filename,
            entry.chains,
            "input structure not found",
        );
    }

    try jsonl_stream.flush();
    if (jsonl_stream.hasError()) return error.JsonlWriteFailed;
    if (ctx.worker_failed.load(.acquire)) return error.SelectionMapWorkerFailed;

    return .{
        .successful = ctx.counter.successful.load(.monotonic),
        .failed = ctx.counter.failed.load(.monotonic),
        .total_sasa_time_ns = ctx.counter.total_sasa_time_ns.load(.monotonic),
        .read_parse_count = ctx.counter.read_parse_count.load(.monotonic),
        .classifier_count = ctx.counter.classifier_count.load(.monotonic),
        .calculation_count = ctx.counter.calculation_count.load(.monotonic),
        .file_threads = file_threads,
        .sasa_threads = sasa_threads,
    };
}

fn effectiveWorkflowSasaThreads(config: BatchConfig) usize {
    return if (config.n_threads == 0)
        std.Thread.getCpuCount() catch 1
    else
        config.n_threads;
}

fn effectiveWorkflowFileThreads(config: BatchConfig, file_count: usize) usize {
    const cpu_count = std.Thread.getCpuCount() catch 1;
    return @min(resolveBatchThreadCount(config.n_threads, cpu_count), file_count);
}

const WorkflowParallelContext = struct {
    files: []const []const u8,
    input_dir: []const u8,
    jobs: []const workflow_manifest.Job,
    runtimes: []WorkflowJobRuntime,
    resource_config: BatchConfig,
    result_allocator: Allocator,
    next_file: std.atomic.Value(usize),
    processed_count: std.atomic.Value(usize),
    lut_f64: ?*const bitmask_lut.BitmaskLut,
    lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    coarse_lut_f64: ?*const bitmask_lut.BitmaskLut,
    fine_lut_f64: ?*const bitmask_lut.BitmaskLut,
    coarse_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    fine_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32),
    io: std.Io,
};

/// A result that only reports `message` for `filename`; the message is not owned.
fn failedFileResult(filename: []const u8, message: []const u8) FileResult {
    return .{
        .filename = filename,
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .err,
        .error_msg = message,
    };
}

/// Count `filename` as failed for one job of the parallel file-first runner
/// and write its JSONL error row, as the other runners do for a failed input.
fn workflowFailJob(io: std.Io, runtime: *WorkflowJobRuntime, arena: Allocator, filename: []const u8, message: []const u8) void {
    _ = runtime.counter.failed.fetchAdd(1, .monotonic);
    runtime.state.failures.record(io, filename, null, message);
    if (runtime.jsonl_stream) |*stream| {
        var result = failedFileResult(filename, message);
        stream.writeResult(arena, &result);
    }
}

/// `workflowFailJob` for every job: the input could not be read, parsed or
/// classified, so no job has a result for it.
fn workflowFailAllJobs(
    ctx: *WorkflowParallelContext,
    arena: Allocator,
    filename: []const u8,
    comptime stage: []const u8,
    err: anyerror,
) void {
    const message = std.fmt.allocPrint(arena, stage ++ " failed: {s}", .{@errorName(err)}) catch stage ++ " failed";
    for (ctx.runtimes) |*runtime| workflowFailJob(ctx.io, runtime, arena, filename, message);
}

fn workflowClassifySourceInput(
    input: *AtomInput,
    inline_ccd: ?*const ccd_parser.ComponentDict,
    input_path: []const u8,
    config: BatchConfig,
) !void {
    if (config.custom_classifier) |custom| {
        if (input.hasClassificationInfo()) {
            try applyCustomClassifier(input, custom, config.quiet);
        }
    } else if (config.classifier_type) |ct| {
        const format = format_detect.detectInputFormat(input_path);
        if (format != .json and input.hasClassificationInfo()) {
            try applyBuiltinClassifier(input, ct, config.sdf_ccd, inline_ccd, config.external_ccd);
        }
    }
}

fn workflowParallelWorker(ctx: *WorkflowParallelContext) void {
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();

    while (true) {
        const file_index = ctx.next_file.fetchAdd(1, .monotonic);
        if (file_index >= ctx.files.len) break;

        const filename = ctx.files[file_index];
        const input_path = std.fs.path.join(arena.allocator(), &.{ ctx.input_dir, filename }) catch |err| {
            workflowFailAllJobs(ctx, arena.allocator(), filename, "path join", err);
            _ = ctx.processed_count.fetchAdd(1, .release);
            _ = arena.reset(.retain_capacity);
            continue;
        };

        var source_config = ctx.resource_config;
        source_config.chain_filter = null;
        var source_parsed = readInputFile(arena.allocator(), ctx.io, input_path, source_config) catch |err| {
            workflowFailAllJobs(ctx, arena.allocator(), filename, "read/parse", err);
            _ = ctx.processed_count.fetchAdd(1, .release);
            _ = arena.reset(.retain_capacity);
            continue;
        };

        workflowClassifySourceInput(&source_parsed.input, source_parsed.inlineCcdPtr(), input_path, source_config) catch |err| {
            workflowFailAllJobs(ctx, arena.allocator(), filename, "classifier", err);
            source_parsed.deinit();
            _ = ctx.processed_count.fetchAdd(1, .release);
            _ = arena.reset(.retain_capacity);
            continue;
        };

        const format = format_detect.detectInputFormat(input_path);
        for (ctx.jobs, 0..) |job, job_index| {
            var runtime = &ctx.runtimes[job_index];
            const selected_chains: ?[]const []const u8 = if (format == .json) null else job.chains;
            var selected_input = copySelectedAtomInput(arena.allocator(), source_parsed.input, selected_chains) catch |err| {
                const message = std.fmt.allocPrint(arena.allocator(), "selection failed: {s}", .{@errorName(err)}) catch "selection failed";
                workflowFailJob(ctx.io, runtime, arena.allocator(), filename, message);
                continue;
            };
            defer selected_input.deinit();

            const sasa_threads: usize = 1;
            var result = switch (runtime.state.config.precision) {
                .f64 => calculatePreparedInputResult(f64, arena.allocator(), ctx.io, ctx.result_allocator, selected_input, runtime.state.output_dir, filename, runtime.state.config, sasa_threads, ctx.lut_f64, ctx.coarse_lut_f64, ctx.fine_lut_f64),
                .f32 => calculatePreparedInputResult(f32, arena.allocator(), ctx.io, ctx.result_allocator, selected_input, runtime.state.output_dir, filename, runtime.state.config, sasa_threads, ctx.lut_f32, ctx.coarse_lut_f32, ctx.fine_lut_f32),
            };
            defer if (result.error_msg) |msg| ctx.result_allocator.free(msg);

            if (result.status == .ok) {
                _ = runtime.counter.successful.fetchAdd(1, .monotonic);
                _ = runtime.counter.total_sasa_time_ns.fetchAdd(result.sasa_time_ns, .monotonic);
            } else {
                _ = runtime.counter.failed.fetchAdd(1, .monotonic);
                runtime.state.failures.record(ctx.io, filename, null, result.error_msg orelse "unknown error");
            }

            if (runtime.jsonl_stream) |*stream| {
                stream.writeResult(arena.allocator(), &result);
            }

            result.atom_areas = null;
            result.residue_map = null;
        }

        source_parsed.deinit();
        _ = ctx.processed_count.fetchAdd(1, .release);
        _ = arena.reset(.retain_capacity);
    }
}

/// Write the JSONL row of one result of a job in the sequential file-first
/// runner; nothing for a job with per-file output.
fn workflowWriteSequentialJsonl(io: std.Io, arena: Allocator, state: *const WorkflowJobState, result: *FileResult) !void {
    if (!batchWritesJsonl(state.config)) return;
    if (state.jsonl_output_path) |path| {
        const file = try openJsonlForAppend(io, path);
        defer file.close(io);
        try appendJsonlResultToFile(io, file, arena, result, jsonlOptions(state.config));
    } else {
        var stdout_write_buf: [64 * 1024]u8 = undefined;
        var stdout_writer = std.Io.File.Writer.initStreaming(std.Io.File.stdout(), io, &stdout_write_buf);
        try writeJsonlResult(&stdout_writer, arena, result, jsonlOptions(state.config));
    }
}

/// `workflowFailJob` for the sequential file-first runner.
fn workflowFailJobSequential(io: std.Io, arena: Allocator, state: *WorkflowJobState, filename: []const u8, message: []const u8) !void {
    state.failed += 1;
    state.failures.record(io, filename, null, message);
    var result = failedFileResult(filename, message);
    try workflowWriteSequentialJsonl(io, arena, state, &result);
}

/// `workflowFailAllJobs` for the sequential file-first runner.
fn workflowFailAllJobsSequential(
    io: std.Io,
    arena: Allocator,
    states: []WorkflowJobState,
    filename: []const u8,
    comptime stage: []const u8,
    err: anyerror,
) !void {
    const message = std.fmt.allocPrint(arena, stage ++ " failed: {s}", .{@errorName(err)}) catch stage ++ " failed";
    for (states) |*state| try workflowFailJobSequential(io, arena, state, filename, message);
}

/// Report a workflow job that failed as a whole, naming the job and the cause.
fn printWorkflowJobError(job_name: []const u8, comptime cause_fmt: []const u8, cause_args: anytype) void {
    std.debug.print("Error running workflow job '{s}': " ++ cause_fmt ++ "\n", .{job_name} ++ cause_args);
}

/// The line that closes a workflow in which whole jobs failed: how many of
/// the jobs, and which. The cause of each was reported when it failed.
/// Caller frees the result.
fn formatFailedWorkflowJobs(allocator: Allocator, failed_jobs: []const []const u8, n_jobs: usize) ![]u8 {
    var aw = std.Io.Writer.Allocating.init(allocator);
    defer aw.deinit();
    const w = &aw.writer;

    try w.print("{d} of {d} job{s} failed:", .{ failed_jobs.len, n_jobs, if (n_jobs == 1) "" else "s" });
    for (failed_jobs, 0..) |name, i| {
        try w.print("{s} {s}", .{ if (i == 0) "" else ",", name });
    }
    try w.writeByte('\n');
    return aw.toOwnedSlice();
}

fn printFailedWorkflowJobs(allocator: Allocator, failed_jobs: []const []const u8, n_jobs: usize) void {
    const text = formatFailedWorkflowJobs(allocator, failed_jobs, n_jobs) catch {
        std.debug.print("{d} of {d} jobs failed\n", .{ failed_jobs.len, n_jobs });
        return;
    };
    defer allocator.free(text);
    std.debug.print("{s}", .{text});
}

fn runWorkflowJobFirst(allocator: Allocator, io: std.Io, args: BatchArgs) !void {
    if (args.chain_filter != null) {
        std.debug.print("Error: --workflow cannot be combined with --chain; --manifest is a compatibility alias for --workflow; use [[jobs]].chains in the workflow\n", .{});
        return error.InvalidArgument;
    }

    const workflow_path = args.workflow_path.?;
    var workflow = parseWorkflowFile(allocator, io, workflow_path) catch |err| {
        printWorkflowReadError(workflow_path, err);
        return err;
    };
    defer workflow.deinit();

    try batchWorkflowFindings(workflow, args).report();
    try workflow.applyInputChainToJobs();

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
    try validateBitmaskCorrectionConfig(resource_config);
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
    // Jobs that failed as a whole. They are not inputs, so they are kept out
    // of the counters above and make the run an error once every job has had
    // its turn.
    var failed_jobs = std.ArrayListUnmanaged([]const u8).empty;
    defer failed_jobs.deinit(allocator);
    // The failed inputs of each job, reported together when the workflow ends.
    var report_arena = std.heap.ArenaAllocator.init(allocator);
    defer report_arena.deinit();
    var reports = std.ArrayListUnmanaged(FailureReport).empty;

    for (workflow.jobs) |job| {
        var loaded_chain_map: ?chain_map.ChainMap = null;
        defer if (loaded_chain_map) |*map| map.deinit();

        var config = BatchConfig{};
        try applyWorkflowToBatchConfig(&config, args, workflow.calculation, workflow.output, workflow.classifier);
        applyCliOverrides(&config, args);

        config.store_atom_areas = batchShouldStoreAtomAreas(config);
        config.external_ccd = if (ext_ccd != null) &ext_ccd.? else null;
        config.sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null;
        if (custom_classifier) |*c| config.custom_classifier = c;
        applyWorkflowJobOverrides(&config, args, job);
        if (job.chain_map) |path| {
            loaded_chain_map = chain_map.loadFile(allocator, io, path) catch |err| {
                std.debug.print("Error loading chain map '{s}' for workflow job '{s}': {s}\n", .{ path, job.name, @errorName(err) });
                return err;
            };
            config.chain_map = &loaded_chain_map.?;
        }
        if (config.jsonl_include_atom_identity and loaded_chain_map == null) {
            std.debug.print("Error: output.jsonl.atom_identity is supported only for chain_map jobs\n", .{});
            return error.AtomIdentityRequiresSelectionMap;
        }
        try validateBitmaskCorrectionConfig(config);
        validateBatchOutputFormat(config.output_format) catch {
            std.debug.print("Error: freesasa and rsa output formats are only supported by the calc command\n", .{});
            return error.InvalidArgument;
        };
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
                if (config.jsonl_metadata == .sidecar) {
                    const metadata_path = try workflowJsonlMetadataPath(allocator, out, job.name);
                    defer allocator.free(metadata_path);
                    try writeWorkflowJsonlMetadata(allocator, io, metadata_path, job.name, config);
                }
            }
        } else if (output_dir) |out| {
            job_output_dir = try workflowPerFileOutputDir(allocator, out, job.name);
        }

        if (!config.quiet) {
            std.debug.print("Workflow job: {s}\n", .{job.name});
        }

        if (loaded_chain_map != null and config.output_format == .jsonl) {
            var selection_failures = FailureLog{ .allocator = allocator };
            defer selection_failures.deinit();
            const stats = runSelectionMapBatch(
                allocator,
                io,
                input_dir,
                config,
                jsonl_output_path,
                &loaded_chain_map.?,
                &selection_failures,
            ) catch |err| {
                std.debug.print("Error running selection-map workflow job '{s}': {s}\n", .{ job.name, @errorName(err) });
                try failed_jobs.append(allocator, job.name);
                continue;
            };
            successful += stats.successful;
            failed += stats.failed;
            if (stats.failed > 0) {
                const report = FailureReport{
                    .job = job.name,
                    .unit = .selection,
                    .total = stats.successful + stats.failed,
                    .failed = stats.failed,
                    .failures = selection_failures.sorted(),
                    .jsonl = JsonlDestination.of(config, jsonl_output_path),
                };
                try reports.append(report_arena.allocator(), try report.dupe(report_arena.allocator()));
            }
            continue;
        }
        if (loaded_chain_map) |*map| {
            if (map.hasMultipleSelections()) {
                std.debug.print("Error: chain maps with multiple selections per file require format = \"jsonl\"\n", .{});
                return error.MultiSelectionMapRequiresJsonl;
            }
        }

        // The steps of `runBatch`, taken one at a time so that a failure can
        // say what went wrong: nothing is created before the inputs are
        // accepted, and a job that cannot write its output calculates nothing.
        var prepared = prepareBatch(allocator, io, input_dir, job_output_dir, config) catch |err| {
            switch (err) {
                error.OutputNameCollision => printWorkflowJobError(job.name, "its inputs share output names (listed above)", .{}),
                error.UnsupportedChainMapInputFormat, error.OutOfMemory => printWorkflowJobError(job.name, "{s}", .{@errorName(err)}),
                else => printWorkflowJobError(job.name, "cannot read input directory '{s}': {s}", .{ input_dir, @errorName(err) }),
            }
            try failed_jobs.append(allocator, job.name);
            continue;
        };
        defer prepared.deinit(allocator);

        if (job_output_dir) |dir| {
            std.Io.Dir.cwd().createDirPath(io, dir) catch |err| {
                printWorkflowJobError(job.name, "cannot create output directory '{s}': {s}", .{ dir, @errorName(err) });
                try failed_jobs.append(allocator, job.name);
                continue;
            };
        }
        if (jsonl_output_path) |path| {
            truncateJsonlOutput(io, path) catch |err| {
                printWorkflowJobError(job.name, "cannot create JSONL output '{s}': {s}", .{ path, @errorName(err) });
                try failed_jobs.append(allocator, job.name);
                continue;
            };
        }

        var result = runPrepared(allocator, io, input_dir, job_output_dir, config, jsonl_output_path, &prepared) catch |err| {
            printWorkflowJobError(job.name, "{s}", .{@errorName(err)});
            try failed_jobs.append(allocator, job.name);
            continue;
        };
        defer result.deinit();

        successful += result.successful;
        failed += result.failed;
        if (result.failed > 0) {
            var buffer: [max_listed_failures]Failure = undefined;
            const report = result.failureReport(&buffer, job.name, JsonlDestination.of(config, jsonl_output_path));
            try reports.append(report_arena.allocator(), try report.dupe(report_arena.allocator()));
        }
    }

    std.debug.print("Workflow complete: {d} successful, {d} failed\n", .{ successful, failed });
    for (reports.items) |report| report.print(allocator);

    if (failed_jobs.items.len > 0) {
        printFailedWorkflowJobs(allocator, failed_jobs.items, workflow.jobs.len);
        return error.WorkflowJobFailed;
    }
}

fn runWorkflowFileFirst(allocator: Allocator, io: std.Io, args: BatchArgs, pre_scanned_files: ?[]const []const u8) !void {
    if (args.chain_filter != null) {
        std.debug.print("Error: --workflow cannot be combined with --chain; --manifest is a compatibility alias for --workflow; use [[jobs]].chains in the workflow\n", .{});
        return error.InvalidArgument;
    }

    const workflow_path = args.workflow_path.?;
    var workflow = parseWorkflowFile(allocator, io, workflow_path) catch |err| {
        printWorkflowReadError(workflow_path, err);
        return err;
    };
    defer workflow.deinit();

    try batchWorkflowFindings(workflow, args).report();
    try workflow.applyInputChainToJobs();

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
    try validateBitmaskCorrectionConfig(resource_config);
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

    resource_config.external_ccd = if (ext_ccd != null) &ext_ccd.? else null;
    resource_config.sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null;
    if (custom_classifier) |*c| resource_config.custom_classifier = c;

    const files = pre_scanned_files orelse try scanDirectory(allocator, io, input_dir);
    defer if (pre_scanned_files == null) freeScannedFiles(allocator, files);

    var states = try allocator.alloc(WorkflowJobState, workflow.jobs.len);
    var states_initialized: usize = 0;
    defer {
        for (states[0..states_initialized]) |*state| state.deinit(allocator);
        allocator.free(states);
    }
    for (workflow.jobs, 0..) |job, i| {
        var config = BatchConfig{};
        try applyWorkflowToBatchConfig(&config, args, workflow.calculation, workflow.output, workflow.classifier);
        applyCliOverrides(&config, args);
        config.store_atom_areas = batchShouldStoreAtomAreas(config);
        config.external_ccd = if (ext_ccd != null) &ext_ccd.? else null;
        config.sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null;
        if (custom_classifier) |*c| config.custom_classifier = c;
        applyWorkflowJobOverrides(&config, args, job);
        if (config.jsonl_include_atom_identity) {
            std.debug.print("Error: output.jsonl.atom_identity is supported only for chain_map jobs\n", .{});
            return error.AtomIdentityRequiresSelectionMap;
        }
        try validateBitmaskCorrectionConfig(config);
        validateBatchOutputFormat(config.output_format) catch {
            std.debug.print("Error: freesasa and rsa output formats are only supported by the calc command\n", .{});
            return error.InvalidArgument;
        };
        validateResidueMapFormat(config.output_format, config.residue_map) catch {
            std.debug.print("Error: residue_map is only supported with format = \"jsonl\"\n", .{});
            return error.InvalidArgument;
        };

        if ((config.classifier_type == .ccd or config.classifier_type == .protor) and config.include_hydrogens and !config.quiet) {
            std.debug.print("Warning: --include-hydrogens with CCD classifier may give inaccurate results\n", .{});
            std.debug.print("         CCD uses united-atom radii that already account for implicit hydrogens\n", .{});
        }

        var job_output_dir: ?[]const u8 = null;
        var jsonl_output_path: ?[]const u8 = null;
        errdefer if (job_output_dir) |path| allocator.free(path);
        errdefer if (jsonl_output_path) |path| allocator.free(path);

        if (config.output_format == .jsonl) {
            if (output_dir) |out| {
                try std.Io.Dir.cwd().createDirPath(io, out);
                jsonl_output_path = try workflowJsonlOutputPath(allocator, out, job.name);
                try truncateJsonlOutput(io, jsonl_output_path.?);
                if (config.jsonl_metadata == .sidecar) {
                    const metadata_path = try workflowJsonlMetadataPath(allocator, out, job.name);
                    defer allocator.free(metadata_path);
                    try writeWorkflowJsonlMetadata(allocator, io, metadata_path, job.name, config);
                }
            }
        } else if (output_dir) |out| {
            try validateUniqueFileOutputNames(allocator, files, out, config);
            job_output_dir = try workflowPerFileOutputDir(allocator, out, job.name);
            try std.Io.Dir.cwd().createDirPath(io, job_output_dir.?);
        }

        states[i] = .{
            .name = job.name,
            .config = config,
            .output_dir = job_output_dir,
            .jsonl_output_path = jsonl_output_path,
            .failures = .{ .allocator = allocator },
        };
        states_initialized += 1;
    }

    for (states) |state| {
        if (!state.config.quiet) std.debug.print("Workflow job: {s}\n", .{state.name});
    }

    var luts = try BatchLuts.init(allocator, resource_config);
    defer luts.deinit();

    const file_threads = effectiveWorkflowFileThreads(resource_config, files.len);
    if (file_threads > 1 and files.len > 1) {
        var runtimes = try allocator.alloc(WorkflowJobRuntime, states.len);
        var runtimes_initialized: usize = 0;
        defer {
            for (runtimes[0..runtimes_initialized]) |*runtime| runtime.close(io);
            allocator.free(runtimes);
        }
        for (states, 0..) |*state, i| {
            runtimes[i] = .{ .state = state };
            if (state.config.output_format == .jsonl) {
                if (state.jsonl_output_path) |path| {
                    const file = try std.Io.Dir.cwd().createFile(io, path, .{});
                    runtimes[i].jsonl_file = file;
                    runtimes[i].jsonl_file_needs_close = true;
                    runtimes[i].jsonl_stream = JsonlStreamWriter.init(file, io, jsonlOptions(state.config), &runtimes[i].jsonl_buffer);
                } else {
                    const file = std.Io.File.stdout();
                    runtimes[i].jsonl_file = file;
                    runtimes[i].jsonl_stream = JsonlStreamWriter.init(file, io, jsonlOptions(state.config), &runtimes[i].jsonl_buffer);
                }
            }
            runtimes_initialized += 1;
        }

        var ctx = WorkflowParallelContext{
            .files = files,
            .input_dir = input_dir,
            .jobs = workflow.jobs,
            .runtimes = runtimes,
            .resource_config = resource_config,
            .result_allocator = allocator,
            .next_file = std.atomic.Value(usize).init(0),
            .processed_count = std.atomic.Value(usize).init(0),
            .lut_f64 = luts.f64Ptr(),
            .lut_f32 = luts.f32Ptr(),
            .coarse_lut_f64 = luts.coarseF64Ptr(),
            .fine_lut_f64 = luts.fineF64Ptr(),
            .coarse_lut_f32 = luts.coarseF32Ptr(),
            .fine_lut_f32 = luts.fineF32Ptr(),
            .io = io,
        };

        const threads = try allocator.alloc(std.Thread, file_threads);
        defer allocator.free(threads);
        {
            // Scoped so a later error cannot join the same threads again.
            var spawned_count: usize = 0;
            errdefer joinSpawnedThreads(threads, spawned_count);
            for (threads) |*thread| {
                thread.* = try std.Thread.spawn(.{}, workflowParallelWorker, .{&ctx});
                spawned_count += 1;
            }
            joinSpawnedThreads(threads, spawned_count);
        }

        for (runtimes) |*runtime| {
            runtime.state.successful = runtime.counter.successful.load(.monotonic);
            runtime.state.failed = runtime.counter.failed.load(.monotonic);
            runtime.state.total_sasa_time_ns = runtime.counter.total_sasa_time_ns.load(.monotonic);
            if (runtime.jsonl_stream) |*stream| {
                try stream.flush();
                if (stream.hasError()) return error.JsonlWriteFailed;
            }
        }
        printWorkflowJobStates(allocator, states);
        return;
    }

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    for (files) |filename| {
        const input_path = try std.fs.path.join(arena.allocator(), &.{ input_dir, filename });

        var source_config = resource_config;
        source_config.chain_filter = null;
        var source_parsed = readInputFile(arena.allocator(), io, input_path, source_config) catch |err| {
            try workflowFailAllJobsSequential(io, arena.allocator(), states, filename, "read/parse", err);
            _ = arena.reset(.retain_capacity);
            continue;
        };

        workflowClassifySourceInput(&source_parsed.input, source_parsed.inlineCcdPtr(), input_path, source_config) catch |err| {
            try workflowFailAllJobsSequential(io, arena.allocator(), states, filename, "classifier", err);
            source_parsed.deinit();
            _ = arena.reset(.retain_capacity);
            continue;
        };

        for (workflow.jobs, 0..) |job, job_index| {
            var state = &states[job_index];

            const selected_chains: ?[]const []const u8 = if (format_detect.detectInputFormat(input_path) == .json) null else job.chains;
            var selected_input = copySelectedAtomInput(arena.allocator(), source_parsed.input, selected_chains) catch |err| {
                const message = std.fmt.allocPrint(arena.allocator(), "selection failed: {s}", .{@errorName(err)}) catch "selection failed";
                try workflowFailJobSequential(io, arena.allocator(), state, filename, message);
                continue;
            };
            defer selected_input.deinit();

            const n_threads = effectiveWorkflowSasaThreads(state.config);
            var result = switch (state.config.precision) {
                .f64 => calculatePreparedInputResult(f64, arena.allocator(), io, allocator, selected_input, state.output_dir, filename, state.config, n_threads, luts.f64Ptr(), luts.coarseF64Ptr(), luts.fineF64Ptr()),
                .f32 => calculatePreparedInputResult(f32, arena.allocator(), io, allocator, selected_input, state.output_dir, filename, state.config, n_threads, luts.f32Ptr(), luts.coarseF32Ptr(), luts.fineF32Ptr()),
            };
            defer if (result.error_msg) |msg| allocator.free(msg);

            if (result.status == .ok) {
                state.successful += 1;
                state.total_sasa_time_ns += result.sasa_time_ns;
            } else {
                state.failed += 1;
                state.failures.record(io, filename, null, result.error_msg orelse "unknown error");
            }

            try workflowWriteSequentialJsonl(io, arena.allocator(), state, &result);

            result.atom_areas = null;
            result.residue_map = null;
        }

        source_parsed.deinit();
        _ = arena.reset(.retain_capacity);
    }

    printWorkflowJobStates(allocator, states);
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

    var chain_filter_slice: ?[]const []const u8 = null;
    if (args.chain_filter) |filter_str| {
        chain_filter_slice = parseBatchChainFilter(allocator, filter_str) catch |err| switch (err) {
            error.EmptyChainFilter => {
                std.debug.print("Error: --chain needs at least one chain ID (for example --chain=A or --chain=A,B), got '{s}'\n", .{filter_str});
                return error.InvalidArgument;
            },
            else => |e| return e,
        };
    }
    defer if (chain_filter_slice) |s| allocator.free(s);

    // Load external CCD dictionary if specified
    var ext_ccd: ?ccd_parser.ComponentDict = null;
    const use_ccd_resources = batchArgsUseCcdResources(args);
    if (use_ccd_resources) {
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
    }
    defer if (ext_ccd) |*d| d.deinit();

    // Load SDF components from --sdf option
    var sdf_ccd: ?ccd_parser.ComponentDict = null;
    if (use_ccd_resources and args.sdf_paths.len > 0) {
        sdf_ccd = loadSdfComponents(allocator, io, args.sdf_paths.constSlice(), args.quiet) catch |err| {
            std.debug.print("Error loading SDF components: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    }
    defer if (sdf_ccd) |*d| d.deinit();

    validateBatchOutputFormat(args.output_format) catch {
        std.debug.print("Error: freesasa and rsa output formats are only supported by the calc command\n", .{});
        return error.InvalidArgument;
    };
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
    if (args.bitmask_correction and !args.use_bitmask) {
        std.debug.print("Error: --bitmask-correction requires --use-bitmask\n", .{});
        return error.InvalidArgument;
    }
    if (args.bitmask_correction and args.adaptive_sr) {
        std.debug.print("Error: --bitmask-correction is not supported with --adaptive-sr\n", .{});
        return error.InvalidArgument;
    }

    // For jsonl, don't pass output_dir to runBatch (no per-file I/O during computation)
    const output_dir: ?[]const u8 = if (args.output_format == .jsonl) null else args.output_path;

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
        .lr_trig = args.lr_trig,
        .probe_radius = args.probe_radius,
        .precision = args.precision,
        .output_format = args.output_format,
        .show_timing = args.show_timing,
        .profile_stages = args.profile_stages,
        .quiet = args.quiet,
        .show_progress = args.show_progress,
        .classifier_type = args.classifier_type,
        .include_hydrogens = args.include_hydrogens,
        .include_hetatm = args.include_hetatm,
        .use_bitmask = args.use_bitmask,
        .bitmask_correction = args.bitmask_correction,
        .bitmask_correction_coeff = args.bitmask_correction_coeff,
        .adaptive_sr = args.adaptive_sr,
        .coarse_points = args.coarse_points,
        .fine_points = args.fine_points,
        .adaptive_low = args.adaptive_low,
        .adaptive_high = args.adaptive_high,
        .store_atom_areas = batchShouldStoreAtomAreas(.{ .output_format = args.output_format }),
        .external_ccd = if (ext_ccd != null) &ext_ccd.? else null,
        .sdf_ccd = if (sdf_ccd != null) &sdf_ccd.? else null,
        .chain_filter = chain_filter_slice,
        .use_auth_chain = args.use_auth_chain,
        .alt_loc_mode = args.alt_loc_mode,
        .alt_loc_id = args.alt_loc_id,
        .af_model_fast = args.af_model_fast,
        .residue_map = args.residue_map,
        .jsonl_decimals = args.jsonl_decimals,
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

    // Print results. Quiet mode drops the summary, not the failed inputs.
    const jsonl_destination = JsonlDestination.of(config, jsonl_output_path);
    if (!args.quiet) {
        result.printSummary(args.show_timing, jsonl_destination);
    } else {
        result.printFailures(jsonl_destination);
    }

    // Always print benchmark output (for script parsing)
    if (args.show_timing) {
        std.debug.print("\n", .{});
        result.printBenchmarkOutput();
        if (args.profile_stages) {
            result.printStageProfile();
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

test "scanDirectory returns supported structure files in sorted order" {
    const allocator = std.testing.allocator;
    const io = std.testing.io;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // Written in a scrambled order on purpose; unsupported names and a
    // directory that merely looks like a structure file must be skipped.
    const names = [_][]const u8{
        "e.sdf",
        "c.json",
        "readme.md",
        "b.pdb",
        "x.xyz",
        "d.bcif.zst",
        "UPPER.PDB",
        "a.cif.gz",
        "structure.pdb.bak",
    };
    for (names) |name| {
        try tmp_dir.dir.writeFile(io, .{ .sub_path = name, .data = "" });
    }
    try tmp_dir.dir.createDir(io, "subdir.pdb", .default_dir);

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(io, &root_buf);

    const files = try scanDirectory(allocator, io, root_buf[0..root_len]);
    defer {
        for (files) |f| allocator.free(f);
        allocator.free(files);
    }

    const expected = [_][]const u8{ "UPPER.PDB", "a.cif.gz", "b.pdb", "c.json", "d.bcif.zst", "e.sdf" };
    try std.testing.expectEqual(expected.len, files.len);
    for (expected, files) |want, got| {
        try std.testing.expectEqualStrings(want, got);
    }
}

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

test "copySelectedAtomInput prefers extended chain IDs over truncated chain IDs" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "ABCD", "ABCD" });
    defer input.deinit();
    const chain_id_full = try allocator.alloc([]const u8, 2);
    chain_id_full[0] = try allocator.dupe(u8, "ABCDE");
    chain_id_full[1] = try allocator.dupe(u8, "ABCD");
    input.chain_id_full = chain_id_full;

    var selected_long = try copySelectedAtomInput(allocator, input, &.{"ABCDE"});
    defer selected_long.deinit();
    try std.testing.expectEqual(@as(usize, 1), selected_long.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), selected_long.x[0], 1e-12);
    try std.testing.expectEqualStrings("ABCDE", selected_long.chain_id_full.?[0]);

    var selected_prefix = try copySelectedAtomInput(allocator, input, &.{"ABCD"});
    defer selected_prefix.deinit();
    try std.testing.expectEqual(@as(usize, 1), selected_prefix.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), selected_prefix.x[0], 1e-12);
    try std.testing.expectEqualStrings("ABCD", selected_prefix.chain_id_full.?[0]);
}

test "copySelectedAtomInput keeps extended chain IDs lossless beyond fixed buffers" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "ABCD", "ABCD" });
    defer input.deinit();
    const chain_id_full = try allocator.alloc([]const u8, 2);
    chain_id_full[0] = try allocator.dupe(u8, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg");
    chain_id_full[1] = try allocator.dupe(u8, "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef");
    input.chain_id_full = chain_id_full;

    var selected_long = try copySelectedAtomInput(allocator, input, &.{"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg"});
    defer selected_long.deinit();
    try std.testing.expectEqual(@as(usize, 1), selected_long.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), selected_long.x[0], 1e-12);
    try std.testing.expectEqualStrings("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg", selected_long.chain_id_full.?[0]);

    var selected_prefix = try copySelectedAtomInput(allocator, input, &.{"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"});
    defer selected_prefix.deinit();
    try std.testing.expectEqual(@as(usize, 1), selected_prefix.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), selected_prefix.x[0], 1e-12);
    try std.testing.expectEqualStrings("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef", selected_prefix.chain_id_full.?[0]);
}

test "calculatePreparedInputResult computes total for selected input" {
    const allocator = std.testing.allocator;
    var input = try makeTestAtomInput(allocator, &.{ "A", "A" });
    defer input.deinit();

    const result = calculatePreparedInputResult(
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

test "BatchConfig default values" {
    const config = BatchConfig{};

    try std.testing.expectEqual(@as(usize, 0), config.n_threads);
    try std.testing.expectEqual(Algorithm.sr, config.algorithm);
    try std.testing.expectEqual(@as(u32, 100), config.n_points);
    try std.testing.expectEqual(@as(f64, 1.4), config.probe_radius);
}

test "BatchResult deinit frees every buffer a file result owns" {
    // A DebugAllocator reports a leak through its return value, which makes the
    // check an explicit assertion.
    var gpa: std.heap.DebugAllocator(.{}) = .init;
    const allocator = gpa.allocator();

    const results = try allocator.alloc(FileResult, 2);
    results[0] = FileResult{
        .filename = try allocator.dupe(u8, "ok.json"),
        .n_atoms = 100,
        .sasa_time_ns = 1000000,
        .total_sasa = 123.45,
        .status = .ok,
        .atom_areas = try allocator.alloc(f64, 3),
    };
    results[1] = FileResult{
        .filename = try allocator.dupe(u8, "bad.json"),
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .err,
        .error_msg = try allocator.dupe(u8, "read/parse failed"),
    };

    var batch_result = BatchResult{
        .total_files = 2,
        .successful = 1,
        .failed = 1,
        .total_sasa_time_ns = 1000000,
        .total_time_ns = 2000000,
        .file_results = results,
        .allocator = allocator,
    };
    batch_result.deinit();

    try std.testing.expectEqual(std.heap.Check.ok, gpa.deinit());
}

test "BatchResult phase timing fields default to zero" {
    const result = BatchResult{
        .total_files = 0,
        .successful = 0,
        .failed = 0,
        .total_sasa_time_ns = 0,
        .total_time_ns = 0,
        .file_results = &.{},
        .allocator = std.testing.allocator,
    };

    try std.testing.expectEqual(@as(u64, 0), result.scan_time_ns);
    try std.testing.expectEqual(@as(u64, 0), result.build_items_time_ns);
    try std.testing.expectEqual(@as(u64, 0), result.process_time_ns);
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
        "zsasa",            "batch",
        "--algorithm=lr",   "--threads=4",
        "--quiet",          "--timing",
        "--profile-stages", "input_dir/",
    };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input_dir/", parsed.input_path.?);
    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqual(true, parsed.quiet);
    try std.testing.expectEqual(false, parsed.show_progress);
    try std.testing.expectEqual(true, parsed.show_timing);
    try std.testing.expectEqual(true, parsed.profile_stages);
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

test "BatchArgs --altloc modes" {
    const none_args = [_][]const u8{ "zsasa", "batch", "--altloc=none", "input_dir/" };
    const all_args = [_][]const u8{ "zsasa", "batch", "--altloc", "all", "input_dir/" };
    const selected_args = [_][]const u8{ "zsasa", "batch", "--altloc=C", "input_dir/" };
    const occupancy_args = [_][]const u8{ "zsasa", "batch", "--altloc=highest-occupancy", "input_dir/" };

    const none = parseArgs(&none_args, 2);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.none, none.alt_loc_mode);

    const all = parseArgs(&all_args, 2);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.all, all.alt_loc_mode);

    const selected = parseArgs(&selected_args, 2);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.selected, selected.alt_loc_mode);
    try std.testing.expectEqual(@as(u8, 'C'), selected.alt_loc_id);

    const occupancy = parseArgs(&occupancy_args, 2);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.highest_occupancy, occupancy.alt_loc_mode);
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

test "batch rejects a --chain value without any chain ID" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    for ([_][]const u8{ "", ",", " , ,", " " }) |value| {
        try std.testing.expectError(error.EmptyChainFilter, parseBatchChainFilter(allocator, value));
    }

    // The command stops before it scans the input directory: every input
    // would otherwise fail to parse with no atom selected.
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root = root_buf[0..try tmp_dir.dir.realPath(std.testing.io, &root_buf)];
    const output_dir = try std.fs.path.join(allocator, &.{ root, "out" });
    defer allocator.free(output_dir);

    for ([_][]const u8{ "--chain=,", "--chain=" }) |flag| {
        const argv = [_][]const u8{ "zsasa", "batch", "-q", flag, root, output_dir };
        try std.testing.expectError(error.InvalidArgument, run(allocator, std.testing.io, parseArgs(&argv, 2)));
    }
    try std.testing.expectError(error.FileNotFound, std.Io.Dir.cwd().access(std.testing.io, output_dir, .{}));
}

test "workflow rejects a job with an empty chains array before running anything" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("tiny.pdb", test_naming_pdb);
    const workflow_path = try sandbox.path("workflow.toml");
    defer allocator.free(workflow_path);
    const output_dir = try sandbox.path("output");
    defer allocator.free(output_dir);

    // With and without the job option that selects the job-first runner: an
    // empty list used to mean every chain in one runner and none in the other.
    inline for (.{ "", "auth_chain = true\n" }) |job_option| {
        const workflow = try std.fmt.allocPrint(allocator,
            \\version = 1
            \\kind = "workflow"
            \\
            \\[input]
            \\dir = "{s}"
            \\
            \\[output]
            \\dir = "{s}"
            \\format = "jsonl"
            \\
            \\[calculation]
            \\quiet = true
            \\
            \\[[jobs]]
            \\name = "nothing"
            \\chains = []
            \\{s}
        , .{ sandbox.input_dir, output_dir, job_option });
        defer allocator.free(workflow);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

        try std.testing.expectError(
            error.EmptyJobChains,
            runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path }),
        );
    }
    try sandbox.expectTree(&.{ "input/", "input/tiny.pdb", "workflow.toml" });
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

test "validateBatchOutputFormat rejects single-calc compatibility formats" {
    try validateBatchOutputFormat(.json);
    try validateBatchOutputFormat(.compact);
    try validateBatchOutputFormat(.csv);
    try validateBatchOutputFormat(.jsonl);
    try std.testing.expectError(error.InvalidArgument, validateBatchOutputFormat(.freesasa));
    try std.testing.expectError(error.InvalidArgument, validateBatchOutputFormat(.rsa));
}

test "public batch runners reject single-calc compatibility formats before scanning" {
    try std.testing.expectError(error.InvalidArgument, runBatch(
        std.testing.allocator,
        std.testing.io,
        "/definitely/missing/zsasa/input",
        null,
        .{ .output_format = .freesasa },
        null,
    ));
    try std.testing.expectError(error.InvalidArgument, runBatchSequential(
        std.testing.allocator,
        std.testing.io,
        "/definitely/missing/zsasa/input",
        null,
        .{ .output_format = .rsa },
        null,
    ));
    try std.testing.expectError(error.InvalidArgument, runBatchParallel(
        std.testing.allocator,
        std.testing.io,
        "/definitely/missing/zsasa/input",
        null,
        .{ .output_format = .freesasa },
        null,
    ));
}

test "findOutputNameCollisions groups inputs that share an output name" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const files = [_][]const u8{
        "1crn.cif.gz",
        "1crn.pdb",
        "1ubq.cif",
        "1ubq.cif.zst",
        "1ubq.pdb",
        "3hhb.cif.gz",
    };
    const items = try plainWorkItems(arena.allocator(), &files);
    const collisions = try findOutputNameCollisions(arena.allocator(), items, .{});

    try std.testing.expectEqual(@as(usize, 5), collisions.len);
    const expected = [_][2][]const u8{
        .{ "1crn.json", "1crn.cif.gz" },
        .{ "1crn.json", "1crn.pdb" },
        .{ "1ubq.json", "1ubq.cif" },
        .{ "1ubq.json", "1ubq.cif.zst" },
        .{ "1ubq.json", "1ubq.pdb" },
    };
    for (expected, collisions) |want, got| {
        try std.testing.expectEqualStrings(want[0], got.output_name);
        try std.testing.expectEqualStrings(want[1], got.filename);
    }
}

test "findOutputNameCollisions uses the output format extension" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const files = [_][]const u8{ "1crn.ent", "1crn.pdb" };
    const items = try plainWorkItems(arena.allocator(), &files);
    const collisions = try findOutputNameCollisions(arena.allocator(), items, .{ .output_format = .csv });

    try std.testing.expectEqual(@as(usize, 2), collisions.len);
    try std.testing.expectEqualStrings("1crn.csv", collisions[0].output_name);
}

test "findOutputNameCollisions accepts distinct stems" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    // "1crn.v2.pdb" keeps its inner dot ("1crn.v2.json"), so it does not
    // collide with "1crn.pdb".
    const files = [_][]const u8{ "1crn.pdb", "1crn.v2.pdb", "1ubq.cif.gz", "3hhb.json" };
    const items = try plainWorkItems(arena.allocator(), &files);
    const collisions = try findOutputNameCollisions(arena.allocator(), items, .{});

    try std.testing.expectEqual(@as(usize, 0), collisions.len);
}

test "findOutputNameCollisions compares names without regard to ASCII case" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const files = [_][]const u8{ "PROT.pdb", "other.pdb", "prot.cif.gz", "x.cif", "x.ent", "X.pdb" };
    const items = try plainWorkItems(arena.allocator(), &files);
    const collisions = try findOutputNameCollisions(arena.allocator(), items, .{});

    try std.testing.expectEqual(@as(usize, 5), collisions.len);
    const expected = [_][2][]const u8{
        .{ "PROT.json", "PROT.pdb" },
        .{ "prot.json", "prot.cif.gz" },
        .{ "X.json", "X.pdb" },
        .{ "x.json", "x.cif" },
        .{ "x.json", "x.ent" },
    };
    for (expected, collisions) |want, got| {
        try std.testing.expectEqualStrings(want[0], got.output_name);
        try std.testing.expectEqualStrings(want[1], got.filename);
    }

    // The message says when case is the only difference.
    const message = try formatOutputNameCollisions(arena.allocator(), collisions);
    try std.testing.expectEqualStrings(
        "Error: 2 output names are shared by more than one input:\n" ++
            "  PROT.json <- PROT.pdb, prot.cif.gz (the output names differ only in case)\n" ++
            "  X.json <- X.pdb, x.cif, x.ent (some of the output names differ only in case)\n" ++
            "Split these inputs into separate directories, or use JSONL output (--format=jsonl) to keep one record per input.\n",
        message,
    );

    const same = [_][]const u8{ "1crn.cif.gz", "1crn.pdb" };
    const same_items = try plainWorkItems(arena.allocator(), &same);
    const same_message = try formatOutputNameCollisions(
        arena.allocator(),
        try findOutputNameCollisions(arena.allocator(), same_items, .{}),
    );
    try std.testing.expectEqualStrings(
        "Error: 1 output name is shared by more than one input:\n" ++
            "  1crn.json <- 1crn.cif.gz, 1crn.pdb\n" ++
            "Split these inputs into separate directories, or use JSONL output (--format=jsonl) to keep one record per input.\n",
        same_message,
    );
}

test "findOutputNameCollisions compares the real names of SDF molecule outputs" {
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();

    const mol = sdf_parser.SdfMolecule{ .name = "", .atoms = &.{}, .bonds = &.{} };
    const sdf_item = struct {
        fn make(molecule: *const sdf_parser.SdfMolecule, filename: []const u8, display_name: []const u8, mol_idx: usize) WorkItem {
            return .{ .filename = filename, .display_name = display_name, .molecule = molecule, .mol_idx = mol_idx };
        }
    }.make;

    // Molecules of SDF files that share a stem only collide when their names do.
    const distinct = [_]WorkItem{
        sdf_item(&mol, "lig.mol", "lig_one", 0),
        .{ .filename = "lig.pdb", .display_name = "lig.pdb" },
        sdf_item(&mol, "lig.sdf", "lig_two", 0),
        sdf_item(&mol, "lig.sdf", "lig_three", 1),
    };
    try std.testing.expectEqual(
        @as(usize, 0),
        (try findOutputNameCollisions(arena.allocator(), &distinct, .{})).len,
    );

    // An SDF molecule output can equal a non-SDF output or the output of a
    // molecule in another SDF file. The separator in "a/b" is replaced.
    const clashing = [_]WorkItem{
        sdf_item(&mol, "lig.mol", "lig_same", 0),
        sdf_item(&mol, "lig.sdf", "lig_1", 0),
        sdf_item(&mol, "lig.sdf", "lig_Same", 1),
        sdf_item(&mol, "lig.sdf", "lig_a/b", 2),
        .{ .filename = "lig_1.pdb", .display_name = "lig_1.pdb" },
        .{ .filename = "lig_a_b.cif", .display_name = "lig_a_b.cif" },
        // An SDF file that could not be loaded writes nothing.
        .{ .filename = "lig_1.sdf", .display_name = "lig_1.sdf", .load_error = error.InvalidCountsLine },
    };
    const collisions = try findOutputNameCollisions(arena.allocator(), &clashing, .{});
    try std.testing.expectEqualStrings(
        "Error: 3 output names are shared by more than one input:\n" ++
            "  lig_1.json <- lig.sdf (molecule 1), lig_1.pdb\n" ++
            "  lig_a_b.json <- lig.sdf (molecule 3), lig_a_b.cif\n" ++
            "  lig_Same.json <- lig.sdf (molecule 2), lig.mol (molecule 1) (the output names differ only in case)\n" ++
            "Split these inputs into separate directories, or use JSONL output (--format=jsonl) to keep one record per input.\n",
        try formatOutputNameCollisions(arena.allocator(), collisions),
    );
}

fn expectSdfDisplayNames(filename: []const u8, titles: []const []const u8, expected: []const []const u8) !void {
    const allocator = std.testing.allocator;
    const molecules = try allocator.alloc(sdf_parser.SdfMolecule, titles.len);
    defer allocator.free(molecules);
    for (titles, molecules) |title, *mol| mol.* = .{ .name = title, .atoms = &.{}, .bonds = &.{} };

    const names = try sdfMoleculeDisplayNames(allocator, filename, molecules);
    defer {
        for (names) |name| allocator.free(name);
        allocator.free(names);
    }

    try std.testing.expectEqual(expected.len, names.len);
    for (expected, names) |want, got| try std.testing.expectEqualStrings(want, got);
}

test "sdfMoleculeDisplayNames keeps the names of molecules that do not clash" {
    // "stem_title", or "stem_N" for a blank title.
    try expectSdfDisplayNames(
        "two_molecules.sdf",
        &.{ "methane", "", "water" },
        &.{ "two_molecules_methane", "two_molecules_2", "two_molecules_water" },
    );
    // Only the format (and compression) extension is removed from the stem,
    // and titles are kept verbatim.
    try expectSdfDisplayNames("lig.v2.sdf.gz", &.{ "a.b", "c/d" }, &.{ "lig.v2_a.b", "lig.v2_c/d" });
    try expectSdfDisplayNames("lig.mol.zst", &.{"x"}, &.{"lig_x"});
    // A title shared by other molecules does not affect a unique one.
    try expectSdfDisplayNames("c.sdf", &.{ "x", "y", "x" }, &.{ "c_x_1", "c_y", "c_x_3" });
}

test "sdfMoleculeDisplayNames appends the position to names that clash" {
    try expectSdfDisplayNames("dup.sdf", &.{ "ethanol", "ethanol" }, &.{ "dup_ethanol_1", "dup_ethanol_2" });

    // The result must not be the name of another molecule: "_N" is appended
    // until it is not, and the molecule with the unique title keeps its name.
    try expectSdfDisplayNames("c.sdf", &.{ "x", "x", "x_2" }, &.{ "c_x_1", "c_x_2_2", "c_x_2" });
    try expectSdfDisplayNames(
        "c.sdf",
        &.{ "x", "x", "x_2", "x_2_2" },
        &.{ "c_x_1", "c_x_2_2_2", "c_x_2", "c_x_2_2" },
    );
    try expectSdfDisplayNames(
        "c.sdf",
        &.{ "x", "x", "x_1", "x_1" },
        &.{ "c_x_1", "c_x_2", "c_x_1_3", "c_x_1_4" },
    );

    // A blank title is named by position, which can be another molecule's title.
    try expectSdfDisplayNames("c.sdf", &.{ "2", "" }, &.{ "c_2_1", "c_2_2" });

    // Names that would be written to one file clash too: the same name in a
    // different case, or with a different unsafe character.
    try expectSdfDisplayNames("c.sdf", &.{ "Ethanol", "ethanol" }, &.{ "c_Ethanol_1", "c_ethanol_2" });
    try expectSdfDisplayNames(
        "c.sdf",
        &.{ "a/b", "a_b", "a\\b", "ab" },
        &.{ "c_a/b_1", "c_a_b_2", "c_a\\b_3", "c_ab" },
    );
}

test "sdfMoleculeOutputName appends the extension and replaces unsafe characters" {
    const allocator = std.testing.allocator;
    const cases = [_][2][]const u8{
        .{ "two_molecules_methane", "two_molecules_methane.json" },
        // Dots in the stem or the title are not an extension.
        .{ "lig.v2_methane", "lig.v2_methane.json" },
        .{ "lig_v1.5", "lig_v1.5.json" },
        .{ "lig_.", "lig_..json" },
        .{ "lig_..", "lig_...json" },
        // Path separators and control characters.
        .{ "lig_a/b", "lig_a_b.json" },
        .{ "lig_c\\d", "lig_c_d.json" },
        .{ "lig_../../escape", "lig_.._.._escape.json" },
        .{ "lig_/etc/passwd", "lig__etc_passwd.json" },
        .{ "lig_a\x00b\tc\x7fd\r", "lig_a_b_c_d_.json" },
        // Never "", "." or "..", whatever the display name is.
        .{ "", ".json" },
        .{ ".", "..json" },
        .{ "..", "...json" },
        .{ "/", "_.json" },
    };
    for (cases) |case| {
        const name = try sdfMoleculeOutputName(allocator, case[0], ".json");
        defer allocator.free(name);
        try std.testing.expectEqualStrings(case[1], name);
        try std.testing.expect(std.mem.findAny(u8, name, "/\\") == null);

        const via_source = try perFileOutputName(allocator, .sdf_molecule, case[0], ".json");
        defer allocator.free(via_source);
        try std.testing.expectEqualStrings(case[1], via_source);
    }

    // Input file names keep their naming: the extension is replaced.
    const file_cases = [_][2][]const u8{
        .{ "1ubq.pdb", "1ubq.csv" },
        .{ "1ubq.cif.gz", "1ubq.csv" },
        .{ "1crn.v2.pdb", "1crn.v2.csv" },
    };
    for (file_cases) |case| {
        const name = try perFileOutputName(allocator, .input_file, case[0], ".csv");
        defer allocator.free(name);
        try std.testing.expectEqualStrings(case[1], name);
    }
}

test "isUnsafeFileNameByte replaces Windows reserved characters only on Windows" {
    for ("/\\\x00\x01\n\t\x1f\x7f") |c| {
        try std.testing.expect(isUnsafeFileNameByte(c, .linux));
        try std.testing.expect(isUnsafeFileNameByte(c, .macos));
        try std.testing.expect(isUnsafeFileNameByte(c, .windows));
    }
    for ("<>:\"|?*") |c| {
        try std.testing.expect(!isUnsafeFileNameByte(c, .linux));
        try std.testing.expect(!isUnsafeFileNameByte(c, .macos));
        try std.testing.expect(isUnsafeFileNameByte(c, .windows));
    }
    for ("aZ09._- ()[]+,=@~\xc3\xa9") |c| {
        try std.testing.expect(!isUnsafeFileNameByte(c, .linux));
        try std.testing.expect(!isUnsafeFileNameByte(c, .windows));
    }
}

/// One-atom V2000 record with the given title line, for output naming tests.
/// A title of a single space is a blank title, like an empty line.
fn testSdfRecord(comptime title: []const u8) []const u8 {
    return title ++ "\n" ++
        "  zsasa\n" ++
        "\n" ++
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n" ++
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
        "M  END\n" ++
        "$$$$\n";
}

const test_naming_pdb =
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n" ++
    "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C\n" ++
    "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00 20.00           C\n" ++
    "END\n";

/// Thread counts for output naming tests: 1 runs the sequential runner, 4
/// the parallel one (for more than one work item).
const test_naming_threads = [_]usize{ 1, 4 };

/// A temporary directory with an "input" directory for output naming tests.
const NamingSandbox = struct {
    tmp: std.testing.TmpDir,
    root: []const u8,
    input_dir: []const u8,

    fn init() !NamingSandbox {
        const allocator = std.testing.allocator;
        var tmp = std.testing.tmpDir(.{ .iterate = true });
        errdefer tmp.cleanup();

        var root_buf: [std.fs.max_path_bytes]u8 = undefined;
        const root_len = try tmp.dir.realPath(std.testing.io, &root_buf);
        const root = try allocator.dupe(u8, root_buf[0..root_len]);
        errdefer allocator.free(root);
        const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
        errdefer allocator.free(input_dir);
        try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

        return .{ .tmp = tmp, .root = root, .input_dir = input_dir };
    }

    fn deinit(self: *NamingSandbox) void {
        std.testing.allocator.free(self.input_dir);
        std.testing.allocator.free(self.root);
        self.tmp.cleanup();
    }

    /// Absolute path of `name` below the sandbox root. Caller frees.
    fn path(self: NamingSandbox, name: []const u8) ![]u8 {
        return std.fs.path.join(std.testing.allocator, &.{ self.root, name });
    }

    fn writeInput(self: NamingSandbox, name: []const u8, data: []const u8) !void {
        const file_path = try std.fs.path.join(std.testing.allocator, &.{ self.input_dir, name });
        defer std.testing.allocator.free(file_path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = file_path, .data = data });
    }

    fn config(n_threads: usize) BatchConfig {
        return .{ .n_threads = n_threads, .n_points = 8, .quiet = true, .show_progress = false };
    }

    /// Run a batch over the input directory, writing per-file JSON output to
    /// the directory `out_name` below the sandbox root.
    fn run(self: NamingSandbox, n_threads: usize, out_name: []const u8) !BatchResult {
        const output_dir = try self.path(out_name);
        defer std.testing.allocator.free(output_dir);
        return runBatch(std.testing.allocator, std.testing.io, self.input_dir, output_dir, config(n_threads), null);
    }

    /// Run a batch over the input directory in JSONL mode and return the
    /// `filename` of every row, sorted. Caller frees with `freeNames`.
    fn runJsonl(self: NamingSandbox, n_threads: usize) ![][]const u8 {
        const allocator = std.testing.allocator;
        const jsonl_path = try self.path("results.jsonl");
        defer allocator.free(jsonl_path);

        var jsonl_config = config(n_threads);
        jsonl_config.output_format = .jsonl;
        jsonl_config.store_atom_areas = true;
        var result = try runBatch(allocator, std.testing.io, self.input_dir, null, jsonl_config, jsonl_path);
        defer result.deinit();

        const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(64 * 1024));
        defer allocator.free(content);
        try std.Io.Dir.cwd().deleteFile(std.testing.io, jsonl_path);

        var names = std.ArrayListUnmanaged([]const u8).empty;
        errdefer freeNames(names.items);
        var lines = std.mem.tokenizeScalar(u8, content, '\n');
        while (lines.next()) |line| {
            const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
            defer parsed.deinit();
            try std.testing.expectEqualStrings("ok", parsed.value.object.get("status").?.string);
            const name = try allocator.dupe(u8, parsed.value.object.get("filename").?.string);
            errdefer allocator.free(name);
            try names.append(allocator, name);
        }
        std.mem.sort([]const u8, names.items, {}, stringLessThan);
        return names.toOwnedSlice(allocator);
    }

    fn freeNames(names: []const []const u8) void {
        for (names) |name| std.testing.allocator.free(name);
        std.testing.allocator.free(names);
    }

    fn stringLessThan(_: void, a: []const u8, b: []const u8) bool {
        return std.mem.lessThan(u8, a, b);
    }

    /// Expect the sandbox to hold exactly `expected`: every file and
    /// directory below the root, in any order, directories with a trailing
    /// '/'. Nothing may have been written anywhere else.
    fn expectTree(self: NamingSandbox, expected: []const []const u8) !void {
        const allocator = std.testing.allocator;

        var actual = std.ArrayListUnmanaged([]const u8).empty;
        defer {
            for (actual.items) |entry| allocator.free(entry);
            actual.deinit(allocator);
        }
        var walker = try self.tmp.dir.walk(allocator);
        defer walker.deinit();
        while (try walker.next(std.testing.io)) |entry| {
            const line = try std.fmt.allocPrint(allocator, "{s}{s}", .{
                entry.path,
                if (entry.kind == .directory) "/" else "",
            });
            errdefer allocator.free(line);
            if (std.fs.path.sep != '/') std.mem.replaceScalar(u8, line, std.fs.path.sep, '/');
            try actual.append(allocator, line);
        }
        std.mem.sort([]const u8, actual.items, {}, stringLessThan);

        const wanted = try allocator.dupe([]const u8, expected);
        defer allocator.free(wanted);
        std.mem.sort([]const u8, wanted, {}, stringLessThan);

        const actual_text = try std.mem.join(allocator, "\n", actual.items);
        defer allocator.free(actual_text);
        const wanted_text = try std.mem.join(allocator, "\n", wanted);
        defer allocator.free(wanted_text);
        try std.testing.expectEqualStrings(wanted_text, actual_text);
    }
};

fn expectResultNames(result: BatchResult, expected: []const []const u8) !void {
    try std.testing.expectEqual(expected.len, result.total_files);
    try std.testing.expectEqual(expected.len, result.successful);
    try std.testing.expectEqual(@as(usize, 0), result.failed);
    for (expected, result.file_results) |want, got| {
        try std.testing.expectEqualStrings(want, got.filename);
    }
}

test "batch keeps the names of an SDF file whose titles are distinct" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("two_molecules.sdf", comptime testSdfRecord("methane") ++ testSdfRecord(" ") ++ testSdfRecord("water"));
    try sandbox.writeInput("tiny.pdb", test_naming_pdb);

    const names = [_][]const u8{ "tiny.pdb", "two_molecules_methane", "two_molecules_2", "two_molecules_water" };
    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try expectResultNames(result, &names);

        const jsonl_names = try sandbox.runJsonl(n_threads);
        defer NamingSandbox.freeNames(jsonl_names);
        try std.testing.expectEqual(names.len, jsonl_names.len);
        for ([_][]const u8{ "tiny.pdb", "two_molecules_2", "two_molecules_methane", "two_molecules_water" }, jsonl_names) |want, got| {
            try std.testing.expectEqualStrings(want, got);
        }
    }
    try sandbox.expectTree(&.{
        "input/",
        "input/tiny.pdb",
        "input/two_molecules.sdf",
        "out1/",
        "out1/tiny.json",
        "out1/two_molecules_2.json",
        "out1/two_molecules_methane.json",
        "out1/two_molecules_water.json",
        "out4/",
        "out4/tiny.json",
        "out4/two_molecules_2.json",
        "out4/two_molecules_methane.json",
        "out4/two_molecules_water.json",
    });
}

test "batch gives SDF molecules with the same title their own output and JSONL name" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("dup.sdf", comptime testSdfRecord("ethanol") ++ testSdfRecord("ethanol"));
    try sandbox.writeInput("c.sdf", comptime testSdfRecord("x") ++ testSdfRecord("x") ++ testSdfRecord("x_2"));

    const names = [_][]const u8{ "c_x_1", "c_x_2_2", "c_x_2", "dup_ethanol_1", "dup_ethanol_2" };
    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try expectResultNames(result, &names);

        const jsonl_names = try sandbox.runJsonl(n_threads);
        defer NamingSandbox.freeNames(jsonl_names);
        try std.testing.expectEqual(names.len, jsonl_names.len);
        for ([_][]const u8{ "c_x_1", "c_x_2", "c_x_2_2", "dup_ethanol_1", "dup_ethanol_2" }, jsonl_names) |want, got| {
            try std.testing.expectEqualStrings(want, got);
        }
    }
    try sandbox.expectTree(&.{
        "input/",
        "input/c.sdf",
        "input/dup.sdf",
        "out1/",
        "out1/c_x_1.json",
        "out1/c_x_2.json",
        "out1/c_x_2_2.json",
        "out1/dup_ethanol_1.json",
        "out1/dup_ethanol_2.json",
        "out4/",
        "out4/c_x_1.json",
        "out4/c_x_2.json",
        "out4/c_x_2_2.json",
        "out4/dup_ethanol_1.json",
        "out4/dup_ethanol_2.json",
    });
}

test "batch keeps dots of SDF stems and titles in output names" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("lig.v2.sdf", comptime testSdfRecord("methane") ++ testSdfRecord("water") ++ testSdfRecord("v1.5"));
    try sandbox.writeInput("lig.sdf", comptime testSdfRecord("v1.5") ++ testSdfRecord("v1.6"));

    const names = [_][]const u8{ "lig_v1.5", "lig_v1.6", "lig.v2_methane", "lig.v2_water", "lig.v2_v1.5" };
    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try expectResultNames(result, &names);
    }
    try sandbox.expectTree(&.{
        "input/",
        "input/lig.sdf",
        "input/lig.v2.sdf",
        "out1/",
        "out1/lig_v1.5.json",
        "out1/lig_v1.6.json",
        "out1/lig.v2_methane.json",
        "out1/lig.v2_water.json",
        "out1/lig.v2_v1.5.json",
        "out4/",
        "out4/lig_v1.5.json",
        "out4/lig_v1.6.json",
        "out4/lig.v2_methane.json",
        "out4/lig.v2_water.json",
        "out4/lig.v2_v1.5.json",
    });
}

test "batch never writes an SDF title with path separators outside the output directory" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("lig.sdf", comptime testSdfRecord("a/b") ++ testSdfRecord("c\\d") ++ testSdfRecord("..") ++
        testSdfRecord("../../escape") ++ testSdfRecord(".") ++ testSdfRecord("/abs") ++ testSdfRecord("../input/lig.sdf"));

    // Results and JSONL rows keep the title as written.
    const names = [_][]const u8{ "lig_a/b", "lig_c\\d", "lig_..", "lig_../../escape", "lig_.", "lig_/abs", "lig_../input/lig.sdf" };
    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try expectResultNames(result, &names);

        const jsonl_names = try sandbox.runJsonl(n_threads);
        defer NamingSandbox.freeNames(jsonl_names);
        try std.testing.expectEqual(names.len, jsonl_names.len);
        for ([_][]const u8{ "lig_.", "lig_..", "lig_../../escape", "lig_../input/lig.sdf", "lig_/abs", "lig_a/b", "lig_c\\d" }, jsonl_names) |want, got| {
            try std.testing.expectEqualStrings(want, got);
        }
    }

    // The input file is intact and every output is a file directly in its
    // output directory.
    const input_path = try sandbox.path("input/lig.sdf");
    defer std.testing.allocator.free(input_path);
    const input_content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, input_path, std.testing.allocator, .limited(64 * 1024));
    defer std.testing.allocator.free(input_content);
    try std.testing.expect(std.mem.startsWith(u8, input_content, "a/b\n"));

    try sandbox.expectTree(&.{
        "input/",
        "input/lig.sdf",
        "out1/",
        "out1/lig_a_b.json",
        "out1/lig_c_d.json",
        "out1/lig_...json",
        "out1/lig_.._.._escape.json",
        "out1/lig_..json",
        "out1/lig__abs.json",
        "out1/lig_.._input_lig.sdf.json",
        "out4/",
        "out4/lig_a_b.json",
        "out4/lig_c_d.json",
        "out4/lig_...json",
        "out4/lig_.._.._escape.json",
        "out4/lig_..json",
        "out4/lig__abs.json",
        "out4/lig_.._input_lig.sdf.json",
    });
}

test "batch rejects output names that differ only in case" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("PROT.pdb", test_naming_pdb);
    try sandbox.writeInput("prot.ent", test_naming_pdb);

    inline for (test_naming_threads) |n_threads| {
        try std.testing.expectError(error.OutputNameCollision, sandbox.run(n_threads, "out"));

        // JSONL output keeps one record per input.
        const jsonl_names = try sandbox.runJsonl(n_threads);
        defer NamingSandbox.freeNames(jsonl_names);
        try std.testing.expectEqual(@as(usize, 2), jsonl_names.len);
        try std.testing.expectEqualStrings("PROT.pdb", jsonl_names[0]);
        try std.testing.expectEqualStrings("prot.ent", jsonl_names[1]);
    }
    try sandbox.expectTree(&.{ "input/", "input/PROT.pdb", "input/prot.ent" });
}

test "batch rejects an SDF molecule output that equals a non-SDF output" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    // The unnamed molecule is written to "lig_1.json", and so is "lig_1.pdb".
    try sandbox.writeInput("lig.sdf", comptime testSdfRecord(" ") ++ testSdfRecord("water"));
    try sandbox.writeInput("lig_1.pdb", test_naming_pdb);

    inline for (test_naming_threads) |n_threads| {
        try std.testing.expectError(error.OutputNameCollision, sandbox.run(n_threads, "out"));
    }
    try sandbox.expectTree(&.{ "input/", "input/lig.sdf", "input/lig_1.pdb" });
}

test "batch accepts SDF files that share a stem unless their molecules clash" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("lig.sdf", comptime testSdfRecord("one") ++ testSdfRecord("two"));
    try sandbox.writeInput("lig.mol", comptime testSdfRecord("three"));

    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try expectResultNames(result, &.{ "lig_three", "lig_one", "lig_two" });
    }
    try sandbox.expectTree(&.{
        "input/",
        "input/lig.mol",
        "input/lig.sdf",
        "out1/",
        "out1/lig_one.json",
        "out1/lig_three.json",
        "out1/lig_two.json",
        "out4/",
        "out4/lig_one.json",
        "out4/lig_three.json",
        "out4/lig_two.json",
    });

    // "lig.mol" now holds a molecule named like one of "lig.sdf".
    try sandbox.writeInput("lig.mol", comptime testSdfRecord("two"));
    inline for (test_naming_threads) |n_threads| {
        try std.testing.expectError(error.OutputNameCollision, sandbox.run(n_threads, "rejected"));
    }
    try std.testing.expectError(
        error.FileNotFound,
        sandbox.tmp.dir.access(std.testing.io, "rejected", .{}),
    );
}

test "batch reports an SDF file that cannot be parsed and lets it claim no output name" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    // "bad.sdf" would be expanded to "bad_1.json" if it could be parsed.
    try sandbox.writeInput("bad.sdf", "broken\n\n\nnot a counts line\n");
    try sandbox.writeInput("bad_1.pdb", test_naming_pdb);
    try sandbox.writeInput("lig.sdf", comptime testSdfRecord("one"));

    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 3), result.total_files);
        try std.testing.expectEqual(@as(usize, 2), result.successful);
        try std.testing.expectEqual(@as(usize, 1), result.failed);
        try std.testing.expectEqualStrings("bad.sdf", result.file_results[0].filename);
        try std.testing.expectEqualStrings("read/parse failed: InvalidCountsLine", result.file_results[0].error_msg.?);
        try std.testing.expectEqualStrings("bad_1.pdb", result.file_results[1].filename);
        try std.testing.expectEqualStrings("lig_one", result.file_results[2].filename);
    }
    try sandbox.expectTree(&.{
        "input/",
        "input/bad.sdf",
        "input/bad_1.pdb",
        "input/lig.sdf",
        "out1/",
        "out1/bad_1.json",
        "out1/lig_one.json",
        "out4/",
        "out4/bad_1.json",
        "out4/lig_one.json",
    });
}

test "workflow job writes per-molecule SDF outputs and rejects collisions" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("dup.sdf", comptime testSdfRecord("x") ++ testSdfRecord("x") ++ testSdfRecord("a/b"));
    try sandbox.writeInput("lig.v2.sdf", comptime testSdfRecord("v1.5"));

    const output_dir = try sandbox.path("output");
    defer allocator.free(output_dir);
    const workflow_path = try sandbox.path("workflow.toml");
    defer allocator.free(workflow_path);
    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "json"
        \\
        \\[calculation]
        \\n_points = 8
        \\quiet = true
        \\
        \\[[jobs]]
        \\name = "all"
        \\
    , .{ sandbox.input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    inline for (test_naming_threads) |n_threads| {
        try runWorkflow(allocator, std.testing.io, .{
            .workflow_path = workflow_path,
            .n_threads = n_threads,
            .threads_explicit = true,
        });
        try sandbox.expectTree(&.{
            "input/",
            "input/dup.sdf",
            "input/lig.v2.sdf",
            "output/",
            "output/all/",
            "output/all/dup_a_b.json",
            "output/all/dup_x_1.json",
            "output/all/dup_x_2.json",
            "output/all/lig.v2_v1.5.json",
            "workflow.toml",
        });
        try sandbox.tmp.dir.deleteTree(std.testing.io, "output");
    }

    // "dup_x_1.pdb" claims the output of the first molecule of "dup.sdf": the
    // job is rejected before its output directory is created, and the run
    // is an error.
    try sandbox.writeInput("dup_x_1.pdb", test_naming_pdb);
    inline for (test_naming_threads) |n_threads| {
        try std.testing.expectError(error.WorkflowJobFailed, runWorkflow(allocator, std.testing.io, .{
            .workflow_path = workflow_path,
            .n_threads = n_threads,
            .threads_explicit = true,
        }));
        try sandbox.expectTree(&.{
            "input/",
            "input/dup.sdf",
            "input/dup_x_1.pdb",
            "input/lig.v2.sdf",
            "workflow.toml",
        });
    }
}

test "workflow rejects per-file output names that differ only in case" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("PROT.pdb", test_naming_pdb);
    try sandbox.writeInput("prot.ent", test_naming_pdb);

    const output_dir = try sandbox.path("output");
    defer allocator.free(output_dir);
    const workflow_path = try sandbox.path("workflow.toml");
    defer allocator.free(workflow_path);
    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "json"
        \\
        \\[calculation]
        \\n_points = 8
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "chain_a"
        \\chains = ["A"]
        \\
    , .{ sandbox.input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try std.testing.expectError(
        error.OutputNameCollision,
        runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path }),
    );
    try sandbox.expectTree(&.{ "input/", "input/PROT.pdb", "input/prot.ent", "workflow.toml" });
}

test "validateUniqueOutputNames only applies to per-file output" {
    const files = [_][]const u8{ "1crn.cif", "1crn.pdb" };
    const items = try plainWorkItems(std.testing.allocator, &files);
    defer std.testing.allocator.free(items);

    // No output directory: nothing is written per file.
    try validateUniqueOutputNames(std.testing.allocator, items, null, .{});
    try validateUniqueFileOutputNames(std.testing.allocator, &files, null, .{});
    // JSONL keeps one record per input.
    const jsonl_config = BatchConfig{ .output_format = .jsonl, .store_atom_areas = true };
    try validateUniqueOutputNames(std.testing.allocator, items, "out", jsonl_config);
    try validateUniqueFileOutputNames(std.testing.allocator, &files, "out", jsonl_config);
    // Also when its rows carry no atom areas.
    const jsonl_totals_config = BatchConfig{ .output_format = .jsonl, .jsonl_include_atom_areas = false };
    try std.testing.expect(!batchShouldStoreAtomAreas(jsonl_totals_config));
    try validateUniqueOutputNames(std.testing.allocator, items, "out", jsonl_totals_config);
    try validateUniqueFileOutputNames(std.testing.allocator, &files, "out", jsonl_totals_config);
}

test "JSONL without atom areas writes rows and no per-file outputs" {
    const allocator = std.testing.allocator;
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    // Two inputs that would share "tiny.jsonl" as a per-file output
    try sandbox.writeInput("tiny.pdb", test_naming_pdb);
    try sandbox.writeInput("tiny.ent", test_naming_pdb);
    try sandbox.writeInput("lig.sdf", comptime testSdfRecord("one"));
    const output_dir = try sandbox.path("out");
    defer allocator.free(output_dir);

    inline for (test_naming_threads) |n_threads| {
        const jsonl_path = try sandbox.path(std.fmt.comptimePrint("rows{d}.jsonl", .{n_threads}));
        defer allocator.free(jsonl_path);

        // What a workflow with `[output.jsonl] atom_areas = false` configures:
        // JSONL output, and no atom areas kept in the results.
        var config = NamingSandbox.config(n_threads);
        config.output_format = .jsonl;
        config.jsonl_include_atom_areas = false;
        config.store_atom_areas = batchShouldStoreAtomAreas(config);
        try std.testing.expect(!config.store_atom_areas);
        try std.testing.expect(batchWritesJsonl(config));

        var result = try runBatch(allocator, std.testing.io, sandbox.input_dir, output_dir, config, jsonl_path);
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 3), result.successful);

        const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(4096));
        defer allocator.free(content);
        try std.testing.expectEqual(@as(usize, 3), std.mem.count(u8, content, "\"status\":\"ok\""));
        try std.testing.expectEqual(@as(usize, 3), std.mem.count(u8, content, "\"total_area\""));
        try std.testing.expectEqual(@as(usize, 0), std.mem.count(u8, content, "\"atom_areas\""));
    }
    // The output directory holds no per-file output.
    try sandbox.expectTree(&.{
        "input/",
        "input/lig.sdf",
        "input/tiny.ent",
        "input/tiny.pdb",
        "out/",
        "rows1.jsonl",
        "rows4.jsonl",
    });
}

test "batch runners reject colliding output names before writing anything" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root_path, "output" });
    defer allocator.free(output_dir);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const pdb_data =
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n" ++
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C\n" ++
        "ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00 20.00           C\n" ++
        "END\n";
    for ([_][]const u8{ "tiny.pdb", "tiny.ent", "other.pdb" }) |filename| {
        const path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = path, .data = pdb_data });
    }

    const config = BatchConfig{
        .n_threads = 2,
        .n_points = 8,
        .quiet = true,
        .show_progress = false,
        .classifier_type = .naccess,
    };

    try std.testing.expectError(
        error.OutputNameCollision,
        runBatchSequential(allocator, std.testing.io, input_dir, output_dir, config, null),
    );
    try std.testing.expectError(
        error.OutputNameCollision,
        runBatchParallel(allocator, std.testing.io, input_dir, output_dir, config, null),
    );
    try std.testing.expectError(
        error.FileNotFound,
        std.Io.Dir.cwd().access(std.testing.io, output_dir, .{}),
    );

    // Without per-file output every input is still processed.
    var result = try runBatchParallel(allocator, std.testing.io, input_dir, null, config, null);
    defer result.deinit();
    try std.testing.expectEqual(@as(usize, 3), result.total_files);
    try std.testing.expectEqual(@as(usize, 3), result.successful);
}

test "runBatchParallel reports sequential fallback errors without freeing results twice" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, "tiny.pdb" });
    defer allocator.free(input_path);
    const jsonl_path = try std.fs.path.join(allocator, &.{ root_path, "missing-dir", "results.jsonl" });
    defer allocator.free(jsonl_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path,
        .data =
        \\ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
        \\ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
        \\END
        \\
        ,
    });

    // A single work item makes the parallel runner hand over to the
    // sequential one, which then fails to create the JSONL file.
    try std.testing.expectError(error.FileNotFound, runBatchParallel(allocator, std.testing.io, input_dir, null, .{
        .n_threads = 4,
        .n_points = 8,
        .quiet = true,
        .show_progress = false,
        .output_format = .jsonl,
        .store_atom_areas = true,
        .classifier_type = .naccess,
    }, jsonl_path));
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
    try std.testing.expect(batchArgsUseCcdResources(.{ .classifier_type = .ccd }));
    try std.testing.expect(!batchArgsUseCcdResources(.{ .classifier_type = .protor }));
    try std.testing.expect(!batchArgsUseCcdResources(.{ .classifier_type = .naccess }));
    try std.testing.expect(!batchArgsUseCcdResources(.{ .classifier_type = .oons }));

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
    try std.testing.expect(resolveWorkflowCcdPath(workflow_only_args, classifier_config, .protor) == null);
    try std.testing.expectEqual(@as(usize, 1), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], .ccd).len);
    try std.testing.expectEqual(@as(usize, 0), resolveWorkflowSdfPaths(workflow_only_args, workflow_sdf_paths[0..], .protor).len);

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

test "batch mmCIF CCD classifier uses inline CCD while ProtOr skips it" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];
    const cif_path = try std.fs.path.join(allocator, &.{ root_path, "ligand.cif" });
    defer allocator.free(cif_path);

    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = cif_path,
        .data =
        \\data_LIGAND
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.pdbx_aromatic_flag
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\LIG C1 C N N
        \\LIG N1 N N N
        \\#
        \\loop_
        \\_chem_comp_bond.comp_id
        \\_chem_comp_bond.atom_id_1
        \\_chem_comp_bond.atom_id_2
        \\_chem_comp_bond.value_order
        \\LIG C1 N1 SING
        \\#
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.group_PDB
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C C1 LIG HETATM 0.0 0.0 0.0
        \\2 N N1 LIG HETATM 2.0 0.0 0.0
        \\#
        \\
        ,
    });

    var ccd_parsed = try readInputFile(allocator, std.testing.io, cif_path, .{
        .classifier_type = .ccd,
        .include_hetatm = true,
    });
    defer ccd_parsed.deinit();
    try std.testing.expect(ccd_parsed.inlineCcdPtr() != null);
    try applyBuiltinClassifier(&ccd_parsed.input, .ccd, null, ccd_parsed.inlineCcdPtr(), null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), ccd_parsed.input.r[1], 0.001);

    var protor_parsed = try readInputFile(allocator, std.testing.io, cif_path, .{
        .classifier_type = .protor,
        .include_hetatm = true,
    });
    defer protor_parsed.deinit();
    try std.testing.expect(protor_parsed.inlineCcdPtr() == null);
    try applyBuiltinClassifier(&protor_parsed.input, .protor, null, protor_parsed.inlineCcdPtr(), null);
    try std.testing.expect(protor_parsed.input.r[1] != 1.64);
}

test "batch NACCESS and OONS take the element of unlisted atoms from the element column" {
    const allocator = std.testing.allocator;

    // Atom names that start with a two-letter element symbol
    const pdb_content =
        \\ATOM      1  HG  SER A   1      10.000   0.000   0.000  1.00 20.00           H
        \\ATOM      2 HG21 VAL A   2      20.000   0.000   0.000  1.00 20.00           H
        \\HETATM    3  NA  HEM A   3      30.000   0.000   0.000  1.00 20.00           N
        \\HETATM    4  PB  ATP A   4      40.000   0.000   0.000  1.00 20.00           P
        \\HETATM    5  CD  PCA A   5      50.000   0.000   0.000  1.00 20.00           C
        \\HETATM    6 ZN    ZN A   6      60.000   0.000   0.000  1.00 20.00          ZN
        \\END
    ;
    // H, H, N, P, C, Zn
    const expected = [_]f64{ 1.10, 1.10, 1.55, 1.80, 1.70, 1.39 };

    for ([_]ClassifierType{ .naccess, .oons }) |ct| {
        var parser = pdb_parser.PdbParser.init(allocator);
        parser.atom_only = false;
        parser.skip_hydrogens = false;
        var input = try parser.parse(pdb_content);
        defer input.deinit();

        try applyBuiltinClassifier(&input, ct, null, null, null);

        try std.testing.expectEqualSlices(f64, &expected, input.r);
    }
}

test "resolveBatchThreadCount allows explicit overcommit for IO-bound runs" {
    try std.testing.expectEqual(@as(usize, 10), resolveBatchThreadCount(0, 10));
    try std.testing.expectEqual(@as(usize, 20), resolveBatchThreadCount(20, 10));
}

fn testCountingThread(counter: *std.atomic.Value(usize)) void {
    _ = counter.fetchAdd(1, .seq_cst);
}

test "joinSpawnedThreads joins only initialized thread slots" {
    var counter = std.atomic.Value(usize).init(0);
    // The third slot is never spawned and stays undefined: joining it would crash.
    var threads: [3]std.Thread = undefined;
    threads[0] = try std.Thread.spawn(.{}, testCountingThread, .{&counter});
    threads[1] = try std.Thread.spawn(.{}, testCountingThread, .{&counter});
    joinSpawnedThreads(threads[0..], 2);
    // Both workers had finished when the join returned.
    try std.testing.expectEqual(@as(usize, 2), counter.load(.seq_cst));
}

test "workflow file-first keeps existing output layout" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const jsonl_output_dir = try std.fs.path.join(allocator, &.{ root_path, "jsonl-output" });
    defer allocator.free(jsonl_output_dir);
    const json_output_dir = try std.fs.path.join(allocator, &.{ root_path, "json-output" });
    defer allocator.free(json_output_dir);
    const jsonl_workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow-jsonl.toml" });
    defer allocator.free(jsonl_workflow_path);
    const json_workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow-json.toml" });
    defer allocator.free(json_workflow_path);
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, "tiny.pdb" });
    defer allocator.free(input_path);
    const input_path_2 = try std.fs.path.join(allocator, &.{ input_dir, "tiny2.pdb" });
    defer allocator.free(input_path_2);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path,
        .data =
        \\ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
        \\ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
        \\ATOM      3  N   GLY B   1       5.000   0.000   0.000  1.00 20.00           N
        \\ATOM      4  CA  GLY B   1       6.500   0.000   0.000  1.00 20.00           C
        \\END
        \\
        ,
    });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path_2,
        .data =
        \\ATOM      1  N   ALA A   1       0.000   1.000   0.000  1.00 20.00           N
        \\ATOM      2  CA  ALA A   1       1.500   1.000   0.000  1.00 20.00           C
        \\ATOM      3  N   GLY B   1       5.000   1.000   0.000  1.00 20.00           N
        \\ATOM      4  CA  GLY B   1       6.500   1.000   0.000  1.00 20.00           C
        \\END
        \\
        ,
    });

    const jsonl_workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 1
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "chain_a"
        \\chains = ["A"]
        \\
        \\[[jobs]]
        \\name = "chain_b"
        \\chains = ["B"]
        \\
    , .{ input_dir, jsonl_output_dir });
    defer allocator.free(jsonl_workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = jsonl_workflow_path, .data = jsonl_workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = jsonl_workflow_path });

    const chain_a_jsonl = try std.fs.path.join(allocator, &.{ jsonl_output_dir, "chain_a.jsonl" });
    defer allocator.free(chain_a_jsonl);
    const chain_b_jsonl = try std.fs.path.join(allocator, &.{ jsonl_output_dir, "chain_b.jsonl" });
    defer allocator.free(chain_b_jsonl);
    const chain_a_jsonl_content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, chain_a_jsonl, allocator, .limited(4096));
    defer allocator.free(chain_a_jsonl_content);
    const chain_b_jsonl_content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, chain_b_jsonl, allocator, .limited(4096));
    defer allocator.free(chain_b_jsonl_content);
    try std.testing.expect(std.mem.indexOf(u8, chain_a_jsonl_content, "\"filename\":\"tiny.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, chain_a_jsonl_content, "\"filename\":\"tiny2.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, chain_b_jsonl_content, "\"filename\":\"tiny.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, chain_b_jsonl_content, "\"filename\":\"tiny2.pdb\"") != null);

    const json_workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "json"
        \\
        \\[calculation]
        \\n_points = 1
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "chain_a"
        \\chains = ["A"]
        \\
    , .{ input_dir, json_output_dir });
    defer allocator.free(json_workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = json_workflow_path, .data = json_workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = json_workflow_path });

    const chain_a_json = try std.fs.path.join(allocator, &.{ json_output_dir, "chain_a", "tiny.json" });
    defer allocator.free(chain_a_json);
    const chain_a_json_2 = try std.fs.path.join(allocator, &.{ json_output_dir, "chain_a", "tiny2.json" });
    defer allocator.free(chain_a_json_2);
    const chain_a_json_content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, chain_a_json, allocator, .limited(4096));
    defer allocator.free(chain_a_json_content);
    const chain_a_json_content_2 = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, chain_a_json_2, allocator, .limited(4096));
    defer allocator.free(chain_a_json_content_2);
    try std.testing.expect(std.mem.indexOf(u8, chain_a_json_content, "\"total_area\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, chain_a_json_content_2, "\"total_area\"") != null);
}

/// The runner a batch workflow is routed to (see `runWorkflow`).
const TestWorkflowRunner = enum {
    /// Every file is parsed once and reused by all jobs.
    file_first,
    /// One `runBatch` per job; chosen here by a job whose `auth_chain`
    /// differs from the shared setting.
    job_first,

    fn jobOption(self: TestWorkflowRunner) []const u8 {
        return switch (self) {
            .file_first => "",
            .job_first => "auth_chain = true\n",
        };
    }
};

const test_workflow_runners = [_]TestWorkflowRunner{ .file_first, .job_first };

/// Two chains, one residue each.
const test_two_chain_pdb =
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n" ++
    "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C\n" ++
    "ATOM      3  N   GLY B   1       5.000   0.000   0.000  1.00 20.00           N\n" ++
    "ATOM      4  CA  GLY B   1       6.500   0.000   0.000  1.00 20.00           C\n" ++
    "END\n";

/// A sandbox whose input directory holds two readable and two unreadable
/// structures, and the workflow files that run two jobs over it.
const FailureSandbox = struct {
    sandbox: NamingSandbox,

    const good_files = [_][]const u8{ "good1.pdb", "good2.pdb" };
    /// Unreadable inputs with the reason every runner reports for them.
    const bad_files = [_][2][]const u8{
        .{ "bad1.pdb", "read/parse failed: NoAtomsFound" },
        .{ "bad2.cif", "read/parse failed: NoAtomSiteLoop" },
    };
    const jobs = [_][]const u8{ "chain_a", "everything" };

    fn init() !FailureSandbox {
        var sandbox = try NamingSandbox.init();
        errdefer sandbox.deinit();
        for (good_files) |name| try sandbox.writeInput(name, test_two_chain_pdb);
        try sandbox.writeInput("bad1.pdb", "not a structure\n");
        try sandbox.writeInput("bad2.cif", "data_bad\nloop_\n_cell.length_a\n1.0\n");
        return .{ .sandbox = sandbox };
    }

    fn deinit(self: *FailureSandbox) void {
        self.sandbox.deinit();
    }

    /// Write "workflow.toml" for `runner` and return its path, allocated
    /// from `arena`. `input_name` and `output_name` are below the sandbox
    /// root; without `output_name` the workflow has no output directory and a
    /// single job. `extra` is inserted before the jobs.
    fn writeWorkflow(
        self: FailureSandbox,
        arena: Allocator,
        runner: TestWorkflowRunner,
        input_name: []const u8,
        output_name: ?[]const u8,
        format: []const u8,
        extra: []const u8,
    ) ![]const u8 {
        const output_dir_line = if (output_name) |name|
            try std.fmt.allocPrint(arena, "dir = \"{s}/{s}\"\n", .{ self.sandbox.root, name })
        else
            "";
        const first_job = if (output_name != null) "[[jobs]]\nname = \"chain_a\"\nchains = [\"A\"]\n\n" else "";
        const workflow = try std.fmt.allocPrint(arena,
            \\version = 1
            \\kind = "workflow"
            \\
            \\[input]
            \\dir = "{s}/{s}"
            \\
            \\[output]
            \\{s}format = "{s}"
            \\
            \\[calculation]
            \\n_points = 8
            \\quiet = true
            \\
            \\[classifier]
            \\type = "naccess"
            \\
            \\{s}
            \\{s}[[jobs]]
            \\name = "everything"
            \\{s}
        , .{ self.sandbox.root, input_name, output_dir_line, format, extra, first_job, runner.jobOption() });
        const workflow_path = try std.fs.path.join(arena, &.{ self.sandbox.root, "workflow.toml" });
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });
        return workflow_path;
    }

    fn args(workflow_path: []const u8, n_threads: usize) BatchArgs {
        return .{ .workflow_path = workflow_path, .n_threads = n_threads, .threads_explicit = true };
    }

    /// The rows of a JSONL file below the sandbox root, sorted, as one string.
    fn sortedRows(self: FailureSandbox, arena: Allocator, sub_path: []const u8) ![]const u8 {
        const path = try std.fs.path.join(arena, &.{ self.sandbox.root, sub_path });
        const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, path, arena, .limited(1 << 20));
        var rows = std.ArrayListUnmanaged([]const u8).empty;
        var lines = std.mem.tokenizeScalar(u8, content, '\n');
        while (lines.next()) |line| try rows.append(arena, line);
        std.mem.sort([]const u8, rows.items, {}, NamingSandbox.stringLessThan);
        return std.mem.join(arena, "\n", rows.items);
    }
};

test "workflow writes JSONL error rows for unreadable inputs in every runner" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var failures = try FailureSandbox.init();
    defer failures.deinit();
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    // Rows of each job as the first runner wrote them; the others must agree.
    var reference: [FailureSandbox.jobs.len]?[]const u8 = @splat(null);

    for (test_workflow_runners) |runner| {
        inline for (test_naming_threads) |n_threads| {
            const output_name = try std.fmt.allocPrint(arena, "out-{s}-{d}", .{ @tagName(runner), n_threads });
            const workflow_path = try failures.writeWorkflow(arena, runner, "input", output_name, "jsonl", "");
            try runWorkflow(std.testing.allocator, std.testing.io, FailureSandbox.args(workflow_path, n_threads));

            for (FailureSandbox.jobs, 0..) |job, job_index| {
                const rows = try failures.sortedRows(arena, try std.fmt.allocPrint(arena, "{s}/{s}.jsonl", .{ output_name, job }));

                // One row per input: the readable ones succeed, and each
                // unreadable one has an error row with the usual wording.
                try std.testing.expectEqual(@as(usize, 4), std.mem.count(u8, rows, "\"status\":"));
                for (FailureSandbox.good_files) |name| {
                    const row_start = try std.fmt.allocPrint(arena, "{{\"status\":\"ok\",\"filename\":\"{s}\",", .{name});
                    try std.testing.expect(std.mem.indexOf(u8, rows, row_start) != null);
                }
                for (FailureSandbox.bad_files) |bad| {
                    const row = try std.fmt.allocPrint(arena, "{{\"status\":\"err\",\"filename\":\"{s}\",\"error\":\"{s}\"}}", .{ bad[0], bad[1] });
                    try std.testing.expect(std.mem.indexOf(u8, rows, row) != null);
                }

                if (reference[job_index]) |expected| {
                    try std.testing.expectEqualStrings(expected, rows);
                } else {
                    reference[job_index] = rows;
                }
            }
        }
    }
}

/// `count` failures named "bad00.pdb", "bad01.pdb", ... Caller frees with
/// `freeTestFailures`.
fn makeTestFailures(allocator: Allocator, count: usize) ![]Failure {
    const failures = try allocator.alloc(Failure, count);
    var made: usize = 0;
    errdefer freeTestFailures(allocator, failures[0..made]);
    for (failures, 0..) |*failure, i| {
        failure.* = .{
            .name = try std.fmt.allocPrint(allocator, "bad{d:0>2}.pdb", .{i}),
            .reason = "read/parse failed: NoAtomsFound",
        };
        made += 1;
    }
    return failures;
}

fn freeTestFailures(allocator: Allocator, failures: []const Failure) void {
    for (failures) |failure| allocator.free(failure.name);
    allocator.free(failures);
}

fn expectFailureReport(expected: []const u8, report: FailureReport) !void {
    const text = try report.format(std.testing.allocator);
    defer std.testing.allocator.free(text);
    try std.testing.expectEqualStrings(expected, text);
}

test "failure report lists the failed inputs with their reasons" {
    const failures = [_]Failure{
        .{ .name = "bad.pdb", .reason = "read/parse failed: NoAtomsFound" },
        .{ .name = "empty.cif", .reason = "read/parse failed: NoAtomSiteLoop" },
    };

    // Nothing failed: no report at all, so quiet mode stays quiet.
    try expectFailureReport("", .{ .total = 4, .failed = 0, .failures = &.{} });

    try expectFailureReport(
        \\2 of 4 inputs failed:
        \\  bad.pdb: read/parse failed: NoAtomsFound
        \\  empty.cif: read/parse failed: NoAtomSiteLoop
        \\
    , .{ .total = 4, .failed = 2, .failures = &failures });

    // The JSONL destination only matters when some failures are not listed.
    try expectFailureReport(
        \\Job 'chain_a': 2 of 2 inputs failed:
        \\  bad.pdb: read/parse failed: NoAtomsFound
        \\  empty.cif: read/parse failed: NoAtomSiteLoop
        \\
    , .{ .job = "chain_a", .total = 2, .failed = 2, .failures = &failures, .jsonl = .{ .file = "out/chain_a.jsonl" } });

    try expectFailureReport(
        \\1 of 1 input failed:
        \\  bad.pdb: read/parse failed: NoAtomsFound
        \\
    , .{ .total = 1, .failed = 1, .failures = failures[0..1] });

    try expectFailureReport(
        \\1 of 3 selections failed:
        \\  bad.pdb: read/parse failed: NoAtomsFound
        \\
    , .{ .unit = .selection, .total = 3, .failed = 1, .failures = failures[0..1] });
    try expectFailureReport(
        \\1 of 1 interface failed:
        \\  bad.pdb: read/parse failed: NoAtomsFound
        \\
    , .{ .unit = .interface, .total = 1, .failed = 1, .failures = failures[0..1] });
}

test "failure report lists at most max_listed_failures inputs and counts the rest" {
    const allocator = std.testing.allocator;
    const failures = try makeTestFailures(allocator, max_listed_failures + 5);
    defer freeTestFailures(allocator, failures);

    var listed = std.Io.Writer.Allocating.init(allocator);
    defer listed.deinit();
    for (failures[0..max_listed_failures]) |failure| {
        try listed.writer.print("  {s}: {s}\n", .{ failure.name, failure.reason });
    }

    const Case = struct { jsonl: JsonlDestination, last_line: []const u8 };
    const cases = [_]Case{
        .{ .jsonl = .none, .last_line = "  ... and 5 more\n" },
        .{ .jsonl = .stdout, .last_line = "  ... and 5 more (every failure is a \"status\":\"err\" row in the JSONL output)\n" },
        .{ .jsonl = .{ .file = "out/results.jsonl" }, .last_line = "  ... and 5 more (every failure is a \"status\":\"err\" row in out/results.jsonl)\n" },
    };
    for (cases) |case| {
        const expected = try std.mem.concat(allocator, u8, &.{ "25 of 300 inputs failed:\n", listed.written(), case.last_line });
        defer allocator.free(expected);
        try expectFailureReport(expected, .{ .total = 300, .failed = 25, .failures = failures, .jsonl = case.jsonl });
    }

    // Exactly the limit: every failure is listed and nothing is left to count.
    const at_limit = try std.mem.concat(allocator, u8, &.{ "20 of 300 inputs failed:\n", listed.written() });
    defer allocator.free(at_limit);
    try expectFailureReport(at_limit, .{ .total = 300, .failed = max_listed_failures, .failures = failures[0..max_listed_failures] });

    // Failures that could not be recorded are counted, not lost.
    try expectFailureReport(
        \\7 of 300 inputs failed:
        \\  bad00.pdb: read/parse failed: NoAtomsFound
        \\  bad01.pdb: read/parse failed: NoAtomsFound
        \\  ... and 5 more
        \\
    , .{ .total = 300, .failed = 7, .failures = failures[0..2] });

    // A copy for later keeps what is listed and owns its strings.
    var arena_state = std.heap.ArenaAllocator.init(allocator);
    defer arena_state.deinit();
    const original = FailureReport{ .job = "all", .total = 300, .failed = 25, .failures = failures, .jsonl = .{ .file = "out/all.jsonl" } };
    const copy = try original.dupe(arena_state.allocator());
    try std.testing.expectEqual(@as(usize, max_listed_failures), copy.failures.len);
    try std.testing.expect(copy.failures[0].name.ptr != failures[0].name.ptr);
    const original_text = try original.format(allocator);
    defer allocator.free(original_text);
    try expectFailureReport(original_text, copy);
}

fn recordTestFailures(log: *FailureLog, first: usize) void {
    var name_buf: [32]u8 = undefined;
    for (0..50) |i| {
        const name = std.fmt.bufPrint(&name_buf, "file{d:0>3}.pdb", .{first + i}) catch unreachable;
        log.record(std.testing.io, name, null, "read/parse failed: NoAtomsFound");
    }
}

test "FailureLog orders failures by name whatever order they were recorded in" {
    var log = FailureLog{ .allocator = std.testing.allocator };
    defer log.deinit();

    log.record(std.testing.io, "b.cif", "b.cif", "read/parse failed: NoAtomSiteLoop");
    log.record(std.testing.io, "a.cif", "pair-2", "partner B chain not found");
    log.record(std.testing.io, "a.cif", "pair-1", "partner A chain not found");
    log.record(std.testing.io, "c.pdb", null, "read/parse failed: NoAtomsFound");

    try expectFailureReport(
        \\4 of 9 interfaces failed:
        \\  a.cif [pair-1]: partner A chain not found
        \\  a.cif [pair-2]: partner B chain not found
        \\  b.cif: read/parse failed: NoAtomSiteLoop
        \\  c.pdb: read/parse failed: NoAtomsFound
        \\
    , .{ .unit = .interface, .total = 9, .failed = 4, .failures = log.sorted() });

    // Workers record concurrently.
    var shared = FailureLog{ .allocator = std.testing.allocator };
    defer shared.deinit();
    var threads: [4]std.Thread = undefined;
    for (&threads, 0..) |*thread, t| {
        thread.* = try std.Thread.spawn(.{}, recordTestFailures, .{ &shared, (threads.len - 1 - t) * 50 });
    }
    for (threads) |thread| thread.join();
    const sorted = shared.sorted();
    try std.testing.expectEqual(@as(usize, 200), sorted.len);
    try std.testing.expectEqualStrings("file000.pdb", sorted[0].name);
    try std.testing.expectEqualStrings("file199.pdb", sorted[199].name);
}

test "batch result reports its failed inputs in input order in both runners" {
    const allocator = std.testing.allocator;
    var failures = try FailureSandbox.init();
    defer failures.deinit();

    inline for (test_naming_threads) |n_threads| {
        // Per-file output: the failed inputs leave no output file behind, so
        // the report is the only trace of them.
        var result = try failures.sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 2), result.failed);

        var buffer: [max_listed_failures]Failure = undefined;
        try expectFailureReport(
            \\2 of 4 inputs failed:
            \\  bad1.pdb: read/parse failed: NoAtomsFound
            \\  bad2.cif: read/parse failed: NoAtomSiteLoop
            \\
        , result.failureReport(&buffer, null, .none));
        try expectFailureReport(
            \\Job 'everything': 2 of 4 inputs failed:
            \\  bad1.pdb: read/parse failed: NoAtomsFound
            \\  bad2.cif: read/parse failed: NoAtomSiteLoop
            \\
        , result.failureReport(&buffer, "everything", .stdout));
    }
    try failures.sandbox.expectTree(&.{
        "input/",
        "input/bad1.pdb",
        "input/bad2.cif",
        "input/good1.pdb",
        "input/good2.pdb",
        "out1/",
        "out1/good1.json",
        "out1/good2.json",
        "out4/",
        "out4/good1.json",
        "out4/good2.json",
    });

    // More failures than the report lists
    var name_buf: [32]u8 = undefined;
    for (0..max_listed_failures + 3) |i| {
        const name = try std.fmt.bufPrint(&name_buf, "worse{d:0>2}.pdb", .{i});
        try failures.sandbox.writeInput(name, "not a structure\n");
    }
    const jsonl_path = try failures.sandbox.path("results.jsonl");
    defer allocator.free(jsonl_path);
    inline for (test_naming_threads) |n_threads| {
        var config = NamingSandbox.config(n_threads);
        config.output_format = .jsonl;
        config.store_atom_areas = true;
        var result = try runBatch(allocator, std.testing.io, failures.sandbox.input_dir, null, config, jsonl_path);
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 27), result.total_files);
        try std.testing.expectEqual(@as(usize, 25), result.failed);

        var buffer: [max_listed_failures]Failure = undefined;
        const report = result.failureReport(&buffer, null, JsonlDestination.of(config, jsonl_path));
        try std.testing.expectEqual(@as(usize, max_listed_failures), report.failures.len);
        try std.testing.expectEqualStrings("bad1.pdb", report.failures[0].name);
        try std.testing.expectEqualStrings("worse17.pdb", report.failures[max_listed_failures - 1].name);

        const text = try report.format(allocator);
        defer allocator.free(text);
        try std.testing.expect(std.mem.startsWith(u8, text, "25 of 27 inputs failed:\n  bad1.pdb: read/parse failed: NoAtomsFound\n"));
        const last_line = try std.fmt.allocPrint(allocator, "\n  ... and 5 more (every failure is a \"status\":\"err\" row in {s})\n", .{jsonl_path});
        defer allocator.free(last_line);
        try std.testing.expect(std.mem.endsWith(u8, text, last_line));
        try std.testing.expectEqual(@as(usize, 2 + max_listed_failures), std.mem.count(u8, text, "\n"));

        // The rows the last line points to
        const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(1 << 16));
        defer allocator.free(content);
        try std.testing.expectEqual(@as(usize, 25), std.mem.count(u8, content, "\"status\":\"err\""));
    }
}

test "JsonlDestination follows the output format" {
    const jsonl = BatchConfig{ .output_format = .jsonl, .jsonl_include_atom_areas = false };
    try std.testing.expectEqualStrings("out.jsonl", JsonlDestination.of(jsonl, "out.jsonl").file);
    try std.testing.expect(JsonlDestination.of(jsonl, null) == .stdout);
    try std.testing.expect(JsonlDestination.of(.{}, null) == .none);
    try std.testing.expect(JsonlDestination.of(.{ .output_format = .csv }, null) == .none);
}

test "workflow in which a whole job fails is an error in every runner" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var failures = try FailureSandbox.init();
    defer failures.deinit();
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    const cwd = std.Io.Dir.cwd();

    for (test_workflow_runners) |runner| {
        inline for (test_naming_threads) |n_threads| {
            const tag = try std.fmt.allocPrint(arena, "{s}-{d}", .{ @tagName(runner), n_threads });

            // The input directory does not exist.
            {
                const output_name = try std.fmt.allocPrint(arena, "missing-{s}", .{tag});
                const workflow_path = try failures.writeWorkflow(arena, runner, "no-such-input", output_name, "jsonl", "");
                const result = runWorkflow(std.testing.allocator, std.testing.io, FailureSandbox.args(workflow_path, n_threads));
                try std.testing.expectError(switch (runner) {
                    .file_first => error.FileNotFound,
                    .job_first => error.WorkflowJobFailed,
                }, result);
            }

            // The JSONL file of the second job cannot be created: a directory
            // is in its place.
            {
                const output_name = try std.fmt.allocPrint(arena, "blocked-{s}", .{tag});
                try cwd.createDirPath(std.testing.io, try std.fs.path.join(arena, &.{ failures.sandbox.root, output_name, "everything.jsonl" }));
                const workflow_path = try failures.writeWorkflow(arena, runner, "input", output_name, "jsonl", "");
                const result = runWorkflow(std.testing.allocator, std.testing.io, FailureSandbox.args(workflow_path, n_threads));
                try std.testing.expectError(switch (runner) {
                    .file_first => error.IsDir,
                    .job_first => error.WorkflowJobFailed,
                }, result);

                // The job-first runner still runs the other job; its failed
                // inputs are rows of that job, not failed jobs.
                if (runner == .job_first) {
                    const rows = try failures.sortedRows(arena, try std.fmt.allocPrint(arena, "{s}/chain_a.jsonl", .{output_name}));
                    try std.testing.expectEqual(@as(usize, 2), std.mem.count(u8, rows, "\"status\":\"ok\""));
                    try std.testing.expectEqual(@as(usize, 2), std.mem.count(u8, rows, "\"status\":\"err\""));
                }
            }

            // Two inputs share a per-file output name.
            {
                const input_name = try std.fmt.allocPrint(arena, "clash-{s}", .{tag});
                const output_name = try std.fmt.allocPrint(arena, "clash-out-{s}", .{tag});
                const input_dir = try std.fs.path.join(arena, &.{ failures.sandbox.root, input_name });
                try cwd.createDirPath(std.testing.io, input_dir);
                for ([_][]const u8{ "same.pdb", "same.ent" }) |name| {
                    try cwd.writeFile(std.testing.io, .{ .sub_path = try std.fs.path.join(arena, &.{ input_dir, name }), .data = test_two_chain_pdb });
                }
                const workflow_path = try failures.writeWorkflow(arena, runner, input_name, output_name, "json", "");
                const result = runWorkflow(std.testing.allocator, std.testing.io, FailureSandbox.args(workflow_path, n_threads));
                try std.testing.expectError(switch (runner) {
                    .file_first => error.OutputNameCollision,
                    .job_first => error.WorkflowJobFailed,
                }, result);
                // Nothing was written for the rejected jobs.
                try std.testing.expectError(
                    error.FileNotFound,
                    cwd.access(std.testing.io, try std.fs.path.join(arena, &.{ failures.sandbox.root, output_name }), .{}),
                );
            }
        }
    }
}

test "formatFailedWorkflowJobs counts and names the failed jobs" {
    const allocator = std.testing.allocator;

    const one = try formatFailedWorkflowJobs(allocator, &.{"chain_a"}, 3);
    defer allocator.free(one);
    try std.testing.expectEqualStrings("1 of 3 jobs failed: chain_a\n", one);

    const all = try formatFailedWorkflowJobs(allocator, &.{ "chain_a", "everything" }, 2);
    defer allocator.free(all);
    try std.testing.expectEqualStrings("2 of 2 jobs failed: chain_a, everything\n", all);

    const single = try formatFailedWorkflowJobs(allocator, &.{"only"}, 1);
    defer allocator.free(single);
    try std.testing.expectEqualStrings("1 of 1 job failed: only\n", single);
}

test "workflow rejects colliding per-file output names before creating output" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root_path, "output" });
    defer allocator.free(output_dir);
    const workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    for ([_][]const u8{ "tiny.pdb", "tiny.ent" }) |filename| {
        const path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{
            .sub_path = path,
            .data =
            \\ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
            \\ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
            \\END
            \\
            ,
        });
    }

    inline for (.{ "json", "jsonl" }) |format| {
        const workflow = try std.fmt.allocPrint(allocator,
            \\version = 1
            \\kind = "workflow"
            \\
            \\[input]
            \\dir = "{s}"
            \\
            \\[output]
            \\dir = "{s}"
            \\format = "{s}"
            \\
            \\[calculation]
            \\n_points = 1
            \\quiet = true
            \\
            \\[classifier]
            \\type = "naccess"
            \\
            \\[[jobs]]
            \\name = "chain_a"
            \\chains = ["A"]
            \\
        , .{ input_dir, output_dir, format });
        defer allocator.free(workflow);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

        if (comptime std.mem.eql(u8, format, "json")) {
            try std.testing.expectError(
                error.OutputNameCollision,
                runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path }),
            );
            try std.testing.expectError(
                error.FileNotFound,
                std.Io.Dir.cwd().access(std.testing.io, output_dir, .{}),
            );
        } else {
            // JSONL keeps one record per input, so the same inputs are accepted.
            try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });
            const jsonl_path = try std.fs.path.join(allocator, &.{ output_dir, "chain_a.jsonl" });
            defer allocator.free(jsonl_path);
            const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(4096));
            defer allocator.free(content);
            try std.testing.expect(std.mem.indexOf(u8, content, "\"filename\":\"tiny.pdb\"") != null);
            try std.testing.expect(std.mem.indexOf(u8, content, "\"filename\":\"tiny.ent\"") != null);
        }
    }
}

test "workflow reports errors raised after job setup without freeing job states twice" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root_path, "output" });
    defer allocator.free(output_dir);
    const workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    for ([_][]const u8{ "tiny.pdb", "tiny2.pdb" }) |filename| {
        const path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{
            .sub_path = path,
            .data =
            \\ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
            \\ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
            \\END
            \\
            ,
        });
    }

    // The bitmask LUT is built after every job state has been set up, and
    // rejects point counts above its supported range.
    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 2000
        \\use_bitmask = true
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "chain_a"
        \\chains = ["A"]
        \\
        \\[[jobs]]
        \\name = "all"
        \\
    , .{ input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try std.testing.expectError(
        error.UnsupportedNPoints,
        runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path }),
    );
}

test "workflow chain map selects per-file PDB and mmCIF chain complexes" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root, "output" });
    defer allocator.free(output_dir);
    const pdb_path = try std.fs.path.join(allocator, &.{ input_dir, "pdb-complex.pdb" });
    defer allocator.free(pdb_path);
    const auth_cif_path = try std.fs.path.join(allocator, &.{ input_dir, "auth-complex.cif" });
    defer allocator.free(auth_cif_path);
    const label_cif_path = try std.fs.path.join(allocator, &.{ input_dir, "label-complex.cif" });
    defer allocator.free(label_cif_path);
    const map_path = try std.fs.path.join(allocator, &.{ root, "chains.csv" });
    defer allocator.free(map_path);
    const workflow_path = try std.fs.path.join(allocator, &.{ root, "workflow.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = pdb_path,
        .data =
        \\ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
        \\ATOM      2  N   ALA B   1       4.000   0.000   0.000  1.00 20.00           N
        \\ATOM      3  N   ALA C   1       8.000   0.000   0.000  1.00 20.00           N
        \\END
        \\
        ,
    });

    const cif_template =
        \\data_CHAIN_MAP
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.auth_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N ALA L1 A 1 0.000 0.000 0.000
        \\ATOM 2 N N ALA L2 B 1 4.000 0.000 0.000
        \\ATOM 3 N N ALA L3 C 1 8.000 0.000 0.000
        \\#
        \\
    ;
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = auth_cif_path, .data = cif_template });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = label_cif_path, .data = cif_template });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = map_path,
        .data =
        \\filename,chains,asym_id_type
        \\pdb-complex.pdb,"A,C",label
        \\auth-complex.cif,"A,C",auth
        \\label-complex.cif,"L1,L3",label
        \\
        ,
    });

    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\threads = 2
        \\n_points = 8
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "selected_complexes"
        \\chain_map = "{s}"
        \\
    , .{ input_dir, output_dir, map_path });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const output_path = try std.fs.path.join(allocator, &.{ output_dir, "selected_complexes.jsonl" });
    defer allocator.free(output_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(16384));
    defer allocator.free(content);

    var found: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqualStrings("ok", object.get("status").?.string);
        try std.testing.expectEqual(@as(usize, 2), object.get("atom_areas").?.array.items.len);
        found += 1;
    }
    try std.testing.expectEqual(@as(usize, 3), found);
}

test "selection map LPT claims highest-cost files first with deterministic ties" {
    const allocator = std.testing.allocator;
    var map = try chain_map.parseCsv(
        allocator,
        "filename,id,chains\n" ++
            "z-heavy.pdb,z1,A\n" ++
            "z-heavy.pdb,z2,B\n" ++
            "z-heavy.pdb,z3,C\n" ++
            "a-heavy.pdb,a1,A\n" ++
            "a-heavy.pdb,a2,B\n" ++
            "a-heavy.pdb,a3,C\n" ++
            "middle.pdb,m1,A\n" ++
            "middle.pdb,m2,B\n" ++
            "light.pdb,light,A\n",
    );
    defer map.deinit();

    var files = [_][]const u8{
        "z-heavy.pdb",
        "unmapped.pdb",
        "light.pdb",
        "middle.pdb",
        "a-heavy.pdb",
    };
    scheduleSelectionMapFiles(files[0..], &map);

    var next_file = std.atomic.Value(usize).init(0);
    const expected = [_][]const u8{
        "a-heavy.pdb",
        "z-heavy.pdb",
        "middle.pdb",
        "light.pdb",
        "unmapped.pdb",
    };
    for (expected) |filename| {
        try std.testing.expectEqualStrings(filename, claimSelectionMapFile(files[0..], &next_file).?);
    }
    try std.testing.expect(claimSelectionMapFile(files[0..], &next_file) == null);
}

test "selection map scheduling leaves legacy one-row file order unchanged" {
    const allocator = std.testing.allocator;
    var map = try chain_map.parseCsv(
        allocator,
        "filename,chains\n" ++
            "charlie.pdb,A\n" ++
            "alpha.pdb,A\n" ++
            "bravo.pdb,A\n",
    );
    defer map.deinit();

    var files = [_][]const u8{ "alpha.pdb", "bravo.pdb", "charlie.pdb" };
    scheduleSelectionMapFiles(files[0..], &map);

    try std.testing.expectEqualStrings("alpha.pdb", files[0]);
    try std.testing.expectEqualStrings("bravo.pdb", files[1]);
    try std.testing.expectEqualStrings("charlie.pdb", files[2]);
}

test "selection map LPT preserves JSONL results while processing a heavy file first" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];
    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_path = try std.fs.path.join(allocator, &.{ root, "lpt.jsonl" });
    defer allocator.free(output_path);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const pdb =
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      2  N   ALA B   2       4.000   0.000   0.000  1.00 20.00           N  \n" ++
        "END\n";
    for ([_][]const u8{ "a-light.pdb", "z-heavy.pdb" }) |filename| {
        const path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = path, .data = pdb });
    }

    var map = try chain_map.parseCsv(
        allocator,
        "filename,id,chains\n" ++
            "a-light.pdb,light,A\n" ++
            "z-heavy.pdb,heavy-a,A\n" ++
            "z-heavy.pdb,heavy-b,B\n" ++
            "z-heavy.pdb,heavy-ab,\"A,B\"\n",
    );
    defer map.deinit();

    const stats = try runSelectionMapBatch(allocator, std.testing.io, input_dir, .{
        .n_threads = 1,
        .n_points = 1,
        .classifier_type = .naccess,
        .output_format = .jsonl,
        .store_atom_areas = true,
        .quiet = true,
    }, output_path, &map, null);

    try std.testing.expectEqual(@as(usize, 4), stats.successful);
    try std.testing.expectEqual(@as(usize, 0), stats.failed);
    try std.testing.expectEqual(@as(usize, 2), stats.read_parse_count);
    try std.testing.expectEqual(@as(usize, 2), stats.classifier_count);
    try std.testing.expectEqual(@as(usize, 4), stats.calculation_count);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(16384));
    defer allocator.free(content);
    const expected_ids = [_][]const u8{ "heavy-a", "heavy-b", "heavy-ab", "light" };
    var row_index: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqualStrings("ok", object.get("status").?.string);
        try std.testing.expectEqualStrings(expected_ids[row_index], object.get("id").?.string);
        try std.testing.expect(object.get("total_area") != null);
        try std.testing.expect(object.get("atom_areas") != null);
        row_index += 1;
    }
    try std.testing.expectEqual(expected_ids.len, row_index);
}

test "selection map parses and classifies once, reuses chain sets, and emits joinable JSONL" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];
    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, "multi.pdb" });
    defer allocator.free(input_path);
    const output_path = try std.fs.path.join(allocator, &.{ root, "selections.jsonl" });
    defer allocator.free(output_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path,
        .data = "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
            "ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C  \n" ++
            "ATOM      3  N   ALA B   2       3.000   0.000   0.000  1.00 20.00           N  \n" ++
            "ATOM      4  CA  ALA B   2       4.500   0.000   0.000  1.00 20.00           C  \n" ++
            "ATOM      5  N   SER C   3       6.000   0.000   0.000  1.00 20.00           N  \n" ++
            "ATOM      6  CA  SER C   3       7.500   0.000   0.000  1.00 20.00           C  \n" ++
            "END\n",
    });

    var map = try chain_map.parseCsv(
        allocator,
        "filename,id,chains,asym_id_type\n" ++
            "multi.pdb,a,A,label\n" ++
            "multi.pdb,bc,\"B,C\",label\n" ++
            "multi.pdb,abc,\"A,B,C\",label\n" ++
            "multi.pdb,bc-copy,\"C,B\",label\n" ++
            "multi.pdb,bad,Z,label\n" ++
            "absent.pdb,absent,A,label\n",
    );
    defer map.deinit();

    var failures = FailureLog{ .allocator = allocator };
    defer failures.deinit();
    const stats = try runSelectionMapBatch(allocator, std.testing.io, input_dir, .{
        .n_threads = 4,
        .n_points = 128,
        .use_bitmask = true,
        .precision = .f64,
        .classifier_type = .ccd,
        .include_hetatm = true,
        .output_format = .jsonl,
        .store_atom_areas = true,
        .residue_map = true,
        .jsonl_include_atom_areas = true,
        .jsonl_include_atom_identity = true,
        .quiet = true,
    }, output_path, &map, &failures);

    try std.testing.expectEqual(@as(usize, 4), stats.successful);
    try std.testing.expectEqual(@as(usize, 2), stats.failed);
    // The failed selections, as the workflow reports them for the job
    try expectFailureReport(
        \\Job 'sel': 2 of 6 selections failed:
        \\  absent.pdb [absent]: input structure not found
        \\  multi.pdb [bad]: selected chain not found: Z
        \\
    , .{
        .job = "sel",
        .unit = .selection,
        .total = stats.successful + stats.failed,
        .failed = stats.failed,
        .failures = failures.sorted(),
    });
    try std.testing.expectEqual(@as(usize, 1), stats.read_parse_count);
    try std.testing.expectEqual(@as(usize, 1), stats.classifier_count);
    try std.testing.expectEqual(@as(usize, 3), stats.calculation_count);
    try std.testing.expectEqual(@as(usize, 1), stats.file_threads);
    try std.testing.expectEqual(@as(usize, 4), stats.sasa_threads);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(65536));
    defer allocator.free(content);
    const expected_ids = [_][]const u8{ "a", "bc", "abc", "bc-copy", "bad", "absent" };
    var row_index: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqualStrings(expected_ids[row_index], object.get("id").?.string);
        try std.testing.expect(object.get("filename") != null);
        try std.testing.expect(object.get("chains") != null);
        if (row_index < 4) {
            try std.testing.expectEqualStrings("ok", object.get("status").?.string);
            try std.testing.expect(object.get("total_area") != null);
            const areas = object.get("atom_areas").?.array.items;
            const source_indices = object.get("source_atom_index").?.array.items;
            try std.testing.expectEqual(areas.len, source_indices.len);
            try std.testing.expectEqual(areas.len, object.get("atom_element").?.array.items.len);
            try std.testing.expect(object.get("residue_sasa") != null);
            if (std.mem.eql(u8, expected_ids[row_index], "a")) {
                try std.testing.expectEqual(@as(i64, 0), source_indices[0].integer);
                try std.testing.expectEqual(@as(i64, 1), source_indices[1].integer);
            } else if (std.mem.eql(u8, expected_ids[row_index], "bc") or
                std.mem.eql(u8, expected_ids[row_index], "bc-copy"))
            {
                try std.testing.expectEqual(@as(i64, 2), source_indices[0].integer);
                try std.testing.expectEqual(@as(i64, 5), source_indices[3].integer);
            } else {
                try std.testing.expectEqual(@as(i64, 0), source_indices[0].integer);
                try std.testing.expectEqual(@as(i64, 5), source_indices[5].integer);
            }
        } else {
            try std.testing.expectEqualStrings("err", object.get("status").?.string);
            try std.testing.expect(object.get("error") != null);
        }
        row_index += 1;
    }
    try std.testing.expectEqual(expected_ids.len, row_index);
}

test "selection map parallelizes files and keeps internal SASA single-threaded" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];
    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_path = try std.fs.path.join(allocator, &.{ root, "parallel.jsonl" });
    defer allocator.free(output_path);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const pdb =
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C  \n" ++
        "END\n";
    for ([_][]const u8{ "one.pdb", "two.pdb" }) |filename| {
        const path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = path, .data = pdb });
    }

    var map = try chain_map.parseCsv(
        allocator,
        "filename,id,chains\n" ++
            "one.pdb,one,A\n" ++
            "two.pdb,two,A\n",
    );
    defer map.deinit();
    const stats = try runSelectionMapBatch(allocator, std.testing.io, input_dir, .{
        .n_threads = 4,
        .n_points = 8,
        .classifier_type = .naccess,
        .output_format = .jsonl,
        .store_atom_areas = true,
        .quiet = true,
    }, output_path, &map, null);

    try std.testing.expectEqual(@as(usize, 2), stats.file_threads);
    try std.testing.expectEqual(@as(usize, 1), stats.sasa_threads);
    try std.testing.expectEqual(@as(usize, 2), stats.read_parse_count);
    try std.testing.expectEqual(@as(usize, 2), stats.classifier_count);
    try std.testing.expectEqual(@as(usize, 2), stats.calculation_count);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(8192));
    defer allocator.free(content);
    var row_count: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        try std.testing.expectEqualStrings("ok", parsed.value.object.get("status").?.string);
        row_count += 1;
    }
    try std.testing.expectEqual(@as(usize, 2), row_count);
}

test "workflow BSA analysis writes analysis JSONL" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root, "output" });
    defer allocator.free(output_dir);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const pdb_path = try std.fs.path.join(allocator, &.{ input_dir, "tiny_ab.pdb" });
    defer allocator.free(pdb_path);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = pdb_path, .data = "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C  \n" ++
        "ATOM      3  N   ALA B   2       3.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      4  CA  ALA B   2       4.500   0.000   0.000  1.00 20.00           C  \n" ++
        "END\n" });

    const workflow_path = try std.fs.path.join(allocator, &.{ root, "bsa.toml" });
    defer allocator.free(workflow_path);
    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 16
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "interface_ab"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\level = "residue"
    , .{ input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const jsonl_path = try std.fs.path.join(allocator, &.{ output_dir, "interface_ab.jsonl" });
    defer allocator.free(jsonl_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(8192));
    defer allocator.free(content);

    try std.testing.expect(std.mem.indexOf(u8, content, "\"analysis\":\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"delta_sasa_total\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"residue_delta_sasa\"") != null);
}

test "workflow BSA analysis writes chain IDs longer than four characters in full" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    const cwd = std.Io.Dir.cwd();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root = root_buf[0..try tmp_dir.dir.realPath(std.testing.io, &root_buf)];
    const input_dir = try std.fs.path.join(arena, &.{ root, "input" });
    const output_dir = try std.fs.path.join(arena, &.{ root, "output" });
    try cwd.createDirPath(std.testing.io, input_dir);
    try cwd.writeFile(std.testing.io, .{
        .sub_path = try std.fs.path.join(arena, &.{ input_dir, "long.cif" }),
        .data =
        \\data_LONG
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N  GLY LONGA 1 0.000 0.000 0.000
        \\ATOM 2 C CA GLY LONGA 1 1.500 0.000 0.000
        \\ATOM 3 N N  ALA LONGB 2 3.000 0.000 0.000
        \\ATOM 4 C CA ALA LONGB 2 4.500 0.000 0.000
        \\#
        \\
        ,
    });

    const workflow_path = try std.fs.path.join(arena, &.{ root, "bsa.toml" });
    try cwd.writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = try std.fmt.allocPrint(arena,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 16
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "long"
        \\partner_a = ["LONGA"]
        \\partner_b = ["LONGB"]
        \\level = "residue"
        \\atom_output = true
        \\
    , .{ input_dir, output_dir }) });

    try runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const content = try cwd.readFileAlloc(std.testing.io, try std.fs.path.join(arena, &.{ output_dir, "long.jsonl" }), arena, .limited(1 << 16));
    const row = (try std.json.parseFromSliceLeaky(std.json.Value, arena, std.mem.trimEnd(u8, content, "\n"), .{})).object;
    try std.testing.expectEqualStrings("ok", row.get("status").?.string);

    // One residue and two atoms per chain; both arrays name the same chains.
    const residue_chain = row.get("residue_chain").?.array.items;
    try std.testing.expectEqual(@as(usize, 2), residue_chain.len);
    try std.testing.expectEqualStrings("LONGA", residue_chain[0].string);
    try std.testing.expectEqualStrings("LONGB", residue_chain[1].string);
    const atom_chain = row.get("atom_chain").?.array.items;
    try std.testing.expectEqual(@as(usize, 4), atom_chain.len);
    try std.testing.expectEqualStrings("LONGA", atom_chain[0].string);
    try std.testing.expectEqualStrings("LONGB", atom_chain[3].string);
}

test "workflow BSA analysis uses per-file multi-chain interface map" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root, "output" });
    defer allocator.free(output_dir);
    const cif_path = try std.fs.path.join(allocator, &.{ input_dir, "abcd.cif" });
    defer allocator.free(cif_path);
    const map_path = try std.fs.path.join(allocator, &.{ root, "interfaces.csv" });
    defer allocator.free(map_path);
    const workflow_path = try std.fs.path.join(allocator, &.{ root, "bsa-map.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = cif_path,
        .data =
        \\data_ABCD
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.auth_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N GLY L1 A 1 0.000 0.000 0.000
        \\ATOM 2 N N GLY L2 B 1 2.000 0.000 0.000
        \\ATOM 3 N N GLY L3 C 1 4.000 0.000 0.000
        \\ATOM 4 N N GLY L4 D 1 6.000 0.000 0.000
        \\#
        \\
        ,
    });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = map_path,
        .data =
        \\filename,partner_a,partner_b,asym_id_type
        \\abcd.cif,"A,B","C,D",auth
        \\
        ,
    });

    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 8
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "interfaces"
        \\chain_map = "{s}"
        \\level = "residue"
        \\
    , .{ input_dir, output_dir, map_path });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const output_path = try std.fs.path.join(allocator, &.{ output_dir, "interfaces.jsonl" });
    defer allocator.free(output_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(16384));
    defer allocator.free(content);
    const line = std.mem.trim(u8, content, " \t\r\n");
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
    defer parsed.deinit();
    const object = parsed.value.object;

    try std.testing.expectEqualStrings("ok", object.get("status").?.string);
    try std.testing.expectEqualStrings("abcd.cif", object.get("id").?.string);
    const partner_a = object.get("partner_a").?.array.items;
    const partner_b = object.get("partner_b").?.array.items;
    try std.testing.expectEqual(@as(usize, 2), partner_a.len);
    try std.testing.expectEqualStrings("A", partner_a[0].string);
    try std.testing.expectEqualStrings("B", partner_a[1].string);
    try std.testing.expectEqual(@as(usize, 2), partner_b.len);
    try std.testing.expectEqualStrings("C", partner_b[0].string);
    try std.testing.expectEqualStrings("D", partner_b[1].string);

    const residue_chains = object.get("residue_chain").?.array.items;
    try std.testing.expectEqual(@as(usize, 4), residue_chains.len);
    try std.testing.expectEqualStrings("A", residue_chains[0].string);
    try std.testing.expectEqualStrings("D", residue_chains[3].string);
    try std.testing.expectEqual(@as(usize, 4), object.get("residue_partner").?.array.items.len);
    try std.testing.expectEqual(@as(usize, 4), object.get("residue_sasa_isolated").?.array.items.len);
    try std.testing.expectEqual(@as(usize, 4), object.get("residue_sasa_complex").?.array.items.len);
}

fn testJsonNumber(value: std.json.Value) f64 {
    return switch (value) {
        .float => |number| number,
        .integer => |number| @floatFromInt(number),
        else => unreachable,
    };
}

fn testJsonNumberArraySum(values: std.json.Array) f64 {
    var total: f64 = 0;
    for (values.items) |value| total += testJsonNumber(value);
    return total;
}

test "workflow BSA analysis emits one detailed row per interface and stable error IDs" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root, "output" });
    defer allocator.free(output_dir);
    const pdb_path = try std.fs.path.join(allocator, &.{ input_dir, "multi.pdb" });
    defer allocator.free(pdb_path);
    const map_path = try std.fs.path.join(allocator, &.{ root, "interfaces.csv" });
    defer allocator.free(map_path);
    const workflow_path = try std.fs.path.join(allocator, &.{ root, "bsa-multi.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = pdb_path,
        .data = "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
            "ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C  \n" ++
            "ATOM      3  N   ALA B   2       3.000   0.000   0.000  1.00 20.00           N  \n" ++
            "ATOM      4  CA  ALA B   2       4.500   0.000   0.000  1.00 20.00           C  \n" ++
            "ATOM      5  N   SER C   3       6.000   0.000   0.000  1.00 20.00           N  \n" ++
            "ATOM      6  CA  SER C   3       7.500   0.000   0.000  1.00 20.00           C  \n" ++
            "END\n",
    });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = map_path,
        .data =
        \\filename,id,partner_a,partner_b,asym_id_type
        \\multi.pdb,interface-ab,A,B,label
        \\multi.pdb,interface-cb,C,B,label
        \\multi.pdb,interface-invalid,Z,B,label
        \\missing.pdb,interface-missing,A,B,label
        \\
        ,
    });

    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 16
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "interfaces"
        \\chain_map = "{s}"
        \\level = "residue"
        \\atom_output = true
        \\
    , .{ input_dir, output_dir, map_path });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const output_path = try std.fs.path.join(allocator, &.{ output_dir, "interfaces.jsonl" });
    defer allocator.free(output_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(65536));
    defer allocator.free(content);

    var row_count: usize = 0;
    var ok_count: usize = 0;
    var error_count: usize = 0;
    var saw_invalid = false;
    var saw_missing = false;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        row_count += 1;
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        const status = object.get("status").?.string;
        const id = object.get("id").?.string;
        if (std.mem.eql(u8, status, "err")) {
            error_count += 1;
            if (std.mem.eql(u8, id, "interface-invalid")) saw_invalid = true;
            if (std.mem.eql(u8, id, "interface-missing")) saw_missing = true;
            try std.testing.expect(object.get("error") != null);
            continue;
        }

        ok_count += 1;
        try std.testing.expect(std.mem.eql(u8, id, "interface-ab") or std.mem.eql(u8, id, "interface-cb"));
        const atom_isolated = object.get("atom_sasa_isolated").?.array;
        const atom_complex = object.get("atom_sasa_complex").?.array;
        const atom_delta = object.get("atom_delta_sasa").?.array;
        try std.testing.expectEqual(atom_isolated.items.len, atom_complex.items.len);
        try std.testing.expectEqual(atom_isolated.items.len, atom_delta.items.len);
        try std.testing.expectEqual(atom_isolated.items.len, object.get("atom_partner").?.array.items.len);
        try std.testing.expectApproxEqAbs(
            testJsonNumberArraySum(atom_isolated) - testJsonNumberArraySum(atom_complex),
            testJsonNumberArraySum(atom_delta),
            1e-9,
        );

        const residue_isolated = object.get("residue_sasa_isolated").?.array;
        const residue_complex = object.get("residue_sasa_complex").?.array;
        const residue_delta = object.get("residue_delta_sasa").?.array;
        try std.testing.expectApproxEqAbs(
            testJsonNumberArraySum(residue_isolated) - testJsonNumberArraySum(residue_complex),
            testJsonNumberArraySum(residue_delta),
            1e-9,
        );
        try std.testing.expectApproxEqAbs(
            testJsonNumber(object.get("delta_sasa_total").?),
            testJsonNumberArraySum(atom_delta),
            1e-9,
        );
        try std.testing.expectApproxEqAbs(
            testJsonNumber(object.get("delta_sasa_total").?) / 2.0,
            testJsonNumber(object.get("bsa").?),
            1e-9,
        );
    }

    try std.testing.expectEqual(@as(usize, 4), row_count);
    try std.testing.expectEqual(@as(usize, 2), ok_count);
    try std.testing.expectEqual(@as(usize, 2), error_count);
    try std.testing.expect(saw_invalid);
    try std.testing.expect(saw_missing);
}

test "workflow BSA reports every interface across parallel success and error paths" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root, "output" });
    defer allocator.free(output_dir);
    const map_path = try std.fs.path.join(allocator, &.{ root, "interfaces.csv" });
    defer allocator.free(map_path);
    const workflow_path = try std.fs.path.join(allocator, &.{ root, "bsa-parallel.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    const pdb =
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C  \n" ++
        "ATOM      3  N   ALA B   2       3.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      4  CA  ALA B   2       4.500   0.000   0.000  1.00 20.00           C  \n" ++
        "ATOM      5  N   SER C   3       6.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      6  CA  SER C   3       7.500   0.000   0.000  1.00 20.00           C  \n" ++
        "END\n";
    for ([_][]const u8{ "alpha.pdb", "beta.pdb" }) |filename| {
        const pdb_path = try std.fs.path.join(allocator, &.{ input_dir, filename });
        defer allocator.free(pdb_path);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = pdb_path, .data = pdb });
    }
    const broken_path = try std.fs.path.join(allocator, &.{ input_dir, "broken.pdb" });
    defer allocator.free(broken_path);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = broken_path, .data = "not a structure\n" });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = map_path,
        .data =
        \\filename,id,partner_a,partner_b,asym_id_type
        \\alpha.pdb,alpha-ab,A,B,label
        \\alpha.pdb,alpha-cb,C,B,label
        \\beta.pdb,beta-ab,A,B,label
        \\beta.pdb,beta-invalid,Z,B,label
        \\broken.pdb,broken-ab,A,B,label
        \\broken.pdb,broken-cb,C,B,label
        \\missing.pdb,missing-ab,A,B,label
        \\
        ,
    });

    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 16
        \\quiet = false
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "interfaces"
        \\chain_map = "{s}"
        \\level = "residue"
        \\
    , .{ input_dir, output_dir, map_path });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    // std.Progress.start may run only once per process, and a test binary run
    // directly (not through `zig build test`) has already started it in the
    // Zig test runner. Progress is therefore off here, and the file loop is
    // exercised with a no-op progress node.
    try runWorkflow(allocator, std.testing.io, .{
        .workflow_path = workflow_path,
        .n_threads = 4,
        .threads_explicit = true,
        .quiet = false,
        .quiet_explicit = true,
        .show_progress = false,
    });

    const output_path = try std.fs.path.join(allocator, &.{ output_dir, "interfaces.jsonl" });
    defer allocator.free(output_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(65536));
    defer allocator.free(content);

    const expected_alpha = [_][]const u8{ "alpha-ab", "alpha-cb" };
    const expected_beta = [_][]const u8{ "beta-ab", "beta-invalid" };
    const expected_broken = [_][]const u8{ "broken-ab", "broken-cb" };
    var alpha_index: usize = 0;
    var beta_index: usize = 0;
    var broken_index: usize = 0;
    var row_count: usize = 0;
    var ok_count: usize = 0;
    var error_count: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        row_count += 1;
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        const filename = object.get("filename").?.string;
        const id = object.get("id").?.string;
        const status = object.get("status").?.string;

        if (std.mem.eql(u8, filename, "alpha.pdb")) {
            try std.testing.expectEqualStrings(expected_alpha[alpha_index], id);
            alpha_index += 1;
        } else if (std.mem.eql(u8, filename, "beta.pdb")) {
            try std.testing.expectEqualStrings(expected_beta[beta_index], id);
            beta_index += 1;
        } else if (std.mem.eql(u8, filename, "broken.pdb")) {
            try std.testing.expectEqualStrings(expected_broken[broken_index], id);
            broken_index += 1;
        }

        if (std.mem.eql(u8, status, "ok")) {
            ok_count += 1;
            try std.testing.expect(object.get("sasa_partner_a") != null);
            try std.testing.expect(object.get("sasa_partner_b") != null);
            try std.testing.expect(object.get("sasa_complex") != null);
            try std.testing.expect(object.get("delta_sasa_total") != null);
            try std.testing.expect(object.get("bsa") != null);
        } else {
            try std.testing.expectEqualStrings("err", status);
            error_count += 1;
        }
    }

    try std.testing.expectEqual(@as(usize, 7), row_count);
    try std.testing.expectEqual(@as(usize, 3), ok_count);
    try std.testing.expectEqual(@as(usize, 4), error_count);
    try std.testing.expectEqual(expected_alpha.len, alpha_index);
    try std.testing.expectEqual(expected_beta.len, beta_index);
    try std.testing.expectEqual(expected_broken.len, broken_index);
}

test "workflow mmCIF chain filters preserve long chain IDs" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root_path, "output" });
    defer allocator.free(output_dir);
    const workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow.toml" });
    defer allocator.free(workflow_path);
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, "long-chain.cif" });
    defer allocator.free(input_path);
    const input_path_2 = try std.fs.path.join(allocator, &.{ input_dir, "long-chain-2.cif" });
    defer allocator.free(input_path_2);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path,
        .data =
        \\data_LONG_CHAIN
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N ALA ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg 1 0.000 0.000 0.000
        \\ATOM 2 C CA ALA ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg 1 1.500 0.000 0.000
        \\#
        \\
        ,
    });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path_2,
        .data =
        \\data_LONG_CHAIN_2
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N ALA ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef 1 0.000 1.000 0.000
        \\ATOM 2 C CA ALA ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef 1 1.500 1.000 0.000
        \\#
        \\
        ,
    });

    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 1
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "long_chain"
        \\chains = ["ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg"]
        \\
        \\[[jobs]]
        \\name = "prefix_chain"
        \\chains = ["ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"]
        \\
    , .{ input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const long_chain_jsonl = try std.fs.path.join(allocator, &.{ output_dir, "long_chain.jsonl" });
    defer allocator.free(long_chain_jsonl);
    const prefix_jsonl = try std.fs.path.join(allocator, &.{ output_dir, "prefix_chain.jsonl" });
    defer allocator.free(prefix_jsonl);
    const long_chain_content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, long_chain_jsonl, allocator, .limited(4096));
    defer allocator.free(long_chain_content);
    const prefix_content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, prefix_jsonl, allocator, .limited(4096));
    defer allocator.free(prefix_content);

    try std.testing.expect(std.mem.indexOf(u8, long_chain_content, "\"status\":\"ok\",\"filename\":\"long-chain.cif\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, long_chain_content, "\"status\":\"err\",\"filename\":\"long-chain-2.cif\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, prefix_content, "\"status\":\"err\",\"filename\":\"long-chain.cif\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, prefix_content, "\"status\":\"ok\",\"filename\":\"long-chain-2.cif\"") != null);
}

test "workflowJsonlOutputPath uses job file under output dir" {
    const path = try workflowJsonlOutputPath(std.testing.allocator, "results", "chain_A");
    defer std.testing.allocator.free(path);
    try std.testing.expectEqualStrings("results/chain_A.jsonl", path);
}

test "workflow JSONL output options control fields and metadata sidecar" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root_path, "output" });
    defer allocator.free(output_dir);
    const workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow.toml" });
    defer allocator.free(workflow_path);
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, "tiny.pdb" });
    defer allocator.free(input_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path,
        .data =
        \\ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
        \\ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C
        \\END
        \\
        ,
    });

    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[output.jsonl]
        \\atom_areas = false
        \\total_area = true
        \\decimals = 2
        \\metadata = "sidecar"
        \\
        \\[calculation]
        \\n_points = 8
        \\quiet = true
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[[jobs]]
        \\name = "all"
        \\
    , .{ input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const jsonl_path = try std.fs.path.join(allocator, &.{ output_dir, "all.jsonl" });
    defer allocator.free(jsonl_path);
    const meta_path = try std.fs.path.join(allocator, &.{ output_dir, "all.meta.json" });
    defer allocator.free(meta_path);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(4096));
    defer allocator.free(content);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"status\":\"ok\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"total_area\":") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"atom_areas\"") == null);

    const meta = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, meta_path, allocator, .limited(4096));
    defer allocator.free(meta);
    const parsed_meta = try std.json.parseFromSlice(std.json.Value, allocator, meta, .{});
    defer parsed_meta.deinit();
    const meta_obj = parsed_meta.value.object;
    try std.testing.expectEqualStrings("all", meta_obj.get("job").?.string);
    try std.testing.expectEqual(false, meta_obj.get("jsonl").?.object.get("atom_areas").?.bool);
    try std.testing.expectEqual(@as(i64, 2), meta_obj.get("jsonl").?.object.get("decimals").?.integer);
}

test "workflowPerFileOutputDir uses job directory under output dir" {
    const path = try workflowPerFileOutputDir(std.testing.allocator, "results", "complex_AB");
    defer std.testing.allocator.free(path);
    try std.testing.expectEqualStrings("results/complex_AB", path);
}

test "workflowFilesContainSdf detects SDF in pre-scanned files" {
    const files = [_][]const u8{ "protein.pdb", "ligand.sdf", "other.cif" };
    try std.testing.expect(workflowFilesContainSdf(files[0..]));
}

test "workflowRequiresJobFirstForLongChainFormats allows chain-filtered mmCIF and BinaryCIF" {
    const mmcif_files = [_][]const u8{"long-chain.cif"};
    const bcif_files = [_][]const u8{"long-chain.bcif"};
    const pdb_json_files = [_][]const u8{ "chain.pdb", "atoms.json" };

    var long_chain_jobs = [_]workflow_manifest.Job{
        .{ .name = "long", .chains = &.{"ABCDE"} },
    };
    const long_chain_workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .jobs = long_chain_jobs[0..],
    };
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(mmcif_files[0..], long_chain_workflow));
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(bcif_files[0..], long_chain_workflow));
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(pdb_json_files[0..], long_chain_workflow));

    var prefix_chain_jobs = [_]workflow_manifest.Job{
        .{ .name = "prefix", .chains = &.{"ABCD"} },
    };
    const prefix_chain_workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .jobs = prefix_chain_jobs[0..],
    };
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(mmcif_files[0..], prefix_chain_workflow));
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(bcif_files[0..], prefix_chain_workflow));

    var unfiltered_jobs = [_]workflow_manifest.Job{
        .{ .name = "all" },
    };
    const unfiltered_workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .jobs = unfiltered_jobs[0..],
    };
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(mmcif_files[0..], unfiltered_workflow));
    try std.testing.expect(!workflowRequiresJobFirstForLongChainFormats(bcif_files[0..], unfiltered_workflow));
}

test "workflowRequiresJobFirstForAuthChain detects job-level auth-chain mismatch" {
    var jobs = [_]workflow_manifest.Job{
        .{ .name = "label", .chains = &.{"A"}, .auth_chain = false },
        .{ .name = "auth", .chains = &.{"A"}, .auth_chain = true },
    };
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .calculation = .{ .auth_chain = false },
        .jobs = jobs[0..],
    };

    try std.testing.expect(workflowRequiresJobFirstForAuthChain(BatchArgs{}, workflow, BatchConfig{ .use_auth_chain = false }));
    try std.testing.expect(!workflowRequiresJobFirstForAuthChain(BatchArgs{ .use_auth_chain = true }, workflow, BatchConfig{ .use_auth_chain = true }));
}

test "effectiveWorkflowSasaThreads honors explicit and auto thread counts" {
    try std.testing.expectEqual(@as(usize, 7), effectiveWorkflowSasaThreads(BatchConfig{ .n_threads = 7 }));
    try std.testing.expect(effectiveWorkflowSasaThreads(BatchConfig{ .n_threads = 0 }) >= 1);
}

test "BSA workflow uses explicit file workers without nested SASA threads" {
    const config = BatchConfig{ .n_threads = 40 };
    const file_threads = effectiveWorkflowFileThreads(config, 100);
    try std.testing.expectEqual(@as(usize, 40), file_threads);
    try std.testing.expectEqual(@as(usize, 1), effectiveBsaSasaThreads(config, file_threads));
    try std.testing.expectEqual(@as(usize, 40), effectiveBsaSasaThreads(config, 1));
}

test "appendJsonlResultToFile appends without truncating existing JSONL content" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const output_path = try std.fs.path.join(allocator, &.{ root_buf[0..root_len], "results.jsonl" });
    defer allocator.free(output_path);

    try truncateJsonlOutput(std.testing.io, output_path);

    const areas = [_]f64{1.0};
    var first = FileResult{
        .filename = "first.pdb",
        .n_atoms = 1,
        .sasa_time_ns = 1,
        .total_sasa = 1.0,
        .status = .ok,
        .atom_areas = areas[0..],
    };
    var second = FileResult{
        .filename = "second.pdb",
        .n_atoms = 1,
        .sasa_time_ns = 1,
        .total_sasa = 1.0,
        .status = .ok,
        .atom_areas = areas[0..],
    };

    {
        const file = try openJsonlForAppend(std.testing.io, output_path);
        defer file.close(std.testing.io);
        try appendJsonlResultToFile(std.testing.io, file, allocator, &first, .{});
    }
    {
        const file = try openJsonlForAppend(std.testing.io, output_path);
        defer file.close(std.testing.io);
        try appendJsonlResultToFile(std.testing.io, file, allocator, &second, .{});
    }

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(4096));
    defer allocator.free(content);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"filename\":\"first.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"filename\":\"second.pdb\"") != null);
    try std.testing.expectEqual(@as(usize, 2), std.mem.count(u8, content, "\n"));
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

test "FileResult JSONL serializes error result without atom areas" {
    const allocator = std.testing.allocator;
    var result = FileResult{
        .filename = "bad.pdb",
        .n_atoms = 0,
        .sasa_time_ns = 0,
        .total_sasa = 0,
        .status = .err,
        .error_msg = "read/parse failed: InvalidFormat",
    };

    const line = try fileResultToJsonlLine(allocator, &result);
    defer allocator.free(line);

    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
    defer parsed.deinit();

    const obj = parsed.value.object;
    try std.testing.expectEqualStrings("err", obj.get("status").?.string);
    try std.testing.expectEqualStrings("bad.pdb", obj.get("filename").?.string);
    try std.testing.expectEqualStrings("read/parse failed: InvalidFormat", obj.get("error").?.string);
    try std.testing.expect(obj.get("atom_areas") == null);
}

test "JsonlStreamWriter writes many parseable JSONL lines" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const output_path = try std.fs.path.join(allocator, &.{ root_buf[0..root_len], "stream.jsonl" });
    defer allocator.free(output_path);

    const file = try std.Io.Dir.cwd().createFile(std.testing.io, output_path, .{});
    var stream_buf: [64 * 1024]u8 = undefined;
    var stream = JsonlStreamWriter.init(file, std.testing.io, .{}, &stream_buf);

    var arena = std.heap.ArenaAllocator.init(std.heap.page_allocator);
    defer arena.deinit();

    const atom_areas = [_]f64{ 1.0, 2.0, 3.0 };
    var name_buf: [32]u8 = undefined;
    for (0..50) |i| {
        _ = arena.reset(.retain_capacity);
        const name = try std.fmt.bufPrint(&name_buf, "file-{d}.pdb", .{i});
        var result = FileResult{
            .filename = name,
            .n_atoms = atom_areas.len,
            .sasa_time_ns = 1,
            .total_sasa = 6.0,
            .status = .ok,
            .atom_areas = atom_areas[0..],
        };
        stream.writeResult(arena.allocator(), &result);
        try std.testing.expect(!stream.hasError());
    }
    try stream.flush();
    file.close(std.testing.io);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(64 * 1024));
    defer allocator.free(content);
    try std.testing.expectEqual(@as(usize, 50), std.mem.count(u8, content, "\n"));

    var lines = std.mem.tokenizeScalar(u8, content, '\n');
    var count: usize = 0;
    while (lines.next()) |line| {
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqualStrings("ok", object.get("status").?.string);
        try std.testing.expect(object.get("filename") != null);
        try std.testing.expect(object.get("atom_areas") != null);
        count += 1;
    }
    try std.testing.expectEqual(@as(usize, 50), count);
}

test "runBatchParallel writes parseable JSONL with multiple threads" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_path = try std.fs.path.join(allocator, &.{ root_path, "results.jsonl" });
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

    var lines = std.mem.tokenizeScalar(u8, content, '\n');
    var count: usize = 0;
    while (lines.next()) |line| {
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        const object = parsed.value.object;
        try std.testing.expectEqualStrings("ok", object.get("status").?.string);
        try std.testing.expect(object.get("filename") != null);
        try std.testing.expect(object.get("total_area") != null);
        try std.testing.expect(object.get("atom_areas") != null);
        try std.testing.expectEqual(@as(usize, 3), object.get("atom_areas").?.array.items.len);
        count += 1;
    }
    try std.testing.expectEqual(@as(usize, 10), count);
}

test "batch and workflow exclude HETATM by default, also with the CCD classifier" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const input_path = try std.fs.path.join(allocator, &.{ input_dir, "hetatm.pdb" });
    defer allocator.free(input_path);
    const workflow_path = try std.fs.path.join(allocator, &.{ root_path, "workflow.toml" });
    defer allocator.free(workflow_path);

    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{
        .sub_path = input_path,
        .data =
        \\ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N
        \\HETATM    2  O   HOH A   2      20.000   0.000   0.000  1.00 20.00           O
        \\END
        \\
        ,
    });

    const Case = struct { name: []const u8, include_hetatm: bool, workflow: bool, expected_atoms: usize };
    const cases = [_]Case{
        .{ .name = "batch-default.jsonl", .include_hetatm = false, .workflow = false, .expected_atoms = 1 },
        .{ .name = "batch-hetatm.jsonl", .include_hetatm = true, .workflow = false, .expected_atoms = 2 },
        .{ .name = "workflow-default", .include_hetatm = false, .workflow = true, .expected_atoms = 1 },
        .{ .name = "workflow-hetatm", .include_hetatm = true, .workflow = true, .expected_atoms = 2 },
    };

    for (cases) |case| {
        const output_path = try std.fs.path.join(allocator, &.{ root_path, case.name });
        defer allocator.free(output_path);

        var jsonl_path: []const u8 = undefined;
        if (case.workflow) {
            const workflow = try std.fmt.allocPrint(allocator,
                \\version = 1
                \\kind = "workflow"
                \\
                \\[input]
                \\dir = "{s}"
                \\
                \\[output]
                \\dir = "{s}"
                \\format = "jsonl"
                \\
                \\[calculation]
                \\n_points = 8
                \\quiet = true
                \\include_hetatm = {}
                \\
                \\[classifier]
                \\type = "ccd"
                \\
                \\[[jobs]]
                \\name = "all"
                \\
            , .{ input_dir, output_path, case.include_hetatm });
            defer allocator.free(workflow);
            try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });
            try run(allocator, std.testing.io, .{ .workflow_path = workflow_path });
            jsonl_path = try std.fs.path.join(allocator, &.{ output_path, "all.jsonl" });
        } else {
            try run(allocator, std.testing.io, .{
                .input_path = input_dir,
                .output_path = output_path,
                .output_format = .jsonl,
                .n_threads = 1,
                .n_points = 8,
                .include_hetatm = case.include_hetatm,
                .quiet = true,
                .show_progress = false,
            });
            jsonl_path = try allocator.dupe(u8, output_path);
        }
        defer allocator.free(jsonl_path);

        const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(4096));
        defer allocator.free(content);
        const line = std.mem.trimEnd(u8, content, "\n");
        const parsed = try std.json.parseFromSlice(std.json.Value, allocator, line, .{});
        defer parsed.deinit();
        try std.testing.expectEqualStrings("ok", parsed.value.object.get("status").?.string);
        try std.testing.expectEqual(case.expected_atoms, parsed.value.object.get("atom_areas").?.array.items.len);
    }
}

/// Chain A has CA as altLoc A (0.30) and B (0.70) and a CB with altLoc A
/// only; chain B is one atom. Atoms kept: 5 with `auto`, 6 with `all`, 4 with
/// `B`.
const altloc_workflow_pdb =
    \\ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00 10.00           N
    \\ATOM      2  CA AALA A   1       2.000   0.000   0.000  0.30 10.00           C
    \\ATOM      3  CA BALA A   1       3.000   0.000   0.000  0.70 10.00           C
    \\ATOM      4  CB AALA A   1       4.000   0.000   0.000  0.30 10.00           C
    \\ATOM      5  C   ALA A   1       5.000   0.000   0.000  1.00 10.00           C
    \\ATOM      6  N   GLY B   1       7.000   0.000   0.000  1.00 10.00           N
    \\END
    \\
;
const altloc_workflow_cif =
    \\data_ALTLOC
    \\loop_
    \\_atom_site.group_PDB
    \\_atom_site.type_symbol
    \\_atom_site.label_atom_id
    \\_atom_site.label_alt_id
    \\_atom_site.label_comp_id
    \\_atom_site.label_asym_id
    \\_atom_site.label_seq_id
    \\_atom_site.Cartn_x
    \\_atom_site.Cartn_y
    \\_atom_site.Cartn_z
    \\_atom_site.occupancy
    \\ATOM N N  . ALA A 1 1.000 0.000 0.000 1.00
    \\ATOM C CA A ALA A 1 2.000 0.000 0.000 0.30
    \\ATOM C CA B ALA A 1 3.000 0.000 0.000 0.70
    \\ATOM C CB A ALA A 1 4.000 0.000 0.000 0.30
    \\ATOM C C  . ALA A 1 5.000 0.000 0.000 1.00
    \\ATOM N N  . GLY B 1 7.000 0.000 0.000 1.00
    \\#
    \\
;

/// One way a workflow reaches the input parser.
const AltLocWorkflowPath = struct {
    /// Workflow text after the [calculation] table
    body: []const u8,
    /// Output file inside the output directory
    output_name: []const u8,
    /// Row field with one entry per atom
    atoms_field: []const u8,
};

/// Run `path` over a directory that holds `altloc_workflow_pdb` and
/// `altloc_workflow_cif`, with an optional `[calculation].altloc` value and an
/// optional `--altloc` flag, and return the atom count of the two output rows.
fn altLocWorkflowAtomCount(
    root: []const u8,
    path: AltLocWorkflowPath,
    workflow_altloc: ?[]const u8,
    cli_flag: ?[]const u8,
) !usize {
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();
    const cwd = std.Io.Dir.cwd();

    const input_dir = try std.fs.path.join(arena, &.{ root, "input" });
    try cwd.createDirPath(std.testing.io, input_dir);
    try cwd.writeFile(std.testing.io, .{ .sub_path = try std.fs.path.join(arena, &.{ input_dir, "altloc.pdb" }), .data = altloc_workflow_pdb });
    try cwd.writeFile(std.testing.io, .{ .sub_path = try std.fs.path.join(arena, &.{ input_dir, "altloc.cif" }), .data = altloc_workflow_cif });
    const map_path = try std.fs.path.join(arena, &.{ root, "chains.csv" });
    try cwd.writeFile(std.testing.io, .{ .sub_path = map_path, .data =
        \\filename,chains,asym_id_type
        \\altloc.pdb,"A,B",label
        \\altloc.cif,"A,B",label
        \\
    });

    const output_dir = try std.fmt.allocPrint(arena, "{s}/out-{s}-{s}-{s}", .{
        root,
        path.output_name,
        workflow_altloc orelse "unset",
        if (cli_flag) |flag| flag["--altloc=".len..] else "unset",
    });
    const altloc_line = if (workflow_altloc) |value|
        try std.fmt.allocPrint(arena, "altloc = \"{s}\"\n", .{value})
    else
        "";
    const body = try std.mem.replaceOwned(u8, arena, path.body, "MAP", map_path);
    const workflow = try std.fmt.allocPrint(arena,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[classifier]
        \\type = "naccess"
        \\
        \\[calculation]
        \\threads = 1
        \\n_points = 8
        \\quiet = true
        \\{s}
        \\{s}
    , .{ input_dir, output_dir, altloc_line, body });
    const workflow_path = try std.fs.path.join(arena, &.{ root, "workflow.toml" });
    try cwd.writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    var argv = std.ArrayListUnmanaged([]const u8).empty;
    try argv.appendSlice(arena, &.{ "zsasa", "batch", "--workflow", workflow_path });
    if (cli_flag) |flag| try argv.append(arena, flag);
    try run(std.testing.allocator, std.testing.io, parseArgs(argv.items, 2));

    const content = try cwd.readFileAlloc(std.testing.io, try std.fs.path.join(arena, &.{ output_dir, path.output_name }), arena, .limited(1 << 20));
    var n_rows: usize = 0;
    var n_atoms: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |line| {
        if (line.len == 0) continue;
        const object = (try std.json.parseFromSliceLeaky(std.json.Value, arena, line, .{})).object;
        try std.testing.expectEqualStrings("ok", object.get("status").?.string);
        const row_atoms = object.get(path.atoms_field).?.array.items.len;
        // The PDB and the mmCIF file hold the same atoms
        if (n_rows > 0) try std.testing.expectEqual(n_atoms, row_atoms);
        n_atoms = row_atoms;
        n_rows += 1;
    }
    try std.testing.expectEqual(@as(usize, 2), n_rows);
    return n_atoms;
}

test "workflow honors --altloc and the calculation altloc key in every path" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root = root_buf[0..try tmp_dir.dir.realPath(std.testing.io, &root_buf)];

    const paths = [_]AltLocWorkflowPath{
        // File-first: jobs without chain maps share one parse of each file
        .{
            .body =
            \\[[jobs]]
            \\name = "file_first"
            \\
            ,
            .output_name = "file_first.jsonl",
            .atoms_field = "atom_areas",
        },
        // Job-first: a job whose auth_chain differs from the shared setting
        .{
            .body =
            \\[[jobs]]
            \\name = "job_first"
            \\auth_chain = true
            \\
            ,
            .output_name = "job_first.jsonl",
            .atoms_field = "atom_areas",
        },
        // Selection map: a job with a chain map and JSONL output
        .{
            .body =
            \\[[jobs]]
            \\name = "selection_map"
            \\chain_map = "MAP"
            \\
            ,
            .output_name = "selection_map.jsonl",
            .atoms_field = "atom_areas",
        },
        // BSA analysis
        .{
            .body =
            \\[analysis]
            \\type = "bsa"
            \\name = "bsa"
            \\partner_a = ["A"]
            \\partner_b = ["B"]
            \\level = "residue"
            \\atom_output = true
            \\
            ,
            .output_name = "bsa.jsonl",
            .atoms_field = "atom_delta_sasa",
        },
    };

    for (paths) |path| {
        // Default: auto
        try std.testing.expectEqual(@as(usize, 5), try altLocWorkflowAtomCount(root, path, null, null));
        // The command-line flag is not dropped
        try std.testing.expectEqual(@as(usize, 4), try altLocWorkflowAtomCount(root, path, null, "--altloc=B"));
        try std.testing.expectEqual(@as(usize, 6), try altLocWorkflowAtomCount(root, path, null, "--altloc=all"));
        // The workflow key
        try std.testing.expectEqual(@as(usize, 6), try altLocWorkflowAtomCount(root, path, "all", null));
        try std.testing.expectEqual(@as(usize, 4), try altLocWorkflowAtomCount(root, path, "B", null));
        // The flag takes precedence over the key, also when it is `auto`
        try std.testing.expectEqual(@as(usize, 4), try altLocWorkflowAtomCount(root, path, "all", "--altloc=B"));
        try std.testing.expectEqual(@as(usize, 5), try altLocWorkflowAtomCount(root, path, "all", "--altloc=auto"));
    }
}

test "workflow altloc applies to the batch config unless --altloc is given" {
    const calculation = workflow_manifest.Calculation{ .altloc = .{ .mode = .selected, .id = 'B' } };

    var from_workflow = BatchConfig{};
    const no_flag = parseArgs(&.{ "zsasa", "batch", "--workflow", "wf.toml" }, 2);
    try applyWorkflowToBatchConfig(&from_workflow, no_flag, calculation, .{}, .{});
    applyCliOverrides(&from_workflow, no_flag);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.selected, from_workflow.alt_loc_mode);
    try std.testing.expectEqual(@as(u8, 'B'), from_workflow.alt_loc_id);

    var from_flag = BatchConfig{};
    const flag = parseArgs(&.{ "zsasa", "batch", "--workflow", "wf.toml", "--altloc=C" }, 2);
    try applyWorkflowToBatchConfig(&from_flag, flag, calculation, .{}, .{});
    applyCliOverrides(&from_flag, flag);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.selected, from_flag.alt_loc_mode);
    try std.testing.expectEqual(@as(u8, 'C'), from_flag.alt_loc_id);

    var flag_only = BatchConfig{};
    const none_flag = parseArgs(&.{ "zsasa", "batch", "--workflow", "wf.toml", "--altloc=none" }, 2);
    try applyWorkflowToBatchConfig(&flag_only, none_flag, .{}, .{}, .{});
    applyCliOverrides(&flag_only, none_flag);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.none, flag_only.alt_loc_mode);
}

test "runBatchParallel writes JSONL error rows for failed files" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root_path = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root_path, "input" });
    defer allocator.free(input_dir);
    const output_path = try std.fs.path.join(allocator, &.{ root_path, "results.jsonl" });
    defer allocator.free(output_path);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const good_pdb =
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n" ++
        "ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 20.00           C\n" ++
        "END\n";
    const good_path = try std.fs.path.join(allocator, &.{ input_dir, "good.pdb" });
    defer allocator.free(good_path);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = good_path, .data = good_pdb });

    const bad_path = try std.fs.path.join(allocator, &.{ input_dir, "bad.pdb" });
    defer allocator.free(bad_path);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = bad_path, .data = "not a pdb file\n" });

    var result = try runBatchParallel(allocator, std.testing.io, input_dir, null, .{
        .n_threads = 2,
        .algorithm = .sr,
        .n_points = 8,
        .quiet = true,
        .show_progress = false,
        .output_format = .jsonl,
        .store_atom_areas = true,
        .classifier_type = .naccess,
    }, output_path);
    defer result.deinit();

    try std.testing.expectEqual(@as(usize, 2), result.total_files);
    try std.testing.expectEqual(@as(usize, 1), result.successful);
    try std.testing.expectEqual(@as(usize, 1), result.failed);

    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, output_path, allocator, .limited(4096));
    defer allocator.free(content);
    try std.testing.expectEqual(@as(usize, 2), std.mem.count(u8, content, "\n"));
    try std.testing.expect(std.mem.indexOf(u8, content, "\"status\":\"ok\",\"filename\":\"good.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"status\":\"err\",\"filename\":\"bad.pdb\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"error\":") != null);
}

test "batch on an empty directory leaves an empty JSONL file whatever the thread count" {
    const allocator = std.testing.allocator;
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    const jsonl_path = try sandbox.path("results.jsonl");
    defer allocator.free(jsonl_path);

    inline for (test_naming_threads) |n_threads| {
        // Left over from an earlier run
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = jsonl_path, .data = "{\"status\":\"ok\",\"filename\":\"stale.pdb\"}\n" });

        var config = NamingSandbox.config(n_threads);
        config.output_format = .jsonl;
        config.store_atom_areas = true;
        var result = try runBatch(allocator, std.testing.io, sandbox.input_dir, null, config, jsonl_path);
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 0), result.total_files);
        try std.testing.expectEqual(@as(usize, 0), result.failed);

        const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(4096));
        defer allocator.free(content);
        try std.testing.expectEqualStrings("", content);
    }

    // Per-file output: both runners create the (empty) output directory.
    inline for (test_naming_threads) |n_threads| {
        var result = try sandbox.run(n_threads, std.fmt.comptimePrint("out{d}", .{n_threads}));
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 0), result.total_files);
    }
    try sandbox.expectTree(&.{ "input/", "out1/", "out4/", "results.jsonl" });
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

test "BatchArgs --lr-trig defaults to exact" {
    const args = [_][]const u8{ "zsasa", "batch", "--algorithm=lr", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(TrigMode.exact, parsed.lr_trig);
    try std.testing.expectEqual(false, parsed.lr_trig_explicit);
}

test "BatchArgs --lr-trig=MODE and --lr-trig MODE" {
    const eq = [_][]const u8{ "zsasa", "batch", "--algorithm=lr", "--lr-trig=fast", "input_dir/" };
    const parsed_eq = parseArgs(&eq, 2);
    try std.testing.expectEqual(TrigMode.fast, parsed_eq.lr_trig);
    try std.testing.expectEqual(true, parsed_eq.lr_trig_explicit);

    const spaced = [_][]const u8{ "zsasa", "batch", "--lr-trig", "fast", "input_dir/" };
    const parsed_spaced = parseArgs(&spaced, 2);
    try std.testing.expectEqual(TrigMode.fast, parsed_spaced.lr_trig);
    try std.testing.expectEqual(true, parsed_spaced.lr_trig_explicit);
    try std.testing.expectEqualStrings("input_dir/", parsed_spaced.input_path.?);

    const exact = [_][]const u8{ "zsasa", "batch", "--lr-trig=exact", "input_dir/" };
    const parsed_exact = parseArgs(&exact, 2);
    try std.testing.expectEqual(TrigMode.exact, parsed_exact.lr_trig);
    try std.testing.expectEqual(true, parsed_exact.lr_trig_explicit);
}

test "workflow lr_trig applies to the batch config unless --lr-trig was given, and rejects unknown values" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const Calculation = @import("workflow_manifest.zig").Calculation;
    const output = @import("workflow_manifest.zig").Output{};
    const classifier_config = @import("workflow_manifest.zig").ClassifierConfig{};

    {
        var config = BatchConfig{};
        try applyWorkflowToBatchConfig(&config, .{}, Calculation{ .lr_trig = "fast" }, output, classifier_config);
        try std.testing.expectEqual(TrigMode.fast, config.lr_trig);
    }
    {
        // The CLI value wins: the workflow value is skipped and the override is applied.
        var config = BatchConfig{};
        const args = BatchArgs{ .lr_trig = .exact, .lr_trig_explicit = true };
        try applyWorkflowToBatchConfig(&config, args, Calculation{ .lr_trig = "fast" }, output, classifier_config);
        applyCliOverrides(&config, args);
        try std.testing.expectEqual(TrigMode.exact, config.lr_trig);
    }
    {
        var config = BatchConfig{};
        const args = BatchArgs{ .lr_trig = .fast, .lr_trig_explicit = true };
        try applyWorkflowToBatchConfig(&config, args, Calculation{}, output, classifier_config);
        applyCliOverrides(&config, args);
        try std.testing.expectEqual(TrigMode.fast, config.lr_trig);
    }
    {
        var config = BatchConfig{};
        try std.testing.expectError(
            error.InvalidArgument,
            applyWorkflowToBatchConfig(&config, .{}, Calculation{ .lr_trig = "approximate" }, output, classifier_config),
        );
    }
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

test "BatchArgs --jsonl-decimals" {
    const args = [_][]const u8{ "zsasa", "batch", "--format=jsonl", "--jsonl-decimals=3", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(?u8, 3), parsed.jsonl_decimals);
    try std.testing.expectEqual(true, parsed.jsonl_decimals_explicit);
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
    try std.testing.expectEqual(@as(f64, 0.10), parsed.adaptive_low);
    try std.testing.expectEqual(@as(f64, 0.90), parsed.adaptive_high);
}

test "BatchArgs --use-bitmask" {
    const args = [_][]const u8{ "zsasa", "batch", "--use-bitmask", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_bitmask);
}

test "BatchArgs --bitmask-correction-coeff" {
    const args = [_][]const u8{ "zsasa", "batch", "--use-bitmask", "--bitmask-correction-coeff=0.2", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.bitmask_correction);
    try std.testing.expectEqual(@as(f64, 0.2), parsed.bitmask_correction_coeff);
}

test "BatchConfig bitmask correction requires bitmask" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    try std.testing.expectError(
        error.InvalidArgument,
        validateBitmaskCorrectionConfig(.{ .bitmask_correction = true }),
    );
}

test "BatchConfig bitmask correction rejects adaptive SR" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    try std.testing.expectError(
        error.InvalidArgument,
        validateBitmaskCorrectionConfig(.{
            .use_bitmask = true,
            .bitmask_correction = true,
            .adaptive_sr = true,
        }),
    );
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

test "BatchArgs --af-model-fast" {
    const args = [_][]const u8{ "zsasa", "batch", "--af-model-fast", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expect(parsed.af_model_fast);
}

test "BatchArgs --input-io" {
    const args = [_][]const u8{ "zsasa", "batch", "--input-io=read", "input_dir/" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(InputIoMode.read, parsed.input_io);
    try std.testing.expect(parsed.input_io_explicit);
}

test "readInputFile preserves AF model atom metadata when fast parsing is requested" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const source =
        \\data_AF_FAST
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_alt_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_ins_code
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM 1 N N  . GLY A 1 1 ? 1.000 2.000 3.000
        \\ATOM 2 C CA . GLY A 1 1 ? 2.000 3.000 4.000
        \\ATOM 3 N N  . ALA B 2 1 ? 4.000 5.000 6.000
        \\ATOM 4 C CA . ALA B 2 1 ? 5.000 6.000 7.000
        \\
    ;

    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "af.cif", .data = source });
    const path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "af.cif", allocator);
    defer allocator.free(path);

    var parsed = try readInputFile(allocator, std.testing.io, path, .{
        .af_model_fast = true,
        .classifier_type = null,
    });
    defer parsed.deinit();

    try std.testing.expectEqual(@as(usize, 4), parsed.input.atomCount());
    try std.testing.expectEqualStrings("GLY", parsed.input.residue.?[0].slice());
    try std.testing.expectEqualStrings("CA", parsed.input.atom_name.?[3].slice());
    try std.testing.expectEqualStrings("A", parsed.input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("B", parsed.input.chain_id.?[2].slice());
    try std.testing.expectEqual(@as(i32, 1), parsed.input.residue_num.?[0]);
}

test "AF model fast fallback is limited to unsupported layouts" {
    try std.testing.expect(shouldFallbackAfModelFastError(error.UnsupportedLayout));
    try std.testing.expect(!shouldFallbackAfModelFastError(error.InvalidCoordinate));
    try std.testing.expect(!shouldFallbackAfModelFastError(error.AccessDenied));
}

test "batch --altloc applies to PDB input as it does to mmCIF input" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "altloc.pdb", .data =
        \\ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00 10.00           N
        \\ATOM      2  CA AALA A   1       2.000   0.000   0.000  0.30 10.00           C
        \\ATOM      3  CA BALA A   1       3.000   0.000   0.000  0.70 10.00           C
        \\ATOM      4  C   ALA A   1       4.000   0.000   0.000  1.00 10.00           C
        \\END
        \\
    });
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "altloc.cif", .data =
        \\data_ALTLOC
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_alt_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\_atom_site.occupancy
        \\ATOM N N  . ALA A 1 1.000 0.000 0.000 1.00
        \\ATOM C CA A ALA A 1 2.000 0.000 0.000 0.30
        \\ATOM C CA B ALA A 1 3.000 0.000 0.000 0.70
        \\ATOM C C  . ALA A 1 4.000 0.000 0.000 1.00
        \\#
        \\
    });
    const pdb_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "altloc.pdb", allocator);
    defer allocator.free(pdb_path);
    const cif_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "altloc.cif", allocator);
    defer allocator.free(cif_path);

    // CA has the alternates A (0.30) and B (0.70)
    const Case = struct { flag: []const u8, x: []const f64 };
    const cases = [_]Case{
        .{ .flag = "--altloc=auto", .x = &.{ 1, 2, 4 } },
        .{ .flag = "--altloc=all", .x = &.{ 1, 2, 3, 4 } },
        .{ .flag = "--altloc=A", .x = &.{ 1, 2, 4 } },
        .{ .flag = "--altloc=B", .x = &.{ 1, 3, 4 } },
        .{ .flag = "--altloc=C", .x = &.{ 1, 4 } },
        .{ .flag = "--altloc=highest-occupancy", .x = &.{ 1, 3, 4 } },
    };
    for ([_][]const u8{ pdb_path, cif_path }) |path| {
        for (cases) |case| {
            const args = parseArgs(&.{ "zsasa", "batch", case.flag, "input_dir/" }, 2);
            var parsed = try readInputFile(allocator, std.testing.io, path, .{
                .alt_loc_mode = args.alt_loc_mode,
                .alt_loc_id = args.alt_loc_id,
            });
            defer parsed.deinit();
            try std.testing.expectEqualSlices(f64, case.x, parsed.input.x);
        }

        // `none` exists to fail fast
        try std.testing.expectError(
            error.UnexpectedAltLoc,
            readInputFile(allocator, std.testing.io, path, .{ .alt_loc_mode = .none }),
        );
    }
}

test "readInputFile falls back to generic mmCIF for unsupported fast layout" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    const source =
        \\data_GENERIC
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM C CA GLY 1.0 2.0 3.0
        \\#
    ;

    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "generic.cif", .data = source });
    const path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "generic.cif", allocator);
    defer allocator.free(path);

    var parsed = try readInputFile(allocator, std.testing.io, path, .{
        .af_model_fast = true,
        .classifier_type = null,
    });
    defer parsed.deinit();

    try std.testing.expectEqual(@as(usize, 1), parsed.input.atomCount());
    try std.testing.expectEqualStrings("GLY", parsed.input.residue.?[0].slice());
    try std.testing.expectEqualStrings("CA", parsed.input.atom_name.?[0].slice());
}

test "AF model fast parser is skipped when chain filtering is requested" {
    const chains = [_][]const u8{"Z"};
    try std.testing.expect(!shouldTryAfModelFastParser(.{
        .af_model_fast = true,
        .chain_filter = chains[0..],
    }));
}

test "AF model fast parser leaves files to the generic parser when HETATM is included" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    // ATOM rows followed by HETATM rows: the fast parser stops at the first
    // row that is not ATOM.
    const source =
        \\data_AF_HETATM
        \\loop_
        \\_atom_site.group_PDB
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_alt_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_entity_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_ins_code
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\ATOM   1 N  N  . GLY A 1 1 ? 1.000 2.000 3.000
        \\ATOM   2 C  CA . GLY A 1 1 ? 2.000 3.000 4.000
        \\HETATM 3 ZN ZN . ZN  B 2 . ? 9.000 9.000 9.000
        \\HETATM 4 O  O  . HOH C 3 . ? 12.00 12.00 12.00
        \\
    ;
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "af.cif", .data = source });
    const path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "af.cif", allocator);
    defer allocator.free(path);

    const Case = struct { af_model_fast: bool, include_hetatm: bool, expected_atoms: usize };
    const cases = [_]Case{
        .{ .af_model_fast = false, .include_hetatm = false, .expected_atoms = 2 },
        .{ .af_model_fast = false, .include_hetatm = true, .expected_atoms = 4 },
        .{ .af_model_fast = true, .include_hetatm = false, .expected_atoms = 2 },
        .{ .af_model_fast = true, .include_hetatm = true, .expected_atoms = 4 },
    };
    for (cases) |case| {
        var parsed = try readInputFile(allocator, std.testing.io, path, .{
            .af_model_fast = case.af_model_fast,
            .include_hetatm = case.include_hetatm,
            .classifier_type = null,
        });
        defer parsed.deinit();
        try std.testing.expectEqual(case.expected_atoms, parsed.input.atomCount());
    }

    try std.testing.expect(shouldTryAfModelFastParser(.{ .af_model_fast = true }));
    try std.testing.expect(!shouldTryAfModelFastParser(.{ .af_model_fast = true, .include_hetatm = true }));
}

test "BatchConfig progress respects show_progress and quiet" {
    try std.testing.expect(shouldShowProgress(.{}));
    try std.testing.expect(!shouldShowProgress(.{ .show_progress = false }));
    try std.testing.expect(!shouldShowProgress(.{ .quiet = true }));
}

const test_topology_header_rest = "  zsasa\n\n";

// Heavy atoms of acetonitrile (radii 1.88, 1.61, 1.64) and of acetaldehyde
// (1.88, 1.76, 1.42), as V2000 records without the title line. Both have
// atoms C1 and C2, with a different radius for C2; the element fallback
// would give 1.70 to every carbon. The atoms are 50 A apart, so each one is
// fully exposed and its area gives its radius.
const test_topology_acetonitrile = test_topology_header_rest ++
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n" ++
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "   50.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  100.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  1  2  1  0  0  0  0\n  2  3  3  0  0  0  0\n" ++
    "M  END\n$$$$\n";
const test_topology_acetaldehyde = test_topology_header_rest ++
    "  3  2  0  0  0  0  0  0  0  0999 V2000\n" ++
    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "   50.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  100.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" ++
    "  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n" ++
    "M  END\n$$$$\n";

/// Expects the per-molecule CSV output `name` below the sandbox root to hold
/// one fully exposed atom for each of `radii`.
fn expectExposedAtomRadii(sandbox: NamingSandbox, name: []const u8, radii: []const f64) !void {
    const allocator = std.testing.allocator;
    const file_path = try sandbox.path(name);
    defer allocator.free(file_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, file_path, allocator, .limited(64 * 1024));
    defer allocator.free(content);

    var lines = std.mem.tokenizeScalar(u8, content, '\n');
    try std.testing.expectEqualStrings("atom_index,area", lines.next().?);
    for (radii) |radius| {
        const line = lines.next().?;
        const area = try std.fmt.parseFloat(f64, line[std.mem.findScalar(u8, line, ',').? + 1 ..]);
        const probe_radius = 1.4;
        const exposed = 4.0 * std.math.pi * (radius + probe_radius) * (radius + probe_radius);
        try std.testing.expectApproxEqAbs(exposed, area, 1e-5);
    }
    try std.testing.expect(std.mem.startsWith(u8, lines.next().?, "total,"));
    try std.testing.expectEqual(@as(?[]const u8, null), lines.next());
}

test "batch classifies every SDF molecule from its own bond topology, with or without a title" {
    const allocator = std.testing.allocator;
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();

    // The same two molecules with titles, without a title, with the same
    // title, with a title that is a residue name of the built-in table, and
    // without a title next to a titled molecule
    try sandbox.writeInput("named.sdf", "nitrile\n" ++ test_topology_acetonitrile ++ "aldehyde\n" ++ test_topology_acetaldehyde);
    try sandbox.writeInput("unnamed.sdf", "\n" ++ test_topology_acetonitrile ++ "\n" ++ test_topology_acetaldehyde);
    try sandbox.writeInput("same.sdf", "lig\n" ++ test_topology_acetonitrile ++ "lig\n" ++ test_topology_acetaldehyde);
    try sandbox.writeInput("residue.sdf", "A\n" ++ test_topology_acetonitrile ++ "ALA\n" ++ test_topology_acetaldehyde);
    try sandbox.writeInput("mixed.sdf", "\n" ++ test_topology_acetonitrile ++ "aldehyde\n" ++ test_topology_acetaldehyde ++ "\n" ++ test_topology_acetonitrile);

    const nitrile_outputs = [_][]const u8{ "named_nitrile.csv", "unnamed_1.csv", "same_lig_1.csv", "residue_A.csv", "mixed_1.csv", "mixed_3.csv" };
    const aldehyde_outputs = [_][]const u8{ "named_aldehyde.csv", "unnamed_2.csv", "same_lig_2.csv", "residue_ALA.csv", "mixed_aldehyde.csv" };

    inline for (test_naming_threads) |n_threads| {
        const out_name = std.fmt.comptimePrint("csv{d}", .{n_threads});
        const output_dir = try sandbox.path(out_name);
        defer allocator.free(output_dir);

        var config = NamingSandbox.config(n_threads);
        config.output_format = .csv;
        var result = try runBatch(allocator, std.testing.io, sandbox.input_dir, output_dir, config, null);
        defer result.deinit();
        try std.testing.expectEqual(@as(usize, 11), result.total_files);
        try std.testing.expectEqual(@as(usize, 11), result.successful);

        // The radii of each molecule's bond table, not the element fallback
        inline for (nitrile_outputs) |name| {
            try expectExposedAtomRadii(sandbox, out_name ++ "/" ++ name, &.{ 1.88, 1.61, 1.64 });
        }
        inline for (aldehyde_outputs) |name| {
            try expectExposedAtomRadii(sandbox, out_name ++ "/" ++ name, &.{ 1.88, 1.76, 1.42 });
        }
    }
}

// -----------------------------------------------------------------------------
// Keys and options a batch workflow does not read
// -----------------------------------------------------------------------------

const test_bsa_analysis = "[analysis]\ntype = \"bsa\"\npartner_a = [\"A\"]\npartner_b = [\"B\"]\n";

/// Write "workflow.toml" in the sandbox: a JSONL workflow over its input
/// directory with the given extra lines in `[input]` and `[calculation]`,
/// followed by `body`. Returns the path, allocated from `arena`.
fn writeKeyWorkflow(
    sandbox: NamingSandbox,
    arena: Allocator,
    input_extra: []const u8,
    calc_extra: []const u8,
    body: []const u8,
) ![]const u8 {
    const workflow = try std.fmt.allocPrint(arena,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\{s}
        \\[output]
        \\dir = "{s}/output"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 8
        \\quiet = true
        \\{s}
        \\[classifier]
        \\type = "naccess"
        \\
        \\{s}
    , .{ sandbox.input_dir, input_extra, sandbox.root, calc_extra, body });
    const workflow_path = try std.fs.path.join(arena, &.{ sandbox.root, "workflow.toml" });
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });
    return workflow_path;
}

fn readSandboxFile(sandbox: NamingSandbox, arena: Allocator, name: []const u8) ![]const u8 {
    const path = try std.fs.path.join(arena, &.{ sandbox.root, name });
    return std.Io.Dir.cwd().readFileAlloc(std.testing.io, path, arena, .limited(1 << 20));
}

test "workflow rejects keys that batch cannot honor before running anything" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("two.pdb", test_two_chain_pdb);
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    const Case = struct { input_extra: []const u8 = "", calc_extra: []const u8 = "" };
    const cases = [_]Case{
        .{ .input_extra = "model = 2\n" },
        .{ .input_extra = "mol = \"1\"\n" },
        .{ .input_extra = "path = \"two.pdb\"\n" },
        .{ .calc_extra = "rsa = true\n" },
        .{ .calc_extra = "per_residue = true\n" },
        .{ .calc_extra = "polar = true\n" },
        .{ .calc_extra = "validate_only = true\n" },
    };
    for (cases) |case| {
        for (test_workflow_runners) |runner| {
            const body = try std.fmt.allocPrint(arena, "[[jobs]]\nname = \"everything\"\n{s}", .{runner.jobOption()});
            const workflow_path = try writeKeyWorkflow(sandbox, arena, case.input_extra, case.calc_extra, body);
            try std.testing.expectError(
                error.InvalidArgument,
                runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path }),
            );
        }
        // The same keys in an [analysis] workflow
        const workflow_path = try writeKeyWorkflow(sandbox, arena, case.input_extra, case.calc_extra, test_bsa_analysis);
        try std.testing.expectError(
            error.InvalidArgument,
            runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path }),
        );
    }
    try sandbox.expectTree(&.{ "input/", "input/two.pdb", "workflow.toml" });
}

test "workflow rejects [input] chain where no job can use it" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("two.pdb", test_two_chain_pdb);
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    const bodies = [_][]const u8{
        test_bsa_analysis,
        "[[jobs]]\nname = \"m\"\nchain_map = \"chains.csv\"\n",
        "[[jobs]]\nname = \"own\"\nchains = [\"A\"]\n",
    };
    for (bodies) |body| {
        const workflow_path = try writeKeyWorkflow(sandbox, arena, "chain = \"B\"\n", "", body);
        try std.testing.expectError(
            error.InvalidArgument,
            runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path }),
        );
    }
    try sandbox.expectTree(&.{ "input/", "input/two.pdb", "workflow.toml" });
}

test "workflow honors [input] chain as the default chains of the jobs that have none" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("two.pdb", test_two_chain_pdb);
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    for (test_workflow_runners) |runner| {
        try sandbox.tmp.dir.deleteTree(std.testing.io, "output");
        // "default" takes chain B from [input]; "explicit" lists it itself,
        // "both" lists A and B.
        const body = try std.fmt.allocPrint(arena,
            \\[[jobs]]
            \\name = "default"
            \\
            \\[[jobs]]
            \\name = "explicit"
            \\chains = ["B"]
            \\{s}
            \\[[jobs]]
            \\name = "both"
            \\chains = ["A", "B"]
            \\
        , .{runner.jobOption()});
        const workflow_path = try writeKeyWorkflow(sandbox, arena, "chain = \"B\"\n", "", body);
        try runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path });

        const default_rows = try readSandboxFile(sandbox, arena, "output/default.jsonl");
        const explicit_rows = try readSandboxFile(sandbox, arena, "output/explicit.jsonl");
        const both_rows = try readSandboxFile(sandbox, arena, "output/both.jsonl");
        try std.testing.expectEqualStrings(explicit_rows, default_rows);
        try std.testing.expect(!std.mem.eql(u8, both_rows, default_rows));
        // Chain B has two atoms, the whole structure four
        try std.testing.expectEqual(@as(usize, 1), std.mem.count(u8, default_rows, "\"atom_areas\":["));
        try std.testing.expect(std.mem.count(u8, default_rows, ",") < std.mem.count(u8, both_rows, ","));
    }
}

test "workflow with timing or an output path still runs and only warns" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("two.pdb", test_two_chain_pdb);
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    for (test_workflow_runners) |runner| {
        try sandbox.tmp.dir.deleteTree(std.testing.io, "output");
        const body = try std.fmt.allocPrint(arena, "[[jobs]]\nname = \"everything\"\n{s}", .{runner.jobOption()});
        const workflow_path = try writeKeyWorkflow(sandbox, arena, "", "timing = true\n", body);
        try runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path });
        try sandbox.expectTree(&.{
            "input/",
            "input/two.pdb",
            "output/",
            "output/everything.jsonl",
            "workflow.toml",
        });
    }
}

fn expectBatchFinding(
    content: []const u8,
    args: BatchArgs,
    severity: workflow_manifest.Severity,
    needle: ?[]const u8,
) !void {
    var workflow = try workflow_manifest.parse(std.testing.allocator, content);
    defer workflow.deinit();
    const findings = batchWorkflowFindings(workflow, args);
    if (needle) |text| {
        for (findings.slice()) |finding| {
            if (finding.severity == severity and std.mem.find(u8, finding.message, text) != null) return;
        }
        std.debug.print("no {s} containing '{s}'\n", .{ @tagName(severity), text });
        return error.TestExpectedFinding;
    }
    try std.testing.expectEqual(@as(usize, 0), findings.len);
}

test "batch workflow findings: timing and stage profiling only warn" {
    const jobs = "version = 1\n[input]\ndir = \"d\"\n[[jobs]]\nname = \"j\"\n";
    const jobs_timing = "version = 1\n[input]\ndir = \"d\"\n[calculation]\ntiming = true\n[[jobs]]\nname = \"j\"\n";
    const analysis_timing = "version = 1\n[input]\ndir = \"d\"\n[calculation]\ntiming = true\n" ++ test_bsa_analysis;

    try expectBatchFinding(jobs, .{}, .warning, null);
    try expectBatchFinding(jobs_timing, .{}, .warning, "[calculation] timing = true");
    try expectBatchFinding(jobs, .{ .show_timing = true, .timing_explicit = true }, .warning, "--timing");
    // The command-line option is named when both are given
    try expectBatchFinding(jobs_timing, .{ .show_timing = true, .timing_explicit = true }, .warning, "--timing");
    // An [analysis] workflow reports its SASA time
    try expectBatchFinding(analysis_timing, .{ .show_timing = true, .timing_explicit = true }, .warning, null);
    for ([_][]const u8{ jobs, analysis_timing }) |content| {
        try expectBatchFinding(content, .{ .profile_stages = true, .profile_stages_explicit = true }, .warning, "--profile-stages");
    }
}

test "batch workflow findings: an [analysis] workflow rejects --residue-map and a format other than jsonl" {
    const analysis = "version = 1\n[input]\ndir = \"d\"\n" ++ test_bsa_analysis;
    try expectBatchFinding(analysis, .{ .residue_map = true }, .err, "--residue-map");
    try expectBatchFinding(analysis, .{ .output_format = .csv, .format_explicit = true }, .err, "--format");
    try expectBatchFinding(analysis, .{ .output_format = .jsonl, .format_explicit = true }, .err, null);
    // Workflows with jobs take both
    const jobs = "version = 1\n[input]\ndir = \"d\"\n[[jobs]]\nname = \"j\"\n";
    try expectBatchFinding(jobs, .{ .residue_map = true, .output_format = .csv, .format_explicit = true }, .err, null);
}

test "workflow rejects [analysis] with [[jobs]] and a [[jobs]] table without a name" {
    var sandbox = try NamingSandbox.init();
    defer sandbox.deinit();
    try sandbox.writeInput("two.pdb", test_two_chain_pdb);
    var arena_state = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena_state.deinit();
    const arena = arena_state.allocator();

    const cases = [_]struct { body: []const u8, expected: anyerror }{
        .{ .body = test_bsa_analysis ++ "[[jobs]]\nname = \"j\"\n", .expected = error.AnalysisWithJobs },
        .{ .body = "[[jobs]]\n[[jobs]]\nname = \"second\"\n", .expected = error.MissingJobName },
        .{ .body = "[[jobs]]\nname = \"first\"\n[[jobs]]\n", .expected = error.MissingJobName },
    };
    for (cases) |case| {
        const workflow_path = try writeKeyWorkflow(sandbox, arena, "", "", case.body);
        try std.testing.expectError(
            case.expected,
            runWorkflow(std.testing.allocator, std.testing.io, .{ .workflow_path = workflow_path }),
        );
    }
    try sandbox.expectTree(&.{ "input/", "input/two.pdb", "workflow.toml" });
}

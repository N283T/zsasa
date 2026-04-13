// Trajectory analysis mode
// Calculates SASA for each frame in an XTC or DCD trajectory
//
// Supports two parallelism strategies:
//   - Sequential (1 thread): processes frames one at a time
//   - Batch parallel (N threads): reads batches of frames, distributes across threads
//
const std = @import("std");
const xtc = @import("zxdrfile").xtc;
const dcd = @import("dcd.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const pdb_parser = @import("pdb_parser.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const classifier = @import("classifier.zig");
const classifier_naccess = @import("classifier_naccess.zig");
// classifier_protor removed — ProtOr is now an alias for CCD
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const gzip = @import("gzip.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
const Config = types.Config;
const Configf32 = types.Configf32;
const Precision = types.Precision;
const ClassifierType = classifier.ClassifierType;

/// Algorithm selection
pub const Algorithm = enum {
    sr, // Shrake-Rupley
    lr, // Lee-Richards
};

/// Trajectory file format
pub const TrajectoryFormat = enum {
    xtc, // GROMACS XTC (compressed, coordinates in nm)
    dcd, // NAMD/CHARMM DCD (uncompressed, coordinates in Angstroms)
};

/// Detect trajectory format from file extension
pub fn detectTrajectoryFormat(path: []const u8) ?TrajectoryFormat {
    if (std.mem.endsWith(u8, path, ".xtc")) {
        return .xtc;
    } else if (std.mem.endsWith(u8, path, ".dcd")) {
        return .dcd;
    }
    return null;
}

/// Common frame type for trajectory readers
const TrajectoryFrame = struct {
    step: i32,
    time: f32,
    coords: []f32, // flat array of x,y,z coordinates (length = natoms * 3)

    /// The underlying frame owns the coords allocation.
    /// Track which reader type produced it so we can free correctly.
    _xtc_frame: ?xtc.XtcFrame = null,
    _dcd_frame: ?dcd.DcdFrame = null,

    fn deinit(self: *TrajectoryFrame, allocator: Allocator) void {
        if (self._xtc_frame) |*f| f.deinit(allocator);
        if (self._dcd_frame) |*f| f.deinit(allocator);
    }
};

/// Unified trajectory reader wrapping either XTC or DCD readers
const TrajectoryReader = struct {
    xtc_reader: ?*xtc.XtcReader = null,
    dcd_reader: ?*dcd.DcdReader = null,
    format: TrajectoryFormat,

    /// Coordinate scale factor: 10.0 for XTC (nm→Å), 1.0 for DCD (already Å)
    fn coordScale(self: *const TrajectoryReader) f64 {
        return switch (self.format) {
            .xtc => 10.0,
            .dcd => 1.0,
        };
    }

    /// Read next frame, returning null-like via EndOfFile error
    fn readFrame(self: *TrajectoryReader) !TrajectoryFrame {
        switch (self.format) {
            .xtc => {
                const frame = try self.xtc_reader.?.readFrame();
                return TrajectoryFrame{
                    .step = frame.step,
                    .time = frame.time,
                    .coords = frame.coords,
                    ._xtc_frame = frame,
                };
            },
            .dcd => {
                const frame = try self.dcd_reader.?.readFrame();
                return TrajectoryFrame{
                    .step = frame.step,
                    .time = frame.time,
                    .coords = frame.coords,
                    ._dcd_frame = frame,
                };
            },
        }
    }

    /// Check if a read error is an EOF condition
    fn isEof(self: *const TrajectoryReader, err: anyerror) bool {
        return switch (self.format) {
            .xtc => err == xtc.XtcError.EndOfFile,
            .dcd => err == dcd.DcdError.EndOfFile,
        };
    }
};

/// Trajectory mode arguments
pub const TrajArgs = struct {
    traj_path: ?[]const u8 = null,
    topology_path: ?[]const u8 = null,
    output_path: []const u8 = "traj_sasa.csv",
    algorithm: Algorithm = .sr,
    n_threads: usize = 0,
    probe_radius: f64 = 1.4,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    precision: Precision = .f32, // Default f32 for trajectory (speed)
    classifier_type: ?ClassifierType = .naccess, // Default: NACCESS for trajectories (supports explicit H)
    ccd_path: ?[]const u8 = null, // External CCD dictionary file (.zsdc or .cif[.gz])
    stride: u32 = 1, // Process every Nth frame
    start_frame: u32 = 0, // Start frame
    end_frame: ?u32 = null, // End frame (null = all)
    include_hydrogens: bool = true, // Include hydrogen atoms (default: include for MD trajectories)
    batch_size: u32 = 0, // Frames per batch for parallel processing (0 = auto)
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 1..1024)
    quiet: bool = false,
    show_help: bool = false,
};

/// Parse trajectory mode arguments
pub fn parseArgs(args: []const []const u8, start_idx: usize) TrajArgs {
    var result = TrajArgs{};
    var i: usize = start_idx;
    var positional_count: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.startsWith(u8, arg, "--")) {
            // Options
            if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
                result.show_help = true;
            } else if (std.mem.startsWith(u8, arg, "--algorithm=")) {
                const value = arg["--algorithm=".len..];
                result.algorithm = parseAlgorithm(value);
            } else if (std.mem.startsWith(u8, arg, "--threads=")) {
                const value = arg["--threads=".len..];
                result.n_threads = std.fmt.parseInt(usize, value, 10) catch {
                    std.debug.print("Error: Invalid thread count: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "--probe-radius=")) {
                const value = arg["--probe-radius=".len..];
                result.probe_radius = std.fmt.parseFloat(f64, value) catch {
                    std.debug.print("Error: Invalid probe radius: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "--n-points=")) {
                const value = arg["--n-points=".len..];
                result.n_points = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid n-points: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "--n-slices=")) {
                const value = arg["--n-slices=".len..];
                result.n_slices = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid n-slices: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "--precision=")) {
                const value = arg["--precision=".len..];
                result.precision = parsePrecision(value);
            } else if (std.mem.startsWith(u8, arg, "--classifier=")) {
                const value = arg["--classifier=".len..];
                result.classifier_type = parseClassifierType(value);
            } else if (std.mem.startsWith(u8, arg, "--ccd=")) {
                const value = arg["--ccd=".len..];
                result.ccd_path = value;
            } else if (std.mem.eql(u8, arg, "--ccd")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --ccd\n", .{});
                    std.process.exit(1);
                }
                result.ccd_path = args[i];
            } else if (std.mem.startsWith(u8, arg, "--stride=")) {
                const value = arg["--stride=".len..];
                result.stride = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid stride: {s}\n", .{value});
                    std.process.exit(1);
                };
                if (result.stride == 0) result.stride = 1;
            } else if (std.mem.startsWith(u8, arg, "--start=")) {
                const value = arg["--start=".len..];
                result.start_frame = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid start frame: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "--end=")) {
                const value = arg["--end=".len..];
                result.end_frame = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid end frame: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "--batch-size=")) {
                const value = arg["--batch-size=".len..];
                const bs = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid batch size: {s}\n", .{value});
                    std.process.exit(1);
                };
                if (bs == 0) {
                    std.debug.print("Error: Batch size must be >= 1 (omit for auto)\n", .{});
                    std.process.exit(1);
                }
                result.batch_size = bs;
            } else if (std.mem.startsWith(u8, arg, "-o=") or std.mem.startsWith(u8, arg, "--output=")) {
                const prefix_len = if (std.mem.startsWith(u8, arg, "-o=")) "-o=".len else "--output=".len;
                result.output_path = arg[prefix_len..];
            } else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
                result.include_hydrogens = true;
            } else if (std.mem.eql(u8, arg, "--no-hydrogens") or std.mem.eql(u8, arg, "--exclude-hydrogens")) {
                result.include_hydrogens = false;
            } else if (std.mem.eql(u8, arg, "--use-bitmask")) {
                result.use_bitmask = true;
            } else if (std.mem.eql(u8, arg, "-q") or std.mem.eql(u8, arg, "--quiet")) {
                result.quiet = true;
            } else {
                std.debug.print("Error: Unknown option: {s}\n", .{arg});
                std.process.exit(1);
            }
        } else if (std.mem.eql(u8, arg, "-h")) {
            result.show_help = true;
        } else if (std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
        } else if (std.mem.startsWith(u8, arg, "-o")) {
            // -o FILE
            if (arg.len > 2) {
                result.output_path = arg[2..];
            } else {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for -o\n", .{});
                    std.process.exit(1);
                }
                result.output_path = args[i];
            }
        } else {
            // Positional arguments
            if (positional_count == 0) {
                result.traj_path = arg;
            } else if (positional_count == 1) {
                result.topology_path = arg;
            } else {
                std.debug.print("Error: Too many positional arguments\n", .{});
                std.process.exit(1);
            }
            positional_count += 1;
        }
    }

    return result;
}

fn parseAlgorithm(value: []const u8) Algorithm {
    if (std.mem.eql(u8, value, "sr") or std.mem.eql(u8, value, "shrake-rupley")) {
        return .sr;
    } else if (std.mem.eql(u8, value, "lr") or std.mem.eql(u8, value, "lee-richards")) {
        return .lr;
    } else {
        std.debug.print("Error: Invalid algorithm: {s}\n", .{value});
        std.process.exit(1);
    }
}

fn parsePrecision(value: []const u8) Precision {
    if (std.mem.eql(u8, value, "f32")) {
        return .f32;
    } else if (std.mem.eql(u8, value, "f64")) {
        return .f64;
    } else {
        std.debug.print("Error: Invalid precision: {s}\n", .{value});
        std.process.exit(1);
    }
}

fn parseClassifierType(value: []const u8) ClassifierType {
    if (ClassifierType.fromString(value)) |ct| {
        return ct;
    } else {
        std.debug.print("Error: Invalid classifier: {s}\n", .{value});
        std.debug.print("Valid classifiers: ccd, protor, naccess, oons\n", .{});
        std.process.exit(1);
    }
}

/// Print help for trajectory mode
pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\Usage: {s} traj <trajectory> <topology> [options]
        \\
        \\Calculate SASA for each frame in a trajectory.
        \\Supported formats: XTC (GROMACS), DCD (NAMD/CHARMM).
        \\Format is auto-detected from file extension.
        \\
        \\ARGUMENTS:
        \\    <trajectory> Trajectory file (.xtc or .dcd)
        \\    <topology>   Topology file (PDB or mmCIF) for atom names and radii
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: ccd, protor, naccess, oons
        \\                       Default: naccess (supports explicit H in MD trajectories)
        \\    --ccd=PATH         External CCD dictionary file (.zsdc or .cif[.gz])
        \\                       Used with --classifier=ccd for non-standard residues
        \\    --threads=N        Number of threads (default: auto-detect)
        \\    --probe-radius=R   Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N       Test points per atom (default: 100, for sr)
        \\    --n-slices=N       Slices per atom diameter (default: 20, for lr)
        \\    --precision=PREC    Floating-point precision: f32, f64 (default: f32)
        \\    --no-hydrogens     Exclude hydrogen atoms (default: included)
        \\    --include-hydrogens
        \\                       Include hydrogen atoms (default, for backward compat)
        \\    --use-bitmask      Use bitmask LUT optimization for SR algorithm
        \\                       (n-points must be 1..1024)
        \\    --stride=N         Process every Nth frame (default: 1)
        \\    --start=N          Start from frame N (default: 0)
        \\    --end=N            End at frame N (default: all)
        \\    --batch-size=N     Frames per batch for parallel processing
        \\                       Default: auto (threads * 2)
        \\    -o, --output=FILE  Output CSV file (default: traj_sasa.csv)
        \\    -q, --quiet        Suppress progress output
        \\    -h, --help         Show this help message
        \\
        \\OUTPUT FORMAT (CSV):
        \\    frame,step,time,total_sasa
        \\    0,1,0.0,12345.67
        \\    1,2,1.0,12340.12
        \\    ...
        \\
        \\EXAMPLES:
        \\    {s} traj trajectory.xtc topology.pdb
        \\    {s} traj trajectory.dcd topology.pdb
        \\    {s} traj trajectory.xtc topology.pdb -o sasa.csv
        \\    {s} traj trajectory.xtc topology.pdb --stride=10
        \\    {s} traj trajectory.xtc topology.pdb --classifier=naccess
        \\    {s} traj trajectory.xtc topology.pdb --algorithm=lr --n-slices=50
        \\    {s} traj trajectory.xtc topology.pdb --threads=8
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
}

/// Detect topology file format
fn detectTopologyFormat(path: []const u8) enum { pdb, mmcif } {
    if (std.mem.endsWith(u8, path, ".cif") or std.mem.endsWith(u8, path, ".mmcif")) {
        return .mmcif;
    }
    return .pdb;
}

// =============================================================================
// Batch Processing Infrastructure
// =============================================================================

/// Frame metadata for batch processing
const FrameData = struct {
    frame_idx: u32,
    step: i32,
    time: f32,
};

/// SASA result for one frame
const FrameResult = struct {
    frame_idx: u32 = 0,
    step: i32 = 0,
    time: f32 = 0,
    total_sasa: f64 = 0,
};

/// Helper to build and hold bitmask LUTs for trajectory processing.
/// Builds the appropriate LUT once based on args, and provides typed pointers.
const TrajLuts = struct {
    lut_f64: ?bitmask_lut.BitmaskLut = null,
    lut_f32: ?bitmask_lut.BitmaskLutGen(f32) = null,

    fn init(allocator: Allocator, args: TrajArgs) !TrajLuts {
        if (!args.use_bitmask) return .{};
        if (args.algorithm != .sr) return error.BitmaskRequiresSR;
        return switch (args.precision) {
            .f64 => .{ .lut_f64 = try bitmask_lut.BitmaskLut.init(allocator, args.n_points) },
            .f32 => .{ .lut_f32 = try bitmask_lut.BitmaskLutGen(f32).init(allocator, args.n_points) },
        };
    }

    fn deinit(self: *TrajLuts) void {
        if (self.lut_f64) |*lut| lut.deinit();
        if (self.lut_f32) |*lut| lut.deinit();
        self.* = .{};
    }

    fn f64Ptr(self: *const TrajLuts) ?*const bitmask_lut.BitmaskLut {
        return if (self.lut_f64 != null) &self.lut_f64.? else null;
    }

    fn f32Ptr(self: *const TrajLuts) ?*const bitmask_lut.BitmaskLutGen(f32) {
        return if (self.lut_f32 != null) &self.lut_f32.? else null;
    }
};

/// Worker arguments for batch frame processing
const BatchWorkerArgs = struct {
    coord_pool: []const f32,
    frame_data: []const FrameData,
    results: []FrameResult,
    batch_count: usize,
    natoms: usize,
    radii: []f64,
    error_flag: *std.atomic.Value(bool),
    error_msg: *[128]u8,
    thread_id: usize,
    n_threads: usize,
    algorithm: Algorithm,
    precision: Precision,
    probe_radius: f64,
    n_points: u32,
    n_slices: u32,
    coord_scale: f64, // 10.0 for XTC (nm→Å), 1.0 for DCD (already Å)
    bitmask_lut_f64: ?*const bitmask_lut.BitmaskLut = null,
    bitmask_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32) = null,
};

/// Set error flag and store a descriptive message (first writer wins).
fn setWorkerError(args: BatchWorkerArgs, comptime fmt: []const u8, fmt_args: anytype) void {
    // Only the first error writes the message
    if (!args.error_flag.load(.acquire)) {
        _ = std.fmt.bufPrint(args.error_msg, fmt, fmt_args) catch {};
    }
    args.error_flag.store(true, .release);
}

/// Worker function for batch frame processing.
/// Each thread processes frames at indices: thread_id, thread_id + n_threads, ...
/// Allocates coordinate buffers per frame; arena reset retains backing memory.
fn batchWorkerFn(args: BatchWorkerArgs) void {
    // Use smp_allocator as arena backing to avoid mmap/munmap syscall contention
    // that page_allocator causes under multi-threaded workloads.
    var arena = std.heap.ArenaAllocator.init(std.heap.smp_allocator);
    defer arena.deinit();
    const thread_alloc = arena.allocator();

    // Process assigned frames (stride distribution across threads)
    var batch_idx = args.thread_id;
    while (batch_idx < args.batch_count) : (batch_idx += args.n_threads) {
        if (args.error_flag.load(.acquire)) return;

        // Allocate coordinate buffers per frame (free after arena reset;
        // re-alloc is ~free with retain_capacity since the backing memory persists)
        const x = thread_alloc.alloc(f64, args.natoms) catch {
            setWorkerError(args, "thread {d}: failed to allocate x buffer", .{args.thread_id});
            return;
        };
        const y = thread_alloc.alloc(f64, args.natoms) catch {
            setWorkerError(args, "thread {d}: failed to allocate y buffer", .{args.thread_id});
            return;
        };
        const z = thread_alloc.alloc(f64, args.natoms) catch {
            setWorkerError(args, "thread {d}: failed to allocate z buffer", .{args.thread_id});
            return;
        };

        const coord_offset = batch_idx * args.natoms * 3;

        // Convert f32 coordinates to f64 with unit scaling (nm→Å for XTC, no-op for DCD)
        const scale = args.coord_scale;
        for (0..args.natoms) |i| {
            x[i] = @as(f64, args.coord_pool[coord_offset + i * 3 + 0]) * scale;
            y[i] = @as(f64, args.coord_pool[coord_offset + i * 3 + 1]) * scale;
            z[i] = @as(f64, args.coord_pool[coord_offset + i * 3 + 2]) * scale;
        }

        const input = AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = args.radii,
            .allocator = thread_alloc,
        };

        // Calculate SASA (single-threaded per frame)
        var total_sasa: f64 = 0;
        const frame_id = args.frame_data[batch_idx].frame_idx;

        if (args.precision == .f32) {
            const config = Configf32{
                .probe_radius = @floatCast(args.probe_radius),
                .n_points = args.n_points,
            };
            if (args.bitmask_lut_f32) |lut| {
                var result = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithLut(thread_alloc, input, config, lut) catch |err| {
                    setWorkerError(args, "frame {d}: bitmask-f32 failed: {s}", .{ frame_id, @errorName(err) });
                    return;
                };
                total_sasa = @floatCast(result.total_area);
                result.deinit();
            } else {
                var result = shrake_rupley.calculateSasaf32(thread_alloc, input, config) catch |err| {
                    setWorkerError(args, "frame {d}: SR-f32 SASA failed: {s}", .{ frame_id, @errorName(err) });
                    return;
                };
                total_sasa = @floatCast(result.total_area);
                result.deinit();
            }
        } else {
            if (args.algorithm == .sr) {
                const config = Config{
                    .probe_radius = args.probe_radius,
                    .n_points = args.n_points,
                };
                if (args.bitmask_lut_f64) |lut| {
                    var result = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithLut(thread_alloc, input, config, lut) catch |err| {
                        setWorkerError(args, "frame {d}: bitmask-f64 failed: {s}", .{ frame_id, @errorName(err) });
                        return;
                    };
                    total_sasa = result.total_area;
                    result.deinit();
                } else {
                    var result = shrake_rupley.calculateSasa(thread_alloc, input, config) catch |err| {
                        setWorkerError(args, "frame {d}: SR-f64 SASA failed: {s}", .{ frame_id, @errorName(err) });
                        return;
                    };
                    total_sasa = result.total_area;
                    result.deinit();
                }
            } else {
                const config = lee_richards.LeeRichardsConfig{
                    .probe_radius = args.probe_radius,
                    .n_slices = args.n_slices,
                };
                var result = lee_richards.calculateSasa(thread_alloc, input, config) catch |err| {
                    setWorkerError(args, "frame {d}: LR SASA failed: {s}", .{ frame_id, @errorName(err) });
                    return;
                };
                total_sasa = result.total_area;
                result.deinit();
            }
        }

        args.results[batch_idx] = .{
            .frame_idx = frame_id,
            .step = args.frame_data[batch_idx].step,
            .time = args.frame_data[batch_idx].time,
            .total_sasa = total_sasa,
        };

        // Reset arena for next frame (retain backing memory to avoid syscalls)
        _ = arena.reset(.retain_capacity);
    }
}

/// Read frames into batch buffer, applying stride/start/end filtering.
/// Returns the number of frames read and whether EOF was reached.
fn readBatch(
    reader: *TrajectoryReader,
    allocator: Allocator,
    coord_pool: []f32,
    frame_data: []FrameData,
    batch_size: usize,
    natoms: usize,
    frame_idx: *u32,
    traj_args: TrajArgs,
) !struct { count: usize, eof: bool } {
    var batch_count: usize = 0;

    while (batch_count < batch_size) {
        var frame = reader.readFrame() catch |err| {
            if (reader.isEof(err)) return .{ .count = batch_count, .eof = true };
            return err;
        };
        defer frame.deinit(allocator);

        // Apply frame range filter
        if (frame_idx.* < traj_args.start_frame) {
            frame_idx.* += 1;
            continue;
        }
        if (traj_args.end_frame) |end| {
            if (frame_idx.* > end) return .{ .count = batch_count, .eof = true };
        }

        // Apply stride filter
        if ((frame_idx.* - traj_args.start_frame) % traj_args.stride != 0) {
            frame_idx.* += 1;
            continue;
        }

        // Copy coordinates to pool
        const offset = batch_count * natoms * 3;
        @memcpy(coord_pool[offset .. offset + natoms * 3], frame.coords[0 .. natoms * 3]);

        frame_data[batch_count] = .{
            .frame_idx = frame_idx.*,
            .step = frame.step,
            .time = frame.time,
        };
        batch_count += 1;
        frame_idx.* += 1;
    }

    return .{ .count = batch_count, .eof = false };
}

// =============================================================================
// Main Entry Point
// =============================================================================

/// Run trajectory analysis
pub fn run(allocator: Allocator, args: TrajArgs) !void {
    // Validate required arguments
    const traj_path = args.traj_path orelse {
        std.debug.print("Error: Missing trajectory file\n", .{});
        return error.MissingArgument;
    };
    const topology_path = args.topology_path orelse {
        std.debug.print("Error: Missing topology file\n", .{});
        return error.MissingArgument;
    };

    // Detect trajectory format
    const traj_format = detectTrajectoryFormat(traj_path) orelse {
        std.debug.print("Error: Unknown trajectory format. Supported: .xtc, .dcd\n", .{});
        return error.UnsupportedFormat;
    };

    // CCD/ProtOr use united-atom radii (implicit H) — warn if explicit H included
    if ((args.classifier_type == .ccd or args.classifier_type == .protor) and args.include_hydrogens and !args.quiet) {
        std.debug.print("Warning: --include-hydrogens with CCD classifier may give inaccurate results\n", .{});
        std.debug.print("         CCD uses united-atom radii that already account for implicit hydrogens\n", .{});
    }

    // Read topology to get atom names and radii
    if (!args.quiet) {
        std.debug.print("Reading topology: {s}\n", .{topology_path});
    }

    var topology = switch (detectTopologyFormat(topology_path)) {
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.skip_hydrogens = !args.include_hydrogens;
            break :blk try parser.parseFile(topology_path);
        },
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.skip_hydrogens = !args.include_hydrogens;
            break :blk try parser.parseFile(topology_path);
        },
    };
    defer topology.deinit();

    const natoms = topology.atomCount();
    if (!args.quiet) {
        std.debug.print("Topology: {d} atoms\n", .{natoms});
    }

    // Load external CCD dictionary if specified
    var ext_ccd: ?ccd_parser.ComponentDict = null;
    if (args.ccd_path) |ccd_path| {
        const ccd_data = if (std.mem.endsWith(u8, ccd_path, ".gz"))
            try gzip.readGzip(allocator, ccd_path)
        else blk: {
            const f = try std.fs.cwd().openFile(ccd_path, .{});
            defer f.close();
            break :blk try f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024);
        };
        defer allocator.free(ccd_data);

        ext_ccd = try ccd_binary.loadDict(allocator, ccd_data);
        if (!args.quiet) {
            std.debug.print("External CCD: loaded {d} components from '{s}'\n", .{ ext_ccd.?.components.count(), ccd_path });
        }
    }
    defer if (ext_ccd) |*d| d.deinit();

    // Apply classifier if specified
    if (args.classifier_type) |ct| {
        const ext_ccd_ptr: ?*const ccd_parser.ComponentDict = if (ext_ccd != null) &ext_ccd.? else null;
        try applyBuiltinClassifier(allocator, &topology, ct, ext_ccd_ptr, args.quiet);
    }

    // Open trajectory file
    if (!args.quiet) {
        const format_name: []const u8 = switch (traj_format) {
            .xtc => "XTC",
            .dcd => "DCD",
        };
        std.debug.print("Opening trajectory ({s}): {s}\n", .{ format_name, traj_path });
    }

    // Open XTC or DCD reader and verify atom count
    var xtc_reader: ?xtc.XtcReader = null;
    var dcd_reader: ?dcd.DcdReader = null;
    const traj_natoms: i32 = switch (traj_format) {
        .xtc => blk: {
            xtc_reader = try xtc.XtcReader.open(allocator, traj_path);
            break :blk xtc_reader.?.natoms;
        },
        .dcd => blk: {
            dcd_reader = try dcd.DcdReader.open(allocator, traj_path);
            break :blk dcd_reader.?.natoms;
        },
    };
    defer {
        if (xtc_reader) |*r| r.close();
        if (dcd_reader) |*r| r.close();
    }

    // Verify atom count matches
    if (traj_natoms != @as(i32, @intCast(natoms))) {
        std.debug.print("Error: Atom count mismatch - trajectory has {d} atoms, topology has {d}\n", .{
            traj_natoms,
            natoms,
        });
        return error.AtomCountMismatch;
    }

    // Open output file with buffered writer
    const output_file = try std.fs.cwd().createFile(args.output_path, .{});
    defer output_file.close();
    var write_buffer: [4096]u8 = undefined;
    var buffered_writer = output_file.writer(&write_buffer);
    const writer = &buffered_writer.interface;

    // Write CSV header
    try writer.writeAll("frame,step,time,total_sasa\n");

    // Resolve thread count
    const n_threads = if (args.n_threads == 0)
        @as(usize, @intCast(std.Thread.getCpuCount() catch 1))
    else
        args.n_threads;

    // Build bitmask LUT once (if enabled) — reused across all frames
    var luts = TrajLuts.init(allocator, args) catch |err| {
        switch (err) {
            error.BitmaskRequiresSR => {
                std.debug.print("Error: --use-bitmask requires --algorithm=sr\n", .{});
                return error.BitmaskRequiresSR;
            },
            error.UnsupportedNPoints => {
                std.debug.print("Error: --use-bitmask requires --n-points=1..1024\n", .{});
                return error.UnsupportedNPoints;
            },
            else => {
                std.debug.print("Error: failed to build bitmask LUT: {s}\n", .{@errorName(err)});
                return err;
            },
        }
    };
    defer luts.deinit();

    if (!args.quiet) {
        std.debug.print("Algorithm: {s}, Threads: {d}, Precision: {s}{s}\n", .{
            if (args.algorithm == .sr) "Shrake-Rupley" else "Lee-Richards",
            n_threads,
            if (args.precision == .f32) "f32" else "f64",
            if (args.use_bitmask) ", Bitmask: enabled" else "",
        });
        std.debug.print("Processing frames", .{});
        if (args.stride > 1) std.debug.print(" (stride={d})", .{args.stride});
        if (args.start_frame > 0) std.debug.print(" from {d}", .{args.start_frame});
        if (args.end_frame) |end| std.debug.print(" to {d}", .{end});
        std.debug.print("...\n\n", .{});
    }

    // Process frames
    var frame_idx: u32 = 0;
    var processed_count: u32 = 0;
    const timer_start = std.time.milliTimestamp();

    var traj_reader = TrajectoryReader{
        .xtc_reader = if (xtc_reader) |*r| r else null,
        .dcd_reader = if (dcd_reader) |*r| r else null,
        .format = traj_format,
    };

    if (n_threads <= 1) {
        // Sequential path: single-threaded SASA per frame
        try runSequential(allocator, &traj_reader, writer, &topology, args, &luts, &frame_idx, &processed_count);
    } else {
        // Batch parallel path: frame-level parallelism across threads
        try runBatchParallel(allocator, &traj_reader, writer, &topology, args, &luts, n_threads, natoms, &frame_idx, &processed_count);
    }

    // Flush buffered output
    try writer.flush();

    const elapsed = std.time.milliTimestamp() - timer_start;

    if (!args.quiet) {
        std.debug.print("\r  Processed {d} frames in {d}ms ({d:.1} frames/sec)\n", .{
            processed_count,
            elapsed,
            if (elapsed > 0) @as(f64, @floatFromInt(processed_count)) * 1000.0 / @as(f64, @floatFromInt(elapsed)) else 0,
        });
        std.debug.print("Output written to: {s}\n", .{args.output_path});
    }
}

/// Sequential frame processing (single-threaded SASA per frame)
fn runSequential(
    allocator: Allocator,
    reader: *TrajectoryReader,
    writer: anytype,
    topology: *AtomInput,
    args: TrajArgs,
    luts: *const TrajLuts,
    frame_idx: *u32,
    processed_count: *u32,
) !void {
    const natoms = topology.atomCount();

    // Allocate mutable coordinate buffers
    const frame_x = try allocator.alloc(f64, natoms);
    defer allocator.free(frame_x);
    const frame_y = try allocator.alloc(f64, natoms);
    defer allocator.free(frame_y);
    const frame_z = try allocator.alloc(f64, natoms);
    defer allocator.free(frame_z);

    const frame_input = AtomInput{
        .x = frame_x,
        .y = frame_y,
        .z = frame_z,
        .r = topology.r,
        .residue = topology.residue,
        .atom_name = topology.atom_name,
        .element = topology.element,
        .chain_id = topology.chain_id,
        .residue_num = topology.residue_num,
        .insertion_code = topology.insertion_code,
        .allocator = allocator,
    };

    const scale = reader.coordScale();

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (reader.isEof(err)) break;
            return err;
        };
        defer frame.deinit(allocator);

        // Check frame range
        if (frame_idx.* < args.start_frame) {
            frame_idx.* += 1;
            continue;
        }
        if (args.end_frame) |end| {
            if (frame_idx.* > end) break;
        }

        // Check stride
        if ((frame_idx.* - args.start_frame) % args.stride != 0) {
            frame_idx.* += 1;
            continue;
        }

        // Update frame coordinates with unit scaling (nm→Å for XTC, no-op for DCD)
        for (0..natoms) |i| {
            frame_x[i] = @as(f64, frame.coords[i * 3 + 0]) * scale;
            frame_y[i] = @as(f64, frame.coords[i * 3 + 1]) * scale;
            frame_z[i] = @as(f64, frame.coords[i * 3 + 2]) * scale;
        }

        // Calculate SASA (single-threaded)
        const total_sasa: f64 = switch (args.precision) {
            .f32 => blk: {
                const config = Configf32{
                    .probe_radius = @floatCast(args.probe_radius),
                    .n_points = args.n_points,
                };
                if (luts.f32Ptr()) |lut| {
                    var result = try shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithLut(allocator, frame_input, config, lut);
                    defer result.deinit();
                    break :blk @floatCast(result.total_area);
                } else {
                    var result = try shrake_rupley.calculateSasaf32(allocator, frame_input, config);
                    defer result.deinit();
                    break :blk @floatCast(result.total_area);
                }
            },
            .f64 => blk: {
                switch (args.algorithm) {
                    .sr => {
                        const config = Config{
                            .probe_radius = args.probe_radius,
                            .n_points = args.n_points,
                        };
                        if (luts.f64Ptr()) |lut| {
                            var result = try shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithLut(allocator, frame_input, config, lut);
                            defer result.deinit();
                            break :blk result.total_area;
                        } else {
                            var result = try shrake_rupley.calculateSasa(allocator, frame_input, config);
                            defer result.deinit();
                            break :blk result.total_area;
                        }
                    },
                    .lr => {
                        const config = lee_richards.LeeRichardsConfig{
                            .probe_radius = args.probe_radius,
                            .n_slices = args.n_slices,
                        };
                        var result = try lee_richards.calculateSasa(allocator, frame_input, config);
                        defer result.deinit();
                        break :blk result.total_area;
                    },
                }
            },
        };

        // Write result
        try writer.print("{d},{d},{d:.3},{d:.2}\n", .{
            frame_idx.*,
            frame.step,
            frame.time,
            total_sasa,
        });

        processed_count.* += 1;
        frame_idx.* += 1;

        // Progress indicator
        if (!args.quiet and processed_count.* % 100 == 0) {
            std.debug.print("\r  Processed {d} frames...", .{processed_count.*});
        }
    }
}

/// Batch parallel frame processing (frame-level parallelism)
fn runBatchParallel(
    allocator: Allocator,
    reader: *TrajectoryReader,
    writer: anytype,
    topology: *AtomInput,
    args: TrajArgs,
    luts: *const TrajLuts,
    n_threads: usize,
    natoms: usize,
    frame_idx: *u32,
    processed_count: *u32,
) !void {
    const batch_size: usize = if (args.batch_size > 0) args.batch_size else n_threads * 2;

    if (!args.quiet) {
        std.debug.print("Batch size: {d} frames\n", .{batch_size});
    }

    // Pre-allocate batch buffers
    const coord_pool = try allocator.alloc(f32, batch_size * natoms * 3);
    defer allocator.free(coord_pool);
    const frame_data = try allocator.alloc(FrameData, batch_size);
    defer allocator.free(frame_data);
    const frame_results = try allocator.alloc(FrameResult, batch_size);
    defer allocator.free(frame_results);

    while (true) {
        // Phase 1: Read batch of frames (sequential, with filtering)
        const read_result = try readBatch(
            reader,
            allocator,
            coord_pool,
            frame_data,
            batch_size,
            natoms,
            frame_idx,
            args,
        );

        if (read_result.count == 0) break;

        // Phase 2: Compute batch (parallel, frame-level distribution)
        const thread_count = @min(n_threads, read_result.count);
        var error_flag = std.atomic.Value(bool).init(false);
        var error_msg: [128]u8 = @splat(0);

        const threads = try allocator.alloc(std.Thread, thread_count);
        defer allocator.free(threads);

        for (0..thread_count) |i| {
            threads[i] = std.Thread.spawn(.{}, batchWorkerFn, .{BatchWorkerArgs{
                .coord_pool = coord_pool,
                .frame_data = frame_data,
                .results = frame_results,
                .batch_count = read_result.count,
                .natoms = natoms,
                .radii = topology.r,
                .error_flag = &error_flag,
                .error_msg = &error_msg,
                .thread_id = i,
                .n_threads = thread_count,
                .algorithm = args.algorithm,
                .precision = args.precision,
                .probe_radius = args.probe_radius,
                .n_points = args.n_points,
                .n_slices = args.n_slices,
                .coord_scale = reader.coordScale(),
                .bitmask_lut_f64 = luts.f64Ptr(),
                .bitmask_lut_f32 = luts.f32Ptr(),
            }}) catch {
                error_flag.store(true, .release);
                for (threads[0..i]) |thread| {
                    thread.join();
                }
                return error.ThreadSpawnFailed;
            };
        }

        // Wait for all threads to complete
        for (threads[0..thread_count]) |thread| {
            thread.join();
        }

        if (error_flag.load(.acquire)) {
            const msg = std.mem.sliceTo(&error_msg, 0);
            if (msg.len > 0) {
                std.debug.print("Error: {s}\n", .{msg});
            }
            return error.BatchCalculationFailed;
        }

        // Phase 3: Write results (sequential, in frame order)
        for (0..read_result.count) |i| {
            try writer.print("{d},{d},{d:.3},{d:.2}\n", .{
                frame_results[i].frame_idx,
                frame_results[i].step,
                frame_results[i].time,
                frame_results[i].total_sasa,
            });
        }

        processed_count.* += @intCast(read_result.count);

        if (!args.quiet and processed_count.* % 100 < @as(u32, @intCast(read_result.count))) {
            std.debug.print("\r  Processed {d} frames...", .{processed_count.*});
        }

        if (read_result.eof) break;
    }
}

/// Apply built-in classifier to topology
fn applyBuiltinClassifier(
    _: Allocator,
    input: *AtomInput,
    ct: ClassifierType,
    external_ccd: ?*const ccd_parser.ComponentDict,
    quiet: bool,
) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    // For CCD/ProtOr: create classifier instance and feed external CCD components
    var ccd_clf: ?classifier_ccd.CcdClassifier = if (ct == .ccd or ct == .protor) classifier_ccd.CcdClassifier.init(input.allocator) else null;
    defer if (ccd_clf) |*c| c.deinit();

    if (ccd_clf != null) {
        if (external_ccd) |dict| {
            // Deduplicate: collect unique non-hardcoded residues
            var needed = std.StringHashMap(void).init(input.allocator);
            defer needed.deinit();
            for (0..n) |i| {
                const res = residues[i].slice();
                if (!classifier_ccd.CcdClassifier.isHardcoded(res)) {
                    try needed.put(res, {});
                }
            }
            var loaded: usize = 0;
            var it = needed.keyIterator();
            while (it.next()) |key_ptr| {
                if (dict.get(key_ptr.*)) |comp| {
                    ccd_clf.?.addComponent(&comp) catch |err| {
                        if (!quiet) {
                            std.debug.print("Warning: Could not derive CCD properties for '{s}': {s}\n", .{ key_ptr.*, @errorName(err) });
                        }
                        continue;
                    };
                    loaded += 1;
                }
            }
            if (!quiet and loaded > 0) {
                std.debug.print("CCD: {d} non-standard components derived from external CCD\n", .{loaded});
            }
        }
    }

    var classified_count: usize = 0;
    var fallback_count: usize = 0;

    for (0..n) |i| {
        const radius_opt: ?f64 = switch (ct) {
            .naccess => classifier_naccess.getRadius(residues[i].slice(), atom_names[i].slice()),
            .protor, .ccd => if (ccd_clf) |*c| c.getRadius(residues[i].slice(), atom_names[i].slice()) else null,
            .oons => classifier_oons.getRadius(residues[i].slice(), atom_names[i].slice()),
        };
        if (radius_opt) |r| {
            input.r[i] = r;
            classified_count += 1;
        } else if (input.element) |elements| {
            if (classifier.guessRadiusFromAtomicNumber(elements[i])) |r| {
                input.r[i] = r;
                fallback_count += 1;
            }
        } else if (classifier.guessRadiusFromAtomName(atom_names[i].slice())) |r| {
            input.r[i] = r;
            fallback_count += 1;
        }
    }

    if (!quiet) {
        std.debug.print("Classifier '{s}': {d} atoms classified, {d} fallback\n", .{
            ct.name(),
            classified_count,
            fallback_count,
        });
    }
}

// Trajectory analysis mode
// Calculates SASA for each frame in an XTC trajectory
//
// Supports two parallelism strategies:
//   - Sequential (1 thread): processes frames one at a time
//   - Batch parallel (N threads): reads batches of frames, distributes across threads
//
const std = @import("std");
const xtc = @import("xtc.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const pdb_parser = @import("pdb_parser.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const classifier = @import("classifier.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_protor = @import("classifier_protor.zig");
const classifier_oons = @import("classifier_oons.zig");

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

/// Trajectory mode arguments
pub const TrajArgs = struct {
    xtc_path: ?[]const u8 = null,
    topology_path: ?[]const u8 = null,
    output_path: []const u8 = "traj_sasa.csv",
    algorithm: Algorithm = .sr,
    n_threads: usize = 0,
    probe_radius: f64 = 1.4,
    n_points: u32 = 100,
    n_slices: u32 = 20,
    precision: Precision = .f32, // Default f32 for trajectory (speed)
    classifier_type: ?ClassifierType = null,
    stride: u32 = 1, // Process every Nth frame
    start_frame: u32 = 0, // Start frame
    end_frame: ?u32 = null, // End frame (null = all)
    include_hydrogens: bool = false, // Include hydrogen atoms (default: exclude)
    batch_size: u32 = 0, // Frames per batch for parallel processing (0 = auto)
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
                result.batch_size = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid batch size: {s}\n", .{value});
                    std.process.exit(1);
                };
            } else if (std.mem.startsWith(u8, arg, "-o=") or std.mem.startsWith(u8, arg, "--output=")) {
                const prefix_len = if (std.mem.startsWith(u8, arg, "-o=")) "-o=".len else "--output=".len;
                result.output_path = arg[prefix_len..];
            } else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
                result.include_hydrogens = true;
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
                result.xtc_path = arg;
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
    if (std.mem.eql(u8, value, "naccess")) {
        return .naccess;
    } else if (std.mem.eql(u8, value, "protor")) {
        return .protor;
    } else if (std.mem.eql(u8, value, "oons")) {
        return .oons;
    } else {
        std.debug.print("Error: Invalid classifier: {s}\n", .{value});
        std.process.exit(1);
    }
}

/// Print help for trajectory mode
pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\Usage: {s} traj <xtc> <topology> [options]
        \\
        \\Calculate SASA for each frame in an XTC trajectory.
        \\
        \\ARGUMENTS:
        \\    <xtc>        XTC trajectory file
        \\    <topology>   Topology file (PDB or mmCIF) for atom names and radii
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: naccess, protor, oons
        \\    --threads=N        Number of threads (default: auto-detect)
        \\    --probe-radius=R   Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N       Test points per atom (default: 100, for sr)
        \\    --n-slices=N       Slices per atom diameter (default: 20, for lr)
        \\    --precision=PREC    Floating-point precision: f32, f64 (default: f32)
        \\    --include-hydrogens
        \\                       Include hydrogen atoms (default: excluded)
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
        \\    frame,time,total_sasa
        \\    0,0.0,12345.67
        \\    1,1.0,12340.12
        \\    ...
        \\
        \\EXAMPLES:
        \\    {s} traj trajectory.xtc topology.pdb
        \\    {s} traj trajectory.xtc topology.pdb -o sasa.csv
        \\    {s} traj trajectory.xtc topology.pdb --stride=10
        \\    {s} traj trajectory.xtc topology.pdb --classifier=naccess
        \\    {s} traj trajectory.xtc topology.pdb --algorithm=lr --n-slices=50
        \\    {s} traj trajectory.xtc topology.pdb --threads=8
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name, program_name });
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

/// Worker arguments for batch frame processing
const BatchWorkerArgs = struct {
    coord_pool: []const f32,
    frame_data: []const FrameData,
    results: []FrameResult,
    batch_count: usize,
    natoms: usize,
    radii: []f64,
    error_flag: *std.atomic.Value(bool),
    thread_id: usize,
    n_threads: usize,
    algorithm: Algorithm,
    precision: Precision,
    probe_radius: f64,
    n_points: u32,
    n_slices: u32,
};

/// Worker function for batch frame processing.
/// Each thread processes frames at indices: thread_id, thread_id + n_threads, ...
/// Allocates per-thread coordinate buffers and reuses them across frames.
fn batchWorkerFn(args: BatchWorkerArgs) void {
    const thread_alloc = std.heap.page_allocator;

    // Pre-allocate coordinate buffers (reused across all frames in this thread)
    const x = thread_alloc.alloc(f64, args.natoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer thread_alloc.free(x);

    const y = thread_alloc.alloc(f64, args.natoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer thread_alloc.free(y);

    const z = thread_alloc.alloc(f64, args.natoms) catch {
        args.error_flag.store(true, .release);
        return;
    };
    defer thread_alloc.free(z);

    // Process assigned frames (stride distribution across threads)
    var batch_idx = args.thread_id;
    while (batch_idx < args.batch_count) : (batch_idx += args.n_threads) {
        if (args.error_flag.load(.acquire)) return;

        const coord_offset = batch_idx * args.natoms * 3;

        // Convert f32 XTC coordinates (nm) to f64 (Angstrom)
        for (0..args.natoms) |i| {
            x[i] = @as(f64, args.coord_pool[coord_offset + i * 3 + 0]) * 10.0;
            y[i] = @as(f64, args.coord_pool[coord_offset + i * 3 + 1]) * 10.0;
            z[i] = @as(f64, args.coord_pool[coord_offset + i * 3 + 2]) * 10.0;
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

        if (args.precision == .f32) {
            const config = Configf32{
                .probe_radius = @floatCast(args.probe_radius),
                .n_points = args.n_points,
            };
            var result = shrake_rupley.calculateSasaf32(thread_alloc, input, config) catch {
                args.error_flag.store(true, .release);
                return;
            };
            total_sasa = @floatCast(result.total_area);
            result.deinit();
        } else {
            if (args.algorithm == .sr) {
                const config = Config{
                    .probe_radius = args.probe_radius,
                    .n_points = args.n_points,
                };
                var result = shrake_rupley.calculateSasa(thread_alloc, input, config) catch {
                    args.error_flag.store(true, .release);
                    return;
                };
                total_sasa = result.total_area;
                result.deinit();
            } else {
                const config = lee_richards.LeeRichardsConfig{
                    .probe_radius = args.probe_radius,
                    .n_slices = args.n_slices,
                };
                var result = lee_richards.calculateSasa(thread_alloc, input, config) catch {
                    args.error_flag.store(true, .release);
                    return;
                };
                total_sasa = result.total_area;
                result.deinit();
            }
        }

        args.results[batch_idx] = .{
            .frame_idx = args.frame_data[batch_idx].frame_idx,
            .step = args.frame_data[batch_idx].step,
            .time = args.frame_data[batch_idx].time,
            .total_sasa = total_sasa,
        };
    }
}

/// Read frames into batch buffer, applying stride/start/end filtering.
/// Returns the number of frames read and whether EOF was reached.
fn readBatch(
    reader: *xtc.XtcReader,
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
            if (err == xtc.XtcError.EndOfFile) return .{ .count = batch_count, .eof = true };
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
    const xtc_path = args.xtc_path orelse {
        std.debug.print("Error: Missing XTC file\n", .{});
        return error.MissingArgument;
    };
    const topology_path = args.topology_path orelse {
        std.debug.print("Error: Missing topology file\n", .{});
        return error.MissingArgument;
    };

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

    // Apply classifier if specified
    if (args.classifier_type) |ct| {
        try applyBuiltinClassifier(allocator, &topology, ct, args.quiet);
    }

    // Open XTC file
    if (!args.quiet) {
        std.debug.print("Opening trajectory: {s}\n", .{xtc_path});
    }

    var reader = try xtc.XtcReader.open(allocator, xtc_path);
    defer reader.close();

    // Verify atom count matches
    if (reader.natoms != @as(i32, @intCast(natoms))) {
        std.debug.print("Error: Atom count mismatch - XTC has {d} atoms, topology has {d}\n", .{
            reader.natoms,
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

    if (!args.quiet) {
        std.debug.print("Algorithm: {s}, Threads: {d}, Precision: {s}\n", .{
            if (args.algorithm == .sr) "Shrake-Rupley" else "Lee-Richards",
            n_threads,
            if (args.precision == .f32) "f32" else "f64",
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

    if (n_threads <= 1) {
        // Sequential path: single-threaded SASA per frame
        try runSequential(allocator, &reader, writer, &topology, args, &frame_idx, &processed_count);
    } else {
        // Batch parallel path: frame-level parallelism across threads
        try runBatchParallel(allocator, &reader, writer, &topology, args, n_threads, natoms, &frame_idx, &processed_count);
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
    reader: *xtc.XtcReader,
    writer: anytype,
    topology: *AtomInput,
    args: TrajArgs,
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

    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == xtc.XtcError.EndOfFile) break;
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

        // Update frame coordinates from XTC (nm -> Angstrom)
        for (0..natoms) |i| {
            frame_x[i] = @as(f64, frame.coords[i * 3 + 0]) * 10.0;
            frame_y[i] = @as(f64, frame.coords[i * 3 + 1]) * 10.0;
            frame_z[i] = @as(f64, frame.coords[i * 3 + 2]) * 10.0;
        }

        // Calculate SASA (single-threaded)
        const total_sasa: f64 = switch (args.precision) {
            .f32 => blk: {
                const config = Configf32{
                    .probe_radius = @floatCast(args.probe_radius),
                    .n_points = args.n_points,
                };
                var result = try shrake_rupley.calculateSasaf32(allocator, frame_input, config);
                defer result.deinit();
                break :blk @floatCast(result.total_area);
            },
            .f64 => blk: {
                switch (args.algorithm) {
                    .sr => {
                        const config = Config{
                            .probe_radius = args.probe_radius,
                            .n_points = args.n_points,
                        };
                        var result = try shrake_rupley.calculateSasa(allocator, frame_input, config);
                        defer result.deinit();
                        break :blk result.total_area;
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
    reader: *xtc.XtcReader,
    writer: anytype,
    topology: *AtomInput,
    args: TrajArgs,
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
                .thread_id = i,
                .n_threads = thread_count,
                .algorithm = args.algorithm,
                .precision = args.precision,
                .probe_radius = args.probe_radius,
                .n_points = args.n_points,
                .n_slices = args.n_slices,
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
    quiet: bool,
) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    var classified_count: usize = 0;
    var fallback_count: usize = 0;

    for (0..n) |i| {
        const radius_opt: ?f64 = switch (ct) {
            .naccess => classifier_naccess.getRadius(residues[i].slice(), atom_names[i].slice()),
            .protor => classifier_protor.getRadius(residues[i].slice(), atom_names[i].slice()),
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

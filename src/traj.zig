// Trajectory analysis mode
// Calculates SASA for each frame in a molecular dynamics trajectory
//
// Supports two parallelism strategies:
//   - Sequential (1 thread): processes frames one at a time
//   - Batch parallel (N threads): reads batches of frames, distributes across threads
//
const std = @import("std");
const ztraj = @import("ztraj");
const xtc = ztraj.io.xtc;
const dcd = ztraj.io.dcd;
const trr = ztraj.io.trr;
const nc = ztraj.io.nc;
const shrake_rupley = @import("shrake_rupley.zig");
const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const lee_richards = @import("lee_richards.zig");
const calc = @import("calc.zig");
const types = @import("types.zig");
const pdb_parser = @import("pdb_parser.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const classifier = @import("classifier.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const sdf_parser = @import("sdf_parser.zig");
const compressed = @import("compressed.zig");

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
    xtc, // GROMACS XTC
    dcd, // NAMD/CHARMM DCD
    trr, // GROMACS TRR
    nc, // AMBER NetCDF
};

/// Detect trajectory format from file extension
pub fn detectTrajectoryFormat(path: []const u8) ?TrajectoryFormat {
    if (std.mem.endsWith(u8, path, ".xtc")) {
        return .xtc;
    } else if (std.mem.endsWith(u8, path, ".dcd")) {
        return .dcd;
    } else if (std.mem.endsWith(u8, path, ".trr")) {
        return .trr;
    } else if (std.mem.endsWith(u8, path, ".nc") or std.mem.endsWith(u8, path, ".ncdf")) {
        return .nc;
    }
    return null;
}

/// Common frame type for trajectory readers
const TrajectoryFrame = struct {
    step: i32,
    time: f32,
    coords: []f32, // flat array of x,y,z coordinates (length = natoms * 3)

    fn deinit(self: *TrajectoryFrame, allocator: Allocator) void {
        allocator.free(self.coords);
    }
};

/// Unified trajectory reader wrapping ztraj trajectory readers
const TrajectoryReader = struct {
    xtc_reader: ?*xtc.XtcReader = null,
    dcd_reader: ?*dcd.DcdReader = null,
    trr_reader: ?*trr.TrrReader = null,
    nc_reader: ?*nc.NcReader = null,
    format: TrajectoryFormat,

    /// Coordinate scale factor. ztraj readers yield coordinates in Å.
    fn coordScale(self: *const TrajectoryReader) f64 {
        _ = self;
        return 1.0;
    }

    /// Read next frame, returning null-like via EndOfFile error
    fn readFrame(self: *TrajectoryReader, allocator: Allocator) !TrajectoryFrame {
        const frame = switch (self.format) {
            .xtc => try self.xtc_reader.?.next() orelse return error.EndOfFile,
            .dcd => try self.dcd_reader.?.next() orelse return error.EndOfFile,
            .trr => try self.trr_reader.?.next() orelse return error.EndOfFile,
            .nc => try self.nc_reader.?.next() orelse return error.EndOfFile,
        };

        const natoms = frame.x.len;
        const coords = try allocator.alloc(f32, natoms * 3);
        errdefer allocator.free(coords);

        for (0..natoms) |i| {
            coords[i * 3 + 0] = frame.x[i];
            coords[i * 3 + 1] = frame.y[i];
            coords[i * 3 + 2] = frame.z[i];
        }

        return .{
            .step = frame.step,
            .time = frame.time,
            .coords = coords,
        };
    }

    /// Check if a read error is an EOF condition
    fn isEof(self: *const TrajectoryReader, err: anyerror) bool {
        _ = self;
        return err == error.EndOfFile;
    }
};

const SdfPathList = sdf_parser.SdfPathList;

const bitmask_lut_cycle_count: usize = 4;

const BitmaskLutMode = enum {
    single,
    per_frame,
    cycle,
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
    ccd_path: ?[]const u8 = null, // External CCD dictionary file (.zsdc or .cif[.gz|.zst])
    sdf_paths: SdfPathList = .{}, // --sdf=PATH (up to 16)
    stride: u32 = 1, // Process every Nth frame
    start_frame: u32 = 0, // Start frame
    end_frame: ?u32 = null, // End frame (null = all)
    include_hydrogens: bool = true, // Include hydrogen atoms (default: include for MD trajectories)
    alt_loc_mode: mmcif_parser.AltLocMode = .auto, // Alternate-location handling for the topology
    alt_loc_id: u8 = 'A',
    batch_size: u32 = 0, // Frames per batch for parallel processing (0 = auto)
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 1..1024)
    bitmask_lut_mode: BitmaskLutMode = .single,
    bitmask_correction: bool = false, // Experimental exposed-fraction correction for bitmask SR
    bitmask_correction_coeff: f64 = shrake_rupley_bitmask.default_bitmask_correction_coeff,
    quiet: bool = false,
    show_progress: bool = true,
    show_help: bool = false,
};

fn shouldShowProgress(args: TrajArgs) bool {
    return args.show_progress and !args.quiet;
}

/// Error returned by the argument parser and the option range checks. The
/// message has already been printed when it is returned.
const ArgError = error{InvalidArgument};

fn parseBitmaskCorrectionCoeff(value: []const u8) ArgError!f64 {
    const coeff = std.fmt.parseFloat(f64, value) catch {
        std.debug.print("Error: Invalid bitmask correction coefficient: {s}\n", .{value});
        return error.InvalidArgument;
    };
    if (!std.math.isFinite(coeff) or coeff < 0.0) {
        std.debug.print("Error: Bitmask correction coefficient must be finite and non-negative: {d}\n", .{coeff});
        return error.InvalidArgument;
    }
    return coeff;
}

fn parseBitmaskLutMode(value: []const u8) ArgError!BitmaskLutMode {
    if (std.mem.eql(u8, value, "single")) return .single;
    if (std.mem.eql(u8, value, "per-frame")) return .per_frame;
    if (std.mem.eql(u8, value, "cycle")) return .cycle;
    std.debug.print("Error: Invalid bitmask LUT mode: {s}\n", .{value});
    return error.InvalidArgument;
}

// Range checks shared by the argument parser and `validateArgs`. The ranges
// are the ones `calc` enforces.

fn checkProbeRadius(radius: f64) ArgError!void {
    _ = calc.validateWorkflowProbeRadius(radius) catch {
        std.debug.print("Error: Probe radius must be between 0 and 10 Angstroms: {d}\n", .{radius});
        return error.InvalidArgument;
    };
}

fn checkNPoints(n: u32) ArgError!void {
    _ = calc.validateWorkflowNPoints(n) catch {
        std.debug.print("Error: n-points must be between 1 and 10000: {d}\n", .{n});
        return error.InvalidArgument;
    };
}

fn checkNSlices(n: u32) ArgError!void {
    _ = calc.validateWorkflowNSlices(n) catch {
        std.debug.print("Error: n-slices must be between 1 and 1000: {d}\n", .{n});
        return error.InvalidArgument;
    };
}

fn checkStride(stride: u32) ArgError!void {
    if (stride == 0) {
        std.debug.print("Error: Stride must be >= 1: {d}\n", .{stride});
        return error.InvalidArgument;
    }
}

fn validateBitmaskLutMode(args: TrajArgs) !void {
    if (args.bitmask_lut_mode != .single and !args.use_bitmask) {
        std.debug.print("Error: --bitmask-lut-mode requires --use-bitmask\n", .{});
        return error.InvalidArgument;
    }
}

/// Check every option value and combination that does not depend on the input
/// files. `run` calls this before it opens or creates any file, so an invalid
/// option can never truncate an existing results file.
fn validateArgs(args: TrajArgs) !void {
    try checkProbeRadius(args.probe_radius);
    try checkNPoints(args.n_points);
    try checkNSlices(args.n_slices);
    try checkStride(args.stride);
    if (args.output_path.len == 0) {
        std.debug.print("Error: Output path is empty\n", .{});
        return error.InvalidArgument;
    }
    if (args.bitmask_correction and !args.use_bitmask) {
        std.debug.print("Error: --bitmask-correction requires --use-bitmask\n", .{});
        return error.InvalidArgument;
    }
    try validateBitmaskLutMode(args);
    if (args.use_bitmask) {
        if (args.algorithm != .sr) {
            std.debug.print("Error: --use-bitmask requires --algorithm=sr\n", .{});
            return error.BitmaskRequiresSR;
        }
        if (!bitmask_lut.isSupportedNPoints(args.n_points)) {
            std.debug.print("Error: --use-bitmask requires --n-points=1..1024\n", .{});
            return error.UnsupportedNPoints;
        }
    }
}

/// Parse trajectory mode arguments. Prints a message and exits on invalid input.
pub fn parseArgs(args: []const []const u8, start_idx: usize) TrajArgs {
    return parseArgsChecked(args, start_idx) catch std.process.exit(1);
}

/// Parse trajectory mode arguments. Invalid input is reported on stderr and
/// returned as an error instead of ending the process.
fn parseArgsChecked(args: []const []const u8, start_idx: usize) ArgError!TrajArgs {
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
                result.algorithm = try parseAlgorithm(value);
            } else if (std.mem.startsWith(u8, arg, "--threads=")) {
                const value = arg["--threads=".len..];
                result.n_threads = std.fmt.parseInt(usize, value, 10) catch {
                    std.debug.print("Error: Invalid thread count: {s}\n", .{value});
                    return error.InvalidArgument;
                };
            } else if (std.mem.startsWith(u8, arg, "--probe-radius=")) {
                const value = arg["--probe-radius=".len..];
                result.probe_radius = std.fmt.parseFloat(f64, value) catch {
                    std.debug.print("Error: Invalid probe radius: {s}\n", .{value});
                    return error.InvalidArgument;
                };
                try checkProbeRadius(result.probe_radius);
            } else if (std.mem.startsWith(u8, arg, "--n-points=")) {
                const value = arg["--n-points=".len..];
                result.n_points = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid n-points: {s}\n", .{value});
                    return error.InvalidArgument;
                };
                try checkNPoints(result.n_points);
            } else if (std.mem.startsWith(u8, arg, "--n-slices=")) {
                const value = arg["--n-slices=".len..];
                result.n_slices = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid n-slices: {s}\n", .{value});
                    return error.InvalidArgument;
                };
                try checkNSlices(result.n_slices);
            } else if (std.mem.startsWith(u8, arg, "--precision=")) {
                const value = arg["--precision=".len..];
                result.precision = try parsePrecision(value);
            } else if (std.mem.startsWith(u8, arg, "--classifier=")) {
                const value = arg["--classifier=".len..];
                result.classifier_type = try parseClassifierType(value);
            } else if (std.mem.startsWith(u8, arg, "--ccd=")) {
                const value = arg["--ccd=".len..];
                result.ccd_path = value;
            } else if (std.mem.eql(u8, arg, "--ccd")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --ccd\n", .{});
                    return error.InvalidArgument;
                }
                result.ccd_path = args[i];
            } else if (std.mem.startsWith(u8, arg, "--sdf=")) {
                const value = arg["--sdf=".len..];
                result.sdf_paths.append(value) catch {
                    std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                    return error.InvalidArgument;
                };
            } else if (std.mem.eql(u8, arg, "--sdf")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --sdf\n", .{});
                    return error.InvalidArgument;
                }
                result.sdf_paths.append(args[i]) catch {
                    std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                    return error.InvalidArgument;
                };
            } else if (std.mem.startsWith(u8, arg, "--stride=")) {
                const value = arg["--stride=".len..];
                result.stride = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid stride: {s}\n", .{value});
                    return error.InvalidArgument;
                };
                try checkStride(result.stride);
            } else if (std.mem.startsWith(u8, arg, "--start=")) {
                const value = arg["--start=".len..];
                result.start_frame = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid start frame: {s}\n", .{value});
                    return error.InvalidArgument;
                };
            } else if (std.mem.startsWith(u8, arg, "--end=")) {
                const value = arg["--end=".len..];
                result.end_frame = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid end frame: {s}\n", .{value});
                    return error.InvalidArgument;
                };
            } else if (std.mem.startsWith(u8, arg, "--batch-size=")) {
                const value = arg["--batch-size=".len..];
                const bs = std.fmt.parseInt(u32, value, 10) catch {
                    std.debug.print("Error: Invalid batch size: {s}\n", .{value});
                    return error.InvalidArgument;
                };
                if (bs == 0) {
                    std.debug.print("Error: Batch size must be >= 1 (omit for auto)\n", .{});
                    return error.InvalidArgument;
                }
                result.batch_size = bs;
            } else if (std.mem.startsWith(u8, arg, "--output=")) {
                result.output_path = arg["--output=".len..];
            } else if (std.mem.eql(u8, arg, "--output")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --output\n", .{});
                    return error.InvalidArgument;
                }
                result.output_path = args[i];
            } else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
                result.include_hydrogens = true;
            } else if (std.mem.eql(u8, arg, "--no-hydrogens") or std.mem.eql(u8, arg, "--exclude-hydrogens")) {
                result.include_hydrogens = false;
            } else if (std.mem.startsWith(u8, arg, "--altloc=")) {
                const setting = try parseAltLocSetting(arg["--altloc=".len..]);
                result.alt_loc_mode = setting.mode;
                result.alt_loc_id = setting.id;
            } else if (std.mem.eql(u8, arg, "--altloc")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --altloc\n", .{});
                    return error.InvalidArgument;
                }
                const setting = try parseAltLocSetting(args[i]);
                result.alt_loc_mode = setting.mode;
                result.alt_loc_id = setting.id;
            } else if (std.mem.eql(u8, arg, "--use-bitmask")) {
                result.use_bitmask = true;
            } else if (std.mem.startsWith(u8, arg, "--bitmask-lut-mode=")) {
                result.bitmask_lut_mode = try parseBitmaskLutMode(arg["--bitmask-lut-mode=".len..]);
            } else if (std.mem.eql(u8, arg, "--bitmask-lut-mode")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --bitmask-lut-mode\n", .{});
                    return error.InvalidArgument;
                }
                result.bitmask_lut_mode = try parseBitmaskLutMode(args[i]);
            } else if (std.mem.eql(u8, arg, "--bitmask-correction")) {
                result.bitmask_correction = true;
            } else if (std.mem.startsWith(u8, arg, "--bitmask-correction-coeff=")) {
                result.bitmask_correction_coeff = try parseBitmaskCorrectionCoeff(arg["--bitmask-correction-coeff=".len..]);
                result.bitmask_correction = true;
            } else if (std.mem.eql(u8, arg, "--bitmask-correction-coeff")) {
                i += 1;
                if (i >= args.len) {
                    std.debug.print("Error: Missing value for --bitmask-correction-coeff\n", .{});
                    return error.InvalidArgument;
                }
                result.bitmask_correction_coeff = try parseBitmaskCorrectionCoeff(args[i]);
                result.bitmask_correction = true;
            } else if (std.mem.eql(u8, arg, "-q") or std.mem.eql(u8, arg, "--quiet")) {
                result.quiet = true;
                result.show_progress = false;
            } else {
                std.debug.print("Error: Unknown option: {s}\n", .{arg});
                return error.InvalidArgument;
            }
        } else if (std.mem.eql(u8, arg, "-h")) {
            result.show_help = true;
        } else if (std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
            result.show_progress = false;
        } else if (std.mem.eql(u8, arg, "-o")) {
            // -o FILE
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for -o\n", .{});
                return error.InvalidArgument;
            }
            result.output_path = args[i];
        } else if (std.mem.startsWith(u8, arg, "-o=")) {
            // -o=FILE
            result.output_path = arg["-o=".len..];
        } else if (std.mem.startsWith(u8, arg, "-")) {
            // Any other dash argument (including -oFILE) is not an option we know
            std.debug.print("Error: Unknown option: {s}\n", .{arg});
            return error.InvalidArgument;
        } else {
            // Positional arguments
            if (positional_count == 0) {
                result.traj_path = arg;
            } else if (positional_count == 1) {
                result.topology_path = arg;
            } else {
                std.debug.print("Error: Too many positional arguments\n", .{});
                return error.InvalidArgument;
            }
            positional_count += 1;
        }
    }

    return result;
}

fn parseAlgorithm(value: []const u8) ArgError!Algorithm {
    if (std.mem.eql(u8, value, "sr") or std.mem.eql(u8, value, "shrake-rupley")) {
        return .sr;
    } else if (std.mem.eql(u8, value, "lr") or std.mem.eql(u8, value, "lee-richards")) {
        return .lr;
    } else {
        std.debug.print("Error: Invalid algorithm: {s}\n", .{value});
        return error.InvalidArgument;
    }
}

fn parsePrecision(value: []const u8) ArgError!Precision {
    if (std.mem.eql(u8, value, "f32")) {
        return .f32;
    } else if (std.mem.eql(u8, value, "f64")) {
        return .f64;
    } else {
        std.debug.print("Error: Invalid precision: {s}\n", .{value});
        return error.InvalidArgument;
    }
}

fn parseAltLocSetting(value: []const u8) ArgError!mmcif_parser.AltLocSetting {
    return mmcif_parser.parseAltLocSetting(value) orelse {
        std.debug.print("Error: Invalid altloc mode: {s}\n", .{value});
        std.debug.print("Valid altloc modes: auto, none, all, highest-occupancy, or a single altLoc ID like A\n", .{});
        return error.InvalidArgument;
    };
}

fn parseClassifierType(value: []const u8) ArgError!ClassifierType {
    if (ClassifierType.fromString(value)) |ct| {
        return ct;
    } else {
        std.debug.print("Error: Invalid classifier: {s}\n", .{value});
        std.debug.print("Valid classifiers: ccd, protor, naccess, oons\n", .{});
        return error.InvalidArgument;
    }
}

/// Print help for trajectory mode
pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\Usage: {s} traj <trajectory> <topology> [options]
        \\
        \\Calculate SASA for each frame in a trajectory.
        \\Supported formats: XTC, TRR (GROMACS), DCD (NAMD/CHARMM), AMBER NetCDF.
        \\Format is auto-detected from file extension.
        \\
        \\ARGUMENTS:
        \\    <trajectory> Trajectory file (.xtc, .trr, .dcd, .nc, or .ncdf)
        \\    <topology>   Topology file (PDB or mmCIF) for atom names and radii.
        \\                 Must list the atoms of the trajectory in the same order:
        \\                 all ATOM and HETATM records of its first model are read.
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: ccd, protor, naccess, oons
        \\                       Default: naccess (supports explicit H in MD trajectories)
        \\    --ccd=PATH         External CCD dictionary file (.zsdc or .cif[.gz|.zst])
        \\                       Used with --classifier=ccd for non-standard residues
        \\    --sdf=PATH         SDF file with bond topology for CCD classifier
        \\                       Can be specified multiple times for multiple ligands
        \\    --threads=N        Number of threads (default: auto-detect)
        \\    --probe-radius=R   Probe radius in Angstroms, 0 < R <= 10 (default: 1.4)
        \\    --n-points=N       Test points per atom, 1..10000 (default: 100, for sr)
        \\    --n-slices=N       Slices per atom diameter, 1..1000 (default: 20, for lr)
        \\    --precision=PREC    Floating-point precision: f32, f64 (default: f32)
        \\    --no-hydrogens     Exclude hydrogen atoms from the calculation; they stay
        \\                       in the topology and trajectory files (default: included)
        \\    --include-hydrogens
        \\                       Include hydrogen atoms (default, for backward compat)
        \\    --altloc=MODE      Topology alternate-location handling (PDB/mmCIF): auto,
        \\                       none, all, highest-occupancy, or a single ID like A
        \\                       (default: auto)
        \\    --use-bitmask      Use bitmask LUT optimization for SR algorithm
        \\                       (n-points must be 1..1024)
        \\    --bitmask-lut-mode=MODE
        \\                       Trajectory bitmask LUT reuse mode:
        \\                       single (default), per-frame, cycle
        \\                       Non-default modes require --use-bitmask
        \\    --bitmask-correction
        \\                       Experimental correction for bitmask quantization bias
        \\                       Requires --use-bitmask
        \\    --bitmask-correction-coeff=V
        \\                       Override correction coefficient (default: 0.020)
        \\    --stride=N         Process every Nth frame, N >= 1 (default: 1)
        \\    --start=N          Start from frame N (default: 0)
        \\    --end=N            End at frame N (default: all)
        \\    --batch-size=N     Frames per batch for parallel processing
        \\                       Default: auto (threads * 2)
        \\    -o FILE, --output=FILE
        \\                       Output CSV file (default: traj_sasa.csv)
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
        \\    {s} traj trajectory.trr topology.pdb
        \\    {s} traj trajectory.nc topology.pdb
        \\    {s} traj trajectory.xtc topology.pdb -o sasa.csv
        \\    {s} traj trajectory.xtc topology.pdb --stride=10
        \\    {s} traj trajectory.xtc topology.pdb --classifier=naccess
        \\    {s} traj trajectory.xtc topology.pdb --algorithm=lr --n-slices=50
        \\    {s} traj trajectory.xtc topology.pdb --threads=8
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
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

fn estimateProcessedFrames(args: TrajArgs) usize {
    const end = args.end_frame orelse return 0;
    if (end < args.start_frame) return 0;
    return @as(usize, @intCast((end - args.start_frame) / args.stride)) + 1;
}

/// SASA result for one frame
const FrameResult = struct {
    frame_idx: u32 = 0,
    step: i32 = 0,
    time: f32 = 0,
    total_sasa: f64 = 0,
    /// Set by the worker once the frame is computed. After a failed batch the
    /// frames before the first unfinished one are still written.
    done: bool = false,
};

/// Helper to build and hold bitmask LUTs for trajectory processing.
/// Builds the appropriate LUT once based on args, and provides typed pointers.
const TrajLuts = struct {
    lut_f64: ?bitmask_lut.BitmaskLut = null,
    lut_f32: ?bitmask_lut.BitmaskLutGen(f32) = null,
    cycle_f64: [bitmask_lut_cycle_count]bitmask_lut.BitmaskLut = undefined,
    cycle_f32: [bitmask_lut_cycle_count]bitmask_lut.BitmaskLutGen(f32) = undefined,
    cycle_len: usize = 0,
    cycle_precision: Precision = .f64,

    fn init(allocator: Allocator, args: TrajArgs) !TrajLuts {
        if (!args.use_bitmask) return .{};
        if (args.algorithm != .sr) return error.BitmaskRequiresSR;
        if (!bitmask_lut.isSupportedNPoints(args.n_points)) return error.UnsupportedNPoints;
        switch (args.bitmask_lut_mode) {
            .per_frame => return .{},
            .single => return switch (args.precision) {
                .f64 => .{ .lut_f64 = try bitmask_lut.BitmaskLut.init(allocator, args.n_points) },
                .f32 => .{ .lut_f32 = try bitmask_lut.BitmaskLutGen(f32).init(allocator, args.n_points) },
            },
            .cycle => {
                var luts = TrajLuts{
                    .cycle_precision = args.precision,
                };
                errdefer luts.deinit();
                switch (args.precision) {
                    .f64 => for (0..bitmask_lut_cycle_count) |i| {
                        luts.cycle_f64[i] = try bitmask_lut.BitmaskLut.initVariant(allocator, args.n_points, i);
                        luts.cycle_len = i + 1;
                    },
                    .f32 => for (0..bitmask_lut_cycle_count) |i| {
                        luts.cycle_f32[i] = try bitmask_lut.BitmaskLutGen(f32).initVariant(allocator, args.n_points, i);
                        luts.cycle_len = i + 1;
                    },
                }
                return luts;
            },
        }
    }

    fn deinit(self: *TrajLuts) void {
        if (self.lut_f64) |*lut| lut.deinit();
        if (self.lut_f32) |*lut| lut.deinit();
        if (self.cycle_len > 0) {
            switch (self.cycle_precision) {
                .f64 => for (self.cycle_f64[0..self.cycle_len]) |*lut| lut.deinit(),
                .f32 => for (self.cycle_f32[0..self.cycle_len]) |*lut| lut.deinit(),
            }
        }
        self.* = .{};
    }

    fn f64Ptr(self: *const TrajLuts, frame_idx: usize) ?*const bitmask_lut.BitmaskLut {
        if (self.lut_f64 != null) return &self.lut_f64.?;
        if (self.cycle_len > 0 and self.cycle_precision == .f64) {
            return &self.cycle_f64[frame_idx % self.cycle_len];
        }
        return null;
    }

    fn f32Ptr(self: *const TrajLuts, frame_idx: usize) ?*const bitmask_lut.BitmaskLutGen(f32) {
        if (self.lut_f32 != null) return &self.lut_f32.?;
        if (self.cycle_len > 0 and self.cycle_precision == .f32) {
            return &self.cycle_f32[frame_idx % self.cycle_len];
        }
        return null;
    }

    fn reusableLutCount(self: *const TrajLuts) usize {
        if (self.cycle_len > 0) return self.cycle_len;
        if (self.lut_f64 != null or self.lut_f32 != null) return 1;
        return 0;
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
    coord_scale: f64, // ztraj readers already yield Å; kept for reader abstraction
    use_bitmask: bool = false,
    bitmask_luts: ?*const TrajLuts = null,
    bitmask_lut_f64: ?*const bitmask_lut.BitmaskLut = null,
    bitmask_lut_f32: ?*const bitmask_lut.BitmaskLutGen(f32) = null,
    bitmask_correction: bool = false,
    bitmask_correction_coeff: f64 = shrake_rupley_bitmask.default_bitmask_correction_coeff,
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

        // Convert f32 coordinates to f64. ztraj readers already yield Å.
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
        const frame_slot: usize = @intCast(frame_id);

        if (args.precision == .f32) {
            if (args.algorithm == .sr) {
                const config = Configf32{
                    .probe_radius = @floatCast(args.probe_radius),
                    .n_points = args.n_points,
                };
                const lut_opt = if (args.bitmask_luts) |luts| luts.f32Ptr(frame_slot) else args.bitmask_lut_f32;
                if (lut_opt) |lut| {
                    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f32){
                        .enabled = args.bitmask_correction,
                        .coeff = @floatCast(args.bitmask_correction_coeff),
                    };
                    var result = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithLutAndCorrection(thread_alloc, input, config, lut, correction) catch |err| {
                        setWorkerError(args, "frame {d}: bitmask-f32 failed: {s}", .{ frame_id, @errorName(err) });
                        return;
                    };
                    total_sasa = @floatCast(result.total_area);
                    result.deinit();
                } else if (args.use_bitmask) {
                    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f32){
                        .enabled = args.bitmask_correction,
                        .coeff = @floatCast(args.bitmask_correction_coeff),
                    };
                    var result = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithCorrection(thread_alloc, input, config, correction) catch |err| {
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
                const config = lee_richards.LeeRichardsConfigf32{
                    .probe_radius = @floatCast(args.probe_radius),
                    .n_slices = args.n_slices,
                };
                var result = lee_richards.calculateSasaf32(thread_alloc, input, config) catch |err| {
                    setWorkerError(args, "frame {d}: LR-f32 SASA failed: {s}", .{ frame_id, @errorName(err) });
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
                const lut_opt = if (args.bitmask_luts) |luts| luts.f64Ptr(frame_slot) else args.bitmask_lut_f64;
                if (lut_opt) |lut| {
                    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
                        .enabled = args.bitmask_correction,
                        .coeff = args.bitmask_correction_coeff,
                    };
                    var result = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithLutAndCorrection(thread_alloc, input, config, lut, correction) catch |err| {
                        setWorkerError(args, "frame {d}: bitmask-f64 failed: {s}", .{ frame_id, @errorName(err) });
                        return;
                    };
                    total_sasa = result.total_area;
                    result.deinit();
                } else if (args.use_bitmask) {
                    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
                        .enabled = args.bitmask_correction,
                        .coeff = args.bitmask_correction_coeff,
                    };
                    var result = shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithCorrection(thread_alloc, input, config, correction) catch |err| {
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
            .done = true,
        };

        // Reset arena for next frame (retain backing memory to avoid syscalls)
        _ = arena.reset(.retain_capacity);
    }
}

/// Copy the coordinates of one frame, keeping only the atoms selected from the
/// topology. `atom_indices` holds the trajectory atom index of each kept atom
/// (null = every atom, in order). `dst` holds `natoms * 3` values.
fn gatherFrameCoords(dst: []f32, src: []const f32, atom_indices: ?[]const u32) void {
    if (atom_indices) |indices| {
        for (indices, 0..) |src_atom, i| {
            dst[i * 3 ..][0..3].* = src[@as(usize, src_atom) * 3 ..][0..3].*;
        }
    } else {
        @memcpy(dst, src[0..dst.len]);
    }
}

/// Outcome of `readBatch`.
const BatchRead = struct {
    /// Frames placed in the batch buffers.
    count: usize,
    /// No further frames will be read (end of file, end of range, or error).
    eof: bool,
    /// Reader failure hit after `count` frames. Those frames are still valid,
    /// so the caller computes and writes them before reporting the error.
    read_error: ?anyerror = null,
};

/// Read frames into batch buffer, applying stride/start/end filtering.
/// Returns the number of frames read and whether EOF was reached.
fn readBatch(
    reader: *TrajectoryReader,
    allocator: Allocator,
    coord_pool: []f32,
    frame_data: []FrameData,
    batch_size: usize,
    natoms: usize,
    atom_indices: ?[]const u32,
    frame_idx: *u32,
    traj_args: TrajArgs,
) BatchRead {
    var batch_count: usize = 0;

    while (batch_count < batch_size) {
        var frame = reader.readFrame(allocator) catch |err| {
            if (reader.isEof(err)) return .{ .count = batch_count, .eof = true };
            return .{ .count = batch_count, .eof = true, .read_error = err };
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
        gatherFrameCoords(coord_pool[offset .. offset + natoms * 3], frame.coords, atom_indices);

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
// Topology
// =============================================================================

/// Topology atoms and how they map onto the atoms of a trajectory frame.
const Topology = struct {
    /// Atoms SASA is calculated for, in file order.
    atoms: AtomInput,
    /// Atoms read from the topology file. The trajectory must have this many.
    n_file_atoms: usize,
    /// Trajectory atom index of each atom in `atoms`, or null when every
    /// topology atom is used (`atoms[i]` is trajectory atom `i`).
    atom_indices: ?[]u32,
    allocator: Allocator,

    fn deinit(self: *Topology) void {
        self.atoms.deinit();
        if (self.atom_indices) |indices| self.allocator.free(indices);
    }
};

/// Copy the entries of `src` selected by `indices`.
fn gatherAtoms(comptime T: type, allocator: Allocator, src: []const T, indices: []const u32) ![]T {
    const out = try allocator.alloc(T, indices.len);
    for (indices, out) |src_idx, *dst| dst.* = src[src_idx];
    return out;
}

/// Build an AtomInput holding only the atoms selected by `indices`. Copies the
/// fields the trajectory pipeline uses (coordinates, radii and the names the
/// classifiers read).
fn selectTopologyAtoms(allocator: Allocator, full: *const AtomInput, indices: []const u32) !AtomInput {
    var out = AtomInput{ .x = &.{}, .y = &.{}, .z = &.{}, .r = &.{}, .allocator = allocator };
    errdefer out.deinit();
    out.x = try gatherAtoms(f64, allocator, full.x, indices);
    out.y = try gatherAtoms(f64, allocator, full.y, indices);
    out.z = try gatherAtoms(f64, allocator, full.z, indices);
    out.r = try gatherAtoms(f64, allocator, full.r, indices);
    if (full.residue) |v| out.residue = try gatherAtoms(types.FixedString5, allocator, v, indices);
    if (full.atom_name) |v| out.atom_name = try gatherAtoms(types.FixedString4, allocator, v, indices);
    if (full.element) |v| out.element = try gatherAtoms(u8, allocator, v, indices);
    return out;
}

/// Read the topology.
///
/// A topology has to map one-to-one onto the trajectory atoms, so every atom
/// record (ATOM and HETATM, hydrogens included) of the first model is read, in
/// file order. With `--no-hydrogens` the hydrogens are removed afterwards and
/// the indices of the remaining atoms are kept, so the same atoms can be
/// picked out of each frame.
fn loadTopology(allocator: Allocator, io: std.Io, path: []const u8, args: TrajArgs) !Topology {
    var hydrogen_flags: std.ArrayListUnmanaged(bool) = .empty;
    defer hydrogen_flags.deinit(allocator);
    const flags_out: ?*std.ArrayListUnmanaged(bool) = if (args.include_hydrogens) null else &hydrogen_flags;

    var full = switch (detectTopologyFormat(path)) {
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.atom_only = false;
            parser.skip_hydrogens = false;
            parser.first_model_only = true;
            parser.hydrogen_flags = flags_out;
            parser.alt_loc_mode = args.alt_loc_mode;
            parser.alt_loc_id = args.alt_loc_id;
            break :blk try parser.parseFile(io, path);
        },
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.atom_only = false;
            parser.skip_hydrogens = false;
            parser.first_model_only = true;
            parser.hydrogen_flags = flags_out;
            parser.alt_loc_mode = args.alt_loc_mode;
            parser.alt_loc_id = args.alt_loc_id;
            parser.parse_inline_ccd = false;
            break :blk try parser.parseFile(io, path);
        },
    };
    errdefer full.deinit();
    const n_file_atoms = full.atomCount();

    const n_hydrogens = std.mem.count(bool, hydrogen_flags.items, &.{true});
    if (n_hydrogens == 0) {
        return .{ .atoms = full, .n_file_atoms = n_file_atoms, .atom_indices = null, .allocator = allocator };
    }
    if (n_hydrogens == n_file_atoms) {
        std.debug.print("Error: --no-hydrogens leaves no atoms: the topology contains only hydrogens\n", .{});
        return error.NoAtomsFound;
    }

    const indices = try allocator.alloc(u32, n_file_atoms - n_hydrogens);
    errdefer allocator.free(indices);
    var n_kept: usize = 0;
    for (hydrogen_flags.items, 0..) |is_hydrogen, i| {
        if (is_hydrogen) continue;
        indices[n_kept] = @intCast(i);
        n_kept += 1;
    }

    const kept = try selectTopologyAtoms(allocator, &full, indices);
    full.deinit();
    return .{ .atoms = kept, .n_file_atoms = n_file_atoms, .atom_indices = indices, .allocator = allocator };
}

// =============================================================================
// Main Entry Point
// =============================================================================

/// Run trajectory analysis
pub fn run(allocator: Allocator, io: std.Io, args: TrajArgs) !void {
    // Validate required arguments
    const traj_path = args.traj_path orelse {
        std.debug.print("Error: Missing trajectory file\n", .{});
        return error.MissingArgument;
    };
    const topology_path = args.topology_path orelse {
        std.debug.print("Error: Missing topology file\n", .{});
        return error.MissingArgument;
    };
    try validateArgs(args);

    // Detect trajectory format
    const traj_format = detectTrajectoryFormat(traj_path) orelse {
        std.debug.print("Error: Unknown trajectory format. Supported: .xtc, .trr, .dcd, .nc, .ncdf\n", .{});
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

    var loaded = try loadTopology(allocator, io, topology_path, args);
    defer loaded.deinit();
    const topology = &loaded.atoms;
    const atom_indices: ?[]const u32 = loaded.atom_indices;

    const natoms = topology.atomCount();
    if (!args.quiet) {
        if (natoms == loaded.n_file_atoms) {
            std.debug.print("Topology: {d} atoms\n", .{natoms});
        } else {
            std.debug.print("Topology: {d} atoms ({d} after excluding hydrogens)\n", .{ loaded.n_file_atoms, natoms });
        }
    }

    // Load external CCD dictionary if specified
    const use_ccd_resources = args.classifier_type == .ccd;
    var ext_ccd: ?ccd_parser.ComponentDict = null;
    if (use_ccd_resources) {
        if (args.ccd_path) |ccd_path| {
            const ccd_data = if (compressed.isCompressed(ccd_path))
                try compressed.read(allocator, ccd_path)
            else blk: {
                const f = try std.Io.Dir.cwd().openFile(io, ccd_path, .{});
                defer f.close(io);
                var read_buf_ccd: [65536]u8 = undefined;
                var file_r_ccd = f.reader(io, &read_buf_ccd);
                break :blk try file_r_ccd.interface.allocRemaining(allocator, .unlimited);
            };
            defer allocator.free(ccd_data);

            ext_ccd = try ccd_binary.loadDict(allocator, ccd_data);
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

    // Apply classifier if specified
    if (args.classifier_type) |ct| {
        const sdf_ccd_ptr: ?*const ccd_parser.ComponentDict = if (sdf_ccd != null) &sdf_ccd.? else null;
        const ext_ccd_ptr: ?*const ccd_parser.ComponentDict = if (ext_ccd != null) &ext_ccd.? else null;
        try applyBuiltinClassifier(allocator, topology, ct, sdf_ccd_ptr, ext_ccd_ptr, args.quiet);
    }

    // Open trajectory file
    if (!args.quiet) {
        const format_name: []const u8 = switch (traj_format) {
            .xtc => "XTC",
            .dcd => "DCD",
            .trr => "TRR",
            .nc => "AMBER NetCDF",
        };
        std.debug.print("Opening trajectory ({s}): {s}\n", .{ format_name, traj_path });
    }

    // Open trajectory reader and verify atom count
    var xtc_reader: ?xtc.XtcReader = null;
    var dcd_reader: ?dcd.DcdReader = null;
    var trr_reader: ?trr.TrrReader = null;
    var nc_reader: ?nc.NcReader = null;
    const traj_natoms: i32 = switch (traj_format) {
        .xtc => blk: {
            xtc_reader = try xtc.XtcReader.open(io, allocator, traj_path);
            break :blk @intCast(xtc_reader.?.nAtoms());
        },
        .dcd => blk: {
            dcd_reader = try dcd.DcdReader.open(io, allocator, traj_path);
            break :blk @intCast(dcd_reader.?.nAtoms());
        },
        .trr => blk: {
            trr_reader = try trr.TrrReader.open(io, allocator, traj_path);
            break :blk @intCast(trr_reader.?.nAtoms());
        },
        .nc => blk: {
            nc_reader = try nc.NcReader.open(io, allocator, traj_path);
            break :blk @intCast(nc_reader.?.nAtoms());
        },
    };
    defer {
        if (xtc_reader) |*r| r.deinit();
        if (dcd_reader) |*r| r.deinit();
        if (trr_reader) |*r| r.deinit();
        if (nc_reader) |*r| r.deinit();
    }

    // Verify atom count matches. The trajectory is compared with the whole
    // topology: atoms removed by --no-hydrogens are skipped in every frame, so
    // they still have to be present in the trajectory.
    if (traj_natoms != @as(i32, @intCast(loaded.n_file_atoms))) {
        std.debug.print("Error: Atom count mismatch - trajectory has {d} atoms, topology has {d}\n", .{
            traj_natoms,
            loaded.n_file_atoms,
        });
        std.debug.print("       The topology must list exactly the atoms of the trajectory, in the same order\n", .{});
        std.debug.print("       (every ATOM and HETATM record of its first model is read).\n", .{});
        if (atom_indices != null and traj_natoms == @as(i32, @intCast(natoms))) {
            std.debug.print("       This trajectory matches the topology without its hydrogens. --no-hydrogens removes\n", .{});
            std.debug.print("       hydrogens from the frames as well, so use a hydrogen-free topology and drop the flag.\n", .{});
        }
        return error.AtomCountMismatch;
    }

    // Resolve thread count
    const n_threads = if (args.n_threads == 0)
        @as(usize, @intCast(std.Thread.getCpuCount() catch 1))
    else
        args.n_threads;

    // Prepare trajectory bitmask LUTs according to the selected reuse mode.
    // The option combinations were checked by validateArgs.
    var luts = TrajLuts.init(allocator, args) catch |err| {
        std.debug.print("Error: failed to build bitmask LUT: {s}\n", .{@errorName(err)});
        return err;
    };
    defer luts.deinit();

    // Open output file with buffered writer. This is the last step before
    // processing: every check that can fail without reading frames has passed,
    // so an existing file is only replaced by a run that is about to produce
    // results.
    const output_file = try std.Io.Dir.cwd().createFile(io, args.output_path, .{});
    defer output_file.close(io);
    var write_buffer: [4096]u8 = undefined;
    var buffered_writer = output_file.writer(io, &write_buffer);
    const writer = &buffered_writer.interface;

    // Write CSV header
    try writer.writeAll("frame,step,time,total_sasa\n");

    if (!args.quiet) {
        std.debug.print("Algorithm: {s}, Threads: {d}, Precision: {s}{s}\n", .{
            if (args.algorithm == .sr) "Shrake-Rupley" else "Lee-Richards",
            n_threads,
            if (args.precision == .f32) "f32" else "f64",
            if (args.use_bitmask) ", Bitmask: enabled" else "",
        });
        if (args.use_bitmask) {
            const mode_name: []const u8 = switch (args.bitmask_lut_mode) {
                .single => "single",
                .per_frame => "per-frame",
                .cycle => "cycle",
            };
            std.debug.print("Bitmask LUT mode: {s}\n", .{mode_name});
        }
        std.debug.print("Processing frames", .{});
        if (args.stride > 1) std.debug.print(" (stride={d})", .{args.stride});
        if (args.start_frame > 0) std.debug.print(" from {d}", .{args.start_frame});
        if (args.end_frame) |end| std.debug.print(" to {d}", .{end});
        std.debug.print("...\n\n", .{});
    }

    // Process frames
    var frame_idx: u32 = 0;
    var processed_count: u32 = 0;
    var timer_start = std.Io.Timestamp.now(io, .awake);

    var traj_reader = TrajectoryReader{
        .xtc_reader = if (xtc_reader) |*r| r else null,
        .dcd_reader = if (dcd_reader) |*r| r else null,
        .trr_reader = if (trr_reader) |*r| r else null,
        .nc_reader = if (nc_reader) |*r| r else null,
        .format = traj_format,
    };

    const process_result = blk: {
        var progress_root: std.Progress.Node = if (shouldShowProgress(args))
            std.Progress.start(io, .{ .root_name = "Processing frames", .estimated_total_items = estimateProcessedFrames(args) })
        else
            .none;
        defer progress_root.end();

        if (n_threads <= 1) {
            // Sequential path: single-threaded SASA per frame
            break :blk runSequential(allocator, &traj_reader, writer, topology, atom_indices, args, &luts, &frame_idx, &processed_count, progress_root);
        } else {
            // Batch parallel path: frame-level parallelism across threads
            break :blk runBatchParallel(allocator, &traj_reader, writer, topology, atom_indices, args, &luts, n_threads, natoms, &frame_idx, &processed_count, progress_root);
        }
    };

    // Flush buffered output. This also runs when a frame failed, so the rows
    // of the frames computed before the failure reach the file.
    const flush_result = writer.flush();
    process_result catch |err| {
        if (processed_count > 0) {
            std.debug.print("Error: stopped after {d} frames; their results were written to: {s}\n", .{ processed_count, args.output_path });
        }
        return err;
    };
    try flush_result;

    const elapsed_ns = timer_start.untilNow(io, .awake).nanoseconds;
    const elapsed_ms = @divTrunc(elapsed_ns, std.time.ns_per_ms);

    if (!args.quiet) {
        std.debug.print("\r  Processed {d} frames in {d}ms ({d:.1} frames/sec)\n", .{
            processed_count,
            elapsed_ms,
            if (elapsed_ms > 0) @as(f64, @floatFromInt(processed_count)) * 1000.0 / @as(f64, @floatFromInt(elapsed_ms)) else 0,
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
    atom_indices: ?[]const u32,
    args: TrajArgs,
    luts: *const TrajLuts,
    frame_idx: *u32,
    processed_count: *u32,
    progress_node: std.Progress.Node,
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
        var frame = reader.readFrame(allocator) catch |err| {
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

        // Update frame coordinates, keeping only the atoms selected from the
        // topology. ztraj readers already yield Å.
        for (0..natoms) |i| {
            const src: usize = if (atom_indices) |indices| indices[i] else i;
            frame_x[i] = @as(f64, frame.coords[src * 3 + 0]) * scale;
            frame_y[i] = @as(f64, frame.coords[src * 3 + 1]) * scale;
            frame_z[i] = @as(f64, frame.coords[src * 3 + 2]) * scale;
        }

        // Calculate SASA (single-threaded)
        const total_sasa: f64 = switch (args.precision) {
            .f32 => blk: {
                switch (args.algorithm) {
                    .sr => {
                        const config = Configf32{
                            .probe_radius = @floatCast(args.probe_radius),
                            .n_points = args.n_points,
                        };
                        if (luts.f32Ptr(@intCast(frame_idx.*))) |lut| {
                            const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f32){
                                .enabled = args.bitmask_correction,
                                .coeff = @floatCast(args.bitmask_correction_coeff),
                            };
                            var result = try shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithLutAndCorrection(allocator, frame_input, config, lut, correction);
                            defer result.deinit();
                            break :blk @floatCast(result.total_area);
                        } else if (args.use_bitmask) {
                            const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f32){
                                .enabled = args.bitmask_correction,
                                .coeff = @floatCast(args.bitmask_correction_coeff),
                            };
                            var result = try shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithCorrection(allocator, frame_input, config, correction);
                            defer result.deinit();
                            break :blk @floatCast(result.total_area);
                        } else {
                            var result = try shrake_rupley.calculateSasaf32(allocator, frame_input, config);
                            defer result.deinit();
                            break :blk @floatCast(result.total_area);
                        }
                    },
                    .lr => {
                        const config = lee_richards.LeeRichardsConfigf32{
                            .probe_radius = @floatCast(args.probe_radius),
                            .n_slices = args.n_slices,
                        };
                        var result = try lee_richards.calculateSasaf32(allocator, frame_input, config);
                        defer result.deinit();
                        break :blk @floatCast(result.total_area);
                    },
                }
            },
            .f64 => blk: {
                switch (args.algorithm) {
                    .sr => {
                        const config = Config{
                            .probe_radius = args.probe_radius,
                            .n_points = args.n_points,
                        };
                        if (luts.f64Ptr(@intCast(frame_idx.*))) |lut| {
                            const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
                                .enabled = args.bitmask_correction,
                                .coeff = args.bitmask_correction_coeff,
                            };
                            var result = try shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithLutAndCorrection(allocator, frame_input, config, lut, correction);
                            defer result.deinit();
                            break :blk result.total_area;
                        } else if (args.use_bitmask) {
                            const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
                                .enabled = args.bitmask_correction,
                                .coeff = args.bitmask_correction_coeff,
                            };
                            var result = try shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithCorrection(allocator, frame_input, config, correction);
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

        progress_node.completeOne();
    }
}

/// Batch parallel frame processing (frame-level parallelism)
fn runBatchParallel(
    allocator: Allocator,
    reader: *TrajectoryReader,
    writer: anytype,
    topology: *AtomInput,
    atom_indices: ?[]const u32,
    args: TrajArgs,
    luts: *const TrajLuts,
    n_threads: usize,
    natoms: usize,
    frame_idx: *u32,
    processed_count: *u32,
    progress_node: std.Progress.Node,
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
        const read_result = readBatch(
            reader,
            allocator,
            coord_pool,
            frame_data,
            batch_size,
            natoms,
            atom_indices,
            frame_idx,
            args,
        );

        if (read_result.count == 0) {
            if (read_result.read_error) |err| return err;
            break;
        }
        @memset(frame_results[0..read_result.count], .{});

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
                .use_bitmask = args.use_bitmask,
                .bitmask_luts = luts,
                .bitmask_correction = args.bitmask_correction,
                .bitmask_correction_coeff = args.bitmask_correction_coeff,
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

        // Phase 3: Write results (sequential, in frame order). When a worker
        // failed, the frames before the first unfinished one are still written
        // so the output stays a gap-free prefix of the trajectory.
        var completed: usize = 0;
        while (completed < read_result.count and frame_results[completed].done) completed += 1;

        for (0..completed) |i| {
            try writer.print("{d},{d},{d:.3},{d:.2}\n", .{
                frame_results[i].frame_idx,
                frame_results[i].step,
                frame_results[i].time,
                frame_results[i].total_sasa,
            });
        }

        processed_count.* += @intCast(completed);

        progress_node.setCompletedItems(processed_count.*);

        if (error_flag.load(.acquire)) {
            const msg = std.mem.sliceTo(&error_msg, 0);
            if (msg.len > 0) {
                std.debug.print("Error: {s}\n", .{msg});
            }
            return error.BatchCalculationFailed;
        }

        // A reader failure ends the run after the frames read before it.
        if (read_result.read_error) |err| return err;

        if (read_result.eof) break;
    }
}

/// Load SDF files and build a ComponentDict from their bond topology.
/// Delegates to sdf_parser.loadSdfComponents.
const loadSdfComponents = sdf_parser.loadSdfComponents;

fn applyBuiltinClassifier(
    _: Allocator,
    input: *AtomInput,
    ct: ClassifierType,
    sdf_ccd: ?*const ccd_parser.ComponentDict,
    external_ccd: ?*const ccd_parser.ComponentDict,
    quiet: bool,
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
            var loaded: usize = 0;
            const dicts: [2]?*const ccd_parser.ComponentDict = .{ sdf_ccd, external_ccd };
            for (dicts) |maybe_dict| {
                if (maybe_dict) |dict| {
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
                }
            }
            if (!quiet and loaded > 0) {
                std.debug.print("CCD: {d} non-standard components derived from CCD data\n", .{loaded});
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
        } else if (classifier.guessFallbackRadius(
            ct,
            if (input.element) |elements| elements[i] else null,
            residues[i].slice(),
            atom_names[i].slice(),
        )) |r| {
            // Fall back to element-based radius, or atom name-based without an element
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

test "traj NACCESS and OONS take the element of unlisted atoms from the element column" {
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

        try applyBuiltinClassifier(allocator, &input, ct, null, null, true);

        try std.testing.expectEqualSlices(f64, &expected, input.r);
    }
}

test "TrajArgs progress defaults to visible" {
    const args = [_][]const u8{ "zsasa", "traj", "traj.xtc", "topology.pdb" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(false, parsed.quiet);
    try std.testing.expectEqual(true, parsed.show_progress);
}

test "TrajArgs quiet disables progress" {
    const args = [_][]const u8{ "zsasa", "traj", "--quiet", "traj.xtc", "topology.pdb" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.quiet);
    try std.testing.expectEqual(false, parsed.show_progress);
}

test "TrajArgs --altloc modes" {
    const none_args = [_][]const u8{ "zsasa", "traj", "--altloc=none", "traj.xtc", "topology.cif" };
    const selected_args = [_][]const u8{ "zsasa", "traj", "--altloc", "D", "traj.xtc", "topology.cif" };

    const none = parseArgs(&none_args, 2);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.none, none.alt_loc_mode);

    const selected = parseArgs(&selected_args, 2);
    try std.testing.expectEqual(mmcif_parser.AltLocMode.selected, selected.alt_loc_mode);
    try std.testing.expectEqual(@as(u8, 'D'), selected.alt_loc_id);
}

test "TrajArgs --bitmask-correction-coeff" {
    const args = [_][]const u8{ "zsasa", "traj", "--use-bitmask", "--bitmask-correction-coeff=0.2", "traj.xtc", "topology.pdb" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.bitmask_correction);
    try std.testing.expectEqual(@as(f64, 0.2), parsed.bitmask_correction_coeff);
}

test "TrajArgs --bitmask-lut-mode parses supported modes" {
    const single_args = [_][]const u8{ "zsasa", "traj", "--use-bitmask", "--bitmask-lut-mode=single", "traj.xtc", "topology.pdb" };
    const per_frame_args = [_][]const u8{ "zsasa", "traj", "--use-bitmask", "--bitmask-lut-mode", "per-frame", "traj.xtc", "topology.pdb" };
    const cycle_args = [_][]const u8{ "zsasa", "traj", "--use-bitmask", "--bitmask-lut-mode=cycle", "traj.xtc", "topology.pdb" };

    try std.testing.expectEqual(BitmaskLutMode.single, parseArgs(&single_args, 2).bitmask_lut_mode);
    try std.testing.expectEqual(BitmaskLutMode.per_frame, parseArgs(&per_frame_args, 2).bitmask_lut_mode);
    try std.testing.expectEqual(BitmaskLutMode.cycle, parseArgs(&cycle_args, 2).bitmask_lut_mode);
}

test "TrajArgs non-default bitmask LUT mode requires bitmask" {
    try std.testing.expectError(
        error.InvalidArgument,
        validateBitmaskLutMode(.{ .bitmask_lut_mode = .per_frame }),
    );
    try validateBitmaskLutMode(.{ .use_bitmask = true, .bitmask_lut_mode = .per_frame });
}

test "TrajLuts honors bitmask LUT reuse modes" {
    const allocator = std.testing.allocator;

    var single = try TrajLuts.init(allocator, .{ .use_bitmask = true, .bitmask_lut_mode = .single, .precision = .f64 });
    defer single.deinit();
    try std.testing.expect(single.f64Ptr(0) != null);
    try std.testing.expectEqual(@as(usize, 1), single.reusableLutCount());

    var per_frame = try TrajLuts.init(allocator, .{ .use_bitmask = true, .bitmask_lut_mode = .per_frame, .precision = .f64 });
    defer per_frame.deinit();
    try std.testing.expect(per_frame.f64Ptr(0) == null);
    try std.testing.expectEqual(@as(usize, 0), per_frame.reusableLutCount());
    try std.testing.expectError(
        error.UnsupportedNPoints,
        TrajLuts.init(allocator, .{ .use_bitmask = true, .bitmask_lut_mode = .per_frame, .precision = .f64, .n_points = 2000 }),
    );

    var cycle = try TrajLuts.init(allocator, .{ .use_bitmask = true, .bitmask_lut_mode = .cycle, .precision = .f64 });
    defer cycle.deinit();
    try std.testing.expect(cycle.f64Ptr(0) != null);
    try std.testing.expect(!std.mem.eql(u64, cycle.f64Ptr(0).?.masks, cycle.f64Ptr(1).?.masks));
    try std.testing.expect(cycle.f64Ptr(bitmask_lut_cycle_count) == cycle.f64Ptr(0));
    try std.testing.expectEqual(bitmask_lut_cycle_count, cycle.reusableLutCount());
}

test "TrajArgs quiet suppresses progress even when show_progress defaults true" {
    const args = TrajArgs{ .quiet = true };
    try std.testing.expectEqual(false, shouldShowProgress(args));
}
test "TrajectoryReader uses angstrom coordinates from ztraj readers" {
    const xtc_reader = TrajectoryReader{ .format = .xtc };
    const dcd_reader = TrajectoryReader{ .format = .dcd };
    const trr_reader = TrajectoryReader{ .format = .trr };
    const nc_reader = TrajectoryReader{ .format = .nc };
    try std.testing.expectEqual(@as(f64, 1.0), xtc_reader.coordScale());
    try std.testing.expectEqual(@as(f64, 1.0), dcd_reader.coordScale());
    try std.testing.expectEqual(@as(f64, 1.0), trr_reader.coordScale());
    try std.testing.expectEqual(@as(f64, 1.0), nc_reader.coordScale());
}

test "detectTrajectoryFormat recognizes ztraj-backed formats" {
    try std.testing.expectEqual(TrajectoryFormat.xtc, detectTrajectoryFormat("traj.xtc").?);
    try std.testing.expectEqual(TrajectoryFormat.dcd, detectTrajectoryFormat("traj.dcd").?);
    try std.testing.expectEqual(TrajectoryFormat.trr, detectTrajectoryFormat("traj.trr").?);
    try std.testing.expectEqual(TrajectoryFormat.nc, detectTrajectoryFormat("traj.nc").?);
    try std.testing.expectEqual(TrajectoryFormat.nc, detectTrajectoryFormat("traj.ncdf").?);
    try std.testing.expectEqual(@as(?TrajectoryFormat, null), detectTrajectoryFormat("traj.gro"));
}

test "parseArgsChecked accepts -o FILE, -o=FILE, --output=FILE and --output FILE" {
    const forms = [_][]const []const u8{
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "-o", "out.csv" },
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "-o=out.csv" },
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "--output=out.csv" },
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "--output", "out.csv" },
    };
    for (forms) |args| {
        const parsed = try parseArgsChecked(args, 2);
        try std.testing.expectEqualStrings("out.csv", parsed.output_path);
        try std.testing.expectEqualStrings("traj.xtc", parsed.traj_path.?);
        try std.testing.expectEqualStrings("topology.pdb", parsed.topology_path.?);
    }
}

test "parseArgsChecked rejects unknown dash arguments and a positional output" {
    const rejected = [_][]const []const u8{
        // -oFILE used to be taken as the output path "XYZ"
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "-oXYZ" },
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "-x" },
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "-o" },
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "--output" },
        // The output file is an option, not a third positional argument
        &.{ "zsasa", "traj", "traj.xtc", "topology.pdb", "out.csv" },
    };
    for (rejected) |args| {
        try std.testing.expectError(error.InvalidArgument, parseArgsChecked(args, 2));
    }
}

test "parseArgsChecked range-checks stride, probe radius, n-points and n-slices" {
    const rejected = [_][]const u8{
        "--stride=0",
        "--probe-radius=0",
        "--probe-radius=-1",
        "--probe-radius=10.5",
        "--probe-radius=1e6",
        "--probe-radius=nan",
        "--n-points=0",
        "--n-points=10001",
        "--n-slices=0",
        "--n-slices=1001",
    };
    for (rejected) |option| {
        const args = [_][]const u8{ "zsasa", "traj", option, "traj.xtc", "topology.pdb" };
        try std.testing.expectError(error.InvalidArgument, parseArgsChecked(&args, 2));
    }

    const accepted = [_][]const u8{ "zsasa", "traj", "--stride=3", "--probe-radius=10", "--n-points=10000", "--n-slices=1000", "traj.xtc", "topology.pdb" };
    const parsed = try parseArgsChecked(&accepted, 2);
    try std.testing.expectEqual(@as(u32, 3), parsed.stride);
    try std.testing.expectEqual(@as(f64, 10.0), parsed.probe_radius);
    try std.testing.expectEqual(@as(u32, 10000), parsed.n_points);
    try std.testing.expectEqual(@as(u32, 1000), parsed.n_slices);
}

// =============================================================================
// End-to-end tests: run the command path on the repository fixtures
// =============================================================================

const test_xtc_path = "test_data/1l2y.xtc"; // 38 frames, 304 atoms
const test_pdb_path = "test_data/1l2y.pdb"; // 38 NMR models, 304 atoms each (150 hydrogens)

/// Temporary directory plus an arena for the paths and file contents of a test.
const TestWorkspace = struct {
    tmp: std.testing.TmpDir,
    arena: std.heap.ArenaAllocator,

    fn init() TestWorkspace {
        return .{ .tmp = std.testing.tmpDir(.{}), .arena = std.heap.ArenaAllocator.init(std.testing.allocator) };
    }

    fn deinit(self: *TestWorkspace) void {
        self.arena.deinit();
        self.tmp.cleanup();
    }

    fn allocator(self: *TestWorkspace) Allocator {
        return self.arena.allocator();
    }

    /// Absolute path of `name` inside the temporary directory.
    fn path(self: *TestWorkspace, name: []const u8) ![]const u8 {
        var root_buf: [std.fs.max_path_bytes]u8 = undefined;
        const root_len = try self.tmp.dir.realPath(std.testing.io, &root_buf);
        return std.fs.path.join(self.allocator(), &.{ root_buf[0..root_len], name });
    }

    /// Write `data` to `name` and return its absolute path.
    fn write(self: *TestWorkspace, name: []const u8, data: []const u8) ![]const u8 {
        const file_path = try self.path(name);
        try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = file_path, .data = data });
        return file_path;
    }

    fn read(self: *TestWorkspace, file_path: []const u8) ![]const u8 {
        return std.Io.Dir.cwd().readFileAlloc(std.testing.io, file_path, self.allocator(), .limited(8 << 20));
    }

    /// Run the trajectory command quietly and return the CSV it wrote.
    fn runTraj(self: *TestWorkspace, base: TrajArgs, out_name: []const u8) ![]const u8 {
        var args = base;
        args.output_path = try self.path(out_name);
        args.quiet = true;
        try run(std.testing.allocator, std.testing.io, args);
        return self.read(args.output_path);
    }

    /// PDB text of the first model of the 1L2Y fixture, up to and including ENDMDL.
    fn firstModelPdb(self: *TestWorkspace) ![]const u8 {
        const data = try self.read(test_pdb_path);
        const record = std.mem.indexOf(u8, data, "ENDMDL") orelse return error.TestUnexpectedResult;
        const line_end = std.mem.indexOfScalarPos(u8, data, record, '\n') orelse return error.TestUnexpectedResult;
        return data[0 .. line_end + 1];
    }

    /// The total_sasa column of a result CSV.
    fn totals(self: *TestWorkspace, csv: []const u8) ![]const f64 {
        var values: std.ArrayListUnmanaged(f64) = .empty;
        var lines = std.mem.tokenizeScalar(u8, csv, '\n');
        try std.testing.expectEqualStrings("frame,step,time,total_sasa", lines.next() orelse "");
        while (lines.next()) |line| {
            const comma = std.mem.lastIndexOfScalar(u8, line, ',') orelse return error.TestUnexpectedResult;
            try values.append(self.allocator(), try std.fmt.parseFloat(f64, line[comma + 1 ..]));
        }
        return values.items;
    }

    /// Copy the first `n_frames` frames of the XTC fixture into two TRR files:
    /// `full_name` gets every atom, `subset_name` only the atoms in `subset`.
    /// Both go through the same conversion, so a kept atom has bit-identical
    /// coordinates in the two files. `nan_frame` poisons one coordinate of
    /// that frame in the full file.
    fn writeTrrPair(
        self: *TestWorkspace,
        full_name: []const u8,
        subset_name: []const u8,
        subset: []const u32,
        n_frames: usize,
        nan_frame: ?usize,
    ) !void {
        const io = std.testing.io;
        const gpa = std.testing.allocator;

        var reader = try xtc.XtcReader.open(io, gpa, test_xtc_path);
        defer reader.deinit();
        const n_atoms: usize = reader.nAtoms();

        var full_writer = try trr.TrrWriter.open(io, gpa, try self.path(full_name), n_atoms);
        defer full_writer.deinit();
        var subset_writer = try trr.TrrWriter.open(io, gpa, try self.path(subset_name), subset.len);
        defer subset_writer.deinit();

        var full_frame = try ztraj.types.Frame.init(gpa, n_atoms);
        defer full_frame.deinit();
        var subset_frame = try ztraj.types.Frame.init(gpa, subset.len);
        defer subset_frame.deinit();

        for (0..n_frames) |frame_no| {
            const frame = (try reader.next()) orelse return error.TestUnexpectedResult;

            @memcpy(full_frame.x, frame.x);
            @memcpy(full_frame.y, frame.y);
            @memcpy(full_frame.z, frame.z);
            if (nan_frame == frame_no) full_frame.x[0] = std.math.nan(f32);
            for (subset, 0..) |src, i| {
                subset_frame.x[i] = frame.x[src];
                subset_frame.y[i] = frame.y[src];
                subset_frame.z[i] = frame.z[src];
            }
            inline for (.{ &full_frame, &subset_frame }) |out| {
                out.box_vectors = frame.box_vectors;
                out.time = frame.time;
                out.step = frame.step;
            }

            try full_writer.writeFrame(full_frame);
            try subset_writer.writeFrame(subset_frame);
        }

        try full_writer.close();
        try subset_writer.close();
    }
};

test "traj run: every precision, algorithm and bitmask combination agrees across paths, and f32 honors --algorithm=lr" {
    var ws = TestWorkspace.init();
    defer ws.deinit();
    const topology = try ws.write("m1.pdb", try ws.firstModelPdb());

    const Combo = struct {
        algorithm: Algorithm,
        use_bitmask: bool = false,
        bitmask_lut_mode: BitmaskLutMode = .single,
    };
    const combos = [_]Combo{
        .{ .algorithm = .sr },
        .{ .algorithm = .sr, .use_bitmask = true },
        .{ .algorithm = .sr, .use_bitmask = true, .bitmask_lut_mode = .per_frame },
        .{ .algorithm = .sr, .use_bitmask = true, .bitmask_lut_mode = .cycle },
        .{ .algorithm = .lr },
    };
    const precisions = [_]Precision{ .f32, .f64 };

    // totals[combo][precision] = total SASA of frames 0 and 1
    var totals: [combos.len][precisions.len][]const f64 = undefined;
    for (combos, 0..) |combo, c| {
        for (precisions, 0..) |precision, p| {
            var args = TrajArgs{
                .traj_path = test_xtc_path,
                .topology_path = topology,
                .algorithm = combo.algorithm,
                .precision = precision,
                .n_points = 64,
                .end_frame = 1,
                .use_bitmask = combo.use_bitmask,
                .bitmask_lut_mode = combo.bitmask_lut_mode,
            };
            args.n_threads = 1;
            const sequential = try ws.runTraj(args, "sequential.csv");
            args.n_threads = 2;
            const parallel = try ws.runTraj(args, "parallel.csv");

            // The sequential and the batch path must select the same kernel.
            try std.testing.expectEqualStrings(sequential, parallel);
            totals[c][p] = try ws.totals(sequential);
            try std.testing.expectEqual(@as(usize, 2), totals[c][p].len);
        }
        // f32 and f64 run the same algorithm, so they agree to f32 accuracy.
        for (totals[c][0], totals[c][1]) |total_f32, total_f64| {
            try std.testing.expectApproxEqRel(total_f64, total_f32, 1e-3);
        }
    }

    // Lee-Richards at f32 is a different algorithm from Shrake-Rupley at f32
    // (before the fix the f32 path always ran Shrake-Rupley).
    const sr_f32 = totals[0][0];
    const lr_f32 = totals[combos.len - 1][0];
    for (sr_f32, lr_f32) |sr_total, lr_total| {
        try std.testing.expect(@abs(sr_total - lr_total) > 5.0);
    }
}

test "traj run: --no-hydrogens equals a run on a hydrogen-free topology and trajectory" {
    var ws = TestWorkspace.init();
    defer ws.deinit();
    const a = ws.allocator();

    // Split the first model into all atoms and heavy atoms using the PDB
    // element column, independently of the parser under test.
    const full_pdb = try ws.firstModelPdb();
    var heavy_pdb: std.ArrayListUnmanaged(u8) = .empty;
    var heavy_indices: std.ArrayListUnmanaged(u32) = .empty;
    var n_atoms: u32 = 0;
    var lines = std.mem.splitScalar(u8, full_pdb, '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "ATOM")) continue;
        defer n_atoms += 1;
        if (std.mem.eql(u8, std.mem.trim(u8, line[76..78], " "), "H")) continue;
        try heavy_pdb.appendSlice(a, line);
        try heavy_pdb.append(a, '\n');
        try heavy_indices.append(a, n_atoms);
    }
    try std.testing.expectEqual(@as(u32, 304), n_atoms);
    try std.testing.expectEqual(@as(usize, 154), heavy_indices.items.len);

    const full_topology = try ws.write("full.pdb", full_pdb);
    const heavy_topology = try ws.write("heavy.pdb", heavy_pdb.items);
    try ws.writeTrrPair("full.trr", "heavy.trr", heavy_indices.items, 3, null);
    const full_traj = try ws.path("full.trr");
    const heavy_traj = try ws.path("heavy.trr");

    for ([_]usize{ 1, 2 }) |n_threads| {
        const stripped = try ws.runTraj(.{
            .traj_path = full_traj,
            .topology_path = full_topology,
            .include_hydrogens = false,
            .n_points = 64,
            .n_threads = n_threads,
        }, "stripped.csv");
        const reference = try ws.runTraj(.{
            .traj_path = heavy_traj,
            .topology_path = heavy_topology,
            .n_points = 64,
            .n_threads = n_threads,
        }, "reference.csv");
        const with_hydrogens = try ws.runTraj(.{
            .traj_path = full_traj,
            .topology_path = full_topology,
            .n_points = 64,
            .n_threads = n_threads,
        }, "with_hydrogens.csv");

        try std.testing.expectEqual(@as(usize, 3), (try ws.totals(stripped)).len);
        try std.testing.expectEqualStrings(reference, stripped);
        try std.testing.expect(!std.mem.eql(u8, with_hydrogens, stripped));
    }

    // The trajectory is validated against the whole topology: a hydrogen-free
    // trajectory does not match a topology that still lists its hydrogens.
    try std.testing.expectError(error.AtomCountMismatch, ws.runTraj(.{
        .traj_path = heavy_traj,
        .topology_path = full_topology,
        .include_hydrogens = false,
        .n_points = 64,
        .n_threads = 1,
    }, "mismatch.csv"));
}

test "traj --altloc applies to PDB topologies as it does to mmCIF topologies" {
    var ws = TestWorkspace.init();
    defer ws.deinit();

    const pdb_path = try ws.write("altloc.pdb",
        \\ATOM      1  N   ALA A   1       1.000   0.000   0.000  1.00 10.00           N
        \\ATOM      2  CA AALA A   1       2.000   0.000   0.000  0.30 10.00           C
        \\ATOM      3  CA BALA A   1       3.000   0.000   0.000  0.70 10.00           C
        \\ATOM      4  C   ALA A   1       4.000   0.000   0.000  1.00 10.00           C
        \\END
        \\
    );
    const cif_path = try ws.write("altloc.cif",
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
    );

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
            const args = parseArgs(&.{ "zsasa", "traj", case.flag, "traj.xtc", path }, 2);
            var topology = try loadTopology(std.testing.allocator, std.testing.io, path, args);
            defer topology.deinit();
            try std.testing.expectEqualSlices(f64, case.x, topology.atoms.x);
            try std.testing.expectEqual(case.x.len, topology.n_file_atoms);
        }

        // `none` exists to fail fast
        const none_args = parseArgs(&.{ "zsasa", "traj", "--altloc=none", "traj.xtc", path }, 2);
        try std.testing.expectError(
            error.UnexpectedAltLoc,
            loadTopology(std.testing.allocator, std.testing.io, path, none_args),
        );
    }
}

test "traj run: multi-model, HETATM and mmCIF topologies map onto the trajectory" {
    var ws = TestWorkspace.init();
    defer ws.deinit();
    const a = ws.allocator();

    const first_model = try ws.firstModelPdb();
    const base = TrajArgs{ .traj_path = test_xtc_path, .n_points = 64, .end_frame = 1, .n_threads = 1 };

    var reference_args = base;
    reference_args.topology_path = try ws.write("m1.pdb", first_model);
    const reference = try ws.runTraj(reference_args, "reference.csv");
    try std.testing.expectEqual(@as(usize, 2), (try ws.totals(reference)).len);

    // The fixture itself: 38 models of 304 atoms. Only the first is the topology.
    var multi_model_args = base;
    multi_model_args.topology_path = test_pdb_path;
    try std.testing.expectEqualStrings(reference, try ws.runTraj(multi_model_args, "multi_model.csv"));

    // The same atoms with the last residue written as HETATM records, as PDB
    // and as a two-model mmCIF.
    var hetatm_pdb: std.ArrayListUnmanaged(u8) = .empty;
    var cif_rows: [2]std.ArrayListUnmanaged(u8) = .{ .empty, .empty };
    var n_hetatm: usize = 0;
    var serial: usize = 0;
    var lines = std.mem.splitScalar(u8, first_model, '\n');
    while (lines.next()) |line| {
        if (!std.mem.startsWith(u8, line, "ATOM")) {
            try hetatm_pdb.appendSlice(a, line);
            try hetatm_pdb.append(a, '\n');
            continue;
        }
        serial += 1;
        const is_last_residue = std.mem.eql(u8, std.mem.trim(u8, line[22..26], " "), "20");
        if (is_last_residue) n_hetatm += 1;
        const group: []const u8 = if (is_last_residue) "HETATM" else "ATOM";

        try hetatm_pdb.appendSlice(a, if (is_last_residue) "HETATM" else "ATOM  ");
        try hetatm_pdb.appendSlice(a, line[6..]);
        try hetatm_pdb.append(a, '\n');

        for (&cif_rows, 1..) |*rows, model| {
            // group id type_symbol atom comp asym seq x y z model
            const row = try std.fmt.allocPrint(a, "{s} {d} {s} {s} {s} {s} {s} {s} {s} {s} {d}\n", .{
                group,
                serial,
                std.mem.trim(u8, line[76..78], " "),
                std.mem.trim(u8, line[12..16], " "),
                std.mem.trim(u8, line[17..20], " "),
                line[21..22],
                std.mem.trim(u8, line[22..26], " "),
                std.mem.trim(u8, line[30..38], " "),
                std.mem.trim(u8, line[38..46], " "),
                std.mem.trim(u8, line[46..54], " "),
                model + 6, // models 7 and 8: the first model is not numbered 1
            });
            try rows.appendSlice(a, row);
        }
    }
    try std.testing.expectEqual(@as(usize, 304), serial);
    try std.testing.expectEqual(@as(usize, 12), n_hetatm);

    var hetatm_args = base;
    hetatm_args.topology_path = try ws.write("hetatm.pdb", hetatm_pdb.items);
    try std.testing.expectEqualStrings(reference, try ws.runTraj(hetatm_args, "hetatm.csv"));

    const cif_header =
        \\data_1L2Y
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
        \\_atom_site.pdbx_PDB_model_num
        \\
    ;
    const cif = try std.mem.concat(a, u8, &.{ cif_header, cif_rows[0].items, cif_rows[1].items, "#\n" });
    var cif_args = base;
    cif_args.topology_path = try ws.write("two_models.cif", cif);
    try std.testing.expectEqualStrings(reference, try ws.runTraj(cif_args, "cif.csv"));
}

test "traj run: invalid options are rejected before the output file is created" {
    var ws = TestWorkspace.init();
    defer ws.deinit();

    const topology = try ws.write("m1.pdb", try ws.firstModelPdb());
    const sentinel = "results from an earlier run\n";
    const output = try ws.write("keep.csv", sentinel);
    const base = TrajArgs{ .traj_path = test_xtc_path, .topology_path = topology, .output_path = output, .quiet = true, .n_threads = 1 };

    const Case = struct { args: TrajArgs, expected: anyerror };
    var cases = [_]Case{
        .{ .args = base, .expected = error.UnsupportedNPoints },
        .{ .args = base, .expected = error.BitmaskRequiresSR },
        .{ .args = base, .expected = error.InvalidArgument },
        .{ .args = base, .expected = error.InvalidArgument },
        .{ .args = base, .expected = error.InvalidArgument },
        .{ .args = base, .expected = error.InvalidArgument },
        .{ .args = base, .expected = error.InvalidArgument },
    };
    cases[0].args.use_bitmask = true;
    cases[0].args.n_points = 2000;
    cases[1].args.use_bitmask = true;
    cases[1].args.algorithm = .lr;
    cases[2].args.probe_radius = 1e6;
    cases[3].args.n_points = 0;
    cases[4].args.n_slices = 0;
    cases[5].args.stride = 0;
    cases[6].args.bitmask_correction = true;

    for (cases) |case| {
        try std.testing.expectError(case.expected, run(std.testing.allocator, std.testing.io, case.args));
        try std.testing.expectEqualStrings(sentinel, try ws.read(output));
    }
}

test "traj run: frames computed before a failure are written to the output" {
    var ws = TestWorkspace.init();
    defer ws.deinit();

    const topology = try ws.write("m1.pdb", try ws.firstModelPdb());
    const base = TrajArgs{ .topology_path = topology, .n_points = 64, .batch_size = 4 };

    // Reference rows from the intact trajectory.
    var full_args = base;
    full_args.traj_path = test_xtc_path;
    full_args.end_frame = 5;
    full_args.n_threads = 1;
    const full = try ws.runTraj(full_args, "full.csv");

    // A trajectory cut in the middle of a frame: the reader fails after the
    // complete frames.
    const xtc_bytes = try ws.read(test_xtc_path);
    const truncated_traj = try ws.write("truncated.xtc", xtc_bytes[0..6000]);
    // One coordinate of frame 2 is NaN: the SASA calculation of that frame fails.
    try ws.writeTrrPair("nan.trr", "unused.trr", &.{0}, 4, 2);
    var clean_args = base;
    clean_args.traj_path = try ws.path("nan.trr");
    clean_args.end_frame = 1;
    clean_args.n_threads = 1;
    const clean_rows = try ws.runTraj(clean_args, "clean.csv");
    try std.testing.expectEqual(@as(usize, 2), (try ws.totals(clean_rows)).len);

    for ([_]usize{ 1, 2 }) |n_threads| {
        var truncated_args = base;
        truncated_args.traj_path = truncated_traj;
        truncated_args.n_threads = n_threads;
        truncated_args.output_path = try ws.path("truncated.csv");
        truncated_args.quiet = true;
        if (run(std.testing.allocator, std.testing.io, truncated_args)) |_| {
            return error.TestExpectedError;
        } else |_| {}
        const truncated_rows = try ws.read(truncated_args.output_path);
        // 6000 bytes hold three complete frames of about 1650 bytes each.
        try std.testing.expectEqual(@as(usize, 3), (try ws.totals(truncated_rows)).len);
        try std.testing.expect(std.mem.startsWith(u8, full, truncated_rows));

        var nan_args = base;
        nan_args.traj_path = try ws.path("nan.trr");
        nan_args.n_threads = n_threads;
        nan_args.output_path = try ws.path("nan.csv");
        nan_args.quiet = true;
        if (run(std.testing.allocator, std.testing.io, nan_args)) |_| {
            return error.TestExpectedError;
        } else |_| {}
        // Frames 0 and 1 precede the failing frame and are kept.
        try std.testing.expectEqualStrings(clean_rows, try ws.read(nan_args.output_path));
    }
}

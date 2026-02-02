const std = @import("std");
const build_options = @import("build_options");
const analysis = @import("analysis.zig");
const batch = @import("batch.zig");
const traj = @import("traj.zig");
const format_detect = @import("format_detect.zig");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const pdb_parser = @import("pdb_parser.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const classifier = @import("classifier.zig");
const classifier_parser = @import("classifier_parser.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_protor = @import("classifier_protor.zig");
const classifier_oons = @import("classifier_oons.zig");

const Config = types.Config;
const Configf32 = types.Configf32;
const Precision = types.Precision;
const OutputFormat = json_writer.OutputFormat;
const LeeRichardsConfig = lee_richards.LeeRichardsConfig;
const LeeRichardsConfigf32 = lee_richards.LeeRichardsConfigf32;
const ClassifierType = classifier.ClassifierType;
const Parallelism = batch.Parallelism;

const version = build_options.version;

/// SASA algorithm selection
const Algorithm = enum {
    sr, // Shrake-Rupley (test point method)
    lr, // Lee-Richards (slice method)
};

/// Parsed command-line arguments
const Args = struct {
    input_path: ?[]const u8 = null,
    output_path: []const u8 = "output.json",
    output_path_explicit: bool = false, // Track if -o was explicitly set
    n_threads: usize = 0, // 0 = auto-detect
    probe_radius: f64 = 1.4,
    n_points: u32 = 100, // For Shrake-Rupley
    n_slices: u32 = 20, // For Lee-Richards
    algorithm: Algorithm = .sr, // Default: Shrake-Rupley
    precision: Precision = .f64, // f32 or f64 (default: f64)
    parallelism: Parallelism = .file, // file, atom, or pipeline (batch mode only)
    output_format: OutputFormat = .json,
    classifier_type: ?ClassifierType = null, // Built-in classifier (naccess/protor/oons)
    config_path: ?[]const u8 = null, // Custom config file path
    chain_filter: ?[]const u8 = null, // Chain filter (e.g., "A" or "A,B,C")
    model_num: ?u32 = null, // Model number for NMR structures
    use_auth_chain: bool = false, // Use auth_asym_id instead of label_asym_id
    per_residue: bool = false, // Show per-residue SASA
    rsa: bool = false, // Show RSA (Relative Solvent Accessibility)
    polar: bool = false, // Show polar/nonpolar SASA summary
    quiet: bool = false,
    validate_only: bool = false,
    show_timing: bool = false, // Show timing breakdown for benchmarking
    show_help: bool = false,
    show_version: bool = false,
};

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

/// Parse and validate output format value
fn parseOutputFormat(value: []const u8) OutputFormat {
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

/// Parse and validate classifier type value
fn parseClassifierType(value: []const u8) ClassifierType {
    if (ClassifierType.fromString(value)) |ct| {
        return ct;
    } else {
        std.debug.print("Error: Invalid classifier: {s}\n", .{value});
        std.debug.print("Valid classifiers: naccess, protor, oons\n", .{});
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

/// Parse and validate parallelism strategy value (batch mode only)
fn parseParallelism(value: []const u8) Parallelism {
    if (std.mem.eql(u8, value, "file")) {
        return .file;
    } else if (std.mem.eql(u8, value, "atom")) {
        return .atom;
    } else if (std.mem.eql(u8, value, "pipeline")) {
        return .pipeline;
    } else {
        std.debug.print("Error: Invalid parallelism strategy: {s}\n", .{value});
        std.debug.print("Valid values: file, atom, pipeline\n", .{});
        std.process.exit(1);
    }
}

/// Parse chain filter string into array of chain IDs
fn parseChainFilter(allocator: std.mem.Allocator, filter_str: []const u8) ![]const []const u8 {
    var chains = std.ArrayListUnmanaged([]const u8){};
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

/// Read input file (auto-detect format)
fn readInputFile(allocator: std.mem.Allocator, path: []const u8, args: Args) !types.AtomInput {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => json_parser.readAtomInputFromFile(allocator, path),
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.model_num = args.model_num;
            parser.use_auth_chain = args.use_auth_chain;

            // Parse chain filter if specified
            var chain_filter_slice: ?[]const []const u8 = null;
            if (args.chain_filter) |filter_str| {
                chain_filter_slice = try parseChainFilter(allocator, filter_str);
                parser.chain_filter = chain_filter_slice;
            }
            defer if (chain_filter_slice) |s| allocator.free(s);

            break :blk parser.parseFile(path);
        },
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.model_num = args.model_num;

            // Parse chain filter if specified
            var chain_filter_slice: ?[]const []const u8 = null;
            if (args.chain_filter) |filter_str| {
                chain_filter_slice = try parseChainFilter(allocator, filter_str);
                parser.chain_filter = chain_filter_slice;
            }
            defer if (chain_filter_slice) |s| allocator.free(s);

            break :blk parser.parseFile(path);
        },
    };
}

/// Parse command-line arguments
fn parseArgs(args: []const []const u8) Args {
    var result = Args{};
    var i: usize = 1;

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
        // --format=FORMAT or --format FORMAT
        else if (std.mem.startsWith(u8, arg, "--format=")) {
            const value = arg["--format=".len..];
            result.output_format = parseOutputFormat(value);
        } else if (std.mem.eql(u8, arg, "--format")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --format\n", .{});
                std.process.exit(1);
            }
            result.output_format = parseOutputFormat(args[i]);
        }
        // --algorithm=ALGO or --algorithm ALGO
        else if (std.mem.startsWith(u8, arg, "--algorithm=")) {
            const value = arg["--algorithm=".len..];
            result.algorithm = parseAlgorithm(value);
        } else if (std.mem.eql(u8, arg, "--algorithm")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --algorithm\n", .{});
                std.process.exit(1);
            }
            result.algorithm = parseAlgorithm(args[i]);
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
        // --classifier=TYPE or --classifier TYPE
        else if (std.mem.startsWith(u8, arg, "--classifier=")) {
            const value = arg["--classifier=".len..];
            result.classifier_type = parseClassifierType(value);
        } else if (std.mem.eql(u8, arg, "--classifier")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --classifier\n", .{});
                std.process.exit(1);
            }
            result.classifier_type = parseClassifierType(args[i]);
        }
        // --precision=PREC or --precision PREC
        else if (std.mem.startsWith(u8, arg, "--precision=")) {
            const value = arg["--precision=".len..];
            result.precision = parsePrecision(value);
        } else if (std.mem.eql(u8, arg, "--precision")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --precision\n", .{});
                std.process.exit(1);
            }
            result.precision = parsePrecision(args[i]);
        }
        // --parallelism=MODE or --parallelism MODE (batch mode only)
        else if (std.mem.startsWith(u8, arg, "--parallelism=")) {
            const value = arg["--parallelism=".len..];
            result.parallelism = parseParallelism(value);
        } else if (std.mem.eql(u8, arg, "--parallelism")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --parallelism\n", .{});
                std.process.exit(1);
            }
            result.parallelism = parseParallelism(args[i]);
        }
        // --config=FILE or --config FILE
        else if (std.mem.startsWith(u8, arg, "--config=")) {
            const value = arg["--config=".len..];
            result.config_path = value;
        } else if (std.mem.eql(u8, arg, "--config")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --config\n", .{});
                std.process.exit(1);
            }
            result.config_path = args[i];
        }
        // --chain=ID or --chain ID (e.g., --chain=A or --chain=A,B,C)
        else if (std.mem.startsWith(u8, arg, "--chain=")) {
            const value = arg["--chain=".len..];
            result.chain_filter = value;
        } else if (std.mem.eql(u8, arg, "--chain")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --chain\n", .{});
                std.process.exit(1);
            }
            result.chain_filter = args[i];
        }
        // --model=N or --model N
        else if (std.mem.startsWith(u8, arg, "--model=")) {
            const value = arg["--model=".len..];
            const model = std.fmt.parseInt(u32, value, 10) catch {
                std.debug.print("Error: Invalid model number: {s}\n", .{value});
                std.process.exit(1);
            };
            if (model == 0) {
                std.debug.print("Error: Model number must be >= 1\n", .{});
                std.process.exit(1);
            }
            result.model_num = model;
        } else if (std.mem.eql(u8, arg, "--model")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --model\n", .{});
                std.process.exit(1);
            }
            const model = std.fmt.parseInt(u32, args[i], 10) catch {
                std.debug.print("Error: Invalid model number: {s}\n", .{args[i]});
                std.process.exit(1);
            };
            if (model == 0) {
                std.debug.print("Error: Model number must be >= 1\n", .{});
                std.process.exit(1);
            }
            result.model_num = model;
        }
        // --auth-chain (use auth_asym_id instead of label_asym_id)
        else if (std.mem.eql(u8, arg, "--auth-chain")) {
            result.use_auth_chain = true;
        }
        // --per-residue (show per-residue SASA)
        else if (std.mem.eql(u8, arg, "--per-residue")) {
            result.per_residue = true;
        }
        // --rsa (show RSA - Relative Solvent Accessibility)
        else if (std.mem.eql(u8, arg, "--rsa")) {
            result.rsa = true;
            result.per_residue = true; // RSA implies per-residue
        }
        // --polar (show polar/nonpolar SASA summary)
        else if (std.mem.eql(u8, arg, "--polar")) {
            result.polar = true;
            result.per_residue = true; // Polar analysis requires per-residue
        }
        // --quiet or -q
        else if (std.mem.eql(u8, arg, "--quiet") or std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
        }
        // --validate
        else if (std.mem.eql(u8, arg, "--validate")) {
            result.validate_only = true;
        }
        // --timing
        else if (std.mem.eql(u8, arg, "--timing")) {
            result.show_timing = true;
        }
        // --help or -h
        else if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            result.show_help = true;
        }
        // --version or -V
        else if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")) {
            result.show_version = true;
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
            std.debug.print("Try '--help' for more information.\n", .{});
            std.process.exit(1);
        }
        // Positional arguments
        else if (result.input_path == null) {
            result.input_path = arg;
        } else if (!result.output_path_explicit) {
            // Only use positional output if -o was not explicitly set
            result.output_path = arg;
            result.output_path_explicit = true;
        }
    }

    return result;
}

fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\zsasa {s} - Solvent Accessible Surface Area calculator
        \\
        \\USAGE:
        \\    {s} [OPTIONS] <input> [output.json]
        \\
        \\ARGUMENTS:
        \\    <input>          Input file (JSON or mmCIF format, auto-detected by extension)
        \\                     Supported: .json, .cif, .mmcif, .pdb, .ent
        \\    [output.json]    Output JSON file (default: output.json)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: naccess, protor, oons
        \\                       Use with residue/atom_name input for auto-radius
        \\    --config=FILE      Custom classifier config file (FreeSASA format)
        \\    --chain=ID         Filter by chain ID (e.g., --chain=A or --chain=A,B,C)
        \\                       Default: label_asym_id (mmCIF standard)
        \\    --auth-chain       Use auth_asym_id instead of label_asym_id
        \\    --model=N          Model number for NMR structures (default: all)
        \\    --per-residue      Show per-residue SASA aggregation
        \\    --rsa              Show RSA (Relative Solvent Accessibility)
        \\    --polar            Show polar/nonpolar SASA summary
        \\    --threads=N        Number of threads (default: auto-detect)
        \\                       Use --threads=1 for single-threaded mode
        \\    --probe-radius=R   Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N       Test points per atom (default: 100, for sr)
        \\    --n-slices=N       Slices per atom diameter (default: 20, for lr)
        \\    --format=FORMAT    Output format: json, compact, csv (default: json)
        \\    --precision=PREC   Floating-point precision: f32, f64 (default: f64)
        \\                       f32 is faster but less accurate
        \\    --parallelism=MODE Batch mode parallelism: file, atom, pipeline (default: file)
        \\                       file: N files in parallel, 1 thread per file
        \\                       atom: 1 file at a time, N threads for SASA
        \\                       pipeline: I/O prefetch + atom-level SASA
        \\    -o, --output=FILE  Output file (alternative to positional argument)
        \\    --validate         Validate input only, do not calculate SASA
        \\    --timing           Show timing breakdown (for benchmarking)
        \\    -q, --quiet        Suppress progress output
        \\    -h, --help         Show this help message
        \\    -V, --version      Show version
        \\
        \\ALGORITHMS:
        \\    sr, shrake-rupley  Test point method (default)
        \\    lr, lee-richards   Slice-based method
        \\
        \\CLASSIFIERS:
        \\    naccess  NACCESS-compatible radii (default if --classifier used)
        \\    protor   ProtOr radii (Tsai et al. 1999)
        \\    oons     OONS radii (Ooi et al.)
        \\
        \\OUTPUT FORMATS:
        \\    json     Pretty-printed JSON with indentation
        \\    compact  Single-line JSON (no whitespace)
        \\    csv      CSV with atom_index,area columns
        \\
        \\EXAMPLES:
        \\    {s} input.json output.json
        \\    {s} structure.cif output.json              # mmCIF input
        \\    {s} structure.pdb output.json              # PDB input
        \\    {s} --algorithm=lr input.json output.json
        \\    {s} --algorithm=lr --n-slices=50 input.json output.json
        \\    {s} --threads=4 input.json output.json
        \\    {s} --probe-radius=1.5 --n-points=200 input.json
        \\    {s} --format=csv input.json output.csv
        \\    {s} --classifier=naccess input.json output.json
        \\    {s} --classifier=naccess structure.cif output.json
        \\    {s} --config=custom.config input.json output.json
        \\
        \\BATCH MODE (directory input):
        \\    {s} input_dir/ output_dir/
        \\    {s} --threads=8 input_dir/ output_dir/           # file-level parallelism
        \\    {s} --parallelism=atom --threads=8 input_dir/    # atom-level parallelism
        \\
    , .{ version, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
}

fn printVersion() void {
    std.debug.print("zsasa {s}\n", .{version});
}

/// Apply a custom classifier to input, replacing radii based on residue/atom names
fn applyClassifier(
    _: std.mem.Allocator,
    input: *types.AtomInput,
    custom_classifier: *const classifier.Classifier,
    quiet: bool,
) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    // Allocate new radii array (use input.allocator for consistency with deinit)
    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    var classified_count: usize = 0;
    var fallback_count: usize = 0;

    for (0..n) |i| {
        // Try classifier lookup
        if (custom_classifier.getRadius(residues[i].slice(), atom_names[i].slice())) |r| {
            new_radii[i] = r;
            classified_count += 1;
        } else if (input.element) |elements| {
            // Fall back to element-based radius
            if (classifier.guessRadiusFromAtomicNumber(elements[i])) |r| {
                new_radii[i] = r;
                fallback_count += 1;
            } else {
                // Keep original radius
                new_radii[i] = input.r[i];
            }
        } else if (classifier.guessRadiusFromAtomName(atom_names[i].slice())) |r| {
            // Fall back to atom name-based radius
            new_radii[i] = r;
            fallback_count += 1;
        } else {
            // Keep original radius
            new_radii[i] = input.r[i];
        }
    }

    // Free old radii and replace with classified radii
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

/// Apply a built-in classifier to input
fn applyBuiltinClassifier(
    _: std.mem.Allocator,
    input: *types.AtomInput,
    ct: ClassifierType,
    quiet: bool,
) !void {
    const n = input.atomCount();
    const residues = input.residue orelse return error.MissingClassificationInfo;
    const atom_names = input.atom_name orelse return error.MissingClassificationInfo;

    // Allocate new radii array (use input.allocator for consistency with deinit)
    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    var classified_count: usize = 0;
    var fallback_count: usize = 0;

    for (0..n) |i| {
        // Try built-in classifier lookup
        const maybe_radius: ?f64 = switch (ct) {
            .naccess => classifier_naccess.getRadius(residues[i].slice(), atom_names[i].slice()),
            .protor => classifier_protor.getRadius(residues[i].slice(), atom_names[i].slice()),
            .oons => classifier_oons.getRadius(residues[i].slice(), atom_names[i].slice()),
        };

        if (maybe_radius) |r| {
            new_radii[i] = r;
            classified_count += 1;
        } else if (input.element) |elements| {
            // Fall back to element-based radius
            if (classifier.guessRadiusFromAtomicNumber(elements[i])) |r| {
                new_radii[i] = r;
                fallback_count += 1;
            } else {
                // Keep original radius
                new_radii[i] = input.r[i];
            }
        } else if (classifier.guessRadiusFromAtomName(atom_names[i].slice())) |r| {
            // Fall back to atom name-based radius
            new_radii[i] = r;
            fallback_count += 1;
        } else {
            // Keep original radius
            new_radii[i] = input.r[i];
        }
    }

    // Free old radii and replace with classified radii
    input.allocator.free(input.r);
    input.r = new_radii;

    if (!quiet) {
        std.debug.print("Classifier '{s}': {d} atoms classified, {d} fallback\n", .{
            ct.name(),
            classified_count,
            fallback_count,
        });
    }
}

/// Print per-chain SASA results
fn printPerChainResults(chain_ids: []const types.FixedString4, atom_areas: []const f64) void {
    // Use a simple approach: iterate through to find unique chains and sum areas
    // For efficiency, we'll use a fixed-size buffer for up to 64 unique chains
    const max_chains = 64;
    var chain_names: [max_chains]types.FixedString4 = undefined;
    var chain_areas: [max_chains]f64 = undefined;
    var chain_counts: [max_chains]usize = undefined;
    var num_chains: usize = 0;
    var warned_overflow = false;

    for (chain_ids, 0..) |chain_id, i| {
        // Find if this chain already exists
        var found_idx: ?usize = null;
        for (0..num_chains) |j| {
            if (std.mem.eql(u8, chain_names[j].slice(), chain_id.slice())) {
                found_idx = j;
                break;
            }
        }

        if (found_idx) |idx| {
            // Add to existing chain
            chain_areas[idx] += atom_areas[i];
            chain_counts[idx] += 1;
        } else if (num_chains < max_chains) {
            // Add new chain
            chain_names[num_chains] = chain_id;
            chain_areas[num_chains] = atom_areas[i];
            chain_counts[num_chains] = 1;
            num_chains += 1;
        } else if (!warned_overflow) {
            // Warn once when limit is exceeded
            std.debug.print("Warning: More than {d} unique chains; some chains omitted from summary\n", .{max_chains});
            warned_overflow = true;
        }
    }

    if (num_chains > 0) {
        std.debug.print("\nPer-chain SASA:\n", .{});
        for (0..num_chains) |i| {
            std.debug.print("  {s}: {d:.2} Å² ({d} atoms)\n", .{
                chain_names[i].slice(),
                chain_areas[i],
                chain_counts[i],
            });
        }
    }
}

/// Run batch processing mode for directory input
fn runBatchMode(allocator: std.mem.Allocator, parsed: Args) !void {
    const input_dir = parsed.input_path.?;
    const output_dir: ?[]const u8 = if (parsed.output_path_explicit)
        parsed.output_path
    else
        null;

    // Build batch config from parsed args
    const config = batch.BatchConfig{
        .n_threads = parsed.n_threads,
        .algorithm = switch (parsed.algorithm) {
            .sr => .sr,
            .lr => .lr,
        },
        .parallelism = parsed.parallelism,
        .n_points = parsed.n_points,
        .n_slices = parsed.n_slices,
        .probe_radius = parsed.probe_radius,
        .precision = parsed.precision,
        .output_format = parsed.output_format,
        .show_timing = parsed.show_timing,
        .quiet = parsed.quiet,
    };

    if (!parsed.quiet) {
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

    var result = try batch.runBatch(allocator, input_dir, output_dir, config);
    defer result.deinit();

    // Print results
    if (!parsed.quiet) {
        result.printSummary(parsed.show_timing);
    }

    // Always print benchmark output (for script parsing)
    if (parsed.show_timing) {
        std.debug.print("\n", .{});
        result.printBenchmarkOutput();
    }
}

pub fn main() !void {
    // Setup allocator
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse command line arguments
    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        printHelp(args[0]);
        std.process.exit(1);
    }

    // Check for subcommands
    if (std.mem.eql(u8, args[1], "traj")) {
        // Trajectory mode
        const traj_args = traj.parseArgs(args, 2);
        if (traj_args.show_help) {
            traj.printHelp(args[0]);
            return;
        }
        traj.run(allocator, traj_args) catch |err| {
            std.debug.print("Error in trajectory mode: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
        return;
    }

    const parsed = parseArgs(args);

    // Handle --help
    if (parsed.show_help) {
        printHelp(args[0]);
        return;
    }

    // Handle --version
    if (parsed.show_version) {
        printVersion();
        return;
    }

    // Validate input path
    if (parsed.input_path == null) {
        std.debug.print("Error: Missing input file\n", .{});
        printHelp(args[0]);
        std.process.exit(1);
    }

    // Check if input is a directory (batch mode)
    const input_stat = std.fs.cwd().statFile(parsed.input_path.?) catch |err| {
        std.debug.print("Error: Cannot access '{s}': {s}\n", .{ parsed.input_path.?, @errorName(err) });
        std.process.exit(1);
    };

    if (input_stat.kind == .directory) {
        // Batch mode: process directory
        runBatchMode(allocator, parsed) catch |err| {
            std.debug.print("Error in batch processing: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
        return;
    }

    // Single file mode continues below

    // Timing variables (in nanoseconds)
    var timer = std.time.Timer.start() catch {
        std.debug.print("Error: Timer not supported\n", .{});
        std.process.exit(1);
    };
    var time_parse: u64 = 0;
    var time_classify: u64 = 0;
    var time_sasa: u64 = 0;
    var time_write: u64 = 0;

    // Read input file (JSON or mmCIF)
    timer.reset();
    var input = readInputFile(allocator, parsed.input_path.?, parsed) catch |err| {
        std.debug.print("Error reading input file '{s}': {s}\n", .{ parsed.input_path.?, @errorName(err) });
        std.process.exit(1);
    };
    defer input.deinit();

    // Validate input data
    var validation = json_parser.validateInput(allocator, input) catch |err| {
        std.debug.print("Error validating input: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer validation.deinit();

    if (!validation.valid) {
        json_parser.printValidationErrors(validation.errors);
        std.process.exit(1);
    }

    // Check for duplicate coordinates (warning only)
    if (!parsed.quiet) {
        _ = json_parser.checkDuplicateCoordinates(allocator, input) catch |err| {
            std.debug.print("Warning: Could not check for duplicate coordinates: {s}\n", .{@errorName(err)});
        };
    }
    time_parse = timer.read();

    // Handle --validate (dry-run)
    if (parsed.validate_only) {
        if (!parsed.quiet) {
            std.debug.print("Input validation passed: {} atoms\n", .{input.atomCount()});
        }
        return;
    }

    // Apply classifier if specified (--config takes precedence over --classifier)
    timer.reset();
    if (parsed.config_path != null or parsed.classifier_type != null) {
        // Warn if both are specified
        if (parsed.config_path != null and parsed.classifier_type != null) {
            if (!parsed.quiet) {
                std.debug.print("Warning: Both --classifier and --config specified; using --config\n", .{});
            }
        }

        // Check if input has classification info
        if (!input.hasClassificationInfo()) {
            std.debug.print("Error: Classifier requires 'residue' and 'atom_name' fields in input JSON\n", .{});
            std.process.exit(1);
        }

        // Load classifier and apply radii
        if (parsed.config_path) |config_path| {
            // Load from custom config file
            var custom_classifier = classifier_parser.parseConfigFile(allocator, config_path) catch |err| {
                std.debug.print("Error loading config file '{s}': {s}\n", .{ config_path, @errorName(err) });
                std.process.exit(1);
            };
            defer custom_classifier.deinit();

            applyClassifier(allocator, &input, &custom_classifier, parsed.quiet) catch |err| {
                std.debug.print("Error applying classifier: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        } else if (parsed.classifier_type) |ct| {
            // Use built-in classifier
            applyBuiltinClassifier(allocator, &input, ct, parsed.quiet) catch |err| {
                std.debug.print("Error applying classifier: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        }
    }
    time_classify = timer.read();

    // Calculate SASA with configured parameters
    timer.reset();
    var result = switch (parsed.precision) {
        .f64 => switch (parsed.algorithm) {
            .sr => blk: {
                const config = Config{
                    .n_points = parsed.n_points,
                    .probe_radius = parsed.probe_radius,
                };
                break :blk if (parsed.n_threads == 1)
                    shrake_rupley.calculateSasa(allocator, input, config) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    }
                else
                    shrake_rupley.calculateSasaParallel(allocator, input, config, parsed.n_threads) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    };
            },
            .lr => blk: {
                const lr_config = LeeRichardsConfig{
                    .n_slices = parsed.n_slices,
                    .probe_radius = parsed.probe_radius,
                };
                break :blk if (parsed.n_threads == 1)
                    lee_richards.calculateSasa(allocator, input, lr_config) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    }
                else
                    lee_richards.calculateSasaParallel(allocator, input, lr_config, parsed.n_threads) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    };
            },
        },
        .f32 => blk: {
            // Calculate in f32, convert to f64 for output
            var result_f32 = switch (parsed.algorithm) {
                .sr => inner: {
                    const config = Configf32{
                        .n_points = parsed.n_points,
                        .probe_radius = @floatCast(parsed.probe_radius),
                    };
                    break :inner if (parsed.n_threads == 1)
                        shrake_rupley.calculateSasaf32(allocator, input, config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        shrake_rupley.calculateSasaParallelf32(allocator, input, config, parsed.n_threads) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        };
                },
                .lr => inner: {
                    const lr_config = LeeRichardsConfigf32{
                        .n_slices = parsed.n_slices,
                        .probe_radius = @floatCast(parsed.probe_radius),
                    };
                    break :inner if (parsed.n_threads == 1)
                        lee_richards.calculateSasaf32(allocator, input, lr_config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        lee_richards.calculateSasaParallelf32(allocator, input, lr_config, parsed.n_threads) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        };
                },
            };
            defer result_f32.deinit();
            // Convert to f64 for output
            break :blk result_f32.toF64(allocator) catch |err| {
                std.debug.print("Error converting result: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        },
    };
    defer result.deinit();
    time_sasa = timer.read();

    // Write output file
    timer.reset();
    json_writer.writeSasaResultWithFormatAndInput(allocator, result, input, parsed.output_path, parsed.output_format) catch |err| {
        std.debug.print("Error writing output file '{s}': {s}\n", .{ parsed.output_path, @errorName(err) });
        std.process.exit(1);
    };
    time_write = timer.read();

    // Calculate total time
    const time_total = time_parse + time_classify + time_sasa + time_write;

    // Print summary (unless quiet mode)
    if (!parsed.quiet) {
        std.debug.print("Calculated SASA for {} atoms\n", .{input.atomCount()});
        std.debug.print("Total area: {d:.2} Å²\n", .{result.total_area});

        // Print per-chain results if chain info is available
        if (input.chain_id) |chain_ids| {
            printPerChainResults(chain_ids, result.atom_areas);
        }

        // Print per-residue results if requested
        if (parsed.per_residue and input.hasResidueInfo()) {
            var residue_result = analysis.aggregateByResidue(allocator, input, result.atom_areas) catch |err| {
                std.debug.print("Error calculating per-residue SASA: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
            defer residue_result.deinit();
            if (parsed.rsa) {
                analysis.printResidueResultsWithRsa(residue_result.residues);
            } else {
                analysis.printResidueResults(residue_result.residues);
            }

            // Print polar/nonpolar summary if requested
            if (parsed.polar) {
                const polar_summary = analysis.calculatePolarSummary(residue_result.residues);
                analysis.printPolarSummary(polar_summary);
            }
        }

        std.debug.print("Output written to {s}\n", .{parsed.output_path});
    }

    // Print timing breakdown if requested
    if (parsed.show_timing) {
        const ns_to_ms = 1_000_000.0;
        std.debug.print("\nTiming breakdown:\n", .{});
        std.debug.print("  Parse + validate: {d:.2} ms\n", .{@as(f64, @floatFromInt(time_parse)) / ns_to_ms});
        if (time_classify > 0) {
            std.debug.print("  Classifier:       {d:.2} ms\n", .{@as(f64, @floatFromInt(time_classify)) / ns_to_ms});
        }
        std.debug.print("  SASA calculation: {d:.2} ms\n", .{@as(f64, @floatFromInt(time_sasa)) / ns_to_ms});
        std.debug.print("  Write output:     {d:.2} ms\n", .{@as(f64, @floatFromInt(time_write)) / ns_to_ms});
        std.debug.print("  Total:            {d:.2} ms\n", .{@as(f64, @floatFromInt(time_total)) / ns_to_ms});
    }
}

// Tests for argument parsing
test "parseArgs defaults" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("output.json", parsed.output_path);
    try std.testing.expectEqual(@as(usize, 0), parsed.n_threads);
    try std.testing.expectEqual(@as(f64, 1.4), parsed.probe_radius);
    try std.testing.expectEqual(@as(u32, 100), parsed.n_points);
    try std.testing.expectEqual(false, parsed.quiet);
    try std.testing.expectEqual(false, parsed.show_help);
    try std.testing.expectEqual(false, parsed.show_version);
}

test "parseArgs with output path" {
    const args = [_][]const u8{ "zsasa", "input.json", "result.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.json", parsed.output_path);
}

test "parseArgs --threads=N" {
    const args = [_][]const u8{ "zsasa", "--threads=4", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs --probe-radius=R" {
    const args = [_][]const u8{ "zsasa", "--probe-radius=1.5", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(f64, 1.5), parsed.probe_radius);
}

test "parseArgs --n-points=N" {
    const args = [_][]const u8{ "zsasa", "--n-points=200", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

test "parseArgs --quiet" {
    const args = [_][]const u8{ "zsasa", "--quiet", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.quiet);
}

test "parseArgs -q" {
    const args = [_][]const u8{ "zsasa", "-q", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.quiet);
}

test "parseArgs --help" {
    const args = [_][]const u8{ "zsasa", "--help" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_help);
}

test "parseArgs -h" {
    const args = [_][]const u8{ "zsasa", "-h" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_help);
}

test "parseArgs --version" {
    const args = [_][]const u8{ "zsasa", "--version" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_version);
}

test "parseArgs -V" {
    const args = [_][]const u8{ "zsasa", "-V" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_version);
}

test "parseArgs multiple options" {
    const args = [_][]const u8{
        "zsasa",
        "--threads=8",
        "--probe-radius=1.6",
        "--n-points=150",
        "--quiet",
        "input.json",
        "output.json",
    };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(usize, 8), parsed.n_threads);
    try std.testing.expectEqual(@as(f64, 1.6), parsed.probe_radius);
    try std.testing.expectEqual(@as(u32, 150), parsed.n_points);
    try std.testing.expectEqual(true, parsed.quiet);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("output.json", parsed.output_path);
}

// Tests for space-separated option syntax
test "parseArgs --threads N (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--threads", "4", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs --probe-radius R (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--probe-radius", "1.5", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(f64, 1.5), parsed.probe_radius);
}

test "parseArgs --n-points N (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--n-points", "200", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

// Tests for --format option
test "parseArgs --format=json" {
    const args = [_][]const u8{ "zsasa", "--format=json", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.json, parsed.output_format);
}

test "parseArgs --format=compact" {
    const args = [_][]const u8{ "zsasa", "--format=compact", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.compact, parsed.output_format);
}

test "parseArgs --format=csv" {
    const args = [_][]const u8{ "zsasa", "--format=csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "parseArgs --format csv (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--format", "csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "parseArgs default format is json" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.json, parsed.output_format);
}

test "parseArgs --validate" {
    const args = [_][]const u8{ "zsasa", "--validate", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.validate_only);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs default validate_only is false" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(false, parsed.validate_only);
}

// Tests for --algorithm option
test "parseArgs default algorithm is sr" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "parseArgs --algorithm=sr" {
    const args = [_][]const u8{ "zsasa", "--algorithm=sr", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "parseArgs --algorithm=lr" {
    const args = [_][]const u8{ "zsasa", "--algorithm=lr", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

test "parseArgs --algorithm=shrake-rupley" {
    const args = [_][]const u8{ "zsasa", "--algorithm=shrake-rupley", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "parseArgs --algorithm=lee-richards" {
    const args = [_][]const u8{ "zsasa", "--algorithm=lee-richards", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

test "parseArgs --algorithm lr (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--algorithm", "lr", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

// Tests for --n-slices option
test "parseArgs default n_slices is 20" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 20), parsed.n_slices);
}

test "parseArgs --n-slices=50" {
    const args = [_][]const u8{ "zsasa", "--n-slices=50", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 50), parsed.n_slices);
}

test "parseArgs --n-slices 100 (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--n-slices", "100", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 100), parsed.n_slices);
}

test "parseArgs combined algorithm and n-slices" {
    const args = [_][]const u8{ "zsasa", "--algorithm=lr", "--n-slices=50", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
    try std.testing.expectEqual(@as(u32, 50), parsed.n_slices);
}

// Tests for --classifier option
test "parseArgs default classifier_type is null" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(?ClassifierType, null), parsed.classifier_type);
}

test "parseArgs --classifier=naccess" {
    const args = [_][]const u8{ "zsasa", "--classifier=naccess", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.naccess, parsed.classifier_type.?);
}

test "parseArgs --classifier=protor" {
    const args = [_][]const u8{ "zsasa", "--classifier=protor", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.protor, parsed.classifier_type.?);
}

test "parseArgs --classifier=oons" {
    const args = [_][]const u8{ "zsasa", "--classifier=oons", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.oons, parsed.classifier_type.?);
}

test "parseArgs --classifier protor (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--classifier", "protor", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.protor, parsed.classifier_type.?);
}

// Tests for --config option
test "parseArgs default config_path is null" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(?[]const u8, null), parsed.config_path);
}

test "parseArgs --config=custom.config" {
    const args = [_][]const u8{ "zsasa", "--config=custom.config", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("custom.config", parsed.config_path.?);
}

test "parseArgs --config custom.config (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--config", "custom.config", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("custom.config", parsed.config_path.?);
}

test "parseArgs combined classifier and config (config takes precedence)" {
    const args = [_][]const u8{ "zsasa", "--classifier=naccess", "--config=custom.config", "input.json" };
    const parsed = parseArgs(&args);

    // Both are set - config should take precedence in main() logic
    try std.testing.expectEqual(ClassifierType.naccess, parsed.classifier_type.?);
    try std.testing.expectEqualStrings("custom.config", parsed.config_path.?);
}

// Tests for --timing option
test "parseArgs default show_timing is false" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(false, parsed.show_timing);
}

test "parseArgs --timing" {
    const args = [_][]const u8{ "zsasa", "--timing", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_timing);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs -o FILE" {
    const args = [_][]const u8{ "zsasa", "-o", "result.csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.csv", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "parseArgs --output=FILE" {
    const args = [_][]const u8{ "zsasa", "--output=result.csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.csv", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "parseArgs --output FILE (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--output", "result.csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.csv", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "parseArgs -o takes precedence over positional output" {
    // -o before positional
    const args1 = [_][]const u8{ "zsasa", "-o", "explicit.json", "input.json", "positional.json" };
    const parsed1 = parseArgs(&args1);
    try std.testing.expectEqualStrings("explicit.json", parsed1.output_path);

    // -o after positional (should still win)
    const args2 = [_][]const u8{ "zsasa", "input.json", "positional.json", "-o", "explicit.json" };
    const parsed2 = parseArgs(&args2);
    try std.testing.expectEqualStrings("explicit.json", parsed2.output_path);
}

// Tests for --precision option
test "parseArgs default precision is f64" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Precision.f64, parsed.precision);
}

test "parseArgs --precision=f32" {
    const args = [_][]const u8{ "zsasa", "--precision=f32", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Precision.f32, parsed.precision);
}

test "parseArgs --precision=f64" {
    const args = [_][]const u8{ "zsasa", "--precision=f64", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Precision.f64, parsed.precision);
}

test "parseArgs --precision f32 (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--precision", "f32", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Precision.f32, parsed.precision);
}

test "parseArgs --precision=single" {
    const args = [_][]const u8{ "zsasa", "--precision=single", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Precision.f32, parsed.precision);
}

test "parseArgs --precision=double" {
    const args = [_][]const u8{ "zsasa", "--precision=double", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Precision.f64, parsed.precision);
}

// Tests for --parallelism option
test "parseArgs --parallelism default is file" {
    const args = [_][]const u8{ "zsasa", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Parallelism.file, parsed.parallelism);
}

test "parseArgs --parallelism=file" {
    const args = [_][]const u8{ "zsasa", "--parallelism=file", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Parallelism.file, parsed.parallelism);
}

test "parseArgs --parallelism=atom" {
    const args = [_][]const u8{ "zsasa", "--parallelism=atom", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Parallelism.atom, parsed.parallelism);
}

test "parseArgs --parallelism atom (space-separated)" {
    const args = [_][]const u8{ "zsasa", "--parallelism", "atom", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Parallelism.atom, parsed.parallelism);
}

test "parseArgs --parallelism=pipeline" {
    const args = [_][]const u8{ "zsasa", "--parallelism=pipeline", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Parallelism.pipeline, parsed.parallelism);
}

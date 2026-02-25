// Single-structure SASA calculation subcommand
// Calculates SASA for a single input file (JSON, PDB, or mmCIF)
//
const std = @import("std");
const json_writer = @import("json_writer.zig");
const classifier = @import("classifier.zig");
const types = @import("types.zig");

const Precision = types.Precision;
const OutputFormat = json_writer.OutputFormat;
const ClassifierType = classifier.ClassifierType;

/// SASA algorithm selection
pub const Algorithm = enum {
    sr, // Shrake-Rupley (test point method)
    lr, // Lee-Richards (slice method)
};

/// Parsed command-line arguments for the calc subcommand
pub const CalcArgs = struct {
    input_path: ?[]const u8 = null,
    output_path: []const u8 = "output.json",
    output_path_explicit: bool = false, // Track if -o was explicitly set
    n_threads: usize = 0, // 0 = auto-detect
    probe_radius: f64 = 1.4,
    n_points: u32 = 100, // For Shrake-Rupley
    n_slices: u32 = 20, // For Lee-Richards
    algorithm: Algorithm = .sr, // Default: Shrake-Rupley
    precision: Precision = .f64, // f32 or f64 (default: f64)
    output_format: OutputFormat = .json,
    classifier_type: ?ClassifierType = null, // Built-in classifier (naccess/protor/oons)
    config_path: ?[]const u8 = null, // Custom config file path
    chain_filter: ?[]const u8 = null, // Chain filter (e.g., "A" or "A,B,C")
    model_num: ?u32 = null, // Model number for NMR structures
    use_auth_chain: bool = false, // Use auth_asym_id instead of label_asym_id
    include_hydrogens: bool = false, // Include hydrogen atoms (default: exclude)
    include_hetatm: bool = false, // Include HETATM records (default: exclude)
    per_residue: bool = false, // Show per-residue SASA
    rsa: bool = false, // Show RSA (Relative Solvent Accessibility)
    polar: bool = false, // Show polar/nonpolar SASA summary
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 64/128/256)
    quiet: bool = false,
    validate_only: bool = false,
    show_timing: bool = false, // Show timing breakdown for benchmarking
    show_help: bool = false,
};

// =============================================================================
// Parse helper functions
// =============================================================================

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

// =============================================================================
// Argument parsing
// =============================================================================

/// Parse calc subcommand arguments
pub fn parseArgs(args: []const []const u8, start_idx: usize) CalcArgs {
    var result = CalcArgs{};
    var i: usize = start_idx;

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
        // --include-hydrogens (include hydrogen atoms, default: exclude)
        else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
            result.include_hydrogens = true;
        }
        // --include-hetatm (include HETATM records, default: exclude)
        else if (std.mem.eql(u8, arg, "--include-hetatm")) {
            result.include_hetatm = true;
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
        // --use-bitmask (bitmask LUT optimization for SR, n-points must be 64/128/256)
        else if (std.mem.eql(u8, arg, "--use-bitmask")) {
            result.use_bitmask = true;
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
            std.debug.print("Try 'zsasa calc --help' for more information.\n", .{});
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

// =============================================================================
// Help
// =============================================================================

/// Print help for the calc subcommand
pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\zsasa calc - Calculate SASA for a single structure file
        \\
        \\USAGE:
        \\    {s} calc [OPTIONS] <input> [output.json]
        \\
        \\ARGUMENTS:
        \\    <input>          Input file (JSON, PDB, or mmCIF format, auto-detected)
        \\                     Supported: .json, .cif, .mmcif, .pdb, .ent
        \\    [output.json]    Output file (default: output.json)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: naccess, protor, oons
        \\                       Default: protor for PDB/mmCIF, none for JSON
        \\    --config=FILE      Custom classifier config file (FreeSASA format)
        \\    --chain=ID         Filter by chain ID (e.g., --chain=A or --chain=A,B,C)
        \\                       Default: label_asym_id (mmCIF standard)
        \\    --auth-chain       Use auth_asym_id instead of label_asym_id
        \\    --include-hydrogens Include hydrogen atoms (default: excluded)
        \\    --include-hetatm   Include HETATM records (default: excluded)
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
        \\    -o, --output=FILE  Output file (alternative to positional argument)
        \\    --use-bitmask      Use bitmask LUT optimization for SR algorithm
        \\                       Faster but approximate; n-points must be 64, 128, or 256
        \\    --validate         Validate input only, do not calculate SASA
        \\    --timing           Show timing breakdown (for benchmarking)
        \\    -q, --quiet        Suppress progress output
        \\    -h, --help         Show this help message
        \\
        \\ALGORITHMS:
        \\    sr, shrake-rupley  Test point method (default)
        \\    lr, lee-richards   Slice-based method
        \\
        \\CLASSIFIERS:
        \\    protor   ProtOr radii (Tsai et al. 1999) - default for PDB/mmCIF
        \\    naccess  NACCESS-compatible radii
        \\    oons     OONS radii (Ooi et al.)
        \\
        \\OUTPUT FORMATS:
        \\    json     Pretty-printed JSON with indentation
        \\    compact  Single-line JSON (no whitespace)
        \\    csv      CSV with atom_index,area columns
        \\
        \\EXAMPLES:
        \\    {s} calc input.json output.json
        \\    {s} calc structure.cif output.json              # mmCIF input
        \\    {s} calc structure.pdb output.json              # PDB input
        \\    {s} calc --algorithm=lr input.json output.json
        \\    {s} calc --threads=4 input.json output.json
        \\    {s} calc --probe-radius=1.5 --n-points=200 input.json
        \\    {s} calc --format=csv input.json output.csv
        \\    {s} calc --classifier=naccess input.json output.json
        \\    {s} calc --config=custom.config input.json output.json
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
}

// =============================================================================
// Run (stub)
// =============================================================================

/// Run the calc subcommand (stub - implementation in a later task)
pub fn run(allocator: std.mem.Allocator, args: CalcArgs) !void {
    _ = allocator;
    _ = args;
}

// =============================================================================
// Tests
// =============================================================================

test "CalcArgs defaults" {
    const args = [_][]const u8{ "zsasa", "calc", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("output.json", parsed.output_path);
    try std.testing.expectEqual(@as(usize, 0), parsed.n_threads);
    try std.testing.expectEqual(@as(f64, 1.4), parsed.probe_radius);
    try std.testing.expectEqual(@as(u32, 100), parsed.n_points);
    try std.testing.expectEqual(false, parsed.quiet);
    try std.testing.expectEqual(false, parsed.show_help);
}

test "CalcArgs with options" {
    const args = [_][]const u8{ "zsasa", "calc", "--threads=4", "--algorithm=lr", "input.pdb", "result.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
    try std.testing.expectEqualStrings("input.pdb", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.json", parsed.output_path);
}

test "CalcArgs --validate flag" {
    const args = [_][]const u8{ "zsasa", "calc", "--validate", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.validate_only);
}

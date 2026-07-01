// Single-structure SASA calculation subcommand
// Calculates SASA for a single input file (JSON, PDB, or mmCIF)
//
const std = @import("std");
const builtin = @import("builtin");
const analysis = @import("analysis.zig");
const format_detect = @import("format_detect.zig");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const bcif_parser = @import("bcif_parser.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const pdb_parser = @import("pdb_parser.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const classifier = @import("classifier.zig");
const classifier_parser = @import("classifier_parser.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const sdf_parser = @import("sdf_parser.zig");
const compressed = @import("compressed.zig");
const workflow_manifest = @import("workflow_manifest.zig");

const Config = types.Config;
const Configf32 = types.Configf32;
const Precision = types.Precision;
const OutputFormat = json_writer.OutputFormat;
const ClassifierType = classifier.ClassifierType;
const LeeRichardsConfig = lee_richards.LeeRichardsConfig;
const LeeRichardsConfigf32 = lee_richards.LeeRichardsConfigf32;

/// SASA algorithm selection
pub const Algorithm = enum {
    sr, // Shrake-Rupley (test point method)
    lr, // Lee-Richards (slice method)
};

const SdfPathList = sdf_parser.SdfPathList;

/// Parsed command-line arguments for the calc subcommand
pub const CalcArgs = struct {
    input_path: ?[]const u8 = null,
    workflow_path: ?[]const u8 = null,
    output_path: []const u8 = "output.json",
    output_path_explicit: bool = false, // Track if -o or positional output was explicitly set
    n_threads: usize = 0, // 0 = auto-detect
    probe_radius: f64 = 1.4,
    n_points: u32 = 100, // For Shrake-Rupley
    n_slices: u32 = 20, // For Lee-Richards
    algorithm: Algorithm = .sr, // Default: Shrake-Rupley
    precision: Precision = .f64, // f32 or f64 (default: f64)
    output_format: OutputFormat = .json,
    classifier_type: ?ClassifierType = null, // Built-in classifier (naccess/protor/oons/ccd)
    config_path: ?[]const u8 = null, // Custom config file path
    chain_filter: ?[]const u8 = null, // Chain filter (e.g., "A" or "A,B,C")
    model_num: ?u32 = null, // Model number for NMR structures
    use_auth_chain: bool = false, // Use auth_asym_id instead of label_asym_id
    include_hydrogens: bool = false, // Include hydrogen atoms (default: exclude)
    include_hetatm: bool = false, // Include HETATM records (default: exclude)
    per_residue: bool = false, // Show per-residue SASA
    rsa: bool = false, // Show RSA (Relative Solvent Accessibility)
    polar: bool = false, // Show polar/nonpolar SASA summary
    use_bitmask: bool = false, // Use bitmask LUT optimization for SR (n_points must be 1..1024)
    bitmask_correction: bool = false, // Experimental exposed-fraction correction for bitmask SR
    bitmask_correction_coeff: f64 = shrake_rupley_bitmask.default_bitmask_correction_coeff,
    ccd_path: ?[]const u8 = null, // External CCD dictionary file (.zsdc or .cif[.gz|.zst])
    sdf_paths: SdfPathList = .{}, // --sdf=PATH (up to 16)
    mol_selector: ?[]const u8 = null, // --mol=NAME or --mol=N (1-based index)
    quiet: bool = false,
    validate_only: bool = false,
    show_timing: bool = false, // Show timing breakdown for benchmarking
    show_help: bool = false,
    threads_explicit: bool = false,
    probe_radius_explicit: bool = false,
    n_points_explicit: bool = false,
    n_slices_explicit: bool = false,
    algorithm_explicit: bool = false,
    precision_explicit: bool = false,
    format_explicit: bool = false,
    classifier_explicit: bool = false,
    config_explicit: bool = false,
    chain_explicit: bool = false,
    model_explicit: bool = false,
    auth_chain_explicit: bool = false,
    include_hydrogens_explicit: bool = false,
    include_hetatm_explicit: bool = false,
    per_residue_explicit: bool = false,
    rsa_explicit: bool = false,
    polar_explicit: bool = false,
    use_bitmask_explicit: bool = false,
    bitmask_correction_explicit: bool = false,
    bitmask_correction_coeff_explicit: bool = false,
    ccd_explicit: bool = false,
    sdf_explicit: bool = false,
    mol_explicit: bool = false,
    quiet_explicit: bool = false,
    validate_explicit: bool = false,
    timing_explicit: bool = false,
};

// =============================================================================
// Parse helper functions
// =============================================================================

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

/// Parse and validate experimental bitmask correction coefficient.
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

/// Parse and validate output format value
fn parseOutputFormat(value: []const u8) OutputFormat {
    if (std.mem.eql(u8, value, "json")) {
        return .json;
    } else if (std.mem.eql(u8, value, "compact")) {
        return .compact;
    } else if (std.mem.eql(u8, value, "csv")) {
        return .csv;
    } else if (std.mem.eql(u8, value, "freesasa")) {
        return .freesasa;
    } else if (std.mem.eql(u8, value, "rsa")) {
        return .rsa;
    } else {
        std.debug.print("Error: Invalid format: {s}\n", .{value});
        std.debug.print("Valid formats: json, compact, csv, freesasa, rsa\n", .{});
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
    return validateWorkflowNSlices(n) catch {
        std.debug.print("Error: n-slices must be between 1 and 1000: {d}\n", .{n});
        std.process.exit(1);
    };
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

/// Parse chain filter string into array of chain IDs
fn parseChainFilter(allocator: std.mem.Allocator, filter_str: []const u8) ![]const []const u8 {
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

// =============================================================================
// Argument parsing
// =============================================================================

/// Parse calc subcommand arguments
pub fn parseArgs(args: []const []const u8, start_idx: usize) CalcArgs {
    var result = CalcArgs{};
    var i: usize = start_idx;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        // --workflow=PATH or --workflow PATH
        if (std.mem.startsWith(u8, arg, "--workflow=")) {
            result.workflow_path = arg["--workflow=".len..];
        } else if (std.mem.eql(u8, arg, "--workflow")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --workflow\n", .{});
                std.process.exit(1);
            }
            result.workflow_path = args[i];
        }
        // --threads=N or --threads N
        else if (std.mem.startsWith(u8, arg, "--threads=")) {
            const value = arg["--threads=".len..];
            result.n_threads = std.fmt.parseInt(usize, value, 10) catch {
                std.debug.print("Error: Invalid thread count: {s}\n", .{value});
                std.process.exit(1);
            };
            result.threads_explicit = true;
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
            result.threads_explicit = true;
        }
        // --probe-radius=R or --probe-radius R
        else if (std.mem.startsWith(u8, arg, "--probe-radius=")) {
            const value = arg["--probe-radius=".len..];
            result.probe_radius = parseProbeRadius(value);
            result.probe_radius_explicit = true;
        } else if (std.mem.eql(u8, arg, "--probe-radius")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --probe-radius\n", .{});
                std.process.exit(1);
            }
            result.probe_radius = parseProbeRadius(args[i]);
            result.probe_radius_explicit = true;
        }
        // --n-points=N or --n-points N
        else if (std.mem.startsWith(u8, arg, "--n-points=")) {
            const value = arg["--n-points=".len..];
            result.n_points = parseNPoints(value);
            result.n_points_explicit = true;
        } else if (std.mem.eql(u8, arg, "--n-points")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --n-points\n", .{});
                std.process.exit(1);
            }
            result.n_points = parseNPoints(args[i]);
            result.n_points_explicit = true;
        }
        // --format=FORMAT or --format FORMAT
        else if (std.mem.startsWith(u8, arg, "--format=")) {
            const value = arg["--format=".len..];
            result.output_format = parseOutputFormat(value);
            result.format_explicit = true;
        } else if (std.mem.eql(u8, arg, "--format")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --format\n", .{});
                std.process.exit(1);
            }
            result.output_format = parseOutputFormat(args[i]);
            result.format_explicit = true;
        }
        // --algorithm=ALGO or --algorithm ALGO
        else if (std.mem.startsWith(u8, arg, "--algorithm=")) {
            const value = arg["--algorithm=".len..];
            result.algorithm = parseAlgorithm(value);
            result.algorithm_explicit = true;
        } else if (std.mem.eql(u8, arg, "--algorithm")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --algorithm\n", .{});
                std.process.exit(1);
            }
            result.algorithm = parseAlgorithm(args[i]);
            result.algorithm_explicit = true;
        }
        // --n-slices=N or --n-slices N (for Lee-Richards)
        else if (std.mem.startsWith(u8, arg, "--n-slices=")) {
            const value = arg["--n-slices=".len..];
            result.n_slices = parseNSlices(value);
            result.n_slices_explicit = true;
        } else if (std.mem.eql(u8, arg, "--n-slices")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --n-slices\n", .{});
                std.process.exit(1);
            }
            result.n_slices = parseNSlices(args[i]);
            result.n_slices_explicit = true;
        }
        // --classifier=TYPE or --classifier TYPE
        else if (std.mem.startsWith(u8, arg, "--classifier=")) {
            const value = arg["--classifier=".len..];
            result.classifier_type = parseClassifierType(value);
            result.classifier_explicit = true;
        } else if (std.mem.eql(u8, arg, "--classifier")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --classifier\n", .{});
                std.process.exit(1);
            }
            result.classifier_type = parseClassifierType(args[i]);
            result.classifier_explicit = true;
        }
        // --precision=PREC or --precision PREC
        else if (std.mem.startsWith(u8, arg, "--precision=")) {
            const value = arg["--precision=".len..];
            result.precision = parsePrecision(value);
            result.precision_explicit = true;
        } else if (std.mem.eql(u8, arg, "--precision")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --precision\n", .{});
                std.process.exit(1);
            }
            result.precision = parsePrecision(args[i]);
            result.precision_explicit = true;
        }
        // --config=FILE or --config FILE
        else if (std.mem.startsWith(u8, arg, "--config=")) {
            const value = arg["--config=".len..];
            result.config_path = value;
            result.config_explicit = true;
        } else if (std.mem.eql(u8, arg, "--config")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --config\n", .{});
                std.process.exit(1);
            }
            result.config_path = args[i];
            result.config_explicit = true;
        }
        // --chain=ID or --chain ID (e.g., --chain=A or --chain=A,B,C)
        else if (std.mem.startsWith(u8, arg, "--chain=")) {
            const value = arg["--chain=".len..];
            result.chain_filter = value;
            result.chain_explicit = true;
        } else if (std.mem.eql(u8, arg, "--chain")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --chain\n", .{});
                std.process.exit(1);
            }
            result.chain_filter = args[i];
            result.chain_explicit = true;
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
            result.model_explicit = true;
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
            result.model_explicit = true;
        }
        // --auth-chain (use auth_asym_id instead of label_asym_id)
        else if (std.mem.eql(u8, arg, "--auth-chain")) {
            result.use_auth_chain = true;
            result.auth_chain_explicit = true;
        }
        // --include-hydrogens (include hydrogen atoms, default: exclude)
        else if (std.mem.eql(u8, arg, "--include-hydrogens")) {
            result.include_hydrogens = true;
            result.include_hydrogens_explicit = true;
        }
        // --include-hetatm (include HETATM records, default: exclude)
        else if (std.mem.eql(u8, arg, "--include-hetatm")) {
            result.include_hetatm = true;
            result.include_hetatm_explicit = true;
        }
        // --per-residue (show per-residue SASA)
        else if (std.mem.eql(u8, arg, "--per-residue")) {
            result.per_residue = true;
            result.per_residue_explicit = true;
        }
        // --rsa (show RSA - Relative Solvent Accessibility)
        else if (std.mem.eql(u8, arg, "--rsa")) {
            result.rsa = true;
            result.rsa_explicit = true;
            result.per_residue = true; // RSA implies per-residue
            result.per_residue_explicit = true;
        }
        // --polar (show polar/nonpolar SASA summary)
        else if (std.mem.eql(u8, arg, "--polar")) {
            result.polar = true;
            result.polar_explicit = true;
            result.per_residue = true; // Polar analysis requires per-residue
            result.per_residue_explicit = true;
        }
        // --quiet or -q
        else if (std.mem.eql(u8, arg, "--quiet") or std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
            result.quiet_explicit = true;
        }
        // --validate
        else if (std.mem.eql(u8, arg, "--validate")) {
            result.validate_only = true;
            result.validate_explicit = true;
        }
        // --timing
        else if (std.mem.eql(u8, arg, "--timing")) {
            result.show_timing = true;
            result.timing_explicit = true;
        }
        // --ccd=PATH or --ccd PATH (external CCD dictionary)
        else if (std.mem.startsWith(u8, arg, "--ccd=")) {
            const value = arg["--ccd=".len..];
            result.ccd_path = value;
            result.ccd_explicit = true;
        } else if (std.mem.eql(u8, arg, "--ccd")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --ccd\n", .{});
                std.process.exit(1);
            }
            result.ccd_path = args[i];
            result.ccd_explicit = true;
        }
        // --sdf=PATH or --sdf PATH (SDF file with bond topology for CCD classifier)
        else if (std.mem.startsWith(u8, arg, "--sdf=")) {
            const value = arg["--sdf=".len..];
            result.sdf_paths.append(value) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
            result.sdf_explicit = true;
        } else if (std.mem.eql(u8, arg, "--sdf")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --sdf\n", .{});
                std.process.exit(1);
            }
            result.sdf_paths.append(args[i]) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
            result.sdf_explicit = true;
        }
        // --mol=NAME|N or --mol NAME|N (select molecule from multi-molecule SDF)
        else if (std.mem.startsWith(u8, arg, "--mol=")) {
            const value = arg["--mol=".len..];
            result.mol_selector = value;
            result.mol_explicit = true;
        } else if (std.mem.eql(u8, arg, "--mol")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --mol\n", .{});
                std.process.exit(1);
            }
            result.mol_selector = args[i];
            result.mol_explicit = true;
        }
        // --use-bitmask (bitmask LUT optimization for SR, n-points must be 1..1024)
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
// Workflow application
// =============================================================================

fn applyWorkflowToCalcArgs(args: *CalcArgs, workflow: workflow_manifest.Workflow) !void {
    if (args.input_path == null) {
        if (workflow.input.path) |path| args.input_path = path;
    }
    if (!args.chain_explicit) {
        if (workflow.input.chain) |chain| args.chain_filter = chain;
    }
    if (!args.model_explicit) {
        if (workflow.input.model) |model| {
            if (model == 0) {
                if (!builtin.is_test) {
                    std.debug.print("Error: workflow model must be >= 1\n", .{});
                }
                return error.InvalidArgument;
            }
            args.model_num = model;
        }
    }
    if (!args.mol_explicit) {
        if (workflow.input.mol) |mol| args.mol_selector = mol;
    }

    if (!args.output_path_explicit) {
        if (workflow.output.path) |path| args.output_path = path;
    }
    if (!args.format_explicit) {
        if (workflow.output.format) |format| args.output_format = parseOutputFormat(format);
    }

    try applyWorkflowCalculationToCalcArgs(args, workflow.calculation);
    try applyWorkflowClassifierToCalcArgs(args, workflow.classifier);
}

fn applyWorkflowCalculationToCalcArgs(args: *CalcArgs, calculation: workflow_manifest.Calculation) !void {
    if (!args.threads_explicit) {
        if (calculation.threads) |threads| args.n_threads = threads;
    }
    if (!args.probe_radius_explicit) {
        if (calculation.probe_radius) |radius| {
            args.probe_radius = validateWorkflowProbeRadius(radius) catch |err| {
                std.debug.print("Error: workflow probe_radius must be finite and between 0 and 10 Angstroms: {d}\n", .{radius});
                return err;
            };
        }
    }
    if (!args.n_points_explicit) {
        if (calculation.n_points) |n_points| {
            args.n_points = validateWorkflowNPoints(n_points) catch |err| {
                std.debug.print("Error: workflow n_points must be between 1 and 10000: {d}\n", .{n_points});
                return err;
            };
        }
    }
    if (!args.n_slices_explicit) {
        if (calculation.n_slices) |n_slices| {
            args.n_slices = validateWorkflowNSlices(n_slices) catch |err| {
                std.debug.print("Error: workflow n_slices must be between 1 and 1000: {d}\n", .{n_slices});
                return err;
            };
        }
    }
    if (!args.algorithm_explicit) {
        if (calculation.algorithm) |algorithm| args.algorithm = parseAlgorithm(algorithm);
    }
    if (!args.precision_explicit) {
        if (calculation.precision) |precision| args.precision = parsePrecision(precision);
    }
    if (!args.include_hydrogens_explicit) {
        if (calculation.include_hydrogens) |include_hydrogens| args.include_hydrogens = include_hydrogens;
    }
    if (!args.include_hetatm_explicit) {
        if (calculation.include_hetatm) |include_hetatm| args.include_hetatm = include_hetatm;
    }
    if (!args.use_bitmask_explicit) {
        if (calculation.use_bitmask) |use_bitmask| args.use_bitmask = use_bitmask;
    }
    if (!args.timing_explicit) {
        if (calculation.timing) |timing| args.show_timing = timing;
    }
    if (!args.quiet_explicit) {
        if (calculation.quiet) |quiet| args.quiet = quiet;
    }
    if (!args.auth_chain_explicit) {
        if (calculation.auth_chain) |auth_chain| args.use_auth_chain = auth_chain;
    }
    if (!args.per_residue_explicit) {
        if (calculation.per_residue) |per_residue| args.per_residue = per_residue;
    }
    if (!args.rsa_explicit) {
        if (calculation.rsa) |rsa| {
            args.rsa = rsa;
            if (rsa) args.per_residue = true;
        }
    }
    if (!args.polar_explicit) {
        if (calculation.polar) |polar| {
            args.polar = polar;
            if (polar) args.per_residue = true;
        }
    }
    if (!args.validate_explicit) {
        if (calculation.validate_only) |validate_only| args.validate_only = validate_only;
    }
}

fn applyWorkflowClassifierToCalcArgs(args: *CalcArgs, classifier_config: workflow_manifest.ClassifierConfig) !void {
    if (!args.classifier_explicit and !args.config_explicit) {
        if (classifier_config.type) |classifier_type| {
            if (std.mem.eql(u8, classifier_type, "custom")) {
                args.config_path = classifier_config.config orelse return error.InvalidArgument;
                args.classifier_type = null;
            } else {
                args.classifier_type = parseClassifierType(classifier_type);
                args.config_path = null;
            }
        }
    }

    if (calcArgsUseCcdResources(args.*)) {
        if (!args.ccd_explicit) {
            if (classifier_config.ccd) |ccd| args.ccd_path = ccd;
        }
        if (!args.sdf_explicit) {
            if (classifier_config.sdf) |sdf_paths| {
                args.sdf_paths = .{};
                for (sdf_paths) |path| {
                    args.sdf_paths.append(path) catch return error.InvalidArgument;
                }
            }
        }
    } else {
        if (!args.ccd_explicit) args.ccd_path = null;
        if (!args.sdf_explicit) args.sdf_paths = .{};
    }
}

fn classifierUsesCcdResources(effective_classifier_type: ?ClassifierType) bool {
    const classifier_type = effective_classifier_type orelse return false;
    return classifier_type == .ccd;
}

fn calcArgsUseCcdResources(args: CalcArgs) bool {
    if (args.config_path != null) return false;
    return classifierUsesCcdResources(args.classifier_type);
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
        \\    {s} calc [OPTIONS] [input] [output.json]
        \\
        \\ARGUMENTS:
        \\    [input]          Input file (JSON, PDB, or mmCIF format, auto-detected)
        \\                     Supported: .json, .cif, .mmcif, .pdb, .ent
        \\                     Optional when --workflow provides [input].path
        \\    [output.json]    Output file (default: output.json)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: ccd, protor, naccess, oons
        \\                       Default: ccd for PDB/mmCIF, none for JSON
        \\                       protor uses static ProtOr-compatible radii only
        \\    --ccd=PATH         External CCD dictionary file (.zsdc or .cif[.gz|.zst])
        \\                       Extends CCD coverage for non-standard residues
        \\    --sdf=PATH         SDF file with bond topology for CCD classifier
        \\                       Can be specified multiple times for multiple ligands
        \\    --mol=NAME|N       Select molecule from multi-molecule SDF by name or
        \\                       1-based index (default: first molecule)
        \\    --config=FILE      Custom classifier config file (TOML format; .toml only)
        \\    --workflow=PATH    TOML workflow file for input, output, calculation, and classifier settings
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
        \\    --format=FORMAT    Output format: json, compact, csv, freesasa, rsa
        \\                       (default: json; freesasa/rsa are single-calc text formats)
        \\    --precision=PREC   Floating-point precision: f32, f64 (default: f64)
        \\                       f32 is faster but less accurate
        \\    -o, --output=FILE  Output file (alternative to positional argument)
        \\    --use-bitmask      Use bitmask LUT optimization for SR algorithm
        \\                       Faster but approximate; n-points must be 1..1024
        \\    --bitmask-correction
        \\                       Experimental correction for bitmask quantization bias
        \\                       Requires --use-bitmask
        \\    --bitmask-correction-coeff=V
        \\                       Override correction coefficient (default: 0.020)
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
        \\    ccd      CCD bond-topology radii (default for PDB/mmCIF)
        \\    protor   Static ProtOr-compatible radii (Tsai et al. 1999)
        \\    naccess  NACCESS-compatible radii
        \\    oons     OONS radii (Ooi et al.)
        \\
        \\OUTPUT FORMATS:
        \\    json     Pretty-printed JSON with indentation
        \\    compact  Single-line JSON (no whitespace)
        \\    csv      CSV with atom_index,area columns
        \\    freesasa FreeSASA-compatible text summary
        \\    rsa      FreeSASA/NACCESS-compatible RSA residue table
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
        \\    {s} calc --config=custom.toml input.json output.json
        \\    {s} calc --workflow sasa.toml
        \\
    , .{ program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
}

// =============================================================================
// Input reading and classification
// =============================================================================

/// Result from readInputFile — includes optional structure parsers for inline
/// CCD access.
const ReadResult = struct {
    input: types.AtomInput,
    bcif: ?bcif_parser.BcifParser = null,
    mmcif: ?mmcif_parser.MmcifParser = null,
    sdf_ccd: ?ccd_parser.ComponentDict = null, // Auto-registered SDF bond topology

    fn deinitCcd(self: *ReadResult) void {
        if (self.bcif) |*p| p.deinitCcd();
        if (self.mmcif) |*p| p.deinitCcd();
        if (self.sdf_ccd) |*d| d.deinit();
    }
};

/// Select a molecule from parsed SDF by --mol value or default to first.
/// Warns if multiple molecules exist and no --mol specified.
fn selectMolecule(molecules: []const sdf_parser.SdfMolecule, mol_selector: ?[]const u8, quiet: bool) usize {
    if (molecules.len == 0) {
        std.debug.print("Error: SDF file contains no molecules\n", .{});
        std.process.exit(1);
    }

    if (mol_selector) |selector| {
        // Try as 1-based index first
        if (std.fmt.parseInt(usize, selector, 10)) |idx| {
            if (idx >= 1 and idx <= molecules.len) {
                return idx - 1;
            }
            std.debug.print("Error: --mol={d} out of range (SDF has {d} molecules, use 1-based index)\n", .{ idx, molecules.len });
            std.process.exit(1);
        } else |_| {
            // Try as molecule name
            for (molecules, 0..) |mol, i| {
                if (std.mem.eql(u8, std.mem.trim(u8, mol.name, " "), selector)) {
                    return i;
                }
            }
            std.debug.print("Error: --mol='{s}' not found in SDF file\n", .{selector});
            std.process.exit(1);
        }
    }

    // No --mol specified: default to first, warn if multiple
    if (molecules.len > 1 and !quiet) {
        std.debug.print("Warning: SDF contains {d} molecules, processing only the first.\n", .{molecules.len});
        std.debug.print("         Use 'batch' for all molecules or --mol to select one.\n", .{});
    }
    return 0;
}

/// Read input file (auto-detect format)
fn readInputFile(allocator: std.mem.Allocator, io: std.Io, path: []const u8, args: CalcArgs) !ReadResult {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => .{ .input = try json_parser.readAtomInputFromFile(allocator, io, path) },
        .bcif => blk: {
            var parser = bcif_parser.BcifParser.init(allocator);
            parser.model_num = args.model_num;
            parser.use_auth_chain = args.use_auth_chain;
            parser.skip_hydrogens = !args.include_hydrogens;
            parser.atom_only = !args.include_hetatm;
            parser.parse_inline_ccd = calcArgsUseCcdResources(args);

            // Parse chain filter if specified
            var chain_filter_slice: ?[]const []const u8 = null;
            if (args.chain_filter) |filter_str| {
                chain_filter_slice = try parseChainFilter(allocator, filter_str);
                parser.chain_filter = chain_filter_slice;
            }
            defer if (chain_filter_slice) |s| allocator.free(s);

            const input = try parser.parseFile(io, path);
            break :blk .{ .input = input, .bcif = parser };
        },
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.model_num = args.model_num;
            parser.use_auth_chain = args.use_auth_chain;
            parser.skip_hydrogens = !args.include_hydrogens;
            parser.atom_only = !args.include_hetatm;
            parser.parse_inline_ccd = calcArgsUseCcdResources(args);

            // Parse chain filter if specified
            var chain_filter_slice: ?[]const []const u8 = null;
            if (args.chain_filter) |filter_str| {
                chain_filter_slice = try parseChainFilter(allocator, filter_str);
                parser.chain_filter = chain_filter_slice;
            }
            defer if (chain_filter_slice) |s| allocator.free(s);

            const input = try parser.parseFile(io, path);
            break :blk .{ .input = input, .mmcif = parser };
        },
        .pdb => blk: {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.model_num = args.model_num;
            parser.skip_hydrogens = !args.include_hydrogens;
            parser.atom_only = !args.include_hetatm;

            // Parse chain filter if specified
            var chain_filter_slice: ?[]const []const u8 = null;
            if (args.chain_filter) |filter_str| {
                chain_filter_slice = try parseChainFilter(allocator, filter_str);
                parser.chain_filter = chain_filter_slice;
            }
            defer if (chain_filter_slice) |s| allocator.free(s);

            break :blk .{ .input = try parser.parseFile(io, path) };
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

            // Select a single molecule (--mol or default to first)
            const mol_idx = selectMolecule(molecules, args.mol_selector, args.quiet);
            const selected = molecules[mol_idx .. mol_idx + 1];

            const input_result = try sdf_parser.toAtomInput(allocator, selected, !args.include_hydrogens);

            if (!calcArgsUseCcdResources(args)) {
                break :blk .{ .input = input_result };
            }

            // Build CCD component dict from selected molecule's bond topology
            var sdf_dict = ccd_parser.ComponentDict.init(allocator);
            errdefer sdf_dict.deinit();
            var has_components = false;
            for (selected) |mol| {
                if (mol.name.len == 0) continue;
                const stored = sdf_parser.toStoredComponent(allocator, &mol) catch continue;
                const comp_id_str = mol.name[0..@min(mol.name.len, 5)];

                // Skip if already registered (avoid StoredComponent leak from duplicate names)
                if (sdf_dict.components.get(comp_id_str) != null) {
                    var s = stored;
                    s.deinit();
                    continue;
                }

                const dict_key = allocator.dupe(u8, comp_id_str) catch {
                    var s = stored;
                    s.deinit();
                    continue;
                };
                sdf_dict.owned_keys.append(allocator, dict_key) catch {
                    allocator.free(dict_key);
                    var s = stored;
                    s.deinit();
                    continue;
                };
                sdf_dict.components.put(allocator, dict_key, stored) catch {
                    var s = stored;
                    s.deinit();
                    continue;
                };
                has_components = true;
            }

            break :blk .{
                .input = input_result,
                .sdf_ccd = if (has_components) sdf_dict else null_blk: {
                    sdf_dict.deinit();
                    break :null_blk null;
                },
            };
        },
    };
}

/// Apply a custom classifier to input, replacing radii based on residue/atom names
fn applyClassifier(
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
    input: *types.AtomInput,
    ct: ClassifierType,
    sdf_ccd: ?*const ccd_parser.ComponentDict,
    inline_ccd: ?*const ccd_parser.ComponentDict,
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

    // Feed CCD components for non-standard residues
    // Only load components that are actually present in the input structure
    if (ct == .ccd and ccd_clf != null) {
        // Collect unique non-hardcoded residue names from input
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
            const dicts: [3]?*const ccd_parser.ComponentDict = .{ sdf_ccd, inline_ccd, external_ccd };
            for (dicts) |maybe_dict| {
                if (maybe_dict) |dict| {
                    var it = needed.keyIterator();
                    while (it.next()) |key_ptr| {
                        const comp_id = key_ptr.*;
                        if (dict.get(comp_id)) |comp| {
                            ccd_clf.?.addComponent(&comp) catch |err| {
                                if (!quiet) {
                                    std.debug.print("Warning: Could not derive CCD properties for '{s}': {s}\n", .{ comp_id, @errorName(err) });
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

    // Allocate new radii array (use input.allocator for consistency with deinit)
    const new_radii = try input.allocator.alloc(f64, n);
    errdefer input.allocator.free(new_radii);

    var classified_count: usize = 0;
    var fallback_count: usize = 0;

    for (0..n) |i| {
        // Try built-in classifier lookup
        const maybe_radius: ?f64 = switch (ct) {
            .naccess => classifier_naccess.getRadius(residues[i].slice(), atom_names[i].slice()),
            .protor, .ccd => if (ccd_clf) |*c| c.getRadius(residues[i].slice(), atom_names[i].slice()) else null,
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

/// Load SDF files and build a ComponentDict from their bond topology.
/// Delegates to sdf_parser.loadSdfComponents.
const loadSdfComponents = sdf_parser.loadSdfComponents;

const ChainSasaSummary = struct {
    chain_id: types.FixedString4,
    chain_id_full: ?[]const u8 = null,
    area: f64,
    atom_count: usize,

    fn label(self: ChainSasaSummary) []const u8 {
        return self.chain_id_full orelse self.chain_id.slice();
    }
};

fn sameChainSummary(summary: ChainSasaSummary, chain_id: types.FixedString4, chain_id_full: ?[]const u8) bool {
    if (chain_id_full) |full| {
        if (summary.chain_id_full) |summary_full| {
            return std.mem.eql(u8, summary_full, full);
        }
        return false;
    }
    if (summary.chain_id_full != null) return false;
    return std.mem.eql(u8, summary.chain_id.slice(), chain_id.slice());
}

const PerChainBuildResult = struct {
    count: usize,
    overflow: bool,
};

fn buildPerChainSummaries(
    chain_ids: []const types.FixedString4,
    chain_ids_full: ?[]const []const u8,
    atom_areas: []const f64,
    summaries: []ChainSasaSummary,
) PerChainBuildResult {
    var num_chains: usize = 0;
    var overflow = false;

    for (chain_ids, 0..) |chain_id, i| {
        const full = if (chain_ids_full) |full_ids| full_ids[i] else null;
        var found_idx: ?usize = null;
        for (summaries[0..num_chains], 0..) |summary, j| {
            if (sameChainSummary(summary, chain_id, full)) {
                found_idx = j;
                break;
            }
        }

        if (found_idx) |idx| {
            summaries[idx].area += atom_areas[i];
            summaries[idx].atom_count += 1;
        } else if (num_chains < summaries.len) {
            summaries[num_chains] = .{
                .chain_id = chain_id,
                .chain_id_full = full,
                .area = atom_areas[i],
                .atom_count = 1,
            };
            num_chains += 1;
        } else {
            overflow = true;
        }
    }

    return .{ .count = num_chains, .overflow = overflow };
}

/// Print per-chain SASA results
fn printPerChainResults(chain_ids: []const types.FixedString4, chain_ids_full: ?[]const []const u8, atom_areas: []const f64) void {
    // Use a simple approach: iterate through to find unique chains and sum areas
    // For efficiency, we'll use a fixed-size buffer for up to 64 unique chains
    const max_chains = 64;
    var summaries: [max_chains]ChainSasaSummary = undefined;
    const build_result = buildPerChainSummaries(chain_ids, chain_ids_full, atom_areas, summaries[0..]);
    const num_chains = build_result.count;

    if (build_result.overflow) {
        std.debug.print("Warning: More than {d} unique chains; some chains omitted from summary\n", .{max_chains});
    }

    if (num_chains > 0) {
        std.debug.print("\nPer-chain SASA:\n", .{});
        for (0..num_chains) |i| {
            std.debug.print("  {s}: {d:.2} Å² ({d} atoms)\n", .{
                summaries[i].label(),
                summaries[i].area,
                summaries[i].atom_count,
            });
        }
    }
}

fn ensureCalcOutputParentDir(io: std.Io, output_path: []const u8) !void {
    const parent = std.fs.path.dirname(output_path) orelse return;
    if (parent.len == 0) return;
    try std.Io.Dir.cwd().createDirPath(io, parent);
}

fn algorithmName(algorithm: Algorithm) []const u8 {
    return switch (algorithm) {
        .sr => "Shrake & Rupley",
        .lr => "Lee & Richards",
    };
}

fn detailLabel(algorithm: Algorithm) []const u8 {
    return switch (algorithm) {
        .sr => "Test-points",
        .lr => "Slices",
    };
}

fn detailCount(args: CalcArgs) u32 {
    return switch (args.algorithm) {
        .sr => args.n_points,
        .lr => args.n_slices,
    };
}

fn classifierName(effective_classifier: ?ClassifierType, config_path: ?[]const u8) []const u8 {
    if (config_path != null) return "custom";
    if (effective_classifier) |ct| return ct.name();
    return "none";
}

// =============================================================================
// Run
// =============================================================================

/// Run the calc subcommand — single-file SASA processing
pub fn run(allocator: std.mem.Allocator, io: std.Io, args: CalcArgs) !void {
    // Timing variables (in nanoseconds)
    var timer = std.Io.Timestamp.now(io, .awake);
    var time_parse: u64 = 0;
    var time_classify: u64 = 0;
    var time_sasa: u64 = 0;
    var time_write: u64 = 0;

    var effective_args = args;
    var workflow: ?workflow_manifest.Workflow = null;
    defer if (workflow) |*loaded_workflow| loaded_workflow.deinit();
    if (effective_args.workflow_path) |workflow_path| {
        workflow = workflow_manifest.parseFile(allocator, io, workflow_path) catch |err| {
            std.debug.print("Error reading workflow file '{s}': {s}\n", .{ workflow_path, @errorName(err) });
            return err;
        };
        try applyWorkflowToCalcArgs(&effective_args, workflow.?);
    }

    // Validate required arguments
    const input_path = effective_args.input_path orelse {
        std.debug.print("Error: Missing input file\n", .{});
        std.debug.print("Usage: zsasa calc [OPTIONS] [input] [output.json]\n", .{});
        return error.MissingArgument;
    };
    const input_format = format_detect.detectInputFormat(input_path);
    const effective_classifier: ?ClassifierType = effective_args.classifier_type orelse
        if (effective_args.config_path == null and input_format != .json) .ccd else null;
    effective_args.classifier_type = if (effective_args.validate_only) null else effective_classifier;

    // CCD classifier implies HETATM inclusion (the whole point is classifying non-standard residues)
    if (effective_classifier) |ct| {
        if (ct == .ccd and !effective_args.include_hetatm) {
            effective_args.include_hetatm = true;
        }
        // CCD/ProtOr use united-atom radii (implicit H) — warn if explicit H included
        if ((ct == .ccd or ct == .protor) and effective_args.include_hydrogens and !effective_args.quiet) {
            std.debug.print("Warning: --include-hydrogens with CCD classifier may give inaccurate results\n", .{});
            std.debug.print("         CCD uses united-atom radii that already account for implicit hydrogens\n", .{});
        }
    }

    // Read input file (JSON, PDB, or mmCIF)
    timer = std.Io.Timestamp.now(io, .awake);
    var read_result = readInputFile(allocator, io, input_path, effective_args) catch |err| {
        std.debug.print("Error reading input file '{s}': {s}\n", .{ input_path, @errorName(err) });
        std.process.exit(1);
    };
    defer read_result.deinitCcd();
    var input = read_result.input;
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
    if (!effective_args.quiet) {
        _ = json_parser.checkDuplicateCoordinates(allocator, input) catch |err| {
            std.debug.print("Warning: Could not check for duplicate coordinates: {s}\n", .{@errorName(err)});
        };
    }
    time_parse = @intCast(timer.untilNow(io, .awake).nanoseconds);

    // Handle --validate (dry-run)
    if (effective_args.validate_only) {
        if (!effective_args.quiet) {
            std.debug.print("Input validation passed: {} atoms\n", .{input.atomCount()});
        }
        return;
    }

    // Validate --use-bitmask constraints
    if (effective_args.use_bitmask) {
        if (effective_args.algorithm != .sr) {
            std.debug.print("Error: --use-bitmask is only supported with the sr (shrake-rupley) algorithm\n", .{});
            std.process.exit(1);
        }
        if (!bitmask_lut.isSupportedNPoints(effective_args.n_points)) {
            std.debug.print("Error: --use-bitmask requires --n-points to be 1..1024 (got {d})\n", .{effective_args.n_points});
            std.process.exit(1);
        }
    }
    if (effective_args.bitmask_correction and !effective_args.use_bitmask) {
        std.debug.print("Error: --bitmask-correction requires --use-bitmask\n", .{});
        std.process.exit(1);
    }

    // Apply classifier (--config takes precedence over --classifier)
    // Default: ccd for PDB/mmCIF input (ProtOr-compatible with CCD extension)
    timer = std.Io.Timestamp.now(io, .awake);

    if (effective_args.config_path != null or effective_classifier != null) {
        // Warn if both are specified
        if (effective_args.config_path != null and effective_args.classifier_type != null) {
            if (!effective_args.quiet) {
                std.debug.print("Warning: Both --classifier and --config specified; using --config\n", .{});
            }
        }

        // Check if input has classification info
        if (!input.hasClassificationInfo()) {
            std.debug.print("Error: Classifier requires 'residue' and 'atom_name' fields in input\n", .{});
            std.process.exit(1);
        }

        // Load classifier and apply radii
        if (effective_args.config_path) |config_path| {
            // Load from custom config file
            var custom_classifier = classifier_parser.parseConfigFile(allocator, io, config_path) catch |err| {
                switch (err) {
                    error.UnsupportedConfigExtension => std.debug.print("Error loading config file '{s}': custom classifier configs are TOML-only; rename or convert the file to .toml\n", .{config_path}),
                    error.UnsupportedLegacyFormat => std.debug.print("Error loading config file '{s}': FreeSASA-style custom classifier configs are no longer supported; convert to TOML [types] and [[atoms]]\n", .{config_path}),
                    else => std.debug.print("Error loading config file '{s}': {s}\n", .{ config_path, @errorName(err) }),
                }
                std.process.exit(1);
            };
            defer custom_classifier.deinit();

            applyClassifier(&input, &custom_classifier, effective_args.quiet) catch |err| {
                std.debug.print("Error applying classifier: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        } else if (effective_classifier) |ct| {
            // Use built-in classifier
            const use_ccd_resources = classifierUsesCcdResources(ct);
            const inline_ccd: ?*const ccd_parser.ComponentDict = if (use_ccd_resources)
                if (read_result.mmcif) |*p|
                    p.getInlineCcd()
                else if (read_result.bcif) |*p|
                    p.getInlineCcd()
                else
                    null
            else
                null;

            // Load external CCD dictionary if specified
            var ext_ccd: ?ccd_parser.ComponentDict = null;
            if (use_ccd_resources) {
                if (effective_args.ccd_path) |ccd_path| {
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
                        var read_buf: [65536]u8 = undefined;
                        var file_r = f.reader(io, &read_buf);
                        break :blk file_r.interface.allocRemaining(allocator, .unlimited) catch |err| {
                            std.debug.print("Error reading CCD file '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                            std.process.exit(1);
                        };
                    };
                    defer allocator.free(ccd_data);

                    ext_ccd = ccd_binary.loadDict(allocator, ccd_data) catch |err| {
                        std.debug.print("Error loading CCD dictionary '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                        std.process.exit(1);
                    };
                    if (!effective_args.quiet) {
                        std.debug.print("External CCD: loaded {d} components from '{s}'\n", .{ ext_ccd.?.components.count(), ccd_path });
                    }
                }
            }
            defer if (ext_ccd) |*d| d.deinit();

            // Load SDF components from --sdf option
            var sdf_ccd: ?ccd_parser.ComponentDict = null;
            if (use_ccd_resources and effective_args.sdf_paths.len > 0) {
                sdf_ccd = loadSdfComponents(allocator, io, effective_args.sdf_paths.constSlice(), effective_args.quiet) catch |err| {
                    std.debug.print("Error loading SDF components: {s}\n", .{@errorName(err)});
                    std.process.exit(1);
                };
            }
            defer if (sdf_ccd) |*d| d.deinit();

            // SDF auto-components (from SDF input file) + explicit --sdf components
            const sdf_auto_ptr: ?*const ccd_parser.ComponentDict = if (read_result.sdf_ccd != null) &read_result.sdf_ccd.? else null;
            const sdf_explicit_ptr: ?*const ccd_parser.ComponentDict = if (sdf_ccd != null) &sdf_ccd.? else null;
            // Explicit --sdf takes priority over auto-detected
            const sdf_ccd_ptr = sdf_explicit_ptr orelse sdf_auto_ptr;

            const ext_ccd_ptr: ?*const ccd_parser.ComponentDict = if (ext_ccd != null) &ext_ccd.? else null;
            applyBuiltinClassifier(&input, ct, sdf_ccd_ptr, inline_ccd, ext_ccd_ptr, effective_args.quiet) catch |err| {
                std.debug.print("Error applying classifier: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        }
    }
    time_classify = @intCast(timer.untilNow(io, .awake).nanoseconds);

    // Calculate SASA with configured parameters
    timer = std.Io.Timestamp.now(io, .awake);
    var result = switch (effective_args.precision) {
        .f64 => switch (effective_args.algorithm) {
            .sr => blk: {
                const config = Config{
                    .n_points = effective_args.n_points,
                    .probe_radius = effective_args.probe_radius,
                };
                break :blk if (effective_args.use_bitmask) bitmask_blk: {
                    const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f64){
                        .enabled = effective_args.bitmask_correction,
                        .coeff = effective_args.bitmask_correction_coeff,
                    };
                    break :bitmask_blk if (effective_args.n_threads == 1)
                        shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaWithCorrection(allocator, input, config, correction) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f64).calculateSasaParallelWithCorrection(allocator, input, config, effective_args.n_threads, correction) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        };
                } else standard_blk: {
                    break :standard_blk if (effective_args.n_threads == 1)
                        shrake_rupley.calculateSasa(allocator, input, config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        shrake_rupley.calculateSasaParallel(allocator, input, config, effective_args.n_threads) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        };
                };
            },
            .lr => blk: {
                const lr_config = LeeRichardsConfig{
                    .n_slices = effective_args.n_slices,
                    .probe_radius = effective_args.probe_radius,
                };
                break :blk if (effective_args.n_threads == 1)
                    lee_richards.calculateSasa(allocator, input, lr_config) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    }
                else
                    lee_richards.calculateSasaParallel(allocator, input, lr_config, effective_args.n_threads) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    };
            },
        },
        .f32 => blk: {
            // Calculate in f32, convert to f64 for output
            var result_f32 = switch (effective_args.algorithm) {
                .sr => inner: {
                    const config = Configf32{
                        .n_points = effective_args.n_points,
                        .probe_radius = @floatCast(effective_args.probe_radius),
                    };
                    break :inner if (effective_args.use_bitmask) bitmask_inner: {
                        const correction = shrake_rupley_bitmask.BitmaskCorrectionGen(f32){
                            .enabled = effective_args.bitmask_correction,
                            .coeff = @floatCast(effective_args.bitmask_correction_coeff),
                        };
                        break :bitmask_inner if (effective_args.n_threads == 1)
                            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaWithCorrection(allocator, input, config, correction) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            }
                        else
                            shrake_rupley_bitmask.ShrakeRupleyBitmaskGen(f32).calculateSasaParallelWithCorrection(allocator, input, config, effective_args.n_threads, correction) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            };
                    } else standard_inner: {
                        break :standard_inner if (effective_args.n_threads == 1)
                            shrake_rupley.calculateSasaf32(allocator, input, config) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            }
                        else
                            shrake_rupley.calculateSasaParallelf32(allocator, input, config, effective_args.n_threads) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            };
                    };
                },
                .lr => inner: {
                    const lr_config = LeeRichardsConfigf32{
                        .n_slices = effective_args.n_slices,
                        .probe_radius = @floatCast(effective_args.probe_radius),
                    };
                    break :inner if (effective_args.n_threads == 1)
                        lee_richards.calculateSasaf32(allocator, input, lr_config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        lee_richards.calculateSasaParallelf32(allocator, input, lr_config, effective_args.n_threads) catch |err| {
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
    time_sasa = @intCast(timer.untilNow(io, .awake).nanoseconds);

    // Write output file
    timer = std.Io.Timestamp.now(io, .awake);
    ensureCalcOutputParentDir(io, effective_args.output_path) catch |err| {
        std.debug.print("Error creating output directory for '{s}': {s}\n", .{ effective_args.output_path, @errorName(err) });
        std.process.exit(1);
    };
    json_writer.writeSasaResultWithFormatAndInputOptions(allocator, io, result, input, effective_args.output_path, effective_args.output_format, .{
        .input_name = input_path,
        .classifier_name = classifierName(effective_classifier, effective_args.config_path),
        .algorithm_name = algorithmName(effective_args.algorithm),
        .probe_radius = effective_args.probe_radius,
        .detail_count = detailCount(effective_args),
        .detail_label = detailLabel(effective_args.algorithm),
    }) catch |err| {
        if (effective_args.output_format == .rsa and err == error.MissingResidueInfo) {
            std.debug.print("Error: --format=rsa requires chain, residue name, residue number, and insertion code metadata; use PDB/mmCIF input or another structural input with residue metadata\n", .{});
        } else {
            std.debug.print("Error writing output file '{s}': {s}\n", .{ effective_args.output_path, @errorName(err) });
        }
        std.process.exit(1);
    };
    time_write = @intCast(timer.untilNow(io, .awake).nanoseconds);

    // Calculate total time
    const time_total = time_parse + time_classify + time_sasa + time_write;

    // Print summary (unless quiet mode)
    if (!effective_args.quiet) {
        std.debug.print("Calculated SASA for {} atoms\n", .{input.atomCount()});
        std.debug.print("Total area: {d:.2} Å²\n", .{result.total_area});

        // Print per-chain results if chain info is available
        if (input.chain_id) |chain_ids| {
            printPerChainResults(chain_ids, input.chain_id_full, result.atom_areas);
        }

        // Print per-residue results if requested
        if (effective_args.per_residue and input.hasResidueInfo()) {
            var residue_result = analysis.aggregateByResidue(allocator, input, result.atom_areas) catch |err| {
                std.debug.print("Error calculating per-residue SASA: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
            defer residue_result.deinit();
            if (effective_args.rsa) {
                analysis.printResidueResultsWithRsa(residue_result.residues);
            } else {
                analysis.printResidueResults(residue_result.residues);
            }

            // Print polar/nonpolar summary if requested
            if (effective_args.polar) {
                const polar_summary = analysis.calculatePolarSummary(residue_result.residues);
                analysis.printPolarSummary(polar_summary);
            }
        }

        std.debug.print("Output written to {s}\n", .{effective_args.output_path});
    }

    // Print timing breakdown if requested
    if (effective_args.show_timing) {
        const ns_to_ms = 1_000_000.0;
        const parse_ms = @as(f64, @floatFromInt(time_parse)) / ns_to_ms;
        const classify_ms = @as(f64, @floatFromInt(time_classify)) / ns_to_ms;
        const sasa_ms = @as(f64, @floatFromInt(time_sasa)) / ns_to_ms;
        const write_ms = @as(f64, @floatFromInt(time_write)) / ns_to_ms;
        const total_ms = @as(f64, @floatFromInt(time_total)) / ns_to_ms;

        // Human-readable breakdown
        std.debug.print("\nTiming breakdown:\n", .{});
        std.debug.print("  Parse + validate: {d:.2} ms\n", .{parse_ms});
        if (time_classify > 0) {
            std.debug.print("  Classifier:       {d:.2} ms\n", .{classify_ms});
        }
        std.debug.print("  SASA calculation: {d:.2} ms\n", .{sasa_ms});
        std.debug.print("  Write output:     {d:.2} ms\n", .{write_ms});
        std.debug.print("  Total:            {d:.2} ms\n", .{total_ms});

        // Machine-readable format (compatible with freesasa/rustsasa forks)
        // PARSE = file read + validation + classification
        // SASA  = algorithm only
        // TOTAL = PARSE + SASA (excludes output write)
        std.debug.print("PARSE_TIME_MS:{d:.2}\n", .{parse_ms + classify_ms});
        std.debug.print("SASA_TIME_MS:{d:.2}\n", .{sasa_ms});
        std.debug.print("TOTAL_TIME_MS:{d:.2}\n", .{total_ms - write_ms});
    }
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

test "CalcArgs with output path" {
    const args = [_][]const u8{ "zsasa", "calc", "input.json", "result.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.json", parsed.output_path);
}

test "CalcArgs --threads=N" {
    const args = [_][]const u8{ "zsasa", "calc", "--threads=4", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "CalcArgs --threads N (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--threads", "4", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "CalcArgs --probe-radius=R" {
    const args = [_][]const u8{ "zsasa", "calc", "--probe-radius=1.5", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(f64, 1.5), parsed.probe_radius);
}

test "CalcArgs --probe-radius R (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--probe-radius", "1.5", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(f64, 1.5), parsed.probe_radius);
}

test "CalcArgs --n-points=N" {
    const args = [_][]const u8{ "zsasa", "calc", "--n-points=200", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

test "CalcArgs --n-points N (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--n-points", "200", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

test "CalcArgs --quiet" {
    const args = [_][]const u8{ "zsasa", "calc", "--quiet", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.quiet);
}

test "CalcArgs -q" {
    const args = [_][]const u8{ "zsasa", "calc", "-q", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.quiet);
}

test "CalcArgs --help" {
    const args = [_][]const u8{ "zsasa", "calc", "--help" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.show_help);
}

test "CalcArgs -h" {
    const args = [_][]const u8{ "zsasa", "calc", "-h" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.show_help);
}

test "CalcArgs --format=json" {
    const args = [_][]const u8{ "zsasa", "calc", "--format=json", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.json, parsed.output_format);
}

test "CalcArgs --format=compact" {
    const args = [_][]const u8{ "zsasa", "calc", "--format=compact", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.compact, parsed.output_format);
}

test "CalcArgs --format=csv" {
    const args = [_][]const u8{ "zsasa", "calc", "--format=csv", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "CalcArgs --format=freesasa" {
    const args = [_][]const u8{ "zsasa", "calc", "--format=freesasa", "input.pdb" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.freesasa, parsed.output_format);
}

test "CalcArgs --format=rsa" {
    const args = [_][]const u8{ "zsasa", "calc", "--format=rsa", "input.pdb" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.rsa, parsed.output_format);
}

test "CalcArgs --format csv (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--format", "csv", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "CalcArgs --algorithm=sr" {
    const args = [_][]const u8{ "zsasa", "calc", "--algorithm=sr", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "CalcArgs --algorithm=lr" {
    const args = [_][]const u8{ "zsasa", "calc", "--algorithm=lr", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

test "CalcArgs --precision=f32" {
    const args = [_][]const u8{ "zsasa", "calc", "--precision=f32", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Precision.f32, parsed.precision);
}

test "CalcArgs --precision=f64" {
    const args = [_][]const u8{ "zsasa", "calc", "--precision=f64", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Precision.f64, parsed.precision);
}

test "CalcArgs -o FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "-o", "out.json", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("out.json", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "CalcArgs --output=FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "--output=out.json", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("out.json", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "CalcArgs --use-bitmask" {
    const args = [_][]const u8{ "zsasa", "calc", "--use-bitmask", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_bitmask);
}

test "CalcArgs --bitmask-correction" {
    const args = [_][]const u8{ "zsasa", "calc", "--use-bitmask", "--bitmask-correction", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_bitmask);
    try std.testing.expectEqual(true, parsed.bitmask_correction);
}

test "CalcArgs --bitmask-correction-coeff implies correction" {
    const args = [_][]const u8{ "zsasa", "calc", "--use-bitmask", "--bitmask-correction-coeff=0.2", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.bitmask_correction);
    try std.testing.expectEqual(@as(f64, 0.2), parsed.bitmask_correction_coeff);
}

test "CalcArgs --timing" {
    const args = [_][]const u8{ "zsasa", "calc", "--timing", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.show_timing);
}

test "CalcArgs multiple options" {
    const args = [_][]const u8{
        "zsasa",
        "calc",
        "--threads=8",
        "--probe-radius=1.6",
        "--n-points=150",
        "--quiet",
        "input.json",
        "output.json",
    };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(usize, 8), parsed.n_threads);
    try std.testing.expectEqual(@as(f64, 1.6), parsed.probe_radius);
    try std.testing.expectEqual(@as(u32, 150), parsed.n_points);
    try std.testing.expectEqual(true, parsed.quiet);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("output.json", parsed.output_path);
}

test "CalcArgs --classifier=naccess" {
    const args = [_][]const u8{ "zsasa", "calc", "--classifier=naccess", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(ClassifierType.naccess, parsed.classifier_type.?);
}

test "CalcArgs --include-hydrogens" {
    const args = [_][]const u8{ "zsasa", "calc", "--include-hydrogens", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.include_hydrogens);
}

test "CalcArgs --include-hetatm" {
    const args = [_][]const u8{ "zsasa", "calc", "--include-hetatm", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.include_hetatm);
}

test "CalcArgs --per-residue" {
    const args = [_][]const u8{ "zsasa", "calc", "--per-residue", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.per_residue);
}

test "CalcArgs --rsa implies --per-residue" {
    const args = [_][]const u8{ "zsasa", "calc", "--rsa", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.rsa);
    try std.testing.expectEqual(true, parsed.per_residue);
}

test "CalcArgs --polar implies --per-residue" {
    const args = [_][]const u8{ "zsasa", "calc", "--polar", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.polar);
    try std.testing.expectEqual(true, parsed.per_residue);
}

test "CalcArgs --algorithm=shrake-rupley (long form)" {
    const args = [_][]const u8{ "zsasa", "calc", "--algorithm=shrake-rupley", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "CalcArgs --algorithm=lee-richards (long form)" {
    const args = [_][]const u8{ "zsasa", "calc", "--algorithm=lee-richards", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

test "CalcArgs -o takes precedence over positional output" {
    const args = [_][]const u8{ "zsasa", "calc", "-o", "explicit.json", "input.json", "positional.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("explicit.json", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "CalcArgs --n-slices=N" {
    const args = [_][]const u8{ "zsasa", "calc", "--n-slices=50", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 50), parsed.n_slices);
}

test "CalcArgs --n-slices N (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--n-slices", "30", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 30), parsed.n_slices);
}

test "CalcArgs --config=FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "--config=my_classifier.toml", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("my_classifier.toml", parsed.config_path.?);
}

test "CalcArgs --config FILE (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--config", "my_classifier.toml", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("my_classifier.toml", parsed.config_path.?);
}

test "CalcArgs --chain=A" {
    const args = [_][]const u8{ "zsasa", "calc", "--chain=A", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("A", parsed.chain_filter.?);
}

test "CalcArgs --model=1" {
    const args = [_][]const u8{ "zsasa", "calc", "--model=1", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(@as(u32, 1), parsed.model_num.?);
}

test "CalcArgs --auth-chain" {
    const args = [_][]const u8{ "zsasa", "calc", "--auth-chain", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqual(true, parsed.use_auth_chain);
}

test "CalcArgs --output FILE (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--output", "result.json", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("result.json", parsed.output_path);
    try std.testing.expectEqual(true, parsed.output_path_explicit);
}

test "CalcArgs --workflow=FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "--workflow=calc-workflow.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("calc-workflow.toml", parsed.workflow_path.?);
}

test "CalcArgs --workflow FILE" {
    const args = [_][]const u8{ "zsasa", "calc", "--workflow", "calc-workflow.toml" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("calc-workflow.toml", parsed.workflow_path.?);
}

test "ensureCalcOutputParentDir creates nested output parent directories" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const output_path = try std.fs.path.join(allocator, &.{ root_buf[0..root_len], "nested", "deeper", "result.json" });
    defer allocator.free(output_path);

    try ensureCalcOutputParentDir(std.testing.io, output_path);
    _ = try tmp_dir.dir.statFile(std.testing.io, "nested/deeper", .{});
}

fn readAtomAreasLenFromJson(allocator: std.mem.Allocator, path: []const u8) !usize {
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, path, allocator, .limited(4096));
    defer allocator.free(content);
    const parsed = try std.json.parseFromSlice(std.json.Value, allocator, content, .{});
    defer parsed.deinit();
    return parsed.value.object.get("atom_areas").?.array.items.len;
}

test "calc default CCD includes HETATM like explicit CCD" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const pdb_path = try std.fs.path.join(allocator, &.{ root, "hetatm.pdb" });
    defer allocator.free(pdb_path);
    const default_out = try std.fs.path.join(allocator, &.{ root, "default.json" });
    defer allocator.free(default_out);
    const explicit_out = try std.fs.path.join(allocator, &.{ root, "explicit.json" });
    defer allocator.free(explicit_out);

    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = pdb_path, .data = "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
        "HETATM    2  O   HOH A   2      20.000   0.000   0.000  1.00 20.00           O  \n" ++
        "END\n" });

    try run(allocator, std.testing.io, .{
        .input_path = pdb_path,
        .output_path = default_out,
        .n_threads = 1,
        .n_points = 8,
        .quiet = true,
    });
    try run(allocator, std.testing.io, .{
        .input_path = pdb_path,
        .output_path = explicit_out,
        .n_threads = 1,
        .n_points = 8,
        .classifier_type = .ccd,
        .quiet = true,
    });

    try std.testing.expectEqual(@as(usize, 2), try readAtomAreasLenFromJson(allocator, default_out));
    try std.testing.expectEqual(@as(usize, 2), try readAtomAreasLenFromJson(allocator, explicit_out));
}

test "per-chain summaries prefer full chain IDs over truncated prefixes" {
    const chain_ids = [_]types.FixedString4{
        types.FixedString4.fromSlice("ABCD"),
        types.FixedString4.fromSlice("ABCD"),
        types.FixedString4.fromSlice("ABCD"),
    };
    const chain_ids_full = [_][]const u8{ "ABCD1", "ABCD2", "ABCD1" };
    const atom_areas = [_]f64{ 10.0, 20.0, 5.0 };
    var summaries: [4]ChainSasaSummary = undefined;

    const build_result = buildPerChainSummaries(chain_ids[0..], chain_ids_full[0..], atom_areas[0..], summaries[0..]);

    try std.testing.expectEqual(@as(usize, 2), build_result.count);
    try std.testing.expect(!build_result.overflow);
    try std.testing.expectEqualStrings("ABCD1", summaries[0].label());
    try std.testing.expectEqual(@as(f64, 15.0), summaries[0].area);
    try std.testing.expectEqual(@as(usize, 2), summaries[0].atom_count);
    try std.testing.expectEqualStrings("ABCD2", summaries[1].label());
    try std.testing.expectEqual(@as(f64, 20.0), summaries[1].area);
    try std.testing.expectEqual(@as(usize, 1), summaries[1].atom_count);
}

test "calc workflow applies fields when CLI did not override" {
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .input = .{
            .path = "structure.cif",
            .chain = "A,B",
            .model = 2,
            .mol = "LIG",
        },
        .output = .{
            .path = "result.csv",
            .format = "csv",
        },
        .calculation = .{
            .algorithm = "lr",
            .threads = 4,
            .probe_radius = 1.6,
            .n_points = 240,
            .n_slices = 32,
            .precision = "f32",
            .include_hydrogens = true,
            .include_hetatm = true,
            .use_bitmask = true,
            .timing = true,
            .quiet = true,
            .auth_chain = true,
            .per_residue = true,
            .rsa = true,
            .polar = true,
            .validate_only = true,
        },
        .classifier = .{
            .type = "custom",
            .config = "custom-radii.toml",
        },
    };
    var args = CalcArgs{};

    try applyWorkflowToCalcArgs(&args, workflow);

    try std.testing.expectEqualStrings("structure.cif", args.input_path.?);
    try std.testing.expectEqualStrings("A,B", args.chain_filter.?);
    try std.testing.expectEqual(@as(u32, 2), args.model_num.?);
    try std.testing.expectEqualStrings("LIG", args.mol_selector.?);
    try std.testing.expectEqualStrings("result.csv", args.output_path);
    try std.testing.expectEqual(OutputFormat.csv, args.output_format);
    try std.testing.expectEqual(Algorithm.lr, args.algorithm);
    try std.testing.expectEqual(@as(usize, 4), args.n_threads);
    try std.testing.expectEqual(@as(f64, 1.6), args.probe_radius);
    try std.testing.expectEqual(@as(u32, 240), args.n_points);
    try std.testing.expectEqual(@as(u32, 32), args.n_slices);
    try std.testing.expectEqual(Precision.f32, args.precision);
    try std.testing.expectEqual(true, args.include_hydrogens);
    try std.testing.expectEqual(true, args.include_hetatm);
    try std.testing.expectEqual(true, args.use_bitmask);
    try std.testing.expectEqual(true, args.show_timing);
    try std.testing.expectEqual(true, args.quiet);
    try std.testing.expectEqual(true, args.use_auth_chain);
    try std.testing.expectEqual(true, args.per_residue);
    try std.testing.expectEqual(true, args.rsa);
    try std.testing.expectEqual(true, args.polar);
    try std.testing.expectEqual(true, args.validate_only);
    try std.testing.expectEqualStrings("custom-radii.toml", args.config_path.?);
    try std.testing.expectEqual(@as(?ClassifierType, null), args.classifier_type);
    try std.testing.expectEqual(@as(?[]const u8, null), args.ccd_path);
    try std.testing.expectEqual(@as(usize, 0), args.sdf_paths.len);
}

test "calc CLI explicit classifier overrides workflow classifier" {
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{
            .type = "custom",
            .config = "custom-radii.toml",
        },
    };
    var args = CalcArgs{
        .classifier_type = ClassifierType.naccess,
        .classifier_explicit = true,
    };

    try applyWorkflowToCalcArgs(&args, workflow);

    try std.testing.expectEqual(ClassifierType.naccess, args.classifier_type.?);
    try std.testing.expectEqual(@as(?[]const u8, null), args.config_path);
}

test "calc workflow applies CCD resources only for CCD classifier" {
    const sdf_paths = [_][]const u8{ "workflow-ligand.sdf", "workflow-cofactor.sdf" };

    const ccd_workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{
            .type = "ccd",
            .ccd = "workflow.zsdc",
            .sdf = sdf_paths[0..],
        },
    };
    var ccd_args = CalcArgs{};

    try applyWorkflowToCalcArgs(&ccd_args, ccd_workflow);

    try std.testing.expectEqual(ClassifierType.ccd, ccd_args.classifier_type.?);
    try std.testing.expectEqualStrings("workflow.zsdc", ccd_args.ccd_path.?);
    try std.testing.expectEqual(@as(usize, 2), ccd_args.sdf_paths.len);
    try std.testing.expectEqualStrings("workflow-ligand.sdf", ccd_args.sdf_paths.constSlice()[0]);
    try std.testing.expectEqualStrings("workflow-cofactor.sdf", ccd_args.sdf_paths.constSlice()[1]);

    const protor_workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{
            .type = "protor",
            .ccd = "workflow.zsdc",
            .sdf = sdf_paths[0..],
        },
    };
    var protor_args = CalcArgs{};

    try applyWorkflowToCalcArgs(&protor_args, protor_workflow);

    try std.testing.expectEqual(ClassifierType.protor, protor_args.classifier_type.?);
    try std.testing.expectEqual(@as(?[]const u8, null), protor_args.ccd_path);
    try std.testing.expectEqual(@as(usize, 0), protor_args.sdf_paths.len);
}

test "calc CLI explicit CCD resources are preserved for CCD workflow classifier" {
    const workflow_sdf_paths = [_][]const u8{"workflow-ligand.sdf"};
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{
            .type = "ccd",
            .ccd = "workflow.zsdc",
            .sdf = workflow_sdf_paths[0..],
        },
    };
    var args = CalcArgs{
        .ccd_explicit = true,
        .ccd_path = "cli.zsdc",
        .sdf_explicit = true,
    };
    try args.sdf_paths.append("cli.sdf");

    try applyWorkflowToCalcArgs(&args, workflow);

    try std.testing.expectEqual(ClassifierType.ccd, args.classifier_type.?);
    try std.testing.expectEqualStrings("cli.zsdc", args.ccd_path.?);
    try std.testing.expectEqual(@as(usize, 1), args.sdf_paths.len);
    try std.testing.expectEqualStrings("cli.sdf", args.sdf_paths.constSlice()[0]);
}

test "calc CLI explicit naccess ignores workflow CCD resources" {
    const sdf_paths = [_][]const u8{ "workflow-ligand.sdf", "workflow-cofactor.sdf" };
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{
            .type = "ccd",
            .ccd = "missing.zsdc",
            .sdf = sdf_paths[0..],
        },
    };
    var args = CalcArgs{
        .classifier_type = ClassifierType.naccess,
        .classifier_explicit = true,
    };

    try applyWorkflowToCalcArgs(&args, workflow);

    try std.testing.expectEqual(ClassifierType.naccess, args.classifier_type.?);
    try std.testing.expectEqual(@as(?[]const u8, null), args.ccd_path);
    try std.testing.expectEqual(@as(usize, 0), args.sdf_paths.len);
}

test "calc CLI explicit config ignores workflow CCD resources even with CCD classifier" {
    const sdf_paths = [_][]const u8{"workflow-ligand.sdf"};
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .classifier = .{
            .type = "ccd",
            .ccd = "missing.zsdc",
            .sdf = sdf_paths[0..],
        },
    };
    var args = CalcArgs{
        .classifier_type = ClassifierType.ccd,
        .classifier_explicit = true,
        .config_path = "custom-radii.toml",
        .config_explicit = true,
    };

    try applyWorkflowToCalcArgs(&args, workflow);

    try std.testing.expectEqual(ClassifierType.ccd, args.classifier_type.?);
    try std.testing.expectEqualStrings("custom-radii.toml", args.config_path.?);
    try std.testing.expectEqual(@as(?[]const u8, null), args.ccd_path);
    try std.testing.expectEqual(@as(usize, 0), args.sdf_paths.len);
}

test "calc workflow rejects zero model number" {
    const workflow = workflow_manifest.Workflow{
        .allocator = std.testing.allocator,
        .content = "",
        .input = .{
            .model = 0,
        },
    };
    var args = CalcArgs{};

    try std.testing.expectError(error.InvalidArgument, applyWorkflowToCalcArgs(&args, workflow));
}

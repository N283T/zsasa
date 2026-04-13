// Single-structure SASA calculation subcommand
// Calculates SASA for a single input file (JSON, PDB, or mmCIF)
//
const std = @import("std");
const analysis = @import("analysis.zig");
const format_detect = @import("format_detect.zig");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
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
// classifier_protor removed — ProtOr is now an alias for CCD
const classifier_oons = @import("classifier_oons.zig");
const classifier_ccd = @import("classifier_ccd.zig");
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const sdf_parser = @import("sdf_parser.zig");
const gzip = @import("gzip.zig");

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

/// Fixed-capacity list for SDF paths (max 16)
const SdfPathList = struct {
    items: [16][]const u8 = undefined,
    len: usize = 0,

    fn append(self: *SdfPathList, value: []const u8) error{Overflow}!void {
        if (self.len >= 16) return error.Overflow;
        self.items[self.len] = value;
        self.len += 1;
    }

    fn constSlice(self: *const SdfPathList) []const []const u8 {
        return self.items[0..self.len];
    }
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
    ccd_path: ?[]const u8 = null, // External CCD dictionary file (.zsdc or .cif[.gz])
    sdf_paths: SdfPathList = .{}, // --sdf=PATH (up to 16)
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
        // --ccd=PATH or --ccd PATH (external CCD dictionary)
        else if (std.mem.startsWith(u8, arg, "--ccd=")) {
            const value = arg["--ccd=".len..];
            result.ccd_path = value;
        } else if (std.mem.eql(u8, arg, "--ccd")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for --ccd\n", .{});
                std.process.exit(1);
            }
            result.ccd_path = args[i];
        }
        // --sdf=PATH or --sdf PATH (SDF file with bond topology for CCD classifier)
        else if (std.mem.startsWith(u8, arg, "--sdf=")) {
            const value = arg["--sdf=".len..];
            result.sdf_paths.append(value) catch {
                std.debug.print("Error: Too many --sdf paths (max 16)\n", .{});
                std.process.exit(1);
            };
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
        }
        // --use-bitmask (bitmask LUT optimization for SR, n-points must be 1..1024)
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
        \\    --classifier=TYPE  Built-in classifier: ccd, protor, naccess, oons
        \\                       Default: ccd for PDB/mmCIF, none for JSON
        \\                       protor is an alias for ccd
        \\    --ccd=PATH         External CCD dictionary file (.zsdc or .cif[.gz])
        \\                       Extends CCD coverage for non-standard residues
        \\    --sdf=PATH         SDF file with bond topology for CCD classifier
        \\                       Can be specified multiple times for multiple ligands
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
        \\                       Faster but approximate; n-points must be 1..1024
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
        \\    protor   Alias for ccd (ProtOr-compatible, Tsai et al. 1999)
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
// Input reading and classification
// =============================================================================

/// Result from readInputFile — includes optional mmCIF parser for inline CCD access.
const ReadResult = struct {
    input: types.AtomInput,
    mmcif: ?mmcif_parser.MmcifParser = null,
    sdf_ccd: ?ccd_parser.ComponentDict = null, // Auto-registered SDF bond topology

    fn deinitCcd(self: *ReadResult) void {
        if (self.mmcif) |*p| p.deinitCcd();
        if (self.sdf_ccd) |*d| d.deinit();
    }
};

/// Read input file (auto-detect format)
fn readInputFile(allocator: std.mem.Allocator, path: []const u8, args: CalcArgs) !ReadResult {
    const format = format_detect.detectInputFormat(path);
    return switch (format) {
        .json => .{ .input = try json_parser.readAtomInputFromFile(allocator, path) },
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
            parser.model_num = args.model_num;
            parser.use_auth_chain = args.use_auth_chain;
            parser.skip_hydrogens = !args.include_hydrogens;
            parser.atom_only = !args.include_hetatm;

            // Parse chain filter if specified
            var chain_filter_slice: ?[]const []const u8 = null;
            if (args.chain_filter) |filter_str| {
                chain_filter_slice = try parseChainFilter(allocator, filter_str);
                parser.chain_filter = chain_filter_slice;
            }
            defer if (chain_filter_slice) |s| allocator.free(s);

            const input = try parser.parseFile(path);
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

            break :blk .{ .input = try parser.parseFile(path) };
        },
        .sdf => blk: {
            const source = if (std.mem.endsWith(u8, path, ".gz"))
                try gzip.readGzip(allocator, path)
            else file_blk: {
                const f = try std.fs.cwd().openFile(path, .{});
                defer f.close();
                break :file_blk try f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024);
            };
            defer allocator.free(source);

            const molecules = try sdf_parser.parse(allocator, source);
            defer sdf_parser.freeMolecules(allocator, molecules);

            const input_result = try sdf_parser.toAtomInput(allocator, molecules, !args.include_hydrogens);

            // Build CCD component dict from SDF bond topology for auto-classification
            var sdf_dict = ccd_parser.ComponentDict.init(allocator);
            var has_components = false;
            for (molecules) |mol| {
                if (mol.name.len == 0) continue;
                const stored = sdf_parser.toStoredComponent(allocator, &mol) catch continue;
                const comp_id_str = mol.name[0..@min(mol.name.len, 5)];
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
                sdf_dict.components.put(dict_key, stored) catch continue;
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

    // For CCD/ProtOr: create classifier instance and feed inline/external CCD components
    // ProtOr is an alias for CCD (same hardcoded radii, plus runtime CCD analysis)
    var ccd_clf: ?classifier_ccd.CcdClassifier = if (ct == .ccd or ct == .protor) classifier_ccd.CcdClassifier.init(input.allocator) else null;
    defer if (ccd_clf) |*c| c.deinit();

    // Feed CCD components for non-standard residues
    // Only load components that are actually present in the input structure
    if (ccd_clf != null) {
        // Collect unique non-hardcoded residue names from input
        var needed = std.StringHashMap(void).init(input.allocator);
        defer needed.deinit();
        for (0..n) |i| {
            const res = residues[i].slice();
            if (!classifier_ccd.CcdClassifier.isHardcoded(res)) {
                try needed.put(res, {});
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

/// Load SDF files and build a ComponentDict from their bond topology
fn loadSdfComponents(
    allocator: std.mem.Allocator,
    sdf_paths: []const []const u8,
    quiet: bool,
) !?ccd_parser.ComponentDict {
    if (sdf_paths.len == 0) return null;

    var dict = ccd_parser.ComponentDict.init(allocator);
    errdefer dict.deinit();

    for (sdf_paths) |sdf_path| {
        const source = if (std.mem.endsWith(u8, sdf_path, ".gz"))
            gzip.readGzip(allocator, sdf_path) catch |err| {
                std.debug.print("Error reading SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            }
        else file_blk: {
            const f = std.fs.cwd().openFile(sdf_path, .{}) catch |err| {
                std.debug.print("Error opening SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            };
            defer f.close();
            break :file_blk f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024) catch |err| {
                std.debug.print("Error reading SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
                std.process.exit(1);
            };
        };
        defer allocator.free(source);

        const molecules = sdf_parser.parse(allocator, source) catch |err| {
            std.debug.print("Error parsing SDF file '{s}': {s}\n", .{ sdf_path, @errorName(err) });
            std.process.exit(1);
        };
        defer sdf_parser.freeMolecules(allocator, molecules);

        for (molecules) |mol| {
            if (mol.name.len == 0) {
                if (!quiet) std.debug.print("Warning: SDF molecule has no name, skipping\n", .{});
                continue;
            }
            const stored = sdf_parser.toStoredComponent(allocator, &mol) catch |err| {
                if (!quiet) std.debug.print("Warning: Could not convert SDF molecule '{s}': {s}\n", .{ mol.name, @errorName(err) });
                continue;
            };
            const comp_id_str = mol.name[0..@min(mol.name.len, 5)];
            const dict_key = allocator.dupe(u8, comp_id_str) catch continue;
            dict.owned_keys.append(allocator, dict_key) catch {
                allocator.free(dict_key);
                continue;
            };
            dict.components.put(dict_key, stored) catch continue;
        }

        if (!quiet) {
            std.debug.print("SDF: loaded from '{s}'\n", .{sdf_path});
        }
    }

    if (dict.components.count() == 0) {
        dict.deinit();
        return null;
    }
    return dict;
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

// =============================================================================
// Run
// =============================================================================

/// Run the calc subcommand — single-file SASA processing
pub fn run(allocator: std.mem.Allocator, args: CalcArgs) !void {
    // Timing variables (in nanoseconds)
    var timer = std.time.Timer.start() catch {
        std.debug.print("Error: Timer not supported\n", .{});
        std.process.exit(1);
    };
    var time_parse: u64 = 0;
    var time_classify: u64 = 0;
    var time_sasa: u64 = 0;
    var time_write: u64 = 0;

    // Validate required arguments
    const input_path = args.input_path orelse {
        std.debug.print("Error: Missing input file\n", .{});
        std.debug.print("Usage: zsasa calc [OPTIONS] <input> [output.json]\n", .{});
        return error.MissingArgument;
    };

    // CCD classifier implies HETATM inclusion (the whole point is classifying non-standard residues)
    var effective_args = args;
    if (effective_args.classifier_type) |ct| {
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
    timer.reset();
    var read_result = readInputFile(allocator, input_path, effective_args) catch |err| {
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
    if (!args.quiet) {
        _ = json_parser.checkDuplicateCoordinates(allocator, input) catch |err| {
            std.debug.print("Warning: Could not check for duplicate coordinates: {s}\n", .{@errorName(err)});
        };
    }
    time_parse = timer.read();

    // Handle --validate (dry-run)
    if (args.validate_only) {
        if (!args.quiet) {
            std.debug.print("Input validation passed: {} atoms\n", .{input.atomCount()});
        }
        return;
    }

    // Validate --use-bitmask constraints
    if (args.use_bitmask) {
        if (args.algorithm != .sr) {
            std.debug.print("Error: --use-bitmask is only supported with the sr (shrake-rupley) algorithm\n", .{});
            std.process.exit(1);
        }
        if (!bitmask_lut.isSupportedNPoints(args.n_points)) {
            std.debug.print("Error: --use-bitmask requires --n-points to be 1..1024 (got {d})\n", .{args.n_points});
            std.process.exit(1);
        }
    }

    // Apply classifier (--config takes precedence over --classifier)
    // Default: ccd for PDB/mmCIF input (ProtOr-compatible with CCD extension)
    timer.reset();
    const input_format = format_detect.detectInputFormat(args.input_path.?);
    const effective_classifier: ?ClassifierType = args.classifier_type orelse
        if (args.config_path == null and input_format != .json) .ccd else null;

    if (args.config_path != null or effective_classifier != null) {
        // Warn if both are specified
        if (args.config_path != null and args.classifier_type != null) {
            if (!args.quiet) {
                std.debug.print("Warning: Both --classifier and --config specified; using --config\n", .{});
            }
        }

        // Check if input has classification info
        if (!input.hasClassificationInfo()) {
            std.debug.print("Error: Classifier requires 'residue' and 'atom_name' fields in input\n", .{});
            std.process.exit(1);
        }

        // Load classifier and apply radii
        if (args.config_path) |config_path| {
            // Load from custom config file
            var custom_classifier = classifier_parser.parseConfigFile(allocator, config_path) catch |err| {
                std.debug.print("Error loading config file '{s}': {s}\n", .{ config_path, @errorName(err) });
                std.process.exit(1);
            };
            defer custom_classifier.deinit();

            applyClassifier(&input, &custom_classifier, args.quiet) catch |err| {
                std.debug.print("Error applying classifier: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        } else if (effective_classifier) |ct| {
            // Use built-in classifier
            const inline_ccd: ?*const ccd_parser.ComponentDict = if (read_result.mmcif) |*p| p.getInlineCcd() else null;

            // Load external CCD dictionary if specified
            var ext_ccd: ?ccd_parser.ComponentDict = null;
            if (args.ccd_path) |ccd_path| {
                const ccd_data = if (std.mem.endsWith(u8, ccd_path, ".gz"))
                    gzip.readGzip(allocator, ccd_path) catch |err| {
                        std.debug.print("Error reading CCD file '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                        std.process.exit(1);
                    }
                else blk: {
                    const f = std.fs.cwd().openFile(ccd_path, .{}) catch |err| {
                        std.debug.print("Error opening CCD file '{s}': {s}\n", .{ ccd_path, @errorName(err) });
                        std.process.exit(1);
                    };
                    defer f.close();
                    break :blk f.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024) catch |err| {
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
            if (effective_args.sdf_paths.len > 0) {
                sdf_ccd = loadSdfComponents(allocator, effective_args.sdf_paths.constSlice(), effective_args.quiet) catch |err| {
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
            applyBuiltinClassifier(&input, ct, sdf_ccd_ptr, inline_ccd, ext_ccd_ptr, args.quiet) catch |err| {
                std.debug.print("Error applying classifier: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
        }
    }
    time_classify = timer.read();

    // Calculate SASA with configured parameters
    timer.reset();
    var result = switch (args.precision) {
        .f64 => switch (args.algorithm) {
            .sr => blk: {
                const config = Config{
                    .n_points = args.n_points,
                    .probe_radius = args.probe_radius,
                };
                break :blk if (args.use_bitmask) bitmask_blk: {
                    break :bitmask_blk if (args.n_threads == 1)
                        shrake_rupley_bitmask.calculateSasa(allocator, input, config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        shrake_rupley_bitmask.calculateSasaParallel(allocator, input, config, args.n_threads) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        };
                } else standard_blk: {
                    break :standard_blk if (args.n_threads == 1)
                        shrake_rupley.calculateSasa(allocator, input, config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        shrake_rupley.calculateSasaParallel(allocator, input, config, args.n_threads) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        };
                };
            },
            .lr => blk: {
                const lr_config = LeeRichardsConfig{
                    .n_slices = args.n_slices,
                    .probe_radius = args.probe_radius,
                };
                break :blk if (args.n_threads == 1)
                    lee_richards.calculateSasa(allocator, input, lr_config) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    }
                else
                    lee_richards.calculateSasaParallel(allocator, input, lr_config, args.n_threads) catch |err| {
                        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                        std.process.exit(1);
                    };
            },
        },
        .f32 => blk: {
            // Calculate in f32, convert to f64 for output
            var result_f32 = switch (args.algorithm) {
                .sr => inner: {
                    const config = Configf32{
                        .n_points = args.n_points,
                        .probe_radius = @floatCast(args.probe_radius),
                    };
                    break :inner if (args.use_bitmask) bitmask_inner: {
                        break :bitmask_inner if (args.n_threads == 1)
                            shrake_rupley_bitmask.calculateSasaf32(allocator, input, config) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            }
                        else
                            shrake_rupley_bitmask.calculateSasaParallelf32(allocator, input, config, args.n_threads) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            };
                    } else standard_inner: {
                        break :standard_inner if (args.n_threads == 1)
                            shrake_rupley.calculateSasaf32(allocator, input, config) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            }
                        else
                            shrake_rupley.calculateSasaParallelf32(allocator, input, config, args.n_threads) catch |err| {
                                std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                                std.process.exit(1);
                            };
                    };
                },
                .lr => inner: {
                    const lr_config = LeeRichardsConfigf32{
                        .n_slices = args.n_slices,
                        .probe_radius = @floatCast(args.probe_radius),
                    };
                    break :inner if (args.n_threads == 1)
                        lee_richards.calculateSasaf32(allocator, input, lr_config) catch |err| {
                            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
                            std.process.exit(1);
                        }
                    else
                        lee_richards.calculateSasaParallelf32(allocator, input, lr_config, args.n_threads) catch |err| {
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
    json_writer.writeSasaResultWithFormatAndInput(allocator, result, input, args.output_path, args.output_format) catch |err| {
        std.debug.print("Error writing output file '{s}': {s}\n", .{ args.output_path, @errorName(err) });
        std.process.exit(1);
    };
    time_write = timer.read();

    // Calculate total time
    const time_total = time_parse + time_classify + time_sasa + time_write;

    // Print summary (unless quiet mode)
    if (!args.quiet) {
        std.debug.print("Calculated SASA for {} atoms\n", .{input.atomCount()});
        std.debug.print("Total area: {d:.2} Å²\n", .{result.total_area});

        // Print per-chain results if chain info is available
        if (input.chain_id) |chain_ids| {
            printPerChainResults(chain_ids, result.atom_areas);
        }

        // Print per-residue results if requested
        if (args.per_residue and input.hasResidueInfo()) {
            var residue_result = analysis.aggregateByResidue(allocator, input, result.atom_areas) catch |err| {
                std.debug.print("Error calculating per-residue SASA: {s}\n", .{@errorName(err)});
                std.process.exit(1);
            };
            defer residue_result.deinit();
            if (args.rsa) {
                analysis.printResidueResultsWithRsa(residue_result.residues);
            } else {
                analysis.printResidueResults(residue_result.residues);
            }

            // Print polar/nonpolar summary if requested
            if (args.polar) {
                const polar_summary = analysis.calculatePolarSummary(residue_result.residues);
                analysis.printPolarSummary(polar_summary);
            }
        }

        std.debug.print("Output written to {s}\n", .{args.output_path});
    }

    // Print timing breakdown if requested
    if (args.show_timing) {
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
    const args = [_][]const u8{ "zsasa", "calc", "--config=my_classifier.conf", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("my_classifier.conf", parsed.config_path.?);
}

test "CalcArgs --config FILE (space-separated)" {
    const args = [_][]const u8{ "zsasa", "calc", "--config", "my_classifier.conf", "input.json" };
    const parsed = parseArgs(&args, 2);
    try std.testing.expectEqualStrings("my_classifier.conf", parsed.config_path.?);
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

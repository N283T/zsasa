const std = @import("std");
const build_options = @import("build_options");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const mmcif_parser = @import("mmcif_parser.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const classifier = @import("classifier.zig");
const classifier_parser = @import("classifier_parser.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_protor = @import("classifier_protor.zig");
const classifier_oons = @import("classifier_oons.zig");

const Config = types.Config;
const OutputFormat = json_writer.OutputFormat;
const LeeRichardsConfig = lee_richards.LeeRichardsConfig;
const ClassifierType = classifier.ClassifierType;

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
    n_threads: usize = 0, // 0 = auto-detect
    probe_radius: f64 = 1.4,
    n_points: u32 = 100, // For Shrake-Rupley
    n_slices: u32 = 20, // For Lee-Richards
    algorithm: Algorithm = .sr, // Default: Shrake-Rupley
    output_format: OutputFormat = .json,
    classifier_type: ?ClassifierType = null, // Built-in classifier (naccess/protor/oons)
    config_path: ?[]const u8 = null, // Custom config file path
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

/// Input file format
const InputFormat = enum {
    json,
    mmcif,
};

/// Detect input file format from extension
fn detectInputFormat(path: []const u8) InputFormat {
    if (std.mem.endsWith(u8, path, ".cif")) return .mmcif;
    if (std.mem.endsWith(u8, path, ".mmcif")) return .mmcif;
    if (std.mem.endsWith(u8, path, ".CIF")) return .mmcif;
    if (std.mem.endsWith(u8, path, ".mmCIF")) return .mmcif;
    return .json;
}

/// Read input file (auto-detect format)
fn readInputFile(allocator: std.mem.Allocator, path: []const u8) !types.AtomInput {
    const format = detectInputFormat(path);
    return switch (format) {
        .json => json_parser.readAtomInputFromFile(allocator, path),
        .mmcif => blk: {
            var parser = mmcif_parser.MmcifParser.init(allocator);
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
        // Unknown option
        else if (std.mem.startsWith(u8, arg, "-")) {
            std.debug.print("Error: Unknown option: {s}\n", .{arg});
            std.debug.print("Try '--help' for more information.\n", .{});
            std.process.exit(1);
        }
        // Positional arguments
        else if (result.input_path == null) {
            result.input_path = arg;
        } else {
            result.output_path = arg;
        }
    }

    return result;
}

fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\freesasa_zig {s} - Solvent Accessible Surface Area calculator
        \\
        \\USAGE:
        \\    {s} [OPTIONS] <input> [output.json]
        \\
        \\ARGUMENTS:
        \\    <input>          Input file (JSON or mmCIF format, auto-detected by extension)
        \\                     Supported: .json, .cif, .mmcif
        \\    [output.json]    Output JSON file (default: output.json)
        \\
        \\OPTIONS:
        \\    --algorithm=ALGO   Algorithm: sr (shrake-rupley), lr (lee-richards)
        \\                       Default: sr
        \\    --classifier=TYPE  Built-in classifier: naccess, protor, oons
        \\                       Use with residue/atom_name input for auto-radius
        \\    --config=FILE      Custom classifier config file (FreeSASA format)
        \\    --threads=N        Number of threads (default: auto-detect)
        \\                       Use --threads=1 for single-threaded mode
        \\    --probe-radius=R   Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N       Test points per atom (default: 100, for sr)
        \\    --n-slices=N       Slices per atom diameter (default: 20, for lr)
        \\    --format=FORMAT    Output format: json, compact, csv (default: json)
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
        \\    {s} --algorithm=lr input.json output.json
        \\    {s} --algorithm=lr --n-slices=50 input.json output.json
        \\    {s} --threads=4 input.json output.json
        \\    {s} --probe-radius=1.5 --n-points=200 input.json
        \\    {s} --format=csv input.json output.csv
        \\    {s} --classifier=naccess input.json output.json
        \\    {s} --classifier=naccess structure.cif output.json
        \\    {s} --config=custom.config input.json output.json
        \\
    , .{ version, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name, program_name });
}

fn printVersion() void {
    std.debug.print("freesasa_zig {s}\n", .{version});
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
        if (custom_classifier.getRadius(residues[i], atom_names[i])) |r| {
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
        } else if (classifier.guessRadiusFromAtomName(atom_names[i])) |r| {
            // Fall back to atom name-based radius
            new_radii[i] = r;
            fallback_count += 1;
        } else {
            // Keep original radius
            new_radii[i] = input.r[i];
        }
    }

    // Free old radii and replace
    // TODO: Refactor AtomInput to have mutable radii field instead of using @constCast
    input.allocator.free(@constCast(input.r));
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
            .naccess => classifier_naccess.getRadius(residues[i], atom_names[i]),
            .protor => classifier_protor.getRadius(residues[i], atom_names[i]),
            .oons => classifier_oons.getRadius(residues[i], atom_names[i]),
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
        } else if (classifier.guessRadiusFromAtomName(atom_names[i])) |r| {
            // Fall back to atom name-based radius
            new_radii[i] = r;
            fallback_count += 1;
        } else {
            // Keep original radius
            new_radii[i] = input.r[i];
        }
    }

    // Free old radii and replace
    // TODO: Refactor AtomInput to have mutable radii field instead of using @constCast
    input.allocator.free(@constCast(input.r));
    input.r = new_radii;

    if (!quiet) {
        std.debug.print("Classifier '{s}': {d} atoms classified, {d} fallback\n", .{
            ct.name(),
            classified_count,
            fallback_count,
        });
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
    var input = readInputFile(allocator, parsed.input_path.?) catch |err| {
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
    var result = switch (parsed.algorithm) {
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
            // FIX: Previously used `if (n_threads > 1)` which made n_threads=0 (auto-detect)
            // fall through to single-threaded mode. Now consistent with SR: use parallel
            // unless explicitly single-threaded (n_threads=1).
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
    };
    defer result.deinit();
    time_sasa = timer.read();

    // Write output file
    timer.reset();
    json_writer.writeSasaResultWithFormat(allocator, result, parsed.output_path, parsed.output_format) catch |err| {
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
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
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
    const args = [_][]const u8{ "freesasa_zig", "input.json", "result.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
    try std.testing.expectEqualStrings("result.json", parsed.output_path);
}

test "parseArgs --threads=N" {
    const args = [_][]const u8{ "freesasa_zig", "--threads=4", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs --probe-radius=R" {
    const args = [_][]const u8{ "freesasa_zig", "--probe-radius=1.5", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(f64, 1.5), parsed.probe_radius);
}

test "parseArgs --n-points=N" {
    const args = [_][]const u8{ "freesasa_zig", "--n-points=200", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

test "parseArgs --quiet" {
    const args = [_][]const u8{ "freesasa_zig", "--quiet", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.quiet);
}

test "parseArgs -q" {
    const args = [_][]const u8{ "freesasa_zig", "-q", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.quiet);
}

test "parseArgs --help" {
    const args = [_][]const u8{ "freesasa_zig", "--help" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_help);
}

test "parseArgs -h" {
    const args = [_][]const u8{ "freesasa_zig", "-h" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_help);
}

test "parseArgs --version" {
    const args = [_][]const u8{ "freesasa_zig", "--version" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_version);
}

test "parseArgs -V" {
    const args = [_][]const u8{ "freesasa_zig", "-V" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_version);
}

test "parseArgs multiple options" {
    const args = [_][]const u8{
        "freesasa_zig",
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
    const args = [_][]const u8{ "freesasa_zig", "--threads", "4", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(usize, 4), parsed.n_threads);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs --probe-radius R (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--probe-radius", "1.5", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(f64, 1.5), parsed.probe_radius);
}

test "parseArgs --n-points N (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--n-points", "200", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 200), parsed.n_points);
}

// Tests for --format option
test "parseArgs --format=json" {
    const args = [_][]const u8{ "freesasa_zig", "--format=json", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.json, parsed.output_format);
}

test "parseArgs --format=compact" {
    const args = [_][]const u8{ "freesasa_zig", "--format=compact", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.compact, parsed.output_format);
}

test "parseArgs --format=csv" {
    const args = [_][]const u8{ "freesasa_zig", "--format=csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "parseArgs --format csv (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--format", "csv", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.csv, parsed.output_format);
}

test "parseArgs default format is json" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(OutputFormat.json, parsed.output_format);
}

test "parseArgs --validate" {
    const args = [_][]const u8{ "freesasa_zig", "--validate", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.validate_only);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

test "parseArgs default validate_only is false" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(false, parsed.validate_only);
}

// Tests for --algorithm option
test "parseArgs default algorithm is sr" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "parseArgs --algorithm=sr" {
    const args = [_][]const u8{ "freesasa_zig", "--algorithm=sr", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "parseArgs --algorithm=lr" {
    const args = [_][]const u8{ "freesasa_zig", "--algorithm=lr", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

test "parseArgs --algorithm=shrake-rupley" {
    const args = [_][]const u8{ "freesasa_zig", "--algorithm=shrake-rupley", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.sr, parsed.algorithm);
}

test "parseArgs --algorithm=lee-richards" {
    const args = [_][]const u8{ "freesasa_zig", "--algorithm=lee-richards", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

test "parseArgs --algorithm lr (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--algorithm", "lr", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
}

// Tests for --n-slices option
test "parseArgs default n_slices is 20" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 20), parsed.n_slices);
}

test "parseArgs --n-slices=50" {
    const args = [_][]const u8{ "freesasa_zig", "--n-slices=50", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 50), parsed.n_slices);
}

test "parseArgs --n-slices 100 (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--n-slices", "100", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(u32, 100), parsed.n_slices);
}

test "parseArgs combined algorithm and n-slices" {
    const args = [_][]const u8{ "freesasa_zig", "--algorithm=lr", "--n-slices=50", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(Algorithm.lr, parsed.algorithm);
    try std.testing.expectEqual(@as(u32, 50), parsed.n_slices);
}

// Tests for --classifier option
test "parseArgs default classifier_type is null" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(?ClassifierType, null), parsed.classifier_type);
}

test "parseArgs --classifier=naccess" {
    const args = [_][]const u8{ "freesasa_zig", "--classifier=naccess", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.naccess, parsed.classifier_type.?);
}

test "parseArgs --classifier=protor" {
    const args = [_][]const u8{ "freesasa_zig", "--classifier=protor", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.protor, parsed.classifier_type.?);
}

test "parseArgs --classifier=oons" {
    const args = [_][]const u8{ "freesasa_zig", "--classifier=oons", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.oons, parsed.classifier_type.?);
}

test "parseArgs --classifier protor (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--classifier", "protor", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(ClassifierType.protor, parsed.classifier_type.?);
}

// Tests for --config option
test "parseArgs default config_path is null" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(@as(?[]const u8, null), parsed.config_path);
}

test "parseArgs --config=custom.config" {
    const args = [_][]const u8{ "freesasa_zig", "--config=custom.config", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("custom.config", parsed.config_path.?);
}

test "parseArgs --config custom.config (space-separated)" {
    const args = [_][]const u8{ "freesasa_zig", "--config", "custom.config", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqualStrings("custom.config", parsed.config_path.?);
}

test "parseArgs combined classifier and config (config takes precedence)" {
    const args = [_][]const u8{ "freesasa_zig", "--classifier=naccess", "--config=custom.config", "input.json" };
    const parsed = parseArgs(&args);

    // Both are set - config should take precedence in main() logic
    try std.testing.expectEqual(ClassifierType.naccess, parsed.classifier_type.?);
    try std.testing.expectEqualStrings("custom.config", parsed.config_path.?);
}

// Tests for --timing option
test "parseArgs default show_timing is false" {
    const args = [_][]const u8{ "freesasa_zig", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(false, parsed.show_timing);
}

test "parseArgs --timing" {
    const args = [_][]const u8{ "freesasa_zig", "--timing", "input.json" };
    const parsed = parseArgs(&args);

    try std.testing.expectEqual(true, parsed.show_timing);
    try std.testing.expectEqualStrings("input.json", parsed.input_path.?);
}

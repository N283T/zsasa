const std = @import("std");
const build_options = @import("build_options");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const types = @import("types.zig");

const Config = types.Config;

const version = build_options.version;

/// Parsed command-line arguments
const Args = struct {
    input_path: ?[]const u8 = null,
    output_path: []const u8 = "output.json",
    n_threads: usize = 0, // 0 = auto-detect
    probe_radius: f64 = 1.4,
    n_points: u32 = 100,
    quiet: bool = false,
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
        // --quiet or -q
        else if (std.mem.eql(u8, arg, "--quiet") or std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
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
        \\    {s} [OPTIONS] <input.json> [output.json]
        \\
        \\ARGUMENTS:
        \\    <input.json>     Input JSON file with atom coordinates and radii
        \\    [output.json]    Output JSON file (default: output.json)
        \\
        \\OPTIONS:
        \\    --threads=N        Number of threads (default: auto-detect)
        \\                       Use --threads=1 for single-threaded mode
        \\    --probe-radius=R   Probe radius in Angstroms (default: 1.4)
        \\    --n-points=N       Test points per atom (default: 100)
        \\    -q, --quiet        Suppress progress output
        \\    -h, --help         Show this help message
        \\    -V, --version      Show version
        \\
        \\EXAMPLES:
        \\    {s} input.json output.json
        \\    {s} --threads=4 input.json output.json
        \\    {s} --probe-radius=1.5 --n-points=200 input.json
        \\    {s} --quiet input.json output.json
        \\
    , .{ version, program_name, program_name, program_name, program_name, program_name });
}

fn printVersion() void {
    std.debug.print("freesasa_zig {s}\n", .{version});
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

    // Read input JSON file
    var input = json_parser.readAtomInputFromFile(allocator, parsed.input_path.?) catch |err| {
        std.debug.print("Error reading input file '{s}': {s}\n", .{ parsed.input_path.?, @errorName(err) });
        std.process.exit(1);
    };
    defer input.deinit();

    // Calculate SASA with configured parameters
    const config = Config{
        .n_points = parsed.n_points,
        .probe_radius = parsed.probe_radius,
    };

    var result = if (parsed.n_threads == 1)
        shrake_rupley.calculateSasa(allocator, input, config) catch |err| {
            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        }
    else
        shrake_rupley.calculateSasaParallel(allocator, input, config, parsed.n_threads) catch |err| {
            std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    defer result.deinit();

    // Write output JSON file
    json_writer.writeSasaResult(allocator, result, parsed.output_path) catch |err| {
        std.debug.print("Error writing output file '{s}': {s}\n", .{ parsed.output_path, @errorName(err) });
        std.process.exit(1);
    };

    // Print summary (unless quiet mode)
    if (!parsed.quiet) {
        std.debug.print("Calculated SASA for {} atoms\n", .{input.atomCount()});
        std.debug.print("Total area: {d:.2} Å²\n", .{result.total_area});
        std.debug.print("Output written to {s}\n", .{parsed.output_path});
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

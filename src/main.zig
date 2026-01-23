const std = @import("std");
const json_parser = @import("json_parser.zig");
const json_writer = @import("json_writer.zig");
const shrake_rupley = @import("shrake_rupley.zig");
const types = @import("types.zig");

const Config = types.Config;

pub fn main() !void {
    // Setup allocator
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse command line arguments
    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        try printUsage(args[0]);
        std.process.exit(1);
    }

    const input_path = args[1];
    const output_path = if (args.len >= 3) args[2] else "output.json";

    // Read input JSON file
    var input = json_parser.readAtomInputFromFile(allocator, input_path) catch |err| {
        std.debug.print("Error reading input file '{s}': {s}\n", .{ input_path, @errorName(err) });
        std.process.exit(1);
    };
    defer input.deinit();

    // Calculate SASA with default config
    const config = Config{};
    var result = shrake_rupley.calculateSasa(allocator, input, config) catch |err| {
        std.debug.print("Error calculating SASA: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer result.deinit();

    // Write output JSON file
    json_writer.writeSasaResult(allocator, result, output_path) catch |err| {
        std.debug.print("Error writing output file '{s}': {s}\n", .{ output_path, @errorName(err) });
        std.process.exit(1);
    };

    // Print summary
    std.debug.print("Calculated SASA for {} atoms\n", .{input.atomCount()});
    std.debug.print("Total area: {d:.2} Ų\n", .{result.total_area});
    std.debug.print("Output written to {s}\n", .{output_path});
}

fn printUsage(program_name: []const u8) !void {
    std.debug.print("Usage: {s} <input.json> [output.json]\n", .{program_name});
    std.debug.print("\n", .{});
    std.debug.print("Arguments:\n", .{});
    std.debug.print("  input.json   Input JSON file with atom coordinates and radii\n", .{});
    std.debug.print("  output.json  Output JSON file (default: output.json)\n", .{});
}

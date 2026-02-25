const std = @import("std");
const build_options = @import("build_options");
const batch = @import("batch.zig");
const calc = @import("calc.zig");
const traj = @import("traj.zig");

const version = build_options.version;

fn isSubcommand(arg: []const u8) bool {
    return std.mem.eql(u8, arg, "calc") or
        std.mem.eql(u8, arg, "batch") or
        std.mem.eql(u8, arg, "traj");
}

fn printUsage(program_name: []const u8) void {
    std.debug.print(
        \\zsasa {s} - Solvent Accessible Surface Area calculator
        \\
        \\USAGE:
        \\    {s} <command> [OPTIONS] <args>
        \\
        \\COMMANDS:
        \\    calc    Calculate SASA for a single structure file
        \\    batch   Calculate SASA for all files in a directory
        \\    traj    Calculate SASA across trajectory frames
        \\
        \\GLOBAL OPTIONS:
        \\    -h, --help       Show this help message
        \\    -V, --version    Show version
        \\
        \\Use '{s} <command> --help' for more information about a command.
        \\
    , .{ version, program_name, program_name });
}

fn printVersion() void {
    std.debug.print("zsasa {s}\n", .{version});
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const args = try std.process.argsAlloc(allocator);
    defer std.process.argsFree(allocator, args);

    if (args.len < 2) {
        printUsage(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];

    // Global flags
    if (std.mem.eql(u8, subcmd, "--help") or std.mem.eql(u8, subcmd, "-h")) {
        printUsage(args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "--version") or std.mem.eql(u8, subcmd, "-V")) {
        printVersion();
        return;
    }

    // Subcommand dispatch
    if (std.mem.eql(u8, subcmd, "calc")) {
        const calc_args = calc.parseArgs(args, 2);
        if (calc_args.show_help) {
            calc.printHelp(args[0]);
            return;
        }
        calc.run(allocator, calc_args) catch |err| {
            std.debug.print("Error: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else if (std.mem.eql(u8, subcmd, "batch")) {
        const batch_args = batch.parseArgs(args, 2);
        if (batch_args.show_help) {
            batch.printHelp(args[0]);
            return;
        }
        batch.run(allocator, batch_args) catch |err| {
            std.debug.print("Error: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else if (std.mem.eql(u8, subcmd, "traj")) {
        const traj_args = traj.parseArgs(args, 2);
        if (traj_args.show_help) {
            traj.printHelp(args[0]);
            return;
        }
        traj.run(allocator, traj_args) catch |err| {
            std.debug.print("Error in trajectory mode: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else {
        std.debug.print("Error: unknown subcommand '{s}'\n", .{subcmd});
        printUsage(args[0]);
        std.process.exit(1);
    }
}

// Tests
test {
    // Discover tests from subcommand modules
    _ = calc;
    _ = batch;
    _ = traj;
}

test "isSubcommand" {
    try std.testing.expect(isSubcommand("calc"));
    try std.testing.expect(isSubcommand("batch"));
    try std.testing.expect(isSubcommand("traj"));
    try std.testing.expect(!isSubcommand("unknown"));
    try std.testing.expect(!isSubcommand("--help"));
}

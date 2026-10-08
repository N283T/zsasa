const std = @import("std");
const build_options = @import("build_options");
const batch = @import("batch.zig");
const calc = @import("calc.zig");
const cli_output = @import("cli_output.zig");
const compile_dict = @import("compile_dict.zig");
const traj = @import("traj.zig");

const version = build_options.version;

fn isSubcommand(arg: []const u8) bool {
    return std.mem.eql(u8, arg, "calc") or
        std.mem.eql(u8, arg, "batch") or
        std.mem.eql(u8, arg, "traj") or
        std.mem.eql(u8, arg, "compile-dict");
}

const usage_format =
    \\zsasa {s} - Solvent Accessible Surface Area calculator
    \\
    \\USAGE:
    \\    {s} <command> [OPTIONS] <args>
    \\
    \\COMMANDS:
    \\    calc          Calculate SASA for a single structure file
    \\    batch         Calculate SASA for all files in a directory
    \\    traj          Calculate SASA across trajectory frames
    \\    compile-dict  Compile CIF dictionary to binary ZSDC format
    \\
    \\GLOBAL OPTIONS:
    \\    -h, --help       Show this help message
    \\    -V, --version    Show version
    \\
    \\Use '{s} <command> --help' for more information about a command.
    \\
;

/// Usage on standard error, next to the error message that precedes it.
fn printUsageError(program_name: []const u8) void {
    std.debug.print(usage_format, .{ version, program_name, program_name });
}

/// Usage on standard output, for `--help`.
fn printUsage(io: std.Io, program_name: []const u8) void {
    cli_output.print(io, usage_format, .{ version, program_name, program_name });
}

fn printVersion(io: std.Io) void {
    cli_output.print(io, "zsasa {s}\n", .{version});
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    // Collect args into a slice using the arena allocator
    const args_z = try init.minimal.args.toSlice(init.arena.allocator());
    // Coerce [:0]const u8 -> []const u8 for each element
    const args: []const []const u8 = @ptrCast(args_z);

    if (args.len < 2) {
        printUsageError(args[0]);
        std.process.exit(1);
    }

    const subcmd = args[1];

    // Global flags
    if (std.mem.eql(u8, subcmd, "--help") or std.mem.eql(u8, subcmd, "-h")) {
        printUsage(io, args[0]);
        return;
    }
    if (std.mem.eql(u8, subcmd, "--version") or std.mem.eql(u8, subcmd, "-V")) {
        printVersion(io);
        return;
    }

    // Subcommand dispatch
    if (std.mem.eql(u8, subcmd, "calc")) {
        const calc_args = calc.parseArgs(args, 2);
        if (calc_args.show_help) {
            calc.printHelp(io, args[0]);
            return;
        }
        calc.run(allocator, io, calc_args) catch |err| {
            std.debug.print("Error: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else if (std.mem.eql(u8, subcmd, "batch")) {
        const batch_args = batch.parseArgs(args, 2);
        if (batch_args.show_help) {
            batch.printHelp(io, args[0]);
            return;
        }
        batch.run(allocator, io, batch_args) catch |err| {
            std.debug.print("Error: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else if (std.mem.eql(u8, subcmd, "traj")) {
        const traj_args = traj.parseArgs(args, 2);
        if (traj_args.show_help) {
            traj.printHelp(io, args[0]);
            return;
        }
        traj.run(allocator, io, traj_args) catch |err| {
            std.debug.print("Error in trajectory mode: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else if (std.mem.eql(u8, subcmd, "compile-dict")) {
        // Check for --help before running
        for (args[2..]) |a| {
            if (std.mem.eql(u8, a, "--help") or std.mem.eql(u8, a, "-h")) {
                compile_dict.printHelp(io, args[0]);
                return;
            }
        }
        compile_dict.run(allocator, io, args[2..]) catch |err| {
            std.debug.print("Error: {s}\n", .{@errorName(err)});
            std.process.exit(1);
        };
    } else {
        std.debug.print("Error: unknown subcommand '{s}'\n", .{subcmd});
        printUsageError(args[0]);
        std.process.exit(1);
    }
}

// Tests
test {
    // Discover tests from subcommand modules
    _ = calc;
    _ = batch;
    _ = traj;
    _ = compile_dict;
}

test "isSubcommand" {
    try std.testing.expect(isSubcommand("calc"));
    try std.testing.expect(isSubcommand("batch"));
    try std.testing.expect(isSubcommand("traj"));
    try std.testing.expect(isSubcommand("compile-dict"));
    try std.testing.expect(!isSubcommand("unknown"));
    try std.testing.expect(!isSubcommand("--help"));
}

//! compile-dict subcommand: convert CIF text CCD data to ZSDC binary format.
//!
//! Usage:
//!   zsasa compile-dict <input.cif[.gz]> -o <output.zsdc>
//!
//! Reads a CIF (Chemical Component Dictionary) file, parses all components,
//! and writes a compact binary dictionary suitable for use with `--ccd`.

const std = @import("std");
const Allocator = std.mem.Allocator;
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const gzip = @import("gzip.zig");

pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\zsasa compile-dict - Compile CIF dictionary to binary ZSDC format
        \\
        \\USAGE:
        \\    {s} compile-dict <input.cif[.gz]> -o <output.zsdc>
        \\
        \\ARGUMENTS:
        \\    <input>          Input CIF file (supports .gz compression)
        \\
        \\OPTIONS:
        \\    -o, --output PATH  Output binary dictionary file (required)
        \\    -h, --help         Show this help message
        \\
        \\EXAMPLES:
        \\    {s} compile-dict components.cif -o components.zsdc
        \\    {s} compile-dict components.cif.gz -o components.zsdc
        \\
    , .{ program_name, program_name, program_name });
}

pub fn run(allocator: Allocator, io: std.Io, args: []const []const u8) !void {
    _ = io;
    var input_path: ?[]const u8 = null;
    var output_path: ?[]const u8 = null;
    var show_help = false;

    var i: usize = 0;
    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            show_help = true;
        } else if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            if (i >= args.len) {
                std.debug.print("Error: Missing value for {s}\n", .{arg});
                std.process.exit(1);
            }
            output_path = args[i];
        } else if (std.mem.startsWith(u8, arg, "--output=")) {
            output_path = arg["--output=".len..];
        } else if (std.mem.startsWith(u8, arg, "-")) {
            std.debug.print("Error: Unknown option: {s}\n", .{arg});
            std.process.exit(1);
        } else if (input_path == null) {
            input_path = arg;
        } else {
            std.debug.print("Error: Unexpected argument: {s}\n", .{arg});
            std.process.exit(1);
        }
    }

    if (show_help) {
        // Caller handles this; but just in case:
        printHelp("zsasa");
        return;
    }

    const in_path = input_path orelse {
        std.debug.print("Error: Missing input file\n", .{});
        std.debug.print("Usage: zsasa compile-dict <input.cif[.gz]> -o <output.zsdc>\n", .{});
        return error.MissingArgument;
    };

    const out_path = output_path orelse {
        std.debug.print("Error: Missing output file. Use -o <output.zsdc>\n", .{});
        return error.MissingArgument;
    };

    // Read input file (handle .gz)
    std.debug.print("Reading '{s}'...\n", .{in_path});
    const source = if (std.mem.endsWith(u8, in_path, ".gz"))
        try gzip.readGzip(allocator, in_path)
    else blk: {
        const file = std.fs.cwd().openFile(in_path, .{}) catch |err| {
            std.debug.print("Error: Could not open '{s}': {s}\n", .{ in_path, @errorName(err) });
            std.process.exit(1);
        };
        defer file.close();
        break :blk file.readToEndAlloc(allocator, 4 * 1024 * 1024 * 1024) catch |err| {
            std.debug.print("Error: Could not read '{s}': {s}\n", .{ in_path, @errorName(err) });
            std.process.exit(1);
        };
    };
    defer allocator.free(source);

    // Parse CIF
    std.debug.print("Parsing CIF data ({d} bytes)...\n", .{source.len});
    var dict = ccd_parser.parseCcdData(allocator, source, null) catch |err| {
        std.debug.print("Error: Failed to parse CIF data: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    defer dict.deinit();

    const comp_count = dict.components.count();
    std.debug.print("Parsed {d} components\n", .{comp_count});

    // Write binary output
    const out_file = std.fs.cwd().createFile(out_path, .{}) catch |err| {
        std.debug.print("Error: Could not create '{s}': {s}\n", .{ out_path, @errorName(err) });
        std.process.exit(1);
    };
    defer out_file.close();

    var write_buf: [64 * 1024]u8 = undefined;
    var buffered = out_file.writer(&write_buf);
    ccd_binary.writeDict(&buffered.interface, &dict) catch |err| {
        std.debug.print("Error: Failed to write binary dict: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };
    buffered.interface.flush() catch |err| {
        std.debug.print("Error: Failed to flush output: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };

    std.debug.print("Compiled {d} components to '{s}'\n", .{ comp_count, out_path });
}

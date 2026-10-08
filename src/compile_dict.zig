//! compile-dict subcommand: convert CIF text CCD data to ZSDC binary format.
//!
//! Usage:
//!   zsasa compile-dict <input.cif[.gz|.zst]> -o <output.zsdc>
//!
//! Reads a CIF (Chemical Component Dictionary) file, parses all components,
//! and writes a compact binary dictionary suitable for use with `--ccd`.

const std = @import("std");
const Allocator = std.mem.Allocator;
const ccd_parser = @import("ccd_parser.zig");
const ccd_binary = @import("ccd_binary.zig");
const compressed = @import("compressed.zig");

pub fn printHelp(program_name: []const u8) void {
    std.debug.print(
        \\zsasa compile-dict - Compile CIF dictionary to binary ZSDC format
        \\
        \\USAGE:
        \\    {s} compile-dict <input.cif[.gz|.zst]> -o <output.zsdc>
        \\
        \\ARGUMENTS:
        \\    <input>          Input CIF file (supports .gz and .zst compression)
        \\
        \\OPTIONS:
        \\    -o, --output PATH  Output binary dictionary file (required)
        \\    -h, --help         Show this help message
        \\
        \\EXAMPLES:
        \\    {s} compile-dict components.cif -o components.zsdc
        \\    {s} compile-dict components.cif.zst -o components.zsdc
        \\
    , .{ program_name, program_name, program_name });
}

pub fn run(allocator: Allocator, io: std.Io, args: []const []const u8) !void {
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
        std.debug.print("Usage: zsasa compile-dict <input.cif[.gz|.zst]> -o <output.zsdc>\n", .{});
        return error.MissingArgument;
    };

    const out_path = output_path orelse {
        std.debug.print("Error: Missing output file. Use -o <output.zsdc>\n", .{});
        return error.MissingArgument;
    };

    // Read input file (handle .gz/.zst)
    std.debug.print("Reading '{s}'...\n", .{in_path});
    const source = if (compressed.isCompressed(in_path))
        try compressed.read(allocator, in_path)
    else blk: {
        const file = std.Io.Dir.cwd().openFile(io, in_path, .{}) catch |err| {
            std.debug.print("Error: Could not open '{s}': {s}\n", .{ in_path, @errorName(err) });
            std.process.exit(1);
        };
        defer file.close(io);
        var read_buf: [65536]u8 = undefined;
        var file_r = file.reader(io, &read_buf);
        break :blk file_r.interface.allocRemaining(allocator, .unlimited) catch |err| {
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
    const out_file = std.Io.Dir.cwd().createFile(io, out_path, .{}) catch |err| {
        std.debug.print("Error: Could not create '{s}': {s}\n", .{ out_path, @errorName(err) });
        std.process.exit(1);
    };
    defer out_file.close(io);

    var write_buf: [64 * 1024]u8 = undefined;
    var buffered = out_file.writer(io, &write_buf);
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

const test_support = @import("test_support.zig");

const test_cif =
    \\data_ALA
    \\#
    \\loop_
    \\_chem_comp_atom.comp_id
    \\_chem_comp_atom.atom_id
    \\_chem_comp_atom.type_symbol
    \\_chem_comp_atom.pdbx_aromatic_flag
    \\_chem_comp_atom.pdbx_leaving_atom_flag
    \\ALA N   N N N
    \\ALA CA  C N N
    \\ALA C   C N N
    \\ALA O   O N N
    \\ALA OXT O N Y
    \\#
    \\loop_
    \\_chem_comp_bond.comp_id
    \\_chem_comp_bond.atom_id_1
    \\_chem_comp_bond.atom_id_2
    \\_chem_comp_bond.value_order
    \\_chem_comp_bond.pdbx_aromatic_flag
    \\ALA N   CA  SING N
    \\ALA CA  C   SING N
    \\ALA C   O   DOUB N
    \\ALA C   OXT SING N
    \\#
    \\data_HOH
    \\#
    \\loop_
    \\_chem_comp_atom.comp_id
    \\_chem_comp_atom.atom_id
    \\_chem_comp_atom.type_symbol
    \\HOH O  O
    \\HOH H1 H
    \\HOH H2 H
    \\#
;

/// Write `test_cif` into a temporary directory and return the paths of the
/// input and of the (not yet existing) output, allocated from `arena`.
fn prepareCompile(tmp: *std.testing.TmpDir, arena: Allocator) ![2][]const u8 {
    try tmp.dir.writeFile(std.testing.io, .{ .sub_path = "components.cif", .data = test_cif });
    const input = try tmp.dir.realPathFileAlloc(std.testing.io, "components.cif", arena);
    const dir = std.fs.path.dirname(input) orelse return error.TestUnexpectedResult;
    return .{ input, try std.fs.path.join(arena, &.{ dir, "components.zsdc" }) };
}

test "compile-dict writes the dictionary that parsing the CIF gives, with -o and --output=" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();
    const input, const output = try prepareCompile(&tmp, arena.allocator());

    // The expected bytes: the same parse, written by the binary writer.
    var expected_dict = try ccd_parser.parseCcdData(allocator, test_cif, null);
    defer expected_dict.deinit();
    var expected_out: std.Io.Writer.Allocating = .init(allocator);
    defer expected_out.deinit();
    try ccd_binary.writeDict(&expected_out.writer, &expected_dict);

    // The three spellings of the output option.
    const output_equals = try std.fmt.allocPrint(arena.allocator(), "--output={s}", .{output});
    const argument_lists = [_][]const []const u8{
        &.{ input, "-o", output },
        &.{ input, "--output", output },
        &.{ output_equals, input },
    };
    for (argument_lists) |arguments| {
        // A stale file from the previous form must not satisfy the check.
        tmp.dir.deleteFile(std.testing.io, "components.zsdc") catch {};
        try run(allocator, std.testing.io, arguments);

        const written = try tmp.dir.readFileAlloc(std.testing.io, "components.zsdc", allocator, .limited(1 << 20));
        defer allocator.free(written);
        try std.testing.expect(ccd_binary.isBinaryDict(written));
        try std.testing.expectEqualSlices(u8, expected_out.written(), written);
    }

    // And the file loads back with what the CIF listed.
    const written = try tmp.dir.readFileAlloc(std.testing.io, "components.zsdc", allocator, .limited(1 << 20));
    defer allocator.free(written);
    var loaded = try ccd_binary.loadDict(allocator, written);
    defer loaded.deinit();
    try std.testing.expectEqual(@as(usize, 2), loaded.count());
    const ala = loaded.get("ALA") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 5), ala.atoms.len);
    try std.testing.expectEqual(@as(usize, 4), ala.bonds.len);
    const hoh = loaded.get("HOH") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 3), hoh.atoms.len);
    try std.testing.expectEqual(@as(usize, 0), hoh.bonds.len);
}

test "compile-dict reports a missing input or output without writing anything" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    const allocator = std.testing.allocator;
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(allocator);
    defer arena.deinit();
    const input, const output = try prepareCompile(&tmp, arena.allocator());

    try std.testing.expectError(error.MissingArgument, run(allocator, std.testing.io, &.{}));
    try std.testing.expectError(error.MissingArgument, run(allocator, std.testing.io, &.{ "-o", output }));
    try std.testing.expectError(error.MissingArgument, run(allocator, std.testing.io, &.{input}));
    try std.testing.expectError(error.FileNotFound, tmp.dir.statFile(std.testing.io, "components.zsdc", .{}));
}

test "compile-dict --help does not need an input and writes nothing" {
    var muted = test_support.muteStderr();
    defer muted.restore();
    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();
    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);
    defer arena.deinit();
    const input, const output = try prepareCompile(&tmp, arena.allocator());

    try run(std.testing.allocator, std.testing.io, &.{"--help"});
    try run(std.testing.allocator, std.testing.io, &.{ "-h", input, "-o", output });
    try std.testing.expectError(error.FileNotFound, tmp.dir.statFile(std.testing.io, "components.zsdc", .{}));
}

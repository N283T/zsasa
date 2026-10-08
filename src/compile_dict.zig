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

    // Parse and serialize completely in memory first, so a failure leaves no
    // output file (and no stale partial one) behind.
    std.debug.print("Parsing CIF data ({d} bytes)...\n", .{source.len});
    var diag: ccd_binary.WriteDiagnostic = .{};
    var comp_count: usize = 0;
    const bytes = compileToBytes(allocator, source, &comp_count, &diag) catch |err| {
        printCompileError(err, in_path, diag);
        std.process.exit(1);
    };
    defer allocator.free(bytes);
    std.debug.print("Parsed {d} components\n", .{comp_count});

    const out_file = std.Io.Dir.cwd().createFile(io, out_path, .{}) catch |err| {
        std.debug.print("Error: Could not create '{s}': {s}\n", .{ out_path, @errorName(err) });
        std.process.exit(1);
    };
    defer out_file.close(io);

    var write_buf: [64 * 1024]u8 = undefined;
    var buffered = out_file.writer(io, &write_buf);
    buffered.interface.writeAll(bytes) catch |err| {
        std.debug.print("Error: Failed to write '{s}': {s}\n", .{ out_path, @errorName(err) });
        std.process.exit(1);
    };
    buffered.interface.flush() catch |err| {
        std.debug.print("Error: Failed to flush output: {s}\n", .{@errorName(err)});
        std.process.exit(1);
    };

    std.debug.print("Compiled {d} components to '{s}'\n", .{ comp_count, out_path });
}

/// Parse CIF text and serialize it to ZSDC bytes. Caller owns the result.
///
/// Fails with `NoComponents` for input that holds no components (an empty
/// dictionary would only produce a file every later `--ccd` run silently
/// ignores). `diag` names the component behind a writer error.
pub fn compileToBytes(
    allocator: Allocator,
    source: []const u8,
    comp_count: *usize,
    diag: *ccd_binary.WriteDiagnostic,
) ![]u8 {
    var dict = try ccd_parser.parseCcdData(allocator, source, null);
    defer dict.deinit();

    comp_count.* = dict.components.count();
    if (comp_count.* == 0) return error.NoComponents;

    var out: std.Io.Writer.Allocating = .init(allocator);
    errdefer out.deinit();
    ccd_binary.writeDictDiag(&out.writer, &dict, diag) catch |err| switch (err) {
        error.WriteFailed => return error.OutOfMemory,
        else => |e| return e,
    };
    return out.toOwnedSlice();
}

fn printCompileError(err: anyerror, in_path: []const u8, diag: ccd_binary.WriteDiagnostic) void {
    switch (err) {
        error.NoComponents => std.debug.print(
            "Error: '{s}' contains no CCD components (expected _chem_comp_atom loops); nothing to compile\n",
            .{in_path},
        ),
        error.TooManyAtoms => if (diag.comp_id.len > 0) std.debug.print(
            "Error: component '{s}' has more than 65535 atoms, which the ZSDC format cannot store\n",
            .{diag.comp_id},
        ) else std.debug.print(
            "Error: '{s}' has a component with more than 65535 atoms, which the ZSDC format cannot store\n",
            .{in_path},
        ),
        error.TooManyBonds => std.debug.print(
            "Error: component '{s}' has more than 65535 bonds, which the ZSDC format cannot store\n",
            .{diag.comp_id},
        ),
        error.InvalidComponentId => std.debug.print(
            "Error: component ID '{s}' is empty or longer than 255 bytes, which the ZSDC format cannot store\n",
            .{diag.comp_id},
        ),
        error.InvalidBondIndex => std.debug.print(
            "Error: component '{s}' has a bond that refers to a missing atom\n",
            .{diag.comp_id},
        ),
        else => std.debug.print("Error: Failed to parse CIF data: {s}\n", .{@errorName(err)}),
    }
}

// =============================================================================
// Tests
// =============================================================================

const test_cif_header =
    \\data_T
    \\loop_
    \\_chem_comp_atom.comp_id
    \\_chem_comp_atom.atom_id
    \\_chem_comp_atom.type_symbol
    \\
;

test "compileToBytes rejects input without components" {
    const allocator = std.testing.allocator;
    var diag: ccd_binary.WriteDiagnostic = .{};
    var count: usize = 99;

    try std.testing.expectError(error.NoComponents, compileToBytes(allocator, "data_x\n", &count, &diag));
    try std.testing.expectEqual(@as(usize, 0), count);
    try std.testing.expectError(error.NoComponents, compileToBytes(allocator, "", &count, &diag));
    // A loop header without rows is still no component.
    try std.testing.expectError(error.NoComponents, compileToBytes(allocator, test_cif_header, &count, &diag));
}

test "compileToBytes produces a loadable image" {
    const allocator = std.testing.allocator;
    var diag: ccd_binary.WriteDiagnostic = .{};
    var count: usize = 0;

    const bytes = try compileToBytes(allocator, test_cif_header ++ "GLY N N\nGLY CA C\n", &count, &diag);
    defer allocator.free(bytes);
    try std.testing.expectEqual(@as(usize, 1), count);

    var dict = try ccd_binary.loadDict(allocator, bytes);
    defer dict.deinit();
    const comp = dict.get("GLY") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 2), comp.atoms.len);
}

test "compileToBytes does not wrap the atom count of a huge component" {
    const allocator = std.testing.allocator;

    var cif: std.Io.Writer.Allocating = .init(allocator);
    defer cif.deinit();
    try cif.writer.writeAll(test_cif_header);
    for (0..70_000) |i| try cif.writer.print("BIG C{d} C\n", .{i});

    var diag: ccd_binary.WriteDiagnostic = .{};
    var count: usize = 0;
    try std.testing.expectError(error.TooManyAtoms, compileToBytes(allocator, cif.written(), &count, &diag));
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

//! Gzip decompression using native std.compress.flate.
//!
//! Restored after #320's C-zlib workaround. The upstream panic
//! (https://github.com/ziglang/zig/issues/25035) is fixed in Zig 0.16.

const std = @import("std");

pub const GzipError = error{ GzipOpenFailed, GzipReadFailed, FileTooLarge, OutOfMemory };

const CHUNK_SIZE = 64 * 1024; // 64 KB read chunks

/// Default max decompressed size: 4 GB.
pub const DEFAULT_MAX_SIZE: usize = 4 * 1024 * 1024 * 1024;

/// Decompress a gzip file. Caller owns the returned slice.
pub fn readGzip(allocator: std.mem.Allocator, path: []const u8) GzipError![]u8 {
    return readGzipLimited(allocator, path, DEFAULT_MAX_SIZE);
}

const ReadOptions = struct {
    log_errors: bool = true,
};

/// Decompress a gzip file with a custom size limit. Caller owns the returned slice.
/// The size limit is checked before each chunk read so a malicious file cannot
/// allocate more than `max_size` bytes before the cap is detected.
pub fn readGzipLimited(allocator: std.mem.Allocator, path: []const u8, max_size: usize) GzipError![]u8 {
    return readGzipLimitedOptions(allocator, path, max_size, .{});
}

fn readGzipLimitedOptions(allocator: std.mem.Allocator, path: []const u8, max_size: usize, options: ReadOptions) GzipError![]u8 {
    // NOTE: gzip.zig deliberately does not take an `io: std.Io` parameter to
    // preserve the existing public API used by all callers (mmcif/pdb/json/sdf
    // parsers + batch + compile_dict + the FFI surface in c_api.zig). The
    // function-local single-threaded Threaded is sufficient for the synchronous
    // file-read + decompress pipeline used here. We use a function-local
    // instance (rather than std.Io.Threaded.global_single_threaded) so the
    // threading guarantee is statically obvious — gzip.zig is reachable from
    // batch.zig worker threads spawned via std.Thread.spawn.
    var threaded: std.Io.Threaded = .init_single_threaded;
    const io = threaded.io();

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch return error.GzipOpenFailed;
    defer file.close(io);

    var file_buf: [CHUNK_SIZE]u8 = undefined;
    var file_reader = file.reader(io, &file_buf);
    const input: *std.Io.Reader = &file_reader.interface;

    // Window buffer for the native gzip decompressor, shared by all members.
    // It must be at least flate.max_window_len (64 KB) to support generic
    // reads that go through the indirect (buffered) vtable used by
    // readSliceShort/readVec.
    // 64 KB on stack — required by the indirect (buffered) Decompress vtable
    // path used by readSliceShort/readVec; see flate.history_len.
    var window_buf: [std.compress.flate.max_window_len]u8 = undefined;

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    errdefer buf.deinit(allocator);

    // A gzip file is a series of members whose decompressed data is
    // concatenated (RFC 1952 section 2.2). `cat a.gz b.gz` produces one, and
    // so does bgzip: a BGZF file is many small members ending with an empty
    // one. flate.Decompress handles a single member and leaves the input
    // positioned right after its trailer, so decode members one at a time
    // until the input is exhausted.
    //
    // Whatever follows a member must be another valid member. Other trailing
    // bytes, including the NUL padding that gzip(1) skips, are an error:
    // there is no way to warn about ignored input from here, and silently
    // ignoring it is how truncated structures went unnoticed before.
    var ctx: MemberContext = .{ .path = path, .member = 1, .log_errors = options.log_errors };
    while (true) : (ctx.member += 1) {
        try readMember(allocator, input, &window_buf, &buf, max_size, ctx);
        if (!try hasNextMember(input, ctx)) break;
    }

    return buf.toOwnedSlice(allocator);
}

const MemberContext = struct {
    path: []const u8,
    /// 1-based index of the member within the file, for diagnostics.
    member: usize,
    log_errors: bool,
};

const gzip_magic = [_]u8{ 0x1f, 0x8b };

/// After member `ctx.member`: returns false at the end of `input` and true if
/// another gzip member starts there. Any other remaining bytes are an error.
/// Only the magic is checked here; flate.Decompress validates the full header.
fn hasNextMember(input: *std.Io.Reader, ctx: MemberContext) error{GzipReadFailed}!bool {
    const next = input.peek(gzip_magic.len) catch |err| switch (err) {
        // Fewer bytes than the magic are left; they stay buffered.
        error.EndOfStream => input.buffered(),
        error.ReadFailed => {
            if (ctx.log_errors) {
                std.log.warn("gzip read failed for {s} after member {d}", .{ ctx.path, ctx.member });
            }
            return error.GzipReadFailed;
        },
    };
    if (next.len == 0) return false;
    if (std.mem.eql(u8, next, &gzip_magic)) return true;
    if (ctx.log_errors) {
        std.log.warn("gzip decode failed for {s}: trailing data after member {d} is not a gzip member", .{ ctx.path, ctx.member });
    }
    return error.GzipReadFailed;
}

/// Decode one gzip member from `input`, append its data to `buf` and verify
/// the member's CRC32 / ISIZE trailer. `max_size` caps `buf.items.len`, i.e.
/// the total across every member decoded so far, not this member alone.
/// On error `buf` may hold part of the member; the caller owns and frees it.
fn readMember(
    allocator: std.mem.Allocator,
    input: *std.Io.Reader,
    window_buf: []u8,
    buf: *std.ArrayListUnmanaged(u8),
    max_size: usize,
    ctx: MemberContext,
) GzipError!void {
    var decompress: std.compress.flate.Decompress = .init(input, .gzip, window_buf);
    const reader: *std.Io.Reader = &decompress.reader;

    // Compute CRC32 over the decompressed bytes so we can verify the gzip
    // trailer ourselves. std.compress.flate parses the trailer fields into
    // decompress.container_metadata.gzip.{crc,count} but does NOT compare them
    // against the actual decoded bytes. Without this verification, a corrupt
    // .cif.gz with a valid frame structure but a bad checksum would silently
    // return wrong bytes — a real correctness regression for the SASA pipeline
    // (the previous C-zlib version verified CRC via gzclose).
    // CRC32 and ISIZE cover one member each, so both start over per member.
    var crc = std.hash.Crc32.init();
    const member_start = buf.items.len;

    while (true) {
        if (buf.items.len == max_size) {
            var extra: [1]u8 = undefined;
            const n = reader.readSliceShort(&extra) catch {
                logDecodeError(&decompress, ctx);
                return error.GzipReadFailed;
            };
            if (n == 0) break;
            return error.FileTooLarge;
        }

        const room = max_size - buf.items.len;
        const want = @min(CHUNK_SIZE, room);
        try buf.ensureUnusedCapacity(allocator, want);
        const dest = buf.unusedCapacitySlice()[0..want];
        const n = reader.readSliceShort(dest) catch {
            logDecodeError(&decompress, ctx);
            return error.GzipReadFailed;
        };
        if (n == 0) break;
        crc.update(dest[0..n]);
        buf.items.len += n;
    }

    // Verify gzip trailer: std.compress.flate does not check CRC / ISIZE for
    // us. Without this, a corrupt .cif.gz can decompress silently to wrong
    // bytes. See the comment above the Crc32.init() for context.
    const meta = decompress.container_metadata.gzip;
    if (crc.final() != meta.crc) {
        if (ctx.log_errors) {
            std.log.warn("gzip CRC mismatch for {s} (member {d})", .{ ctx.path, ctx.member });
        }
        return error.GzipReadFailed;
    }
    const truncated_size: u32 = @truncate(buf.items.len - member_start);
    if (truncated_size != meta.count) {
        if (ctx.log_errors) {
            std.log.warn("gzip ISIZE mismatch for {s} (member {d})", .{ ctx.path, ctx.member });
        }
        return error.GzipReadFailed;
    }
}

/// The underlying decompressor parks the real cause (BadGzipHeader,
/// InvalidCode, WrongStoredBlockNlen, raw I/O error, etc.) on decompress.err.
/// Surface it via the log so debugging corrupt archives doesn't require a
/// debugger.
fn logDecodeError(decompress: *const std.compress.flate.Decompress, ctx: MemberContext) void {
    if (!ctx.log_errors) return;
    const inner = decompress.err orelse return;
    std.log.warn("gzip decode failed for {s} (member {d}): {s}", .{ ctx.path, ctx.member, @errorName(inner) });
}

// -- Tests --

test "readGzip decompresses gzip store block" {
    const allocator = std.testing.allocator;

    // Minimal gzip containing "Hello world\n" (store block, no compression)
    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    const content = try readGzip(allocator, tmp_path);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\n", content);
}

test "readGzip returns GzipOpenFailed for nonexistent file" {
    const allocator = std.testing.allocator;
    const result = readGzip(allocator, "/nonexistent/path/file.gz");
    try std.testing.expectError(error.GzipOpenFailed, result);
}

test "readGzipLimited accepts exact size limit" {
    const allocator = std.testing.allocator;

    // Minimal gzip containing "Hello world\n" (12 bytes decompressed)
    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    const content = try readGzipLimited(allocator, tmp_path, 12);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\n", content);
}

test "readGzipLimited returns FileTooLarge when limit exceeded" {
    const allocator = std.testing.allocator;

    // Minimal gzip containing "Hello world\n" (12 bytes decompressed)
    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    // Limit to 5 bytes — "Hello world\n" is 12 bytes, should fail
    const result = readGzipLimited(allocator, tmp_path, 5);
    try std.testing.expectError(error.FileTooLarge, result);
}

test "readGzip rejects gzip with corrupted CRC" {
    const allocator = std.testing.allocator;

    // Same "Hello world\n" gzip as above, but with one byte of the CRC32 flipped.
    var gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };
    gz_data[27] ^= 0xff; // corrupt first byte of CRC32 trailer

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    const result = readGzipLimitedOptions(allocator, tmp_path, DEFAULT_MAX_SIZE, .{ .log_errors = false });
    try std.testing.expectError(error.GzipReadFailed, result);
}

// -- Multi-member tests (RFC 1952: a gzip file is a series of members) --

// "Hello world\n" as a stored block, the same member the tests above use.
const test_hello_member = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
    0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
    0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
    0xb7, 0x0c, 0x00, 0x00, 0x00,
};

// `printf 'second member\n' | gzip -nc` (fixed Huffman block).
const test_second_member = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
    0x2b, 0x4e, 0x4d, 0xce, 0xcf, 0x4b, 0x51, 0xc8, 0x4d, 0xcd,
    0x4d, 0x4a, 0x2d, 0xe2, 0x02, 0x00, 0x36, 0x18, 0x4b, 0x0e,
    0x0e, 0x00, 0x00, 0x00,
};

// `printf '' | gzip -nc`: a member with a zero-length payload.
const test_empty_member = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
    0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
};

// BGZF blocks: gzip members carrying the "BC" extra subfield (FEXTRA) that
// holds the block size. Payload: "ATOM      1  N   MET A   1\n".
const test_bgzf_block_1 = [_]u8{
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x2f, 0x00, 0x73, 0x0c,
    0xf1, 0xf7, 0x55, 0x00, 0x03, 0x43, 0x05, 0x05, 0x3f, 0x20,
    0xe5, 0xeb, 0x1a, 0xa2, 0xe0, 0x08, 0xe2, 0x72, 0x01, 0x00,
    0x90, 0xc3, 0xbc, 0x77, 0x1b, 0x00, 0x00, 0x00,
};

// Payload: "ATOM      2  CA  MET A   1\n".
const test_bgzf_block_2 = [_]u8{
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x31, 0x00, 0x73, 0x0c,
    0xf1, 0xf7, 0x55, 0x00, 0x03, 0x23, 0x05, 0x05, 0x67, 0x47,
    0x05, 0x05, 0x5f, 0xd7, 0x10, 0x05, 0x20, 0xa5, 0x60, 0xc8,
    0x05, 0x00, 0x94, 0xeb, 0xb9, 0x40, 0x1b, 0x00, 0x00, 0x00,
};

// The 28-byte empty block that terminates every BGZF file.
const test_bgzf_eof = [_]u8{
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff,
    0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
};

/// Write `gz_data` to a temporary file and decompress it without logging.
/// Caller owns the returned slice.
fn readTestGzip(allocator: std.mem.Allocator, gz_data: []const u8, max_size: usize) ![]u8 {
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    return readGzipLimitedOptions(allocator, tmp_path, max_size, .{ .log_errors = false });
}

test "readGzip concatenates two gzip members" {
    const allocator = std.testing.allocator;

    const content = try readTestGzip(allocator, &(test_hello_member ++ test_second_member), DEFAULT_MAX_SIZE);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\nsecond member\n", content);
}

test "readGzip accepts empty members between and after data members" {
    const allocator = std.testing.allocator;

    const gz_data = test_hello_member ++ test_empty_member ++ test_second_member ++ test_empty_member;
    const content = try readTestGzip(allocator, &gz_data, DEFAULT_MAX_SIZE);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\nsecond member\n", content);
}

test "readGzip accepts a file holding only an empty member" {
    const allocator = std.testing.allocator;

    const content = try readTestGzip(allocator, &test_empty_member, DEFAULT_MAX_SIZE);
    defer allocator.free(content);

    try std.testing.expectEqual(@as(usize, 0), content.len);
}

test "readGzip decodes BGZF blocks followed by the EOF block" {
    const allocator = std.testing.allocator;

    const gz_data = test_bgzf_block_1 ++ test_bgzf_block_2 ++ test_bgzf_eof;
    const content = try readTestGzip(allocator, &gz_data, DEFAULT_MAX_SIZE);
    defer allocator.free(content);

    try std.testing.expectEqualStrings(
        "ATOM      1  N   MET A   1\nATOM      2  CA  MET A   1\n",
        content,
    );
}

test "readGzip rejects garbage after a gzip member" {
    const allocator = std.testing.allocator;

    // Longer than a gzip header, so the bad magic is what gets rejected.
    const long = test_hello_member ++ "this is not a gzip member".*;
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &long, DEFAULT_MAX_SIZE));

    // Shorter than a gzip header.
    const short = test_hello_member ++ "x".*;
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &short, DEFAULT_MAX_SIZE));

    // Garbage after the last of several members.
    const after_two = test_hello_member ++ test_second_member ++ "garbage".*;
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &after_two, DEFAULT_MAX_SIZE));

    // The gzip magic followed by an invalid header (compression method 0).
    const bad_header = test_hello_member ++ [_]u8{ 0x1f, 0x8b } ++ [_]u8{0} ** 16;
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &bad_header, DEFAULT_MAX_SIZE));
}

test "readGzip rejects zero padding after a gzip member" {
    const allocator = std.testing.allocator;

    const one = test_hello_member ++ [_]u8{0};
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &one, DEFAULT_MAX_SIZE));

    const block = test_hello_member ++ [_]u8{0} ** 512;
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &block, DEFAULT_MAX_SIZE));
}

test "readGzip rejects a corrupted CRC in the second member" {
    const allocator = std.testing.allocator;

    var gz_data = test_hello_member ++ test_second_member;
    gz_data[gz_data.len - 8] ^= 0xff; // first byte of the second member's CRC32
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &gz_data, DEFAULT_MAX_SIZE));
}

test "readGzip rejects a wrong ISIZE in the second member" {
    const allocator = std.testing.allocator;

    var gz_data = test_hello_member ++ test_second_member;
    gz_data[gz_data.len - 4] ^= 0xff; // first byte of the second member's ISIZE
    try std.testing.expectError(error.GzipReadFailed, readTestGzip(allocator, &gz_data, DEFAULT_MAX_SIZE));
}

test "readGzip rejects a truncated second member" {
    const allocator = std.testing.allocator;

    const gz_data = test_hello_member ++ test_second_member;
    // Cut inside the second member's trailer, its deflate data and its header.
    for ([_]usize{ 1, 4, 8, 12, test_second_member.len - 5, test_second_member.len - 1 }) |cut| {
        const result = readTestGzip(allocator, gz_data[0 .. gz_data.len - cut], DEFAULT_MAX_SIZE);
        try std.testing.expectError(error.GzipReadFailed, result);
    }
}

test "readGzipLimited applies the size limit to the total across members" {
    const allocator = std.testing.allocator;

    // 12 + 14 decompressed bytes.
    const gz_data = test_hello_member ++ test_second_member;

    const content = try readTestGzip(allocator, &gz_data, 26);
    defer allocator.free(content);
    try std.testing.expectEqualStrings("Hello world\nsecond member\n", content);

    // Each member fits on its own, the total does not.
    try std.testing.expectError(error.FileTooLarge, readTestGzip(allocator, &gz_data, 25));
    try std.testing.expectError(error.FileTooLarge, readTestGzip(allocator, &gz_data, 14));
    // The first member fills the limit exactly; the second must not be appended.
    try std.testing.expectError(error.FileTooLarge, readTestGzip(allocator, &gz_data, 12));
}

test "readGzipLimited accepts an empty member once the limit is reached" {
    const allocator = std.testing.allocator;

    const gz_data = test_hello_member ++ test_bgzf_eof;
    const content = try readTestGzip(allocator, &gz_data, 12);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\n", content);
}

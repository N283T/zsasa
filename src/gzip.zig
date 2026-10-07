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

    // The file reader is unbuffered: PaddedInput reads from it straight into
    // its own buffer, which is the one the decompressor consumes from.
    var file_reader = file.reader(io, &.{});
    var input_buf: [CHUNK_SIZE]u8 = undefined;
    var input: PaddedInput = .init(&file_reader.interface, &input_buf);

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
        try readMember(allocator, &input, &window_buf, &buf, max_size, ctx);
        if (!try hasNextMember(&input, ctx)) break;
    }

    return buf.toOwnedSlice(allocator);
}

/// The compressed file as the decompressor sees it: the file's bytes followed
/// by a few zero bytes of padding, with a record of whether any of the
/// padding was consumed. The end of the file is never reported as
/// EndOfStream; reading on after consuming padding fails with ReadFailed.
///
/// Workaround for std.compress.flate.Decompress in Zig 0.16.0, which crashes
/// (integer overflow panic in safe builds, out-of-bounds read otherwise) when
/// its input ends inside a dynamic Huffman block: `tossBitsShort` adds
/// `consumed_bits` to the number of buffered bits where it should subtract it,
/// so it tosses bits that are not there, and the next `peekBitsEnding` then
/// underflows. Those code paths are only taken once the input reports
/// EndOfStream, so this input never does. A decoder that runs past the end of
/// the file gets padding instead, and `pastEnd` tells the caller that the
/// file is truncated.
///
/// Fixed upstream after 0.16.0 (ziglang/zig commit 12815e2223, "flate: Correct
/// math when calculating number of available bits"). Once zsasa requires a
/// Zig release with that fix, this type can go and the decompressor can read
/// the file reader directly again.
const PaddedInput = struct {
    source: *std.Io.Reader,
    /// Set once `source` has reported the end of the file.
    source_ended: bool = false,
    /// Number of padding bytes supplied so far. Padding starts at the end of
    /// the file, so these are always the last bytes that were supplied.
    padding_len: usize = 0,
    interface: std.Io.Reader,

    /// Padding bytes supplied per read once the file is exhausted. More than
    /// the decompressor asks for in one go (the 10-byte gzip header), so it
    /// cannot ask twice without having consumed some of them.
    const padding_chunk = 64;

    /// `buffer` is what the decompressor peeks into: flate.Decompress asserts
    /// that it holds at least 10 bytes.
    fn init(source: *std.Io.Reader, buffer: []u8) PaddedInput {
        return .{
            .source = source,
            .interface = .{
                .vtable = &.{ .stream = stream },
                .buffer = buffer,
                .seek = 0,
                .end = 0,
            },
        };
    }

    fn stream(r: *std.Io.Reader, w: *std.Io.Writer, limit: std.Io.Limit) std.Io.Reader.StreamError!usize {
        const self: *PaddedInput = @alignCast(@fieldParentPtr("interface", r));
        const dest = limit.slice(try w.writableSliceGreedy(1));
        if (dest.len == 0) return 0;

        if (!self.source_ended) {
            // readSliceShort fills `dest` unless the file ends first.
            const n = try self.source.readSliceShort(dest);
            if (n < dest.len) self.source_ended = true;
            if (n != 0) {
                w.advance(n);
                return n;
            }
        }

        // The file is known to be truncated once padding has been consumed.
        // Stop there instead of feeding the decoder zeros for as long as it
        // manages to decode them.
        if (self.pastEnd()) return error.ReadFailed;

        const padding = dest[0..@min(dest.len, padding_chunk)];
        @memset(padding, 0);
        w.advance(padding.len);
        self.padding_len += padding.len;
        return padding.len;
    }

    /// True once bytes from beyond the end of the file have been consumed.
    /// Padding that is merely buffered (peeked at, not consumed) does not
    /// count.
    fn pastEnd(self: *const PaddedInput) bool {
        return self.padding_len > self.interface.bufferedLen();
    }

    /// The next `n` bytes of the file without consuming them, or all that is
    /// left when the file ends before that. Must not be called once `pastEnd`.
    fn peekFile(self: *PaddedInput, n: usize) error{ReadFailed}![]const u8 {
        // Never EndOfStream: padding makes up for the bytes the file lacks.
        _ = self.interface.peek(n) catch return error.ReadFailed;
        const buffered = self.interface.buffered();
        const file_bytes = buffered[0 .. buffered.len - @min(buffered.len, self.padding_len)];
        return file_bytes[0..@min(n, file_bytes.len)];
    }
};

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
fn hasNextMember(input: *PaddedInput, ctx: MemberContext) error{GzipReadFailed}!bool {
    // Shorter than the magic when fewer bytes than that are left in the file.
    const next = input.peekFile(gzip_magic.len) catch {
        if (ctx.log_errors) {
            std.log.warn("gzip read failed for {s} after member {d}", .{ ctx.path, ctx.member });
        }
        return error.GzipReadFailed;
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
    input: *PaddedInput,
    window_buf: []u8,
    buf: *std.ArrayListUnmanaged(u8),
    max_size: usize,
    ctx: MemberContext,
) GzipError!void {
    var decompress: std.compress.flate.Decompress = .init(&input.interface, .gzip, window_buf);

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
            const n = try readDecoded(&decompress, input, &extra, ctx);
            if (n == 0) break;
            return error.FileTooLarge;
        }

        const room = max_size - buf.items.len;
        const want = @min(CHUNK_SIZE, room);
        try buf.ensureUnusedCapacity(allocator, want);
        const dest = buf.unusedCapacitySlice()[0..want];
        const n = try readDecoded(&decompress, input, dest, ctx);
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

/// Read the next decompressed bytes of the current member into `dest`.
/// Returns how many were read; 0 means the member, trailer included, has been
/// decoded completely. Fails if the decompressor reports an error or if it
/// needed bytes from beyond the end of the file, i.e. the file is truncated.
fn readDecoded(
    decompress: *std.compress.flate.Decompress,
    input: *const PaddedInput,
    dest: []u8,
    ctx: MemberContext,
) error{GzipReadFailed}!usize {
    const n = decompress.reader.readSliceShort(dest) catch {
        logDecodeError(decompress, input, ctx);
        return error.GzipReadFailed;
    };
    // Checked on every read, also a successful one: bytes decoded from the
    // padding must not reach the caller, and a member whose trailer was read
    // from the padding ends without an error.
    if (input.pastEnd()) {
        logDecodeError(decompress, input, ctx);
        return error.GzipReadFailed;
    }
    return n;
}

/// The underlying decompressor parks the real cause (BadGzipHeader,
/// InvalidCode, WrongStoredBlockNlen, raw I/O error, etc.) on decompress.err.
/// Surface it via the log so debugging corrupt archives doesn't require a
/// debugger. Whatever it reports after running past the end of the file is an
/// artifact of decoding the padding, so that case is reported as truncation.
fn logDecodeError(decompress: *const std.compress.flate.Decompress, input: *const PaddedInput, ctx: MemberContext) void {
    if (!ctx.log_errors) return;
    if (input.pastEnd()) {
        std.log.warn("gzip decode failed for {s} (member {d}): unexpected end of file", .{ ctx.path, ctx.member });
        return;
    }
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

// -- Truncated input tests (#442) --

// `grep '^ATOM' examples/1ubq.pdb | head -n 32 | gzip -nc`: 2592 bytes of ATOM
// records in one dynamic Huffman block. The stored and fixed Huffman members
// above do not reach the bit reader code that mishandles the end of the input.
const test_dynamic_member = [_]u8{
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
    0x7d, 0x96, 0x4b, 0x72, 0xdb, 0x30, 0x10, 0x44, 0xf7, 0x39,
    0x05, 0x4f, 0x30, 0x85, 0xf9, 0xe0, 0xb7, 0x54, 0x64, 0x95,
    0x9c, 0x2a, 0xdb, 0x4c, 0x95, 0x7d, 0xff, 0xb3, 0x64, 0x66,
    0x20, 0x0a, 0xa4, 0x40, 0x44, 0x1b, 0x9a, 0x92, 0xf9, 0xd8,
    0x98, 0xee, 0x06, 0x79, 0xf9, 0x59, 0x3f, 0x17, 0xff, 0xe0,
    0xb2, 0x7c, 0xe9, 0xe1, 0xf3, 0xf6, 0xb3, 0x5c, 0xda, 0xa9,
    0x7d, 0x28, 0x03, 0x4b, 0xd0, 0xa3, 0x80, 0x70, 0xb0, 0x2f,
    0x20, 0xa1, 0xe8, 0xcf, 0x10, 0xf4, 0xac, 0x42, 0xca, 0x4b,
    0xff, 0xe8, 0xf5, 0xbf, 0x2e, 0x4f, 0x20, 0x2d, 0xcb, 0xf5,
    0x32, 0x02, 0x13, 0x50, 0x4a, 0x7a, 0x8c, 0x20, 0xc8, 0x0e,
    0x2c, 0x42, 0x0f, 0x20, 0x06, 0xe0, 0xb2, 0x03, 0x5e, 0x0f,
    0x40, 0xf6, 0x2f, 0x4e, 0x80, 0xd5, 0x40, 0x7a, 0x4c, 0x5c,
    0xed, 0xdf, 0x20, 0x32, 0xee, 0x14, 0xd2, 0x14, 0xa8, 0xeb,
    0x58, 0x4f, 0x97, 0x5c, 0x4a, 0x72, 0xa0, 0x24, 0x53, 0x28,
    0xaa, 0x98, 0x27, 0xc0, 0xf5, 0x00, 0x8c, 0x7a, 0x87, 0xdf,
    0x23, 0x30, 0x02, 0x22, 0xf9, 0x0c, 0x4b, 0x09, 0xae, 0x30,
    0x49, 0xdd, 0x96, 0xcc, 0x90, 0xf7, 0x33, 0x3c, 0x2a, 0x54,
    0x19, 0xd7, 0xfb, 0x19, 0x90, 0x23, 0x37, 0x60, 0x32, 0xa0,
    0xde, 0x80, 0x37, 0x53, 0x50, 0x07, 0x5c, 0xa7, 0x40, 0xbd,
    0xd5, 0xf7, 0xdb, 0x08, 0x64, 0xa8, 0xe6, 0xae, 0x1d, 0x63,
    0x75, 0x60, 0x0d, 0x4f, 0x60, 0x06, 0xdc, 0x2b, 0xfc, 0x3e,
    0x00, 0xd5, 0xae, 0xeb, 0x6d, 0x04, 0x6a, 0x5c, 0x24, 0x37,
    0x60, 0x11, 0xbb, 0xaf, 0x8e, 0x2d, 0x74, 0x85, 0x88, 0x53,
    0x85, 0xb5, 0xe5, 0xf0, 0xfe, 0xf1, 0xe5, 0x40, 0x7a, 0xba,
    0xcc, 0x1c, 0xdd, 0x9c, 0x9c, 0xdb, 0x0c, 0x29, 0x96, 0x6e,
    0x0a, 0x4d, 0x73, 0x88, 0xa1, 0xe5, 0x70, 0x04, 0x96, 0x68,
    0x4b, 0xae, 0x10, 0x08, 0x1d, 0x58, 0xea, 0x0e, 0x18, 0xa6,
    0xa6, 0x98, 0xf6, 0xeb, 0xa9, 0x42, 0x0c, 0x0d, 0x48, 0x66,
    0x8e, 0xce, 0x90, 0xc2, 0x16, 0xec, 0xa5, 0x40, 0x9e, 0xe6,
    0xd0, 0xb2, 0xb1, 0x9e, 0x00, 0xcd, 0xdd, 0xf8, 0x50, 0x28,
    0x0e, 0x64, 0x0e, 0x1d, 0x48, 0xd3, 0x1c, 0x5a, 0x1d, 0x2c,
    0x87, 0xa3, 0xc2, 0xcc, 0xfa, 0x13, 0x07, 0x40, 0xb1, 0x9a,
    0x91, 0xba, 0x1c, 0x37, 0x53, 0xd4, 0xb1, 0x34, 0x55, 0x28,
    0x2d, 0x87, 0x27, 0x33, 0x2c, 0xfa, 0x37, 0x23, 0x44, 0xbf,
    0x98, 0x41, 0x42, 0xed, 0xb1, 0x09, 0x53, 0x97, 0xd1, 0x9a,
    0xf2, 0x76, 0xaa, 0xd0, 0xaa, 0xc7, 0x04, 0xd1, 0x6b, 0x46,
    0x6a, 0xec, 0xb6, 0x64, 0x52, 0xd9, 0x61, 0x0a, 0xd4, 0xab,
    0xd6, 0x1b, 0x0e, 0x40, 0x8d, 0x4b, 0xb1, 0x25, 0x33, 0xa0,
    0x37, 0x05, 0xa1, 0x74, 0xa0, 0x9e, 0xec, 0x9b, 0x72, 0x9c,
    0xa1, 0xfa, 0xff, 0x75, 0xa3, 0x01, 0x18, 0x9b, 0x32, 0x55,
    0xe8, 0xb3, 0x34, 0x60, 0x48, 0xdb, 0x92, 0x2b, 0xc8, 0x1e,
    0x78, 0xcc, 0x61, 0x69, 0xc1, 0xfe, 0xf3, 0x71, 0x73, 0x20,
    0xf7, 0x19, 0xda, 0x45, 0xea, 0x72, 0x8a, 0x36, 0x43, 0xad,
    0xaf, 0xdd, 0xbb, 0xb9, 0x1c, 0x55, 0xee, 0x14, 0x58, 0x5b,
    0xb0, 0x47, 0x20, 0x59, 0x53, 0xd4, 0xe5, 0x60, 0x0d, 0xd1,
    0xea, 0x49, 0xdd, 0x01, 0xe7, 0xc1, 0xb6, 0x82, 0x5e, 0xcf,
    0x15, 0x3e, 0x5c, 0x16, 0x6a, 0xc0, 0x92, 0x9e, 0xc1, 0x96,
    0xff, 0xb8, 0x6c, 0xbd, 0x5a, 0x4f, 0x80, 0x59, 0x73, 0x97,
    0x1c, 0x98, 0x7d, 0x23, 0xc8, 0xba, 0xc1, 0x4a, 0x07, 0xa6,
    0x3d, 0xf0, 0x60, 0x8a, 0x45, 0xde, 0x82, 0x3d, 0x2a, 0x64,
    0x91, 0xd6, 0x14, 0xeb, 0xb4, 0x96, 0x23, 0xc9, 0x16, 0x6c,
    0x1d, 0x68, 0x8c, 0x53, 0x85, 0xd6, 0x94, 0x3b, 0x9e, 0x28,
    0x2c, 0x96, 0x35, 0xd2, 0xda, 0x7a, 0x53, 0x0a, 0xd4, 0x5a,
    0xbb, 0xc2, 0x79, 0x97, 0xc9, 0x9b, 0x42, 0x03, 0x50, 0x9f,
    0x76, 0x15, 0x1f, 0xdb, 0x17, 0x3a, 0x90, 0xca, 0xce, 0x94,
    0x38, 0x7d, 0xea, 0x91, 0x37, 0xe5, 0x4c, 0x61, 0xb5, 0x67,
    0xaf, 0x2a, 0x0c, 0x0e, 0x0a, 0xfa, 0x38, 0xdd, 0x80, 0x7a,
    0x52, 0x78, 0x0a, 0x4c, 0x2d, 0x87, 0x7f, 0xdf, 0x1b, 0x50,
    0x7a, 0x6c, 0xac, 0xe6, 0x1a, 0xec, 0x60, 0x71, 0xb1, 0x1d,
    0x2b, 0x63, 0x5f, 0xf2, 0x61, 0x86, 0x87, 0x1c, 0xda, 0xd6,
    0x6b, 0x39, 0x1c, 0x81, 0xd9, 0xe6, 0xa4, 0xd5, 0x13, 0xb6,
    0x60, 0x57, 0xc0, 0x9e, 0x43, 0x75, 0x79, 0xbe, 0xe4, 0xd2,
    0x72, 0xf8, 0x0a, 0xd4, 0xe7, 0x50, 0x44, 0x07, 0xb2, 0xe7,
    0x2f, 0x68, 0x63, 0x42, 0x9f, 0x21, 0x4f, 0x37, 0x07, 0x7b,
    0x1e, 0xae, 0xa7, 0x4b, 0x66, 0x03, 0x58, 0x97, 0xb3, 0x6d,
    0xfd, 0x08, 0x5c, 0xe3, 0xc4, 0x94, 0x43, 0x0e, 0xed, 0x56,
    0x96, 0xc3, 0x57, 0x60, 0x84, 0x64, 0x00, 0xd6, 0xbd, 0xb4,
    0x3e, 0x62, 0x23, 0xa9, 0x2f, 0x79, 0x6e, 0x8a, 0xbd, 0xb0,
    0xd8, 0x06, 0x3b, 0x02, 0xa9, 0x14, 0x07, 0x26, 0xdb, 0x58,
    0xbd, 0xcb, 0xcf, 0x1c, 0x6a, 0x04, 0xa6, 0x5d, 0x66, 0x6a,
    0xb1, 0x79, 0x05, 0x8a, 0x6e, 0xfd, 0xd9, 0x67, 0x58, 0xed,
    0x2d, 0x4c, 0x19, 0x81, 0x4b, 0x6f, 0x4a, 0x7a, 0x7d, 0x73,
    0xf8, 0x07, 0xa3, 0x2b, 0xfc, 0x96, 0x20, 0x0a, 0x00, 0x00,
};
const test_dynamic_size = 2592;
const test_dynamic_first_line = "ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  9.67           N  \n";
const test_dynamic_last_line = "ATOM     32  CD1 PHE A   4      24.147  33.966   7.038  1.00  6.69           C  \n";

/// Decompress every proper prefix of `gz_data` (1 to `gz_data.len - 1` bytes)
/// and require each one to fail cleanly with GzipReadFailed. The exception is
/// a prefix of `complete_len` bytes: it ends on a member boundary, so it is a
/// complete gzip file and must decode to `complete_size` bytes.
fn expectTruncationsRejected(gz_data: []const u8, complete_len: ?usize, complete_size: usize) !void {
    const allocator = std.testing.allocator;

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    // Shrink the one file a byte at a time; much cheaper than rewriting it.
    const file = try tmp_dir.dir.openFile(std.testing.io, "test.gz", .{ .mode = .write_only });
    defer file.close(std.testing.io);

    var len = gz_data.len - 1;
    while (len > 0) : (len -= 1) {
        errdefer std.debug.print("gzip input truncated to {d} of {d} bytes\n", .{ len, gz_data.len });
        try file.setLength(std.testing.io, len);

        const result = readGzipLimitedOptions(allocator, tmp_path, DEFAULT_MAX_SIZE, .{ .log_errors = false });
        if (result) |content| {
            defer allocator.free(content);
            try std.testing.expectEqual(complete_len, len);
            try std.testing.expectEqual(complete_size, content.len);
        } else |err| {
            try std.testing.expectEqual(error.GzipReadFailed, err);
            try std.testing.expect(complete_len != len);
        }
    }
}

test "readGzip decompresses a dynamic Huffman block" {
    const allocator = std.testing.allocator;

    const content = try readTestGzip(allocator, &test_dynamic_member, DEFAULT_MAX_SIZE);
    defer allocator.free(content);

    try std.testing.expectEqual(@as(usize, test_dynamic_size), content.len);
    try std.testing.expectStringStartsWith(content, test_dynamic_first_line);
    try std.testing.expectStringEndsWith(content, test_dynamic_last_line);
}

test "readGzip rejects a dynamic Huffman member truncated at any length" {
    try expectTruncationsRejected(&test_dynamic_member, null, 0);
}

test "readGzip rejects a two-member file truncated at any length" {
    // Every cut but one falls inside a member. Cut exactly between the two,
    // the file is the complete first member ("Hello world\n").
    const gz_data = test_hello_member ++ test_dynamic_member;
    try expectTruncationsRejected(&gz_data, test_hello_member.len, 12);
}

test "readGzip decodes members that straddle the input buffer" {
    const allocator = std.testing.allocator;

    // More compressed data than the input buffer holds, so that members, and
    // the look-ahead for the next member, span refills of that buffer.
    const count = 100;
    const gz_data = test_dynamic_member ** count;
    comptime std.debug.assert(gz_data.len > CHUNK_SIZE);

    const content = try readTestGzip(allocator, &gz_data, DEFAULT_MAX_SIZE);
    defer allocator.free(content);

    try std.testing.expectEqual(@as(usize, count * test_dynamic_size), content.len);
    try std.testing.expectStringStartsWith(content, test_dynamic_first_line);
    try std.testing.expectStringEndsWith(content, test_dynamic_last_line);

    // Truncated inside the last member, after the buffer has been refilled.
    const result = readTestGzip(allocator, gz_data[0 .. gz_data.len - 300], DEFAULT_MAX_SIZE);
    try std.testing.expectError(error.GzipReadFailed, result);
}

/// Wrap `data` in a gzip member made of stored (uncompressed) blocks.
/// Caller owns the returned slice.
fn storedGzipMember(allocator: std.mem.Allocator, data: []const u8) ![]u8 {
    var out: std.ArrayListUnmanaged(u8) = .empty;
    errdefer out.deinit(allocator);

    try out.appendSlice(allocator, test_hello_member[0..10]); // gzip header
    var rest = data;
    while (true) {
        const len: u16 = @intCast(@min(rest.len, 0xffff));
        const is_final = rest.len == len;
        var block_header: [5]u8 = undefined;
        block_header[0] = @intFromBool(is_final); // BFINAL, BTYPE = 00 (stored)
        std.mem.writeInt(u16, block_header[1..3], len, .little);
        std.mem.writeInt(u16, block_header[3..5], ~len, .little);
        try out.appendSlice(allocator, &block_header);
        try out.appendSlice(allocator, rest[0..len]);
        rest = rest[len..];
        if (is_final) break;
    }
    var trailer: [8]u8 = undefined;
    std.mem.writeInt(u32, trailer[0..4], std.hash.Crc32.hash(data), .little);
    std.mem.writeInt(u32, trailer[4..8], @truncate(data.len), .little);
    try out.appendSlice(allocator, &trailer);

    return out.toOwnedSlice(allocator);
}

test "readGzip reads stored blocks larger than the input buffer" {
    const allocator = std.testing.allocator;

    // The decompressor copies stored blocks from the input straight into its
    // window instead of going through the input buffer.
    const data = try allocator.alloc(u8, 3 * CHUNK_SIZE - 123);
    defer allocator.free(data);
    for (data, 0..) |*byte, i| byte.* = @truncate(i *% 31 +% (i >> 8));

    const gz_data = try storedGzipMember(allocator, data);
    defer allocator.free(gz_data);

    const content = try readTestGzip(allocator, gz_data, DEFAULT_MAX_SIZE);
    defer allocator.free(content);
    try std.testing.expect(std.mem.eql(u8, data, content));

    // Cut inside each of the three blocks, at the input buffer size and
    // inside the trailer.
    const cuts = [_]usize{ 100, CHUNK_SIZE, CHUNK_SIZE + 1, 2 * CHUNK_SIZE + 17, gz_data.len - 9, gz_data.len - 1 };
    for (cuts) |len| {
        const result = readTestGzip(allocator, gz_data[0..len], DEFAULT_MAX_SIZE);
        try std.testing.expectError(error.GzipReadFailed, result);
    }
}

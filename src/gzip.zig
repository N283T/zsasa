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

/// Decompress a gzip file with a custom size limit. Caller owns the returned slice.
/// The size limit is checked before each chunk read so a malicious file cannot
/// allocate more than `max_size` bytes before the cap is detected.
pub fn readGzipLimited(allocator: std.mem.Allocator, path: []const u8, max_size: usize) GzipError![]u8 {
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

    // Native gzip decompressor. The window buffer must be at least
    // flate.max_window_len (64 KB) to support generic reads that go through
    // the indirect (buffered) vtable used by readSliceShort/readVec.
    // 64 KB on stack — required by the indirect (buffered) Decompress vtable
    // path used by readSliceShort/readVec; see flate.history_len.
    var window_buf: [std.compress.flate.max_window_len]u8 = undefined;
    var decompress: std.compress.flate.Decompress = .init(&file_reader.interface, .gzip, &window_buf);
    const reader: *std.Io.Reader = &decompress.reader;

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    errdefer buf.deinit(allocator);

    // Compute CRC32 over the decompressed bytes so we can verify the gzip
    // trailer ourselves. std.compress.flate parses the trailer fields into
    // decompress.container_metadata.gzip.{crc,count} but does NOT compare them
    // against the actual decoded bytes. Without this verification, a corrupt
    // .cif.gz with a valid frame structure but a bad checksum would silently
    // return wrong bytes — a real correctness regression for the SASA pipeline
    // (the previous C-zlib version verified CRC via gzclose).
    var crc = std.hash.Crc32.init();

    while (true) {
        if (buf.items.len >= max_size) return error.FileTooLarge;
        const room = max_size - buf.items.len;
        const want = @min(CHUNK_SIZE, room);
        try buf.ensureUnusedCapacity(allocator, want);
        const dest = buf.unusedCapacitySlice()[0..want];
        const n = reader.readSliceShort(dest) catch {
            // The underlying decompressor parks the real cause (BadGzipHeader,
            // InvalidCode, WrongStoredBlockNlen, raw I/O error, etc.) on
            // decompress.err. Surface it via the log so debugging corrupt
            // archives doesn't require a debugger.
            if (decompress.err) |inner| {
                std.log.warn("gzip decode failed for {s}: {s}", .{ path, @errorName(inner) });
            }
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
        std.log.warn("gzip CRC mismatch for {s}", .{path});
        return error.GzipReadFailed;
    }
    const truncated_size: u32 = @truncate(buf.items.len);
    if (truncated_size != meta.count) {
        std.log.warn("gzip ISIZE mismatch for {s}", .{path});
        return error.GzipReadFailed;
    }

    return buf.toOwnedSlice(allocator);
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

    const result = readGzip(allocator, tmp_path);
    try std.testing.expectError(error.GzipReadFailed, result);
}

//! Gzip decompression using C zlib.
//!
//! Workaround for Zig 0.15 flate bug (unreachable in Writer.rebase).
//! See: https://github.com/ziglang/zig/issues/25035
//! TODO: revert to native flate when Zig fixes the bug (test with 2oxd.cif.gz)

const std = @import("std");
const c = @import("zlib_c");

pub const GzipError = error{ GzipOpenFailed, GzipReadFailed, FileTooLarge, OutOfMemory };

const CHUNK_SIZE = 64 * 1024; // 64 KB read chunks

/// Default max decompressed size: 4 GB
pub const DEFAULT_MAX_SIZE: usize = 4 * 1024 * 1024 * 1024;

/// Decompress a gzip file using C zlib. Caller owns the returned slice.
pub fn readGzip(allocator: std.mem.Allocator, path: []const u8) GzipError![]u8 {
    return readGzipLimited(allocator, path, DEFAULT_MAX_SIZE);
}

/// Decompress a gzip file with a custom size limit. Caller owns the returned slice.
pub fn readGzipLimited(allocator: std.mem.Allocator, path: []const u8, max_size: usize) GzipError![]u8 {
    const c_path = try allocator.dupeZ(u8, path);
    defer allocator.free(c_path);

    const gz = c.gzopen(c_path.ptr, "rb") orelse return error.GzipOpenFailed;
    var gz_closed = false;
    errdefer if (!gz_closed) {
        _ = c.gzclose(gz);
    };

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    errdefer buf.deinit(allocator);

    while (true) {
        if (buf.items.len >= max_size) return error.FileTooLarge;
        const to_read = @min(CHUNK_SIZE, max_size - buf.items.len);
        try buf.ensureUnusedCapacity(allocator, to_read);
        const dest = buf.unusedCapacitySlice()[0..to_read];
        const n = c.gzread(gz, dest.ptr, @intCast(dest.len));
        if (n < 0) return error.GzipReadFailed;
        if (n == 0) break;
        buf.items.len += @intCast(n);
    }

    // gzclose checks CRC — detect corrupt data before returning
    const close_result = c.gzclose(gz);
    gz_closed = true;
    if (close_result != c.Z_OK) {
        buf.deinit(allocator);
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

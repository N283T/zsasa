//! Zstandard decompression using native std.compress.zstd.

const std = @import("std");

pub const ZstdError = error{ ZstdOpenFailed, ZstdReadFailed, FileTooLarge, OutOfMemory };

const CHUNK_SIZE = 64 * 1024;

/// Default max decompressed size: 4 GB.
pub const DEFAULT_MAX_SIZE: usize = 4 * 1024 * 1024 * 1024;

/// Decompress a zstd file. Caller owns the returned slice.
pub fn readZstd(allocator: std.mem.Allocator, path: []const u8) ZstdError![]u8 {
    return readZstdLimited(allocator, path, DEFAULT_MAX_SIZE);
}

/// Decompress a zstd file with a custom size limit. Caller owns the returned slice.
pub fn readZstdLimited(allocator: std.mem.Allocator, path: []const u8, max_size: usize) ZstdError![]u8 {
    var threaded: std.Io.Threaded = .init_single_threaded;
    const io = threaded.io();

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch return error.ZstdOpenFailed;
    defer file.close(io);

    var file_buf: [CHUNK_SIZE]u8 = undefined;
    var file_reader = file.reader(io, &file_buf);
    return readZstdFromReader(allocator, path, &file_reader.interface, max_size);
}

fn readZstdFromReader(allocator: std.mem.Allocator, path: []const u8, input: *std.Io.Reader, max_size: usize) ZstdError![]u8 {
    const window_len = std.compress.zstd.default_window_len;
    const zstd_buf = try allocator.alloc(u8, window_len + std.compress.zstd.block_size_max);
    defer allocator.free(zstd_buf);

    var decompress: std.compress.zstd.Decompress = .init(input, zstd_buf, .{
        .window_len = window_len,
    });
    const reader: *std.Io.Reader = &decompress.reader;

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    errdefer buf.deinit(allocator);

    while (true) {
        if (buf.items.len >= max_size) return error.FileTooLarge;
        const room = max_size - buf.items.len;
        const want = @min(CHUNK_SIZE, room);
        try buf.ensureUnusedCapacity(allocator, want);
        const dest = buf.unusedCapacitySlice()[0..want];
        const n = reader.readSliceShort(dest) catch {
            if (decompress.err) |inner| {
                std.log.warn("zstd decode failed for {s}: {s}", .{ path, @errorName(inner) });
            }
            return error.ZstdReadFailed;
        };
        if (n == 0) break;
        buf.items.len += n;
    }

    return buf.toOwnedSlice(allocator);
}

test "readZstd decompresses zstd frame" {
    const allocator = std.testing.allocator;

    const zst_data = [_]u8{
        0x28, 0xb5, 0x2f, 0xfd, 0x24, 0x0c, 0x61, 0x00, 0x00, 0x48, 0x65, 0x6c,
        0x6c, 0x6f, 0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0x8c, 0xa0, 0x38,
        0x9a,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.zst", .data = &zst_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.zst", allocator);
    defer allocator.free(tmp_path);

    const content = try readZstd(allocator, tmp_path);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\n", content);
}

test "readZstdLimited returns FileTooLarge when limit exceeded" {
    const allocator = std.testing.allocator;

    const zst_data = [_]u8{
        0x28, 0xb5, 0x2f, 0xfd, 0x24, 0x0c, 0x61, 0x00, 0x00, 0x48, 0x65, 0x6c,
        0x6c, 0x6f, 0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0x8c, 0xa0, 0x38,
        0x9a,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.zst", .data = &zst_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.zst", allocator);
    defer allocator.free(tmp_path);

    const result = readZstdLimited(allocator, tmp_path, 5);
    try std.testing.expectError(error.FileTooLarge, result);
}

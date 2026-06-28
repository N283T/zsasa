//! Shared helpers for transparent compressed input formats.

const std = @import("std");
const gzip = @import("gzip.zig");
const zstd = @import("zstd.zig");

pub const ReadError = gzip.GzipError || zstd.ZstdError || error{NotCompressed};

pub fn isGzip(path: []const u8) bool {
    return std.ascii.endsWithIgnoreCase(path, ".gz");
}

pub fn isZstd(path: []const u8) bool {
    return std.ascii.endsWithIgnoreCase(path, ".zst");
}

pub fn isCompressed(path: []const u8) bool {
    return isGzip(path) or isZstd(path);
}

/// Read a supported compressed file. Caller owns the returned slice.
pub fn read(allocator: std.mem.Allocator, path: []const u8) ReadError![]u8 {
    if (isGzip(path)) return gzip.readGzip(allocator, path);
    if (isZstd(path)) return zstd.readZstd(allocator, path);
    return error.NotCompressed;
}

test "isCompressed recognizes gzip and zstd" {
    try std.testing.expect(isCompressed("file.cif.gz"));
    try std.testing.expect(isCompressed("file.cif.GZ"));
    try std.testing.expect(isCompressed("file.cif.zst"));
    try std.testing.expect(isCompressed("file.cif.ZST"));
    try std.testing.expect(!isCompressed("file.cif"));
}

test "read accepts uppercase gzip extension" {
    const allocator = std.testing.allocator;

    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.PDB.GZ", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.PDB.GZ", allocator);
    defer allocator.free(tmp_path);

    const content = try read(allocator, tmp_path);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\n", content);
}

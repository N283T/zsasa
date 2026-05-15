//! Shared helpers for transparent compressed input formats.

const std = @import("std");
const gzip = @import("gzip.zig");
const zstd = @import("zstd.zig");

pub const ReadError = gzip.GzipError || zstd.ZstdError || error{NotCompressed};

pub fn isGzip(path: []const u8) bool {
    return std.mem.endsWith(u8, path, ".gz");
}

pub fn isZstd(path: []const u8) bool {
    return std.mem.endsWith(u8, path, ".zst");
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
    try std.testing.expect(isCompressed("file.cif.zst"));
    try std.testing.expect(!isCompressed("file.cif"));
}

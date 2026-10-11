const std = @import("std");
const builtin = @import("builtin");

pub const is_windows = builtin.os.tag == .windows;

/// A read-only file region backed by mmap (POSIX) or heap allocation (Windows).
///
/// Ownership: single-owner, non-copyable by convention. Exactly one `deinit`
/// call must occur per `mmapFile` call.
/// Typical usage: `const mapped = try mmapFile(alloc, path); defer mapped.deinit();`
pub const MappedFile = struct {
    data: if (is_windows) []const u8 else []align(std.heap.page_size_min) const u8,
    allocator: if (is_windows) std.mem.Allocator else void,

    pub fn deinit(self: MappedFile) void {
        if (self.data.len == 0) return;
        if (is_windows) {
            self.allocator.free(self.data);
        } else {
            std.posix.munmap(@alignCast(self.data));
        }
    }
};

pub fn mmapFile(allocator: std.mem.Allocator, io: std.Io, path: []const u8) !MappedFile {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    const stat = try file.stat(io);
    const size: usize = std.math.cast(usize, stat.size) orelse return error.FileTooBig;
    if (size == 0) return .{ .data = &.{}, .allocator = if (is_windows) allocator else {} };

    if (is_windows) {
        var read_buf: [65536]u8 = undefined;
        var r = file.reader(io, &read_buf);
        // Use .unlimited rather than .limited64(size). On Windows the file
        // size reported by stat occasionally diverges from the byte count
        // produced by the buffered reader (text-mode line-ending handling
        // and similar quirks), causing spurious StreamTooLong errors on
        // benign reads of input fixtures (#352 / v0.2.10 publish failure).
        // For trusted local input files this is acceptable; callers of
        // mmapFile already trust the path.
        const data = try r.interface.allocRemaining(allocator, .unlimited);
        return .{ .data = data, .allocator = allocator };
    } else {
        const mapped = try std.posix.mmap(
            null,
            size,
            .{ .READ = true },
            .{ .TYPE = .PRIVATE },
            file.handle,
            0,
        );
        return .{ .data = mapped, .allocator = {} };
    }
}

/// Reads a whole file into memory, where possible with a single read of its
/// stat size. Caller frees the result with `allocator`.
///
/// One read of the whole file is one disk request for a file that is not
/// cached, where page faults on a mapping or a loop of small reads issue
/// several.
pub fn readFileOnce(allocator: std.mem.Allocator, io: std.Io, path: []const u8) ![]u8 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    // Not on Windows, where the stat size can differ from the bytes read.
    if (!is_windows) sized: {
        const stat = file.stat(io) catch break :sized;
        const size = std.math.cast(usize, stat.size) orelse return error.FileTooBig;
        // One byte more than the size, so that a file that has grown is noticed.
        const buf = try allocator.alloc(u8, size + 1);
        const n = file.readPositionalAll(io, buf, 0) catch {
            allocator.free(buf);
            break :sized;
        };
        if (n > size) {
            allocator.free(buf);
            break :sized;
        }
        if (allocator.resize(buf, n)) return buf[0..n];
        defer allocator.free(buf);
        return allocator.dupe(u8, buf[0..n]);
    }

    var read_buf: [64 * 1024]u8 = undefined;
    var reader = file.reader(io, &read_buf);
    return reader.interface.allocRemaining(allocator, .unlimited);
}

test "mmapFile reads valid PDB" {
    const mapped = try mmapFile(std.testing.allocator, std.testing.io, "test_data/1l2y.pdb");
    defer mapped.deinit();
    try std.testing.expect(mapped.data.len > 0);
    try std.testing.expect(std.mem.startsWith(u8, mapped.data, "REMARK") or
        std.mem.startsWith(u8, mapped.data, "HEADER") or
        std.mem.startsWith(u8, mapped.data, "ATOM"));
}

test "readFileOnce returns the same bytes as mmapFile" {
    const allocator = std.testing.allocator;
    const mapped = try mmapFile(allocator, std.testing.io, "test_data/1l2y.pdb");
    defer mapped.deinit();
    const data = try readFileOnce(allocator, std.testing.io, "test_data/1l2y.pdb");
    defer allocator.free(data);
    try std.testing.expectEqualSlices(u8, mapped.data, data);
}

test "readFileOnce reads an empty file" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "empty.pdb", .data = "" });
    const path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "empty.pdb", allocator);
    defer allocator.free(path);
    const data = try readFileOnce(allocator, std.testing.io, path);
    defer allocator.free(data);
    try std.testing.expectEqual(@as(usize, 0), data.len);
}

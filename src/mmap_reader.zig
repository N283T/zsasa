const std = @import("std");
const builtin = @import("builtin");

pub const is_windows = builtin.os.tag == .windows;

/// A read-only file region backed by mmap (POSIX) or heap allocation (Windows).
///
/// Ownership: single-owner, non-copyable by convention. Exactly one `deinit`
/// call must occur per `mmapFile` call.
/// Typical usage: `const mapped = try mmapFile(alloc, path); defer mapped.deinit(alloc);`
pub const MappedFile = struct {
    data: if (is_windows) []const u8 else []align(std.heap.page_size_min) const u8,

    pub fn deinit(self: MappedFile, allocator: if (is_windows) std.mem.Allocator else void) void {
        if (self.data.len == 0) return;
        if (is_windows) {
            allocator.free(self.data);
        } else {
            std.posix.munmap(@alignCast(self.data));
        }
    }
};

pub fn mmapFile(allocator: if (is_windows) std.mem.Allocator else void, path: []const u8) !MappedFile {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();

    if (is_windows) {
        const data = try file.readToEndAlloc(allocator, std.math.maxInt(u32));
        return .{ .data = data };
    } else {
        const stat = try file.stat();
        const size: usize = std.math.cast(usize, stat.size) orelse return error.FileTooBig;
        if (size == 0) return .{ .data = &.{} };

        const mapped = try std.posix.mmap(
            null,
            size,
            std.posix.PROT.READ,
            .{ .TYPE = .PRIVATE },
            file.handle,
            0,
        );
        return .{ .data = mapped };
    }
}

test "mmapFile reads valid PDB" {
    const alloc_arg = if (is_windows) std.testing.allocator else {};
    const mapped = try mmapFile(alloc_arg, "test_data/1l2y.pdb");
    defer mapped.deinit(alloc_arg);
    try std.testing.expect(mapped.data.len > 0);
    try std.testing.expect(std.mem.startsWith(u8, mapped.data, "REMARK") or
        std.mem.startsWith(u8, mapped.data, "HEADER") or
        std.mem.startsWith(u8, mapped.data, "ATOM"));
}

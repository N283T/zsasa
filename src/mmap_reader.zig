const std = @import("std");

/// A read-only memory-mapped file region.
///
/// Ownership: single-owner, non-copyable by convention. Exactly one `deinit`
/// call must occur per `mmapFile` call.
/// Typical usage: `const mapped = try mmapFile(path); defer mapped.deinit();`
pub const MappedFile = struct {
    data: []align(std.heap.page_size_min) const u8,

    pub fn deinit(self: MappedFile) void {
        if (self.data.len == 0) return;
        std.posix.munmap(self.data);
    }
};

pub fn mmapFile(path: []const u8) !MappedFile {
    const file = try std.fs.cwd().openFile(path, .{});
    defer file.close();

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

test "mmapFile reads valid PDB" {
    const mapped = try mmapFile("test_data/1l2y.pdb");
    defer mapped.deinit();
    try std.testing.expect(mapped.data.len > 0);
    try std.testing.expect(std.mem.startsWith(u8, mapped.data, "REMARK") or
        std.mem.startsWith(u8, mapped.data, "HEADER") or
        std.mem.startsWith(u8, mapped.data, "ATOM"));
}

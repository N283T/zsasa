const std = @import("std");

pub const InputIoMode = enum {
    /// Let the caller choose its established default.
    auto,
    /// Parse directly from a memory-mapped file.
    mmap,
    /// Read the complete file into an allocated buffer before parsing.
    read,

    pub fn parse(value: []const u8) ?InputIoMode {
        if (std.mem.eql(u8, value, "auto")) return .auto;
        if (std.mem.eql(u8, value, "mmap")) return .mmap;
        if (std.mem.eql(u8, value, "read")) return .read;
        return null;
    }

    pub fn resolve(self: InputIoMode, auto_mode: InputIoMode) InputIoMode {
        std.debug.assert(auto_mode != .auto);
        return if (self == .auto) auto_mode else self;
    }
};

test "InputIoMode parses CLI values" {
    try std.testing.expectEqual(InputIoMode.auto, InputIoMode.parse("auto").?);
    try std.testing.expectEqual(InputIoMode.mmap, InputIoMode.parse("mmap").?);
    try std.testing.expectEqual(InputIoMode.read, InputIoMode.parse("read").?);
    try std.testing.expect(InputIoMode.parse("invalid") == null);
}

test "InputIoMode resolves auto without overriding explicit modes" {
    try std.testing.expectEqual(InputIoMode.read, InputIoMode.auto.resolve(.read));
    try std.testing.expectEqual(InputIoMode.mmap, InputIoMode.mmap.resolve(.read));
}

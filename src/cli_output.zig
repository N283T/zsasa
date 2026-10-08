//! Standard output of the CLI for text the user asked to read (`--help`,
//! `--version`), so that it can be piped and captured. Progress, summaries,
//! warnings and error messages stay on standard error (`std.debug.print`).

const std = @import("std");
const builtin = @import("builtin");

/// Formats to standard output.
///
/// A write error (a closed pipe, a full disk) is ignored: the caller is about
/// to exit and has nothing better to do with help text it cannot write.
/// Prints nothing in a test binary, where standard output carries the
/// protocol between the test runner and `zig build`.
pub fn print(io: std.Io, comptime fmt: []const u8, args: anytype) void {
    if (builtin.is_test) return;
    var buffer: [4096]u8 = undefined;
    var writer = std.Io.File.Writer.initStreaming(std.Io.File.stdout(), io, &buffer);
    writer.interface.print(fmt, args) catch return;
    writer.interface.flush() catch return;
}

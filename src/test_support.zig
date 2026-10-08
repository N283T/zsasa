//! Helpers shared by the Zig tests. Imported from `test` blocks only, so none
//! of this reaches a release build.

const std = @import("std");
const builtin = @import("builtin");

/// Environment variable that keeps `muteStderr` from muting anything, for
/// reading the diagnostics a test prints while debugging it.
pub const show_stderr_env = "ZSASA_TEST_STDERR";

/// Stderr of the process while it is muted; see `muteStderr`.
pub const MutedStderr = struct {
    /// A duplicate of the original stderr descriptor, or null when nothing was muted.
    saved: ?i32 = null,

    pub fn restore(self: *MutedStderr) void {
        const saved = self.saved orelse return;
        self.saved = null;
        _ = posixDup2(saved, 2);
        posixClose(saved);
    }
};

/// Discards what the code under test prints with `std.debug.print` (the
/// "Error: ..." and "Workflow complete: ..." lines of the CLI paths), so a
/// passing `zig build test` prints nothing from the tests. Those messages
/// have no quiet flag or writer to route them through, because they are the
/// output of the commands the tests run. Use it as
///
///     var muted = test_support.muteStderr();
///     defer muted.restore();
///
/// Everything is printed again when ZSASA_TEST_STDERR is set, and on Windows,
/// where nothing is muted.
pub fn muteStderr() MutedStderr {
    if (builtin.os.tag == .windows or builtin.os.tag == .wasi) return .{};
    if (std.testing.environ.getPosix(show_stderr_env) != null) return .{};

    const saved = posixDup(2) orelse return .{};
    const null_file = std.Io.Dir.openFileAbsolute(std.testing.io, "/dev/null", .{ .mode = .write_only }) catch {
        posixClose(saved);
        return .{};
    };
    defer null_file.close(std.testing.io);
    if (!posixDup2(null_file.handle, 2)) {
        posixClose(saved);
        return .{};
    }
    return .{ .saved = saved };
}

fn posixDup(fd: i32) ?i32 {
    const rc = std.posix.system.dup(fd);
    if (@TypeOf(rc) == usize) {
        if (std.os.linux.errno(rc) != .SUCCESS) return null;
        return @intCast(rc);
    }
    return if (rc < 0) null else @intCast(rc);
}

fn posixDup2(old: i32, new: i32) bool {
    const rc = std.posix.system.dup2(old, new);
    if (@TypeOf(rc) == usize) return std.os.linux.errno(rc) == .SUCCESS;
    return rc >= 0;
}

fn posixClose(fd: i32) void {
    _ = std.posix.system.close(fd);
}

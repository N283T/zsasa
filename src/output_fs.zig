//! Creating output directories.

const std = @import("std");

pub const CreateDirPathError = std.Io.Dir.CreateDirError || std.Io.Dir.StatFileError;

/// Create `path` and any missing parent directories, relative to the current
/// working directory.
///
/// Unlike `std.Io.Dir.createDirPath`, an existing component that is a symbolic
/// link to a directory is accepted (macOS `/tmp` links to `/private/tmp`). A
/// component that exists but is not a directory, or a symlink to one that is
/// not, fails with `error.NotDir`; a dangling symlink fails with
/// `error.FileNotFound`.
pub fn createDirPath(io: std.Io, path: []const u8) CreateDirPathError!void {
    const cwd = std.Io.Dir.cwd();
    if (path.len == 0) return;

    // Fast path, and the only check needed when the directory already exists.
    if (cwd.statFile(io, path, .{})) |stat| {
        if (stat.kind != .directory) return error.NotDir;
        return;
    } else |err| switch (err) {
        error.FileNotFound => {},
        else => |e| return e,
    }

    // A missing path whose parent is missing as well: create the parent first.
    if (std.fs.path.dirname(path)) |parent| {
        try createDirPath(io, parent);
    }

    cwd.createDir(io, path, .default_dir) catch |err| switch (err) {
        // Created concurrently, or a dangling symlink: check what is there now.
        error.PathAlreadyExists => {
            const stat = try cwd.statFile(io, path, .{});
            if (stat.kind != .directory) return error.NotDir;
        },
        else => |e| return e,
    };
}

fn realTmpPath(tmp_dir: std.testing.TmpDir, buf: []u8) ![]const u8 {
    const len = try tmp_dir.dir.realPath(std.testing.io, buf);
    return buf[0..len];
}

test "createDirPath creates nested directories" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root = try realTmpPath(tmp_dir, &root_buf);

    const path = try std.fs.path.join(allocator, &.{ root, "a", "b", "c" });
    defer allocator.free(path);
    try createDirPath(std.testing.io, path);
    const stat = try tmp_dir.dir.statFile(std.testing.io, "a/b/c", .{});
    try std.testing.expectEqual(std.Io.File.Kind.directory, stat.kind);

    // Calling it again on the existing path succeeds.
    try createDirPath(std.testing.io, path);
}

test "createDirPath accepts a symlink to a directory and creates under it" {
    if (@import("builtin").os.tag == .windows) return error.SkipZigTest;
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root = try realTmpPath(tmp_dir, &root_buf);

    try tmp_dir.dir.createDir(std.testing.io, "real", .default_dir);
    try tmp_dir.dir.symLink(std.testing.io, "real", "link", .{ .is_directory = true });

    // The symlink itself, as for `calc ... /tmp/out.json` on macOS.
    const link_path = try std.fs.path.join(allocator, &.{ root, "link" });
    defer allocator.free(link_path);
    try createDirPath(std.testing.io, link_path);

    // A nested directory below the symlink.
    const nested = try std.fs.path.join(allocator, &.{ root, "link", "x", "y" });
    defer allocator.free(nested);
    try createDirPath(std.testing.io, nested);
    const stat = try tmp_dir.dir.statFile(std.testing.io, "real/x/y", .{});
    try std.testing.expectEqual(std.Io.File.Kind.directory, stat.kind);
}

test "createDirPath rejects a file and a symlink to a file" {
    if (@import("builtin").os.tag == .windows) return error.SkipZigTest;
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root = try realTmpPath(tmp_dir, &root_buf);

    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "file", .data = "x" });
    try tmp_dir.dir.symLink(std.testing.io, "file", "file_link", .{});

    const file_path = try std.fs.path.join(allocator, &.{ root, "file" });
    defer allocator.free(file_path);
    try std.testing.expectError(error.NotDir, createDirPath(std.testing.io, file_path));

    const link_path = try std.fs.path.join(allocator, &.{ root, "file_link" });
    defer allocator.free(link_path);
    try std.testing.expectError(error.NotDir, createDirPath(std.testing.io, link_path));

    // A directory below a regular file cannot be created either.
    const below = try std.fs.path.join(allocator, &.{ root, "file", "sub" });
    defer allocator.free(below);
    try std.testing.expect(std.meta.isError(createDirPath(std.testing.io, below)));
}

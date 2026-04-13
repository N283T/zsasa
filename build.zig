const std = @import("std");

const version = "0.2.8";

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // zlib dependency (built from source via allyourcodebase/zlib)
    const zlib_dep = b.dependency("zlib", .{
        .target = target,
        .optimize = optimize,
    });
    const zlib_artifact = zlib_dep.artifact("z");

    // PIC-enabled zlib for shared library builds
    const zlib_pic_dep = b.dependency("zlib", .{
        .target = target,
        .optimize = optimize,
        .pie = true,
    });
    const zlib_pic_artifact = zlib_pic_dep.artifact("z");

    // Library module (exposed to package consumers via zig fetch)
    const mod = b.addModule("zsasa", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });
    mod.linkLibrary(zlib_artifact);

    const zxdrfile_dep = b.dependency("zxdrfile", .{
        .target = target,
        .optimize = optimize,
    });
    const zxdrfile_mod = zxdrfile_dep.module("zxdrfile");

    // CLI executable
    const options = b.addOptions();
    options.addOption([]const u8, "version", version);

    const exe_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zsasa", .module = mod },
            .{ .name = "build_options", .module = options.createModule() },
            .{ .name = "zxdrfile", .module = zxdrfile_mod },
        },
    });
    exe_module.linkLibrary(zlib_artifact);

    const exe = b.addExecutable(.{
        .name = "zsasa",
        .root_module = exe_module,
    });
    b.installArtifact(exe);

    // Shared library for C API / Python bindings
    const lib_module = b.createModule(.{
        .root_source_file = b.path("src/c_api.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zxdrfile", .module = zxdrfile_mod },
        },
    });
    lib_module.linkLibrary(zlib_pic_artifact);

    const lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "zsasa",
        .root_module = lib_module,
    });
    lib.linkLibC();
    b.installArtifact(lib);

    // Run step
    const run_step = b.step("run", "Run the app");
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // Test step
    const mod_tests = b.addTest(.{ .root_module = mod });
    const exe_tests = b.addTest(.{ .root_module = exe.root_module });
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&b.addRunArtifact(mod_tests).step);
    test_step.dependOn(&b.addRunArtifact(exe_tests).step);

    // Docs step (zig autodoc)
    const docs_lib = b.addLibrary(.{
        .linkage = .static,
        .name = "zsasa",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/root.zig"),
            .target = target,
        }),
    });
    const install_docs = b.addInstallDirectory(.{
        .source_dir = docs_lib.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });
    const docs_step = b.step("docs", "Emit autodoc to zig-out/docs");
    docs_step.dependOn(&install_docs.step);
}

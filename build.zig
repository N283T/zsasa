const std = @import("std");

const version = "0.1.1";

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Library module (exposed to package consumers via zig fetch)
    const mod = b.addModule("zsasa", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });

    // CLI executable
    const options = b.addOptions();
    options.addOption([]const u8, "version", version);

    const exe = b.addExecutable(.{
        .name = "zsasa",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zsasa", .module = mod },
                .{ .name = "build_options", .module = options.createModule() },
            },
        }),
    });
    b.installArtifact(exe);

    // Shared library for C API / Python bindings
    const lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "zsasa",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/c_api.zig"),
            .target = target,
            .optimize = optimize,
        }),
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
}

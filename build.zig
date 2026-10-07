const std = @import("std");

const version = "0.9.1";

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Library module (exposed to package consumers via zig fetch)
    const mod = b.addModule("zsasa", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });

    const ztraj_dep = b.dependency("ztraj", .{
        .target = target,
        .optimize = optimize,
    });
    const ztraj_mod = ztraj_dep.module("ztraj");

    // CLI executable
    const options = b.addOptions();
    options.addOption([]const u8, "version", version);

    const options_mod = options.createModule();

    const exe_module = b.createModule(.{
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{
            .{ .name = "zsasa", .module = mod },
            .{ .name = "build_options", .module = options_mod },
            .{ .name = "ztraj", .module = ztraj_mod },
        },
    });

    const exe = b.addExecutable(.{
        .name = "zsasa",
        .root_module = exe_module,
    });
    b.installArtifact(exe);

    // Shared library for C API / Python bindings.
    // libc is required because c_api.zig uses std.heap.c_allocator (so the
    // FFI surface uses C's malloc/free rather than a per-call GeneralPurposeAllocator,
    // matching Python ctypes' lifetime expectations). On macOS the dylib auto-links
    // libSystem so this is implicit, but Linux needs the explicit flag.
    const lib_module = b.createModule(.{
        .root_source_file = b.path("src/c_api.zig"),
        .target = target,
        .optimize = optimize,
        .link_libc = true,
        .imports = &.{
            .{ .name = "ztraj", .module = ztraj_mod },
        },
    });

    const lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "zsasa",
        .root_module = lib_module,
    });
    b.installArtifact(lib);

    // Run step
    const run_step = b.step("run", "Run the app");
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    // Test step.
    //
    // Zig runs the tests of every file reachable from a test root, and the
    // three roots overlap heavily (c_api.zig alone reaches almost every file).
    // Each test should run exactly once, so the library root stays unfiltered
    // and the other two roots only run the tests no earlier root reaches:
    //   - the zsasa module (src/root.zig): dcd.zig and root.zig itself
    //   - the executable (src/main.zig): calc.zig, traj.zig and main.zig itself
    // Filters are substring matches on the full test name.
    // scripts/check_test_partition.py verifies that every test runs in exactly
    // one artifact, so a new file reachable from only one root cannot silently
    // lose its tests or run them twice.
    const all_tests = b.option(
        bool,
        "all-tests",
        "Disable the per-artifact test filters (used by scripts/check_test_partition.py)",
    ) orelse false;
    const mod_tests = b.addTest(.{
        .root_module = mod,
        .filters = if (all_tests) &.{} else &.{ "dcd.test.", "root.test" },
    });
    const exe_tests = b.addTest(.{
        .root_module = exe.root_module,
        .filters = if (all_tests) &.{} else &.{ "calc.test.", "traj.test.", "main.test" },
    });
    const lib_tests = b.addTest(.{ .root_module = lib.root_module });
    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&b.addRunArtifact(mod_tests).step);
    test_step.dependOn(&b.addRunArtifact(exe_tests).step);
    test_step.dependOn(&b.addRunArtifact(lib_tests).step);

    // Install the test executables without running them, so
    // scripts/check_test_partition.py can list the tests each one contains.
    const test_bins_step = b.step("test-bins", "Install the test executables to <prefix>/test-bin");
    const test_bins = [_]struct { name: []const u8, artifact: *std.Build.Step.Compile }{
        .{ .name = "mod-tests", .artifact = mod_tests },
        .{ .name = "exe-tests", .artifact = exe_tests },
        .{ .name = "lib-tests", .artifact = lib_tests },
    };
    for (test_bins) |test_bin| {
        const install = b.addInstallArtifact(test_bin.artifact, .{
            .dest_dir = .{ .override = .{ .custom = "test-bin" } },
            .dest_sub_path = test_bin.name,
        });
        test_bins_step.dependOn(&install.step);
    }

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

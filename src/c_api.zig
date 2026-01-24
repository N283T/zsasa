//! C ABI interface for freesasa-zig library.
//!
//! This module provides C-compatible functions that can be called from
//! other languages (Python, C, etc.) via FFI/ctypes.

const std = @import("std");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");

const AtomInput = types.AtomInput;
const Config = types.Config;

/// No error - calculation completed successfully
pub const FREESASA_OK: c_int = 0;
/// Invalid input parameters (n_atoms=0, n_points=0, n_slices=0, or invalid probe_radius)
/// Note: Passing NULL pointers results in undefined behavior
pub const FREESASA_ERROR_INVALID_INPUT: c_int = -1;
/// Memory allocation failed during calculation
pub const FREESASA_ERROR_OUT_OF_MEMORY: c_int = -2;
/// Internal calculation error
pub const FREESASA_ERROR_CALCULATION: c_int = -3;

// Version string
const VERSION = "0.0.5";

/// Thread-safe allocator for C API (uses C allocator for simplicity)
const c_allocator = std.heap.c_allocator;

/// Get library version string.
export fn freesasa_version() callconv(.c) [*:0]const u8 {
    return VERSION;
}

/// Calculate SASA using Shrake-Rupley algorithm.
///
/// Parameters:
///   x, y, z: Atom coordinates (arrays of n_atoms elements)
///   radii: Atom radii (array of n_atoms elements)
///   n_atoms: Number of atoms
///   n_points: Number of test points per atom (e.g., 100)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (must be pre-allocated, n_atoms elements)
///   total_area: Output pointer for total SASA
///
/// Returns:
///   FREESASA_OK (0) on success, negative error code on failure.
export fn freesasa_calc_sr(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
    n_points: u32,
    probe_radius: f64,
    n_threads: usize,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.c) c_int {
    // Validate input
    if (n_atoms == 0 or n_points == 0 or probe_radius <= 0.0) {
        return FREESASA_ERROR_INVALID_INPUT;
    }

    // Create AtomInput from raw arrays
    const input = AtomInput{
        .x = x[0..n_atoms],
        .y = y[0..n_atoms],
        .z = z[0..n_atoms],
        .r = radii[0..n_atoms],
        .allocator = c_allocator,
    };

    const config = Config{
        .n_points = n_points,
        .probe_radius = probe_radius,
    };

    // Calculate SASA
    const result = if (n_threads == 1)
        shrake_rupley.calculateSasa(c_allocator, input, config) catch {
            return FREESASA_ERROR_CALCULATION;
        }
    else
        shrake_rupley.calculateSasaParallel(c_allocator, input, config, n_threads) catch {
            return FREESASA_ERROR_CALCULATION;
        };
    defer {
        // Free the result's internal allocation
        c_allocator.free(result.atom_areas);
    }

    // Copy results to output buffers
    @memcpy(atom_areas[0..n_atoms], result.atom_areas);
    total_area.* = result.total_area;

    return FREESASA_OK;
}

/// Calculate SASA using Lee-Richards algorithm.
///
/// Parameters:
///   x, y, z: Atom coordinates (arrays of n_atoms elements)
///   radii: Atom radii (array of n_atoms elements)
///   n_atoms: Number of atoms
///   n_slices: Number of slices per atom (e.g., 20)
///   probe_radius: Water probe radius in Angstroms (e.g., 1.4)
///   n_threads: Number of threads (0 = auto-detect)
///   atom_areas: Output buffer for per-atom SASA (must be pre-allocated, n_atoms elements)
///   total_area: Output pointer for total SASA
///
/// Returns:
///   FREESASA_OK (0) on success, negative error code on failure.
export fn freesasa_calc_lr(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    radii: [*]const f64,
    n_atoms: usize,
    n_slices: u32,
    probe_radius: f64,
    n_threads: usize,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.c) c_int {
    // Validate input
    if (n_atoms == 0 or n_slices == 0 or probe_radius <= 0.0) {
        return FREESASA_ERROR_INVALID_INPUT;
    }

    // Create AtomInput from raw arrays
    const input = AtomInput{
        .x = x[0..n_atoms],
        .y = y[0..n_atoms],
        .z = z[0..n_atoms],
        .r = radii[0..n_atoms],
        .allocator = c_allocator,
    };

    const config = lee_richards.LeeRichardsConfig{
        .n_slices = n_slices,
        .probe_radius = probe_radius,
    };

    // Calculate SASA
    const result = if (n_threads == 1)
        lee_richards.calculateSasa(c_allocator, input, config) catch {
            return FREESASA_ERROR_CALCULATION;
        }
    else
        lee_richards.calculateSasaParallel(c_allocator, input, config, n_threads) catch {
            return FREESASA_ERROR_CALCULATION;
        };
    defer {
        // Free the result's internal allocation
        c_allocator.free(result.atom_areas);
    }

    // Copy results to output buffers
    @memcpy(atom_areas[0..n_atoms], result.atom_areas);
    total_area.* = result.total_area;

    return FREESASA_OK;
}

// Tests
test "freesasa_version returns valid string" {
    const version = freesasa_version();
    try std.testing.expect(version[0] != 0);
}

test "freesasa_calc_sr with empty input returns error" {
    var total_area: f64 = 0.0;
    const result = freesasa_calc_sr(
        undefined,
        undefined,
        undefined,
        undefined,
        0, // n_atoms = 0
        100,
        1.4,
        1,
        undefined,
        &total_area,
    );
    try std.testing.expectEqual(FREESASA_ERROR_INVALID_INPUT, result);
}

test "freesasa_calc_sr single atom" {
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    const result = freesasa_calc_sr(
        &x,
        &y,
        &z,
        &radii,
        1,
        100,
        1.4,
        1,
        &atom_areas,
        &total_area,
    );

    try std.testing.expectEqual(FREESASA_OK, result);
    // Expected: 4π * (1.5 + 1.4)² ≈ 105.68 Ų
    try std.testing.expect(total_area > 100.0 and total_area < 110.0);
    try std.testing.expectEqual(total_area, atom_areas[0]);
}

test "freesasa_calc_lr single atom" {
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const radii = [_]f64{1.5};
    var atom_areas = [_]f64{0.0};
    var total_area: f64 = 0.0;

    const result = freesasa_calc_lr(
        &x,
        &y,
        &z,
        &radii,
        1,
        20,
        1.4,
        1,
        &atom_areas,
        &total_area,
    );

    try std.testing.expectEqual(FREESASA_OK, result);
    // Expected: 4π * (1.5 + 1.4)² ≈ 105.68 Ų
    try std.testing.expect(total_area > 100.0 and total_area < 110.0);
    try std.testing.expectEqual(total_area, atom_areas[0]);
}

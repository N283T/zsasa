//! C ABI interface for freesasa-zig library.
//!
//! This module provides C-compatible functions that can be called from
//! other languages (Python, C, etc.) via FFI/ctypes.

const std = @import("std");
const shrake_rupley = @import("shrake_rupley.zig");
const lee_richards = @import("lee_richards.zig");
const types = @import("types.zig");
const classifier = @import("classifier.zig");
const classifier_naccess = @import("classifier_naccess.zig");
const classifier_protor = @import("classifier_protor.zig");
const classifier_oons = @import("classifier_oons.zig");

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

// =============================================================================
// Classifier Types
// =============================================================================

/// NACCESS-compatible radii (default)
pub const FREESASA_CLASSIFIER_NACCESS: c_int = 0;
/// ProtOr radii based on hybridization state
pub const FREESASA_CLASSIFIER_PROTOR: c_int = 1;
/// OONS radii (older FreeSASA default)
pub const FREESASA_CLASSIFIER_OONS: c_int = 2;

// =============================================================================
// Atom Classes
// =============================================================================

/// Polar atom class
pub const FREESASA_ATOM_CLASS_POLAR: c_int = 0;
/// Apolar atom class
pub const FREESASA_ATOM_CLASS_APOLAR: c_int = 1;
/// Unknown atom class
pub const FREESASA_ATOM_CLASS_UNKNOWN: c_int = 2;

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

// =============================================================================
// Classifier Functions
// =============================================================================

/// Get van der Waals radius for an atom using the specified classifier.
///
/// Parameters:
///   classifier_type: Classifier to use (FREESASA_CLASSIFIER_NACCESS, etc.)
///   residue: Residue name (e.g., "ALA", "GLY") - null-terminated
///   atom: Atom name (e.g., "CA", "CB") - null-terminated
///
/// Returns:
///   Radius in Angstroms, or NaN if atom is not found in classifier.
export fn freesasa_classifier_get_radius(
    classifier_type: c_int,
    residue: [*:0]const u8,
    atom: [*:0]const u8,
) callconv(.c) f64 {
    const residue_slice = std.mem.span(residue);
    const atom_slice = std.mem.span(atom);

    const radius: ?f64 = switch (classifier_type) {
        FREESASA_CLASSIFIER_NACCESS => classifier_naccess.getRadius(residue_slice, atom_slice),
        FREESASA_CLASSIFIER_PROTOR => classifier_protor.getRadius(residue_slice, atom_slice),
        FREESASA_CLASSIFIER_OONS => classifier_oons.getRadius(residue_slice, atom_slice),
        else => null,
    };

    return radius orelse std.math.nan(f64);
}

/// Get atom polarity class using the specified classifier.
///
/// Parameters:
///   classifier_type: Classifier to use (FREESASA_CLASSIFIER_NACCESS, etc.)
///   residue: Residue name (e.g., "ALA", "GLY") - null-terminated
///   atom: Atom name (e.g., "CA", "CB") - null-terminated
///
/// Returns:
///   Atom class (FREESASA_ATOM_CLASS_POLAR, FREESASA_ATOM_CLASS_APOLAR, or FREESASA_ATOM_CLASS_UNKNOWN)
export fn freesasa_classifier_get_class(
    classifier_type: c_int,
    residue: [*:0]const u8,
    atom: [*:0]const u8,
) callconv(.c) c_int {
    const residue_slice = std.mem.span(residue);
    const atom_slice = std.mem.span(atom);

    const atom_class = switch (classifier_type) {
        FREESASA_CLASSIFIER_NACCESS => classifier_naccess.getClass(residue_slice, atom_slice),
        FREESASA_CLASSIFIER_PROTOR => classifier_protor.getClass(residue_slice, atom_slice),
        FREESASA_CLASSIFIER_OONS => classifier_oons.getClass(residue_slice, atom_slice),
        else => classifier.AtomClass.unknown,
    };

    return switch (atom_class) {
        .polar => FREESASA_ATOM_CLASS_POLAR,
        .apolar => FREESASA_ATOM_CLASS_APOLAR,
        .unknown => FREESASA_ATOM_CLASS_UNKNOWN,
    };
}

/// Guess van der Waals radius from element symbol.
///
/// Parameters:
///   element: Element symbol (e.g., "C", "N", "FE") - null-terminated
///            Case-insensitive, whitespace is trimmed.
///
/// Returns:
///   Radius in Angstroms, or NaN if element is not recognized.
export fn freesasa_guess_radius(
    element: [*:0]const u8,
) callconv(.c) f64 {
    const element_slice = std.mem.span(element);
    return classifier.guessRadius(element_slice) orelse std.math.nan(f64);
}

/// Guess van der Waals radius from PDB-style atom name.
/// Extracts element symbol from atom name and returns corresponding radius.
///
/// Parameters:
///   atom_name: PDB-style atom name (e.g., " CA ", "FE  ") - null-terminated
///              Following PDB conventions:
///              - Leading space indicates single-char element (e.g., " CA " = Carbon alpha)
///              - No leading space may indicate 2-char element (e.g., "FE  " = Iron)
///
/// Returns:
///   Radius in Angstroms, or NaN if element cannot be determined.
export fn freesasa_guess_radius_from_atom_name(
    atom_name: [*:0]const u8,
) callconv(.c) f64 {
    const atom_slice = std.mem.span(atom_name);
    return classifier.guessRadiusFromAtomName(atom_slice) orelse std.math.nan(f64);
}

/// Classify multiple atoms at once (batch operation).
///
/// This is more efficient than calling freesasa_classifier_get_radius
/// for each atom individually.
///
/// Parameters:
///   classifier_type: Classifier to use (FREESASA_CLASSIFIER_NACCESS, etc.)
///   residues: Array of residue name pointers (null-terminated strings)
///   atoms: Array of atom name pointers (null-terminated strings)
///   n_atoms: Number of atoms
///   radii_out: Output buffer for radii (pre-allocated, n_atoms elements)
///              NaN is written for atoms not found in classifier
///   classes_out: Output buffer for classes (pre-allocated, n_atoms elements)
///                Can be NULL if classes are not needed
///
/// Returns:
///   FREESASA_OK on success, negative error code on failure.
export fn freesasa_classify_atoms(
    classifier_type: c_int,
    residues: [*]const [*:0]const u8,
    atoms: [*]const [*:0]const u8,
    n_atoms: usize,
    radii_out: [*]f64,
    classes_out: ?[*]c_int,
) callconv(.c) c_int {
    if (n_atoms == 0) {
        return FREESASA_OK; // Nothing to do
    }

    // Validate classifier type
    if (classifier_type < FREESASA_CLASSIFIER_NACCESS or classifier_type > FREESASA_CLASSIFIER_OONS) {
        return FREESASA_ERROR_INVALID_INPUT;
    }

    for (0..n_atoms) |i| {
        const residue_slice = std.mem.span(residues[i]);
        const atom_slice = std.mem.span(atoms[i]);

        // Get radius
        const radius: ?f64 = switch (classifier_type) {
            FREESASA_CLASSIFIER_NACCESS => classifier_naccess.getRadius(residue_slice, atom_slice),
            FREESASA_CLASSIFIER_PROTOR => classifier_protor.getRadius(residue_slice, atom_slice),
            FREESASA_CLASSIFIER_OONS => classifier_oons.getRadius(residue_slice, atom_slice),
            else => null,
        };
        radii_out[i] = radius orelse std.math.nan(f64);

        // Get class if requested
        if (classes_out) |classes| {
            const atom_class = switch (classifier_type) {
                FREESASA_CLASSIFIER_NACCESS => classifier_naccess.getClass(residue_slice, atom_slice),
                FREESASA_CLASSIFIER_PROTOR => classifier_protor.getClass(residue_slice, atom_slice),
                FREESASA_CLASSIFIER_OONS => classifier_oons.getClass(residue_slice, atom_slice),
                else => classifier.AtomClass.unknown,
            };
            classes[i] = switch (atom_class) {
                .polar => FREESASA_ATOM_CLASS_POLAR,
                .apolar => FREESASA_ATOM_CLASS_APOLAR,
                .unknown => FREESASA_ATOM_CLASS_UNKNOWN,
            };
        }
    }

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

// =============================================================================
// Classifier Tests
// =============================================================================

test "freesasa_classifier_get_radius NACCESS" {
    // Standard backbone atoms
    const ca_radius = freesasa_classifier_get_radius(FREESASA_CLASSIFIER_NACCESS, "ALA", "CA");
    try std.testing.expectApproxEqAbs(1.87, ca_radius, 0.01);

    const n_radius = freesasa_classifier_get_radius(FREESASA_CLASSIFIER_NACCESS, "ALA", "N");
    try std.testing.expectApproxEqAbs(1.65, n_radius, 0.01);

    const o_radius = freesasa_classifier_get_radius(FREESASA_CLASSIFIER_NACCESS, "ALA", "O");
    try std.testing.expectApproxEqAbs(1.40, o_radius, 0.01);

    // Unknown atom should return NaN
    const unknown = freesasa_classifier_get_radius(FREESASA_CLASSIFIER_NACCESS, "ALA", "XX");
    try std.testing.expect(std.math.isNan(unknown));
}

test "freesasa_classifier_get_radius ProtoR" {
    const ca_radius = freesasa_classifier_get_radius(FREESASA_CLASSIFIER_PROTOR, "ALA", "CA");
    try std.testing.expect(ca_radius > 1.0 and ca_radius < 3.0);
}

test "freesasa_classifier_get_radius OONS" {
    const ca_radius = freesasa_classifier_get_radius(FREESASA_CLASSIFIER_OONS, "ALA", "CA");
    try std.testing.expect(ca_radius > 1.0 and ca_radius < 3.0);
}

test "freesasa_classifier_get_radius invalid classifier" {
    const radius = freesasa_classifier_get_radius(99, "ALA", "CA");
    try std.testing.expect(std.math.isNan(radius));
}

test "freesasa_classifier_get_class" {
    // Carbon atoms are apolar
    const ca_class = freesasa_classifier_get_class(FREESASA_CLASSIFIER_NACCESS, "ALA", "CA");
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_APOLAR, ca_class);

    // Oxygen atoms are polar
    const o_class = freesasa_classifier_get_class(FREESASA_CLASSIFIER_NACCESS, "ALA", "O");
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_POLAR, o_class);

    // Nitrogen atoms are polar
    const n_class = freesasa_classifier_get_class(FREESASA_CLASSIFIER_NACCESS, "ALA", "N");
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_POLAR, n_class);

    // Unknown atoms
    const unknown_class = freesasa_classifier_get_class(FREESASA_CLASSIFIER_NACCESS, "ALA", "XX");
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_UNKNOWN, unknown_class);
}

test "freesasa_guess_radius" {
    // Common elements
    try std.testing.expectApproxEqAbs(1.70, freesasa_guess_radius("C"), 0.01);
    try std.testing.expectApproxEqAbs(1.55, freesasa_guess_radius("N"), 0.01);
    try std.testing.expectApproxEqAbs(1.52, freesasa_guess_radius("O"), 0.01);
    try std.testing.expectApproxEqAbs(1.80, freesasa_guess_radius("S"), 0.01);

    // Case insensitive
    try std.testing.expectApproxEqAbs(1.70, freesasa_guess_radius("c"), 0.01);

    // Two-character elements
    try std.testing.expectApproxEqAbs(1.26, freesasa_guess_radius("FE"), 0.01);
    try std.testing.expectApproxEqAbs(1.39, freesasa_guess_radius("ZN"), 0.01);

    // Unknown element returns NaN
    const unknown = freesasa_guess_radius("XX");
    try std.testing.expect(std.math.isNan(unknown));
}

test "freesasa_guess_radius_from_atom_name" {
    // Standard PDB atom names (leading space = single-char element)
    try std.testing.expectApproxEqAbs(1.70, freesasa_guess_radius_from_atom_name(" CA "), 0.01);
    try std.testing.expectApproxEqAbs(1.55, freesasa_guess_radius_from_atom_name(" N  "), 0.01);
    try std.testing.expectApproxEqAbs(1.52, freesasa_guess_radius_from_atom_name(" O  "), 0.01);

    // Metal atoms (no leading space = 2-char element)
    try std.testing.expectApproxEqAbs(1.26, freesasa_guess_radius_from_atom_name("FE  "), 0.01);
    try std.testing.expectApproxEqAbs(1.39, freesasa_guess_radius_from_atom_name("ZN  "), 0.01);
}

test "freesasa_classify_atoms batch" {
    const residues = [_][*:0]const u8{ "ALA", "ALA", "GLY" };
    const atoms = [_][*:0]const u8{ "CA", "O", "N" };
    var radii: [3]f64 = undefined;
    var classes: [3]c_int = undefined;

    const result = freesasa_classify_atoms(
        FREESASA_CLASSIFIER_NACCESS,
        &residues,
        &atoms,
        3,
        &radii,
        &classes,
    );

    try std.testing.expectEqual(FREESASA_OK, result);

    // Check radii
    try std.testing.expectApproxEqAbs(1.87, radii[0], 0.01); // CA
    try std.testing.expectApproxEqAbs(1.40, radii[1], 0.01); // O
    try std.testing.expectApproxEqAbs(1.65, radii[2], 0.01); // N

    // Check classes
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_APOLAR, classes[0]); // CA
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_POLAR, classes[1]); // O
    try std.testing.expectEqual(FREESASA_ATOM_CLASS_POLAR, classes[2]); // N
}

test "freesasa_classify_atoms without classes" {
    const residues = [_][*:0]const u8{ "ALA", "ALA" };
    const atoms = [_][*:0]const u8{ "CA", "CB" };
    var radii: [2]f64 = undefined;

    const result = freesasa_classify_atoms(
        FREESASA_CLASSIFIER_NACCESS,
        &residues,
        &atoms,
        2,
        &radii,
        null, // classes_out is null
    );

    try std.testing.expectEqual(FREESASA_OK, result);
    try std.testing.expectApproxEqAbs(1.87, radii[0], 0.01);
    try std.testing.expectApproxEqAbs(1.87, radii[1], 0.01);
}

test "freesasa_classify_atoms empty" {
    var radii: [0]f64 = undefined;
    const residues: [*]const [*:0]const u8 = undefined;
    const atoms: [*]const [*:0]const u8 = undefined;

    const result = freesasa_classify_atoms(
        FREESASA_CLASSIFIER_NACCESS,
        residues,
        atoms,
        0,
        &radii,
        null,
    );

    try std.testing.expectEqual(FREESASA_OK, result);
}

test "freesasa_classify_atoms invalid classifier" {
    const residues = [_][*:0]const u8{"ALA"};
    const atoms = [_][*:0]const u8{"CA"};
    var radii: [1]f64 = undefined;

    const result = freesasa_classify_atoms(
        99, // Invalid classifier type
        &residues,
        &atoms,
        1,
        &radii,
        null,
    );

    try std.testing.expectEqual(FREESASA_ERROR_INVALID_INPUT, result);
}

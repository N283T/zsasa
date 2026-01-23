const std = @import("std");
const types = @import("types.zig");
const test_points = @import("test_points.zig");

const Vec3 = types.Vec3;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const Config = types.Config;
const Allocator = std.mem.Allocator;

/// Calculate SASA for a single atom using the Shrake-Rupley algorithm.
///
/// Algorithm (from FreeSASA sasa_sr.c sr_atom_area, lines 276-338):
/// 1. Scale test points by (r[i] + probe_radius)
/// 2. Translate to atom position (x[i], y[i], z[i])
/// 3. For each test point:
///    - Check if it's inside any other atom j (distance² < (r[j] + probe_radius)²)
///    - If not inside any atom, count as exposed
/// 4. SASA[i] = 4 * PI * (r[i] + probe_radius)² * (exposed_points / total_points)
///
/// # Parameters
/// - `atom_idx`: Index of the atom to calculate SASA for
/// - `positions`: Array of atom positions (Vec3)
/// - `radii`: Array of atom radii
/// - `test_points`: Array of unit sphere test points
/// - `probe_radius`: Water probe radius in Angstroms
///
/// # Returns
/// SASA value for the atom in Ų
pub fn atomSasa(
    atom_idx: usize,
    positions: []const Vec3,
    radii: []const f64,
    test_points_array: []const Vec3,
    probe_radius: f64,
) f64 {
    const n_atoms = positions.len;
    const n_points = test_points_array.len;

    // Get atom properties
    const atom_pos = positions[atom_idx];
    const atom_radius = radii[atom_idx];
    const atom_radius_probe = atom_radius + probe_radius;

    // Count exposed test points
    var n_exposed: usize = 0;

    for (test_points_array) |test_point| {
        // Scale and translate test point to atom surface
        const scaled = test_point.scale(atom_radius_probe);
        const point = atom_pos.add(scaled);

        // Check if this point is inside any other atom
        var is_buried = false;
        for (0..n_atoms) |j| {
            if (j == atom_idx) continue; // Skip self

            const other_pos = positions[j];
            const other_radius_probe = radii[j] + probe_radius;

            // Calculate squared distance
            const dx = point.x - other_pos.x;
            const dy = point.y - other_pos.y;
            const dz = point.z - other_pos.z;
            const dist_sq = dx * dx + dy * dy + dz * dz;

            // Check if point is inside this atom
            const radius_sq = other_radius_probe * other_radius_probe;
            if (dist_sq < radius_sq) {
                is_buried = true;
                break;
            }
        }

        if (!is_buried) {
            n_exposed += 1;
        }
    }

    // Calculate SASA: 4π * r² * (exposed / total)
    const surface_area = 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
    const exposed_fraction = @as(f64, @floatFromInt(n_exposed)) / @as(f64, @floatFromInt(n_points));
    return surface_area * exposed_fraction;
}

/// Calculate SASA for all atoms in the system.
///
/// # Parameters
/// - `allocator`: Memory allocator for result arrays
/// - `input`: Atom input data (positions and radii)
/// - `config`: Configuration parameters (n_points, probe_radius)
///
/// # Returns
/// SasaResult containing total_area and per-atom areas. Caller must call deinit().
pub fn calculateSasa(
    allocator: Allocator,
    input: AtomInput,
    config: Config,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) return error.NoAtoms;

    // Generate test points
    const test_points_array = try test_points.generateTestPoints(allocator, config.n_points);
    defer allocator.free(test_points_array);

    // Convert input to Vec3 positions
    const positions = try allocator.alloc(Vec3, n_atoms);
    defer allocator.free(positions);

    for (0..n_atoms) |i| {
        positions[i] = Vec3{
            .x = input.x[i],
            .y = input.y[i],
            .z = input.z[i],
        };
    }

    // Allocate result arrays
    const atom_areas = try allocator.alloc(f64, n_atoms);
    errdefer allocator.free(atom_areas);

    // Calculate SASA for each atom
    var total_area: f64 = 0.0;
    for (0..n_atoms) |i| {
        const area = atomSasa(i, positions, input.r, test_points_array, config.probe_radius);
        atom_areas[i] = area;
        total_area += area;
    }

    return SasaResult{
        .total_area = total_area,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}

// Tests

test "atomSasa - single isolated atom" {
    const allocator = std.testing.allocator;

    // Single atom at origin
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
    };
    const radii = &[_]f64{1.0};
    const probe_radius = 1.4;

    // Generate test points
    const test_points_array = try test_points.generateTestPoints(allocator, 1000);
    defer allocator.free(test_points_array);

    const sasa = atomSasa(0, positions, radii, test_points_array, probe_radius);

    // Expected: 4π * (r + probe)²
    const expected_radius = radii[0] + probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    // Should be very close with 1000 test points
    try std.testing.expectApproxEqRel(expected_area, sasa, 0.01);
}

test "atomSasa - two far apart atoms" {
    const allocator = std.testing.allocator;

    // Two atoms far apart (no overlap)
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 100.0, .y = 0.0, .z = 0.0 }, // Very far away
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    const test_points_array = try test_points.generateTestPoints(allocator, 1000);
    defer allocator.free(test_points_array);

    // Each atom should have full SASA (no burial)
    const sasa1 = atomSasa(0, positions, radii, test_points_array, probe_radius);
    const sasa2 = atomSasa(1, positions, radii, test_points_array, probe_radius);

    const expected_radius = radii[0] + probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    try std.testing.expectApproxEqRel(expected_area, sasa1, 0.01);
    try std.testing.expectApproxEqRel(expected_area, sasa2, 0.01);
}

test "atomSasa - two overlapping atoms" {
    const allocator = std.testing.allocator;

    // Two atoms close together (overlapping)
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 2.0, .y = 0.0, .z = 0.0 }, // Close enough to bury some surface
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    const test_points_array = try test_points.generateTestPoints(allocator, 1000);
    defer allocator.free(test_points_array);

    const sasa1 = atomSasa(0, positions, radii, test_points_array, probe_radius);

    // Expected: less than full SASA due to burial
    const expected_radius = radii[0] + probe_radius;
    const full_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    // SASA should be less than full area
    try std.testing.expect(sasa1 < full_area);
    // But still positive
    try std.testing.expect(sasa1 > 0.0);
}

test "calculateSasa - single atom" {
    const allocator = std.testing.allocator;

    // Create input for single atom
    const x = try allocator.alloc(f64, 1);
    const y = try allocator.alloc(f64, 1);
    const z = try allocator.alloc(f64, 1);
    const r = try allocator.alloc(f64, 1);
    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;

    var input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input.deinit();

    const config = Config{
        .n_points = 1000,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    // Expected: 4π * (r + probe)²
    const expected_radius = r[0] + config.probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    try std.testing.expectEqual(@as(usize, 1), result.atom_areas.len);
    try std.testing.expectApproxEqRel(expected_area, result.total_area, 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[0], 0.01);
}

test "calculateSasa - two far apart atoms" {
    const allocator = std.testing.allocator;

    // Create input for two atoms
    const x = try allocator.alloc(f64, 2);
    const y = try allocator.alloc(f64, 2);
    const z = try allocator.alloc(f64, 2);
    const r = try allocator.alloc(f64, 2);
    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;
    x[1] = 100.0;
    y[1] = 0.0;
    z[1] = 0.0;
    r[1] = 1.0;

    var input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input.deinit();

    const config = Config{
        .n_points = 1000,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    // Each atom should have full SASA
    const expected_radius = r[0] + config.probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;
    const expected_total = expected_area * 2.0;

    try std.testing.expectEqual(@as(usize, 2), result.atom_areas.len);
    try std.testing.expectApproxEqRel(expected_total, result.total_area, 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[0], 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[1], 0.01);
}

test "calculateSasa - two overlapping atoms" {
    const allocator = std.testing.allocator;

    // Create input for two overlapping atoms
    const x = try allocator.alloc(f64, 2);
    const y = try allocator.alloc(f64, 2);
    const z = try allocator.alloc(f64, 2);
    const r = try allocator.alloc(f64, 2);
    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;
    x[1] = 2.0;
    y[1] = 0.0;
    z[1] = 0.0;
    r[1] = 1.0;

    var input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input.deinit();

    const config = Config{
        .n_points = 1000,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    // Total should be less than 2 * full_area due to burial
    const expected_radius = r[0] + config.probe_radius;
    const full_area = 4.0 * std.math.pi * expected_radius * expected_radius;
    const full_total = full_area * 2.0;

    try std.testing.expectEqual(@as(usize, 2), result.atom_areas.len);
    try std.testing.expect(result.total_area < full_total);
    try std.testing.expect(result.total_area > 0.0);
    // Each atom should have reduced SASA
    try std.testing.expect(result.atom_areas[0] < full_area);
    try std.testing.expect(result.atom_areas[1] < full_area);
}

test "calculateSasa - no atoms error" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 0);
    const y = try allocator.alloc(f64, 0);
    const z = try allocator.alloc(f64, 0);
    const r = try allocator.alloc(f64, 0);

    var input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input.deinit();

    const config = Config{};

    const result = calculateSasa(allocator, input, config);
    try std.testing.expectError(error.NoAtoms, result);
}

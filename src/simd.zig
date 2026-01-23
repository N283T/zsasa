const std = @import("std");
const types = @import("types.zig");
const Vec3 = types.Vec3;

/// SIMD-optimized batch distance squared calculation.
/// Process 4 positions simultaneously using @Vector(4, f64).
///
/// # Parameters
/// - `point`: The test point to measure distances from
/// - `positions`: Array of 4 Vec3 positions to measure to
///
/// # Returns
/// Array of 4 squared distances
pub fn distanceSquaredBatch4(
    point: Vec3,
    positions: [4]Vec3,
) [4]f64 {
    // Splat point coordinates into vectors
    const px: @Vector(4, f64) = @splat(point.x);
    const py: @Vector(4, f64) = @splat(point.y);
    const pz: @Vector(4, f64) = @splat(point.z);

    // Load other positions into vectors
    const ox = @Vector(4, f64){ positions[0].x, positions[1].x, positions[2].x, positions[3].x };
    const oy = @Vector(4, f64){ positions[0].y, positions[1].y, positions[2].y, positions[3].y };
    const oz = @Vector(4, f64){ positions[0].z, positions[1].z, positions[2].z, positions[3].z };

    // Calculate differences
    const dx = px - ox;
    const dy = py - oy;
    const dz = pz - oz;

    // Calculate squared distances: dx² + dy² + dz²
    const dist_sq = dx * dx + dy * dy + dz * dz;

    return dist_sq;
}

/// Check if point is buried by any of 4 atoms.
///
/// # Parameters
/// - `point`: The test point to check
/// - `positions`: Array of 4 atom positions
/// - `radii_sq`: Array of 4 pre-computed (radius + probe)² values
///
/// # Returns
/// true if point is inside any of the 4 atoms, false otherwise
pub fn isPointBuriedBatch4(
    point: Vec3,
    positions: [4]Vec3,
    radii_sq: [4]f64,
) bool {
    const dist_sq = distanceSquaredBatch4(point, positions);
    const radii_v: @Vector(4, f64) = radii_sq;
    const dist_v: @Vector(4, f64) = dist_sq;

    // Check if any distance < radius (point inside atom)
    const inside = dist_v < radii_v;
    return @reduce(.Or, inside);
}

// Tests

test "distanceSquaredBatch4 - correctness" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 1, .y = 0, .z = 0 }, // dist² = 1
        Vec3{ .x = 0, .y = 2, .z = 0 }, // dist² = 4
        Vec3{ .x = 0, .y = 0, .z = 3 }, // dist² = 9
        Vec3{ .x = 1, .y = 1, .z = 1 }, // dist² = 3
    };

    const result = distanceSquaredBatch4(point, positions);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), result[3], 1e-10);
}

test "distanceSquaredBatch4 - non-origin point" {
    const point = Vec3{ .x = 1, .y = 2, .z = 3 };
    const positions = [4]Vec3{
        Vec3{ .x = 1, .y = 2, .z = 3 }, // dist² = 0
        Vec3{ .x = 2, .y = 2, .z = 3 }, // dist² = 1
        Vec3{ .x = 1, .y = 4, .z = 3 }, // dist² = 4
        Vec3{ .x = 4, .y = 6, .z = 3 }, // dist² = 9 + 16 = 25
    };

    const result = distanceSquaredBatch4(point, positions);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 25.0), result[3], 1e-10);
}

test "isPointBuriedBatch4 - none inside" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 10, .y = 0, .z = 0 }, // dist² = 100
        Vec3{ .x = 0, .y = 10, .z = 0 }, // dist² = 100
        Vec3{ .x = 0, .y = 0, .z = 10 }, // dist² = 100
        Vec3{ .x = 10, .y = 10, .z = 10 }, // dist² = 300
    };
    const radii_sq = [4]f64{ 1.0, 1.0, 1.0, 1.0 }; // All radii² = 1

    try std.testing.expect(!isPointBuriedBatch4(point, positions, radii_sq));
}

test "isPointBuriedBatch4 - one inside (first)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 0, .y = 0, .z = 0 }, // dist² = 0.25 < 1.0
        Vec3{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3{ .x = 10, .y = 0, .z = 0 }, // far
    };
    const radii_sq = [4]f64{ 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch4(point, positions, radii_sq));
}

test "isPointBuriedBatch4 - one inside (last)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3{ .x = 0, .y = 0, .z = 0 }, // dist² = 0.25 < 1.0
    };
    const radii_sq = [4]f64{ 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch4(point, positions, radii_sq));
}

test "isPointBuriedBatch4 - boundary case (exactly on radius)" {
    const point = Vec3{ .x = 1.0, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 0, .y = 0, .z = 0 }, // dist² = 1.0 == radius²
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
    };
    const radii_sq = [4]f64{ 1.0, 1.0, 1.0, 1.0 };

    // dist² == radius² means NOT inside (we use < not <=)
    try std.testing.expect(!isPointBuriedBatch4(point, positions, radii_sq));
}

test "isPointBuriedBatch4 - all inside" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 0.1, .y = 0, .z = 0 },
        Vec3{ .x = 0, .y = 0.1, .z = 0 },
        Vec3{ .x = 0, .y = 0, .z = 0.1 },
        Vec3{ .x = 0.05, .y = 0.05, .z = 0 },
    };
    const radii_sq = [4]f64{ 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch4(point, positions, radii_sq));
}

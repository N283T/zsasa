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

/// SIMD-optimized batch distance squared calculation for 8 atoms.
/// Process 8 positions simultaneously using @Vector(8, f64).
///
/// # Parameters
/// - `point`: The test point to measure distances from
/// - `positions`: Array of 8 Vec3 positions to measure to
///
/// # Returns
/// Array of 8 squared distances
pub fn distanceSquaredBatch8(
    point: Vec3,
    positions: [8]Vec3,
) [8]f64 {
    // Splat point coordinates into vectors
    const px: @Vector(8, f64) = @splat(point.x);
    const py: @Vector(8, f64) = @splat(point.y);
    const pz: @Vector(8, f64) = @splat(point.z);

    // Load other positions into vectors
    const ox = @Vector(8, f64){
        positions[0].x, positions[1].x, positions[2].x, positions[3].x,
        positions[4].x, positions[5].x, positions[6].x, positions[7].x,
    };
    const oy = @Vector(8, f64){
        positions[0].y, positions[1].y, positions[2].y, positions[3].y,
        positions[4].y, positions[5].y, positions[6].y, positions[7].y,
    };
    const oz = @Vector(8, f64){
        positions[0].z, positions[1].z, positions[2].z, positions[3].z,
        positions[4].z, positions[5].z, positions[6].z, positions[7].z,
    };

    // Calculate differences
    const dx = px - ox;
    const dy = py - oy;
    const dz = pz - oz;

    // Calculate squared distances: dx² + dy² + dz²
    const dist_sq = dx * dx + dy * dy + dz * dz;

    return dist_sq;
}

/// Check if point is buried by any of 8 atoms.
///
/// # Parameters
/// - `point`: The test point to check
/// - `positions`: Array of 8 atom positions
/// - `radii_sq`: Array of 8 pre-computed (radius + probe)² values
///
/// # Returns
/// true if point is inside any of the 8 atoms, false otherwise
pub fn isPointBuriedBatch8(
    point: Vec3,
    positions: [8]Vec3,
    radii_sq: [8]f64,
) bool {
    const dist_sq = distanceSquaredBatch8(point, positions);
    const radii_v: @Vector(8, f64) = radii_sq;
    const dist_v: @Vector(8, f64) = dist_sq;

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

test "distanceSquaredBatch8 - correctness" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [8]Vec3{
        Vec3{ .x = 1, .y = 0, .z = 0 }, // dist² = 1
        Vec3{ .x = 0, .y = 2, .z = 0 }, // dist² = 4
        Vec3{ .x = 0, .y = 0, .z = 3 }, // dist² = 9
        Vec3{ .x = 1, .y = 1, .z = 1 }, // dist² = 3
        Vec3{ .x = 2, .y = 0, .z = 0 }, // dist² = 4
        Vec3{ .x = 0, .y = 3, .z = 0 }, // dist² = 9
        Vec3{ .x = 0, .y = 0, .z = 4 }, // dist² = 16
        Vec3{ .x = 1, .y = 1, .z = 0 }, // dist² = 2
    };

    const result = distanceSquaredBatch8(point, positions);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), result[3], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), result[4], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[5], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 16.0), result[6], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), result[7], 1e-10);
}

test "isPointBuriedBatch8 - none inside" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [8]Vec3{
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 0, .y = 10, .z = 0 },
        Vec3{ .x = 0, .y = 0, .z = 10 },
        Vec3{ .x = 10, .y = 10, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 10 },
        Vec3{ .x = 0, .y = 10, .z = 10 },
        Vec3{ .x = 10, .y = 10, .z = 10 },
        Vec3{ .x = 5, .y = 5, .z = 5 },
    };
    const radii_sq = [8]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(!isPointBuriedBatch8(point, positions, radii_sq));
}

test "isPointBuriedBatch8 - one inside (first)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [8]Vec3{
        Vec3{ .x = 0, .y = 0, .z = 0 }, // dist² = 0.25 < 1.0
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
    };
    const radii_sq = [8]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch8(point, positions, radii_sq));
}

test "isPointBuriedBatch8 - one inside (last)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [8]Vec3{
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 0, .y = 0, .z = 0 }, // dist² = 0.25 < 1.0
    };
    const radii_sq = [8]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch8(point, positions, radii_sq));
}

test "isPointBuriedBatch8 - boundary case (exactly on radius)" {
    const point = Vec3{ .x = 1.0, .y = 0, .z = 0 };
    const positions = [8]Vec3{
        Vec3{ .x = 0, .y = 0, .z = 0 }, // dist² = 1.0 == radius²
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
        Vec3{ .x = 10, .y = 0, .z = 0 },
    };
    const radii_sq = [8]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    // dist² == radius² means NOT inside (we use < not <=)
    try std.testing.expect(!isPointBuriedBatch8(point, positions, radii_sq));
}

// ============================================================================
// Lee-Richards SIMD helpers
// ============================================================================

/// SIMD-optimized batch xy-distance calculation for Lee-Richards.
/// Computes sqrt(dx² + dy²) for 4 neighbors simultaneously.
///
/// # Parameters
/// - `xi`, `yi`: Coordinates of the reference atom
/// - `x_neighbors`, `y_neighbors`: Arrays of 4 neighbor x/y coordinates
///
/// # Returns
/// Array of 4 xy-distances
pub fn xyDistanceBatch4(
    xi: f64,
    yi: f64,
    x_neighbors: [4]f64,
    y_neighbors: [4]f64,
) [4]f64 {
    const px: @Vector(4, f64) = @splat(xi);
    const py: @Vector(4, f64) = @splat(yi);

    const nx: @Vector(4, f64) = x_neighbors;
    const ny: @Vector(4, f64) = y_neighbors;

    const dx = nx - px;
    const dy = ny - py;

    const dist_sq = dx * dx + dy * dy;
    const dist = @sqrt(dist_sq);

    return dist;
}

/// Compute slice radii (Rj' = sqrt(Rj² - dj²)) for 4 neighbors.
/// Returns 0 for neighbors that don't intersect the slice (dj >= Rj).
///
/// # Parameters
/// - `slice_z`: Z-coordinate of the current slice
/// - `z_neighbors`: Array of 4 neighbor z-coordinates
/// - `radii`: Array of 4 neighbor radii
///
/// # Returns
/// Array of 4 slice radii (0 if no intersection)
pub fn sliceRadiiBatch4(
    slice_z: f64,
    z_neighbors: [4]f64,
    radii: [4]f64,
) [4]f64 {
    const sz: @Vector(4, f64) = @splat(slice_z);
    const zn: @Vector(4, f64) = z_neighbors;
    const rn: @Vector(4, f64) = radii;

    const dz = zn - sz;
    const dz_sq = dz * dz;
    const r_sq = rn * rn;

    // Rj_prime² = Rj² - dj²
    const rp_sq = r_sq - dz_sq;

    // Clamp negative values to 0 (no intersection)
    const zero: @Vector(4, f64) = @splat(0.0);
    const rp_sq_clamped = @max(rp_sq, zero);

    return @sqrt(rp_sq_clamped);
}

/// Check if circles overlap (dij < Ri' + Rj') for 4 neighbors.
///
/// # Parameters
/// - `dij`: Array of 4 xy-distances
/// - `ri_prime`: Slice radius of reference atom
/// - `rj_primes`: Array of 4 neighbor slice radii
///
/// # Returns
/// Bitmask where bit i is set if circles overlap
pub fn circlesOverlapBatch4(
    dij: [4]f64,
    ri_prime: f64,
    rj_primes: [4]f64,
) u4 {
    const d: @Vector(4, f64) = dij;
    const ri: @Vector(4, f64) = @splat(ri_prime);
    const rj: @Vector(4, f64) = rj_primes;

    const sum_radii = ri + rj;
    const overlaps = d < sum_radii;

    return @bitCast(overlaps);
}

/// SIMD-optimized batch xy-distance calculation for Lee-Richards (8-wide).
/// Computes sqrt(dx² + dy²) for 8 neighbors simultaneously.
///
/// # Parameters
/// - `xi`, `yi`: Coordinates of the reference atom
/// - `x_neighbors`, `y_neighbors`: Arrays of 8 neighbor x/y coordinates
///
/// # Returns
/// Array of 8 xy-distances
pub fn xyDistanceBatch8(
    xi: f64,
    yi: f64,
    x_neighbors: [8]f64,
    y_neighbors: [8]f64,
) [8]f64 {
    const px: @Vector(8, f64) = @splat(xi);
    const py: @Vector(8, f64) = @splat(yi);

    const nx: @Vector(8, f64) = x_neighbors;
    const ny: @Vector(8, f64) = y_neighbors;

    const dx = nx - px;
    const dy = ny - py;

    const dist_sq = dx * dx + dy * dy;
    const dist = @sqrt(dist_sq);

    return dist;
}

/// Compute slice radii (Rj' = sqrt(Rj² - dj²)) for 8 neighbors.
/// Returns 0 for neighbors that don't intersect the slice (dj >= Rj).
///
/// # Parameters
/// - `slice_z`: Z-coordinate of the current slice
/// - `z_neighbors`: Array of 8 neighbor z-coordinates
/// - `radii`: Array of 8 neighbor radii
///
/// # Returns
/// Array of 8 slice radii (0 if no intersection)
pub fn sliceRadiiBatch8(
    slice_z: f64,
    z_neighbors: [8]f64,
    radii: [8]f64,
) [8]f64 {
    const sz: @Vector(8, f64) = @splat(slice_z);
    const zn: @Vector(8, f64) = z_neighbors;
    const rn: @Vector(8, f64) = radii;

    const dz = zn - sz;
    const dz_sq = dz * dz;
    const r_sq = rn * rn;

    // Rj_prime² = Rj² - dj²
    const rp_sq = r_sq - dz_sq;

    // Clamp negative values to 0 (no intersection)
    const zero: @Vector(8, f64) = @splat(0.0);
    const rp_sq_clamped = @max(rp_sq, zero);

    return @sqrt(rp_sq_clamped);
}

/// Check if circles overlap (dij < Ri' + Rj') for 8 neighbors.
///
/// # Parameters
/// - `dij`: Array of 8 xy-distances
/// - `ri_prime`: Slice radius of reference atom
/// - `rj_primes`: Array of 8 neighbor slice radii
///
/// # Returns
/// Bitmask where bit i is set if circles overlap
pub fn circlesOverlapBatch8(
    dij: [8]f64,
    ri_prime: f64,
    rj_primes: [8]f64,
) u8 {
    const d: @Vector(8, f64) = dij;
    const ri: @Vector(8, f64) = @splat(ri_prime);
    const rj: @Vector(8, f64) = rj_primes;

    const sum_radii = ri + rj;
    const overlaps = d < sum_radii;

    return @bitCast(overlaps);
}

// Lee-Richards SIMD tests

test "xyDistanceBatch4 - correctness" {
    const result = xyDistanceBatch4(
        0.0,
        0.0,
        [4]f64{ 3.0, 0.0, 1.0, 3.0 },
        [4]f64{ 4.0, 5.0, 0.0, 4.0 },
    );

    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[0], 1e-10); // 3-4-5 triangle
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[3], 1e-10);
}

test "sliceRadiiBatch4 - correctness" {
    const result = sliceRadiiBatch4(
        0.0, // slice at z=0
        [4]f64{ 0.0, 0.6, 0.8, 2.0 }, // z-coords
        [4]f64{ 1.0, 1.0, 1.0, 1.0 }, // radii (R=1)
    );

    // R' = sqrt(R² - d²)
    // Neighbor 0: sqrt(1 - 0) = 1.0
    // Neighbor 1: sqrt(1 - 0.36) = sqrt(0.64) = 0.8
    // Neighbor 2: sqrt(1 - 0.64) = sqrt(0.36) = 0.6
    // Neighbor 3: sqrt(1 - 4) -> clamped to 0
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result[3], 1e-10);
}

test "circlesOverlapBatch4 - mixed" {
    const result = circlesOverlapBatch4(
        [4]f64{ 1.0, 3.0, 0.5, 2.0 }, // xy-distances
        1.0, // Ri'
        [4]f64{ 1.0, 1.0, 1.0, 1.0 }, // Rj' values
    );

    // Ri' + Rj' = 2.0 for all
    // dij < 2.0?
    // Neighbor 0: 1.0 < 2.0 -> overlaps
    // Neighbor 1: 3.0 < 2.0 -> no
    // Neighbor 2: 0.5 < 2.0 -> overlaps
    // Neighbor 3: 2.0 < 2.0 -> no (equal, not less)
    try std.testing.expectEqual(@as(u4, 0b0101), result);
}

// Lee-Richards 8-wide SIMD tests

test "xyDistanceBatch8 - correctness" {
    const result = xyDistanceBatch8(
        0.0,
        0.0,
        [8]f64{ 3.0, 0.0, 1.0, 3.0, 4.0, 0.0, 5.0, 6.0 },
        [8]f64{ 4.0, 5.0, 0.0, 4.0, 3.0, 12.0, 12.0, 8.0 },
    );

    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[0], 1e-10); // 3-4-5 triangle
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[3], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), result[4], 1e-10); // 4-3-5
    try std.testing.expectApproxEqAbs(@as(f64, 12.0), result[5], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 13.0), result[6], 1e-10); // 5-12-13
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), result[7], 1e-10); // 6-8-10
}

test "sliceRadiiBatch8 - correctness" {
    const result = sliceRadiiBatch8(
        0.0, // slice at z=0
        [8]f64{ 0.0, 0.6, 0.8, 2.0, 0.0, 0.3, 0.4, 1.5 }, // z-coords
        [8]f64{ 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 1.0 }, // radii
    );

    // R' = sqrt(R² - d²)
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[0], 1e-10); // sqrt(1 - 0)
    try std.testing.expectApproxEqAbs(@as(f64, 0.8), result[1], 1e-10); // sqrt(1 - 0.36)
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), result[2], 1e-10); // sqrt(1 - 0.64)
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result[3], 1e-10); // clamped
    try std.testing.expectApproxEqAbs(@as(f64, 0.5), result[4], 1e-10); // sqrt(0.25 - 0)
    try std.testing.expectApproxEqAbs(@as(f64, 0.4), result[5], 1e-10); // sqrt(0.25 - 0.09)
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), result[6], 1e-10); // sqrt(0.25 - 0.16)
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), result[7], 1e-10); // clamped
}

test "circlesOverlapBatch8 - mixed" {
    const result = circlesOverlapBatch8(
        [8]f64{ 1.0, 3.0, 0.5, 2.0, 1.5, 2.5, 0.1, 1.9 }, // xy-distances
        1.0, // Ri'
        [8]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }, // Rj' values
    );

    // Ri' + Rj' = 2.0 for all
    // dij < 2.0?
    // 0: 1.0 < 2.0 -> yes (bit 0)
    // 1: 3.0 < 2.0 -> no
    // 2: 0.5 < 2.0 -> yes (bit 2)
    // 3: 2.0 < 2.0 -> no
    // 4: 1.5 < 2.0 -> yes (bit 4)
    // 5: 2.5 < 2.0 -> no
    // 6: 0.1 < 2.0 -> yes (bit 6)
    // 7: 1.9 < 2.0 -> yes (bit 7)
    try std.testing.expectEqual(@as(u8, 0b11010101), result);
}

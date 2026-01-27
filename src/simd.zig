const std = @import("std");
const types = @import("types.zig");
const Vec3 = types.Vec3;

// ============================================================================
// Generic SIMD types (inspired by astroz)
// ============================================================================

/// Generic N-wide f64 vector type
pub fn VecN(comptime N: usize) type {
    return @Vector(N, f64);
}

/// Generic N-wide bool vector type (for masks)
pub fn VecNBool(comptime N: usize) type {
    return @Vector(N, bool);
}

// ============================================================================
// Fast approximate trigonometric functions
// ============================================================================

/// Fast approximate acos using polynomial approximation.
/// Based on Handbook of Mathematical Functions (Abramowitz & Stegun).
/// Max error: ~0.0003 radians (~0.02 degrees)
///
/// # Parameters
/// - `x`: Input value in range [-1, 1]
///
/// # Returns
/// Approximate acos(x) in radians [0, π]
pub fn fastAcos(x: f64) f64 {
    // Clamp input to valid range
    const clamped = std.math.clamp(x, -1.0, 1.0);
    const abs_x = @abs(clamped);

    // Polynomial approximation for acos(x) when x >= 0
    // acos(x) ≈ sqrt(1-x) * (a0 + a1*x + a2*x² + a3*x³)
    const a0: f64 = 1.5707963267948966; // π/2
    const a1: f64 = -0.2145988016038123;
    const a2: f64 = 0.0889789874093553;
    const a3: f64 = -0.0501743046129726;

    const sqrt_term = @sqrt(1.0 - abs_x);
    const poly = a0 + abs_x * (a1 + abs_x * (a2 + abs_x * a3));
    const result = sqrt_term * poly;

    // For negative x: acos(-x) = π - acos(x)
    return if (clamped < 0) std.math.pi - result else result;
}

/// Fast approximate atan2 using polynomial approximation.
/// Based on approximation with max error ~0.0015 radians (~0.09 degrees)
///
/// # Parameters
/// - `y`: Y coordinate
/// - `x`: X coordinate
///
/// # Returns
/// Approximate atan2(y, x) in radians [-π, π]
pub fn fastAtan2(y: f64, x: f64) f64 {
    const abs_x = @abs(x);
    const abs_y = @abs(y);

    // Handle special cases
    if (abs_x < 1e-10 and abs_y < 1e-10) {
        return 0.0;
    }

    // Use the smaller ratio for better accuracy
    const swap = abs_y > abs_x;
    const ratio = if (swap) abs_x / abs_y else abs_y / abs_x;

    // Polynomial approximation for atan(r) where r = min(|y/x|, |x/y|)
    // atan(r) ≈ r * (c0 + r² * (c1 + r² * c2))
    const c0: f64 = 0.9998660373;
    const c1: f64 = -0.3302994844;
    const c2: f64 = 0.1801410321;

    const r2 = ratio * ratio;
    var atan_r = ratio * (c0 + r2 * (c1 + r2 * c2));

    // Adjust for the octant
    if (swap) {
        atan_r = std.math.pi / 2.0 - atan_r;
    }
    if (x < 0) {
        atan_r = std.math.pi - atan_r;
    }
    if (y < 0) {
        atan_r = -atan_r;
    }

    return atan_r;
}

// ============================================================================
// Generic SIMD distance calculations with FMA optimization
// ============================================================================

/// SIMD-optimized batch distance squared calculation (generic N-wide).
/// Uses FMA (fused multiply-add) for better precision and potential speedup.
///
/// # Parameters
/// - `point`: The test point to measure distances from
/// - `positions`: Array of N Vec3 positions to measure to
///
/// # Returns
/// Array of N squared distances
pub fn distanceSquaredBatchN(comptime N: usize, point: Vec3, positions: [N]Vec3) [N]f64 {
    const Vec = VecN(N);

    // Splat point coordinates into vectors
    const px: Vec = @splat(point.x);
    const py: Vec = @splat(point.y);
    const pz: Vec = @splat(point.z);

    // Load other positions into vectors (unrolled at comptime)
    var ox_arr: [N]f64 = undefined;
    var oy_arr: [N]f64 = undefined;
    var oz_arr: [N]f64 = undefined;
    inline for (0..N) |i| {
        ox_arr[i] = positions[i].x;
        oy_arr[i] = positions[i].y;
        oz_arr[i] = positions[i].z;
    }
    const ox: Vec = ox_arr;
    const oy: Vec = oy_arr;
    const oz: Vec = oz_arr;

    // Calculate differences
    const dx = px - ox;
    const dy = py - oy;
    const dz = pz - oz;

    // Calculate squared distances using FMA: dx² + dy² + dz²
    // FMA: a*b + c in one operation with better precision
    const dz_sq = dz * dz;
    const dist_sq = @mulAdd(Vec, dx, dx, @mulAdd(Vec, dy, dy, dz_sq));

    return dist_sq;
}

/// Check if point is buried by any of N atoms (generic N-wide).
///
/// # Parameters
/// - `point`: The test point to check
/// - `positions`: Array of N atom positions
/// - `radii_sq`: Array of N pre-computed (radius + probe)² values
///
/// # Returns
/// true if point is inside any of the N atoms, false otherwise
pub fn isPointBuriedBatchN(comptime N: usize, point: Vec3, positions: [N]Vec3, radii_sq: [N]f64) bool {
    const Vec = VecN(N);
    const dist_sq = distanceSquaredBatchN(N, point, positions);
    const radii_v: Vec = radii_sq;
    const dist_v: Vec = dist_sq;

    // Check if any distance < radius (point inside atom)
    const inside = dist_v < radii_v;
    return @reduce(.Or, inside);
}

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
    return distanceSquaredBatchN(4, point, positions);
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
    return isPointBuriedBatchN(4, point, positions, radii_sq);
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
    return distanceSquaredBatchN(8, point, positions);
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
    return isPointBuriedBatchN(8, point, positions, radii_sq);
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
// Lee-Richards SIMD helpers (generic N-wide with FMA)
// ============================================================================

/// SIMD-optimized batch xy-distance calculation for Lee-Richards (generic N-wide).
/// Computes sqrt(dx² + dy²) for N neighbors simultaneously using FMA.
///
/// # Parameters
/// - `xi`, `yi`: Coordinates of the reference atom
/// - `x_neighbors`, `y_neighbors`: Arrays of N neighbor x/y coordinates
///
/// # Returns
/// Array of N xy-distances
pub fn xyDistanceBatchN(comptime N: usize, xi: f64, yi: f64, x_neighbors: [N]f64, y_neighbors: [N]f64) [N]f64 {
    const Vec = VecN(N);

    const px: Vec = @splat(xi);
    const py: Vec = @splat(yi);

    const nx: Vec = x_neighbors;
    const ny: Vec = y_neighbors;

    const dx = nx - px;
    const dy = ny - py;

    // Use FMA for better precision: dx² + dy²
    const dist_sq = @mulAdd(Vec, dx, dx, dy * dy);
    const dist = @sqrt(dist_sq);

    return dist;
}

/// Compute slice radii (Rj' = sqrt(Rj² - dj²)) for N neighbors (generic N-wide).
/// Returns 0 for neighbors that don't intersect the slice (dj >= Rj).
/// Uses FMA for precision.
///
/// # Parameters
/// - `slice_z`: Z-coordinate of the current slice
/// - `z_neighbors`: Array of N neighbor z-coordinates
/// - `radii`: Array of N neighbor radii
///
/// # Returns
/// Array of N slice radii (0 if no intersection)
pub fn sliceRadiiBatchN(comptime N: usize, slice_z: f64, z_neighbors: [N]f64, radii: [N]f64) [N]f64 {
    const Vec = VecN(N);

    const sz: Vec = @splat(slice_z);
    const zn: Vec = z_neighbors;
    const rn: Vec = radii;

    const dz = zn - sz;

    // Rj_prime² = Rj² - dj² using FMA: r*r + (-dz*dz)
    const neg_dz_sq = -(dz * dz);
    const rp_sq = @mulAdd(Vec, rn, rn, neg_dz_sq);

    // Clamp negative values to 0 (no intersection)
    const zero: Vec = @splat(0.0);
    const rp_sq_clamped = @max(rp_sq, zero);

    return @sqrt(rp_sq_clamped);
}

/// Check if circles overlap (dij < Ri' + Rj') for N neighbors.
/// Returns bitmask where bit i is set if circles overlap.
pub fn circlesOverlapBatchN(comptime N: usize, comptime MaskType: type, dij: [N]f64, ri_prime: f64, rj_primes: [N]f64) MaskType {
    const Vec = VecN(N);

    const d: Vec = dij;
    const ri: Vec = @splat(ri_prime);
    const rj: Vec = rj_primes;

    const sum_radii = ri + rj;
    const overlaps = d < sum_radii;

    return @bitCast(overlaps);
}

// Convenience wrappers for 4-wide and 8-wide

pub fn xyDistanceBatch4(xi: f64, yi: f64, x_neighbors: [4]f64, y_neighbors: [4]f64) [4]f64 {
    return xyDistanceBatchN(4, xi, yi, x_neighbors, y_neighbors);
}

pub fn sliceRadiiBatch4(slice_z: f64, z_neighbors: [4]f64, radii: [4]f64) [4]f64 {
    return sliceRadiiBatchN(4, slice_z, z_neighbors, radii);
}

pub fn circlesOverlapBatch4(dij: [4]f64, ri_prime: f64, rj_primes: [4]f64) u4 {
    return circlesOverlapBatchN(4, u4, dij, ri_prime, rj_primes);
}

pub fn xyDistanceBatch8(xi: f64, yi: f64, x_neighbors: [8]f64, y_neighbors: [8]f64) [8]f64 {
    return xyDistanceBatchN(8, xi, yi, x_neighbors, y_neighbors);
}

pub fn sliceRadiiBatch8(slice_z: f64, z_neighbors: [8]f64, radii: [8]f64) [8]f64 {
    return sliceRadiiBatchN(8, slice_z, z_neighbors, radii);
}

pub fn circlesOverlapBatch8(dij: [8]f64, ri_prime: f64, rj_primes: [8]f64) u8 {
    return circlesOverlapBatchN(8, u8, dij, ri_prime, rj_primes);
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

// Fast approximation tests

test "fastAcos - accuracy" {
    // Test various values and compare with std.math.acos
    const test_values = [_]f64{ -1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0 };
    const tolerance = 0.001; // ~0.06 degrees

    for (test_values) |x| {
        const expected = std.math.acos(x);
        const actual = fastAcos(x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }
}

test "fastAcos - edge cases" {
    // Values outside [-1, 1] should be clamped
    try std.testing.expectApproxEqAbs(std.math.pi, fastAcos(-1.5), 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), fastAcos(1.5), 0.001);
}

test "fastAtan2 - accuracy" {
    // Test various angles
    const angles = [_]f64{ 0.0, 0.25, 0.5, 1.0, 2.0, 3.0 };
    const tolerance = 0.002; // ~0.1 degrees

    for (angles) |angle| {
        const y = @sin(angle);
        const x = @cos(angle);
        const expected = std.math.atan2(y, x);
        const actual = fastAtan2(y, x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }

    // Test negative quadrants
    try std.testing.expectApproxEqAbs(std.math.atan2(-1.0, 1.0), fastAtan2(-1.0, 1.0), tolerance);
    try std.testing.expectApproxEqAbs(std.math.atan2(-1.0, -1.0), fastAtan2(-1.0, -1.0), tolerance);
    try std.testing.expectApproxEqAbs(std.math.atan2(1.0, -1.0), fastAtan2(1.0, -1.0), tolerance);
}

test "fastAtan2 - edge cases" {
    // Origin should return 0
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), fastAtan2(0.0, 0.0), 0.001);

    // Axis-aligned cases
    const tolerance = 0.002;
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), fastAtan2(0.0, 1.0), tolerance);
    try std.testing.expectApproxEqAbs(std.math.pi / 2.0, fastAtan2(1.0, 0.0), tolerance);
    try std.testing.expectApproxEqAbs(std.math.pi, fastAtan2(0.0, -1.0), tolerance);
    try std.testing.expectApproxEqAbs(-std.math.pi / 2.0, fastAtan2(-1.0, 0.0), tolerance);
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

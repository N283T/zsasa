const std = @import("std");
const builtin = @import("builtin");
const types = @import("types.zig");
const Vec3 = types.Vec3;
const Vec3f32 = types.Vec3f32;
const Vec3Gen = types.Vec3Gen;

// ============================================================================
// Compile-time CPU feature detection
// ============================================================================

/// Compile-time CPU feature detection for SIMD optimization.
pub const cpu_features = struct {
    /// AVX-512F support (512-bit vectors)
    pub const has_avx512f = if (builtin.cpu.arch == .x86_64)
        std.Target.x86.featureSetHas(builtin.cpu.features, .avx512f)
    else
        false;

    /// AVX2 support (256-bit vectors)
    pub const has_avx2 = if (builtin.cpu.arch == .x86_64)
        std.Target.x86.featureSetHas(builtin.cpu.features, .avx2)
    else
        false;

    /// ARM NEON support (128-bit vectors, available on all aarch64)
    pub const has_neon = builtin.cpu.arch == .aarch64;
};

/// Optimal vector widths based on detected CPU features.
/// - AVX-512: 16 f32s or 8 f64s per vector
/// - AVX2/NEON: 8 f32s or 4 f64s per vector
/// - Fallback: 4 f32s or 2 f64s per vector
pub const optimal_vector_width = struct {
    /// Optimal f32 vector width for this CPU
    pub const f32_width: comptime_int = if (cpu_features.has_avx512f) 16 else if (cpu_features.has_avx2 or cpu_features.has_neon) 8 else 4;

    /// Optimal f64 vector width for this CPU
    pub const f64_width: comptime_int = if (cpu_features.has_avx512f) 8 else if (cpu_features.has_avx2) 4 else 2;
};

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

/// SIMD-optimized batch distance squared calculation for 16 atoms.
/// Process 16 positions simultaneously using @Vector(16, f64).
/// Best performance on AVX-512 enabled CPUs.
///
/// # Parameters
/// - `point`: The test point to measure distances from
/// - `positions`: Array of 16 Vec3 positions to measure to
///
/// # Returns
/// Array of 16 squared distances
pub fn distanceSquaredBatch16(
    point: Vec3,
    positions: [16]Vec3,
) [16]f64 {
    // Splat point coordinates into vectors
    const px: @Vector(16, f64) = @splat(point.x);
    const py: @Vector(16, f64) = @splat(point.y);
    const pz: @Vector(16, f64) = @splat(point.z);

    // Load other positions into vectors
    const ox = @Vector(16, f64){
        positions[0].x,  positions[1].x,  positions[2].x,  positions[3].x,
        positions[4].x,  positions[5].x,  positions[6].x,  positions[7].x,
        positions[8].x,  positions[9].x,  positions[10].x, positions[11].x,
        positions[12].x, positions[13].x, positions[14].x, positions[15].x,
    };
    const oy = @Vector(16, f64){
        positions[0].y,  positions[1].y,  positions[2].y,  positions[3].y,
        positions[4].y,  positions[5].y,  positions[6].y,  positions[7].y,
        positions[8].y,  positions[9].y,  positions[10].y, positions[11].y,
        positions[12].y, positions[13].y, positions[14].y, positions[15].y,
    };
    const oz = @Vector(16, f64){
        positions[0].z,  positions[1].z,  positions[2].z,  positions[3].z,
        positions[4].z,  positions[5].z,  positions[6].z,  positions[7].z,
        positions[8].z,  positions[9].z,  positions[10].z, positions[11].z,
        positions[12].z, positions[13].z, positions[14].z, positions[15].z,
    };

    // Calculate differences
    const dx = px - ox;
    const dy = py - oy;
    const dz = pz - oz;

    // Calculate squared distances: dx² + dy² + dz²
    const dist_sq = dx * dx + dy * dy + dz * dz;

    return dist_sq;
}

/// Check if point is buried by any of 16 atoms.
///
/// # Parameters
/// - `point`: The test point to check
/// - `positions`: Array of 16 atom positions
/// - `radii_sq`: Array of 16 pre-computed (radius + probe)² values
///
/// # Returns
/// true if point is inside any of the 16 atoms, false otherwise
pub fn isPointBuriedBatch16(
    point: Vec3,
    positions: [16]Vec3,
    radii_sq: [16]f64,
) bool {
    const dist_sq = distanceSquaredBatch16(point, positions);
    const radii_v: @Vector(16, f64) = radii_sq;
    const dist_v: @Vector(16, f64) = dist_sq;

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

test "distanceSquaredBatch16 - correctness" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [16]Vec3{
        Vec3{ .x = 1, .y = 0, .z = 0 }, // dist² = 1
        Vec3{ .x = 0, .y = 2, .z = 0 }, // dist² = 4
        Vec3{ .x = 0, .y = 0, .z = 3 }, // dist² = 9
        Vec3{ .x = 1, .y = 1, .z = 1 }, // dist² = 3
        Vec3{ .x = 2, .y = 0, .z = 0 }, // dist² = 4
        Vec3{ .x = 0, .y = 3, .z = 0 }, // dist² = 9
        Vec3{ .x = 0, .y = 0, .z = 4 }, // dist² = 16
        Vec3{ .x = 1, .y = 1, .z = 0 }, // dist² = 2
        Vec3{ .x = 3, .y = 0, .z = 0 }, // dist² = 9
        Vec3{ .x = 0, .y = 4, .z = 0 }, // dist² = 16
        Vec3{ .x = 0, .y = 0, .z = 5 }, // dist² = 25
        Vec3{ .x = 2, .y = 2, .z = 0 }, // dist² = 8
        Vec3{ .x = 1, .y = 2, .z = 2 }, // dist² = 9
        Vec3{ .x = 2, .y = 1, .z = 2 }, // dist² = 9
        Vec3{ .x = 2, .y = 2, .z = 1 }, // dist² = 9
        Vec3{ .x = 1, .y = 0, .z = 1 }, // dist² = 2
    };

    const result = distanceSquaredBatch16(point, positions);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), result[0], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), result[1], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 3.0), result[3], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 4.0), result[4], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[5], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 16.0), result[6], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), result[7], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[8], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 16.0), result[9], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 25.0), result[10], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 8.0), result[11], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[12], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[13], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 9.0), result[14], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 2.0), result[15], 1e-10);
}

test "isPointBuriedBatch16 - none inside" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    var positions: [16]Vec3 = undefined;
    for (0..16) |i| {
        positions[i] = Vec3{
            .x = 10.0 + @as(f64, @floatFromInt(i)),
            .y = 0,
            .z = 0,
        };
    }
    const radii_sq = [16]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(!isPointBuriedBatch16(point, positions, radii_sq));
}

test "isPointBuriedBatch16 - one inside (first)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    var positions: [16]Vec3 = undefined;
    positions[0] = Vec3{ .x = 0, .y = 0, .z = 0 }; // dist² = 0.25 < 1.0
    for (1..16) |i| {
        positions[i] = Vec3{ .x = 10, .y = 0, .z = 0 }; // far
    }
    const radii_sq = [16]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch16(point, positions, radii_sq));
}

test "isPointBuriedBatch16 - one inside (last)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    var positions: [16]Vec3 = undefined;
    for (0..15) |i| {
        positions[i] = Vec3{ .x = 10, .y = 0, .z = 0 }; // far
    }
    positions[15] = Vec3{ .x = 0, .y = 0, .z = 0 }; // dist² = 0.25 < 1.0
    const radii_sq = [16]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch16(point, positions, radii_sq));
}

test "isPointBuriedBatch16 - one inside (middle)" {
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    var positions: [16]Vec3 = undefined;
    for (0..16) |i| {
        if (i == 8) {
            positions[i] = Vec3{ .x = 0, .y = 0, .z = 0 }; // dist² = 0.25 < 1.0
        } else {
            positions[i] = Vec3{ .x = 10, .y = 0, .z = 0 }; // far
        }
    }
    const radii_sq = [16]f64{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch16(point, positions, radii_sq));
}

// ============================================================================
// Generic SIMD functions (f32/f64)
// ============================================================================

/// Generic SIMD-optimized batch distance squared calculation (4-wide).
/// Works with both f32 and f64.
pub fn distanceSquaredBatch4Gen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [4]Vec) [4]T {
            const px: @Vector(4, T) = @splat(point.x);
            const py: @Vector(4, T) = @splat(point.y);
            const pz: @Vector(4, T) = @splat(point.z);

            const ox = @Vector(4, T){ positions[0].x, positions[1].x, positions[2].x, positions[3].x };
            const oy = @Vector(4, T){ positions[0].y, positions[1].y, positions[2].y, positions[3].y };
            const oz = @Vector(4, T){ positions[0].z, positions[1].z, positions[2].z, positions[3].z };

            const dx = px - ox;
            const dy = py - oy;
            const dz = pz - oz;

            return dx * dx + dy * dy + dz * dz;
        }
    };
}

/// Generic check if point is buried by any of 4 atoms.
pub fn isPointBuriedBatch4Gen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [4]Vec, radii_sq: [4]T) bool {
            const dist_sq = distanceSquaredBatch4Gen(T).compute(point, positions);
            const radii_v: @Vector(4, T) = radii_sq;
            const dist_v: @Vector(4, T) = dist_sq;

            const inside = dist_v < radii_v;
            return @reduce(.Or, inside);
        }
    };
}

/// Generic SIMD-optimized batch distance squared calculation (8-wide).
/// Works with both f32 and f64.
pub fn distanceSquaredBatch8Gen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [8]Vec) [8]T {
            const px: @Vector(8, T) = @splat(point.x);
            const py: @Vector(8, T) = @splat(point.y);
            const pz: @Vector(8, T) = @splat(point.z);

            const ox = @Vector(8, T){
                positions[0].x, positions[1].x, positions[2].x, positions[3].x,
                positions[4].x, positions[5].x, positions[6].x, positions[7].x,
            };
            const oy = @Vector(8, T){
                positions[0].y, positions[1].y, positions[2].y, positions[3].y,
                positions[4].y, positions[5].y, positions[6].y, positions[7].y,
            };
            const oz = @Vector(8, T){
                positions[0].z, positions[1].z, positions[2].z, positions[3].z,
                positions[4].z, positions[5].z, positions[6].z, positions[7].z,
            };

            const dx = px - ox;
            const dy = py - oy;
            const dz = pz - oz;

            return dx * dx + dy * dy + dz * dz;
        }
    };
}

/// Generic check if point is buried by any of 8 atoms.
pub fn isPointBuriedBatch8Gen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [8]Vec, radii_sq: [8]T) bool {
            const dist_sq = distanceSquaredBatch8Gen(T).compute(point, positions);
            const radii_v: @Vector(8, T) = radii_sq;
            const dist_v: @Vector(8, T) = dist_sq;

            const inside = dist_v < radii_v;
            return @reduce(.Or, inside);
        }
    };
}

/// Generic SIMD-optimized batch distance squared calculation (16-wide).
/// Works with both f32 and f64. Best performance on AVX-512 enabled CPUs.
pub fn distanceSquaredBatch16Gen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [16]Vec) [16]T {
            const px: @Vector(16, T) = @splat(point.x);
            const py: @Vector(16, T) = @splat(point.y);
            const pz: @Vector(16, T) = @splat(point.z);

            const ox = @Vector(16, T){
                positions[0].x,  positions[1].x,  positions[2].x,  positions[3].x,
                positions[4].x,  positions[5].x,  positions[6].x,  positions[7].x,
                positions[8].x,  positions[9].x,  positions[10].x, positions[11].x,
                positions[12].x, positions[13].x, positions[14].x, positions[15].x,
            };
            const oy = @Vector(16, T){
                positions[0].y,  positions[1].y,  positions[2].y,  positions[3].y,
                positions[4].y,  positions[5].y,  positions[6].y,  positions[7].y,
                positions[8].y,  positions[9].y,  positions[10].y, positions[11].y,
                positions[12].y, positions[13].y, positions[14].y, positions[15].y,
            };
            const oz = @Vector(16, T){
                positions[0].z,  positions[1].z,  positions[2].z,  positions[3].z,
                positions[4].z,  positions[5].z,  positions[6].z,  positions[7].z,
                positions[8].z,  positions[9].z,  positions[10].z, positions[11].z,
                positions[12].z, positions[13].z, positions[14].z, positions[15].z,
            };

            const dx = px - ox;
            const dy = py - oy;
            const dz = pz - oz;

            return dx * dx + dy * dy + dz * dz;
        }
    };
}

/// Generic check if point is buried by any of 16 atoms.
pub fn isPointBuriedBatch16Gen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [16]Vec, radii_sq: [16]T) bool {
            const dist_sq = distanceSquaredBatch16Gen(T).compute(point, positions);
            const radii_v: @Vector(16, T) = radii_sq;
            const dist_v: @Vector(16, T) = dist_sq;

            const inside = dist_v < radii_v;
            return @reduce(.Or, inside);
        }
    };
}

/// Generic fast approximate acos using polynomial approximation.
pub fn fastAcosGen(comptime T: type) type {
    return struct {
        pub fn compute(x: T) T {
            const clamped = std.math.clamp(x, -1.0, 1.0);
            const abs_x = @abs(clamped);

            const a0: T = 1.5707963267948966; // π/2
            const a1: T = -0.2145988016038123;
            const a2: T = 0.0889789874093553;
            const a3: T = -0.0501743046129726;

            const sqrt_term = @sqrt(1.0 - abs_x);
            const poly = a0 + abs_x * (a1 + abs_x * (a2 + abs_x * a3));
            const result = sqrt_term * poly;

            return if (clamped < 0) std.math.pi - result else result;
        }
    };
}

/// Generic fast approximate atan2 using polynomial approximation.
pub fn fastAtan2Gen(comptime T: type) type {
    return struct {
        pub fn compute(y: T, x: T) T {
            const abs_x = @abs(x);
            const abs_y = @abs(y);

            const epsilon = types.Epsilon(T).trig;
            if (abs_x < epsilon and abs_y < epsilon) {
                return 0.0;
            }

            const swap = abs_y > abs_x;
            const ratio = if (swap) abs_x / abs_y else abs_y / abs_x;

            const c0: T = 0.9998660373;
            const c1: T = -0.3302994844;
            const c2: T = 0.1801410321;

            const r2 = ratio * ratio;
            var atan_r = ratio * (c0 + r2 * (c1 + r2 * c2));

            const half_pi: T = @as(T, std.math.pi) / 2.0;
            const pi: T = std.math.pi;

            if (swap) {
                atan_r = half_pi - atan_r;
            }
            if (x < 0) {
                atan_r = pi - atan_r;
            }
            if (y < 0) {
                atan_r = -atan_r;
            }

            return atan_r;
        }
    };
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

// ============================================================================
// Generic SIMD tests (f32)
// ============================================================================

test "distanceSquaredBatch4Gen f32 - correctness" {
    const point = Vec3f32{ .x = 0, .y = 0, .z = 0 };
    const positions = [4]Vec3f32{
        Vec3f32{ .x = 1, .y = 0, .z = 0 }, // dist² = 1
        Vec3f32{ .x = 0, .y = 2, .z = 0 }, // dist² = 4
        Vec3f32{ .x = 0, .y = 0, .z = 3 }, // dist² = 9
        Vec3f32{ .x = 1, .y = 1, .z = 1 }, // dist² = 3
    };

    const result = distanceSquaredBatch4Gen(f32).compute(point, positions);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), result[0], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 4.0), result[1], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 9.0), result[2], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 3.0), result[3], 1e-5);
}

test "isPointBuriedBatch4Gen f32 - one inside" {
    const point = Vec3f32{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [4]Vec3f32{
        Vec3f32{ .x = 0, .y = 0, .z = 0 }, // dist² = 0.25 < 1.0
        Vec3f32{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3f32{ .x = 10, .y = 0, .z = 0 }, // far
        Vec3f32{ .x = 10, .y = 0, .z = 0 }, // far
    };
    const radii_sq = [4]f32{ 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch4Gen(f32).compute(point, positions, radii_sq));
}

test "isPointBuriedBatch8Gen f32 - one inside" {
    const point = Vec3f32{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [8]Vec3f32{
        Vec3f32{ .x = 0, .y = 0, .z = 0 }, // dist² = 0.25 < 1.0
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
        Vec3f32{ .x = 10, .y = 0, .z = 0 },
    };
    const radii_sq = [8]f32{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch8Gen(f32).compute(point, positions, radii_sq));
}

test "distanceSquaredBatch16Gen f32 - correctness" {
    const point = Vec3f32{ .x = 0, .y = 0, .z = 0 };
    var positions: [16]Vec3f32 = undefined;
    positions[0] = Vec3f32{ .x = 1, .y = 0, .z = 0 }; // dist² = 1
    positions[1] = Vec3f32{ .x = 0, .y = 2, .z = 0 }; // dist² = 4
    positions[2] = Vec3f32{ .x = 0, .y = 0, .z = 3 }; // dist² = 9
    positions[3] = Vec3f32{ .x = 1, .y = 1, .z = 1 }; // dist² = 3
    for (4..16) |i| {
        const fi: f32 = @floatFromInt(i);
        positions[i] = Vec3f32{ .x = fi, .y = 0, .z = 0 }; // dist² = i²
    }

    const result = distanceSquaredBatch16Gen(f32).compute(point, positions);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), result[0], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 4.0), result[1], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 9.0), result[2], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 3.0), result[3], 1e-5);
    for (4..16) |i| {
        const fi: f32 = @floatFromInt(i);
        try std.testing.expectApproxEqAbs(fi * fi, result[i], 1e-5);
    }
}

test "isPointBuriedBatch16Gen f32 - one inside" {
    const point = Vec3f32{ .x = 0.5, .y = 0, .z = 0 };
    var positions: [16]Vec3f32 = undefined;
    positions[0] = Vec3f32{ .x = 0, .y = 0, .z = 0 }; // dist² = 0.25 < 1.0
    for (1..16) |i| {
        positions[i] = Vec3f32{ .x = 10, .y = 0, .z = 0 }; // far
    }
    const radii_sq = [16]f32{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch16Gen(f32).compute(point, positions, radii_sq));
}

test "isPointBuriedBatch16Gen f32 - none inside" {
    const point = Vec3f32{ .x = 0, .y = 0, .z = 0 };
    var positions: [16]Vec3f32 = undefined;
    for (0..16) |i| {
        positions[i] = Vec3f32{
            .x = 10.0 + @as(f32, @floatFromInt(i)),
            .y = 0,
            .z = 0,
        };
    }
    const radii_sq = [16]f32{ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(!isPointBuriedBatch16Gen(f32).compute(point, positions, radii_sq));
}

// ============================================================================
// CPU feature detection tests
// ============================================================================

test "cpu_features - compile-time detection works" {
    // These tests verify that the CPU feature detection compiles and
    // produces consistent values at compile time
    const has_avx512 = cpu_features.has_avx512f;
    const has_avx2 = cpu_features.has_avx2;
    const has_neon = cpu_features.has_neon;

    // On x86_64, if AVX-512 is available, AVX2 should also be available
    if (builtin.cpu.arch == .x86_64 and has_avx512) {
        try std.testing.expect(has_avx2);
    }

    // On aarch64, NEON should always be available
    if (builtin.cpu.arch == .aarch64) {
        try std.testing.expect(has_neon);
    }

    // Vector widths should be consistent with features
    if (has_avx512) {
        try std.testing.expectEqual(@as(comptime_int, 16), optimal_vector_width.f32_width);
        try std.testing.expectEqual(@as(comptime_int, 8), optimal_vector_width.f64_width);
    } else if (has_avx2 or has_neon) {
        try std.testing.expectEqual(@as(comptime_int, 8), optimal_vector_width.f32_width);
        if (has_avx2) {
            try std.testing.expectEqual(@as(comptime_int, 4), optimal_vector_width.f64_width);
        }
    }
}

test "fastAcosGen f32 - accuracy" {
    const test_values = [_]f32{ -1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0 };
    const tolerance: f32 = 0.002; // Slightly larger tolerance for f32

    for (test_values) |x| {
        const expected = std.math.acos(x);
        const actual = fastAcosGen(f32).compute(x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }
}

test "fastAtan2Gen f32 - accuracy" {
    const angles = [_]f32{ 0.0, 0.25, 0.5, 1.0, 2.0, 3.0 };
    const tolerance: f32 = 0.003;

    for (angles) |angle| {
        const y = @sin(angle);
        const x = @cos(angle);
        const expected = std.math.atan2(y, x);
        const actual = fastAtan2Gen(f32).compute(y, x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }
}

// ============================================================================
// Generic Lee-Richards SIMD helpers (f32/f64)
// ============================================================================

/// Generic SIMD-optimized batch xy-distance calculation for Lee-Richards (4-wide).
pub fn xyDistanceBatch4Gen(comptime T: type) type {
    return struct {
        pub fn compute(
            xi: T,
            yi: T,
            x_neighbors: [4]T,
            y_neighbors: [4]T,
        ) [4]T {
            const px: @Vector(4, T) = @splat(xi);
            const py: @Vector(4, T) = @splat(yi);

            const nx: @Vector(4, T) = x_neighbors;
            const ny: @Vector(4, T) = y_neighbors;

            const dx = nx - px;
            const dy = ny - py;

            const dist_sq = dx * dx + dy * dy;
            return @sqrt(dist_sq);
        }
    };
}

/// Generic compute slice radii (Rj' = sqrt(Rj² - dj²)) for 4 neighbors.
pub fn sliceRadiiBatch4Gen(comptime T: type) type {
    return struct {
        pub fn compute(
            slice_z: T,
            z_neighbors: [4]T,
            radii: [4]T,
        ) [4]T {
            const sz: @Vector(4, T) = @splat(slice_z);
            const zn: @Vector(4, T) = z_neighbors;
            const rn: @Vector(4, T) = radii;

            const dz = zn - sz;
            const dz_sq = dz * dz;
            const r_sq = rn * rn;

            const rp_sq = r_sq - dz_sq;

            const zero: @Vector(4, T) = @splat(0.0);
            const rp_sq_clamped = @max(rp_sq, zero);

            return @sqrt(rp_sq_clamped);
        }
    };
}

/// Generic check if circles overlap (dij < Ri' + Rj') for 4 neighbors.
pub fn circlesOverlapBatch4Gen(comptime T: type) type {
    return struct {
        pub fn compute(
            dij: [4]T,
            ri_prime: T,
            rj_primes: [4]T,
        ) u4 {
            const d: @Vector(4, T) = dij;
            const ri: @Vector(4, T) = @splat(ri_prime);
            const rj: @Vector(4, T) = rj_primes;

            const sum_radii = ri + rj;
            const overlaps = d < sum_radii;

            return @bitCast(overlaps);
        }
    };
}

/// Generic SIMD-optimized batch xy-distance calculation for Lee-Richards (8-wide).
pub fn xyDistanceBatch8Gen(comptime T: type) type {
    return struct {
        pub fn compute(
            xi: T,
            yi: T,
            x_neighbors: [8]T,
            y_neighbors: [8]T,
        ) [8]T {
            const px: @Vector(8, T) = @splat(xi);
            const py: @Vector(8, T) = @splat(yi);

            const nx: @Vector(8, T) = x_neighbors;
            const ny: @Vector(8, T) = y_neighbors;

            const dx = nx - px;
            const dy = ny - py;

            const dist_sq = dx * dx + dy * dy;
            return @sqrt(dist_sq);
        }
    };
}

/// Generic compute slice radii (Rj' = sqrt(Rj² - dj²)) for 8 neighbors.
pub fn sliceRadiiBatch8Gen(comptime T: type) type {
    return struct {
        pub fn compute(
            slice_z: T,
            z_neighbors: [8]T,
            radii: [8]T,
        ) [8]T {
            const sz: @Vector(8, T) = @splat(slice_z);
            const zn: @Vector(8, T) = z_neighbors;
            const rn: @Vector(8, T) = radii;

            const dz = zn - sz;
            const dz_sq = dz * dz;
            const r_sq = rn * rn;

            const rp_sq = r_sq - dz_sq;

            const zero: @Vector(8, T) = @splat(0.0);
            const rp_sq_clamped = @max(rp_sq, zero);

            return @sqrt(rp_sq_clamped);
        }
    };
}

/// Generic check if circles overlap (dij < Ri' + Rj') for 8 neighbors.
pub fn circlesOverlapBatch8Gen(comptime T: type) type {
    return struct {
        pub fn compute(
            dij: [8]T,
            ri_prime: T,
            rj_primes: [8]T,
        ) u8 {
            const d: @Vector(8, T) = dij;
            const ri: @Vector(8, T) = @splat(ri_prime);
            const rj: @Vector(8, T) = rj_primes;

            const sum_radii = ri + rj;
            const overlaps = d < sum_radii;

            return @bitCast(overlaps);
        }
    };
}

// Generic Lee-Richards SIMD tests

test "xyDistanceBatch4Gen f32 - correctness" {
    const result = xyDistanceBatch4Gen(f32).compute(
        0.0,
        0.0,
        [4]f32{ 3.0, 0.0, 1.0, 3.0 },
        [4]f32{ 4.0, 5.0, 0.0, 4.0 },
    );

    try std.testing.expectApproxEqAbs(@as(f32, 5.0), result[0], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 5.0), result[1], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), result[2], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 5.0), result[3], 1e-5);
}

test "sliceRadiiBatch4Gen f32 - correctness" {
    const result = sliceRadiiBatch4Gen(f32).compute(
        0.0,
        [4]f32{ 0.0, 0.6, 0.8, 2.0 },
        [4]f32{ 1.0, 1.0, 1.0, 1.0 },
    );

    try std.testing.expectApproxEqAbs(@as(f32, 1.0), result[0], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.8), result[1], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.6), result[2], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.0), result[3], 1e-5);
}

test "sliceRadiiBatch8Gen f32 - correctness" {
    const result = sliceRadiiBatch8Gen(f32).compute(
        0.0,
        [8]f32{ 0.0, 0.6, 0.8, 2.0, 0.0, 0.3, 0.4, 1.5 },
        [8]f32{ 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 1.0 },
    );

    try std.testing.expectApproxEqAbs(@as(f32, 1.0), result[0], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.8), result[1], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.6), result[2], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.0), result[3], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.5), result[4], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.4), result[5], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.3), result[6], 1e-5);
    try std.testing.expectApproxEqAbs(@as(f32, 0.0), result[7], 1e-5);
}

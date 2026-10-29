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
/// - AVX-512: 16 f32s or 8 f64s per vector (512-bit)
/// - AVX2: 8 f32s or 4 f64s per vector (256-bit)
/// - NEON: 8 f32s or 2 f64s per vector (128-bit, but efficient for f32)
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
/// The coefficients are the first four of the eight-term formula 4.4.46 in the
/// Handbook of Mathematical Functions (Abramowitz & Stegun); the four terms
/// left out cost most of its accuracy.
/// Max error, measured over [-1, 1]: 0.0042 radians (0.24 degrees), reached at
/// |x| = 0.88. The result is too small for x > 0 and too large for x < 0.
///
/// Lee-Richards uses it only in `TrigMode.fast`; see `lee_richards.TrigMode`.
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
/// The coefficients are the first three of the five-term formula 4.4.49 in the
/// Handbook of Mathematical Functions (Abramowitz & Stegun); the two terms left
/// out cost most of its accuracy.
/// Max error, measured over all directions: 0.064 radians (3.7 degrees),
/// reached on the diagonals |y| = |x|. The error grows with r = min(|y/x|, |x/y|):
/// below 0.001 rad for r < 0.5, 0.010 rad at r = 0.75, 0.033 rad at r = 0.9.
/// The result is not continuous: it jumps by 0.129 rad across each diagonal
/// (0.850 just below y = x, 0.721 just above, against π/4 = 0.785).
///
/// Lee-Richards uses it only in `TrigMode.fast`; see `lee_richards.TrigMode`.
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
// Shrake-Rupley SIMD helpers
// ============================================================================
//
// Each helper is written once, generic over the float type `T` (f32 or f64)
// and the vector width `N`. The width-specific names further down
// (`distanceSquaredBatch4`, `isPointBuriedBatch8Gen`, ...) are aliases of it.

/// SIMD batch squared-distance calculation: dx² + dy² + dz² from `point` to
/// each of `N` positions, computed in one @Vector(N, T).
pub fn DistanceSquaredBatch(comptime T: type, comptime N: usize) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [N]Vec) [N]T {
            // Splat point coordinates into vectors
            const px: @Vector(N, T) = @splat(point.x);
            const py: @Vector(N, T) = @splat(point.y);
            const pz: @Vector(N, T) = @splat(point.z);

            // Load other positions into vectors
            var ox: @Vector(N, T) = undefined;
            var oy: @Vector(N, T) = undefined;
            var oz: @Vector(N, T) = undefined;
            inline for (0..N) |i| {
                ox[i] = positions[i].x;
                oy[i] = positions[i].y;
                oz[i] = positions[i].z;
            }

            // Calculate differences
            const dx = px - ox;
            const dy = py - oy;
            const dz = pz - oz;

            // Calculate squared distances: dx² + dy² + dz²
            return dx * dx + dy * dy + dz * dz;
        }
    };
}

/// Check if a point is buried by any of `N` atoms.
///
/// `compute(point, positions, radii_sq)` returns true if the squared distance
/// from `point` to some position is below its pre-computed (radius + probe)²
/// in `radii_sq`. A point exactly on the radius is not buried.
pub fn IsPointBuriedBatch(comptime T: type, comptime N: usize) type {
    const Vec = Vec3Gen(T);
    return struct {
        pub fn compute(point: Vec, positions: [N]Vec, radii_sq: [N]T) bool {
            const dist_sq = DistanceSquaredBatch(T, N).compute(point, positions);
            const radii_v: @Vector(N, T) = radii_sq;
            const dist_v: @Vector(N, T) = dist_sq;

            // Check if any distance < radius (point inside atom)
            const inside = dist_v < radii_v;
            return @reduce(.Or, inside);
        }
    };
}

/// Squared distances from `point` to 4 positions (f64).
pub const distanceSquaredBatch4 = DistanceSquaredBatch(f64, 4).compute;
/// Squared distances from `point` to 8 positions (f64).
pub const distanceSquaredBatch8 = DistanceSquaredBatch(f64, 8).compute;
/// Squared distances from `point` to 16 positions (f64). Best on AVX-512.
pub const distanceSquaredBatch16 = DistanceSquaredBatch(f64, 16).compute;

/// Whether `point` is buried by any of 4 atoms (f64).
pub const isPointBuriedBatch4 = IsPointBuriedBatch(f64, 4).compute;
/// Whether `point` is buried by any of 8 atoms (f64).
pub const isPointBuriedBatch8 = IsPointBuriedBatch(f64, 8).compute;
/// Whether `point` is buried by any of 16 atoms (f64). Best on AVX-512.
pub const isPointBuriedBatch16 = IsPointBuriedBatch(f64, 16).compute;

/// `DistanceSquaredBatch(T, 4)`.
pub fn distanceSquaredBatch4Gen(comptime T: type) type {
    return DistanceSquaredBatch(T, 4);
}
/// `DistanceSquaredBatch(T, 8)`.
pub fn distanceSquaredBatch8Gen(comptime T: type) type {
    return DistanceSquaredBatch(T, 8);
}
/// `DistanceSquaredBatch(T, 16)`.
pub fn distanceSquaredBatch16Gen(comptime T: type) type {
    return DistanceSquaredBatch(T, 16);
}

/// `IsPointBuriedBatch(T, 4)`.
pub fn isPointBuriedBatch4Gen(comptime T: type) type {
    return IsPointBuriedBatch(T, 4);
}
/// `IsPointBuriedBatch(T, 8)`.
pub fn isPointBuriedBatch8Gen(comptime T: type) type {
    return IsPointBuriedBatch(T, 8);
}
/// `IsPointBuriedBatch(T, 16)`.
pub fn isPointBuriedBatch16Gen(comptime T: type) type {
    return IsPointBuriedBatch(T, 16);
}

// ============================================================================
// Generic fast approximate trigonometric functions (f32/f64)
// ============================================================================

/// Generic fast approximate acos using polynomial approximation.
/// Same polynomial as `fastAcos`. Max error, measured over [-1, 1] for f32 and
/// f64: 0.0042 radians (0.24 degrees).
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
/// Same polynomial as `fastAtan2`. Max error, measured over all directions for
/// f32 and f64: 0.064 radians (3.7 degrees), with the same jump of 0.129 rad
/// across the diagonals |y| = |x|.
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
//
// Written once, generic over the float type `T` and the vector width `N`, as
// above; the width-specific names further down are aliases.

/// SIMD batch xy-distance calculation for Lee-Richards.
///
/// `compute(xi, yi, x_neighbors, y_neighbors)` returns sqrt(dx² + dy²) from
/// the reference atom at (`xi`, `yi`) to each of `N` neighbors.
pub fn XyDistanceBatch(comptime T: type, comptime N: usize) type {
    return struct {
        pub fn compute(
            xi: T,
            yi: T,
            x_neighbors: [N]T,
            y_neighbors: [N]T,
        ) [N]T {
            const px: @Vector(N, T) = @splat(xi);
            const py: @Vector(N, T) = @splat(yi);

            const nx: @Vector(N, T) = x_neighbors;
            const ny: @Vector(N, T) = y_neighbors;

            const dx = nx - px;
            const dy = ny - py;

            const dist_sq = dx * dx + dy * dy;
            const result: @Vector(N, T) = @sqrt(dist_sq);
            return result;
        }
    };
}

/// Slice radii (Rj' = sqrt(Rj² - dj²)) for `N` neighbors.
///
/// `compute(slice_z, z_neighbors, radii)` returns 0 for a neighbor that does
/// not intersect the slice at `slice_z` (dj >= Rj).
pub fn SliceRadiiBatch(comptime T: type, comptime N: usize) type {
    return struct {
        pub fn compute(
            slice_z: T,
            z_neighbors: [N]T,
            radii: [N]T,
        ) [N]T {
            const sz: @Vector(N, T) = @splat(slice_z);
            const zn: @Vector(N, T) = z_neighbors;
            const rn: @Vector(N, T) = radii;

            const dz = zn - sz;
            const dz_sq = dz * dz;
            const r_sq = rn * rn;

            // Rj_prime² = Rj² - dj²
            const rp_sq = r_sq - dz_sq;

            // Clamp negative values to 0 (no intersection)
            const zero: @Vector(N, T) = @splat(0.0);
            const rp_sq_clamped = @max(rp_sq, zero);

            const result: @Vector(N, T) = @sqrt(rp_sq_clamped);
            return result;
        }
    };
}

/// Check if circles overlap (dij < Ri' + Rj') for `N` neighbors.
///
/// `compute(dij, ri_prime, rj_primes)` returns a bitmask in which bit i is set
/// if the circle of neighbor i overlaps the one of the reference atom.
pub fn CirclesOverlapBatch(comptime T: type, comptime N: usize) type {
    return struct {
        pub fn compute(
            dij: [N]T,
            ri_prime: T,
            rj_primes: [N]T,
        ) std.meta.Int(.unsigned, N) {
            const d: @Vector(N, T) = dij;
            const ri: @Vector(N, T) = @splat(ri_prime);
            const rj: @Vector(N, T) = rj_primes;

            const sum_radii = ri + rj;
            const overlaps = d < sum_radii;

            // Cross-platform packed bit layout: @Vector(N, bool) layout is target-defined
            // in Zig 0.16+. Forcing through @Vector(N, u1) before @bitCast guarantees a
            // packed uN result on all platforms (Linux x86_64 regression caught in CI).
            return @bitCast(@as(@Vector(N, u1), @intFromBool(overlaps)));
        }
    };
}

/// xy-distances to 4 neighbors (f64).
pub const xyDistanceBatch4 = XyDistanceBatch(f64, 4).compute;
/// xy-distances to 8 neighbors (f64).
pub const xyDistanceBatch8 = XyDistanceBatch(f64, 8).compute;
/// Slice radii of 4 neighbors (f64).
pub const sliceRadiiBatch4 = SliceRadiiBatch(f64, 4).compute;
/// Slice radii of 8 neighbors (f64).
pub const sliceRadiiBatch8 = SliceRadiiBatch(f64, 8).compute;
/// Overlap bitmask (u4) of 4 neighbors (f64).
pub const circlesOverlapBatch4 = CirclesOverlapBatch(f64, 4).compute;
/// Overlap bitmask (u8) of 8 neighbors (f64).
pub const circlesOverlapBatch8 = CirclesOverlapBatch(f64, 8).compute;

/// `XyDistanceBatch(T, 4)`.
pub fn xyDistanceBatch4Gen(comptime T: type) type {
    return XyDistanceBatch(T, 4);
}
/// `XyDistanceBatch(T, 8)`.
pub fn xyDistanceBatch8Gen(comptime T: type) type {
    return XyDistanceBatch(T, 8);
}
/// `SliceRadiiBatch(T, 4)`.
pub fn sliceRadiiBatch4Gen(comptime T: type) type {
    return SliceRadiiBatch(T, 4);
}
/// `SliceRadiiBatch(T, 8)`.
pub fn sliceRadiiBatch8Gen(comptime T: type) type {
    return SliceRadiiBatch(T, 8);
}
/// `CirclesOverlapBatch(T, 4)`.
pub fn circlesOverlapBatch4Gen(comptime T: type) type {
    return CirclesOverlapBatch(T, 4);
}
/// `CirclesOverlapBatch(T, 8)`.
pub fn circlesOverlapBatch8Gen(comptime T: type) type {
    return CirclesOverlapBatch(T, 8);
}

// ============================================================================
// Batch helper tests, run for every width and both float types
// ============================================================================

const float_types = .{ f64, f32 };
const sr_widths = .{ 4, 8, 16 };
const lr_widths = .{ 4, 8 };

fn testTolerance(comptime T: type) T {
    return if (T == f64) 1e-10 else 1e-5;
}

fn testVec(comptime T: type, x: f64, y: f64, z: f64) Vec3Gen(T) {
    return .{ .x = @floatCast(x), .y = @floatCast(y), .z = @floatCast(z) };
}

/// The first `n` values of `values`, as `T`.
fn testArray(comptime T: type, comptime n: usize, values: anytype) [n]T {
    var result: [n]T = undefined;
    for (&result, values[0..n]) |*r, v| r.* = @floatCast(v);
    return result;
}

/// A batch of `n` atoms at (10, 0, 0), far from points near the origin, all
/// with radius² 1.
fn FarBatch(comptime T: type, comptime n: usize) type {
    return struct {
        positions: [n]Vec3Gen(T) = @splat(.{ .x = 10, .y = 0, .z = 0 }),
        radii_sq: [n]T = @splat(1.0),
    };
}

// Positions around the origin and their squared distances to it. A width-n
// test uses the first n.
const test_dist_sq_positions = [16][3]f64{
    .{ 1, 0, 0 }, .{ 0, 2, 0 }, .{ 0, 0, 3 }, .{ 1, 1, 1 },
    .{ 2, 0, 0 }, .{ 0, 3, 0 }, .{ 0, 0, 4 }, .{ 1, 1, 0 },
    .{ 3, 0, 0 }, .{ 0, 4, 0 }, .{ 0, 0, 5 }, .{ 2, 2, 0 },
    .{ 1, 2, 2 }, .{ 2, 1, 2 }, .{ 2, 2, 1 }, .{ 1, 0, 1 },
};
const test_dist_sq_expected = [16]f64{ 1, 4, 9, 3, 4, 9, 16, 2, 9, 16, 25, 8, 9, 9, 9, 2 };

test "DistanceSquaredBatch - correctness" {
    inline for (float_types) |T| {
        inline for (sr_widths) |n| {
            var positions: [n]Vec3Gen(T) = undefined;
            for (&positions, test_dist_sq_positions[0..n]) |*p, c| p.* = testVec(T, c[0], c[1], c[2]);

            const result = DistanceSquaredBatch(T, n).compute(testVec(T, 0, 0, 0), positions);
            for (result, test_dist_sq_expected[0..n]) |actual, expected| {
                try std.testing.expectApproxEqAbs(@as(T, @floatCast(expected)), actual, testTolerance(T));
            }
        }
    }
}

test "DistanceSquaredBatch - non-origin point" {
    inline for (float_types) |T| {
        const positions = [4]Vec3Gen(T){
            testVec(T, 1, 2, 3), // dist² = 0
            testVec(T, 2, 2, 3), // dist² = 1
            testVec(T, 1, 4, 3), // dist² = 4
            testVec(T, 4, 6, 3), // dist² = 9 + 16 = 25
        };

        const result = distanceSquaredBatch4Gen(T).compute(testVec(T, 1, 2, 3), positions);
        const expected = [4]T{ 0, 1, 4, 25 };
        for (result, expected) |a, e| try std.testing.expectApproxEqAbs(e, a, testTolerance(T));
    }
}

test "IsPointBuriedBatch - none inside" {
    inline for (float_types) |T| {
        inline for (sr_widths) |n| {
            var b: FarBatch(T, n) = .{};
            for (&b.positions, 0..) |*p, i| p.* = testVec(T, 10.0 + @as(f64, @floatFromInt(i)), 0, 0);

            try std.testing.expect(!IsPointBuriedBatch(T, n).compute(testVec(T, 0, 0, 0), b.positions, b.radii_sq));
        }
    }
}

test "IsPointBuriedBatch - one inside, in each lane" {
    inline for (float_types) |T| {
        inline for (sr_widths) |n| {
            for (0..n) |lane| {
                var b: FarBatch(T, n) = .{};
                b.positions[lane] = testVec(T, 0, 0, 0); // dist² = 0.25 < 1.0

                try std.testing.expect(IsPointBuriedBatch(T, n).compute(testVec(T, 0.5, 0, 0), b.positions, b.radii_sq));
            }
        }
    }
}

test "IsPointBuriedBatch - boundary case (exactly on radius)" {
    inline for (float_types) |T| {
        inline for (sr_widths) |n| {
            for (0..n) |lane| {
                var b: FarBatch(T, n) = .{};
                b.positions[lane] = testVec(T, 0, 0, 0); // dist² = 1.0 == radius²

                // dist² == radius² means NOT inside (we use < not <=)
                try std.testing.expect(!IsPointBuriedBatch(T, n).compute(testVec(T, 1, 0, 0), b.positions, b.radii_sq));
            }
        }
    }
}

test "IsPointBuriedBatch - all inside" {
    inline for (float_types) |T| {
        inline for (sr_widths) |n| {
            var b: FarBatch(T, n) = .{};
            for (&b.positions, 0..) |*p, i| p.* = switch (i % 3) {
                0 => testVec(T, 0.1, 0, 0),
                1 => testVec(T, 0, 0.1, 0),
                else => testVec(T, 0.05, 0.05, 0),
            };

            try std.testing.expect(IsPointBuriedBatch(T, n).compute(testVec(T, 0, 0, 0), b.positions, b.radii_sq));
        }
    }
}

test "width-specific aliases" {
    const d4 = distanceSquaredBatch4(Vec3{ .x = 1, .y = 2, .z = 3 }, @splat(Vec3{ .x = 4, .y = 6, .z = 3 }));
    const d8 = distanceSquaredBatch8(Vec3{ .x = 1, .y = 2, .z = 3 }, @splat(Vec3{ .x = 4, .y = 6, .z = 3 }));
    const d16 = distanceSquaredBatch16(Vec3{ .x = 1, .y = 2, .z = 3 }, @splat(Vec3{ .x = 4, .y = 6, .z = 3 }));
    try std.testing.expectApproxEqAbs(@as(f64, 25.0), d4[3], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 25.0), d8[7], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 25.0), d16[15], 1e-10);

    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    var b4: FarBatch(f64, 4) = .{};
    var b8: FarBatch(f64, 8) = .{};
    var b16: FarBatch(f64, 16) = .{};
    try std.testing.expect(!isPointBuriedBatch4(point, b4.positions, b4.radii_sq));
    try std.testing.expect(!isPointBuriedBatch8(point, b8.positions, b8.radii_sq));
    try std.testing.expect(!isPointBuriedBatch16(point, b16.positions, b16.radii_sq));
    b4.positions[3] = Vec3{ .x = 0, .y = 0, .z = 0 };
    b8.positions[7] = Vec3{ .x = 0, .y = 0, .z = 0 };
    b16.positions[15] = Vec3{ .x = 0, .y = 0, .z = 0 };
    try std.testing.expect(isPointBuriedBatch4(point, b4.positions, b4.radii_sq));
    try std.testing.expect(isPointBuriedBatch8(point, b8.positions, b8.radii_sq));
    try std.testing.expect(isPointBuriedBatch16(point, b16.positions, b16.radii_sq));
}

// Lee-Richards inputs for 8 neighbors and the expected results. A width-n test
// uses the first n.
const test_xy_x = [8]f64{ 3.0, 0.0, 1.0, 3.0, 4.0, 0.0, 5.0, 6.0 };
const test_xy_y = [8]f64{ 4.0, 5.0, 0.0, 4.0, 3.0, 12.0, 12.0, 8.0 };
const test_xy_expected = [8]f64{ 5.0, 5.0, 1.0, 5.0, 5.0, 12.0, 13.0, 10.0 }; // 3-4-5, 5-12-13, 6-8-10

const test_slice_z = [8]f64{ 0.0, 0.6, 0.8, 2.0, 0.0, 0.3, 0.4, 1.5 };
const test_slice_r = [8]f64{ 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 1.0 };
// R' = sqrt(R² - d²), clamped to 0 where d >= R (lanes 3 and 7)
const test_slice_expected = [8]f64{ 1.0, 0.8, 0.6, 0.0, 0.5, 0.4, 0.3, 0.0 };

// With Ri' = Rj' = 1, a circle overlaps if dij < 2.0; 2.0 itself does not.
const test_overlap_dij = [8]f64{ 1.0, 3.0, 0.5, 2.0, 1.5, 2.5, 0.1, 1.9 };
const test_overlap_expected: u8 = 0b11010101;

test "XyDistanceBatch - correctness" {
    inline for (float_types) |T| {
        inline for (lr_widths) |n| {
            const result = XyDistanceBatch(T, n).compute(0.0, 0.0, testArray(T, n, test_xy_x), testArray(T, n, test_xy_y));
            for (result, testArray(T, n, test_xy_expected)) |a, e| try std.testing.expectApproxEqAbs(e, a, testTolerance(T));
        }
    }
    const r4 = xyDistanceBatch4(0.0, 0.0, test_xy_x[0..4].*, test_xy_y[0..4].*);
    const r8 = xyDistanceBatch8(0.0, 0.0, test_xy_x, test_xy_y);
    try std.testing.expectApproxEqAbs(@as(f64, 5.0), r4[3], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), r8[7], 1e-10);
}

test "SliceRadiiBatch - correctness" {
    inline for (float_types) |T| {
        inline for (lr_widths) |n| {
            const result = SliceRadiiBatch(T, n).compute(0.0, testArray(T, n, test_slice_z), testArray(T, n, test_slice_r));
            for (result, testArray(T, n, test_slice_expected)) |a, e| try std.testing.expectApproxEqAbs(e, a, testTolerance(T));
        }
    }
    const r4 = sliceRadiiBatch4(0.0, test_slice_z[0..4].*, test_slice_r[0..4].*);
    const r8 = sliceRadiiBatch8(0.0, test_slice_z, test_slice_r);
    try std.testing.expectApproxEqAbs(@as(f64, 0.6), r4[2], 1e-10);
    try std.testing.expectApproxEqAbs(@as(f64, 0.3), r8[6], 1e-10);
}

test "CirclesOverlapBatch - mixed" {
    inline for (float_types) |T| {
        inline for (lr_widths) |n| {
            const Mask = std.meta.Int(.unsigned, n);
            const result = CirclesOverlapBatch(T, n).compute(testArray(T, n, test_overlap_dij), 1.0, @splat(1.0));
            try std.testing.expectEqual(@as(Mask, @truncate(test_overlap_expected)), result);
        }
    }
    const ones: [8]f64 = @splat(1.0);
    try std.testing.expectEqual(@as(u4, 0b0101), circlesOverlapBatch4(test_overlap_dij[0..4].*, 1.0, ones[0..4].*));
    try std.testing.expectEqual(@as(u8, 0b11010101), circlesOverlapBatch8(test_overlap_dij, 1.0, ones));
}

// Fast approximation tests

test "fastAcos - accuracy" {
    // Test various values and compare with std.math.acos
    const test_values = [_]f64{ -1.0, -0.9, -0.5, 0.0, 0.5, 0.9, 1.0 };
    const tolerance = 0.005; // ~0.3 degrees, just above the maximum error of 0.0042

    for (test_values) |x| {
        const expected = std.math.acos(x);
        const actual = fastAcos(x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }
}

test "fastAcos - maximum error is the documented 0.0042 rad, for f64 and f32" {
    var max_err: f64 = 0.0;
    var max_err_gen: f64 = 0.0;
    var max_err_f32: f64 = 0.0;
    const n = 20000;
    for (0..n + 1) |i| {
        const x = -1.0 + 2.0 * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        max_err = @max(max_err, @abs(fastAcos(x) - std.math.acos(x)));
        max_err_gen = @max(max_err_gen, @abs(fastAcosGen(f64).compute(x) - std.math.acos(x)));
        const x32: f32 = @floatCast(x);
        max_err_f32 = @max(max_err_f32, @abs(@as(f64, fastAcosGen(f32).compute(x32)) - std.math.acos(@as(f64, x32))));
    }
    for ([_]f64{ max_err, max_err_gen, max_err_f32 }) |err| {
        try std.testing.expect(err > 0.0041);
        try std.testing.expect(err < 0.0043);
    }
}

test "fastAcos - edge cases" {
    // Values outside [-1, 1] should be clamped
    try std.testing.expectApproxEqAbs(std.math.pi, fastAcos(-1.5), 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 0.0), fastAcos(1.5), 0.001);
}

test "fastAtan2 - accuracy" {
    // Test various angles, all at least 0.2 rad away from the diagonals
    const angles = [_]f64{ 0.0, 0.25, 0.5, 1.0, 2.0, 3.0 };
    const tolerance = 0.005; // ~0.3 degrees; the error is larger near the diagonals

    for (angles) |angle| {
        const y = @sin(angle);
        const x = @cos(angle);
        const expected = std.math.atan2(y, x);
        const actual = fastAtan2(y, x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }

    // Test negative quadrants on the diagonals |y| = |x|, where the error is
    // at its maximum of 0.064 rad
    const neg_tolerance = 0.07;
    try std.testing.expectApproxEqAbs(std.math.atan2(@as(f64, -1.0), @as(f64, 1.0)), fastAtan2(-1.0, 1.0), neg_tolerance);
    try std.testing.expectApproxEqAbs(std.math.atan2(@as(f64, -1.0), @as(f64, -1.0)), fastAtan2(-1.0, -1.0), neg_tolerance);
    try std.testing.expectApproxEqAbs(std.math.atan2(@as(f64, 1.0), @as(f64, -1.0)), fastAtan2(1.0, -1.0), neg_tolerance);
}

test "fastAtan2 - maximum error is the documented 0.064 rad on the diagonals, for f64 and f32" {
    var max_err: f64 = 0.0;
    var max_err_gen: f64 = 0.0;
    var max_err_f32: f64 = 0.0;
    const n = 20000; // a multiple of 8, so the diagonals are sampled
    for (1..n) |i| {
        const theta = -std.math.pi + 2.0 * std.math.pi * @as(f64, @floatFromInt(i)) / @as(f64, @floatFromInt(n));
        const y = @sin(theta);
        const x = @cos(theta);
        const expected = std.math.atan2(y, x);
        max_err = @max(max_err, @abs(fastAtan2(y, x) - expected));
        max_err_gen = @max(max_err_gen, @abs(fastAtan2Gen(f64).compute(y, x) - expected));
        const y32: f32 = @floatCast(y);
        const x32: f32 = @floatCast(x);
        const expected32 = std.math.atan2(@as(f64, y32), @as(f64, x32));
        max_err_f32 = @max(max_err_f32, @abs(@as(f64, fastAtan2Gen(f32).compute(y32, x32)) - expected32));
    }
    for ([_]f64{ max_err, max_err_gen, max_err_f32 }) |err| {
        try std.testing.expect(err > 0.064);
        try std.testing.expect(err < 0.065);
    }

    // The polynomial is evaluated on min(|y/x|, |x/y|), so the result jumps
    // across a diagonal instead of passing through π/4.
    const below = fastAtan2(1.0, 1.0);
    const above = fastAtan2(1.0 + 1e-12, 1.0);
    try std.testing.expectApproxEqAbs(@as(f64, 0.8497), below, 1e-4);
    try std.testing.expectApproxEqAbs(@as(f64, 0.7211), above, 1e-4);
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
    const tolerance: f32 = 0.005; // Just above the maximum error of 0.0042

    for (test_values) |x| {
        const expected = std.math.acos(x);
        const actual = fastAcosGen(f32).compute(x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }
}

test "fastAtan2Gen f32 - accuracy" {
    const angles = [_]f32{ 0.0, 0.25, 0.5, 1.0, 2.0, 3.0 };
    const tolerance: f32 = 0.005;

    for (angles) |angle| {
        const y = @sin(angle);
        const x = @cos(angle);
        const expected = std.math.atan2(y, x);
        const actual = fastAtan2Gen(f32).compute(y, x);
        try std.testing.expectApproxEqAbs(expected, actual, tolerance);
    }
}

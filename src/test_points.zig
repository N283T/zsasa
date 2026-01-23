const std = @import("std");
const types = @import("types.zig");
const Vec3 = types.Vec3;
const Allocator = std.mem.Allocator;

/// Generate evenly distributed test points on a unit sphere using the Golden Section Spiral algorithm.
/// This follows the algorithm from FreeSASA sasa_sr.c lines 56-90.
///
/// The Golden Section Spiral distributes points evenly by using the golden angle
/// to separate successive points in longitude, and evenly spacing them in z-coordinate.
///
/// # Parameters
/// - `allocator`: Memory allocator for the points array
/// - `n`: Number of points to generate
///
/// # Returns
/// A slice of Vec3 points on the unit sphere. Caller must free the returned slice.
pub fn generateTestPoints(allocator: Allocator, n: u32) ![]Vec3 {
    if (n == 0) return error.InvalidPointCount;

    const points = try allocator.alloc(Vec3, n);
    errdefer allocator.free(points);

    // Golden angle: 2π * (3 - √5) / 2 = π * (3 - √5)
    const dlong: f64 = std.math.pi * (3.0 - @sqrt(5.0));
    const dz: f64 = 2.0 / @as(f64, @floatFromInt(n));

    var longitude: f64 = 0.0;
    var z: f64 = 1.0 - dz / 2.0;

    for (points) |*point| {
        const r = @sqrt(1.0 - z * z);
        point.* = Vec3{
            .x = @cos(longitude) * r,
            .y = @sin(longitude) * r,
            .z = z,
        };
        z -= dz;
        longitude += dlong;
    }

    return points;
}

// Tests
test "generateTestPoints - point count" {
    const allocator = std.testing.allocator;

    const points = try generateTestPoints(allocator, 100);
    defer allocator.free(points);

    try std.testing.expectEqual(@as(usize, 100), points.len);
}

test "generateTestPoints - points on unit sphere" {
    const allocator = std.testing.allocator;
    const tolerance = 1e-10;

    const points = try generateTestPoints(allocator, 100);
    defer allocator.free(points);

    // Verify all points are on the unit sphere (length ≈ 1.0)
    for (points) |point| {
        const len = point.length();
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), len, tolerance);
    }
}

test "generateTestPoints - n=1000" {
    const allocator = std.testing.allocator;
    const tolerance = 1e-10;

    const points = try generateTestPoints(allocator, 1000);
    defer allocator.free(points);

    try std.testing.expectEqual(@as(usize, 1000), points.len);

    // Verify all points are on the unit sphere
    for (points) |point| {
        const len = point.length();
        try std.testing.expectApproxEqAbs(@as(f64, 1.0), len, tolerance);
    }
}

test "generateTestPoints - zero points error" {
    const allocator = std.testing.allocator;

    const result = generateTestPoints(allocator, 0);
    try std.testing.expectError(error.InvalidPointCount, result);
}

test "generateTestPoints - first and last points" {
    const allocator = std.testing.allocator;

    const points = try generateTestPoints(allocator, 100);
    defer allocator.free(points);

    // First point should be near the north pole (z ≈ 1)
    const first = points[0];
    try std.testing.expect(first.z > 0.98);

    // Last point should be near the south pole (z ≈ -1)
    const last = points[points.len - 1];
    try std.testing.expect(last.z < -0.98);
}

const std = @import("std");
const types = @import("types.zig");
const test_points = @import("test_points.zig");

const Allocator = std.mem.Allocator;

/// Supported n_points values for bitmask mode.
pub const supported_n_points = [_]u32{ 64, 128, 256 };
/// Maximum number of u64 words needed for the largest supported n_points.
pub const max_words: usize = (supported_n_points[supported_n_points.len - 1] + 63) / 64;

comptime {
    // Ensure max_words fits in the fixed-size visibility array used by shrake_rupley_bitmask.
    std.debug.assert(max_words <= 4);
}

/// Check if n_points is supported for bitmask mode.
pub fn isSupportedNPoints(n: u32) bool {
    for (supported_n_points) |s| {
        if (n == s) return true;
    }
    return false;
}

/// Generic bitmask lookup table for occlusion-based SASA calculation.
/// Precomputes, for every (direction, angle) pair, which test points on a unit
/// sphere are occluded. At runtime each neighbor's occlusion is applied with a
/// single `visibility &= ~mask` operation per u64 word.
pub fn BitmaskLutGen(comptime T: type) type {
    const Vec = types.Vec3Gen(T);

    return struct {
        const Self = @This();

        /// Number of direction bins per axis (total dirs = dir_resolution²).
        pub const dir_resolution: usize = 24;
        /// Total direction bins.
        pub const n_dirs: usize = dir_resolution * dir_resolution;
        /// Number of cosine-threshold quantization levels.
        pub const angle_bins: usize = 256;

        /// Flat array of precomputed masks: [n_dirs][angle_bins][words].
        masks: []u64,
        /// Direction unit vectors for each bin.
        dir_vectors: []Vec,
        /// Number of u64 words per mask (ceil(n_points / 64)).
        words: usize,
        /// Number of test points.
        n_points: u32,
        /// Bitmask for valid bits in the last word (handles n_points not multiple of 64).
        last_word_mask: u64,
        allocator: Allocator,

        /// Build the lookup table for the given n_points.
        /// Returns error.UnsupportedNPoints if n_points is not 64, 128, or 256.
        pub fn init(allocator: Allocator, n_points_val: u32) !Self {
            if (!isSupportedNPoints(n_points_val)) return error.UnsupportedNPoints;

            const words: usize = (n_points_val + 63) / 64;
            const last_bits = n_points_val % 64;
            const last_word_mask: u64 = if (last_bits == 0)
                std.math.maxInt(u64)
            else
                (@as(u64, 1) << @intCast(last_bits)) - 1;

            // Generate direction vectors using golden spiral
            const dir_vectors = try generateDirections(allocator, n_dirs);
            errdefer allocator.free(dir_vectors);

            // Generate test points on unit sphere
            const TestPointsGen = test_points.generateTestPointsGen(T);
            const tp = try TestPointsGen.generate(allocator, n_points_val);
            defer allocator.free(tp);

            // Allocate masks: n_dirs * angle_bins * words
            const total_masks = n_dirs * angle_bins * words;
            const masks = try allocator.alloc(u64, total_masks);
            errdefer allocator.free(masks);
            @memset(masks, 0);

            // Build all masks
            for (0..n_dirs) |dir_idx| {
                const dir = dir_vectors[dir_idx];
                for (0..angle_bins) |angle_idx| {
                    // Convert angle bin back to cos threshold
                    const cos_threshold = angleBinToCos(angle_idx);
                    const mask_offset = (dir_idx * angle_bins + angle_idx) * words;

                    for (0..n_points_val) |k| {
                        const dot = tp[k].x * dir.x + tp[k].y * dir.y + tp[k].z * dir.z;
                        // Point is occluded if dot >= cos_threshold
                        if (dot >= cos_threshold) {
                            const word_idx = k / 64;
                            const bit_idx: u6 = @intCast(k % 64);
                            masks[mask_offset + word_idx] |= @as(u64, 1) << bit_idx;
                        }
                    }
                }
            }

            return Self{
                .masks = masks,
                .dir_vectors = dir_vectors,
                .words = words,
                .n_points = n_points_val,
                .last_word_mask = last_word_mask,
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.masks);
            self.allocator.free(self.dir_vectors);
        }

        /// Get the occlusion mask for a given (direction, angle) pair.
        pub fn getMask(self: Self, dir_idx: usize, angle_idx: usize) []const u64 {
            const offset = (dir_idx * angle_bins + angle_idx) * self.words;
            return self.masks[offset..][0..self.words];
        }

        /// Find the closest direction bin for a unit vector.
        pub fn dirBinFromUnit(self: Self, x: T, y: T, z: T) usize {
            var best_idx: usize = 0;
            var best_dot: T = -2.0;
            for (self.dir_vectors, 0..) |dv, i| {
                const dot = x * dv.x + y * dv.y + z * dv.z;
                if (dot > best_dot) {
                    best_dot = dot;
                    best_idx = i;
                }
            }
            return best_idx;
        }

        /// Map a cosine value from [-1, 1] to angle bin [0, 255].
        pub fn angleBinFromCos(cos_value: T) u8 {
            // Map [-1, 1] → [0, 255]
            const normalized = (cos_value + 1.0) * 0.5;
            const clamped = @max(@as(T, 0.0), @min(@as(T, 1.0), normalized));
            const bin_f = clamped * 255.0;
            const bin_rounded = @round(bin_f);
            return @intFromFloat(@max(@as(T, 0.0), @min(@as(T, 255.0), bin_rounded)));
        }

        /// Convert an angle bin index back to its cosine threshold value.
        fn angleBinToCos(angle_idx: usize) T {
            return @as(T, @floatFromInt(angle_idx)) / 255.0 * 2.0 - 1.0;
        }

        /// Generate evenly distributed direction vectors using golden spiral.
        fn generateDirections(allocator: Allocator, n: usize) ![]Vec {
            const dirs = try allocator.alloc(Vec, n);
            errdefer allocator.free(dirs);

            const dlong: f64 = std.math.pi * (3.0 - @sqrt(5.0));
            const dz: f64 = 2.0 / @as(f64, @floatFromInt(n));

            var longitude: f64 = 0.0;
            var z: f64 = 1.0 - dz / 2.0;

            for (dirs) |*dir| {
                const r = @sqrt(1.0 - z * z);
                dir.* = Vec{
                    .x = @floatCast(@cos(longitude) * r),
                    .y = @floatCast(@sin(longitude) * r),
                    .z = @floatCast(z),
                };
                z -= dz;
                longitude += dlong;
            }

            return dirs;
        }
    };
}

/// Default f64 bitmask LUT.
pub const BitmaskLut = BitmaskLutGen(f64);
/// Single precision bitmask LUT.
pub const BitmaskLutf32 = BitmaskLutGen(f32);

// =============================================================================
// Tests
// =============================================================================

test "BitmaskLut init and deinit - 64 points" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 64);
    defer lut.deinit();

    try std.testing.expectEqual(@as(usize, 1), lut.words);
    try std.testing.expectEqual(@as(u32, 64), lut.n_points);
    try std.testing.expectEqual(std.math.maxInt(u64), lut.last_word_mask);
}

test "BitmaskLut init and deinit - 128 points" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 128);
    defer lut.deinit();

    try std.testing.expectEqual(@as(usize, 2), lut.words);
    try std.testing.expectEqual(@as(u32, 128), lut.n_points);
    try std.testing.expectEqual(std.math.maxInt(u64), lut.last_word_mask);
}

test "BitmaskLut init and deinit - 256 points" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 256);
    defer lut.deinit();

    try std.testing.expectEqual(@as(usize, 4), lut.words);
    try std.testing.expectEqual(@as(u32, 256), lut.n_points);
    try std.testing.expectEqual(std.math.maxInt(u64), lut.last_word_mask);
}

test "BitmaskLut invalid n_points" {
    const allocator = std.testing.allocator;
    const result = BitmaskLut.init(allocator, 100);
    try std.testing.expectError(error.UnsupportedNPoints, result);
}

test "BitmaskLut angleBinFromCos edges" {
    // cos = -1.0 → bin 0
    try std.testing.expectEqual(@as(u8, 0), BitmaskLut.angleBinFromCos(-1.0));
    // cos = 1.0 → bin 255
    try std.testing.expectEqual(@as(u8, 255), BitmaskLut.angleBinFromCos(1.0));
    // cos = 0.0 → bin ~128
    const mid = BitmaskLut.angleBinFromCos(0.0);
    try std.testing.expect(mid >= 127 and mid <= 128);
}

test "BitmaskLut dirBinFromUnit - different directions" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 64);
    defer lut.deinit();

    // North pole vs south pole should map to different bins
    const north = lut.dirBinFromUnit(0.0, 0.0, 1.0);
    const south = lut.dirBinFromUnit(0.0, 0.0, -1.0);
    const east = lut.dirBinFromUnit(1.0, 0.0, 0.0);

    try std.testing.expect(north != south);
    try std.testing.expect(north != east);
    try std.testing.expect(south != east);
}

test "BitmaskLut mask patterns - north pole occlusion" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 64);
    defer lut.deinit();

    // Direction pointing at north pole, cos_threshold = 0.0 → should occlude ~half
    const dir_idx = lut.dirBinFromUnit(0.0, 0.0, 1.0);
    const angle_idx = BitmaskLut.angleBinFromCos(0.0);
    const mask = lut.getMask(dir_idx, angle_idx);

    const count = @popCount(mask[0]);
    // Should occlude roughly half the 64 points (allow 20-44 range)
    try std.testing.expect(count >= 20 and count <= 44);
}

test "BitmaskLut mask patterns - full occlusion" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 64);
    defer lut.deinit();

    // cos_threshold = -1.0 → all points occluded
    const dir_idx: usize = 0;
    const angle_idx = BitmaskLut.angleBinFromCos(-1.0);
    const mask = lut.getMask(dir_idx, angle_idx);

    const count = @popCount(mask[0]);
    try std.testing.expectEqual(@as(u7, 64), count);
}

test "BitmaskLut mask patterns - no occlusion" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 64);
    defer lut.deinit();

    // cos_threshold = 1.0 → only points exactly at the pole (very few)
    const dir_idx: usize = 0;
    const angle_idx = BitmaskLut.angleBinFromCos(1.0);
    const mask = lut.getMask(dir_idx, angle_idx);

    const count = @popCount(mask[0]);
    // Very few or no points should be exactly at cos = 1.0
    try std.testing.expect(count <= 5);
}

test "BitmaskLut last_word_mask correctness" {
    // 64 points → 1 word, all bits valid
    try std.testing.expectEqual(std.math.maxInt(u64), blk: {
        const allocator = std.testing.allocator;
        var lut = try BitmaskLut.init(allocator, 64);
        defer lut.deinit();
        break :blk lut.last_word_mask;
    });

    // 128 points → 2 words, all bits valid
    try std.testing.expectEqual(std.math.maxInt(u64), blk: {
        const allocator = std.testing.allocator;
        var lut = try BitmaskLut.init(allocator, 128);
        defer lut.deinit();
        break :blk lut.last_word_mask;
    });
}

test "BitmaskLutf32 init and deinit" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLutf32.init(allocator, 128);
    defer lut.deinit();

    try std.testing.expectEqual(@as(usize, 2), lut.words);
    try std.testing.expectEqual(@as(u32, 128), lut.n_points);
}

test "isSupportedNPoints" {
    try std.testing.expect(isSupportedNPoints(64));
    try std.testing.expect(isSupportedNPoints(128));
    try std.testing.expect(isSupportedNPoints(256));
    try std.testing.expect(!isSupportedNPoints(100));
    try std.testing.expect(!isSupportedNPoints(0));
    try std.testing.expect(!isSupportedNPoints(512));
}

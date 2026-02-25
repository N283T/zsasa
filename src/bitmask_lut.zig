const std = @import("std");
const types = @import("types.zig");
const test_points = @import("test_points.zig");

const Allocator = std.mem.Allocator;

/// Maximum n_points supported by bitmask mode (limited by fixed-size visibility array).
pub const max_n_points: u32 = 1024;
/// Maximum number of u64 words needed for the largest supported n_points.
pub const max_words: usize = (max_n_points + 63) / 64;

comptime {
    std.debug.assert(max_words <= 16);
}

/// Check if n_points is valid for bitmask mode (1..256).
pub fn isSupportedNPoints(n: u32) bool {
    return n >= 1 and n <= max_n_points;
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
        /// Returns error.UnsupportedNPoints if n_points is not in 1..256.
        pub fn init(allocator: Allocator, n_points_val: u32) !Self {
            if (!isSupportedNPoints(n_points_val)) return error.UnsupportedNPoints;

            const words: usize = (n_points_val + 63) / 64;
            const last_bits = n_points_val % 64;
            const last_word_mask: u64 = if (last_bits == 0)
                std.math.maxInt(u64)
            else
                (@as(u64, 1) << @intCast(last_bits)) - 1;

            // Generate direction vectors from octahedral grid cell centers
            const dir_vectors = try generateOctaDirections(allocator);
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

        /// Find the direction bin for a unit vector using octahedral encoding (O(1)).
        pub fn dirBinFromUnit(_: Self, x: T, y: T, z: T) usize {
            return octaDirBin(x, y, z);
        }

        /// O(1) octahedral encoding: maps a unit vector to a direction bin index.
        fn octaDirBin(x: T, y: T, z: T) usize {
            // Project onto L1 octahedron
            const inv_l1 = 1.0 / (@abs(x) + @abs(y) + @abs(z));
            var px = x * inv_l1;
            var py = y * inv_l1;
            const pz = z * inv_l1;

            // Fold lower hemisphere onto upper
            if (pz < 0.0) {
                const px_sign: T = if (px >= 0.0) 1.0 else -1.0;
                const py_sign: T = if (py >= 0.0) 1.0 else -1.0;
                const px_abs = @abs(px);
                const py_abs = @abs(py);
                px = (1.0 - py_abs) * px_sign;
                py = (1.0 - px_abs) * py_sign;
            }

            // Map [-1, 1] → [0, resolution)
            const res_f: T = @floatFromInt(dir_resolution);
            const fx = (px + 1.0) * 0.5 * res_f;
            const fy = (py + 1.0) * 0.5 * res_f;
            const ix: usize = @intFromFloat(@min(@max(fx, 0.0), res_f - 1.0));
            const iy: usize = @intFromFloat(@min(@max(fy, 0.0), res_f - 1.0));
            return iy * dir_resolution + ix;
        }

        /// Decode octahedral grid cell center to unit vector.
        fn octaDecode(ix: usize, iy: usize) Vec {
            const res_f: T = @floatFromInt(dir_resolution);
            var u = (@as(T, @floatFromInt(ix)) + 0.5) / res_f * 2.0 - 1.0;
            var v = (@as(T, @floatFromInt(iy)) + 0.5) / res_f * 2.0 - 1.0;
            const zz: T = 1.0 - @abs(u) - @abs(v);

            if (zz < 0.0) {
                const u_sign: T = if (u >= 0.0) 1.0 else -1.0;
                const v_sign: T = if (v >= 0.0) 1.0 else -1.0;
                const u_abs = @abs(u);
                const v_abs = @abs(v);
                u = (1.0 - v_abs) * u_sign;
                v = (1.0 - u_abs) * v_sign;
            }

            const norm = @sqrt(u * u + v * v + zz * zz);
            if (norm <= 0.0) return Vec{ .x = 0.0, .y = 0.0, .z = 1.0 };
            const inv_norm = 1.0 / norm;
            return Vec{ .x = u * inv_norm, .y = v * inv_norm, .z = zz * inv_norm };
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

        /// Generate direction vectors from octahedral grid cell centers.
        /// Each direction corresponds to the center of its octahedral grid cell,
        /// ensuring consistency between dirBinFromUnit lookup and mask construction.
        fn generateOctaDirections(allocator: Allocator) ![]Vec {
            const dirs = try allocator.alloc(Vec, n_dirs);
            errdefer allocator.free(dirs);

            for (0..dir_resolution) |iy| {
                for (0..dir_resolution) |ix| {
                    dirs[iy * dir_resolution + ix] = octaDecode(ix, iy);
                }
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
    try std.testing.expectError(error.UnsupportedNPoints, BitmaskLut.init(allocator, 0));
    try std.testing.expectError(error.UnsupportedNPoints, BitmaskLut.init(allocator, 1025));
    try std.testing.expectError(error.UnsupportedNPoints, BitmaskLut.init(allocator, 2048));
}

test "BitmaskLut arbitrary n_points - 100 points" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 100);
    defer lut.deinit();

    try std.testing.expectEqual(@as(usize, 2), lut.words);
    try std.testing.expectEqual(@as(u32, 100), lut.n_points);
    // 100 % 64 = 36 → lower 36 bits set
    try std.testing.expectEqual((@as(u64, 1) << 36) - 1, lut.last_word_mask);
}

test "BitmaskLut arbitrary n_points - 200 points" {
    const allocator = std.testing.allocator;
    var lut = try BitmaskLut.init(allocator, 200);
    defer lut.deinit();

    try std.testing.expectEqual(@as(usize, 4), lut.words);
    try std.testing.expectEqual(@as(u32, 200), lut.n_points);
    // 200 % 64 = 8 → lower 8 bits set
    try std.testing.expectEqual((@as(u64, 1) << 8) - 1, lut.last_word_mask);
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
    try std.testing.expect(!isSupportedNPoints(0));
    try std.testing.expect(isSupportedNPoints(1));
    try std.testing.expect(isSupportedNPoints(64));
    try std.testing.expect(isSupportedNPoints(100));
    try std.testing.expect(isSupportedNPoints(128));
    try std.testing.expect(isSupportedNPoints(200));
    try std.testing.expect(isSupportedNPoints(256));
    try std.testing.expect(isSupportedNPoints(512));
    try std.testing.expect(isSupportedNPoints(1024));
    try std.testing.expect(!isSupportedNPoints(1025));
    try std.testing.expect(!isSupportedNPoints(2048));
}

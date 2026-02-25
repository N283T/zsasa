const std = @import("std");
const types = @import("types.zig");
const test_points = @import("test_points.zig");
const neighbor_list = @import("neighbor_list.zig");
const bitmask_lut = @import("bitmask_lut.zig");
const thread_pool = @import("thread_pool.zig");

const Vec3 = types.Vec3;
const Vec3Gen = types.Vec3Gen;
const NeighborList = neighbor_list.NeighborList;
const NeighborListGen = neighbor_list.NeighborListGen;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const SasaResultGen = types.SasaResultGen;
const Config = types.Config;
const ConfigGen = types.ConfigGen;
const BitmaskLut = bitmask_lut.BitmaskLut;
const BitmaskLutGen = bitmask_lut.BitmaskLutGen;
const Allocator = std.mem.Allocator;

/// Calculate SASA for a single atom using bitmask-based occlusion.
///
/// Instead of checking every test point against every neighbor, this approach:
/// 1. Initializes a visibility bitmask with all test points marked as exposed
/// 2. For each neighbor, looks up a precomputed occlusion mask and applies it
/// 3. Counts remaining exposed points via popcount
///
/// This eliminates the O(n_points) inner loop, replacing it with O(words) bitwise ops.
fn atomSasaBitmask(
    atom_idx: usize,
    positions: []const Vec3,
    radii_with_probe_sq: []const f64,
    atom_radius_probe: f64,
    neighbors: []const u32,
    lut: *const BitmaskLut,
) f64 {
    const words = lut.words;
    comptime std.debug.assert(bitmask_lut.max_words <= 4);

    // Early exit: no neighbors means full surface is exposed
    if (neighbors.len == 0) {
        return 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
    }

    const atom_pos = positions[atom_idx];
    const r_i = atom_radius_probe;
    const r_i_sq = r_i * r_i;

    // Initialize visibility: all points exposed
    var visibility: [4]u64 = .{ std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64) };
    // Mask off invalid bits in the last word
    visibility[words - 1] &= lut.last_word_mask;

    for (neighbors) |j| {
        // Compute direction from atom i to neighbor j
        const dx = positions[j].x - atom_pos.x;
        const dy = positions[j].y - atom_pos.y;
        const dz = positions[j].z - atom_pos.z;
        const dist_sq = dx * dx + dy * dy + dz * dz;
        const dist = @sqrt(dist_sq);

        if (dist < 1e-10) continue;

        // Compute cosine threshold for the occlusion cone (law of cosines)
        // cos_threshold = (r_i² + dist² - r_j²) / (2 * r_i * dist)
        // A test point at angle θ from i→j direction is buried if cos(θ) > cos_threshold
        const r_j_sq = radii_with_probe_sq[j];
        const cos_threshold = (r_i_sq + dist_sq - r_j_sq) / (2.0 * r_i * dist);

        // cos_threshold >= 1.0: neighbor too far to occlude anything
        if (cos_threshold >= 1.0) continue;

        // cos_threshold <= -1.0: neighbor completely covers this atom
        if (cos_threshold <= -1.0) return 0.0;

        // Normalize direction
        const inv_dist = 1.0 / dist;
        const nx = dx * inv_dist;
        const ny = dy * inv_dist;
        const nz = dz * inv_dist;

        // Look up precomputed mask
        const dir_idx = lut.dirBinFromUnit(nx, ny, nz);
        const angle_idx = BitmaskLut.angleBinFromCos(cos_threshold);
        const mask = lut.getMask(dir_idx, angle_idx);

        // Apply occlusion: clear bits for occluded points
        for (0..words) |w| {
            visibility[w] &= ~mask[w];
        }

        // Early exit: all points buried
        var any_visible = false;
        for (0..words) |w| {
            if (visibility[w] != 0) {
                any_visible = true;
                break;
            }
        }
        if (!any_visible) return 0.0;
    }

    // Count exposed points
    var exposed: u32 = 0;
    for (0..words) |w| {
        exposed += @popCount(visibility[w]);
    }

    // Calculate SASA: 4π * r² * (exposed / n_points)
    const surface_area = 4.0 * std.math.pi * r_i * r_i;
    const exposed_fraction = @as(f64, @floatFromInt(exposed)) / @as(f64, @floatFromInt(lut.n_points));
    return surface_area * exposed_fraction;
}

/// Calculate SASA for all atoms using bitmask-optimized Shrake-Rupley (single-threaded).
pub fn calculateSasa(
    allocator: Allocator,
    input: AtomInput,
    config: Config,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) return error.NoAtoms;
    if (!bitmask_lut.isSupportedNPoints(config.n_points)) return error.UnsupportedNPoints;

    // Build bitmask LUT
    var lut = try BitmaskLut.init(allocator, config.n_points);
    defer lut.deinit();

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

    // Build neighbor list
    var neighbor_list_data = try NeighborList.init(allocator, positions, input.r, config.probe_radius);
    defer neighbor_list_data.deinit();

    // Pre-compute (r[j] + probe_radius)² for each atom
    const radii_with_probe_sq = try allocator.alloc(f64, n_atoms);
    defer allocator.free(radii_with_probe_sq);

    for (0..n_atoms) |i| {
        const r_probe = input.r[i] + config.probe_radius;
        radii_with_probe_sq[i] = r_probe * r_probe;
    }

    // Allocate result arrays
    const atom_areas = try allocator.alloc(f64, n_atoms);
    errdefer allocator.free(atom_areas);

    // Calculate SASA for each atom
    var total_area: f64 = 0.0;
    for (0..n_atoms) |i| {
        const atom_radius_probe = input.r[i] + config.probe_radius;
        const neighbors = neighbor_list_data.getNeighbors(i);
        const area = atomSasaBitmask(
            i,
            positions,
            radii_with_probe_sq,
            atom_radius_probe,
            neighbors,
            &lut,
        );
        atom_areas[i] = area;
        total_area += area;
    }

    return SasaResult{
        .total_area = total_area,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}

/// Context for parallel bitmask SASA calculation workers.
const ParallelContext = struct {
    positions: []const Vec3,
    radii: []const f64,
    radii_with_probe_sq: []const f64,
    neighbor_list_data: *const NeighborList,
    probe_radius: f64,
    atom_areas: []f64,
    lut: *const BitmaskLut,
};

/// Worker function for parallel bitmask SASA calculation.
fn parallelSasaWorker(ctx: ParallelContext, chunk_start: usize, chunk_end: usize) f64 {
    var chunk_total: f64 = 0.0;

    for (chunk_start..chunk_end) |i| {
        const atom_radius_probe = ctx.radii[i] + ctx.probe_radius;
        const neighbors = ctx.neighbor_list_data.getNeighbors(i);
        const area = atomSasaBitmask(
            i,
            ctx.positions,
            ctx.radii_with_probe_sq,
            atom_radius_probe,
            neighbors,
            ctx.lut,
        );
        ctx.atom_areas[i] = area;
        chunk_total += area;
    }

    return chunk_total;
}

/// Reduce function to sum all chunk totals.
fn sumReducer(results: []const f64) f64 {
    var total: f64 = 0.0;
    for (results) |r| {
        total += r;
    }
    return total;
}

/// Calculate SASA for all atoms using bitmask-optimized Shrake-Rupley (parallel).
pub fn calculateSasaParallel(
    allocator: Allocator,
    input: AtomInput,
    config: Config,
    n_threads: usize,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) return error.NoAtoms;
    if (!bitmask_lut.isSupportedNPoints(config.n_points)) return error.UnsupportedNPoints;

    const actual_threads = if (n_threads == 0)
        try std.Thread.getCpuCount()
    else
        n_threads;

    // Build bitmask LUT
    var lut = try BitmaskLut.init(allocator, config.n_points);
    defer lut.deinit();

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

    // Build neighbor list
    var neighbor_list_data = try NeighborList.init(allocator, positions, input.r, config.probe_radius);
    defer neighbor_list_data.deinit();

    // Pre-compute (r[j] + probe_radius)² for each atom
    const radii_with_probe_sq = try allocator.alloc(f64, n_atoms);
    defer allocator.free(radii_with_probe_sq);

    for (0..n_atoms) |i| {
        const r_probe = input.r[i] + config.probe_radius;
        radii_with_probe_sq[i] = r_probe * r_probe;
    }

    // Allocate result arrays
    const atom_areas = try allocator.alloc(f64, n_atoms);
    errdefer allocator.free(atom_areas);

    // Create parallel context
    const ctx = ParallelContext{
        .positions = positions,
        .radii = input.r,
        .radii_with_probe_sq = radii_with_probe_sq,
        .neighbor_list_data = &neighbor_list_data,
        .probe_radius = config.probe_radius,
        .atom_areas = atom_areas,
        .lut = &lut,
    };

    const chunk_size = @max(64, n_atoms / (actual_threads * 4));

    const total_area = try thread_pool.parallelFor(
        ParallelContext,
        f64,
        allocator,
        actual_threads,
        parallelSasaWorker,
        ctx,
        n_atoms,
        chunk_size,
        sumReducer,
    );

    return SasaResult{
        .total_area = total_area,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}

// =============================================================================
// Generic f32/f64 implementations
// =============================================================================

/// Generic bitmask Shrake-Rupley implementation supporting both f32 and f64.
pub fn ShrakeRupleyBitmaskGen(comptime T: type) type {
    const Vec = types.Vec3Gen(T);
    const Result = SasaResultGen(T);
    const Cfg = ConfigGen(T);
    const NList = NeighborListGen(T);
    const Lut = BitmaskLutGen(T);

    return struct {
        const Self = @This();

        /// Calculate SASA for a single atom using bitmask occlusion.
        fn atomSasaBitmask(
            atom_idx: usize,
            positions: []const Vec,
            radii_with_probe_sq: []const T,
            atom_radius_probe: T,
            neighbors: []const u32,
            lut: *const Lut,
        ) T {
            const words = lut.words;
            comptime std.debug.assert(bitmask_lut.max_words <= 4);

            if (neighbors.len == 0) {
                return 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
            }

            const atom_pos = positions[atom_idx];
            const r_i = atom_radius_probe;
            const r_i_sq = r_i * r_i;

            var visibility: [4]u64 = .{ std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64) };
            visibility[words - 1] &= lut.last_word_mask;

            for (neighbors) |j| {
                const dx = positions[j].x - atom_pos.x;
                const dy = positions[j].y - atom_pos.y;
                const dz = positions[j].z - atom_pos.z;
                const dist_sq = dx * dx + dy * dy + dz * dz;
                const dist = @sqrt(dist_sq);

                if (dist < 1e-10) continue;

                const r_j_sq = radii_with_probe_sq[j];
                const cos_threshold = (r_i_sq + dist_sq - r_j_sq) / (2.0 * r_i * dist);

                if (cos_threshold >= 1.0) continue;
                if (cos_threshold <= -1.0) return 0.0;

                const inv_dist: T = 1.0 / dist;
                const nx = dx * inv_dist;
                const ny = dy * inv_dist;
                const nz = dz * inv_dist;

                const dir_idx = lut.dirBinFromUnit(nx, ny, nz);
                const angle_idx = Lut.angleBinFromCos(cos_threshold);
                const mask = lut.getMask(dir_idx, angle_idx);

                for (0..words) |w| {
                    visibility[w] &= ~mask[w];
                }

                var any_visible = false;
                for (0..words) |w| {
                    if (visibility[w] != 0) {
                        any_visible = true;
                        break;
                    }
                }
                if (!any_visible) return 0.0;
            }

            var exposed: u32 = 0;
            for (0..words) |w| {
                exposed += @popCount(visibility[w]);
            }

            const surface_area = 4.0 * std.math.pi * r_i * r_i;
            const exposed_fraction = @as(T, @floatFromInt(exposed)) / @as(T, @floatFromInt(lut.n_points));
            return surface_area * exposed_fraction;
        }

        pub const ParallelContext = struct {
            positions: []const Vec,
            radii: []const T,
            radii_with_probe_sq: []const T,
            neighbor_list_data: *const NList,
            probe_radius: T,
            atom_areas: []T,
            lut: *const Lut,
        };

        fn parallelSasaWorker(ctx: Self.ParallelContext, chunk_start: usize, chunk_end: usize) T {
            var chunk_total: T = 0.0;

            for (chunk_start..chunk_end) |i| {
                const atom_radius_probe = ctx.radii[i] + ctx.probe_radius;
                const neighbors = ctx.neighbor_list_data.getNeighbors(i);
                const area = Self.atomSasaBitmask(
                    i,
                    ctx.positions,
                    ctx.radii_with_probe_sq,
                    atom_radius_probe,
                    neighbors,
                    ctx.lut,
                );
                ctx.atom_areas[i] = area;
                chunk_total += area;
            }

            return chunk_total;
        }

        fn sumReducer(results: []const T) T {
            var total: T = 0.0;
            for (results) |r| {
                total += r;
            }
            return total;
        }

        /// Calculate SASA (single-threaded, generic precision).
        pub fn calculateSasa(
            allocator: Allocator,
            input: AtomInput,
            config: Cfg,
        ) !Result {
            const n_atoms = input.atomCount();
            if (n_atoms == 0) return error.NoAtoms;
            if (!bitmask_lut.isSupportedNPoints(config.n_points)) return error.UnsupportedNPoints;

            var lut = try Lut.init(allocator, config.n_points);
            defer lut.deinit();

            const positions = try allocator.alloc(Vec, n_atoms);
            defer allocator.free(positions);

            for (0..n_atoms) |i| {
                positions[i] = Vec{
                    .x = @floatCast(input.x[i]),
                    .y = @floatCast(input.y[i]),
                    .z = @floatCast(input.z[i]),
                };
            }

            const radii = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii);
            for (0..n_atoms) |i| {
                radii[i] = @floatCast(input.r[i]);
            }

            var neighbor_list_data = try NList.init(allocator, positions, radii, config.probe_radius);
            defer neighbor_list_data.deinit();

            const radii_with_probe_sq = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii_with_probe_sq);

            for (0..n_atoms) |i| {
                const r_probe = radii[i] + config.probe_radius;
                radii_with_probe_sq[i] = r_probe * r_probe;
            }

            const atom_areas = try allocator.alloc(T, n_atoms);
            errdefer allocator.free(atom_areas);

            var total_area: T = 0.0;
            for (0..n_atoms) |i| {
                const atom_radius_probe = radii[i] + config.probe_radius;
                const neighbors = neighbor_list_data.getNeighbors(i);
                const area = Self.atomSasaBitmask(
                    i,
                    positions,
                    radii_with_probe_sq,
                    atom_radius_probe,
                    neighbors,
                    &lut,
                );
                atom_areas[i] = area;
                total_area += area;
            }

            return Result{
                .total_area = total_area,
                .atom_areas = atom_areas,
                .allocator = allocator,
            };
        }

        /// Calculate SASA (parallel, generic precision).
        pub fn calculateSasaParallel(
            allocator: Allocator,
            input: AtomInput,
            config: Cfg,
            n_threads: usize,
        ) !Result {
            const n_atoms = input.atomCount();
            if (n_atoms == 0) return error.NoAtoms;
            if (!bitmask_lut.isSupportedNPoints(config.n_points)) return error.UnsupportedNPoints;

            const actual_threads = if (n_threads == 0)
                try std.Thread.getCpuCount()
            else
                n_threads;

            var lut = try Lut.init(allocator, config.n_points);
            defer lut.deinit();

            const positions = try allocator.alloc(Vec, n_atoms);
            defer allocator.free(positions);

            for (0..n_atoms) |i| {
                positions[i] = Vec{
                    .x = @floatCast(input.x[i]),
                    .y = @floatCast(input.y[i]),
                    .z = @floatCast(input.z[i]),
                };
            }

            const radii = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii);
            for (0..n_atoms) |i| {
                radii[i] = @floatCast(input.r[i]);
            }

            var neighbor_list_data = try NList.init(allocator, positions, radii, config.probe_radius);
            defer neighbor_list_data.deinit();

            const radii_with_probe_sq = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii_with_probe_sq);

            for (0..n_atoms) |i| {
                const r_probe = radii[i] + config.probe_radius;
                radii_with_probe_sq[i] = r_probe * r_probe;
            }

            const atom_areas = try allocator.alloc(T, n_atoms);
            errdefer allocator.free(atom_areas);

            const ctx = Self.ParallelContext{
                .positions = positions,
                .radii = radii,
                .radii_with_probe_sq = radii_with_probe_sq,
                .neighbor_list_data = &neighbor_list_data,
                .probe_radius = config.probe_radius,
                .atom_areas = atom_areas,
                .lut = &lut,
            };

            const chunk_size = @max(64, n_atoms / (actual_threads * 4));

            const total_area = try thread_pool.parallelFor(
                Self.ParallelContext,
                T,
                allocator,
                actual_threads,
                Self.parallelSasaWorker,
                ctx,
                n_atoms,
                chunk_size,
                Self.sumReducer,
            );

            return Result{
                .total_area = total_area,
                .atom_areas = atom_areas,
                .allocator = allocator,
            };
        }
    };
}

/// f32 precision bitmask Shrake-Rupley implementation.
pub const ShrakeRupleyBitmaskf32 = ShrakeRupleyBitmaskGen(f32);

/// Convenience function: Calculate bitmask SASA using f32 precision (single-threaded).
pub fn calculateSasaf32(allocator: Allocator, input: AtomInput, config: ConfigGen(f32)) !SasaResultGen(f32) {
    return ShrakeRupleyBitmaskf32.calculateSasa(allocator, input, config);
}

/// Convenience function: Calculate bitmask SASA using f32 precision (parallel).
pub fn calculateSasaParallelf32(allocator: Allocator, input: AtomInput, config: ConfigGen(f32), n_threads: usize) !SasaResultGen(f32) {
    return ShrakeRupleyBitmaskf32.calculateSasaParallel(allocator, input, config, n_threads);
}

// =============================================================================
// Tests
// =============================================================================

test "bitmask calculateSasa - single isolated atom" {
    const allocator = std.testing.allocator;

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
        .n_points = 128,
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

test "bitmask calculateSasa - two far apart atoms" {
    const allocator = std.testing.allocator;

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
        .n_points = 128,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    const expected_radius = r[0] + config.probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;
    const expected_total = expected_area * 2.0;

    try std.testing.expectEqual(@as(usize, 2), result.atom_areas.len);
    try std.testing.expectApproxEqRel(expected_total, result.total_area, 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[0], 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[1], 0.01);
}

test "bitmask calculateSasa - two overlapping atoms" {
    const allocator = std.testing.allocator;

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
        .n_points = 128,
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
    try std.testing.expect(result.atom_areas[0] < full_area);
    try std.testing.expect(result.atom_areas[1] < full_area);
}

test "bitmask calculateSasa - unsupported n_points" {
    const allocator = std.testing.allocator;

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
        .n_points = 100,
        .probe_radius = 1.4,
    };

    const result = calculateSasa(allocator, input, config);
    try std.testing.expectError(error.UnsupportedNPoints, result);
}

test "bitmask vs standard SR - equivalence within tolerance" {
    const allocator = std.testing.allocator;
    const shrake_rupley = @import("shrake_rupley.zig");

    // Set up a 4-atom system
    const x = try allocator.alloc(f64, 4);
    const y = try allocator.alloc(f64, 4);
    const z = try allocator.alloc(f64, 4);
    const r = try allocator.alloc(f64, 4);

    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.5;
    x[1] = 2.5;
    y[1] = 0.0;
    z[1] = 0.0;
    r[1] = 1.0;
    x[2] = 0.0;
    y[2] = 2.5;
    z[2] = 0.0;
    r[2] = 1.2;
    x[3] = 1.5;
    y[3] = 1.5;
    z[3] = 1.5;
    r[3] = 1.0;

    // Need two separate inputs since deinit will free the same arrays
    const x2 = try allocator.alloc(f64, 4);
    const y2 = try allocator.alloc(f64, 4);
    const z2 = try allocator.alloc(f64, 4);
    const r2 = try allocator.alloc(f64, 4);
    @memcpy(x2, x);
    @memcpy(y2, y);
    @memcpy(z2, z);
    @memcpy(r2, r);

    var input_bitmask = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input_bitmask.deinit();

    var input_standard = AtomInput{
        .x = x2,
        .y = y2,
        .z = z2,
        .r = r2,
        .allocator = allocator,
    };
    defer input_standard.deinit();

    const config = Config{
        .n_points = 128,
        .probe_radius = 1.4,
    };

    var result_bitmask = try calculateSasa(allocator, input_bitmask, config);
    defer result_bitmask.deinit();

    var result_standard = try shrake_rupley.calculateSasa(allocator, input_standard, config);
    defer result_standard.deinit();

    // Total areas should be within 5% of each other (bitmask is approximate)
    try std.testing.expectApproxEqRel(result_standard.total_area, result_bitmask.total_area, 0.05);

    // Per-atom areas should be within 10% tolerance
    for (0..4) |i| {
        if (result_standard.atom_areas[i] > 0.1) {
            try std.testing.expectApproxEqRel(
                result_standard.atom_areas[i],
                result_bitmask.atom_areas[i],
                0.10,
            );
        }
    }
}

test "bitmask calculateSasa - NoAtoms error" {
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

    const config = Config{
        .n_points = 128,
        .probe_radius = 1.4,
    };

    const result = calculateSasa(allocator, input, config);
    try std.testing.expectError(error.NoAtoms, result);
}

test "bitmask calculateSasa - completely buried atom" {
    const allocator = std.testing.allocator;

    // Place a small atom surrounded by large overlapping neighbors
    const n = 7;
    const x = try allocator.alloc(f64, n);
    const y = try allocator.alloc(f64, n);
    const z = try allocator.alloc(f64, n);
    const r = try allocator.alloc(f64, n);

    // Central small atom
    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 0.5;

    // 6 large surrounding atoms at ±1 along each axis
    const offsets = [_][3]f64{ .{ 1.0, 0, 0 }, .{ -1.0, 0, 0 }, .{ 0, 1.0, 0 }, .{ 0, -1.0, 0 }, .{ 0, 0, 1.0 }, .{ 0, 0, -1.0 } };
    for (offsets, 0..) |off, idx| {
        x[idx + 1] = off[0];
        y[idx + 1] = off[1];
        z[idx + 1] = off[2];
        r[idx + 1] = 2.0;
    }

    var input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input.deinit();

    const config = Config{
        .n_points = 128,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    // Central atom should be nearly or completely buried
    try std.testing.expect(result.atom_areas[0] < 1.0);
}

test "bitmask calculateSasa - n_points 256" {
    const allocator = std.testing.allocator;

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
        .n_points = 256,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    const expected_radius = r[0] + config.probe_radius;
    const full_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    // Both should be partially buried
    try std.testing.expect(result.atom_areas[0] < full_area);
    try std.testing.expect(result.atom_areas[0] > 0.0);
    try std.testing.expect(result.total_area > 0.0);
    try std.testing.expect(result.total_area < full_area * 2.0);
}

test "bitmask calculateSasaf32 - single atom" {
    const allocator = std.testing.allocator;

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

    const config = types.Configf32{
        .n_points = 128,
        .probe_radius = 1.4,
    };

    var result = try calculateSasaf32(allocator, input, config);
    defer result.deinit();

    const expected_radius: f32 = 1.0 + 1.4;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    try std.testing.expectEqual(@as(usize, 1), result.atom_areas.len);
    try std.testing.expectApproxEqRel(expected_area, result.total_area, 0.02);
}

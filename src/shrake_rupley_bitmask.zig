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
/// O(1) octahedral direction bin (standalone, no struct needed).
fn octaDirBinScalar(comptime T: type, x: T, y: T, z: T, resolution: usize) usize {
    const inv_l1 = 1.0 / (@abs(x) + @abs(y) + @abs(z));
    var px = x * inv_l1;
    var py = y * inv_l1;
    const pz = z * inv_l1;

    if (pz < 0.0) {
        const px_sign: T = if (px >= 0.0) 1.0 else -1.0;
        const py_sign: T = if (py >= 0.0) 1.0 else -1.0;
        const px_abs = @abs(px);
        const py_abs = @abs(py);
        px = (1.0 - py_abs) * px_sign;
        py = (1.0 - px_abs) * py_sign;
    }

    const res_f: T = @floatFromInt(resolution);
    const fx = (px + 1.0) * 0.5 * res_f;
    const fy = (py + 1.0) * 0.5 * res_f;
    const ix: usize = @intFromFloat(@min(@max(fx, 0.0), res_f - 1.0));
    const iy: usize = @intFromFloat(@min(@max(fy, 0.0), res_f - 1.0));
    return iy * resolution + ix;
}

/// O(1) angle bin from cosine (standalone, no struct needed).
fn angleBinFromCosScalar(comptime T: type, cos_value: T) usize {
    const normalized = (cos_value + 1.0) * 0.5;
    const clamped = @max(@as(T, 0.0), @min(@as(T, 1.0), normalized));
    const bin_f = clamped * 255.0;
    const bin_rounded = @round(bin_f);
    return @intFromFloat(@max(@as(T, 0.0), @min(@as(T, 255.0), bin_rounded)));
}

/// Branchless octahedral direction bin for a batch of 4 direction vectors.
/// Uses @select instead of if/else for hemisphere fold — maps to blendv (x86) or bsl (NEON).
fn octaDirBinBranchless(comptime T: type, nx: @Vector(4, T), ny: @Vector(4, T), nz: @Vector(4, T), resolution: usize) [4]usize {
    const V4 = @Vector(4, T);
    const ones: V4 = @splat(1.0);
    const zeros: V4 = @splat(0.0);
    const neg_ones: V4 = @splat(-1.0);
    const half: V4 = @splat(0.5);

    const abs_x = @abs(nx);
    const abs_y = @abs(ny);
    const abs_z = @abs(nz);
    const inv_l1 = ones / (abs_x + abs_y + abs_z);
    const px_orig = nx * inv_l1;
    const py_orig = ny * inv_l1;
    const pz = nz * inv_l1;

    // Branchless hemisphere fold for pz < 0
    const neg_z = pz < zeros;
    const px_sign: V4 = @select(T, px_orig >= zeros, ones, neg_ones);
    const py_sign: V4 = @select(T, py_orig >= zeros, ones, neg_ones);
    const folded_px = (ones - @abs(py_orig)) * px_sign;
    const folded_py = (ones - @abs(px_orig)) * py_sign;
    const px = @select(T, neg_z, folded_px, px_orig);
    const py = @select(T, neg_z, folded_py, py_orig);

    // Grid coordinates to bin indices
    const res_f: V4 = @splat(@as(T, @floatFromInt(resolution)));
    const res_m1 = res_f - ones;
    const clamped_fx = @min(@max((px + ones) * half * res_f, zeros), res_m1);
    const clamped_fy = @min(@max((py + ones) * half * res_f, zeros), res_m1);

    var result: [4]usize = undefined;
    inline for (0..4) |i| {
        result[i] = @as(usize, @intFromFloat(clamped_fx[i])) + @as(usize, @intFromFloat(clamped_fy[i])) * resolution;
    }
    return result;
}

/// Batch angle bin computation for 4 cosine values.
fn angleBinFromCosBatch(comptime T: type, cos_vals: @Vector(4, T)) [4]usize {
    const V4 = @Vector(4, T);
    const ones: V4 = @splat(1.0);
    const zeros: V4 = @splat(0.0);
    const max_bin: V4 = @splat(255.0);

    const clamped = @max(zeros, @min(ones, (cos_vals + ones) * @as(V4, @splat(@as(T, 0.5)))));
    const rounded = @round(clamped * max_bin);
    const final_bins = @max(zeros, @min(max_bin, rounded));

    var result: [4]usize = undefined;
    inline for (0..4) |i| {
        result[i] = @intFromFloat(final_bins[i]);
    }
    return result;
}

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

    if (neighbors.len == 0) {
        return 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
    }

    const atom_pos = positions[atom_idx];
    const r_i = atom_radius_probe;
    const r_i_sq = r_i * r_i;
    const two_r_i = 2.0 * r_i;
    const dir_resolution = BitmaskLut.dir_resolution;
    const angle_bins_c = BitmaskLut.angle_bins;

    var visibility: [4]u64 = .{ std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64) };
    visibility[words - 1] &= lut.last_word_mask;

    // SIMD 4-neighbor batching
    const batch_size = 4;
    const V4 = @Vector(batch_size, f64);
    const cx_splat: V4 = @splat(atom_pos.x);
    const cy_splat: V4 = @splat(atom_pos.y);
    const cz_splat: V4 = @splat(atom_pos.z);
    const r_i_sq_splat: V4 = @splat(r_i_sq);
    const inv_two_r_splat: V4 = @splat(@as(f64, 1.0) / two_r_i);
    const ones_v: V4 = @splat(@as(f64, 1.0));
    const neg_ones_v: V4 = @splat(@as(f64, -1.0));

    var idx: usize = 0;

    while (idx + batch_size <= neighbors.len) : (idx += batch_size) {
        var nbr_x: [batch_size]f64 = undefined;
        var nbr_y: [batch_size]f64 = undefined;
        var nbr_z: [batch_size]f64 = undefined;
        var nbr_rsq: [batch_size]f64 = undefined;
        inline for (0..batch_size) |b| {
            const j = neighbors[idx + b];
            nbr_x[b] = positions[j].x;
            nbr_y[b] = positions[j].y;
            nbr_z[b] = positions[j].z;
            nbr_rsq[b] = radii_with_probe_sq[j];
        }

        const vx: V4 = @as(V4, nbr_x) - cx_splat;
        const vy: V4 = @as(V4, nbr_y) - cy_splat;
        const vz: V4 = @as(V4, nbr_z) - cz_splat;
        const dist_sq = vx * vx + vy * vy + vz * vz;
        const inv_dist = ones_v / @sqrt(dist_sq);
        const cos_threshold = (r_i_sq_splat + dist_sq - @as(V4, nbr_rsq)) * inv_two_r_splat * inv_dist;

        if (@reduce(.Or, cos_threshold <= neg_ones_v)) return 0.0;

        const valid = cos_threshold < ones_v;
        if (@reduce(.Or, valid)) {
            const nx = vx * inv_dist;
            const ny = vy * inv_dist;
            const nz = vz * inv_dist;
            const dir_indices = octaDirBinBranchless(f64, nx, ny, nz, dir_resolution);
            const angle_indices = angleBinFromCosBatch(f64, cos_threshold);

            var combined: [4]u64 = .{ 0, 0, 0, 0 };
            inline for (0..batch_size) |b| {
                if (valid[b]) {
                    const offset = (dir_indices[b] * angle_bins_c + angle_indices[b]) * words;
                    for (0..words) |w| {
                        combined[w] |= lut.masks[offset + w];
                    }
                }
            }

            var any_visible: u64 = 0;
            for (0..words) |w| {
                visibility[w] &= ~combined[w];
                any_visible |= visibility[w];
            }
            if (any_visible == 0) return 0.0;
        }
    }

    // Scalar tail for remaining < 4 neighbors
    while (idx < neighbors.len) : (idx += 1) {
        const j = neighbors[idx];
        const dx = positions[j].x - atom_pos.x;
        const dy = positions[j].y - atom_pos.y;
        const dz = positions[j].z - atom_pos.z;
        const dist_sq = dx * dx + dy * dy + dz * dz;
        const dist = @sqrt(dist_sq);

        if (dist < 1e-10) continue;

        const r_j_sq = radii_with_probe_sq[j];
        const cos_threshold = (r_i_sq + dist_sq - r_j_sq) / (two_r_i * dist);

        if (cos_threshold >= 1.0) continue;
        if (cos_threshold <= -1.0) return 0.0;

        const inv_dist = 1.0 / dist;
        const dir_idx = octaDirBinScalar(f64, dx * inv_dist, dy * inv_dist, dz * inv_dist, dir_resolution);
        const angle_idx = angleBinFromCosScalar(f64, cos_threshold);
        const mask_offset = (dir_idx * angle_bins_c + angle_idx) * words;

        var combined: u64 = 0;
        for (0..words) |w| {
            visibility[w] &= ~lut.masks[mask_offset + w];
            combined |= visibility[w];
        }
        if (combined == 0) return 0.0;
    }

    var exposed: u32 = 0;
    for (0..words) |w| {
        exposed += @popCount(visibility[w]);
    }

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

            const r_i = atom_radius_probe;
            const r_i_sq = r_i * r_i;
            const two_r_i = 2.0 * r_i;
            const atom_pos = positions[atom_idx];
            const dir_resolution = Lut.dir_resolution;
            const angle_bins_val = Lut.angle_bins;

            var visibility: [4]u64 = .{ std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64), std.math.maxInt(u64) };
            visibility[words - 1] &= lut.last_word_mask;

            // SIMD 4-neighbor batching
            const batch_size = 4;
            const V4 = @Vector(batch_size, T);
            const cx_splat: V4 = @splat(atom_pos.x);
            const cy_splat: V4 = @splat(atom_pos.y);
            const cz_splat: V4 = @splat(atom_pos.z);
            const r_i_sq_splat: V4 = @splat(r_i_sq);
            const inv_two_r_splat: V4 = @splat(@as(T, 1.0) / two_r_i);
            const ones_v: V4 = @splat(@as(T, 1.0));
            const neg_ones_v: V4 = @splat(@as(T, -1.0));

            var idx: usize = 0;

            while (idx + batch_size <= neighbors.len) : (idx += batch_size) {
                var nbr_x: [batch_size]T = undefined;
                var nbr_y: [batch_size]T = undefined;
                var nbr_z: [batch_size]T = undefined;
                var nbr_rsq: [batch_size]T = undefined;
                inline for (0..batch_size) |b| {
                    const j = neighbors[idx + b];
                    nbr_x[b] = positions[j].x;
                    nbr_y[b] = positions[j].y;
                    nbr_z[b] = positions[j].z;
                    nbr_rsq[b] = radii_with_probe_sq[j];
                }

                const vx: V4 = @as(V4, nbr_x) - cx_splat;
                const vy: V4 = @as(V4, nbr_y) - cy_splat;
                const vz: V4 = @as(V4, nbr_z) - cz_splat;
                const dist_sq = vx * vx + vy * vy + vz * vz;
                const inv_dist = ones_v / @sqrt(dist_sq);
                const cos_threshold = (r_i_sq_splat + dist_sq - @as(V4, nbr_rsq)) * inv_two_r_splat * inv_dist;

                if (@reduce(.Or, cos_threshold <= neg_ones_v)) return 0.0;

                const valid = cos_threshold < ones_v;
                if (@reduce(.Or, valid)) {
                    const nx = vx * inv_dist;
                    const ny = vy * inv_dist;
                    const nz = vz * inv_dist;
                    const dir_indices = octaDirBinBranchless(T, nx, ny, nz, dir_resolution);
                    const angle_indices = angleBinFromCosBatch(T, cos_threshold);

                    var combined: [4]u64 = .{ 0, 0, 0, 0 };
                    inline for (0..batch_size) |b| {
                        if (valid[b]) {
                            const offset = (dir_indices[b] * angle_bins_val + angle_indices[b]) * words;
                            for (0..words) |w| {
                                combined[w] |= lut.masks[offset + w];
                            }
                        }
                    }

                    var any_visible: u64 = 0;
                    for (0..words) |w| {
                        visibility[w] &= ~combined[w];
                        any_visible |= visibility[w];
                    }
                    if (any_visible == 0) return 0.0;
                }
            }

            // Scalar tail for remaining < 4 neighbors
            while (idx < neighbors.len) : (idx += 1) {
                const j = neighbors[idx];
                const dx = positions[j].x - atom_pos.x;
                const dy = positions[j].y - atom_pos.y;
                const dz = positions[j].z - atom_pos.z;
                const dist_sq = dx * dx + dy * dy + dz * dz;
                const dist = @sqrt(dist_sq);

                if (dist < 1e-10) continue;

                const r_j_sq = radii_with_probe_sq[j];
                const cos_threshold = (r_i_sq + dist_sq - r_j_sq) / (two_r_i * dist);

                if (cos_threshold >= 1.0) continue;
                if (cos_threshold <= -1.0) return 0.0;

                const inv_dist: T = 1.0 / dist;
                const dir_idx = octaDirBinScalar(T, dx * inv_dist, dy * inv_dist, dz * inv_dist, dir_resolution);
                const angle_idx = angleBinFromCosScalar(T, cos_threshold);
                const mask_offset = (dir_idx * angle_bins_val + angle_idx) * words;

                var combined: u64 = 0;
                for (0..words) |w| {
                    visibility[w] &= ~lut.masks[mask_offset + w];
                    combined |= visibility[w];
                }
                if (combined == 0) return 0.0;
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

    // Total areas should be within 10% of each other (bitmask with octahedral encoding is approximate)
    try std.testing.expectApproxEqRel(result_standard.total_area, result_bitmask.total_area, 0.10);

    // Per-atom areas should be within 15% tolerance
    // Octahedral grid quantization can cause larger per-atom deviations
    for (0..4) |i| {
        if (result_standard.atom_areas[i] > 0.1) {
            try std.testing.expectApproxEqRel(
                result_standard.atom_areas[i],
                result_bitmask.atom_areas[i],
                0.15,
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

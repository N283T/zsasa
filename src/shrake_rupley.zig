const std = @import("std");
const types = @import("types.zig");
const test_points = @import("test_points.zig");
const neighbor_list = @import("neighbor_list.zig");
const simd = @import("simd.zig");
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
const Allocator = std.mem.Allocator;

/// Calculate SASA for a single atom using the Shrake-Rupley algorithm.
///
/// Algorithm (from FreeSASA sasa_sr.c sr_atom_area, lines 276-338):
/// 1. Scale test points by (r[i] + probe_radius)
/// 2. Translate to atom position (x[i], y[i], z[i])
/// 3. For each test point:
///    - Check if it's inside any other atom j (distance² < (r[j] + probe_radius)²)
///    - If not inside any atom, count as exposed
/// 4. SASA[i] = 4 * PI * (r[i] + probe_radius)² * (exposed_points / total_points)
///
/// # Parameters
/// - `atom_idx`: Index of the atom to calculate SASA for
/// - `positions`: Array of atom positions (Vec3)
/// - `radii`: Array of atom radii
/// - `test_points`: Array of unit sphere test points
/// - `probe_radius`: Water probe radius in Angstroms
///
/// # Returns
/// SASA value for the atom in Ų
pub fn atomSasa(
    atom_idx: usize,
    positions: []const Vec3,
    radii: []const f64,
    test_points_array: []const Vec3,
    probe_radius: f64,
) f64 {
    const n_atoms = positions.len;
    const n_points = test_points_array.len;

    // Get atom properties
    const atom_pos = positions[atom_idx];
    const atom_radius = radii[atom_idx];
    const atom_radius_probe = atom_radius + probe_radius;

    // Count exposed test points
    var n_exposed: usize = 0;

    for (test_points_array) |test_point| {
        // Scale and translate test point to atom surface
        const scaled = test_point.scale(atom_radius_probe);
        const point = atom_pos.add(scaled);

        // Check if this point is inside any other atom
        var is_buried = false;
        for (0..n_atoms) |j| {
            if (j == atom_idx) continue; // Skip self

            const other_pos = positions[j];
            const other_radius_probe = radii[j] + probe_radius;

            // Calculate squared distance
            const dx = point.x - other_pos.x;
            const dy = point.y - other_pos.y;
            const dz = point.z - other_pos.z;
            const dist_sq = dx * dx + dy * dy + dz * dz;

            // Check if point is inside this atom
            const radius_sq = other_radius_probe * other_radius_probe;
            // use <= to mark boundary points as buried
            if (dist_sq < radius_sq) {
                is_buried = true;
                break;
            }
        }

        if (!is_buried) {
            n_exposed += 1;
        }
    }

    // Calculate SASA: 4π * r² * (exposed / total)
    const surface_area = 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
    const exposed_fraction = @as(f64, @floatFromInt(n_exposed)) / @as(f64, @floatFromInt(n_points));
    return surface_area * exposed_fraction;
}

/// Calculate SASA for a single atom using pre-computed neighbor list (O(k) instead of O(N))
///
/// This is the optimized version that only checks neighbors instead of all atoms.
/// The neighbor list must be pre-computed using NeighborList.init().
/// Uses compile-time CPU feature detection to select optimal SIMD width:
/// - AVX-512: 16 → 8 → 4 → scalar
/// - AVX2/NEON: 8 → 4 → scalar
/// - Fallback: 4 → scalar
///
/// # Parameters
/// - `atom_idx`: Index of the atom to calculate SASA for
/// - `positions`: Array of atom positions (Vec3)
/// - `radii_with_probe_sq`: Pre-computed (r[j] + probe_radius)² for each atom
/// - `test_points_array`: Array of unit sphere test points
/// - `atom_radius_probe`: Pre-computed (r[i] + probe_radius) for this atom
/// - `neighbors`: Pre-computed list of neighbor indices for this atom
///
/// # Returns
/// SASA value for the atom in Ų
fn atomSasaWithNeighbors(
    atom_idx: usize,
    positions: []const Vec3,
    radii_with_probe_sq: []const f64,
    test_points_array: []const Vec3,
    atom_radius_probe: f64,
    neighbors: []const u32,
) f64 {
    const n_points = test_points_array.len;

    // Early exit: no neighbors means full surface is exposed
    if (neighbors.len == 0) {
        return 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
    }

    // Get atom position
    const atom_pos = positions[atom_idx];

    // Count exposed test points
    var n_exposed: usize = 0;

    for (test_points_array) |test_point| {
        // Scale and translate test point to atom surface
        const scaled = test_point.scale(atom_radius_probe);
        const point = atom_pos.add(scaled);

        // Check if this point is inside any neighbor atom using SIMD
        var is_buried = false;
        var i: usize = 0;

        // AVX-512: Process 16 neighbors at a time with 16-wide SIMD
        if (comptime simd.cpu_features.has_avx512f) {
            while (i + 16 <= neighbors.len) : (i += 16) {
                const batch_positions = [16]Vec3{
                    positions[neighbors[i]],
                    positions[neighbors[i + 1]],
                    positions[neighbors[i + 2]],
                    positions[neighbors[i + 3]],
                    positions[neighbors[i + 4]],
                    positions[neighbors[i + 5]],
                    positions[neighbors[i + 6]],
                    positions[neighbors[i + 7]],
                    positions[neighbors[i + 8]],
                    positions[neighbors[i + 9]],
                    positions[neighbors[i + 10]],
                    positions[neighbors[i + 11]],
                    positions[neighbors[i + 12]],
                    positions[neighbors[i + 13]],
                    positions[neighbors[i + 14]],
                    positions[neighbors[i + 15]],
                };
                const batch_radii = [16]f64{
                    radii_with_probe_sq[neighbors[i]],
                    radii_with_probe_sq[neighbors[i + 1]],
                    radii_with_probe_sq[neighbors[i + 2]],
                    radii_with_probe_sq[neighbors[i + 3]],
                    radii_with_probe_sq[neighbors[i + 4]],
                    radii_with_probe_sq[neighbors[i + 5]],
                    radii_with_probe_sq[neighbors[i + 6]],
                    radii_with_probe_sq[neighbors[i + 7]],
                    radii_with_probe_sq[neighbors[i + 8]],
                    radii_with_probe_sq[neighbors[i + 9]],
                    radii_with_probe_sq[neighbors[i + 10]],
                    radii_with_probe_sq[neighbors[i + 11]],
                    radii_with_probe_sq[neighbors[i + 12]],
                    radii_with_probe_sq[neighbors[i + 13]],
                    radii_with_probe_sq[neighbors[i + 14]],
                    radii_with_probe_sq[neighbors[i + 15]],
                };

                if (simd.isPointBuriedBatch16(point, batch_positions, batch_radii)) {
                    is_buried = true;
                    break;
                }
            }
        }

        // Process 8 neighbors at a time with 8-wide SIMD
        if (!is_buried) {
            while (i + 8 <= neighbors.len) : (i += 8) {
                const batch_positions = [8]Vec3{
                    positions[neighbors[i]],
                    positions[neighbors[i + 1]],
                    positions[neighbors[i + 2]],
                    positions[neighbors[i + 3]],
                    positions[neighbors[i + 4]],
                    positions[neighbors[i + 5]],
                    positions[neighbors[i + 6]],
                    positions[neighbors[i + 7]],
                };
                const batch_radii = [8]f64{
                    radii_with_probe_sq[neighbors[i]],
                    radii_with_probe_sq[neighbors[i + 1]],
                    radii_with_probe_sq[neighbors[i + 2]],
                    radii_with_probe_sq[neighbors[i + 3]],
                    radii_with_probe_sq[neighbors[i + 4]],
                    radii_with_probe_sq[neighbors[i + 5]],
                    radii_with_probe_sq[neighbors[i + 6]],
                    radii_with_probe_sq[neighbors[i + 7]],
                };

                if (simd.isPointBuriedBatch8(point, batch_positions, batch_radii)) {
                    is_buried = true;
                    break;
                }
            }
        }

        // Process remaining 4-7 neighbors with 4-wide SIMD
        if (!is_buried and i + 4 <= neighbors.len) {
            const batch_positions = [4]Vec3{
                positions[neighbors[i]],
                positions[neighbors[i + 1]],
                positions[neighbors[i + 2]],
                positions[neighbors[i + 3]],
            };
            const batch_radii = [4]f64{
                radii_with_probe_sq[neighbors[i]],
                radii_with_probe_sq[neighbors[i + 1]],
                radii_with_probe_sq[neighbors[i + 2]],
                radii_with_probe_sq[neighbors[i + 3]],
            };

            if (simd.isPointBuriedBatch4(point, batch_positions, batch_radii)) {
                is_buried = true;
            }
            i += 4;
        }

        // Handle remaining neighbors (0-3) with scalar fallback
        if (!is_buried) {
            while (i < neighbors.len) : (i += 1) {
                const j = neighbors[i];
                const other_pos = positions[j];

                // Calculate squared distance
                const dx = point.x - other_pos.x;
                const dy = point.y - other_pos.y;
                const dz = point.z - other_pos.z;
                const dist_sq = dx * dx + dy * dy + dz * dz;

                // Check if point is inside this atom (strict less-than for open boundary)
                if (dist_sq < radii_with_probe_sq[j]) {
                    is_buried = true;
                    break;
                }
            }
        }

        if (!is_buried) {
            n_exposed += 1;
        }
    }

    // Calculate SASA: 4π * r² * (exposed / total)
    const surface_area = 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
    const exposed_fraction = @as(f64, @floatFromInt(n_exposed)) / @as(f64, @floatFromInt(n_points));
    return surface_area * exposed_fraction;
}

/// Calculate SASA for all atoms in the system.
///
/// # Parameters
/// - `allocator`: Memory allocator for result arrays
/// - `input`: Atom input data (positions and radii)
/// - `config`: Configuration parameters (n_points, probe_radius)
///
/// # Returns
/// SasaResult containing total_area and per-atom areas. Caller must call deinit().
pub fn calculateSasa(
    allocator: Allocator,
    input: AtomInput,
    config: Config,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) return error.NoAtoms;

    // Generate test points
    const test_points_array = try test_points.generateTestPoints(allocator, config.n_points);
    defer allocator.free(test_points_array);

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

    // Build neighbor list for O(N) instead of O(N²) neighbor checking
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

    // Calculate SASA for each atom using neighbor list
    var total_area: f64 = 0.0;
    for (0..n_atoms) |i| {
        const atom_radius_probe = input.r[i] + config.probe_radius;
        const neighbors = neighbor_list_data.getNeighbors(i);
        const area = atomSasaWithNeighbors(
            i,
            positions,
            radii_with_probe_sq,
            test_points_array,
            atom_radius_probe,
            neighbors,
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

/// Context for parallel SASA calculation workers.
/// Thread safety: All fields are read-only except `atom_areas` which has
/// disjoint write access (each thread writes to different indices).
const ParallelContext = struct {
    positions: []const Vec3,
    radii: []const f64,
    radii_with_probe_sq: []const f64,
    test_points_array: []const Vec3,
    neighbor_list_data: *const NeighborList,
    probe_radius: f64,
    atom_areas: []f64,
};

/// Worker function for parallel SASA calculation.
/// Processes atoms from chunk_start to chunk_end.
fn parallelSasaWorker(ctx: ParallelContext, chunk_start: usize, chunk_end: usize) f64 {
    var chunk_total: f64 = 0.0;

    for (chunk_start..chunk_end) |i| {
        const atom_radius_probe = ctx.radii[i] + ctx.probe_radius;
        const neighbors = ctx.neighbor_list_data.getNeighbors(i);
        const area = atomSasaWithNeighbors(
            i,
            ctx.positions,
            ctx.radii_with_probe_sq,
            ctx.test_points_array,
            atom_radius_probe,
            neighbors,
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

/// Calculate SASA for all atoms using parallel processing.
///
/// # Parameters
/// - `allocator`: Memory allocator for result arrays
/// - `input`: Atom input data (positions and radii)
/// - `config`: Configuration parameters (n_points, probe_radius)
/// - `n_threads`: Number of worker threads (0 = auto-detect)
///
/// # Returns
/// SasaResult containing total_area and per-atom areas. Caller must call deinit().
pub fn calculateSasaParallel(
    allocator: Allocator,
    input: AtomInput,
    config: Config,
    n_threads: usize,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) return error.NoAtoms;

    // Auto-detect thread count if 0
    const actual_threads = if (n_threads == 0)
        try std.Thread.getCpuCount()
    else
        n_threads;

    // Generate test points
    const test_points_array = try test_points.generateTestPoints(allocator, config.n_points);
    defer allocator.free(test_points_array);

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

    // Build neighbor list for O(N) instead of O(N²) neighbor checking
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
        .test_points_array = test_points_array,
        .neighbor_list_data = &neighbor_list_data,
        .probe_radius = config.probe_radius,
        .atom_areas = atom_areas,
    };

    // Chunk size heuristic:
    // - Minimum 64 atoms per chunk to amortize thread overhead
    // - Target 4 chunks per thread for load balancing (work stealing)
    const chunk_size = @max(64, n_atoms / (actual_threads * 4));

    // Run parallel calculation
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
// Generic implementations for f32/f64 precision support
// =============================================================================

/// Generic Shrake-Rupley algorithm implementation supporting both f32 and f64 precision.
/// Uses compile-time CPU feature detection to select optimal SIMD width.
pub fn ShrakeRupleyGen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    const Result = SasaResultGen(T);
    const Cfg = ConfigGen(T);
    const NList = NeighborListGen(T);
    const TestPointsGen = test_points.generateTestPointsGen(T);
    const IsPointBuriedBatch16 = simd.isPointBuriedBatch16Gen(T);
    const IsPointBuriedBatch8 = simd.isPointBuriedBatch8Gen(T);
    const IsPointBuriedBatch4 = simd.isPointBuriedBatch4Gen(T);

    return struct {
        const Self = @This();

        /// Calculate SASA for a single atom using pre-computed neighbor list.
        /// Uses compile-time CPU feature detection to select optimal SIMD width:
        /// - AVX-512: 16 → 8 → 4 → scalar
        /// - AVX2/NEON: 8 → 4 → scalar
        /// - Fallback: 4 → scalar
        pub fn atomSasaWithNeighbors(
            atom_idx: usize,
            positions: []const Vec,
            radii_with_probe_sq: []const T,
            test_points_array: []const Vec,
            atom_radius_probe: T,
            neighbors: []const u32,
        ) T {
            const n_points = test_points_array.len;

            // Early exit: no neighbors means full surface is exposed
            if (neighbors.len == 0) {
                return 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
            }

            // Get atom position
            const atom_pos = positions[atom_idx];

            // Count exposed test points
            var n_exposed: usize = 0;

            for (test_points_array) |test_point| {
                // Scale and translate test point to atom surface
                const scaled = test_point.scale(atom_radius_probe);
                const point = atom_pos.add(scaled);

                // Check if this point is inside any neighbor atom using SIMD
                var is_buried = false;
                var i: usize = 0;

                // AVX-512: Process 16 neighbors at a time with 16-wide SIMD
                if (comptime simd.cpu_features.has_avx512f) {
                    while (i + 16 <= neighbors.len) : (i += 16) {
                        const batch_positions = [16]Vec{
                            positions[neighbors[i]],
                            positions[neighbors[i + 1]],
                            positions[neighbors[i + 2]],
                            positions[neighbors[i + 3]],
                            positions[neighbors[i + 4]],
                            positions[neighbors[i + 5]],
                            positions[neighbors[i + 6]],
                            positions[neighbors[i + 7]],
                            positions[neighbors[i + 8]],
                            positions[neighbors[i + 9]],
                            positions[neighbors[i + 10]],
                            positions[neighbors[i + 11]],
                            positions[neighbors[i + 12]],
                            positions[neighbors[i + 13]],
                            positions[neighbors[i + 14]],
                            positions[neighbors[i + 15]],
                        };
                        const batch_radii = [16]T{
                            radii_with_probe_sq[neighbors[i]],
                            radii_with_probe_sq[neighbors[i + 1]],
                            radii_with_probe_sq[neighbors[i + 2]],
                            radii_with_probe_sq[neighbors[i + 3]],
                            radii_with_probe_sq[neighbors[i + 4]],
                            radii_with_probe_sq[neighbors[i + 5]],
                            radii_with_probe_sq[neighbors[i + 6]],
                            radii_with_probe_sq[neighbors[i + 7]],
                            radii_with_probe_sq[neighbors[i + 8]],
                            radii_with_probe_sq[neighbors[i + 9]],
                            radii_with_probe_sq[neighbors[i + 10]],
                            radii_with_probe_sq[neighbors[i + 11]],
                            radii_with_probe_sq[neighbors[i + 12]],
                            radii_with_probe_sq[neighbors[i + 13]],
                            radii_with_probe_sq[neighbors[i + 14]],
                            radii_with_probe_sq[neighbors[i + 15]],
                        };

                        if (IsPointBuriedBatch16.compute(point, batch_positions, batch_radii)) {
                            is_buried = true;
                            break;
                        }
                    }
                }

                // Process 8 neighbors at a time with 8-wide SIMD
                if (!is_buried) {
                    while (i + 8 <= neighbors.len) : (i += 8) {
                        const batch_positions = [8]Vec{
                            positions[neighbors[i]],
                            positions[neighbors[i + 1]],
                            positions[neighbors[i + 2]],
                            positions[neighbors[i + 3]],
                            positions[neighbors[i + 4]],
                            positions[neighbors[i + 5]],
                            positions[neighbors[i + 6]],
                            positions[neighbors[i + 7]],
                        };
                        const batch_radii = [8]T{
                            radii_with_probe_sq[neighbors[i]],
                            radii_with_probe_sq[neighbors[i + 1]],
                            radii_with_probe_sq[neighbors[i + 2]],
                            radii_with_probe_sq[neighbors[i + 3]],
                            radii_with_probe_sq[neighbors[i + 4]],
                            radii_with_probe_sq[neighbors[i + 5]],
                            radii_with_probe_sq[neighbors[i + 6]],
                            radii_with_probe_sq[neighbors[i + 7]],
                        };

                        if (IsPointBuriedBatch8.compute(point, batch_positions, batch_radii)) {
                            is_buried = true;
                            break;
                        }
                    }
                }

                // Process remaining 4-7 neighbors with 4-wide SIMD
                if (!is_buried and i + 4 <= neighbors.len) {
                    const batch_positions = [4]Vec{
                        positions[neighbors[i]],
                        positions[neighbors[i + 1]],
                        positions[neighbors[i + 2]],
                        positions[neighbors[i + 3]],
                    };
                    const batch_radii = [4]T{
                        radii_with_probe_sq[neighbors[i]],
                        radii_with_probe_sq[neighbors[i + 1]],
                        radii_with_probe_sq[neighbors[i + 2]],
                        radii_with_probe_sq[neighbors[i + 3]],
                    };

                    if (IsPointBuriedBatch4.compute(point, batch_positions, batch_radii)) {
                        is_buried = true;
                    }
                    i += 4;
                }

                // Handle remaining neighbors (0-3) with scalar fallback
                if (!is_buried) {
                    while (i < neighbors.len) : (i += 1) {
                        const j = neighbors[i];
                        const other_pos = positions[j];

                        // Calculate squared distance
                        const dx = point.x - other_pos.x;
                        const dy = point.y - other_pos.y;
                        const dz = point.z - other_pos.z;
                        const dist_sq = dx * dx + dy * dy + dz * dz;

                        // Check if point is inside this atom (strict less-than for open boundary)
                        if (dist_sq < radii_with_probe_sq[j]) {
                            is_buried = true;
                            break;
                        }
                    }
                }

                if (!is_buried) {
                    n_exposed += 1;
                }
            }

            // Calculate SASA: 4π * r² * (exposed / total)
            const surface_area = 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;
            const exposed_fraction = @as(T, @floatFromInt(n_exposed)) / @as(T, @floatFromInt(n_points));
            return surface_area * exposed_fraction;
        }

        /// Context for parallel SASA calculation workers.
        pub const ParallelContext = struct {
            positions: []const Vec,
            radii: []const T,
            radii_with_probe_sq: []const T,
            test_points_array: []const Vec,
            neighbor_list_data: *const NList,
            probe_radius: T,
            atom_areas: []T,
        };

        /// Worker function for parallel SASA calculation.
        fn parallelSasaWorker(ctx: Self.ParallelContext, chunk_start: usize, chunk_end: usize) T {
            var chunk_total: T = 0.0;

            for (chunk_start..chunk_end) |i| {
                const atom_radius_probe = ctx.radii[i] + ctx.probe_radius;
                const neighbors = ctx.neighbor_list_data.getNeighbors(i);
                const area = Self.atomSasaWithNeighbors(
                    i,
                    ctx.positions,
                    ctx.radii_with_probe_sq,
                    ctx.test_points_array,
                    atom_radius_probe,
                    neighbors,
                );
                ctx.atom_areas[i] = area;
                chunk_total += area;
            }

            return chunk_total;
        }

        /// Reduce function to sum all chunk totals.
        fn sumReducer(results: []const T) T {
            var total: T = 0.0;
            for (results) |r| {
                total += r;
            }
            return total;
        }

        /// Calculate SASA for all atoms in the system (single-threaded).
        pub fn calculateSasa(
            allocator: Allocator,
            input: AtomInput,
            config: Cfg,
        ) !Result {
            const n_atoms = input.atomCount();
            if (n_atoms == 0) return error.NoAtoms;

            // Generate test points
            const test_points_array = try TestPointsGen.generate(allocator, config.n_points);
            defer allocator.free(test_points_array);

            // Convert input to Vec positions (cast from f64)
            const positions = try allocator.alloc(Vec, n_atoms);
            defer allocator.free(positions);

            for (0..n_atoms) |i| {
                positions[i] = Vec{
                    .x = @floatCast(input.x[i]),
                    .y = @floatCast(input.y[i]),
                    .z = @floatCast(input.z[i]),
                };
            }

            // Convert radii (cast from f64)
            const radii = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii);
            for (0..n_atoms) |i| {
                radii[i] = @floatCast(input.r[i]);
            }

            // Build neighbor list
            var neighbor_list_data = try NList.init(allocator, positions, radii, config.probe_radius);
            defer neighbor_list_data.deinit();

            // Pre-compute (r[j] + probe_radius)² for each atom
            const radii_with_probe_sq = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii_with_probe_sq);

            for (0..n_atoms) |i| {
                const r_probe = radii[i] + config.probe_radius;
                radii_with_probe_sq[i] = r_probe * r_probe;
            }

            // Allocate result arrays
            const atom_areas = try allocator.alloc(T, n_atoms);
            errdefer allocator.free(atom_areas);

            // Calculate SASA for each atom using neighbor list
            var total_area: T = 0.0;
            for (0..n_atoms) |i| {
                const atom_radius_probe = radii[i] + config.probe_radius;
                const neighbors = neighbor_list_data.getNeighbors(i);
                const area = Self.atomSasaWithNeighbors(
                    i,
                    positions,
                    radii_with_probe_sq,
                    test_points_array,
                    atom_radius_probe,
                    neighbors,
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

        /// Calculate SASA for all atoms using parallel processing.
        pub fn calculateSasaParallel(
            allocator: Allocator,
            input: AtomInput,
            config: Cfg,
            n_threads: usize,
        ) !Result {
            const n_atoms = input.atomCount();
            if (n_atoms == 0) return error.NoAtoms;

            // Auto-detect thread count if 0
            const actual_threads = if (n_threads == 0)
                try std.Thread.getCpuCount()
            else
                n_threads;

            // Generate test points
            const test_points_array = try TestPointsGen.generate(allocator, config.n_points);
            defer allocator.free(test_points_array);

            // Convert input to Vec positions (cast from f64)
            const positions = try allocator.alloc(Vec, n_atoms);
            defer allocator.free(positions);

            for (0..n_atoms) |i| {
                positions[i] = Vec{
                    .x = @floatCast(input.x[i]),
                    .y = @floatCast(input.y[i]),
                    .z = @floatCast(input.z[i]),
                };
            }

            // Convert radii (cast from f64)
            const radii = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii);
            for (0..n_atoms) |i| {
                radii[i] = @floatCast(input.r[i]);
            }

            // Build neighbor list
            var neighbor_list_data = try NList.init(allocator, positions, radii, config.probe_radius);
            defer neighbor_list_data.deinit();

            // Pre-compute (r[j] + probe_radius)² for each atom
            const radii_with_probe_sq = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii_with_probe_sq);

            for (0..n_atoms) |i| {
                const r_probe = radii[i] + config.probe_radius;
                radii_with_probe_sq[i] = r_probe * r_probe;
            }

            // Allocate result arrays
            const atom_areas = try allocator.alloc(T, n_atoms);
            errdefer allocator.free(atom_areas);

            // Create parallel context
            const ctx = Self.ParallelContext{
                .positions = positions,
                .radii = radii,
                .radii_with_probe_sq = radii_with_probe_sq,
                .test_points_array = test_points_array,
                .neighbor_list_data = &neighbor_list_data,
                .probe_radius = config.probe_radius,
                .atom_areas = atom_areas,
            };

            // Chunk size heuristic:
            // - Minimum 64 atoms per chunk to amortize thread overhead
            // - Target 4 chunks per thread for load balancing
            const chunk_size = @max(64, n_atoms / (actual_threads * 4));

            // Run parallel calculation
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

/// f32 precision Shrake-Rupley implementation
pub const ShrakeRupleyf32 = ShrakeRupleyGen(f32);

/// Convenience function: Calculate SASA using f32 precision (single-threaded)
pub fn calculateSasaf32(allocator: Allocator, input: AtomInput, config: ConfigGen(f32)) !SasaResultGen(f32) {
    return ShrakeRupleyf32.calculateSasa(allocator, input, config);
}

/// Convenience function: Calculate SASA using f32 precision (parallel)
pub fn calculateSasaParallelf32(allocator: Allocator, input: AtomInput, config: ConfigGen(f32), n_threads: usize) !SasaResultGen(f32) {
    return ShrakeRupleyf32.calculateSasaParallel(allocator, input, config, n_threads);
}

// Tests

test "atomSasa - single isolated atom" {
    const allocator = std.testing.allocator;

    // Single atom at origin
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
    };
    const radii = &[_]f64{1.0};
    const probe_radius = 1.4;

    // Generate test points
    const test_points_array = try test_points.generateTestPoints(allocator, 1000);
    defer allocator.free(test_points_array);

    const sasa = atomSasa(0, positions, radii, test_points_array, probe_radius);

    // Expected: 4π * (r + probe)²
    const expected_radius = radii[0] + probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    // Should be very close with 1000 test points
    try std.testing.expectApproxEqRel(expected_area, sasa, 0.01);
}

test "atomSasa - two far apart atoms" {
    const allocator = std.testing.allocator;

    // Two atoms far apart (no overlap)
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 100.0, .y = 0.0, .z = 0.0 }, // Very far away
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    const test_points_array = try test_points.generateTestPoints(allocator, 1000);
    defer allocator.free(test_points_array);

    // Each atom should have full SASA (no burial)
    const sasa1 = atomSasa(0, positions, radii, test_points_array, probe_radius);
    const sasa2 = atomSasa(1, positions, radii, test_points_array, probe_radius);

    const expected_radius = radii[0] + probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    try std.testing.expectApproxEqRel(expected_area, sasa1, 0.01);
    try std.testing.expectApproxEqRel(expected_area, sasa2, 0.01);
}

test "atomSasa - two overlapping atoms" {
    const allocator = std.testing.allocator;

    // Two atoms close together (overlapping)
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 2.0, .y = 0.0, .z = 0.0 }, // Close enough to bury some surface
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    const test_points_array = try test_points.generateTestPoints(allocator, 1000);
    defer allocator.free(test_points_array);

    const sasa1 = atomSasa(0, positions, radii, test_points_array, probe_radius);

    // Expected: less than full SASA due to burial
    const expected_radius = radii[0] + probe_radius;
    const full_area = 4.0 * std.math.pi * expected_radius * expected_radius;

    // SASA should be less than full area
    try std.testing.expect(sasa1 < full_area);
    // But still positive
    try std.testing.expect(sasa1 > 0.0);
}

test "calculateSasa - single atom" {
    const allocator = std.testing.allocator;

    // Create input for single atom
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
        .n_points = 1000,
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

test "calculateSasa - two far apart atoms" {
    const allocator = std.testing.allocator;

    // Create input for two atoms
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
        .n_points = 1000,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    // Each atom should have full SASA
    const expected_radius = r[0] + config.probe_radius;
    const expected_area = 4.0 * std.math.pi * expected_radius * expected_radius;
    const expected_total = expected_area * 2.0;

    try std.testing.expectEqual(@as(usize, 2), result.atom_areas.len);
    try std.testing.expectApproxEqRel(expected_total, result.total_area, 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[0], 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[1], 0.01);
}

test "calculateSasa - two overlapping atoms" {
    const allocator = std.testing.allocator;

    // Create input for two overlapping atoms
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
        .n_points = 1000,
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
    // Each atom should have reduced SASA
    try std.testing.expect(result.atom_areas[0] < full_area);
    try std.testing.expect(result.atom_areas[1] < full_area);
}

test "calculateSasa - no atoms error" {
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

    const config = Config{};

    const result = calculateSasa(allocator, input, config);
    try std.testing.expectError(error.NoAtoms, result);
}

test "optimized vs original - same results" {
    // Verify the neighbor-list optimized implementation produces
    // identical results to the O(N²) implementation
    const allocator = std.testing.allocator;

    // Create a cluster of 10 atoms
    const n_atoms = 10;
    const positions = try allocator.alloc(Vec3, n_atoms);
    defer allocator.free(positions);
    const radii = try allocator.alloc(f64, n_atoms);
    defer allocator.free(radii);

    // Random-ish positions in a 10 Å cube
    positions[0] = Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 };
    positions[1] = Vec3{ .x = 2.5, .y = 0.0, .z = 0.0 };
    positions[2] = Vec3{ .x = 5.0, .y = 1.0, .z = 0.0 };
    positions[3] = Vec3{ .x = 1.0, .y = 3.0, .z = 0.5 };
    positions[4] = Vec3{ .x = 3.5, .y = 2.5, .z = 1.0 };
    positions[5] = Vec3{ .x = 0.5, .y = 1.5, .z = 3.0 };
    positions[6] = Vec3{ .x = 4.0, .y = 0.5, .z = 2.5 };
    positions[7] = Vec3{ .x = 2.0, .y = 4.0, .z = 2.0 };
    positions[8] = Vec3{ .x = 6.0, .y = 3.0, .z = 1.5 };
    positions[9] = Vec3{ .x = 7.0, .y = 1.0, .z = 3.0 };

    // Different radii to test the cutoff calculation
    radii[0] = 1.0;
    radii[1] = 1.2;
    radii[2] = 0.8;
    radii[3] = 1.5;
    radii[4] = 1.0;
    radii[5] = 1.1;
    radii[6] = 0.9;
    radii[7] = 1.3;
    radii[8] = 1.0;
    radii[9] = 1.2;

    const probe_radius = 1.4;
    const test_points_array = try test_points.generateTestPoints(allocator, 500);
    defer allocator.free(test_points_array);

    // Build neighbor list for optimized implementation
    var neighbor_list_data = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list_data.deinit();

    // Pre-compute radii_with_probe_sq
    const radii_with_probe_sq = try allocator.alloc(f64, n_atoms);
    defer allocator.free(radii_with_probe_sq);
    for (0..n_atoms) |i| {
        const r_probe = radii[i] + probe_radius;
        radii_with_probe_sq[i] = r_probe * r_probe;
    }

    // Compare results for each atom
    for (0..n_atoms) |i| {
        const atom_radius_probe = radii[i] + probe_radius;
        const neighbors = neighbor_list_data.getNeighbors(i);

        // Original O(N²) implementation
        const original_sasa = atomSasa(i, positions, radii, test_points_array, probe_radius);

        // Optimized O(k) implementation
        const optimized_sasa = atomSasaWithNeighbors(
            i,
            positions,
            radii_with_probe_sq,
            test_points_array,
            atom_radius_probe,
            neighbors,
        );

        // Should be exactly equal (same algorithm, just checking fewer atoms)
        try std.testing.expectApproxEqAbs(original_sasa, optimized_sasa, 1e-10);
    }
}

test "atomSasaWithNeighbors - handles 0-3 neighbors correctly (scalar fallback)" {
    // Test edge cases where SIMD loop is skipped and only scalar fallback runs
    const allocator = std.testing.allocator;

    const positions = [_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 }, // atom 0 (isolated)
        Vec3{ .x = 100.0, .y = 0.0, .z = 0.0 }, // atom 1 (far away)
        Vec3{ .x = 100.0, .y = 100.0, .z = 0.0 }, // atom 2 (far away)
        Vec3{ .x = 100.0, .y = 100.0, .z = 100.0 }, // atom 3 (far away)
    };
    const radii = [_]f64{ 1.0, 1.0, 1.0, 1.0 };
    const probe_radius = 1.4;

    const test_points_array = try test_points.generateTestPoints(allocator, 100);
    defer allocator.free(test_points_array);

    // Pre-compute radii_with_probe_sq
    var radii_with_probe_sq: [4]f64 = undefined;
    for (0..4) |i| {
        const r_probe = radii[i] + probe_radius;
        radii_with_probe_sq[i] = r_probe * r_probe;
    }

    const atom_radius_probe = radii[0] + probe_radius;
    const expected_full_sasa = 4.0 * std.math.pi * atom_radius_probe * atom_radius_probe;

    // Test with 0 neighbors (completely isolated)
    const neighbors_0 = [_]u32{};
    const sasa_0 = atomSasaWithNeighbors(0, &positions, &radii_with_probe_sq, test_points_array, atom_radius_probe, &neighbors_0);
    try std.testing.expectApproxEqRel(expected_full_sasa, sasa_0, 0.01);

    // Test with 1 neighbor (far away, should still have full SASA)
    const neighbors_1 = [_]u32{1};
    const sasa_1 = atomSasaWithNeighbors(0, &positions, &radii_with_probe_sq, test_points_array, atom_radius_probe, &neighbors_1);
    try std.testing.expectApproxEqRel(expected_full_sasa, sasa_1, 0.01);

    // Test with 2 neighbors (far away)
    const neighbors_2 = [_]u32{ 1, 2 };
    const sasa_2 = atomSasaWithNeighbors(0, &positions, &radii_with_probe_sq, test_points_array, atom_radius_probe, &neighbors_2);
    try std.testing.expectApproxEqRel(expected_full_sasa, sasa_2, 0.01);

    // Test with 3 neighbors (far away)
    const neighbors_3 = [_]u32{ 1, 2, 3 };
    const sasa_3 = atomSasaWithNeighbors(0, &positions, &radii_with_probe_sq, test_points_array, atom_radius_probe, &neighbors_3);
    try std.testing.expectApproxEqRel(expected_full_sasa, sasa_3, 0.01);
}

test "calculateSasaParallel - same results as sequential" {
    const allocator = std.testing.allocator;

    // Create a cluster of 20 atoms for more meaningful parallel test
    const n_atoms = 20;
    const x = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r);

    // Create a grid of atoms
    var idx: usize = 0;
    for (0..4) |ix| {
        for (0..5) |iy| {
            x[idx] = @as(f64, @floatFromInt(ix)) * 3.0;
            y[idx] = @as(f64, @floatFromInt(iy)) * 3.0;
            z[idx] = 0.0;
            r[idx] = 1.0 + @as(f64, @floatFromInt(idx % 5)) * 0.1;
            idx += 1;
        }
    }

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const config = Config{
        .n_points = 100,
        .probe_radius = 1.4,
    };

    // Calculate sequential
    var sequential_result = try calculateSasa(allocator, input, config);
    defer sequential_result.deinit();

    // Calculate parallel with 2 threads
    var parallel_result = try calculateSasaParallel(allocator, input, config, 2);
    defer parallel_result.deinit();

    // Verify total area matches
    try std.testing.expectApproxEqAbs(sequential_result.total_area, parallel_result.total_area, 1e-10);

    // Verify each atom area matches
    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqAbs(
            sequential_result.atom_areas[i],
            parallel_result.atom_areas[i],
            1e-10,
        );
    }
}

test "calculateSasaParallel - auto thread count" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 2);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 2);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 2);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 2);
    defer allocator.free(r);

    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;
    x[1] = 100.0;
    y[1] = 0.0;
    z[1] = 0.0;
    r[1] = 1.0;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const config = Config{
        .n_points = 100,
        .probe_radius = 1.4,
    };

    // n_threads = 0 should auto-detect
    var result = try calculateSasaParallel(allocator, input, config, 0);
    defer result.deinit();

    // Just verify it runs and produces valid output
    try std.testing.expect(result.total_area > 0);
    try std.testing.expectEqual(@as(usize, 2), result.atom_areas.len);
}

test "calculateSasaParallel - no atoms error" {
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

    const config = Config{};
    const result = calculateSasaParallel(allocator, input, config, 2);
    try std.testing.expectError(error.NoAtoms, result);
}

test "calculateSasaParallel - small workload uses single-thread fallback" {
    // When n_atoms <= chunk_size, parallelFor falls back to single-threaded
    const allocator = std.testing.allocator;

    // Create only 3 atoms (less than minimum chunk size of 64)
    const n_atoms = 3;
    const x = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r);

    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;
    x[1] = 5.0;
    y[1] = 0.0;
    z[1] = 0.0;
    r[1] = 1.0;
    x[2] = 10.0;
    y[2] = 0.0;
    z[2] = 0.0;
    r[2] = 1.0;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const config = Config{
        .n_points = 100,
        .probe_radius = 1.4,
    };

    // Request 8 threads but with only 3 atoms, should use single-thread fallback
    var parallel_result = try calculateSasaParallel(allocator, input, config, 8);
    defer parallel_result.deinit();

    // Compare with sequential
    var sequential_result = try calculateSasa(allocator, input, config);
    defer sequential_result.deinit();

    // Results should be identical
    try std.testing.expectApproxEqAbs(sequential_result.total_area, parallel_result.total_area, 1e-10);
    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqAbs(
            sequential_result.atom_areas[i],
            parallel_result.atom_areas[i],
            1e-10,
        );
    }
}

// =============================================================================
// Tests for f32 precision implementation
// =============================================================================

test "calculateSasaf32 - single atom" {
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

    const config = ConfigGen(f32){
        .n_points = 1000,
        .probe_radius = 1.4,
    };

    var result = try calculateSasaf32(allocator, input, config);
    defer result.deinit();

    // Expected: 4π * (r + probe)²
    const expected_radius: f32 = 1.0 + config.probe_radius;
    const expected_area: f32 = 4.0 * std.math.pi * expected_radius * expected_radius;

    try std.testing.expectEqual(@as(usize, 1), result.atom_areas.len);
    try std.testing.expectApproxEqRel(expected_area, result.total_area, 0.01);
    try std.testing.expectApproxEqRel(expected_area, result.atom_areas[0], 0.01);
}

test "calculateSasaf32 - two overlapping atoms" {
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

    const config = ConfigGen(f32){
        .n_points = 1000,
        .probe_radius = 1.4,
    };

    var result = try calculateSasaf32(allocator, input, config);
    defer result.deinit();

    // Total should be less than 2 * full_area due to burial
    const expected_radius: f32 = 1.0 + config.probe_radius;
    const full_area: f32 = 4.0 * std.math.pi * expected_radius * expected_radius;
    const full_total = full_area * 2.0;

    try std.testing.expectEqual(@as(usize, 2), result.atom_areas.len);
    try std.testing.expect(result.total_area < full_total);
    try std.testing.expect(result.total_area > 0.0);
}

test "calculateSasaf32 vs f64 - similar results" {
    // f32 and f64 should produce similar results (within tolerance)
    const allocator = std.testing.allocator;

    const n_atoms = 10;
    const x = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r);

    // Create a cluster of atoms
    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;
    x[1] = 2.5;
    y[1] = 0.0;
    z[1] = 0.0;
    r[1] = 1.2;
    x[2] = 5.0;
    y[2] = 1.0;
    z[2] = 0.0;
    r[2] = 0.8;
    x[3] = 1.0;
    y[3] = 3.0;
    z[3] = 0.5;
    r[3] = 1.5;
    x[4] = 3.5;
    y[4] = 2.5;
    z[4] = 1.0;
    r[4] = 1.0;
    x[5] = 0.5;
    y[5] = 1.5;
    z[5] = 3.0;
    r[5] = 1.1;
    x[6] = 4.0;
    y[6] = 0.5;
    z[6] = 2.5;
    r[6] = 0.9;
    x[7] = 2.0;
    y[7] = 4.0;
    z[7] = 2.0;
    r[7] = 1.3;
    x[8] = 6.0;
    y[8] = 3.0;
    z[8] = 1.5;
    r[8] = 1.0;
    x[9] = 7.0;
    y[9] = 1.0;
    z[9] = 3.0;
    r[9] = 1.2;

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    // Calculate with f64
    const config64 = Config{
        .n_points = 500,
        .probe_radius = 1.4,
    };
    var result64 = try calculateSasa(allocator, input, config64);
    defer result64.deinit();

    // Calculate with f32
    const config32 = ConfigGen(f32){
        .n_points = 500,
        .probe_radius = 1.4,
    };
    var result32 = try calculateSasaf32(allocator, input, config32);
    defer result32.deinit();

    // f32 and f64 should produce similar results
    // Tolerance: 0.1% relative error (SASA calculation is stable)
    try std.testing.expectApproxEqRel(@as(f32, @floatCast(result64.total_area)), result32.total_area, 0.001);

    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqRel(
            @as(f32, @floatCast(result64.atom_areas[i])),
            result32.atom_areas[i],
            0.001,
        );
    }
}

test "calculateSasaParallelf32 - same as sequential f32" {
    const allocator = std.testing.allocator;

    const n_atoms = 20;
    const x = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r);

    // Create a grid of atoms
    var idx: usize = 0;
    for (0..4) |ix| {
        for (0..5) |iy| {
            x[idx] = @as(f64, @floatFromInt(ix)) * 3.0;
            y[idx] = @as(f64, @floatFromInt(iy)) * 3.0;
            z[idx] = 0.0;
            r[idx] = 1.0 + @as(f64, @floatFromInt(idx % 5)) * 0.1;
            idx += 1;
        }
    }

    const input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };

    const config = ConfigGen(f32){
        .n_points = 100,
        .probe_radius = 1.4,
    };

    // Calculate sequential
    var sequential_result = try calculateSasaf32(allocator, input, config);
    defer sequential_result.deinit();

    // Calculate parallel with 2 threads
    var parallel_result = try calculateSasaParallelf32(allocator, input, config, 2);
    defer parallel_result.deinit();

    // Verify total area matches
    try std.testing.expectApproxEqAbs(sequential_result.total_area, parallel_result.total_area, 1e-6);

    // Verify each atom area matches
    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqAbs(
            sequential_result.atom_areas[i],
            parallel_result.atom_areas[i],
            1e-6,
        );
    }
}

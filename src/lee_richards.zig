const std = @import("std");
const types = @import("types.zig");
const neighbor_list_mod = @import("neighbor_list.zig");
const thread_pool = @import("thread_pool.zig");
const simd = @import("simd.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const SasaResultGen = types.SasaResultGen;
const Vec3 = types.Vec3;
const Vec3Gen = types.Vec3Gen;
const NeighborList = neighbor_list_mod.NeighborList;
const NeighborListGen = neighbor_list_mod.NeighborListGen;

const TWOPI: f64 = 2.0 * std.math.pi;

/// How the angles of the arc that a neighbor covers on a slice circle are computed.
pub const TrigMode = enum {
    /// `std.math.acos` and `std.math.atan2` for every neighbor.
    exact,
    /// The polynomial approximations `simd.fastAcos` and `simd.fastAtan2` for
    /// the neighbors that fall in the 8-wide and 4-wide batches, exact
    /// trigonometry for the 0 to 3 neighbors left over. This is how
    /// Lee-Richards was computed up to zsasa 0.9.1. The approximations are off
    /// by up to 0.064 rad, which biases the total upwards by a few tenths of a
    /// percent, independently of `n_slices`, and makes the area of an atom
    /// depend on the order of its neighbors.
    fast,

    /// Parse the value of `--lr-trig` or of the workflow key `lr_trig`.
    pub fn fromString(value: []const u8) ?TrigMode {
        return std.meta.stringToEnum(TrigMode, value);
    }
};

/// Configuration for Lee-Richards algorithm
pub const LeeRichardsConfig = struct {
    /// Number of slices per atom diameter
    n_slices: u32 = 20,
    /// Water probe radius in Angstroms
    probe_radius: f64 = 1.4,
    /// Trigonometry used for the arc angles
    trig: TrigMode = .exact,
};

/// Arc interval representing a buried portion of a circle
///
/// A neighbor circle j covers the arc `[beta - alpha, beta + alpha]` of circle
/// i, with `cos(alpha) = (Ri'^2 + dij^2 - Rj'^2) / (2 Ri' dij)`. The angles are
/// then brought into [0, 2 pi] and an arc that crosses 0 is split in two. Two
/// tangent cases must not reach that step:
///
/// - `dij + Ri' = Rj'`: circle i touches circle j from inside and is buried.
///   Here `cos(alpha) = -1` and `alpha = pi`, the whole circle, but once both
///   ends are wrapped into [0, 2 pi] they can coincide (`start == end`), and
///   the neighbor that covers everything covers nothing.
/// - `dij + Rj' = Ri'`: circle j touches circle i from inside and covers
///   nothing. Here `cos(alpha) = 1` and `alpha = 0`, and with `beta = 2 pi`
///   the wrap turns the empty arc `[2 pi, 2 pi]` into `[0, 2 pi]`: full
///   burial from a neighbor that covers nothing.
///
/// `atomArea` therefore treats `dij + Ri' <= Rj'` and `cos(alpha) <= -1` as
/// buried, `dij + Rj' <= Ri'` and `cos(alpha) >= 1` as no arc (the tests on
/// `cos(alpha)` catch the pairs that rounding lets past the tests on the
/// radii), and drops any arc whose two ends are equal before they are wrapped.
const Arc = struct {
    start: f64, // Start angle (radians)
    end: f64, // End angle (radians)
};

/// Calculate SASA using Lee-Richards algorithm
pub fn calculateSasa(
    allocator: Allocator,
    input: AtomInput,
    config: LeeRichardsConfig,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) {
        return SasaResult{
            .total_area = 0.0,
            .atom_areas = try allocator.alloc(f64, 0),
            .allocator = allocator,
        };
    }
    if (config.n_slices == 0) return error.InvalidInput;
    try types.validateProbeRadius(f64, config.probe_radius);
    try input.validateFiniteAndRange();

    // Convert to Vec3 positions for neighbor list
    const positions = try allocator.alloc(Vec3, n_atoms);
    defer allocator.free(positions);
    for (0..n_atoms) |i| {
        positions[i] = Vec3{ .x = input.x[i], .y = input.y[i], .z = input.z[i] };
    }

    // Pre-compute effective radii (atom radius + probe radius)
    const radii = try allocator.alloc(f64, n_atoms);
    defer allocator.free(radii);
    for (0..n_atoms) |i| {
        radii[i] = input.r[i] + config.probe_radius;
    }

    // Build neighbor list with effective radii
    // Note: pass probe_radius=0 since we already added it to radii
    var neighbor_list = try NeighborList.init(allocator, positions, radii, 0.0);
    defer neighbor_list.deinit();

    // Calculate SASA for each atom
    const atom_areas = try allocator.alloc(f64, n_atoms);
    errdefer allocator.free(atom_areas);
    var total_area: f64 = 0.0;

    // Estimate max neighbors for arc buffer allocation
    var max_neighbors: usize = 0;
    for (0..n_atoms) |i| {
        max_neighbors = @max(max_neighbors, neighbor_list.getNeighbors(i).len);
    }
    // Each neighbor can create up to 2 arcs (when crossing 0)
    const arc_buffer = try allocator.alloc(Arc, (max_neighbors + 1) * 2);
    defer allocator.free(arc_buffer);

    for (0..n_atoms) |i| {
        atom_areas[i] = atomArea(
            i,
            input.x,
            input.y,
            input.z,
            radii,
            &neighbor_list,
            config.n_slices,
            config.trig,
            arc_buffer,
        );
        total_area += atom_areas[i];
    }

    return SasaResult{
        .total_area = total_area,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}

/// Half-angle of the arc that neighbor j covers on the slice circle of atom i,
/// for a neighbor of the 8-wide and 4-wide batches. The neighbors left over
/// after the batches use `std.math.acos` in both modes.
inline fn batchHalfAngle(trig: TrigMode, cos_alpha: f64) f64 {
    return switch (trig) {
        .exact => std.math.acos(std.math.clamp(cos_alpha, -1.0, 1.0)),
        .fast => simd.fastAcos(cos_alpha),
    };
}

/// Direction of neighbor j seen from atom i in the slice plane, for a neighbor
/// of the 8-wide and 4-wide batches. The neighbors left over after the batches
/// use `std.math.atan2` in both modes.
inline fn batchDirection(trig: TrigMode, dy: f64, dx: f64) f64 {
    return switch (trig) {
        .exact => std.math.atan2(dy, dx),
        .fast => simd.fastAtan2(dy, dx),
    };
}

/// Calculate SASA for a single atom using slice-based method with SIMD optimization
fn atomArea(
    atom_idx: usize,
    x: []const f64,
    y: []const f64,
    z: []const f64,
    radii: []const f64,
    neighbor_list: *const NeighborList,
    n_slices: u32,
    trig: TrigMode,
    arc_buffer: []Arc,
) f64 {
    const xi = x[atom_idx];
    const yi = y[atom_idx];
    const zi = z[atom_idx];
    const Ri = radii[atom_idx];

    const neighbors = neighbor_list.getNeighbors(atom_idx);
    if (neighbors.len == 0) {
        // No neighbors, full sphere exposed
        return 4.0 * std.math.pi * Ri * Ri;
    }

    const delta = 2.0 * Ri / @as(f64, @floatFromInt(n_slices));
    var sasa: f64 = 0.0;

    // Iterate over slices
    var slice_idx: u32 = 0;
    while (slice_idx < n_slices) : (slice_idx += 1) {
        // z-coordinate of this slice (center of slice)
        const slice_z = zi - Ri + delta * (@as(f64, @floatFromInt(slice_idx)) + 0.5);

        // Distance from atom center to slice
        const di = @abs(zi - slice_z);

        // Radius of atom i's cross-section at this slice
        const Ri_prime2 = Ri * Ri - di * di;
        if (Ri_prime2 <= 0) continue; // Round-off protection
        const Ri_prime = @sqrt(Ri_prime2);
        if (Ri_prime <= 0) continue;

        // Find buried arcs from neighbors
        var n_arcs: usize = 0;
        var is_buried = false;

        // Process neighbors in batches of 8 using SIMD
        var i: usize = 0;
        while (i + 8 <= neighbors.len and !is_buried) : (i += 8) {
            // Load batch of 8 neighbors
            const batch_x = [8]f64{
                x[neighbors[i]],
                x[neighbors[i + 1]],
                x[neighbors[i + 2]],
                x[neighbors[i + 3]],
                x[neighbors[i + 4]],
                x[neighbors[i + 5]],
                x[neighbors[i + 6]],
                x[neighbors[i + 7]],
            };
            const batch_y = [8]f64{
                y[neighbors[i]],
                y[neighbors[i + 1]],
                y[neighbors[i + 2]],
                y[neighbors[i + 3]],
                y[neighbors[i + 4]],
                y[neighbors[i + 5]],
                y[neighbors[i + 6]],
                y[neighbors[i + 7]],
            };
            const batch_z = [8]f64{
                z[neighbors[i]],
                z[neighbors[i + 1]],
                z[neighbors[i + 2]],
                z[neighbors[i + 3]],
                z[neighbors[i + 4]],
                z[neighbors[i + 5]],
                z[neighbors[i + 6]],
                z[neighbors[i + 7]],
            };
            const batch_r = [8]f64{
                radii[neighbors[i]],
                radii[neighbors[i + 1]],
                radii[neighbors[i + 2]],
                radii[neighbors[i + 3]],
                radii[neighbors[i + 4]],
                radii[neighbors[i + 5]],
                radii[neighbors[i + 6]],
                radii[neighbors[i + 7]],
            };

            // SIMD: Calculate slice radii and xy-distances
            const rj_primes = simd.sliceRadiiBatch8(slice_z, batch_z, batch_r);
            const dij_batch = simd.xyDistanceBatch8(xi, yi, batch_x, batch_y);

            // SIMD: Check which circles overlap
            const overlap_mask = simd.circlesOverlapBatch8(dij_batch, Ri_prime, rj_primes);

            // Process overlapping neighbors
            for (0..8) |k| {
                if ((overlap_mask >> @intCast(k)) & 1 == 0) continue;

                const Rj_prime = rj_primes[k];
                if (Rj_prime <= 0) continue; // No slice intersection

                const dij = dij_batch[k];
                const dx = batch_x[k] - xi;
                const dy = batch_y[k] - yi;
                const Rj_prime2 = Rj_prime * Rj_prime;

                // Handle near-zero distance
                if (dij < 1e-10) {
                    if (Rj_prime > Ri_prime) {
                        is_buried = true;
                        break;
                    }
                    continue;
                }

                // Check if circle i is completely inside circle j (or touches it from inside)
                if (dij + Ri_prime <= Rj_prime) {
                    is_buried = true;
                    break;
                }
                // Check if circle j is completely inside circle i (or touches it from inside)
                if (dij + Rj_prime <= Ri_prime) {
                    continue;
                }

                // Calculate arc
                const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
                // Tangent circles that rounding let past the tests above (see `Arc`)
                if (cos_alpha <= -1.0) {
                    is_buried = true;
                    break;
                }
                if (cos_alpha >= 1.0) continue;
                const alpha = batchHalfAngle(trig, cos_alpha);
                const beta = batchDirection(trig, dy, dx) + std.math.pi;

                var inf = beta - alpha;
                var sup = beta + alpha;
                // An arc without width covers nothing (see `Arc`)
                if (!(inf < sup)) continue;

                while (inf < 0) inf += TWOPI;
                while (inf >= TWOPI) inf -= TWOPI;
                while (sup <= 0) sup += TWOPI;
                while (sup > TWOPI) sup -= TWOPI;

                // Bounds check before writing (defensive)
                if (n_arcs + 2 > arc_buffer.len) break;

                if (sup < inf) {
                    arc_buffer[n_arcs] = Arc{ .start = 0, .end = sup };
                    n_arcs += 1;
                    arc_buffer[n_arcs] = Arc{ .start = inf, .end = TWOPI };
                    n_arcs += 1;
                } else {
                    arc_buffer[n_arcs] = Arc{ .start = inf, .end = sup };
                    n_arcs += 1;
                }
            }
        }

        // Process remaining neighbors in batches of 4 using SIMD
        while (i + 4 <= neighbors.len and !is_buried) : (i += 4) {
            // Load batch of 4 neighbors
            const batch_x = [4]f64{
                x[neighbors[i]],
                x[neighbors[i + 1]],
                x[neighbors[i + 2]],
                x[neighbors[i + 3]],
            };
            const batch_y = [4]f64{
                y[neighbors[i]],
                y[neighbors[i + 1]],
                y[neighbors[i + 2]],
                y[neighbors[i + 3]],
            };
            const batch_z = [4]f64{
                z[neighbors[i]],
                z[neighbors[i + 1]],
                z[neighbors[i + 2]],
                z[neighbors[i + 3]],
            };
            const batch_r = [4]f64{
                radii[neighbors[i]],
                radii[neighbors[i + 1]],
                radii[neighbors[i + 2]],
                radii[neighbors[i + 3]],
            };

            // SIMD: Calculate slice radii and xy-distances
            const rj_primes = simd.sliceRadiiBatch4(slice_z, batch_z, batch_r);
            const dij_batch = simd.xyDistanceBatch4(xi, yi, batch_x, batch_y);

            // SIMD: Check which circles overlap
            const overlap_mask = simd.circlesOverlapBatch4(dij_batch, Ri_prime, rj_primes);

            // Process overlapping neighbors
            for (0..4) |k| {
                if ((overlap_mask >> @intCast(k)) & 1 == 0) continue;

                const Rj_prime = rj_primes[k];
                if (Rj_prime <= 0) continue; // No slice intersection

                const dij = dij_batch[k];
                const dx = batch_x[k] - xi;
                const dy = batch_y[k] - yi;
                const Rj_prime2 = Rj_prime * Rj_prime;

                // Handle near-zero distance
                if (dij < 1e-10) {
                    if (Rj_prime > Ri_prime) {
                        is_buried = true;
                        break;
                    }
                    continue;
                }

                // Check if circle i is completely inside circle j (or touches it from inside)
                if (dij + Ri_prime <= Rj_prime) {
                    is_buried = true;
                    break;
                }
                // Check if circle j is completely inside circle i (or touches it from inside)
                if (dij + Rj_prime <= Ri_prime) {
                    continue;
                }

                // Calculate arc
                const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
                // Tangent circles that rounding let past the tests above (see `Arc`)
                if (cos_alpha <= -1.0) {
                    is_buried = true;
                    break;
                }
                if (cos_alpha >= 1.0) continue;
                const alpha = batchHalfAngle(trig, cos_alpha);
                const beta = batchDirection(trig, dy, dx) + std.math.pi;

                var inf = beta - alpha;
                var sup = beta + alpha;
                // An arc without width covers nothing (see `Arc`)
                if (!(inf < sup)) continue;

                while (inf < 0) inf += TWOPI;
                while (inf >= TWOPI) inf -= TWOPI;
                while (sup <= 0) sup += TWOPI;
                while (sup > TWOPI) sup -= TWOPI;

                // Bounds check before writing (defensive)
                if (n_arcs + 2 > arc_buffer.len) break;

                if (sup < inf) {
                    arc_buffer[n_arcs] = Arc{ .start = 0, .end = sup };
                    n_arcs += 1;
                    arc_buffer[n_arcs] = Arc{ .start = inf, .end = TWOPI };
                    n_arcs += 1;
                } else {
                    arc_buffer[n_arcs] = Arc{ .start = inf, .end = sup };
                    n_arcs += 1;
                }
            }
        }

        // Process remaining neighbors (scalar, exact trigonometry in both modes)
        while (i < neighbors.len and !is_buried) : (i += 1) {
            const j = neighbors[i];
            const zj = z[j];
            const dj = @abs(zj - slice_z);
            const Rj = radii[j];

            if (dj >= Rj) continue;

            const Rj_prime2 = Rj * Rj - dj * dj;
            const Rj_prime = @sqrt(Rj_prime2);

            const dx = x[j] - xi;
            const dy = y[j] - yi;
            const dij = @sqrt(dx * dx + dy * dy);

            if (dij >= Ri_prime + Rj_prime) continue;

            if (dij < 1e-10) {
                if (Rj_prime > Ri_prime) {
                    is_buried = true;
                    break;
                }
                continue;
            }

            if (dij + Ri_prime <= Rj_prime) {
                is_buried = true;
                break;
            }
            if (dij + Rj_prime <= Ri_prime) continue;

            const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
            // Tangent circles that rounding let past the tests above (see `Arc`)
            if (cos_alpha <= -1.0) {
                is_buried = true;
                break;
            }
            if (cos_alpha >= 1.0) continue;
            const alpha = std.math.acos(std.math.clamp(cos_alpha, -1.0, 1.0));
            const beta = std.math.atan2(dy, dx) + std.math.pi;

            var inf = beta - alpha;
            var sup = beta + alpha;
            // An arc without width covers nothing (see `Arc`)
            if (!(inf < sup)) continue;

            while (inf < 0) inf += TWOPI;
            while (inf >= TWOPI) inf -= TWOPI;
            while (sup <= 0) sup += TWOPI;
            while (sup > TWOPI) sup -= TWOPI;

            // Bounds check before writing (defensive)
            if (n_arcs + 2 > arc_buffer.len) break;

            if (sup < inf) {
                arc_buffer[n_arcs] = Arc{ .start = 0, .end = sup };
                n_arcs += 1;
                arc_buffer[n_arcs] = Arc{ .start = inf, .end = TWOPI };
                n_arcs += 1;
            } else {
                arc_buffer[n_arcs] = Arc{ .start = inf, .end = sup };
                n_arcs += 1;
            }
        }

        if (!is_buried) {
            const exposed = exposedArcLength(arc_buffer[0..n_arcs]);
            sasa += delta * Ri * exposed;
        }
    }

    return sasa;
}

/// Calculate total exposed arc length from list of buried arcs
fn exposedArcLength(arcs: []Arc) f64 {
    if (arcs.len == 0) return TWOPI;

    // Sort arcs by start angle (insertion sort - efficient for small arrays)
    sortArcs(arcs);

    // Merge overlapping arcs and calculate exposed length
    var sum: f64 = arcs[0].start; // Exposed before first arc
    var sup: f64 = arcs[0].end; // Current coverage end

    for (arcs[1..]) |arc| {
        if (sup < arc.start) {
            // Gap between arcs - add exposed portion
            sum += arc.start - sup;
        }
        if (arc.end > sup) {
            sup = arc.end;
        }
    }

    // Add exposed portion after last arc
    sum += TWOPI - sup;

    return sum;
}

/// Sort arcs by start angle using insertion sort
fn sortArcs(arcs: []Arc) void {
    if (arcs.len <= 1) return;

    for (1..arcs.len) |i| {
        const tmp = arcs[i];
        var j = i;
        while (j > 0 and arcs[j - 1].start > tmp.start) {
            arcs[j] = arcs[j - 1];
            j -= 1;
        }
        arcs[j] = tmp;
    }
}

/// Context for parallel Lee-Richards calculation workers.
/// Thread safety: All fields are read-only except `atom_areas` and
/// `arc_buffers`, which have disjoint write access (each chunk writes to its
/// own atom indices and to its own slice of `arc_buffers`).
const ParallelContext = struct {
    x: []const f64,
    y: []const f64,
    z: []const f64,
    radii: []const f64,
    neighbor_list: *const NeighborList,
    n_slices: u32,
    trig: TrigMode,
    max_arc_buffer_size: usize,
    /// Scratch space for all chunks: `max_arc_buffer_size` arcs per chunk.
    arc_buffers: []Arc,
    /// Chunk size handed to the thread pool, used to find a chunk's slice.
    chunk_size: usize,
    atom_areas: []f64,
};

/// Allocate the arc scratch space of a parallel run: one buffer per chunk the
/// thread pool will hand out. Allocating it up front leaves the workers with
/// nothing to allocate, so they cannot fail part-way through the atoms.
fn allocArcBuffers(
    comptime ArcT: type,
    allocator: Allocator,
    n_atoms: usize,
    chunk_size: usize,
    max_arc_buffer_size: usize,
) Allocator.Error![]ArcT {
    const n_chunks = thread_pool.chunkCount(n_atoms, chunk_size);
    const total = std.math.mul(usize, n_chunks, max_arc_buffer_size) catch return error.OutOfMemory;
    return allocator.alloc(ArcT, total);
}

/// The slice of `arc_buffers` that belongs to the chunk starting at `chunk_start`.
fn chunkArcBuffer(
    comptime ArcT: type,
    arc_buffers: []ArcT,
    chunk_start: usize,
    chunk_size: usize,
    max_arc_buffer_size: usize,
) []ArcT {
    const offset = (chunk_start / chunk_size) * max_arc_buffer_size;
    return arc_buffers[offset..][0..max_arc_buffer_size];
}

/// Worker function for parallel Lee-Richards calculation.
/// Processes atoms from chunk_start to chunk_end.
fn parallelLeeRichardsWorker(ctx: ParallelContext, chunk_start: usize, chunk_end: usize) f64 {
    const arc_buffer = chunkArcBuffer(Arc, ctx.arc_buffers, chunk_start, ctx.chunk_size, ctx.max_arc_buffer_size);

    var chunk_total: f64 = 0.0;

    for (chunk_start..chunk_end) |i| {
        const area = atomArea(
            i,
            ctx.x,
            ctx.y,
            ctx.z,
            ctx.radii,
            ctx.neighbor_list,
            ctx.n_slices,
            ctx.trig,
            arc_buffer,
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

/// Calculate SASA using Lee-Richards algorithm with parallel processing.
///
/// # Parameters
/// - `allocator`: Memory allocator for result arrays
/// - `input`: Atom input data (positions and radii)
/// - `config`: Configuration parameters (n_slices, probe_radius, trig)
/// - `n_threads`: Number of worker threads (0 = auto-detect)
///
/// # Returns
/// SasaResult containing total_area and per-atom areas. Caller must call deinit().
pub fn calculateSasaParallel(
    allocator: Allocator,
    input: AtomInput,
    config: LeeRichardsConfig,
    n_threads: usize,
) !SasaResult {
    const n_atoms = input.atomCount();
    if (n_atoms == 0) {
        return SasaResult{
            .total_area = 0.0,
            .atom_areas = try allocator.alloc(f64, 0),
            .allocator = allocator,
        };
    }
    if (config.n_slices == 0) return error.InvalidInput;
    try types.validateProbeRadius(f64, config.probe_radius);
    try input.validateFiniteAndRange();

    // Auto-detect thread count if 0
    const actual_threads = if (n_threads == 0)
        try std.Thread.getCpuCount()
    else
        n_threads;

    // Convert to Vec3 positions for neighbor list
    const positions = try allocator.alloc(Vec3, n_atoms);
    defer allocator.free(positions);
    for (0..n_atoms) |i| {
        positions[i] = Vec3{ .x = input.x[i], .y = input.y[i], .z = input.z[i] };
    }

    // Pre-compute effective radii (atom radius + probe radius)
    const radii = try allocator.alloc(f64, n_atoms);
    defer allocator.free(radii);
    for (0..n_atoms) |i| {
        radii[i] = input.r[i] + config.probe_radius;
    }

    // Build neighbor list with effective radii
    // Note: pass probe_radius=0 since we already added it to radii
    var neighbor_list = try NeighborList.init(allocator, positions, radii, 0.0);
    defer neighbor_list.deinit();

    // Estimate max neighbors for arc buffer allocation
    var max_neighbors: usize = 0;
    for (0..n_atoms) |i| {
        max_neighbors = @max(max_neighbors, neighbor_list.getNeighbors(i).len);
    }
    // Each neighbor can create up to 2 arcs (when crossing 0)
    const max_arc_buffer_size = (max_neighbors + 1) * 2;

    // Chunk size heuristic:
    // - Minimum 64 atoms per chunk to amortize thread overhead
    // - Target 4 chunks per thread for load balancing
    const chunk_size = @max(64, n_atoms / (actual_threads * 4));

    // Allocate the arc buffers of all chunks before any worker starts
    const arc_buffers = try allocArcBuffers(Arc, allocator, n_atoms, chunk_size, max_arc_buffer_size);
    defer allocator.free(arc_buffers);

    // Allocate result arrays
    const atom_areas = try allocator.alloc(f64, n_atoms);
    errdefer allocator.free(atom_areas);

    // Create parallel context
    const ctx = ParallelContext{
        .x = input.x,
        .y = input.y,
        .z = input.z,
        .radii = radii,
        .neighbor_list = &neighbor_list,
        .n_slices = config.n_slices,
        .trig = config.trig,
        .max_arc_buffer_size = max_arc_buffer_size,
        .arc_buffers = arc_buffers,
        .chunk_size = chunk_size,
        .atom_areas = atom_areas,
    };

    // Run parallel calculation
    const total_area = try thread_pool.parallelFor(
        ParallelContext,
        f64,
        allocator,
        actual_threads,
        parallelLeeRichardsWorker,
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

/// Generic configuration for Lee-Richards algorithm
pub fn LeeRichardsConfigGen(comptime T: type) type {
    return struct {
        /// Number of slices per atom diameter
        n_slices: u32 = 20,
        /// Water probe radius in Angstroms
        probe_radius: T = 1.4,
        /// Trigonometry used for the arc angles
        trig: TrigMode = .exact,
    };
}

/// Generic Lee-Richards algorithm implementation supporting both f32 and f64 precision.
pub fn LeeRichardsGen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    const Result = SasaResultGen(T);
    const Cfg = LeeRichardsConfigGen(T);
    const NList = NeighborListGen(T);

    return struct {
        const Self = @This();

        const FastAcos = simd.fastAcosGen(T);
        const FastAtan2 = simd.fastAtan2Gen(T);
        const SliceRadiiBatch8 = simd.sliceRadiiBatch8Gen(T);
        const SliceRadiiBatch4 = simd.sliceRadiiBatch4Gen(T);
        const XyDistanceBatch8 = simd.xyDistanceBatch8Gen(T);
        const XyDistanceBatch4 = simd.xyDistanceBatch4Gen(T);
        const CirclesOverlapBatch8 = simd.circlesOverlapBatch8Gen(T);
        const CirclesOverlapBatch4 = simd.circlesOverlapBatch4Gen(T);

        const TWOPI_T: T = 2.0 * std.math.pi;

        /// Arc interval representing a buried portion of a circle.
        /// See the non-generic `Arc` for the handling of tangent circles.
        pub const Arc = struct {
            start: T, // Start angle (radians)
            end: T, // End angle (radians)
        };

        /// Sort arcs by start angle using insertion sort
        fn sortArcs(arcs: []Self.Arc) void {
            if (arcs.len <= 1) return;

            for (1..arcs.len) |i| {
                const tmp = arcs[i];
                var j = i;
                while (j > 0 and arcs[j - 1].start > tmp.start) {
                    arcs[j] = arcs[j - 1];
                    j -= 1;
                }
                arcs[j] = tmp;
            }
        }

        /// Calculate total exposed arc length from list of buried arcs
        fn exposedArcLength(arcs: []Self.Arc) T {
            if (arcs.len == 0) return Self.TWOPI_T;

            // Sort arcs by start angle
            Self.sortArcs(arcs);

            // Merge overlapping arcs and calculate exposed length
            var sum: T = arcs[0].start; // Exposed before first arc
            var sup: T = arcs[0].end; // Current coverage end

            for (arcs[1..]) |arc| {
                if (sup < arc.start) {
                    // Gap between arcs - add exposed portion
                    sum += arc.start - sup;
                }
                if (arc.end > sup) {
                    sup = arc.end;
                }
            }

            // Add exposed portion after last arc
            sum += Self.TWOPI_T - sup;

            return sum;
        }

        /// Half-angle of the arc that neighbor j covers on the slice circle of
        /// atom i, for a neighbor of the 8-wide and 4-wide batches. The
        /// neighbors left over after the batches use `std.math.acos` in both
        /// modes.
        inline fn batchHalfAngle(trig: TrigMode, cos_alpha: T) T {
            return switch (trig) {
                .exact => std.math.acos(std.math.clamp(cos_alpha, @as(T, -1.0), @as(T, 1.0))),
                .fast => FastAcos.compute(cos_alpha),
            };
        }

        /// Direction of neighbor j seen from atom i in the slice plane, for a
        /// neighbor of the 8-wide and 4-wide batches. The neighbors left over
        /// after the batches use `std.math.atan2` in both modes.
        inline fn batchDirection(trig: TrigMode, dy: T, dx: T) T {
            return switch (trig) {
                .exact => std.math.atan2(dy, dx),
                .fast => FastAtan2.compute(dy, dx),
            };
        }

        /// Calculate SASA for a single atom using slice-based method with SIMD optimization
        fn atomArea(
            atom_idx: usize,
            x: []const T,
            y: []const T,
            z: []const T,
            radii: []const T,
            neighbor_list: *const NList,
            n_slices: u32,
            trig: TrigMode,
            arc_buffer: []Self.Arc,
        ) T {
            const xi = x[atom_idx];
            const yi = y[atom_idx];
            const zi = z[atom_idx];
            const Ri = radii[atom_idx];

            const neighbors = neighbor_list.getNeighbors(atom_idx);
            if (neighbors.len == 0) {
                // No neighbors, full sphere exposed
                return 4.0 * std.math.pi * Ri * Ri;
            }

            const delta = 2.0 * Ri / @as(T, @floatFromInt(n_slices));
            var sasa: T = 0.0;

            // Iterate over slices
            var slice_idx: u32 = 0;
            while (slice_idx < n_slices) : (slice_idx += 1) {
                // z-coordinate of this slice (center of slice)
                const slice_z = zi - Ri + delta * (@as(T, @floatFromInt(slice_idx)) + 0.5);

                // Distance from atom center to slice
                const di = @abs(zi - slice_z);

                // Radius of atom i's cross-section at this slice
                const Ri_prime2 = Ri * Ri - di * di;
                if (Ri_prime2 <= 0) continue; // Round-off protection
                const Ri_prime = @sqrt(Ri_prime2);
                if (Ri_prime <= 0) continue;

                // Find buried arcs from neighbors
                var n_arcs: usize = 0;
                var is_buried = false;

                // Process neighbors in batches of 8 using SIMD
                var i: usize = 0;
                while (i + 8 <= neighbors.len and !is_buried) : (i += 8) {
                    // Load batch of 8 neighbors
                    const batch_x = [8]T{
                        x[neighbors[i]],
                        x[neighbors[i + 1]],
                        x[neighbors[i + 2]],
                        x[neighbors[i + 3]],
                        x[neighbors[i + 4]],
                        x[neighbors[i + 5]],
                        x[neighbors[i + 6]],
                        x[neighbors[i + 7]],
                    };
                    const batch_y = [8]T{
                        y[neighbors[i]],
                        y[neighbors[i + 1]],
                        y[neighbors[i + 2]],
                        y[neighbors[i + 3]],
                        y[neighbors[i + 4]],
                        y[neighbors[i + 5]],
                        y[neighbors[i + 6]],
                        y[neighbors[i + 7]],
                    };
                    const batch_z = [8]T{
                        z[neighbors[i]],
                        z[neighbors[i + 1]],
                        z[neighbors[i + 2]],
                        z[neighbors[i + 3]],
                        z[neighbors[i + 4]],
                        z[neighbors[i + 5]],
                        z[neighbors[i + 6]],
                        z[neighbors[i + 7]],
                    };
                    const batch_r = [8]T{
                        radii[neighbors[i]],
                        radii[neighbors[i + 1]],
                        radii[neighbors[i + 2]],
                        radii[neighbors[i + 3]],
                        radii[neighbors[i + 4]],
                        radii[neighbors[i + 5]],
                        radii[neighbors[i + 6]],
                        radii[neighbors[i + 7]],
                    };

                    // SIMD: Calculate slice radii and xy-distances
                    const rj_primes = SliceRadiiBatch8.compute(slice_z, batch_z, batch_r);
                    const dij_batch = XyDistanceBatch8.compute(xi, yi, batch_x, batch_y);

                    // SIMD: Check which circles overlap
                    const overlap_mask = CirclesOverlapBatch8.compute(dij_batch, Ri_prime, rj_primes);

                    // Process overlapping neighbors
                    for (0..8) |k| {
                        if ((overlap_mask >> @intCast(k)) & 1 == 0) continue;

                        const Rj_prime = rj_primes[k];
                        if (Rj_prime <= 0) continue; // No slice intersection

                        const dij = dij_batch[k];
                        const dx = batch_x[k] - xi;
                        const dy = batch_y[k] - yi;
                        const Rj_prime2 = Rj_prime * Rj_prime;

                        // Handle near-zero distance
                        const epsilon = types.Epsilon(T).default;
                        if (dij < epsilon) {
                            if (Rj_prime > Ri_prime) {
                                is_buried = true;
                                break;
                            }
                            continue;
                        }

                        // Check if circle i is completely inside circle j (or touches it from inside)
                        if (dij + Ri_prime <= Rj_prime) {
                            is_buried = true;
                            break;
                        }
                        // Check if circle j is completely inside circle i (or touches it from inside)
                        if (dij + Rj_prime <= Ri_prime) {
                            continue;
                        }

                        // Calculate arc
                        const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
                        // Tangent circles that rounding let past the tests above (see `Arc`)
                        if (cos_alpha <= -1.0) {
                            is_buried = true;
                            break;
                        }
                        if (cos_alpha >= 1.0) continue;
                        const alpha = Self.batchHalfAngle(trig, cos_alpha);
                        const beta = Self.batchDirection(trig, dy, dx) + std.math.pi;

                        var inf = beta - alpha;
                        var sup = beta + alpha;
                        // An arc without width covers nothing (see `Arc`)
                        if (!(inf < sup)) continue;

                        while (inf < 0) inf += Self.TWOPI_T;
                        while (inf >= Self.TWOPI_T) inf -= Self.TWOPI_T;
                        while (sup <= 0) sup += Self.TWOPI_T;
                        while (sup > Self.TWOPI_T) sup -= Self.TWOPI_T;

                        // Bounds check before writing (defensive)
                        if (n_arcs + 2 > arc_buffer.len) break;

                        if (sup < inf) {
                            arc_buffer[n_arcs] = Self.Arc{ .start = 0, .end = sup };
                            n_arcs += 1;
                            arc_buffer[n_arcs] = Self.Arc{ .start = inf, .end = Self.TWOPI_T };
                            n_arcs += 1;
                        } else {
                            arc_buffer[n_arcs] = Self.Arc{ .start = inf, .end = sup };
                            n_arcs += 1;
                        }
                    }
                }

                // Process remaining neighbors in batches of 4 using SIMD
                while (i + 4 <= neighbors.len and !is_buried) : (i += 4) {
                    const batch_x = [4]T{
                        x[neighbors[i]],
                        x[neighbors[i + 1]],
                        x[neighbors[i + 2]],
                        x[neighbors[i + 3]],
                    };
                    const batch_y = [4]T{
                        y[neighbors[i]],
                        y[neighbors[i + 1]],
                        y[neighbors[i + 2]],
                        y[neighbors[i + 3]],
                    };
                    const batch_z = [4]T{
                        z[neighbors[i]],
                        z[neighbors[i + 1]],
                        z[neighbors[i + 2]],
                        z[neighbors[i + 3]],
                    };
                    const batch_r = [4]T{
                        radii[neighbors[i]],
                        radii[neighbors[i + 1]],
                        radii[neighbors[i + 2]],
                        radii[neighbors[i + 3]],
                    };

                    const rj_primes = SliceRadiiBatch4.compute(slice_z, batch_z, batch_r);
                    const dij_batch = XyDistanceBatch4.compute(xi, yi, batch_x, batch_y);
                    const overlap_mask = CirclesOverlapBatch4.compute(dij_batch, Ri_prime, rj_primes);

                    for (0..4) |k| {
                        if ((overlap_mask >> @intCast(k)) & 1 == 0) continue;

                        const Rj_prime = rj_primes[k];
                        if (Rj_prime <= 0) continue;

                        const dij = dij_batch[k];
                        const dx = batch_x[k] - xi;
                        const dy = batch_y[k] - yi;
                        const Rj_prime2 = Rj_prime * Rj_prime;

                        const epsilon = types.Epsilon(T).default;
                        if (dij < epsilon) {
                            if (Rj_prime > Ri_prime) {
                                is_buried = true;
                                break;
                            }
                            continue;
                        }

                        if (dij + Ri_prime <= Rj_prime) {
                            is_buried = true;
                            break;
                        }
                        if (dij + Rj_prime <= Ri_prime) continue;

                        const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
                        // Tangent circles that rounding let past the tests above (see `Arc`)
                        if (cos_alpha <= -1.0) {
                            is_buried = true;
                            break;
                        }
                        if (cos_alpha >= 1.0) continue;
                        const alpha = Self.batchHalfAngle(trig, cos_alpha);
                        const beta = Self.batchDirection(trig, dy, dx) + std.math.pi;

                        var inf = beta - alpha;
                        var sup = beta + alpha;
                        // An arc without width covers nothing (see `Arc`)
                        if (!(inf < sup)) continue;

                        while (inf < 0) inf += Self.TWOPI_T;
                        while (inf >= Self.TWOPI_T) inf -= Self.TWOPI_T;
                        while (sup <= 0) sup += Self.TWOPI_T;
                        while (sup > Self.TWOPI_T) sup -= Self.TWOPI_T;

                        // Bounds check before writing (defensive)
                        if (n_arcs + 2 > arc_buffer.len) break;

                        if (sup < inf) {
                            arc_buffer[n_arcs] = Self.Arc{ .start = 0, .end = sup };
                            n_arcs += 1;
                            arc_buffer[n_arcs] = Self.Arc{ .start = inf, .end = Self.TWOPI_T };
                            n_arcs += 1;
                        } else {
                            arc_buffer[n_arcs] = Self.Arc{ .start = inf, .end = sup };
                            n_arcs += 1;
                        }
                    }
                }

                // Process remaining neighbors (scalar, exact trigonometry in both modes)
                while (i < neighbors.len and !is_buried) : (i += 1) {
                    const j = neighbors[i];
                    const zj = z[j];
                    const dj = @abs(zj - slice_z);
                    const Rj = radii[j];

                    if (dj >= Rj) continue;

                    const Rj_prime2 = Rj * Rj - dj * dj;
                    const Rj_prime = @sqrt(Rj_prime2);

                    const dx = x[j] - xi;
                    const dy = y[j] - yi;
                    const dij = @sqrt(dx * dx + dy * dy);

                    if (dij >= Ri_prime + Rj_prime) continue;

                    const epsilon = types.Epsilon(T).default;
                    if (dij < epsilon) {
                        if (Rj_prime > Ri_prime) {
                            is_buried = true;
                            break;
                        }
                        continue;
                    }

                    if (dij + Ri_prime <= Rj_prime) {
                        is_buried = true;
                        break;
                    }
                    if (dij + Rj_prime <= Ri_prime) continue;

                    const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
                    // Tangent circles that rounding let past the tests above (see `Arc`)
                    if (cos_alpha <= -1.0) {
                        is_buried = true;
                        break;
                    }
                    if (cos_alpha >= 1.0) continue;
                    const alpha = std.math.acos(std.math.clamp(cos_alpha, @as(T, -1.0), @as(T, 1.0)));
                    const beta = std.math.atan2(dy, dx) + std.math.pi;

                    var inf = beta - alpha;
                    var sup = beta + alpha;
                    // An arc without width covers nothing (see `Arc`)
                    if (!(inf < sup)) continue;

                    while (inf < 0) inf += Self.TWOPI_T;
                    while (inf >= Self.TWOPI_T) inf -= Self.TWOPI_T;
                    while (sup <= 0) sup += Self.TWOPI_T;
                    while (sup > Self.TWOPI_T) sup -= Self.TWOPI_T;

                    // Bounds check before writing (defensive)
                    if (n_arcs + 2 > arc_buffer.len) break;

                    if (sup < inf) {
                        arc_buffer[n_arcs] = Self.Arc{ .start = 0, .end = sup };
                        n_arcs += 1;
                        arc_buffer[n_arcs] = Self.Arc{ .start = inf, .end = Self.TWOPI_T };
                        n_arcs += 1;
                    } else {
                        arc_buffer[n_arcs] = Self.Arc{ .start = inf, .end = sup };
                        n_arcs += 1;
                    }
                }

                if (!is_buried) {
                    const exposed = Self.exposedArcLength(arc_buffer[0..n_arcs]);
                    sasa += delta * Ri * exposed;
                }
            }

            return sasa;
        }

        /// Context for parallel Lee-Richards calculation workers.
        pub const ParallelContext = struct {
            x: []const T,
            y: []const T,
            z: []const T,
            radii: []const T,
            neighbor_list: *const NList,
            n_slices: u32,
            trig: TrigMode,
            max_arc_buffer_size: usize,
            /// Scratch space for all chunks: `max_arc_buffer_size` arcs per chunk.
            arc_buffers: []Self.Arc,
            /// Chunk size handed to the thread pool, used to find a chunk's slice.
            chunk_size: usize,
            atom_areas: []T,
        };

        /// Worker function for parallel Lee-Richards calculation.
        fn parallelLeeRichardsWorker(ctx: Self.ParallelContext, chunk_start: usize, chunk_end: usize) T {
            const arc_buffer = chunkArcBuffer(Self.Arc, ctx.arc_buffers, chunk_start, ctx.chunk_size, ctx.max_arc_buffer_size);

            var chunk_total: T = 0.0;

            for (chunk_start..chunk_end) |i| {
                const area = Self.atomArea(
                    i,
                    ctx.x,
                    ctx.y,
                    ctx.z,
                    ctx.radii,
                    ctx.neighbor_list,
                    ctx.n_slices,
                    ctx.trig,
                    arc_buffer,
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

        /// Calculate SASA using Lee-Richards algorithm (single-threaded)
        pub fn calculateSasa(
            allocator: Allocator,
            input: AtomInput,
            config: Cfg,
        ) !Result {
            const n_atoms = input.atomCount();
            if (n_atoms == 0) {
                return Result{
                    .total_area = 0.0,
                    .atom_areas = try allocator.alloc(T, 0),
                    .allocator = allocator,
                };
            }
            if (config.n_slices == 0) return error.InvalidInput;
            try types.validateProbeRadius(T, config.probe_radius);
            try input.validateFiniteAndRange();

            // Convert to Vec positions (cast from f64)
            const positions = try allocator.alloc(Vec, n_atoms);
            defer allocator.free(positions);
            for (0..n_atoms) |i| {
                positions[i] = Vec{
                    .x = @floatCast(input.x[i]),
                    .y = @floatCast(input.y[i]),
                    .z = @floatCast(input.z[i]),
                };
            }

            // Pre-compute effective radii (atom radius + probe radius)
            const radii = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii);
            for (0..n_atoms) |i| {
                radii[i] = @as(T, @floatCast(input.r[i])) + config.probe_radius;
            }

            // Convert x, y, z to T arrays for atomArea
            const x_t = try allocator.alloc(T, n_atoms);
            defer allocator.free(x_t);
            const y_t = try allocator.alloc(T, n_atoms);
            defer allocator.free(y_t);
            const z_t = try allocator.alloc(T, n_atoms);
            defer allocator.free(z_t);
            for (0..n_atoms) |i| {
                x_t[i] = @floatCast(input.x[i]);
                y_t[i] = @floatCast(input.y[i]);
                z_t[i] = @floatCast(input.z[i]);
            }

            // Build neighbor list with effective radii
            var neighbor_list = try NList.init(allocator, positions, radii, 0.0);
            defer neighbor_list.deinit();

            // Calculate SASA for each atom
            const atom_areas = try allocator.alloc(T, n_atoms);
            errdefer allocator.free(atom_areas);
            var total_area: T = 0.0;

            // Estimate max neighbors for arc buffer allocation
            var max_neighbors: usize = 0;
            for (0..n_atoms) |i| {
                max_neighbors = @max(max_neighbors, neighbor_list.getNeighbors(i).len);
            }
            const arc_buffer = try allocator.alloc(Self.Arc, (max_neighbors + 1) * 2);
            defer allocator.free(arc_buffer);

            for (0..n_atoms) |i| {
                atom_areas[i] = Self.atomArea(
                    i,
                    x_t,
                    y_t,
                    z_t,
                    radii,
                    &neighbor_list,
                    config.n_slices,
                    config.trig,
                    arc_buffer,
                );
                total_area += atom_areas[i];
            }

            return Result{
                .total_area = total_area,
                .atom_areas = atom_areas,
                .allocator = allocator,
            };
        }

        /// Calculate SASA using Lee-Richards algorithm with parallel processing.
        pub fn calculateSasaParallel(
            allocator: Allocator,
            input: AtomInput,
            config: Cfg,
            n_threads: usize,
        ) !Result {
            const n_atoms = input.atomCount();
            if (n_atoms == 0) {
                return Result{
                    .total_area = 0.0,
                    .atom_areas = try allocator.alloc(T, 0),
                    .allocator = allocator,
                };
            }
            if (config.n_slices == 0) return error.InvalidInput;
            try types.validateProbeRadius(T, config.probe_radius);
            try input.validateFiniteAndRange();

            // Auto-detect thread count if 0
            const actual_threads = if (n_threads == 0)
                try std.Thread.getCpuCount()
            else
                n_threads;

            // Convert to Vec positions (cast from f64)
            const positions = try allocator.alloc(Vec, n_atoms);
            defer allocator.free(positions);
            for (0..n_atoms) |i| {
                positions[i] = Vec{
                    .x = @floatCast(input.x[i]),
                    .y = @floatCast(input.y[i]),
                    .z = @floatCast(input.z[i]),
                };
            }

            // Pre-compute effective radii (atom radius + probe radius)
            const radii = try allocator.alloc(T, n_atoms);
            defer allocator.free(radii);
            for (0..n_atoms) |i| {
                radii[i] = @as(T, @floatCast(input.r[i])) + config.probe_radius;
            }

            // Convert x, y, z to T arrays
            const x_t = try allocator.alloc(T, n_atoms);
            defer allocator.free(x_t);
            const y_t = try allocator.alloc(T, n_atoms);
            defer allocator.free(y_t);
            const z_t = try allocator.alloc(T, n_atoms);
            defer allocator.free(z_t);
            for (0..n_atoms) |i| {
                x_t[i] = @floatCast(input.x[i]);
                y_t[i] = @floatCast(input.y[i]);
                z_t[i] = @floatCast(input.z[i]);
            }

            // Build neighbor list with effective radii
            var neighbor_list = try NList.init(allocator, positions, radii, 0.0);
            defer neighbor_list.deinit();

            // Estimate max neighbors for arc buffer allocation
            var max_neighbors: usize = 0;
            for (0..n_atoms) |i| {
                max_neighbors = @max(max_neighbors, neighbor_list.getNeighbors(i).len);
            }
            const max_arc_buffer_size = (max_neighbors + 1) * 2;

            // Chunk size heuristic
            const chunk_size = @max(64, n_atoms / (actual_threads * 4));

            // Allocate the arc buffers of all chunks before any worker starts
            const arc_buffers = try allocArcBuffers(Self.Arc, allocator, n_atoms, chunk_size, max_arc_buffer_size);
            defer allocator.free(arc_buffers);

            // Allocate result arrays
            const atom_areas = try allocator.alloc(T, n_atoms);
            errdefer allocator.free(atom_areas);

            // Create parallel context
            const ctx = Self.ParallelContext{
                .x = x_t,
                .y = y_t,
                .z = z_t,
                .radii = radii,
                .neighbor_list = &neighbor_list,
                .n_slices = config.n_slices,
                .trig = config.trig,
                .max_arc_buffer_size = max_arc_buffer_size,
                .arc_buffers = arc_buffers,
                .chunk_size = chunk_size,
                .atom_areas = atom_areas,
            };

            // Run parallel calculation
            const total_area = try thread_pool.parallelFor(
                Self.ParallelContext,
                T,
                allocator,
                actual_threads,
                Self.parallelLeeRichardsWorker,
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

/// f32 precision Lee-Richards implementation
pub const LeeRichardsf32 = LeeRichardsGen(f32);
/// f32 precision Lee-Richards configuration
pub const LeeRichardsConfigf32 = LeeRichardsConfigGen(f32);

/// Convenience function: Calculate SASA using Lee-Richards with f32 precision (single-threaded)
pub fn calculateSasaf32(allocator: Allocator, input: AtomInput, config: LeeRichardsConfigGen(f32)) !SasaResultGen(f32) {
    return LeeRichardsf32.calculateSasa(allocator, input, config);
}

/// Convenience function: Calculate SASA using Lee-Richards with f32 precision (parallel)
pub fn calculateSasaParallelf32(allocator: Allocator, input: AtomInput, config: LeeRichardsConfigGen(f32), n_threads: usize) !SasaResultGen(f32) {
    return LeeRichardsf32.calculateSasaParallel(allocator, input, config, n_threads);
}

// Tests
test "exposedArcLength empty" {
    var arcs: [0]Arc = .{};
    const result = exposedArcLength(&arcs);
    try std.testing.expectApproxEqAbs(TWOPI, result, 1e-10);
}

test "exposedArcLength full coverage" {
    var arcs = [_]Arc{
        Arc{ .start = 0, .end = TWOPI },
    };
    const result = exposedArcLength(&arcs);
    try std.testing.expectApproxEqAbs(0.0, result, 1e-10);
}

test "exposedArcLength partial coverage" {
    // Two arcs covering 20% total (10% each)
    var arcs = [_]Arc{
        Arc{ .start = 0, .end = 0.1 * TWOPI },
        Arc{ .start = 0.9 * TWOPI, .end = TWOPI },
    };
    const result = exposedArcLength(&arcs);
    try std.testing.expectApproxEqAbs(0.8 * TWOPI, result, 1e-10);
}

test "exposedArcLength overlapping arcs" {
    var arcs = [_]Arc{
        Arc{ .start = 0.1 * TWOPI, .end = 0.3 * TWOPI },
        Arc{ .start = 0.15 * TWOPI, .end = 0.2 * TWOPI }, // Inside first
    };
    const result = exposedArcLength(&arcs);
    try std.testing.expectApproxEqAbs(0.8 * TWOPI, result, 1e-10);
}

test "sortArcs" {
    var arcs = [_]Arc{
        Arc{ .start = 0.5, .end = 0.6 },
        Arc{ .start = 0.1, .end = 0.2 },
        Arc{ .start = 0.3, .end = 0.4 },
    };
    sortArcs(&arcs);

    try std.testing.expectApproxEqAbs(0.1, arcs[0].start, 1e-10);
    try std.testing.expectApproxEqAbs(0.3, arcs[1].start, 1e-10);
    try std.testing.expectApproxEqAbs(0.5, arcs[2].start, 1e-10);
}

test "single atom SASA" {
    const allocator = std.testing.allocator;

    const x_arr = try allocator.alloc(f64, 1);
    defer allocator.free(x_arr);
    const y_arr = try allocator.alloc(f64, 1);
    defer allocator.free(y_arr);
    const z_arr = try allocator.alloc(f64, 1);
    defer allocator.free(z_arr);
    const r_arr = try allocator.alloc(f64, 1);
    defer allocator.free(r_arr);

    x_arr[0] = 0.0;
    y_arr[0] = 0.0;
    z_arr[0] = 0.0;
    r_arr[0] = 1.5;

    const input = AtomInput{
        .x = x_arr,
        .y = y_arr,
        .z = z_arr,
        .r = r_arr,
        .allocator = allocator,
    };

    const config = LeeRichardsConfig{
        .n_slices = 50,
        .probe_radius = 1.4,
    };

    var result = try calculateSasa(allocator, input, config);
    defer result.deinit();

    // Expected: 4π(1.5 + 1.4)² = 4π(2.9)² ≈ 105.68
    const expected = 4.0 * std.math.pi * 2.9 * 2.9;
    try std.testing.expectApproxEqRel(expected, result.total_area, 0.01);
}

test "calculateSasa rejects non-finite inputs" {
    const allocator = std.testing.allocator;
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};

    {
        const x = [_]f64{std.math.nan(f64)};
        const r = [_]f64{1.5};
        const input = AtomInput{
            .x = &x,
            .y = &y,
            .z = &z,
            .r = @constCast(&r),
            .allocator = allocator,
        };
        try std.testing.expectError(error.InvalidInput, calculateSasa(allocator, input, .{}));
    }

    {
        const x = [_]f64{0.0};
        const r = [_]f64{std.math.inf(f64)};
        const input = AtomInput{
            .x = &x,
            .y = &y,
            .z = &z,
            .r = @constCast(&r),
            .allocator = allocator,
        };
        try std.testing.expectError(error.InvalidInput, calculateSasa(allocator, input, .{}));
    }

    {
        const x = [_]f64{0.0};
        const r = [_]f64{1.5};
        const input = AtomInput{
            .x = &x,
            .y = &y,
            .z = &z,
            .r = @constCast(&r),
            .allocator = allocator,
        };
        try std.testing.expectError(
            error.InvalidInput,
            calculateSasa(allocator, input, .{ .probe_radius = std.math.inf(f64) }),
        );
    }
}

test "calculateSasaParallel rejects non-finite inputs" {
    const allocator = std.testing.allocator;
    const x = [_]f64{0.0};
    const y = [_]f64{0.0};
    const z = [_]f64{0.0};
    const r = [_]f64{1.5};
    const input = AtomInput{
        .x = &x,
        .y = &y,
        .z = &z,
        .r = @constCast(&r),
        .allocator = allocator,
    };

    try std.testing.expectError(
        error.InvalidInput,
        calculateSasaParallel(allocator, input, .{ .probe_radius = std.math.nan(f64) }, 2),
    );
}

test "parallel calculation matches serial" {
    const allocator = std.testing.allocator;

    // Create a small multi-atom system for testing
    const n_atoms = 100;
    const x_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x_arr);
    const y_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y_arr);
    const z_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z_arr);
    const r_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r_arr);

    // Create a grid of atoms
    for (0..n_atoms) |i| {
        const fi: f64 = @floatFromInt(i);
        x_arr[i] = @mod(fi, 10.0) * 3.0;
        y_arr[i] = @mod(@floor(fi / 10.0), 10.0) * 3.0;
        z_arr[i] = @floor(fi / 100.0) * 3.0;
        r_arr[i] = 1.5;
    }

    const input = AtomInput{
        .x = x_arr,
        .y = y_arr,
        .z = z_arr,
        .r = r_arr,
        .allocator = allocator,
    };

    const config = LeeRichardsConfig{
        .n_slices = 20,
        .probe_radius = 1.4,
    };

    // Calculate using serial version
    var serial_result = try calculateSasa(allocator, input, config);
    defer serial_result.deinit();

    // Calculate using parallel version (2 threads)
    var parallel_result = try calculateSasaParallel(allocator, input, config, 2);
    defer parallel_result.deinit();

    // Total area should match
    try std.testing.expectApproxEqRel(serial_result.total_area, parallel_result.total_area, 1e-10);

    // Per-atom areas should match
    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqRel(serial_result.atom_areas[i], parallel_result.atom_areas[i], 1e-10);
    }
}

// =============================================================================
// Tests for f32 precision Lee-Richards implementation
// =============================================================================

test "calculateSasaf32 - single atom" {
    const allocator = std.testing.allocator;

    const x_arr = try allocator.alloc(f64, 1);
    defer allocator.free(x_arr);
    const y_arr = try allocator.alloc(f64, 1);
    defer allocator.free(y_arr);
    const z_arr = try allocator.alloc(f64, 1);
    defer allocator.free(z_arr);
    const r_arr = try allocator.alloc(f64, 1);
    defer allocator.free(r_arr);

    x_arr[0] = 0.0;
    y_arr[0] = 0.0;
    z_arr[0] = 0.0;
    r_arr[0] = 1.5;

    const input = AtomInput{
        .x = x_arr,
        .y = y_arr,
        .z = z_arr,
        .r = r_arr,
        .allocator = allocator,
    };

    const config = LeeRichardsConfigGen(f32){
        .n_slices = 50,
        .probe_radius = 1.4,
    };

    var result = try calculateSasaf32(allocator, input, config);
    defer result.deinit();

    // Expected: 4π(1.5 + 1.4)² = 4π(2.9)² ≈ 105.68
    const expected: f32 = 4.0 * std.math.pi * 2.9 * 2.9;
    try std.testing.expectApproxEqRel(expected, result.total_area, 0.01);
}

test "calculateSasaf32 vs f64 - similar results" {
    const allocator = std.testing.allocator;

    const n_atoms = 10;
    const x_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x_arr);
    const y_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y_arr);
    const z_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z_arr);
    const r_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r_arr);

    // Create a cluster of atoms
    x_arr[0] = 0.0;
    y_arr[0] = 0.0;
    z_arr[0] = 0.0;
    r_arr[0] = 1.0;
    x_arr[1] = 2.5;
    y_arr[1] = 0.0;
    z_arr[1] = 0.0;
    r_arr[1] = 1.2;
    x_arr[2] = 5.0;
    y_arr[2] = 1.0;
    z_arr[2] = 0.0;
    r_arr[2] = 0.8;
    x_arr[3] = 1.0;
    y_arr[3] = 3.0;
    z_arr[3] = 0.5;
    r_arr[3] = 1.5;
    x_arr[4] = 3.5;
    y_arr[4] = 2.5;
    z_arr[4] = 1.0;
    r_arr[4] = 1.0;
    x_arr[5] = 0.5;
    y_arr[5] = 1.5;
    z_arr[5] = 3.0;
    r_arr[5] = 1.1;
    x_arr[6] = 4.0;
    y_arr[6] = 0.5;
    z_arr[6] = 2.5;
    r_arr[6] = 0.9;
    x_arr[7] = 2.0;
    y_arr[7] = 4.0;
    z_arr[7] = 2.0;
    r_arr[7] = 1.3;
    x_arr[8] = 6.0;
    y_arr[8] = 3.0;
    z_arr[8] = 1.5;
    r_arr[8] = 1.0;
    x_arr[9] = 7.0;
    y_arr[9] = 1.0;
    z_arr[9] = 3.0;
    r_arr[9] = 1.2;

    const input = AtomInput{
        .x = x_arr,
        .y = y_arr,
        .z = z_arr,
        .r = r_arr,
        .allocator = allocator,
    };

    // Calculate with f64
    const config64 = LeeRichardsConfig{
        .n_slices = 30,
        .probe_radius = 1.4,
    };
    var result64 = try calculateSasa(allocator, input, config64);
    defer result64.deinit();

    // Calculate with f32
    const config32 = LeeRichardsConfigGen(f32){
        .n_slices = 30,
        .probe_radius = 1.4,
    };
    var result32 = try calculateSasaf32(allocator, input, config32);
    defer result32.deinit();

    // f32 and f64 should produce similar results (within 0.5% tolerance)
    try std.testing.expectApproxEqRel(@as(f32, @floatCast(result64.total_area)), result32.total_area, 0.005);

    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqRel(
            @as(f32, @floatCast(result64.atom_areas[i])),
            result32.atom_areas[i],
            0.005,
        );
    }
}

test "calculateSasaParallelf32 - same as sequential f32" {
    const allocator = std.testing.allocator;

    const n_atoms = 50;
    const x_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(x_arr);
    const y_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(y_arr);
    const z_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(z_arr);
    const r_arr = try allocator.alloc(f64, n_atoms);
    defer allocator.free(r_arr);

    // Create a grid of atoms
    for (0..n_atoms) |i| {
        const fi: f64 = @floatFromInt(i);
        x_arr[i] = @mod(fi, 10.0) * 3.0;
        y_arr[i] = @mod(@floor(fi / 10.0), 10.0) * 3.0;
        z_arr[i] = @floor(fi / 100.0) * 3.0;
        r_arr[i] = 1.5;
    }

    const input = AtomInput{
        .x = x_arr,
        .y = y_arr,
        .z = z_arr,
        .r = r_arr,
        .allocator = allocator,
    };

    const config = LeeRichardsConfigGen(f32){
        .n_slices = 20,
        .probe_radius = 1.4,
    };

    // Calculate sequential
    var sequential_result = try calculateSasaf32(allocator, input, config);
    defer sequential_result.deinit();

    // Calculate parallel with 2 threads
    var parallel_result = try calculateSasaParallelf32(allocator, input, config, 2);
    defer parallel_result.deinit();

    // Verify total area matches
    try std.testing.expectApproxEqAbs(sequential_result.total_area, parallel_result.total_area, 1e-4);

    // Verify each atom area matches
    for (0..n_atoms) |i| {
        try std.testing.expectApproxEqAbs(
            sequential_result.atom_areas[i],
            parallel_result.atom_areas[i],
            1e-4,
        );
    }
}

// =============================================================================
// Error paths: thread spawn failure and allocation failure
// =============================================================================

/// Number of atoms in `fillErrorPathGrid`: enough for calculateSasaParallel to
/// split the work into several chunks (the minimum chunk size is 64).
const error_path_n_atoms = 400;

/// Fill a grid of overlapping atoms.
fn fillErrorPathGrid(x: []f64, y: []f64, z: []f64, r: []f64) void {
    for (x, y, z, r, 0..) |*xi, *yi, *zi, *ri, i| {
        xi.* = @as(f64, @floatFromInt(i % 8)) * 3.0;
        yi.* = @as(f64, @floatFromInt((i / 8) % 8)) * 3.0;
        zi.* = @as(f64, @floatFromInt(i / 64)) * 3.0;
        ri.* = 1.2 + @as(f64, @floatFromInt(i % 5)) * 0.1;
    }
}

test "chunkArcBuffer gives every chunk of the pool its own slice" {
    const allocator = std.testing.allocator;
    const buffer_size = 6;

    const cases = [_]struct { n_atoms: usize, chunk_size: usize }{
        .{ .n_atoms = 1, .chunk_size = 64 },
        .{ .n_atoms = 64, .chunk_size = 64 },
        .{ .n_atoms = 65, .chunk_size = 64 },
        .{ .n_atoms = 400, .chunk_size = 64 },
        .{ .n_atoms = 1000, .chunk_size = 100 },
    };

    for (cases) |case| {
        const arc_buffers = try allocArcBuffers(Arc, allocator, case.n_atoms, case.chunk_size, buffer_size);
        defer allocator.free(arc_buffers);

        // Walk the chunks the way ThreadPool.workerLoop hands them out.
        var n_chunks: usize = 0;
        var chunk_start: usize = 0;
        while (chunk_start < case.n_atoms) : (chunk_start += case.chunk_size) {
            const arc_buffer = chunkArcBuffer(Arc, arc_buffers, chunk_start, case.chunk_size, buffer_size);
            try std.testing.expectEqual(@as(usize, buffer_size), arc_buffer.len);
            try std.testing.expectEqual(arc_buffers.ptr + n_chunks * buffer_size, arc_buffer.ptr);
            n_chunks += 1;
        }
        try std.testing.expectEqual(n_chunks * buffer_size, arc_buffers.len);
    }
}

test "calculateSasaParallel - spawn failure after K workers matches serial" {
    const allocator = std.testing.allocator;

    var x: [error_path_n_atoms]f64 = undefined;
    var y: [error_path_n_atoms]f64 = undefined;
    var z: [error_path_n_atoms]f64 = undefined;
    var r: [error_path_n_atoms]f64 = undefined;
    fillErrorPathGrid(&x, &y, &z, &r);
    const input = AtomInput{ .x = &x, .y = &y, .z = &z, .r = &r, .allocator = allocator };

    const config = LeeRichardsConfig{ .n_slices = 20, .probe_radius = 1.4 };

    var serial = try calculateSasa(allocator, input, config);
    defer serial.deinit();

    // K == n_threads is the run without a failure.
    const n_threads = 4;
    for (0..n_threads + 1) |k| {
        thread_pool.testing.spawns_until_failure = k;
        defer thread_pool.testing.spawns_until_failure = null;

        var parallel = try calculateSasaParallel(allocator, input, config, n_threads);
        defer parallel.deinit();

        // No worker may outlive the call: its buffers are already freed.
        try std.testing.expectEqual(@as(usize, 0), thread_pool.testing.live_workers.load(.monotonic));
        try std.testing.expectEqualSlices(f64, serial.atom_areas, parallel.atom_areas);
        try std.testing.expectApproxEqRel(serial.total_area, parallel.total_area, 1e-12);
    }
}

test "calculateSasaParallelf32 - spawn failure after K workers matches serial" {
    const allocator = std.testing.allocator;

    var x: [error_path_n_atoms]f64 = undefined;
    var y: [error_path_n_atoms]f64 = undefined;
    var z: [error_path_n_atoms]f64 = undefined;
    var r: [error_path_n_atoms]f64 = undefined;
    fillErrorPathGrid(&x, &y, &z, &r);
    const input = AtomInput{ .x = &x, .y = &y, .z = &z, .r = &r, .allocator = allocator };

    const config = LeeRichardsConfigGen(f32){ .n_slices = 20, .probe_radius = 1.4 };

    var serial = try calculateSasaf32(allocator, input, config);
    defer serial.deinit();

    const n_threads = 4;
    for (0..n_threads + 1) |k| {
        thread_pool.testing.spawns_until_failure = k;
        defer thread_pool.testing.spawns_until_failure = null;

        var parallel = try calculateSasaParallelf32(allocator, input, config, n_threads);
        defer parallel.deinit();

        try std.testing.expectEqual(@as(usize, 0), thread_pool.testing.live_workers.load(.monotonic));
        try std.testing.expectEqualSlices(f32, serial.atom_areas, parallel.atom_areas);
        try std.testing.expectApproxEqRel(serial.total_area, parallel.total_area, 1e-5);
    }
}

/// Test bodies for `std.testing.checkAllAllocationFailures`. `n_threads` selects
/// the entry point: null for calculateSasa, a count for calculateSasaParallel.
/// A run that succeeds must return the reference areas, so an allocation
/// failure can neither leak nor be swallowed into a partly computed result.
const AllocationFailure = struct {
    fn runF64(allocator: Allocator, input: AtomInput, n_threads: ?usize, expected: []const f64) !void {
        const config = LeeRichardsConfig{ .n_slices = 20, .probe_radius = 1.4 };
        var result = if (n_threads) |n|
            try calculateSasaParallel(allocator, input, config, n)
        else
            try calculateSasa(allocator, input, config);
        defer result.deinit();
        try std.testing.expectEqualSlices(f64, expected, result.atom_areas);
    }

    fn runF32(allocator: Allocator, input: AtomInput, n_threads: ?usize, expected: []const f32) !void {
        const config = LeeRichardsConfigGen(f32){ .n_slices = 20, .probe_radius = 1.4 };
        var result = if (n_threads) |n|
            try calculateSasaParallelf32(allocator, input, config, n)
        else
            try calculateSasaf32(allocator, input, config);
        defer result.deinit();
        try std.testing.expectEqualSlices(f32, expected, result.atom_areas);
    }
};

test "calculateSasa and calculateSasaParallel - every allocation failure is reported without a leak" {
    const allocator = std.testing.allocator;

    var x: [error_path_n_atoms]f64 = undefined;
    var y: [error_path_n_atoms]f64 = undefined;
    var z: [error_path_n_atoms]f64 = undefined;
    var r: [error_path_n_atoms]f64 = undefined;
    fillErrorPathGrid(&x, &y, &z, &r);
    const input = AtomInput{ .x = &x, .y = &y, .z = &z, .r = &r, .allocator = allocator };

    var expected = try calculateSasa(allocator, input, .{ .n_slices = 20, .probe_radius = 1.4 });
    defer expected.deinit();
    const areas: []const f64 = expected.atom_areas;

    try std.testing.checkAllAllocationFailures(allocator, AllocationFailure.runF64, .{ input, null, areas });
    // One thread takes the direct path of parallelFor, four go through the pool.
    try std.testing.checkAllAllocationFailures(allocator, AllocationFailure.runF64, .{ input, 1, areas });
    try std.testing.checkAllAllocationFailures(allocator, AllocationFailure.runF64, .{ input, 4, areas });
}

test "calculateSasaf32 and calculateSasaParallelf32 - every allocation failure is reported without a leak" {
    const allocator = std.testing.allocator;

    var x: [error_path_n_atoms]f64 = undefined;
    var y: [error_path_n_atoms]f64 = undefined;
    var z: [error_path_n_atoms]f64 = undefined;
    var r: [error_path_n_atoms]f64 = undefined;
    fillErrorPathGrid(&x, &y, &z, &r);
    const input = AtomInput{ .x = &x, .y = &y, .z = &z, .r = &r, .allocator = allocator };

    var expected = try calculateSasaf32(allocator, input, .{ .n_slices = 20, .probe_radius = 1.4 });
    defer expected.deinit();
    const areas: []const f32 = expected.atom_areas;

    try std.testing.checkAllAllocationFailures(allocator, AllocationFailure.runF32, .{ input, null, areas });
    try std.testing.checkAllAllocationFailures(allocator, AllocationFailure.runF32, .{ input, 1, areas });
    try std.testing.checkAllAllocationFailures(allocator, AllocationFailure.runF32, .{ input, 4, areas });
}

// Far-apart atoms (issue #428): the neighbor grid used to be sized by the bounding box.

test "calculateSasa - far-apart atoms are isolated spheres" {
    const allocator = std.testing.allocator;

    // A second atom far from one at the origin. The dense grid of the first case had more
    // cells than a usize can count; the second needed hundreds of MB.
    const cases = [_]struct { position: [3]f64, radius: f64 }{
        .{ .position = .{ 34359738352.0, 34359738352.0, 0.0 }, .radius = 2.6 },
        .{ .position = .{ 2000.0, 2000.0, 2000.0 }, .radius = 1.7 },
    };
    const config = LeeRichardsConfig{ .n_slices = 20, .probe_radius = 1.4 };
    const config_f32 = LeeRichardsConfigGen(f32){ .n_slices = 20, .probe_radius = 1.4 };

    for (cases) |case| {
        var x = [_]f64{ 0.0, case.position[0] };
        var y = [_]f64{ 0.0, case.position[1] };
        var z = [_]f64{ 0.0, case.position[2] };
        var r = [_]f64{ case.radius, case.radius };
        const input = AtomInput{
            .x = &x,
            .y = &y,
            .z = &z,
            .r = &r,
            .allocator = allocator,
        };
        const expected = 4.0 * std.math.pi * (case.radius + 1.4) * (case.radius + 1.4);
        const expected_f32: f32 = @floatCast(expected);

        var sequential = try calculateSasa(allocator, input, config);
        defer sequential.deinit();
        var parallel = try calculateSasaParallel(allocator, input, config, 4);
        defer parallel.deinit();
        var sequential_f32 = try calculateSasaf32(allocator, input, config_f32);
        defer sequential_f32.deinit();
        var parallel_f32 = try calculateSasaParallelf32(allocator, input, config_f32, 4);
        defer parallel_f32.deinit();

        for (0..2) |i| {
            try std.testing.expectApproxEqRel(expected, sequential.atom_areas[i], 1e-12);
            try std.testing.expectApproxEqRel(expected, parallel.atom_areas[i], 1e-12);
            try std.testing.expectApproxEqRel(expected_f32, sequential_f32.atom_areas[i], 1e-6);
            try std.testing.expectApproxEqRel(expected_f32, parallel_f32.atom_areas[i], 1e-6);
        }
    }
}

test "calculateSasa - extreme finite coordinates" {
    const allocator = std.testing.allocator;

    const config = LeeRichardsConfig{ .n_slices = 20, .probe_radius = 1.4 };
    const config_f32 = LeeRichardsConfigGen(f32){ .n_slices = 20, .probe_radius = 1.4 };
    const expected = 4.0 * std.math.pi * 3.1 * 3.1;

    var x = [_]f64{ 0.0, 1.0e30 };
    var y = [_]f64{ 0.0, -1.0e30 };
    var z = [_]f64{ 0.0, 1.0e30 };
    var r = [_]f64{ 1.7, 1.7 };
    const input = AtomInput{
        .x = &x,
        .y = &y,
        .z = &z,
        .r = &r,
        .allocator = allocator,
    };

    // 1e30 is representable in both precisions: two isolated spheres
    {
        var result = try calculateSasa(allocator, input, config);
        defer result.deinit();
        var result_f32 = try calculateSasaf32(allocator, input, config_f32);
        defer result_f32.deinit();
        for (0..2) |i| {
            try std.testing.expectApproxEqRel(expected, result.atom_areas[i], 1e-12);
            try std.testing.expectApproxEqRel(@as(f32, @floatCast(expected)), result_f32.atom_areas[i], 1e-6);
        }
    }

    // 1e300 is representable in f64 only: f32 reports an error instead of building a grid
    // from infinite coordinates
    x[1] = 1.0e300;
    y[1] = -1.0e300;
    z[1] = 1.0e300;
    {
        var result = try calculateSasa(allocator, input, config);
        defer result.deinit();
        for (0..2) |i| {
            try std.testing.expectApproxEqRel(expected, result.atom_areas[i], 1e-12);
        }
        try std.testing.expectError(error.CoordinateRangeTooLarge, calculateSasaf32(allocator, input, config_f32));
        try std.testing.expectError(error.CoordinateRangeTooLarge, calculateSasaParallelf32(allocator, input, config_f32, 4));
    }
}

test "calculateSasa - a stray distant atom does not change the other areas" {
    const allocator = std.testing.allocator;

    // A 3 x 3 x 3 block of overlapping atoms, then one atom at a sentinel coordinate
    const n_compact = 27;
    var x: [n_compact + 1]f64 = undefined;
    var y: [n_compact + 1]f64 = undefined;
    var z: [n_compact + 1]f64 = undefined;
    var r: [n_compact + 1]f64 = undefined;
    for (0..n_compact) |i| {
        x[i] = @as(f64, @floatFromInt(i % 3)) * 2.9;
        y[i] = @as(f64, @floatFromInt((i / 3) % 3)) * 3.1;
        z[i] = @as(f64, @floatFromInt(i / 9)) * 2.7;
        r[i] = 1.2 + @as(f64, @floatFromInt(i % 5)) * 0.15;
    }
    x[n_compact] = 9999.999;
    y[n_compact] = 9999.999;
    z[n_compact] = 9999.999;
    r[n_compact] = 1.7;

    const compact_input = AtomInput{
        .x = x[0..n_compact],
        .y = y[0..n_compact],
        .z = z[0..n_compact],
        .r = r[0..n_compact],
        .allocator = allocator,
    };
    const stray_input = AtomInput{
        .x = &x,
        .y = &y,
        .z = &z,
        .r = &r,
        .allocator = allocator,
    };
    const config = LeeRichardsConfig{ .n_slices = 20, .probe_radius = 1.4 };
    const config_f32 = LeeRichardsConfigGen(f32){ .n_slices = 20, .probe_radius = 1.4 };
    const isolated = 4.0 * std.math.pi * 3.1 * 3.1;

    var compact = try calculateSasa(allocator, compact_input, config);
    defer compact.deinit();
    var stray = try calculateSasa(allocator, stray_input, config);
    defer stray.deinit();
    var compact_f32 = try calculateSasaf32(allocator, compact_input, config_f32);
    defer compact_f32.deinit();
    var stray_f32 = try calculateSasaf32(allocator, stray_input, config_f32);
    defer stray_f32.deinit();

    // The atom in the middle of the block is buried, so the block is a real test
    try std.testing.expect(compact.atom_areas[13] < compact.atom_areas[0]);

    // The stray atom extends the neighbor grid without changing its cells, so every other
    // atom keeps the same neighbor list and gets exactly the same area.
    try std.testing.expectEqualSlices(f64, compact.atom_areas, stray.atom_areas[0..n_compact]);
    try std.testing.expectEqualSlices(f32, compact_f32.atom_areas, stray_f32.atom_areas[0..n_compact]);
    try std.testing.expectApproxEqRel(isolated, stray.atom_areas[n_compact], 1e-12);
    try std.testing.expectApproxEqRel(@as(f32, @floatCast(isolated)), stray_f32.atom_areas[n_compact], 1e-6);
}

// =============================================================================
// Trig mode: exact arc angles against an independent reference
// =============================================================================

/// Fixtures and an independent Lee-Richards reference for the trig mode tests.
const trig_testing = struct {
    /// Linear congruential generator. The fixtures use it instead of
    /// `std.Random` so that their coordinates, and with them the pinned value
    /// of the fast mode, never change with the standard library.
    const Lcg = struct {
        state: u64,

        /// Uniform in [0, 1).
        fn next(self: *Lcg) f64 {
            self.state = self.state *% 6364136223846793005 +% 1442695040888963407;
            return @as(f64, @floatFromInt(self.state >> 11)) / 9007199254740992.0;
        }
    };

    const Fixture = struct {
        allocator: Allocator,
        x: []f64,
        y: []f64,
        z: []f64,
        r: []f64,

        fn init(allocator: Allocator, n_atoms: usize) !Fixture {
            const x = try allocator.alloc(f64, n_atoms);
            errdefer allocator.free(x);
            const y = try allocator.alloc(f64, n_atoms);
            errdefer allocator.free(y);
            const z = try allocator.alloc(f64, n_atoms);
            errdefer allocator.free(z);
            const r = try allocator.alloc(f64, n_atoms);
            return .{ .allocator = allocator, .x = x, .y = y, .z = z, .r = r };
        }

        fn deinit(self: *Fixture) void {
            self.allocator.free(self.x);
            self.allocator.free(self.y);
            self.allocator.free(self.z);
            self.allocator.free(self.r);
        }

        fn input(self: Fixture) AtomInput {
            return .{ .x = self.x, .y = self.y, .z = self.z, .r = self.r, .allocator = self.allocator };
        }

        /// The same atoms in the order `order[0], order[1], ...`.
        fn permuted(self: Fixture, order: []const usize) !Fixture {
            const result = try Fixture.init(self.allocator, order.len);
            for (order, 0..) |from, to| {
                result.x[to] = self.x[from];
                result.y[to] = self.y[from];
                result.z[to] = self.z[from];
                result.r[to] = self.r[from];
            }
            return result;
        }
    };

    /// `side`^3 atoms on a cubic grid with the given spacing, each coordinate
    /// moved by up to `jitter`, radii between 1.2 and 2.0 A.
    fn jitteredGrid(allocator: Allocator, side: usize, spacing: f64, jitter: f64, seed: u64) !Fixture {
        const fixture = try Fixture.init(allocator, side * side * side);
        var lcg = Lcg{ .state = seed };
        for (0..fixture.x.len) |i| {
            const ix: f64 = @floatFromInt(i % side);
            const iy: f64 = @floatFromInt((i / side) % side);
            const iz: f64 = @floatFromInt(i / (side * side));
            fixture.x[i] = spacing * ix + jitter * (2.0 * lcg.next() - 1.0);
            fixture.y[i] = spacing * iy + jitter * (2.0 * lcg.next() - 1.0);
            fixture.z[i] = spacing * iz + jitter * (2.0 * lcg.next() - 1.0);
            fixture.r[i] = 1.2 + 0.8 * lcg.next();
        }
        return fixture;
    }

    /// 125 atoms at the density of a protein interior (one atom per 12 A^3):
    /// about 60 neighbors per atom, so most neighbors go through the batches.
    fn proteinLikeCluster(allocator: Allocator) !Fixture {
        return jitteredGrid(allocator, 5, 2.3, 0.5, 430);
    }

    /// 216 strongly overlapping atoms (one atom per 2.2 A^3), where every atom
    /// is a neighbor of almost every other one.
    fn denseBlob(allocator: Allocator) !Fixture {
        return jitteredGrid(allocator, 6, 1.3, 0.25, 431);
    }

    /// A fixed pseudo-random order of 0..n (Fisher-Yates).
    fn shuffledOrder(allocator: Allocator, n: usize, seed: u64) ![]usize {
        const order = try allocator.alloc(usize, n);
        for (order, 0..) |*o, i| o.* = i;
        var lcg = Lcg{ .state = seed };
        var i = n;
        while (i > 1) : (i -= 1) {
            const j: usize = @intFromFloat(lcg.next() * @as(f64, @floatFromInt(i)));
            std.mem.swap(usize, &order[i - 1], &order[j]);
        }
        return order;
    }

    /// Length of the union of `intervals`, each `{ start, end }` with
    /// `start <= end`. Sorts `intervals`.
    fn unionLength(intervals: [][2]f64) f64 {
        std.mem.sort([2]f64, intervals, {}, struct {
            fn lessThan(_: void, a: [2]f64, b: [2]f64) bool {
                return a[0] < b[0];
            }
        }.lessThan);
        var covered: f64 = 0.0;
        var reach: f64 = 0.0;
        for (intervals) |interval| {
            const start = @max(interval[0], reach);
            if (interval[1] > start) {
                covered += interval[1] - start;
                reach = interval[1];
            }
        }
        return covered;
    }

    /// Independent Lee-Richards reference: the slicing of `atomArea`, but
    /// `std.math.acos` and `std.math.atan2` for every neighbor, every other
    /// atom tested as a neighbor (no neighbor list, no batches, no early
    /// distance cut-off) and its own interval union. Returns the per-atom
    /// areas, which the caller frees.
    fn reference(allocator: Allocator, input: AtomInput, n_slices: u32, probe_radius: f64) ![]f64 {
        const n_atoms = input.atomCount();
        const areas = try allocator.alloc(f64, n_atoms);
        errdefer allocator.free(areas);
        // An arc that crosses the angle 0 is split in two.
        const intervals = try allocator.alloc([2]f64, 2 * n_atoms);
        defer allocator.free(intervals);

        for (areas, 0..) |*area, i| {
            const ri = input.r[i] + probe_radius;
            const delta = 2.0 * ri / @as(f64, @floatFromInt(n_slices));
            var exposed_angle: f64 = 0.0; // summed over the slices

            for (0..n_slices) |k| {
                const slice_z = input.z[i] - ri + delta * (@as(f64, @floatFromInt(k)) + 0.5);
                const hi = slice_z - input.z[i];
                const ci2 = ri * ri - hi * hi; // squared radius of circle i
                if (ci2 <= 0) continue;
                const ci = @sqrt(ci2);

                var n_intervals: usize = 0;
                var buried = false;
                for (0..n_atoms) |j| {
                    if (j == i) continue;
                    const rj = input.r[j] + probe_radius;
                    const hj = slice_z - input.z[j];
                    const cj2 = rj * rj - hj * hj; // squared radius of circle j
                    if (cj2 <= 0) continue;
                    const cj = @sqrt(cj2);

                    const dx = input.x[j] - input.x[i];
                    const dy = input.y[j] - input.y[i];
                    const d = @sqrt(dx * dx + dy * dy);
                    if (d >= ci + cj) continue; // apart
                    if (d + ci <= cj) { // circle i inside circle j
                        buried = true;
                        break;
                    }
                    if (d + cj <= ci) continue; // circle j inside circle i

                    const cos_half = (ci2 + d * d - cj2) / (2.0 * ci * d);
                    const half = std.math.acos(std.math.clamp(cos_half, -1.0, 1.0));
                    // The covered arc is centered on the direction of j.
                    const start = @mod(std.math.atan2(dy, dx) - half, TWOPI);
                    const end = start + 2.0 * half;
                    if (end > TWOPI) {
                        intervals[n_intervals] = .{ 0.0, end - TWOPI };
                        intervals[n_intervals + 1] = .{ start, TWOPI };
                        n_intervals += 2;
                    } else {
                        intervals[n_intervals] = .{ start, end };
                        n_intervals += 1;
                    }
                }
                if (!buried) exposed_angle += TWOPI - unionLength(intervals[0..n_intervals]);
            }
            area.* = ri * delta * exposed_angle;
        }
        return areas;
    }

    /// Largest absolute per-atom difference between `expected` and `actual`.
    fn maxAbsDiff(comptime T: type, expected: []const f64, actual: []const T) !f64 {
        try std.testing.expectEqual(expected.len, actual.len);
        var max_diff: f64 = 0.0;
        for (expected, actual) |e, a| {
            max_diff = @max(max_diff, @abs(e - @as(f64, a)));
        }
        return max_diff;
    }

    fn sum(comptime T: type, areas: []const T) f64 {
        var total: f64 = 0.0;
        for (areas) |area| total += area;
        return total;
    }

    /// Per-atom tolerance in A^2 for f64 results against the reference: both
    /// use exact trigonometry, so only rounding is left (observed: 2e-14).
    const f64_tolerance = 1e-10;
    /// Per-atom tolerance in A^2 for f32 results against the f64 reference
    /// (observed: 2e-5).
    const f32_tolerance = 5e-4;

    /// Run every Lee-Richards entry point in exact mode on `fixture` and
    /// compare the per-atom areas with the reference.
    fn expectExactMatchesReference(fixture: Fixture, n_slices: u32) !void {
        const allocator = fixture.allocator;
        const input = fixture.input();
        const probe_radius = 1.4;

        const expected = try reference(allocator, input, n_slices, probe_radius);
        defer allocator.free(expected);

        const config = LeeRichardsConfig{ .n_slices = n_slices, .probe_radius = probe_radius, .trig = .exact };
        const config_f32 = LeeRichardsConfigGen(f32){ .n_slices = n_slices, .probe_radius = probe_radius, .trig = .exact };
        const config_gen = LeeRichardsConfigGen(f64){ .n_slices = n_slices, .probe_radius = probe_radius, .trig = .exact };

        {
            var result = try calculateSasa(allocator, input, config);
            defer result.deinit();
            try std.testing.expect(try maxAbsDiff(f64, expected, result.atom_areas) < f64_tolerance);
        }
        {
            var result = try calculateSasaParallel(allocator, input, config, 4);
            defer result.deinit();
            try std.testing.expect(try maxAbsDiff(f64, expected, result.atom_areas) < f64_tolerance);
        }
        {
            var result = try LeeRichardsGen(f64).calculateSasa(allocator, input, config_gen);
            defer result.deinit();
            try std.testing.expect(try maxAbsDiff(f64, expected, result.atom_areas) < f64_tolerance);
        }
        {
            var result = try LeeRichardsGen(f64).calculateSasaParallel(allocator, input, config_gen, 4);
            defer result.deinit();
            try std.testing.expect(try maxAbsDiff(f64, expected, result.atom_areas) < f64_tolerance);
        }
        {
            var result = try calculateSasaf32(allocator, input, config_f32);
            defer result.deinit();
            try std.testing.expect(try maxAbsDiff(f32, expected, result.atom_areas) < f32_tolerance);
        }
        {
            var result = try calculateSasaParallelf32(allocator, input, config_f32, 4);
            defer result.deinit();
            try std.testing.expect(try maxAbsDiff(f32, expected, result.atom_areas) < f32_tolerance);
        }
    }
};

test "TrigMode.fromString accepts exact and fast only" {
    try std.testing.expectEqual(TrigMode.exact, TrigMode.fromString("exact").?);
    try std.testing.expectEqual(TrigMode.fast, TrigMode.fromString("fast").?);
    try std.testing.expect(TrigMode.fromString("") == null);
    try std.testing.expect(TrigMode.fromString("Exact") == null);
    try std.testing.expect(TrigMode.fromString("approximate") == null);
}

test "exact trigonometry is the default of every Lee-Richards configuration" {
    try std.testing.expectEqual(TrigMode.exact, (LeeRichardsConfig{}).trig);
    try std.testing.expectEqual(TrigMode.exact, (LeeRichardsConfigGen(f32){}).trig);
    try std.testing.expectEqual(TrigMode.exact, (LeeRichardsConfigGen(f64){}).trig);

    const allocator = std.testing.allocator;
    var fixture = try trig_testing.proteinLikeCluster(allocator);
    defer fixture.deinit();
    const input = fixture.input();

    var by_default = try calculateSasa(allocator, input, .{});
    defer by_default.deinit();
    var exact = try calculateSasa(allocator, input, .{ .trig = .exact });
    defer exact.deinit();
    try std.testing.expectEqualSlices(f64, exact.atom_areas, by_default.atom_areas);
}

test "exact mode matches the independent reference on a protein-like cluster" {
    var fixture = try trig_testing.proteinLikeCluster(std.testing.allocator);
    defer fixture.deinit();
    try trig_testing.expectExactMatchesReference(fixture, 20);
    // An odd slice count puts a slice through every atom center.
    try trig_testing.expectExactMatchesReference(fixture, 7);
}

test "exact mode matches the independent reference on a dense blob" {
    var fixture = try trig_testing.denseBlob(std.testing.allocator);
    defer fixture.deinit();
    try trig_testing.expectExactMatchesReference(fixture, 20);
}

test "exact mode does not depend on the order of the atoms, fast mode does" {
    const allocator = std.testing.allocator;
    var fixture = try trig_testing.proteinLikeCluster(allocator);
    defer fixture.deinit();
    const n_atoms = fixture.x.len;

    const order = try trig_testing.shuffledOrder(allocator, n_atoms, 432);
    defer allocator.free(order);
    var shuffled = try fixture.permuted(order);
    defer shuffled.deinit();

    // Largest difference between the area of an atom in the original order
    // and the area of the same atom after shuffling.
    const Shuffle = struct {
        fn maxDiff(comptime T: type, original: []const T, reordered: []const T, new_order: []const usize) f64 {
            var max_diff: f64 = 0.0;
            for (new_order, 0..) |from, to| {
                max_diff = @max(max_diff, @abs(@as(f64, original[from]) - @as(f64, reordered[to])));
            }
            return max_diff;
        }
    };

    inline for (.{ TrigMode.exact, TrigMode.fast }) |trig| {
        var original = try calculateSasa(allocator, fixture.input(), .{ .trig = trig });
        defer original.deinit();
        var reordered = try calculateSasa(allocator, shuffled.input(), .{ .trig = trig });
        defer reordered.deinit();
        const diff = Shuffle.maxDiff(f64, original.atom_areas, reordered.atom_areas, order);

        var original_f32 = try calculateSasaf32(allocator, fixture.input(), .{ .trig = trig });
        defer original_f32.deinit();
        var reordered_f32 = try calculateSasaf32(allocator, shuffled.input(), .{ .trig = trig });
        defer reordered_f32.deinit();
        const diff_f32 = Shuffle.maxDiff(f32, original_f32.atom_areas, reordered_f32.atom_areas, order);

        switch (trig) {
            // Only rounding is left: a neighbor's slice radius is squared from
            // its square root in the batches and computed directly otherwise.
            .exact => {
                try std.testing.expect(diff < trig_testing.f64_tolerance);
                try std.testing.expect(diff_f32 < trig_testing.f32_tolerance);
            },
            // The order decides which neighbors get the approximation
            // (observed: 2e-2 A^2).
            .fast => {
                try std.testing.expect(diff > 1e-3);
                try std.testing.expect(diff_f32 > 1e-3);
            },
        }
    }
}

test "fast mode is selectable, keeps the values of zsasa 0.9.1 and is biased on a dense system" {
    const allocator = std.testing.allocator;
    var fixture = try trig_testing.denseBlob(allocator);
    defer fixture.deinit();
    const input = fixture.input();

    const expected = try trig_testing.reference(allocator, input, 20, 1.4);
    defer allocator.free(expected);
    const reference_total = trig_testing.sum(f64, expected);

    var exact = try calculateSasa(allocator, input, .{ .trig = .exact });
    defer exact.deinit();
    var fast = try calculateSasa(allocator, input, .{ .trig = .fast });
    defer fast.deinit();

    try std.testing.expectApproxEqRel(reference_total, trig_testing.sum(f64, exact.atom_areas), 1e-12);

    // Total of `zsasa calc --algorithm=lr --threads=1` from zsasa 0.9.1 (commit
    // 83a4a4c, before the trig mode existed) on this fixture written to JSON.
    const fast_total_0_9_1 = 776.668240633079;
    const fast_total = trig_testing.sum(f64, fast.atom_areas);
    try std.testing.expectApproxEqRel(fast_total_0_9_1, fast_total, 1e-13);

    // The approximation overestimates the exposed area (+0.12% here, 0.12 A^2
    // on the worst atom).
    try std.testing.expect(fast_total > reference_total * 1.001);
    try std.testing.expect(try trig_testing.maxAbsDiff(f64, expected, fast.atom_areas) > 0.1);
    try std.testing.expect(try trig_testing.maxAbsDiff(f64, expected, exact.atom_areas) < trig_testing.f64_tolerance);

    // Every entry point honors the mode and agrees with the others.
    {
        var parallel = try calculateSasaParallel(allocator, input, .{ .trig = .fast }, 4);
        defer parallel.deinit();
        try std.testing.expectEqualSlices(f64, fast.atom_areas, parallel.atom_areas);
    }
    {
        var generic = try LeeRichardsGen(f64).calculateSasa(allocator, input, .{ .trig = .fast });
        defer generic.deinit();
        try std.testing.expectEqualSlices(f64, fast.atom_areas, generic.atom_areas);
    }
    {
        var generic = try LeeRichardsGen(f64).calculateSasaParallel(allocator, input, .{ .trig = .fast }, 4);
        defer generic.deinit();
        try std.testing.expectEqualSlices(f64, fast.atom_areas, generic.atom_areas);
    }
    {
        var fast_f32 = try calculateSasaf32(allocator, input, .{ .trig = .fast });
        defer fast_f32.deinit();
        var parallel_f32 = try calculateSasaParallelf32(allocator, input, .{ .trig = .fast }, 4);
        defer parallel_f32.deinit();
        try std.testing.expectEqualSlices(f32, fast_f32.atom_areas, parallel_f32.atom_areas);
        // f32 fast follows f64 fast, not the exact reference.
        try std.testing.expect(try trig_testing.maxAbsDiff(f32, fast.atom_areas, fast_f32.atom_areas) < trig_testing.f32_tolerance);
        try std.testing.expect(trig_testing.sum(f32, fast_f32.atom_areas) > reference_total * 1.001);
    }
}

// =============================================================================
// Tangent circles
// =============================================================================

/// The tangent pair of issue #430: with the probe, a sphere of radius 2 at the
/// origin touches a sphere of radius 3 at (x_large, 0, 0) from inside when
/// |x_large| = 1. The small atom is buried and the large one fully exposed.
/// An odd slice count puts a slice through both centers, where the two slice
/// circles are tangent; in every other slice the small circle lies strictly
/// inside the large one.
const tangent_testing = struct {
    const small_radius = 0.6;
    const large_radius = 1.6;
    const probe_radius = 1.4;
    /// 4 pi (1.6 + 1.4)^2
    const large_area = 4.0 * std.math.pi * 9.0;
    const max_padding = 12;

    /// The pair plus `n_padding` atoms of radius 0.1 on the z axis, at
    /// |z| >= 1.6. With the probe they are neighbors of both atoms of the pair
    /// but do not reach the plane z = 0, so with one slice they change no area.
    /// They only move the tangent neighbor between the 8-wide batch, the
    /// 4-wide batch and the scalar remainder.
    const System = struct {
        x: [2 + max_padding]f64 = @splat(0.0),
        y: [2 + max_padding]f64 = @splat(0.0),
        z: [2 + max_padding]f64 = @splat(0.0),
        r: [2 + max_padding]f64 = @splat(0.1),
        n_atoms: usize,

        fn init(x_large: f64, n_padding: usize) System {
            std.debug.assert(n_padding <= max_padding);
            var system = System{ .n_atoms = 2 + n_padding };
            system.r[0] = small_radius;
            system.x[1] = x_large;
            system.r[1] = large_radius;
            for (0..n_padding) |k| {
                const height = 1.6 + 0.15 * @as(f64, @floatFromInt(k / 2));
                system.z[2 + k] = if (k % 2 == 0) height else -height;
            }
            return system;
        }

        fn input(self: *System, allocator: Allocator) AtomInput {
            const n = self.n_atoms;
            return .{ .x = self.x[0..n], .y = self.y[0..n], .z = self.z[0..n], .r = self.r[0..n], .allocator = allocator };
        }
    };

    /// Areas of the small and the large atom, as f64.
    fn pairAreas(
        comptime T: type,
        allocator: Allocator,
        input: AtomInput,
        n_slices: u32,
        trig: TrigMode,
        n_threads: ?usize,
    ) ![2]f64 {
        const config = LeeRichardsConfigGen(T){ .n_slices = n_slices, .probe_radius = probe_radius, .trig = trig };
        var result = if (n_threads) |n|
            try LeeRichardsGen(T).calculateSasaParallel(allocator, input, config, n)
        else
            try LeeRichardsGen(T).calculateSasa(allocator, input, config);
        defer result.deinit();
        return .{ result.atom_areas[0], result.atom_areas[1] };
    }

    /// The same through the non-generic f64 entry points.
    fn pairAreasNonGeneric(allocator: Allocator, input: AtomInput, n_slices: u32, trig: TrigMode, n_threads: ?usize) ![2]f64 {
        const config = LeeRichardsConfig{ .n_slices = n_slices, .probe_radius = probe_radius, .trig = trig };
        var result = if (n_threads) |n|
            try calculateSasaParallel(allocator, input, config, n)
        else
            try calculateSasa(allocator, input, config);
        defer result.deinit();
        return .{ result.atom_areas[0], result.atom_areas[1] };
    }

    fn expectBuriedAndExposed(comptime T: type, areas: [2]f64) !void {
        try std.testing.expectEqual(@as(f64, 0.0), areas[0]);
        try std.testing.expectApproxEqRel(@as(f64, large_area), areas[1], if (T == f32) 1e-6 else 1e-13);
    }
};

test "tangent circles: an atom touching a larger one from inside is buried, the larger one stays exposed" {
    const allocator = std.testing.allocator;

    inline for (.{ f64, f32 }) |T| {
        for ([_]f64{ 1.0, -1.0 }) |x_large| {
            var system = tangent_testing.System.init(x_large, 0);
            const input = system.input(allocator);
            // Odd counts have a slice through the centers, where the circles are tangent.
            for ([_]u32{ 1, 2, 20, 21 }) |n_slices| {
                for ([_]TrigMode{ .exact, .fast }) |trig| {
                    for ([_]?usize{ null, 2 }) |n_threads| {
                        const areas = try tangent_testing.pairAreas(T, allocator, input, n_slices, trig, n_threads);
                        try tangent_testing.expectBuriedAndExposed(T, areas);
                    }
                }
            }
        }
    }

    // The non-generic f64 copy
    for ([_]f64{ 1.0, -1.0 }) |x_large| {
        var system = tangent_testing.System.init(x_large, 0);
        const input = system.input(allocator);
        for ([_]u32{ 1, 2, 20, 21 }) |n_slices| {
            for ([_]TrigMode{ .exact, .fast }) |trig| {
                for ([_]?usize{ null, 2 }) |n_threads| {
                    const areas = try tangent_testing.pairAreasNonGeneric(allocator, input, n_slices, trig, n_threads);
                    try tangent_testing.expectBuriedAndExposed(f64, areas);
                }
            }
        }
    }
}

test "tangent circles: the tangent neighbor is handled the same in the batches and in the scalar remainder" {
    const allocator = std.testing.allocator;

    // 1 to 13 neighbors per atom: the tangent neighbor falls in the scalar
    // remainder, in the 4-wide batch and in the 8-wide batch.
    for (0..tangent_testing.max_padding + 1) |n_padding| {
        for ([_]f64{ 1.0, -1.0 }) |x_large| {
            var system = tangent_testing.System.init(x_large, n_padding);
            const input = system.input(allocator);
            for ([_]TrigMode{ .exact, .fast }) |trig| {
                inline for (.{ f64, f32 }) |T| {
                    const areas = try tangent_testing.pairAreas(T, allocator, input, 1, trig, null);
                    try tangent_testing.expectBuriedAndExposed(T, areas);
                }
                const areas = try tangent_testing.pairAreasNonGeneric(allocator, input, 1, trig, null);
                try tangent_testing.expectBuriedAndExposed(f64, areas);
            }
        }
    }
}

test "tangent circles: areas are continuous across tangency" {
    const allocator = std.testing.allocator;

    inline for (.{ f64, f32 }) |T| {
        // Center distances next to 1 in the working precision, and a little further away.
        const one: T = 1.0;
        const distances = [_]f64{
            std.math.nextAfter(T, one, 0.0),
            std.math.nextAfter(T, one, 2.0),
            @as(T, 1.0 - 1e-6),
            @as(T, 1.0 + 1e-6),
        };
        for (distances) |distance| {
            // With one slice (thickness 4, radius 2), a small circle that
            // sticks out by eps exposes an arc of 2 sqrt(3 eps), and covers an
            // arc of 2 sqrt(4 eps / 3) of the large circle (thickness 6, radius 3).
            const eps = @abs(distance - 1.0);
            const small_bound = 2.0 * 16.0 * @sqrt(3.0 * eps) + 1e-5;
            const large_bound = 2.0 * 36.0 * @sqrt(4.0 * eps / 3.0) + 1e-4;

            for ([_]f64{ 1.0, -1.0 }) |sign| {
                var system = tangent_testing.System.init(sign * distance, 0);
                const input = system.input(allocator);
                for ([_]u32{ 1, 21 }) |n_slices| {
                    for ([_]TrigMode{ .exact, .fast }) |trig| {
                        const areas = try tangent_testing.pairAreas(T, allocator, input, n_slices, trig, null);
                        try std.testing.expect(areas[0] >= 0.0);
                        try std.testing.expect(areas[0] < small_bound);
                        try std.testing.expectApproxEqAbs(@as(f64, tangent_testing.large_area), areas[1], large_bound);
                    }
                }
            }
        }
    }
}

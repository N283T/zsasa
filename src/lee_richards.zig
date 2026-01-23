const std = @import("std");
const types = @import("types.zig");
const neighbor_list_mod = @import("neighbor_list.zig");
const thread_pool = @import("thread_pool.zig");

const Allocator = std.mem.Allocator;
const AtomInput = types.AtomInput;
const SasaResult = types.SasaResult;
const Vec3 = types.Vec3;
const NeighborList = neighbor_list_mod.NeighborList;

const TWOPI: f64 = 2.0 * std.math.pi;

/// Thread-safe allocator for parallel workers
const thread_safe_allocator = std.heap.page_allocator;

/// Configuration for Lee-Richards algorithm
pub const LeeRichardsConfig = struct {
    /// Number of slices per atom diameter
    n_slices: u32 = 20,
    /// Water probe radius in Angstroms
    probe_radius: f64 = 1.4,
};

/// Arc interval representing a buried portion of a circle
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

/// Calculate SASA for a single atom using slice-based method
fn atomArea(
    atom_idx: usize,
    x: []const f64,
    y: []const f64,
    z: []const f64,
    radii: []const f64,
    neighbor_list: *const NeighborList,
    n_slices: u32,
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

        for (neighbors) |j| {
            const zj = z[j];
            const dj = @abs(zj - slice_z);
            const Rj = radii[j];

            // Check if neighbor j intersects this slice
            if (dj >= Rj) continue;

            const Rj_prime2 = Rj * Rj - dj * dj;
            const Rj_prime = @sqrt(Rj_prime2);

            // Distance between atom centers in xy-plane
            const dx = x[j] - xi;
            const dy = y[j] - yi;
            const dij = @sqrt(dx * dx + dy * dy);

            // Check circle interactions
            if (dij >= Ri_prime + Rj_prime) {
                // Circles don't touch
                continue;
            }

            // Handle near-zero distance (atoms co-located in xy-plane)
            if (dij < 1e-10) {
                if (Rj_prime > Ri_prime) {
                    // Circle i is inside circle j
                    is_buried = true;
                    break;
                }
                // Circle j is inside circle i, skip
                continue;
            }

            if (dij + Ri_prime < Rj_prime) {
                // Circle i completely inside circle j
                is_buried = true;
                break;
            }
            if (dij + Rj_prime < Ri_prime) {
                // Circle j completely inside circle i
                continue;
            }

            // Calculate arc blocked by neighbor j
            // alpha = half-angle of blocked arc
            const cos_alpha = (Ri_prime2 + dij * dij - Rj_prime2) / (2.0 * Ri_prime * dij);
            const alpha = std.math.acos(std.math.clamp(cos_alpha, -1.0, 1.0));

            // beta = angle to center of blocked arc
            const beta = std.math.atan2(dy, dx) + std.math.pi;

            var inf = beta - alpha;
            var sup = beta + alpha;

            // Normalize angles to [0, 2π)
            // Use while loops for robustness against numerical edge cases
            while (inf < 0) inf += TWOPI;
            while (inf >= TWOPI) inf -= TWOPI;
            while (sup <= 0) sup += TWOPI;
            while (sup > TWOPI) sup -= TWOPI;

            // Store arc(s)
            if (sup < inf) {
                // Arc crosses 0, split into two arcs
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
/// Thread safety: All fields are read-only except `atom_areas` which has
/// disjoint write access (each thread writes to different indices).
const ParallelContext = struct {
    x: []const f64,
    y: []const f64,
    z: []const f64,
    radii: []const f64,
    neighbor_list: *const NeighborList,
    n_slices: u32,
    max_arc_buffer_size: usize,
    atom_areas: []f64,
};

/// Worker function for parallel Lee-Richards calculation.
/// Processes atoms from chunk_start to chunk_end.
fn parallelLeeRichardsWorker(ctx: ParallelContext, chunk_start: usize, chunk_end: usize) f64 {
    // Allocate arc buffer for this chunk (thread-safe allocator)
    const arc_buffer = thread_safe_allocator.alloc(Arc, ctx.max_arc_buffer_size) catch {
        // Log error - this is very unlikely with page_allocator but shouldn't fail silently
        std.log.err("Lee-Richards worker: allocation failed for chunk {d}-{d}", .{ chunk_start, chunk_end });
        return 0.0;
    };
    defer thread_safe_allocator.free(arc_buffer);

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
/// - `config`: Configuration parameters (n_slices, probe_radius)
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
        .max_arc_buffer_size = max_arc_buffer_size,
        .atom_areas = atom_areas,
    };

    // Chunk size heuristic:
    // - Minimum 64 atoms per chunk to amortize thread overhead
    // - Target 4 chunks per thread for load balancing
    const chunk_size = @max(64, n_atoms / (actual_threads * 4));

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

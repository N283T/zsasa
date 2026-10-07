const std = @import("std");
const types = @import("types.zig");

const Vec3 = types.Vec3;
const Vec3Gen = types.Vec3Gen;
const Allocator = std.mem.Allocator;

/// Smallest cell budget of a grid, whatever the atom count. A small or sparse selection (a few
/// ions in a simulation box, a ligand, two domains far apart) legitimately has far more cells
/// than atoms, so the budget must not shrink with the atom count. 2^21 cells is 16 MB of
/// transient memory and covers a box of about 845 Å per side at a typical cell size of 6.6 Å.
const min_grid_cells: usize = 1 << 21;

/// Additional cell budget per atom for large systems. A compact structure needs about one cell
/// per atom or fewer, so this only limits systems whose bounding box is mostly empty. At 8 bytes
/// per cell the grid then stays smaller than the neighbor list built from it.
const grid_cells_per_atom: usize = 16;

/// Maximum number of grid cells for `n_atoms` atoms.
fn maxGridCells(n_atoms: usize) usize {
    return @max(min_grid_cells, n_atoms *| grid_cells_per_atom);
}

/// Cell size and per-axis cell counts of a grid.
fn GridShape(comptime T: type) type {
    return struct {
        cell_size: T,
        nx: usize,
        ny: usize,
        nz: usize,
    };
}

/// Choose the cell size and cell counts for a bounding box with the given extents.
///
/// When `ceil(extent / min_cell_size)` cells per axis fit in `max_cells`, that grid is returned
/// unchanged. Otherwise the cell size is increased until the grid fits. A cell only has to be
/// at least as large as the interaction cutoff, so a larger cell finds exactly the same
/// neighbor pairs; it only puts more candidates in each cell.
///
/// The returned counts satisfy `nx * ny * nz <= max_cells`, so the product cannot overflow.
/// Returns `error.CoordinateRangeTooLarge` when an extent is infinite or NaN, which is how
/// non-finite coordinates and coordinate ranges wider than `T` can represent show up here.
fn fitGrid(
    comptime T: type,
    extent: [3]T,
    min_cell_size: T,
    max_cells: usize,
) error{CoordinateRangeTooLarge}!GridShape(T) {
    var max_extent: T = 0.0;
    for (extent) |e| {
        if (!std.math.isFinite(e)) return error.CoordinateRangeTooLarge;
        max_extent = @max(max_extent, e);
    }

    // Integers up to 2^52 are exact in f64, so the comparisons and casts below are exact.
    const budget: f64 = @min(@as(f64, @floatFromInt(@max(1, max_cells))), 0x1p52);

    // No single axis can have more cells than the whole budget.
    const axis_floor = max_extent / @as(T, @floatCast(budget));

    var cell_size = min_cell_size;
    while (true) {
        var counts: [3]f64 = undefined;
        var total: f64 = 1.0;
        for (extent, &counts) |e, *count| {
            count.* = @max(1.0, @ceil(e / cell_size));
            total *= count.*;
        }
        if (total <= budget) {
            return .{
                .cell_size = cell_size,
                .nx = @intFromFloat(counts[0]),
                .ny = @intFromFloat(counts[1]),
                .nz = @intFromFloat(counts[2]),
            };
        }

        // Jump to the per-axis bound first. From there on every count is finite, even when
        // the requested cell size is tiny compared with the extent.
        if (cell_size < axis_floor) {
            cell_size = axis_floor;
            continue;
        }

        // Grow by the cube root of the excess. A flat or linear box needs several passes,
        // because an axis that is already one cell cannot shrink; the 1% floor guarantees
        // progress, and at `max_extent` the grid is a single cell.
        const growth: T = @floatCast(@max(std.math.cbrt(total / budget), 1.01));
        cell_size = @min(cell_size * growth, max_extent);
    }
}

/// Compute cell index from position coordinates.
///
/// The position must lie inside the bounding box the grid was built from, so that each
/// quotient is at most the cell count of its axis and the float-to-integer casts are in range.
fn computeCellIndex(
    comptime T: type,
    x: T,
    y: T,
    z: T,
    x_min: T,
    y_min: T,
    z_min: T,
    cell_size: T,
    nx: usize,
    ny: usize,
    nz: usize,
) usize {
    const ix = @as(usize, @intFromFloat(@max(@as(T, 0.0), (x - x_min) / cell_size)));
    const iy = @as(usize, @intFromFloat(@max(@as(T, 0.0), (y - y_min) / cell_size)));
    const iz = @as(usize, @intFromFloat(@max(@as(T, 0.0), (z - z_min) / cell_size)));
    return @min(iz, nz - 1) * nx * ny + @min(iy, ny - 1) * nx + @min(ix, nx - 1);
}

/// Generic spatial hash grid with flat storage (counting-sort)
pub fn CellListGen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        const Self = @This();

        atom_indices: []u32,
        cell_offsets: []u32, // length = n_cells + 1
        nx: usize,
        ny: usize,
        nz: usize,
        cell_size: T,
        x_min: T,
        y_min: T,
        z_min: T,
        allocator: Allocator,

        /// Build spatial grid from atom positions
        /// cell_size should be >= 2 * (max_radius + probe_radius) for correctness
        ///
        /// The grid uses `cell_size` unless the bounding box of the atoms would then need more
        /// cells than `maxGridCells` allows; in that case the cells are made larger (see
        /// `fitGrid`), and the `cell_size` field holds the size actually used.
        pub fn init(
            allocator: Allocator,
            positions: []const Vec,
            cell_size: T,
        ) !Self {
            return initWithMaxCells(allocator, positions, cell_size, maxGridCells(positions.len));
        }

        /// Same as `init` with an explicit cell budget.
        fn initWithMaxCells(
            allocator: Allocator,
            positions: []const Vec,
            min_cell_size: T,
            max_cells: usize,
        ) !Self {
            if (positions.len == 0) return error.NoAtoms;
            if (min_cell_size <= 0.0) return error.InvalidCellSize;

            // Compute bounding box
            var x_min = positions[0].x;
            var x_max = positions[0].x;
            var y_min = positions[0].y;
            var y_max = positions[0].y;
            var z_min = positions[0].z;
            var z_max = positions[0].z;

            for (positions) |pos| {
                x_min = @min(x_min, pos.x);
                x_max = @max(x_max, pos.x);
                y_min = @min(y_min, pos.y);
                y_max = @max(y_max, pos.y);
                z_min = @min(z_min, pos.z);
                z_max = @max(z_max, pos.z);
            }

            // Add padding to avoid edge cases
            x_min -= min_cell_size;
            y_min -= min_cell_size;
            z_min -= min_cell_size;
            x_max += min_cell_size;
            y_max += min_cell_size;
            z_max += min_cell_size;

            // Calculate grid dimensions (minimum 1 cell), bounded by the cell budget
            const shape = try fitGrid(
                T,
                .{ x_max - x_min, y_max - y_min, z_max - z_min },
                min_cell_size,
                max_cells,
            );
            const cell_size = shape.cell_size;
            const nx = shape.nx;
            const ny = shape.ny;
            const nz = shape.nz;
            const n_cells = nx * ny * nz;

            // Pass 1: count atoms per cell
            const counts = try allocator.alloc(u32, n_cells);
            defer allocator.free(counts);
            @memset(counts, 0);

            for (positions) |pos| {
                const idx = computeCellIndex(T, pos.x, pos.y, pos.z, x_min, y_min, z_min, cell_size, nx, ny, nz);
                counts[idx] += 1;
            }

            // Build prefix sum into cell_offsets
            const cell_offsets = try allocator.alloc(u32, n_cells + 1);
            errdefer allocator.free(cell_offsets);
            cell_offsets[0] = 0;
            for (0..n_cells) |i| {
                cell_offsets[i + 1] = cell_offsets[i] + counts[i];
            }

            // Pass 2: place atoms (reuse counts as write cursors)
            const atom_indices = try allocator.alloc(u32, positions.len);
            errdefer allocator.free(atom_indices);
            @memset(counts, 0);

            for (positions, 0..) |pos, i| {
                const idx = computeCellIndex(T, pos.x, pos.y, pos.z, x_min, y_min, z_min, cell_size, nx, ny, nz);
                atom_indices[cell_offsets[idx] + counts[idx]] = @intCast(i);
                counts[idx] += 1;
            }

            return Self{
                .atom_indices = atom_indices,
                .cell_offsets = cell_offsets,
                .nx = nx,
                .ny = ny,
                .nz = nz,
                .cell_size = cell_size,
                .x_min = x_min,
                .y_min = y_min,
                .z_min = z_min,
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.atom_indices);
            self.allocator.free(self.cell_offsets);
        }

        /// Get atom indices in a cell
        pub fn getCellAtoms(self: Self, cell_idx: usize) []const u32 {
            return self.atom_indices[self.cell_offsets[cell_idx]..self.cell_offsets[cell_idx + 1]];
        }

        fn getCellIndex(self: Self, pos: Vec) usize {
            return computeCellIndex(
                T,
                pos.x,
                pos.y,
                pos.z,
                self.x_min,
                self.y_min,
                self.z_min,
                self.cell_size,
                self.nx,
                self.ny,
                self.nz,
            );
        }

        /// Get cell coordinates from index
        pub fn getCellCoords(self: Self, idx: usize) struct { ix: usize, iy: usize, iz: usize } {
            const iz = idx / (self.nx * self.ny);
            const remainder = idx % (self.nx * self.ny);
            const iy = remainder / self.nx;
            const ix = remainder % self.nx;
            return .{ .ix = ix, .iy = iy, .iz = iz };
        }

        /// Get cell index from coordinates (returns null if out of bounds)
        pub fn getCellIndexFromCoords(self: Self, ix: i64, iy: i64, iz: i64) ?usize {
            if (ix < 0 or iy < 0 or iz < 0) return null;
            const uix = @as(usize, @intCast(ix));
            const uiy = @as(usize, @intCast(iy));
            const uiz = @as(usize, @intCast(iz));
            if (uix >= self.nx or uiy >= self.ny or uiz >= self.nz) return null;
            return uiz * self.nx * self.ny + uiy * self.nx + uix;
        }
    };
}

const IterMode = enum { count, fill };

/// Shared iteration over all neighbor pairs using cell list.
/// In count mode: increments counts[i] and counts[j] for each pair.
/// In fill mode: writes indices into neighbor_indices using offsets + counts as cursors.
///
/// SAFETY: The count and fill passes MUST iterate in identical order. The `neighbor_indices`
/// buffer must be sized exactly as computed by a prior count-mode pass (via prefix-sum offsets).
/// The `counts` array is reused as write cursors in fill mode and must be zeroed beforehand.
fn processNeighborPairs(
    comptime T: type,
    comptime mode: IterMode,
    cell_list: anytype,
    positions: []const Vec3Gen(T),
    radii: []const T,
    probe_radius: T,
    counts: []u32,
    neighbor_indices: []u32,
    offsets: []const u32,
) void {
    const n_cells = cell_list.nx * cell_list.ny * cell_list.nz;
    for (0..n_cells) |cell_idx| {
        const cell1_atoms = cell_list.getCellAtoms(cell_idx);
        if (cell1_atoms.len == 0) continue;

        const coords = cell_list.getCellCoords(cell_idx);
        const cix = @as(i64, @intCast(coords.ix));
        const ciy = @as(i64, @intCast(coords.iy));
        const ciz = @as(i64, @intCast(coords.iz));

        // Check all 27 neighboring cells (including self)
        var cdz: i64 = -1;
        while (cdz <= 1) : (cdz += 1) {
            var cdy: i64 = -1;
            while (cdy <= 1) : (cdy += 1) {
                var cdx: i64 = -1;
                while (cdx <= 1) : (cdx += 1) {
                    const ncell = cell_list.getCellIndexFromCoords(cix + cdx, ciy + cdy, ciz + cdz);
                    if (ncell) |nidx| {
                        // Only process cell pairs where cell_idx <= nidx to avoid duplicates
                        if (nidx < cell_idx) continue;
                        const same_cell = cell_idx == nidx;
                        const cell2_atoms = cell_list.getCellAtoms(nidx);

                        for (cell1_atoms) |ai| {
                            for (cell2_atoms) |aj| {
                                if (ai == aj) continue;
                                if (same_cell and aj <= ai) continue;

                                const pi = positions[ai];
                                const pj = positions[aj];
                                const dx = pi.x - pj.x;
                                const dy = pi.y - pj.y;
                                const dz = pi.z - pj.z;
                                const dist_sq = dx * dx + dy * dy + dz * dz;

                                const cutoff = radii[ai] + radii[aj] + 2.0 * probe_radius;

                                if (dist_sq < cutoff * cutoff) {
                                    if (mode == .count) {
                                        counts[ai] += 1;
                                        counts[aj] += 1;
                                    } else {
                                        std.debug.assert(offsets[ai] + counts[ai] < offsets[ai + 1]);
                                        neighbor_indices[offsets[ai] + counts[ai]] = aj;
                                        counts[ai] += 1;
                                        std.debug.assert(offsets[aj] + counts[aj] < offsets[aj + 1]);
                                        neighbor_indices[offsets[aj] + counts[aj]] = ai;
                                        counts[aj] += 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Generic pre-computed neighbor list with flat storage (two-pass)
pub fn NeighborListGen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    const CellListT = CellListGen(T);
    return struct {
        const Self = @This();

        neighbor_indices: []u32,
        offsets: []u32, // length = n_atoms + 1
        allocator: Allocator,

        /// Build neighbor list from atom positions and radii
        /// Two atoms i, j are neighbors if distance(i, j) < r[i] + r[j] + 2*probe_radius
        pub fn init(
            allocator: Allocator,
            positions: []const Vec,
            radii: []const T,
            probe_radius: T,
        ) !Self {
            return initWithMaxCells(allocator, positions, radii, probe_radius, maxGridCells(positions.len));
        }

        /// Same as `init` with an explicit cell budget for the underlying grid.
        fn initWithMaxCells(
            allocator: Allocator,
            positions: []const Vec,
            radii: []const T,
            probe_radius: T,
            max_cells: usize,
        ) !Self {
            const n_atoms = positions.len;
            if (n_atoms == 0) return error.NoAtoms;
            std.debug.assert(radii.len == n_atoms);

            // Find maximum radius and validate
            var max_radius: T = 0.0;
            for (radii) |r| {
                if (r < 0.0) return error.InvalidRadius;
                max_radius = @max(max_radius, r);
            }

            const cell_size = 2.0 * (max_radius + probe_radius);

            var cell_list = try CellListT.initWithMaxCells(allocator, positions, cell_size, max_cells);
            defer cell_list.deinit();

            // Pass 1: count neighbor pairs
            const counts = try allocator.alloc(u32, n_atoms);
            defer allocator.free(counts);
            @memset(counts, 0);

            processNeighborPairs(T, .count, cell_list, positions, radii, probe_radius, counts, counts[0..0], counts[0..0]);

            // Build prefix sum
            const offsets = try allocator.alloc(u32, n_atoms + 1);
            errdefer allocator.free(offsets);
            offsets[0] = 0;
            for (0..n_atoms) |i| {
                offsets[i + 1] = offsets[i] + counts[i];
            }

            // Allocate flat neighbor buffer
            const total = offsets[n_atoms];
            const neighbor_indices = try allocator.alloc(u32, total);
            errdefer allocator.free(neighbor_indices);

            // Pass 2: fill neighbors (reuse counts as write cursors)
            @memset(counts, 0);
            processNeighborPairs(T, .fill, cell_list, positions, radii, probe_radius, counts, neighbor_indices, offsets);

            return Self{
                .neighbor_indices = neighbor_indices,
                .offsets = offsets,
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.neighbor_indices);
            self.allocator.free(self.offsets);
        }

        /// Get neighbors for atom i
        pub fn getNeighbors(self: Self, atom_idx: usize) []const u32 {
            return self.neighbor_indices[self.offsets[atom_idx]..self.offsets[atom_idx + 1]];
        }
    };
}

/// Type aliases
pub const CellList = CellListGen(f64);
pub const NeighborList = NeighborListGen(f64);
pub const CellListf32 = CellListGen(f32);
pub const NeighborListf32 = NeighborListGen(f32);

// Tests

test "CellList - single atom" {
    const allocator = std.testing.allocator;

    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
    };

    var cell_list = try CellList.init(allocator, positions, 5.0);
    defer cell_list.deinit();

    // Should have at least 1 cell
    try std.testing.expect(cell_list.nx >= 1);
    try std.testing.expect(cell_list.ny >= 1);
    try std.testing.expect(cell_list.nz >= 1);
}

test "CellList - atoms in different cells" {
    const allocator = std.testing.allocator;

    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 20.0, .y = 0.0, .z = 0.0 }, // Far apart
    };

    var cell_list = try CellList.init(allocator, positions, 5.0);
    defer cell_list.deinit();

    // Should have multiple cells in x direction
    try std.testing.expect(cell_list.nx > 1);
}

test "NeighborList - two far atoms have no neighbors" {
    const allocator = std.testing.allocator;

    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 100.0, .y = 0.0, .z = 0.0 }, // Very far apart
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Neither should be neighbors
    try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(0).len);
    try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(1).len);
}

test "NeighborList - two touching atoms are neighbors" {
    const allocator = std.testing.allocator;

    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 2.0, .y = 0.0, .z = 0.0 }, // Close together
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Distance = 2, cutoff = 1 + 1 + 2*1.4 = 4.8 → neighbors
    try std.testing.expectEqual(@as(usize, 1), neighbor_list.getNeighbors(0).len);
    try std.testing.expectEqual(@as(usize, 1), neighbor_list.getNeighbors(1).len);
    try std.testing.expectEqual(@as(u32, 1), neighbor_list.getNeighbors(0)[0]);
    try std.testing.expectEqual(@as(u32, 0), neighbor_list.getNeighbors(1)[0]);
}

test "NeighborList - symmetry (j in neighbors[i] iff i in neighbors[j])" {
    const allocator = std.testing.allocator;

    // Three atoms in a row
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 3.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 6.0, .y = 0.0, .z = 0.0 },
    };
    const radii = &[_]f64{ 1.0, 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Check symmetry
    for (0..3) |i| {
        for (neighbor_list.getNeighbors(i)) |j| {
            // j should have i as neighbor
            var found = false;
            for (neighbor_list.getNeighbors(j)) |k| {
                if (k == @as(u32, @intCast(i))) {
                    found = true;
                    break;
                }
            }
            try std.testing.expect(found);
        }
    }
}

test "NeighborList - cluster of close atoms" {
    const allocator = std.testing.allocator;

    // 4 atoms at corners of a small tetrahedron
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 2.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 1.0, .y = 1.732, .z = 0.0 },
        Vec3{ .x = 1.0, .y = 0.577, .z = 1.633 },
    };
    const radii = &[_]f64{ 1.0, 1.0, 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // All atoms should be neighbors of each other (distances are ~2 Å)
    for (0..4) |i| {
        try std.testing.expectEqual(@as(usize, 3), neighbor_list.getNeighbors(i).len);
    }
}

test "NeighborList - boundary atoms exactly at cutoff" {
    const allocator = std.testing.allocator;

    // Two atoms exactly at cutoff distance
    // cutoff = r1 + r2 + 2*probe = 1 + 1 + 2*1.4 = 4.8
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 4.79, .y = 0.0, .z = 0.0 }, // Just inside cutoff
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Should be neighbors (just inside cutoff)
    try std.testing.expectEqual(@as(usize, 1), neighbor_list.getNeighbors(0).len);
}

test "NeighborList - boundary atoms outside cutoff" {
    const allocator = std.testing.allocator;

    // Two atoms just outside cutoff distance
    // cutoff = r1 + r2 + 2*probe = 1 + 1 + 2*1.4 = 4.8
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 4.81, .y = 0.0, .z = 0.0 }, // Just outside cutoff
    };
    const radii = &[_]f64{ 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Should NOT be neighbors (just outside cutoff)
    try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(0).len);
}

test "NeighborList - different radii" {
    const allocator = std.testing.allocator;

    // Two atoms with different radii
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 5.0, .y = 0.0, .z = 0.0 },
    };
    const radii = &[_]f64{ 2.0, 1.5 }; // cutoff = 2.0 + 1.5 + 2*1.4 = 6.3
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Distance 5.0 < cutoff 6.3 → neighbors
    try std.testing.expectEqual(@as(usize, 1), neighbor_list.getNeighbors(0).len);
    try std.testing.expectEqual(@as(usize, 1), neighbor_list.getNeighbors(1).len);
}

test "NeighborList - all atoms in same cell" {
    const allocator = std.testing.allocator;

    // 5 atoms very close together, all should fall in same cell
    // cell_size = 2 * (max_radius + probe) = 2 * (1.0 + 1.4) = 4.8
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 0.5, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 0.0, .y = 0.5, .z = 0.0 },
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.5 },
        Vec3{ .x = 0.5, .y = 0.5, .z = 0.5 },
    };
    const radii = &[_]f64{ 1.0, 1.0, 1.0, 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // All atoms should be neighbors of each other (4 neighbors each)
    for (0..5) |i| {
        try std.testing.expectEqual(@as(usize, 4), neighbor_list.getNeighbors(i).len);
    }
}

test "NeighborList - no duplicate entries" {
    const allocator = std.testing.allocator;

    // Create atoms that span multiple cells to test duplicate prevention
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 3.0, .y = 0.0, .z = 0.0 }, // Different cell but within cutoff
        Vec3{ .x = 6.0, .y = 0.0, .z = 0.0 }, // Different cell but within cutoff of atom 1
        Vec3{ .x = 0.0, .y = 3.0, .z = 0.0 }, // Different cell
    };
    const radii = &[_]f64{ 1.0, 1.0, 1.0, 1.0 };
    const probe_radius = 1.4;

    var neighbor_list = try NeighborList.init(allocator, positions, radii, probe_radius);
    defer neighbor_list.deinit();

    // Check for duplicates in each neighbor list
    for (0..4) |i| {
        const neighbors = neighbor_list.getNeighbors(i);
        // Check all pairs for duplicates
        for (0..neighbors.len) |j| {
            for (j + 1..neighbors.len) |k| {
                try std.testing.expect(neighbors[j] != neighbors[k]);
            }
        }
    }
}

test "NeighborList - invalid negative radius" {
    const allocator = std.testing.allocator;

    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
    };
    const radii = &[_]f64{-1.0}; // Invalid negative radius
    const probe_radius = 1.4;

    const result = NeighborList.init(allocator, positions, radii, probe_radius);
    try std.testing.expectError(error.InvalidRadius, result);
}

test "CellList - invalid cell_size" {
    const allocator = std.testing.allocator;

    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
    };

    // Zero cell_size
    const result1 = CellList.init(allocator, positions, 0.0);
    try std.testing.expectError(error.InvalidCellSize, result1);

    // Negative cell_size
    const result2 = CellList.init(allocator, positions, -1.0);
    try std.testing.expectError(error.InvalidCellSize, result2);
}

// Grid bounds (issue #428)

/// Check a neighbor list against an all-pairs search: same neighbor sets, no self entries
/// and no duplicates. The all-pairs relation is symmetric, so this also checks symmetry.
fn expectMatchesBruteForce(
    comptime T: type,
    neighbor_list: NeighborListGen(T),
    positions: []const Vec3Gen(T),
    radii: []const T,
    probe_radius: T,
) !void {
    const allocator = std.testing.allocator;
    const n_atoms = positions.len;
    const listed = try allocator.alloc(bool, n_atoms);
    defer allocator.free(listed);

    for (0..n_atoms) |i| {
        @memset(listed, false);
        for (neighbor_list.getNeighbors(i)) |j| {
            try std.testing.expect(j != i);
            try std.testing.expect(!listed[j]);
            listed[j] = true;
        }
        for (0..n_atoms) |j| {
            const dx = positions[i].x - positions[j].x;
            const dy = positions[i].y - positions[j].y;
            const dz = positions[i].z - positions[j].z;
            const cutoff = radii[i] + radii[j] + 2.0 * probe_radius;
            const expected = i != j and dx * dx + dy * dy + dz * dz < cutoff * cutoff;
            try std.testing.expectEqual(expected, listed[j]);
        }
    }
}

/// Fill `positions` with clusters of atoms around the given centers and `radii` with
/// protein-like radii. Deterministic for a given seed.
fn fillClusters(
    comptime T: type,
    seed: u64,
    centers: []const Vec3Gen(T),
    spread: T,
    positions: []Vec3Gen(T),
    radii: []T,
) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const random = prng.random();
    for (positions, radii, 0..) |*pos, *r, i| {
        const center = centers[i % centers.len];
        pos.* = .{
            .x = center.x + spread * random.float(T),
            .y = center.y + spread * random.float(T),
            .z = center.z + spread * random.float(T),
        };
        r.* = 1.2 + 0.8 * random.float(T);
    }
}

test "fitGrid - a grid within the budget is returned unchanged" {
    inline for (.{ f32, f64 }) |T| {
        const shape = try fitGrid(T, .{ 64.5, 33.0, 6.6 }, 6.6, maxGridCells(1));
        try std.testing.expectEqual(@as(T, 6.6), shape.cell_size);
        try std.testing.expectEqual(@as(usize, 10), shape.nx);
        try std.testing.expectEqual(@as(usize, 5), shape.ny);
        try std.testing.expectEqual(@as(usize, 1), shape.nz);

        // A degenerate box still gets one cell per axis
        const point = try fitGrid(T, .{ 0.0, 0.0, 0.0 }, 6.6, 1);
        try std.testing.expectEqual(@as(usize, 1), point.nx * point.ny * point.nz);
    }
}

test "fitGrid - an oversized grid gets larger cells that fit the budget" {
    const extents = [_][3]f64{
        .{ 2000.0, 2000.0, 2000.0 }, // cube
        .{ 3.4e10, 3.4e10, 10.4 }, // flat: issue #428 overflow input
        .{ 1.0e7, 12.4, 12.4 }, // linear
        .{ 1.0e30, 1.0e30, 1.0e30 },
        .{ 1.0e300, 5.0, 1.0e-300 },
        .{ 123.0, 4567.0, 89012.0 },
    };
    const budgets = [_]usize{ 1, 2, 7, 64, 1000, min_grid_cells, std.math.maxInt(usize) };

    inline for (.{ f32, f64 }) |T| {
        for (extents) |extent_f64| {
            if (T == f32 and extent_f64[0] > std.math.floatMax(f32)) continue;
            const extent = [3]T{ @floatCast(extent_f64[0]), @floatCast(extent_f64[1]), @floatCast(extent_f64[2]) };
            for (budgets) |budget| {
                for ([_]T{ 6.2, 1.0e-30 }) |min_cell_size| {
                    const shape = try fitGrid(T, extent, min_cell_size, budget);
                    try std.testing.expect(std.math.isFinite(shape.cell_size));
                    try std.testing.expect(shape.cell_size >= min_cell_size);
                    try std.testing.expect(shape.nx >= 1 and shape.ny >= 1 and shape.nz >= 1);

                    const n_cells = try std.math.mul(usize, try std.math.mul(usize, shape.nx, shape.ny), shape.nz);
                    try std.testing.expect(n_cells <= budget);

                    // The cells still cover the whole box
                    try std.testing.expect(@as(T, @floatFromInt(shape.nx)) * shape.cell_size >= extent[0] * 0.999);
                    try std.testing.expect(@as(T, @floatFromInt(shape.ny)) * shape.cell_size >= extent[1] * 0.999);
                    try std.testing.expect(@as(T, @floatFromInt(shape.nz)) * shape.cell_size >= extent[2] * 0.999);
                }
            }
        }
    }
}

test "fitGrid - cells are not enlarged much more than the budget requires" {
    // 100^3 cells requested, 1000 allowed: the ideal answer is 10 cells per axis.
    const shape = try fitGrid(f64, .{ 600.0, 600.0, 600.0 }, 6.0, 1000);
    const n_cells = shape.nx * shape.ny * shape.nz;
    try std.testing.expect(n_cells <= 1000);
    try std.testing.expect(n_cells >= 700);
}

test "fitGrid - non-finite extent is an error" {
    inline for (.{ f32, f64 }) |T| {
        const inf = std.math.inf(T);
        const nan = std.math.nan(T);
        try std.testing.expectError(error.CoordinateRangeTooLarge, fitGrid(T, .{ inf, 1.0, 1.0 }, 6.2, 1000));
        try std.testing.expectError(error.CoordinateRangeTooLarge, fitGrid(T, .{ 1.0, nan, 1.0 }, 6.2, 1000));
        try std.testing.expectError(error.CoordinateRangeTooLarge, fitGrid(T, .{ 1.0, 1.0, -inf }, 6.2, 1000));
    }
}

test "CellList - grid size is bounded for far-apart atoms" {
    const allocator = std.testing.allocator;

    // Issue #428: 2^32 x 2^32 x 3 cells of 8 Å; the product used to wrap around.
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 34359738352.0, .y = 34359738352.0, .z = 0.0 },
    };

    var cell_list = try CellList.init(allocator, positions, 8.0);
    defer cell_list.deinit();

    const n_cells = cell_list.nx * cell_list.ny * cell_list.nz;
    try std.testing.expect(n_cells <= maxGridCells(positions.len));
    try std.testing.expectEqual(n_cells + 1, cell_list.cell_offsets.len);
    try std.testing.expect(cell_list.cell_size > 8.0);
    try std.testing.expectEqual(@as(u32, 2), cell_list.cell_offsets[n_cells]);
}

test "CellList - compact structure keeps the requested cell size" {
    const allocator = std.testing.allocator;

    var positions: [200]Vec3 = undefined;
    var radii: [200]f64 = undefined;
    fillClusters(f64, 1, &.{.{ .x = 10.0, .y = -5.0, .z = 3.0 }}, 40.0, &positions, &radii);

    var cell_list = try CellList.init(allocator, &positions, 6.6);
    defer cell_list.deinit();

    try std.testing.expectEqual(@as(f64, 6.6), cell_list.cell_size);
}

test "NeighborList - far-apart atoms have no neighbors" {
    const allocator = std.testing.allocator;

    inline for (.{ f32, f64 }) |T| {
        const cases = [_][2]Vec3Gen(T){
            // Issue #428: cell count overflow
            .{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 34359738352.0, .y = 34359738352.0, .z = 0.0 } },
            // Issue #428: hundreds of MB for two atoms
            .{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 2000.0, .y = 2000.0, .z = 2000.0 } },
            .{ .{ .x = -1.0e30, .y = 0.0, .z = 0.0 }, .{ .x = 1.0e30, .y = 1.0e30, .z = 1.0e30 } },
        };
        for (cases) |positions| {
            var neighbor_list = try NeighborListGen(T).init(allocator, &positions, &.{ 2.6, 2.6 }, 1.4);
            defer neighbor_list.deinit();

            try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(0).len);
            try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(1).len);
        }
    }
}

test "NeighborList - extreme coordinates" {
    const allocator = std.testing.allocator;

    // f64 holds 1e300 and the distance to it: two isolated atoms.
    {
        const positions = &[_]Vec3{
            Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
            Vec3{ .x = 1.0e300, .y = -1.0e300, .z = 1.0e300 },
        };
        var neighbor_list = try NeighborList.init(allocator, positions, &.{ 1.7, 1.7 }, 1.4);
        defer neighbor_list.deinit();
        try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(0).len);
        try std.testing.expectEqual(@as(usize, 0), neighbor_list.getNeighbors(1).len);
    }

    // The width of this range overflows f64.
    {
        const max = std.math.floatMax(f64);
        const positions = &[_]Vec3{
            Vec3{ .x = -max, .y = 0.0, .z = 0.0 },
            Vec3{ .x = max, .y = 0.0, .z = 0.0 },
        };
        try std.testing.expectError(
            error.CoordinateRangeTooLarge,
            NeighborList.init(allocator, positions, &.{ 1.7, 1.7 }, 1.4),
        );
    }

    // 1e300 becomes infinite when the f32 paths cast it.
    {
        const Vec3f32 = Vec3Gen(f32);
        const positions = &[_]Vec3f32{
            Vec3f32{ .x = 0.0, .y = 0.0, .z = 0.0 },
            Vec3f32{ .x = @floatCast(@as(f64, 1.0e300)), .y = 0.0, .z = 0.0 },
        };
        try std.testing.expectError(
            error.CoordinateRangeTooLarge,
            NeighborListf32.init(allocator, positions, &.{ 1.7, 1.7 }, 1.4),
        );
    }

    // A radius too large for the cell size to be finite.
    {
        const positions = &[_]Vec3{
            Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
            Vec3{ .x = 5.0, .y = 0.0, .z = 0.0 },
        };
        try std.testing.expectError(
            error.CoordinateRangeTooLarge,
            NeighborList.init(allocator, positions, &.{ std.math.floatMax(f64), 1.7 }, 1.4),
        );
    }
}

test "NeighborList - sparse clusters match brute force" {
    const allocator = std.testing.allocator;

    inline for (.{ f32, f64 }) |T| {
        const Vec = Vec3Gen(T);
        // Clusters far apart along every axis, plus two that touch each other.
        const centers = [_]Vec{
            .{ .x = 0.0, .y = 0.0, .z = 0.0 },
            .{ .x = 5000.0, .y = 0.0, .z = 0.0 },
            .{ .x = -3000.0, .y = 4000.0, .z = 8000.0 },
            .{ .x = 5008.0, .y = 3.0, .z = -4.0 },
            .{ .x = 9999.999, .y = 9999.999, .z = 9999.999 },
        };
        var positions: [300]Vec = undefined;
        var radii: [300]T = undefined;
        fillClusters(T, 428, &centers, 14.0, &positions, &radii);

        var neighbor_list = try NeighborListGen(T).init(allocator, &positions, &radii, 1.4);
        defer neighbor_list.deinit();
        try expectMatchesBruteForce(T, neighbor_list, &positions, &radii, 1.4);

        // The clusters are dense enough for this to be a real test.
        var total: usize = 0;
        for (0..positions.len) |i| total += neighbor_list.getNeighbors(i).len;
        try std.testing.expect(total > positions.len * 4);
    }
}

test "NeighborList - any cell budget gives the brute-force neighbors" {
    const allocator = std.testing.allocator;

    inline for (.{ f32, f64 }) |T| {
        const Vec = Vec3Gen(T);
        const layouts = [_][]const Vec{
            // Compact
            &.{.{ .x = 0.0, .y = 0.0, .z = 0.0 }},
            // Linear, flat and diagonal arrangements of clusters
            &.{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 60.0, .y = 0.0, .z = 0.0 }, .{ .x = 130.0, .y = 0.0, .z = 0.0 } },
            &.{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 0.0, .y = 70.0, .z = 90.0 }, .{ .x = 0.0, .y = -80.0, .z = 20.0 } },
            &.{ .{ .x = -50.0, .y = -50.0, .z = -50.0 }, .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 75.0, .y = 75.0, .z = 75.0 } },
        };
        var positions: [160]Vec = undefined;
        var radii: [160]T = undefined;

        for (layouts, 0..) |centers, layout_idx| {
            fillClusters(T, 1000 + layout_idx, centers, 18.0, &positions, &radii);
            for ([_]usize{ 1, 2, 5, 27, 100, 999, 20_000, min_grid_cells }) |max_cells| {
                var neighbor_list = try NeighborListGen(T).initWithMaxCells(allocator, &positions, &radii, 1.4, max_cells);
                defer neighbor_list.deinit();
                try expectMatchesBruteForce(T, neighbor_list, &positions, &radii, 1.4);
            }
        }
    }
}

test "NeighborList - symmetry with a bounded grid" {
    const allocator = std.testing.allocator;

    // Two clusters far enough apart that the grid cells are enlarged.
    const centers = [_]Vec3{
        .{ .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .x = 4000.0, .y = -4000.0, .z = 4000.0 },
    };
    var positions: [120]Vec3 = undefined;
    var radii: [120]f64 = undefined;
    fillClusters(f64, 7, &centers, 12.0, &positions, &radii);

    var neighbor_list = try NeighborList.init(allocator, &positions, &radii, 1.4);
    defer neighbor_list.deinit();

    var n_pairs: usize = 0;
    for (0..positions.len) |i| {
        for (neighbor_list.getNeighbors(i)) |j| {
            try std.testing.expect(std.mem.indexOfScalar(u32, neighbor_list.getNeighbors(j), @intCast(i)) != null);
            n_pairs += 1;
        }
    }
    try std.testing.expect(n_pairs > 0);
}

test "NeighborList - a stray distant atom does not change the other neighbor sets" {
    const allocator = std.testing.allocator;

    inline for (.{ f32, f64 }) |T| {
        const Vec = Vec3Gen(T);
        const n_compact = 150;
        var positions: [n_compact + 1]Vec = undefined;
        var radii: [n_compact + 1]T = undefined;
        fillClusters(T, 99, &.{.{ .x = 12.0, .y = 34.0, .z = 56.0 }}, 22.0, positions[0..n_compact], radii[0..n_compact]);
        positions[n_compact] = .{ .x = 9999.999, .y = 9999.999, .z = 9999.999 };
        radii[n_compact] = 1.7;

        var compact = try NeighborListGen(T).init(allocator, positions[0..n_compact], radii[0..n_compact], 1.4);
        defer compact.deinit();
        var with_stray = try NeighborListGen(T).init(allocator, &positions, &radii, 1.4);
        defer with_stray.deinit();

        try std.testing.expectEqual(@as(usize, 0), with_stray.getNeighbors(n_compact).len);

        const sorted_a = try allocator.alloc(u32, n_compact);
        defer allocator.free(sorted_a);
        const sorted_b = try allocator.alloc(u32, n_compact);
        defer allocator.free(sorted_b);
        for (0..n_compact) |i| {
            const a = compact.getNeighbors(i);
            const b = with_stray.getNeighbors(i);
            try std.testing.expectEqual(a.len, b.len);
            @memcpy(sorted_a[0..a.len], a);
            @memcpy(sorted_b[0..b.len], b);
            std.mem.sort(u32, sorted_a[0..a.len], {}, std.sort.asc(u32));
            std.mem.sort(u32, sorted_b[0..b.len], {}, std.sort.asc(u32));
            try std.testing.expectEqualSlices(u32, sorted_a[0..a.len], sorted_b[0..b.len]);
        }
    }
}

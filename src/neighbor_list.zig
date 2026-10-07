const std = @import("std");
const types = @import("types.zig");

const Vec3 = types.Vec3;
const Vec3Gen = types.Vec3Gen;
const Allocator = std.mem.Allocator;

/// Smallest budget of a dense grid, whatever the atom count. A small selection (a ligand, a
/// few ions) has more cells than atoms even when it is compact; 2^16 cells take 0.5 MB.
const min_dense_cells: usize = 1 << 16;

/// Dense cells allowed per atom. A compact structure has about one cell per atom or fewer, so
/// only systems whose bounding box is mostly empty exceed this. Around 64 cells per atom,
/// scanning the empty cells of a dense grid takes as long as sorting the atoms for the sparse
/// layout; a dense grid of that size needs 512 bytes per atom (8 bytes per cell).
const dense_cells_per_atom: usize = 64;

/// Largest dense grid for `n_atoms` atoms. A grid with more cells is stored sparsely.
fn maxDenseCells(n_atoms: usize) usize {
    return @max(min_dense_cells, n_atoms *| dense_cells_per_atom);
}

/// Bits of a sparse cell key that hold the cell coordinate along one axis.
const sparse_axis_bits = 21;

/// Largest number of cells along one axis. 2^21 cells of 6.6 Å span 1.4e7 Å.
const max_axis_cells = 1 << sparse_axis_bits;

/// Cell size, per-axis cell counts and storage layout of a grid.
fn GridShape(comptime T: type) type {
    return struct {
        cell_size: T,
        nx: usize,
        ny: usize,
        nz: usize,
        /// Whether one slot per cell fits the dense budget
        dense: bool,
    };
}

/// Choose the cell size, cell counts and layout for a bounding box with the given extents.
///
/// The cell size is `min_cell_size` and each axis has `ceil(extent / min_cell_size)` cells,
/// as long as no axis needs more than `max_axis_cells` cells. Beyond that the cell size is
/// increased so that the widest axis has `max_axis_cells` cells. A cell only has to be at
/// least as large as the interaction cutoff, so larger cells find the same neighbor pairs.
///
/// Every count is at most `max_axis_cells` (2^21), so the casts are in range and the product
/// of the three fits 64 bits. The grid is dense when that product is at most
/// `max_dense_cells`.
///
/// Returns `error.CoordinateRangeTooLarge` when an extent is infinite or NaN, which is how
/// non-finite coordinates and coordinate ranges wider than `T` can represent show up here.
fn fitGrid(
    comptime T: type,
    extent: [3]T,
    min_cell_size: T,
    max_dense_cells: usize,
) error{CoordinateRangeTooLarge}!GridShape(T) {
    var max_extent: T = 0.0;
    for (extent) |e| {
        if (!std.math.isFinite(e)) return error.CoordinateRangeTooLarge;
        max_extent = @max(max_extent, e);
    }

    const cell_size = @max(min_cell_size, max_extent / @as(T, max_axis_cells));

    var counts: [3]usize = undefined;
    var n_cells: u64 = 1;
    for (extent, &counts) |e, *count| {
        // Clamped as a float, before the cast. Atoms past the last cell of an axis are put
        // into it (see `axisCell`), which only adds candidate pairs.
        count.* = @intFromFloat(@min(@max(1.0, @ceil(e / cell_size)), @as(T, max_axis_cells)));
        n_cells *= count.*;
    }

    return .{
        .cell_size = cell_size,
        .nx = counts[0],
        .ny = counts[1],
        .nz = counts[2],
        .dense = n_cells <= max_dense_cells,
    };
}

/// Cell coordinate of a position along one axis, clamped to the last cell.
///
/// The position must lie inside the bounding box the grid was built from: the quotient is
/// then at most about the cell count of the axis, and the float-to-integer cast is in range.
inline fn axisCell(comptime T: type, v: T, v_min: T, cell_size: T, n: usize) usize {
    return @min(@as(usize, @intFromFloat(@max(@as(T, 0.0), (v - v_min) / cell_size))), n - 1);
}

/// Compute cell index from position coordinates (dense layout)
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
    const ix = axisCell(T, x, x_min, cell_size, nx);
    const iy = axisCell(T, y, y_min, cell_size, ny);
    const iz = axisCell(T, z, z_min, cell_size, nz);
    return iz * nx * ny + iy * nx + ix;
}

/// Key of the cell at the given cell coordinates (sparse layout).
///
/// Keys order cells by z, then y, then x, which is the order of the dense cell index.
fn cellKey(ix: usize, iy: usize, iz: usize) u64 {
    return (@as(u64, iz) << (2 * sparse_axis_bits)) | (@as(u64, iy) << sparse_axis_bits) | @as(u64, ix);
}

/// Generic spatial hash grid with flat storage (counting-sort)
///
/// Two layouts share the storage. The dense layout has one slot per cell of the grid. The
/// sparse layout has one slot per occupied cell and is used when the bounding box of the
/// atoms has more cells than `maxDenseCells` allows (atoms far apart). Both keep the cells in
/// the same order and the atoms of a cell in ascending index order, so the neighbor lists
/// built from them are identical.
pub fn CellListGen(comptime T: type) type {
    const Vec = Vec3Gen(T);
    return struct {
        const Self = @This();

        atom_indices: []u32,
        cell_offsets: []u32, // length = number of slots + 1
        /// Sparse layout: key of each occupied cell, ascending. Empty in the dense layout.
        cell_keys: []u64,
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
        /// The `cell_size` field holds the size actually used, which is larger only when an
        /// axis would otherwise have more than `max_axis_cells` cells (see `fitGrid`).
        pub fn init(
            allocator: Allocator,
            positions: []const Vec,
            cell_size: T,
        ) !Self {
            return initWithDenseLimit(allocator, positions, cell_size, maxDenseCells(positions.len));
        }

        /// Same as `init` with an explicit limit on the size of a dense grid.
        fn initWithDenseLimit(
            allocator: Allocator,
            positions: []const Vec,
            min_cell_size: T,
            max_dense_cells: usize,
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

            // Calculate grid dimensions (minimum 1 cell) and choose the layout
            const shape = try fitGrid(
                T,
                .{ x_max - x_min, y_max - y_min, z_max - z_min },
                min_cell_size,
                max_dense_cells,
            );
            if (!shape.dense) return initSparse(allocator, positions, shape, x_min, y_min, z_min);

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
                .cell_keys = &.{},
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

        /// Build the sparse layout: memory and time depend on the number of atoms, not on
        /// the volume of the bounding box.
        fn initSparse(
            allocator: Allocator,
            positions: []const Vec,
            shape: GridShape(T),
            x_min: T,
            y_min: T,
            z_min: T,
        ) !Self {
            // One entry per atom: cell key in the high bits, atom index in the low 32 bits.
            // Sorting the entries orders the cells by key and the atoms of a cell by index,
            // which is the order the counting sort of the dense layout produces.
            const entries = try allocator.alloc(u128, positions.len);
            defer allocator.free(entries);

            for (positions, entries, 0..) |pos, *entry, i| {
                const key = cellKey(
                    axisCell(T, pos.x, x_min, shape.cell_size, shape.nx),
                    axisCell(T, pos.y, y_min, shape.cell_size, shape.ny),
                    axisCell(T, pos.z, z_min, shape.cell_size, shape.nz),
                );
                entry.* = (@as(u128, key) << 32) | @as(u32, @intCast(i));
            }
            std.mem.sortUnstable(u128, entries, {}, std.sort.asc(u128));

            var n_occupied: usize = 1;
            for (entries[1..], entries[0 .. entries.len - 1]) |entry, previous| {
                if (entry >> 32 != previous >> 32) n_occupied += 1;
            }

            const cell_keys = try allocator.alloc(u64, n_occupied);
            errdefer allocator.free(cell_keys);
            const cell_offsets = try allocator.alloc(u32, n_occupied + 1);
            errdefer allocator.free(cell_offsets);
            const atom_indices = try allocator.alloc(u32, positions.len);
            errdefer allocator.free(atom_indices);

            var slot: usize = 0;
            for (entries, atom_indices, 0..) |entry, *atom_index, i| {
                const key: u64 = @intCast(entry >> 32);
                if (i == 0 or key != cell_keys[slot - 1]) {
                    cell_keys[slot] = key;
                    cell_offsets[slot] = @intCast(i);
                    slot += 1;
                }
                atom_index.* = @truncate(entry);
            }
            cell_offsets[n_occupied] = @intCast(positions.len);

            return Self{
                .atom_indices = atom_indices,
                .cell_offsets = cell_offsets,
                .cell_keys = cell_keys,
                .nx = shape.nx,
                .ny = shape.ny,
                .nz = shape.nz,
                .cell_size = shape.cell_size,
                .x_min = x_min,
                .y_min = y_min,
                .z_min = z_min,
                .allocator = allocator,
            };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.atom_indices);
            self.allocator.free(self.cell_offsets);
            self.allocator.free(self.cell_keys);
        }

        /// Whether only the occupied cells are stored
        pub fn isSparse(self: Self) bool {
            return self.cell_keys.len != 0;
        }

        /// Get atom indices in a cell (dense layout) or in an occupied cell (sparse layout)
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

        /// Get cell coordinates of an occupied cell (sparse layout)
        pub fn getSparseCellCoords(self: Self, slot: usize) struct { ix: usize, iy: usize, iz: usize } {
            const key = self.cell_keys[slot];
            const axis_mask = max_axis_cells - 1;
            return .{
                .ix = @intCast(key & axis_mask),
                .iy = @intCast((key >> sparse_axis_bits) & axis_mask),
                .iz = @intCast(key >> (2 * sparse_axis_bits)),
            };
        }

        /// Get the slot of the occupied cell at the given coordinates, searching the slots
        /// from `first_slot` on (sparse layout). Returns null if the cell is out of bounds,
        /// holds no atoms, or comes before `first_slot`.
        pub fn findSparseCell(self: Self, ix: i64, iy: i64, iz: i64, first_slot: usize) ?usize {
            if (ix < 0 or iy < 0 or iz < 0) return null;
            const uix = @as(usize, @intCast(ix));
            const uiy = @as(usize, @intCast(iy));
            const uiz = @as(usize, @intCast(iz));
            if (uix >= self.nx or uiy >= self.ny or uiz >= self.nz) return null;

            const keys = self.cell_keys;
            const key = cellKey(uix, uiy, uiz);
            if (key < keys[first_slot]) return null;

            // The cells next to a cell are usually a few slots ahead of it, so bracket the
            // key with doubling steps before bisecting. `keys[low] <= key` holds throughout.
            var low = first_slot;
            var step: usize = 1;
            while (low + step < keys.len and keys[low + step] <= key) {
                low += step;
                step *= 2;
            }
            var high = @min(low + step, keys.len);
            while (high - low > 1) {
                const mid = low + (high - low) / 2;
                if (keys[mid] <= key) {
                    low = mid;
                } else {
                    high = mid;
                }
            }
            return if (keys[low] == key) low else null;
        }
    };
}

const IterMode = enum { count, fill };
const GridLayout = enum { dense, sparse };

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
    if (cell_list.isSparse()) {
        processNeighborPairsIn(T, mode, .sparse, cell_list, positions, radii, probe_radius, counts, neighbor_indices, offsets);
    } else {
        processNeighborPairsIn(T, mode, .dense, cell_list, positions, radii, probe_radius, counts, neighbor_indices, offsets);
    }
}

/// `processNeighborPairs` for one grid layout.
///
/// The sparse layout visits only the occupied cells, in ascending key order, and looks the
/// neighboring cells up by key. That is the order in which the dense layout reaches the same
/// cells, so both layouts report the pairs in the same order.
fn processNeighborPairsIn(
    comptime T: type,
    comptime mode: IterMode,
    comptime layout: GridLayout,
    cell_list: anytype,
    positions: []const Vec3Gen(T),
    radii: []const T,
    probe_radius: T,
    counts: []u32,
    neighbor_indices: []u32,
    offsets: []const u32,
) void {
    const n_cells = if (layout == .sparse)
        cell_list.cell_keys.len
    else
        cell_list.nx * cell_list.ny * cell_list.nz;
    for (0..n_cells) |cell_idx| {
        const cell1_atoms = cell_list.getCellAtoms(cell_idx);
        if (cell1_atoms.len == 0) continue;

        const coords = if (layout == .sparse)
            cell_list.getSparseCellCoords(cell_idx)
        else
            cell_list.getCellCoords(cell_idx);
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
                    const ncell = if (layout == .sparse)
                        cell_list.findSparseCell(cix + cdx, ciy + cdy, ciz + cdz, cell_idx)
                    else
                        cell_list.getCellIndexFromCoords(cix + cdx, ciy + cdy, ciz + cdz);
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
            return initWithDenseLimit(allocator, positions, radii, probe_radius, maxDenseCells(positions.len));
        }

        /// Same as `init` with an explicit limit on the size of a dense grid.
        /// The limit selects the layout of the grid only; the result does not depend on it.
        fn initWithDenseLimit(
            allocator: Allocator,
            positions: []const Vec,
            radii: []const T,
            probe_radius: T,
            max_dense_cells: usize,
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

            var cell_list = try CellListT.initWithDenseLimit(allocator, positions, cell_size, max_dense_cells);
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

// Grid bounds and sparse layout (issue #428)

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

test "fitGrid - a compact box is dense with the requested cell size" {
    inline for (.{ f32, f64 }) |T| {
        const shape = try fitGrid(T, .{ 64.5, 33.0, 6.6 }, 6.6, maxDenseCells(1));
        try std.testing.expect(shape.dense);
        try std.testing.expectEqual(@as(T, 6.6), shape.cell_size);
        try std.testing.expectEqual(@as(usize, 10), shape.nx);
        try std.testing.expectEqual(@as(usize, 5), shape.ny);
        try std.testing.expectEqual(@as(usize, 1), shape.nz);

        // A degenerate box still gets one cell per axis
        const point = try fitGrid(T, .{ 0.0, 0.0, 0.0 }, 6.6, 1);
        try std.testing.expect(point.dense);
        try std.testing.expectEqual(@as(usize, 1), point.nx * point.ny * point.nz);
    }
}

test "fitGrid - a box beyond the dense limit is sparse with the same cells" {
    inline for (.{ f32, f64 }) |T| {
        // Two atoms 2000 Å apart on the diagonal: 325^3 cells, 275 MB when stored densely
        const shape = try fitGrid(T, .{ 2012.4, 2012.4, 2012.4 }, 6.2, maxDenseCells(2));
        try std.testing.expect(!shape.dense);
        try std.testing.expectEqual(@as(T, 6.2), shape.cell_size);
        try std.testing.expectEqual(@as(usize, 325), shape.nx);
        try std.testing.expectEqual(@as(usize, 325), shape.ny);
        try std.testing.expectEqual(@as(usize, 325), shape.nz);

        // The limit decides the layout, nothing else
        const dense = try fitGrid(T, .{ 2012.4, 2012.4, 2012.4 }, 6.2, 325 * 325 * 325);
        try std.testing.expect(dense.dense);
        try std.testing.expectEqual(shape.cell_size, dense.cell_size);
        try std.testing.expectEqual(shape.nx, dense.nx);
    }
}

test "fitGrid - cells grow only when an axis exceeds the key range" {
    const extents = [_][3]f64{
        .{ 3.4e10, 3.4e10, 10.4 }, // flat: issue #428 overflow input, 2^32 cells per axis
        .{ 1.0e9, 12.4, 12.4 }, // linear
        .{ 1.0e30, 1.0e30, 1.0e30 },
        .{ 1.0e300, 5.0, 1.0e-300 },
        .{ 123.0, 4567.0, 8.9e12 },
    };

    inline for (.{ f32, f64 }) |T| {
        for (extents) |extent_f64| {
            if (T == f32 and extent_f64[0] > std.math.floatMax(f32)) continue;
            const extent = [3]T{ @floatCast(extent_f64[0]), @floatCast(extent_f64[1]), @floatCast(extent_f64[2]) };
            for ([_]T{ 6.2, 1.0e-30 }) |min_cell_size| {
                const shape = try fitGrid(T, extent, min_cell_size, maxDenseCells(2));
                try std.testing.expect(!shape.dense);
                try std.testing.expect(std.math.isFinite(shape.cell_size));
                try std.testing.expect(shape.cell_size > min_cell_size);

                // The widest axis uses the whole key range; every axis covers its extent
                const counts = [3]usize{ shape.nx, shape.ny, shape.nz };
                try std.testing.expectEqual(@as(usize, max_axis_cells), @max(shape.nx, shape.ny, shape.nz));
                for (counts, extent) |count, e| {
                    try std.testing.expect(count >= 1 and count <= max_axis_cells);
                    try std.testing.expect(@as(T, @floatFromInt(count)) * shape.cell_size >= e);
                }
            }
        }

        // 2^21 cells on an axis still fit: the cell size is kept
        const at_limit = try fitGrid(T, .{ 4.0 * max_axis_cells, 8.0, 8.0 }, 4.0, maxDenseCells(2));
        try std.testing.expectEqual(@as(T, 4.0), at_limit.cell_size);
        try std.testing.expectEqual(@as(usize, max_axis_cells), at_limit.nx);
    }
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

test "CellList - far-apart atoms use the sparse layout" {
    const allocator = std.testing.allocator;

    // Issue #428: 2^32 x 2^32 x 3 cells of 8 Å; the product used to wrap around.
    const positions = &[_]Vec3{
        Vec3{ .x = 0.0, .y = 0.0, .z = 0.0 },
        Vec3{ .x = 34359738352.0, .y = 34359738352.0, .z = 0.0 },
    };

    var cell_list = try CellList.init(allocator, positions, 8.0);
    defer cell_list.deinit();

    try std.testing.expect(cell_list.isSparse());
    try std.testing.expectEqual(@as(usize, 2), cell_list.cell_keys.len);
    try std.testing.expect(cell_list.cell_keys[0] < cell_list.cell_keys[1]);
    try std.testing.expectEqualSlices(u32, &.{ 0, 1, 2 }, cell_list.cell_offsets);
    try std.testing.expectEqualSlices(u32, &.{ 0, 1 }, cell_list.atom_indices);
    try std.testing.expect(cell_list.nx <= max_axis_cells and cell_list.ny <= max_axis_cells);

    // Occupied cells are found by their coordinates, empty and out-of-range cells are not
    for (0..2) |slot| {
        const coords = cell_list.getSparseCellCoords(slot);
        const ix: i64 = @intCast(coords.ix);
        const iy: i64 = @intCast(coords.iy);
        const iz: i64 = @intCast(coords.iz);
        try std.testing.expectEqual(@as(?usize, slot), cell_list.findSparseCell(ix, iy, iz, 0));
        try std.testing.expectEqual(@as(?usize, slot), cell_list.findSparseCell(ix, iy, iz, slot));
        // A cell before the first slot searched is not reported
        if (slot == 0) try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(ix, iy, iz, 1));
    }
    try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(5, 5, 0, 0));
    try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(-1, 0, 0, 0));
    try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(0, 0, @intCast(cell_list.nz), 0));
}

test "CellList - sparse lookup finds every occupied cell from every earlier slot" {
    const allocator = std.testing.allocator;

    var positions: [400]Vec3 = undefined;
    var radii: [400]f64 = undefined;
    const centers = [_]Vec3{
        .{ .x = 0.0, .y = 0.0, .z = 0.0 },
        .{ .x = 300.0, .y = 20.0, .z = -150.0 },
        .{ .x = -40.0, .y = 500.0, .z = 90.0 },
    };
    fillClusters(f64, 3, &centers, 45.0, &positions, &radii);

    var cell_list = try CellList.initWithDenseLimit(allocator, &positions, 6.6, 0);
    defer cell_list.deinit();
    const n_slots = cell_list.cell_keys.len;
    try std.testing.expect(n_slots > 100);

    for (0..n_slots) |slot| {
        const coords = cell_list.getSparseCellCoords(slot);
        const ix: i64 = @intCast(coords.ix);
        const iy: i64 = @intCast(coords.iy);
        const iz: i64 = @intCast(coords.iz);
        // Starting points at every distance exercise the doubling steps and the bisection
        var first_slot: usize = 0;
        while (first_slot <= slot) : (first_slot += 1 + first_slot / 3) {
            try std.testing.expectEqual(@as(?usize, slot), cell_list.findSparseCell(ix, iy, iz, first_slot));
        }
        try std.testing.expectEqual(@as(?usize, slot), cell_list.findSparseCell(ix, iy, iz, slot));
        if (slot + 1 < n_slots) {
            try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(ix, iy, iz, slot + 1));
        }
    }

    // A cell between two occupied ones that holds no atoms is not found
    var n_empty: usize = 0;
    for (0..n_slots - 1) |slot| {
        if (cell_list.cell_keys[slot] + 1 == cell_list.cell_keys[slot + 1]) continue;
        const coords = cell_list.getSparseCellCoords(slot);
        if (coords.ix + 1 >= cell_list.nx) continue;
        const ix: i64 = @intCast(coords.ix + 1);
        try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(ix, @intCast(coords.iy), @intCast(coords.iz), 0));
        try std.testing.expectEqual(@as(?usize, null), cell_list.findSparseCell(ix, @intCast(coords.iy), @intCast(coords.iz), slot));
        n_empty += 1;
    }
    try std.testing.expect(n_empty > 10);
}

test "CellList - compact structure uses the dense layout" {
    const allocator = std.testing.allocator;

    var positions: [200]Vec3 = undefined;
    var radii: [200]f64 = undefined;
    fillClusters(f64, 1, &.{.{ .x = 10.0, .y = -5.0, .z = 3.0 }}, 40.0, &positions, &radii);

    var cell_list = try CellList.init(allocator, &positions, 6.6);
    defer cell_list.deinit();

    try std.testing.expect(!cell_list.isSparse());
    try std.testing.expectEqual(@as(f64, 6.6), cell_list.cell_size);
    try std.testing.expectEqual(cell_list.nx * cell_list.ny * cell_list.nz + 1, cell_list.cell_offsets.len);
}

test "CellList - sparse layout stores the occupied cells of the dense layout" {
    const allocator = std.testing.allocator;

    inline for (.{ f32, f64 }) |T| {
        var positions: [200]Vec3Gen(T) = undefined;
        var radii: [200]T = undefined;
        fillClusters(T, 2, &.{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 90.0, .y = -40.0, .z = 25.0 } }, 30.0, &positions, &radii);

        var dense = try CellListGen(T).initWithDenseLimit(allocator, &positions, 6.6, std.math.maxInt(usize));
        defer dense.deinit();
        var sparse = try CellListGen(T).initWithDenseLimit(allocator, &positions, 6.6, 0);
        defer sparse.deinit();
        try std.testing.expect(!dense.isSparse());
        try std.testing.expect(sparse.isSparse());

        // Same cells in the same order, each with the same atoms in the same order
        try std.testing.expectEqualSlices(u32, dense.atom_indices, sparse.atom_indices);
        var slot: usize = 0;
        for (0..dense.nx * dense.ny * dense.nz) |cell_idx| {
            const atoms = dense.getCellAtoms(cell_idx);
            if (atoms.len == 0) continue;
            try std.testing.expectEqualSlices(u32, atoms, sparse.getCellAtoms(slot));
            const dense_coords = dense.getCellCoords(cell_idx);
            const sparse_coords = sparse.getSparseCellCoords(slot);
            try std.testing.expectEqual(dense_coords.ix, sparse_coords.ix);
            try std.testing.expectEqual(dense_coords.iy, sparse_coords.iy);
            try std.testing.expectEqual(dense_coords.iz, sparse_coords.iz);
            slot += 1;
        }
        try std.testing.expectEqual(sparse.cell_keys.len, slot);
    }
}

test "NeighborList - far-apart atoms have no neighbors" {
    const allocator = std.testing.allocator;

    inline for (.{ f32, f64 }) |T| {
        const cases = [_][2]Vec3Gen(T){
            // Issue #428: cell count overflow
            .{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 34359738352.0, .y = 34359738352.0, .z = 0.0 } },
            // Issue #428: hundreds of MB for two atoms
            .{ .{ .x = 0.0, .y = 0.0, .z = 0.0 }, .{ .x = 2000.0, .y = 2000.0, .z = 2000.0 } },
            // Sentinel coordinate
            .{ .{ .x = 1.0, .y = 2.0, .z = 3.0 }, .{ .x = 9999.999, .y = 9999.999, .z = 9999.999 } },
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

    // Atoms that touch stay neighbors however far from the origin they are
    {
        const positions = &[_]Vec3{
            Vec3{ .x = 1.0e12, .y = 1.0e12, .z = 1.0e12 },
            Vec3{ .x = 1.0e12 + 3.0, .y = 1.0e12, .z = 1.0e12 },
            Vec3{ .x = -1.0e12, .y = 0.0, .z = 0.0 },
        };
        var neighbor_list = try NeighborList.init(allocator, positions, &.{ 1.7, 1.7, 1.7 }, 1.4);
        defer neighbor_list.deinit();
        try expectMatchesBruteForce(f64, neighbor_list, positions, &.{ 1.7, 1.7, 1.7 }, 1.4);
        try std.testing.expectEqualSlices(u32, &.{1}, neighbor_list.getNeighbors(0));
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

test "NeighborList - sparse and dense layouts give identical lists" {
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
        var positions: [240]Vec = undefined;
        var radii: [240]T = undefined;

        for (layouts, 0..) |centers, layout_idx| {
            fillClusters(T, 1000 + layout_idx, centers, 18.0, &positions, &radii);

            var dense = try NeighborListGen(T).initWithDenseLimit(allocator, &positions, &radii, 1.4, std.math.maxInt(usize));
            defer dense.deinit();
            try expectMatchesBruteForce(T, dense, &positions, &radii, 1.4);

            // Same neighbors in the same order, whatever the limit
            for ([_]usize{ 0, 1, 27, 1000 }) |max_dense_cells| {
                var neighbor_list = try NeighborListGen(T).initWithDenseLimit(allocator, &positions, &radii, 1.4, max_dense_cells);
                defer neighbor_list.deinit();
                try std.testing.expectEqualSlices(u32, dense.offsets, neighbor_list.offsets);
                try std.testing.expectEqualSlices(u32, dense.neighbor_indices, neighbor_list.neighbor_indices);
            }
        }
    }
}

test "NeighborList - symmetry with the sparse layout" {
    const allocator = std.testing.allocator;

    // Two clusters far enough apart for the sparse layout
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

test "NeighborList - a stray distant atom does not change the other neighbor lists" {
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

        // The stray atom extends the grid without moving its origin or changing the cell
        // size, so every other atom keeps its neighbors in the same order.
        try std.testing.expectEqual(@as(usize, 0), with_stray.getNeighbors(n_compact).len);
        for (0..n_compact) |i| {
            try std.testing.expectEqualSlices(u32, compact.getNeighbors(i), with_stray.getNeighbors(i));
        }
    }
}

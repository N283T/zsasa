const std = @import("std");
const types = @import("types.zig");

const Vec3 = types.Vec3;
const Vec3Gen = types.Vec3Gen;
const Allocator = std.mem.Allocator;

/// Compute cell index from position coordinates
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
        pub fn init(
            allocator: Allocator,
            positions: []const Vec,
            cell_size: T,
        ) !Self {
            if (positions.len == 0) return error.NoAtoms;
            if (cell_size <= 0.0) return error.InvalidCellSize;

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
            x_min -= cell_size;
            y_min -= cell_size;
            z_min -= cell_size;
            x_max += cell_size;
            y_max += cell_size;
            z_max += cell_size;

            // Calculate grid dimensions (minimum 1 cell)
            const nx = @max(1, @as(usize, @intFromFloat(@ceil((x_max - x_min) / cell_size))));
            const ny = @max(1, @as(usize, @intFromFloat(@ceil((y_max - y_min) / cell_size))));
            const nz = @max(1, @as(usize, @intFromFloat(@ceil((z_max - z_min) / cell_size))));
            const n_cells = nx * ny * nz;

            // Pass 1: count atoms per cell
            const counts = try allocator.alloc(u32, n_cells);
            defer allocator.free(counts);
            @memset(counts, 0);

            for (positions) |pos| {
                counts[computeCellIndex(T, pos.x, pos.y, pos.z, x_min, y_min, z_min, cell_size, nx, ny, nz)] += 1;
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
            return computeCellIndex(T, pos.x, pos.y, pos.z, self.x_min, self.y_min, self.z_min, self.cell_size, self.nx, self.ny, self.nz);
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
        var dz: i64 = -1;
        while (dz <= 1) : (dz += 1) {
            var dy: i64 = -1;
            while (dy <= 1) : (dy += 1) {
                var dx: i64 = -1;
                while (dx <= 1) : (dx += 1) {
                    const ncell = cell_list.getCellIndexFromCoords(cix + dx, ciy + dy, ciz + dz);
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
                                const ddx = pi.x - pj.x;
                                const ddy = pi.y - pj.y;
                                const ddz = pi.z - pj.z;
                                const dist_sq = ddx * ddx + ddy * ddy + ddz * ddz;

                                const cutoff = radii[ai] + radii[aj] + 2.0 * probe_radius;

                                if (dist_sq < cutoff * cutoff) {
                                    if (mode == .count) {
                                        counts[ai] += 1;
                                        counts[aj] += 1;
                                    } else {
                                        neighbor_indices[offsets[ai] + counts[ai]] = aj;
                                        counts[ai] += 1;
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
            const n_atoms = positions.len;
            if (n_atoms == 0) return error.NoAtoms;

            // Find maximum radius and validate
            var max_radius: T = 0.0;
            for (radii) |r| {
                if (r < 0.0) return error.InvalidRadius;
                max_radius = @max(max_radius, r);
            }

            const cell_size = 2.0 * (max_radius + probe_radius);

            var cell_list = try CellListT.init(allocator, positions, cell_size);
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

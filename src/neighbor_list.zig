const std = @import("std");
const types = @import("types.zig");

const Vec3 = types.Vec3;
const Allocator = std.mem.Allocator;

/// Cell in spatial hash grid containing atom indices
pub const Cell = struct {
    atoms: std.ArrayListUnmanaged(u32),

    pub fn init() Cell {
        return Cell{
            .atoms = .{},
        };
    }

    pub fn deinit(self: *Cell, allocator: Allocator) void {
        self.atoms.deinit(allocator);
    }

    pub fn append(self: *Cell, allocator: Allocator, value: u32) !void {
        try self.atoms.append(allocator, value);
    }
};

/// Spatial hash grid for efficient neighbor lookups
pub const CellList = struct {
    cells: []Cell,
    nx: usize,
    ny: usize,
    nz: usize,
    cell_size: f64,
    x_min: f64,
    y_min: f64,
    z_min: f64,
    allocator: Allocator,

    /// Build spatial grid from atom positions
    /// cell_size should be >= 2 * (max_radius + probe_radius) for correctness
    pub fn init(
        allocator: Allocator,
        positions: []const Vec3,
        cell_size: f64,
    ) !CellList {
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

        // Add small padding to avoid edge cases
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

        // Allocate cells
        const cells = try allocator.alloc(Cell, n_cells);
        errdefer allocator.free(cells);

        for (cells) |*cell| {
            cell.* = Cell.init();
        }

        var cell_list = CellList{
            .cells = cells,
            .nx = nx,
            .ny = ny,
            .nz = nz,
            .cell_size = cell_size,
            .x_min = x_min,
            .y_min = y_min,
            .z_min = z_min,
            .allocator = allocator,
        };

        // Assign atoms to cells
        for (positions, 0..) |pos, i| {
            const cell_idx = cell_list.getCellIndex(pos);
            try cell_list.cells[cell_idx].append(allocator, @intCast(i));
        }

        return cell_list;
    }

    pub fn deinit(self: *CellList) void {
        for (self.cells) |*cell| {
            cell.deinit(self.allocator);
        }
        self.allocator.free(self.cells);
    }

    /// Get cell index for a position
    fn getCellIndex(self: CellList, pos: Vec3) usize {
        const ix = @as(usize, @intFromFloat(@max(0.0, (pos.x - self.x_min) / self.cell_size)));
        const iy = @as(usize, @intFromFloat(@max(0.0, (pos.y - self.y_min) / self.cell_size)));
        const iz = @as(usize, @intFromFloat(@max(0.0, (pos.z - self.z_min) / self.cell_size)));

        // Clamp to valid range
        const cix = @min(ix, self.nx - 1);
        const ciy = @min(iy, self.ny - 1);
        const ciz = @min(iz, self.nz - 1);

        return ciz * self.nx * self.ny + ciy * self.nx + cix;
    }

    /// Get cell coordinates from index
    fn getCellCoords(self: CellList, idx: usize) struct { ix: usize, iy: usize, iz: usize } {
        const iz = idx / (self.nx * self.ny);
        const remainder = idx % (self.nx * self.ny);
        const iy = remainder / self.nx;
        const ix = remainder % self.nx;
        return .{ .ix = ix, .iy = iy, .iz = iz };
    }

    /// Get cell index from coordinates (returns null if out of bounds)
    fn getCellIndexFromCoords(self: CellList, ix: i64, iy: i64, iz: i64) ?usize {
        if (ix < 0 or iy < 0 or iz < 0) return null;
        const uix = @as(usize, @intCast(ix));
        const uiy = @as(usize, @intCast(iy));
        const uiz = @as(usize, @intCast(iz));
        if (uix >= self.nx or uiy >= self.ny or uiz >= self.nz) return null;
        return uiz * self.nx * self.ny + uiy * self.nx + uix;
    }
};

/// Pre-computed neighbor list for efficient SASA calculation
pub const NeighborList = struct {
    /// neighbors[i] contains indices of atoms that could potentially bury atom i's surface
    neighbors: []std.ArrayListUnmanaged(u32),
    allocator: Allocator,

    /// Build neighbor list from atom positions and radii
    /// Two atoms i, j are neighbors if distance(i, j) < r[i] + r[j] + 2*probe_radius
    pub fn init(
        allocator: Allocator,
        positions: []const Vec3,
        radii: []const f64,
        probe_radius: f64,
    ) !NeighborList {
        const n_atoms = positions.len;
        if (n_atoms == 0) return error.NoAtoms;

        // Find maximum radius and validate
        var max_radius: f64 = 0.0;
        for (radii) |r| {
            if (r < 0.0) return error.InvalidRadius;
            max_radius = @max(max_radius, r);
        }

        // Cell size: 2 * (max_radius + probe_radius)
        // This ensures that neighbors can only be in adjacent cells
        const cell_size = 2.0 * (max_radius + probe_radius);

        // Build cell list
        var cell_list = try CellList.init(allocator, positions, cell_size);
        defer cell_list.deinit();

        // Allocate neighbor lists
        const neighbors = try allocator.alloc(std.ArrayListUnmanaged(u32), n_atoms);
        errdefer {
            for (neighbors) |*list| {
                list.deinit(allocator);
            }
            allocator.free(neighbors);
        }

        for (neighbors) |*list| {
            list.* = .{};
        }

        // For each cell, check atoms in same cell and 26 neighboring cells
        const n_cells = cell_list.nx * cell_list.ny * cell_list.nz;
        for (0..n_cells) |cell_idx| {
            const coords = cell_list.getCellCoords(cell_idx);
            const ix = @as(i64, @intCast(coords.ix));
            const iy = @as(i64, @intCast(coords.iy));
            const iz = @as(i64, @intCast(coords.iz));

            // Check all 27 neighboring cells (including self)
            var dz: i64 = -1;
            while (dz <= 1) : (dz += 1) {
                var dy: i64 = -1;
                while (dy <= 1) : (dy += 1) {
                    var dx: i64 = -1;
                    while (dx <= 1) : (dx += 1) {
                        const neighbor_idx = cell_list.getCellIndexFromCoords(ix + dx, iy + dy, iz + dz);
                        if (neighbor_idx) |nidx| {
                            // Only process cell pairs where cell_idx <= nidx to avoid duplicates
                            // This ensures (A,B) is processed but not (B,A)
                            if (nidx < cell_idx) continue;

                            // Check all pairs between this cell and neighbor cell
                            try addNeighborPairs(
                                allocator,
                                &cell_list.cells[cell_idx],
                                &cell_list.cells[nidx],
                                positions,
                                radii,
                                probe_radius,
                                neighbors,
                                cell_idx == nidx, // same_cell flag
                            );
                        }
                    }
                }
            }
        }

        return NeighborList{
            .neighbors = neighbors,
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *NeighborList) void {
        for (self.neighbors) |*list| {
            list.deinit(self.allocator);
        }
        self.allocator.free(self.neighbors);
    }

    /// Get neighbors for atom i
    pub fn getNeighbors(self: NeighborList, atom_idx: usize) []const u32 {
        return self.neighbors[atom_idx].items;
    }
};

/// Add neighbor pairs between atoms in two cells
fn addNeighborPairs(
    allocator: Allocator,
    cell1: *const Cell,
    cell2: *const Cell,
    positions: []const Vec3,
    radii: []const f64,
    probe_radius: f64,
    neighbors: []std.ArrayListUnmanaged(u32),
    same_cell: bool,
) !void {
    for (cell1.atoms.items) |i| {
        for (cell2.atoms.items) |j| {
            if (i == j) continue;

            // Check if already added (for same cell, avoid duplicate pairs)
            if (same_cell and j <= i) continue;

            const pos_i = positions[i];
            const pos_j = positions[j];
            const r_i = radii[i];
            const r_j = radii[j];

            // Calculate distance squared
            const dx = pos_i.x - pos_j.x;
            const dy = pos_i.y - pos_j.y;
            const dz = pos_i.z - pos_j.z;
            const dist_sq = dx * dx + dy * dy + dz * dz;

            // Cutoff: r[i] + r[j] + 2*probe_radius
            const cutoff = r_i + r_j + 2.0 * probe_radius;
            const cutoff_sq = cutoff * cutoff;

            if (dist_sq < cutoff_sq) {
                // Add symmetric neighbors
                try neighbors[i].append(allocator, j);
                try neighbors[j].append(allocator, i);
            }
        }
    }
}

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

const std = @import("std");

/// 3D vector/point representation
pub const Vec3 = struct {
    x: f64,
    y: f64,
    z: f64,

    /// Add two vectors
    pub fn add(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x + other.x,
            .y = self.y + other.y,
            .z = self.z + other.z,
        };
    }

    /// Subtract two vectors
    pub fn sub(self: Vec3, other: Vec3) Vec3 {
        return Vec3{
            .x = self.x - other.x,
            .y = self.y - other.y,
            .z = self.z - other.z,
        };
    }

    /// Scale vector by scalar
    pub fn scale(self: Vec3, scalar: f64) Vec3 {
        return Vec3{
            .x = self.x * scalar,
            .y = self.y * scalar,
            .z = self.z * scalar,
        };
    }

    /// Calculate vector length (magnitude)
    pub fn length(self: Vec3) f64 {
        return @sqrt(self.x * self.x + self.y * self.y + self.z * self.z);
    }

    /// Calculate dot product
    pub fn dot(self: Vec3, other: Vec3) f64 {
        return self.x * other.x + self.y * other.y + self.z * other.z;
    }
};

/// Input data structure for SASA calculation
pub const AtomInput = struct {
    x: []const f64,
    y: []const f64,
    z: []const f64,
    r: []const f64,
    /// Residue names (e.g., "ALA", "GLY") - optional, for classifier
    residue: ?[]const []const u8 = null,
    /// Atom names (e.g., "CA", "CB") - optional, for classifier
    atom_name: ?[]const []const u8 = null,
    allocator: std.mem.Allocator,

    /// Get number of atoms
    pub fn atomCount(self: AtomInput) usize {
        return self.x.len;
    }

    /// Check if residue/atom names are available for classification
    pub fn hasClassificationInfo(self: AtomInput) bool {
        return self.residue != null and self.atom_name != null;
    }

    /// Free allocated memory
    pub fn deinit(self: *AtomInput) void {
        self.allocator.free(self.x);
        self.allocator.free(self.y);
        self.allocator.free(self.z);
        self.allocator.free(self.r);
        if (self.residue) |res| {
            for (res) |s| {
                self.allocator.free(s);
            }
            self.allocator.free(res);
        }
        if (self.atom_name) |names| {
            for (names) |s| {
                self.allocator.free(s);
            }
            self.allocator.free(names);
        }
    }
};

/// Output data structure for SASA calculation
pub const SasaResult = struct {
    total_area: f64,
    atom_areas: []f64,
    allocator: std.mem.Allocator,

    /// Free allocated memory
    pub fn deinit(self: *SasaResult) void {
        self.allocator.free(self.atom_areas);
    }
};

/// Configuration parameters for SASA calculation
pub const Config = struct {
    /// Number of test points per atom
    n_points: u32 = 100,
    /// Water probe radius in Angstroms
    probe_radius: f64 = 1.4,
};

// Tests
test "Vec3 add" {
    const v1 = Vec3{ .x = 1.0, .y = 2.0, .z = 3.0 };
    const v2 = Vec3{ .x = 4.0, .y = 5.0, .z = 6.0 };
    const result = v1.add(v2);

    try std.testing.expectEqual(@as(f64, 5.0), result.x);
    try std.testing.expectEqual(@as(f64, 7.0), result.y);
    try std.testing.expectEqual(@as(f64, 9.0), result.z);
}

test "Vec3 sub" {
    const v1 = Vec3{ .x = 4.0, .y = 5.0, .z = 6.0 };
    const v2 = Vec3{ .x = 1.0, .y = 2.0, .z = 3.0 };
    const result = v1.sub(v2);

    try std.testing.expectEqual(@as(f64, 3.0), result.x);
    try std.testing.expectEqual(@as(f64, 3.0), result.y);
    try std.testing.expectEqual(@as(f64, 3.0), result.z);
}

test "Vec3 scale" {
    const v = Vec3{ .x = 1.0, .y = 2.0, .z = 3.0 };
    const result = v.scale(2.0);

    try std.testing.expectEqual(@as(f64, 2.0), result.x);
    try std.testing.expectEqual(@as(f64, 4.0), result.y);
    try std.testing.expectEqual(@as(f64, 6.0), result.z);
}

test "Vec3 length" {
    const v = Vec3{ .x = 3.0, .y = 4.0, .z = 0.0 };
    const result = v.length();

    try std.testing.expectEqual(@as(f64, 5.0), result);
}

test "Vec3 dot product" {
    const v1 = Vec3{ .x = 1.0, .y = 2.0, .z = 3.0 };
    const v2 = Vec3{ .x = 4.0, .y = 5.0, .z = 6.0 };
    const result = v1.dot(v2);

    // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    try std.testing.expectEqual(@as(f64, 32.0), result);
}

test "AtomInput atomCount" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 3);
    const y = try allocator.alloc(f64, 3);
    const z = try allocator.alloc(f64, 3);
    const r = try allocator.alloc(f64, 3);

    var input = AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .allocator = allocator,
    };
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
}

test "SasaResult deinit" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f64, 5);

    var result = SasaResult{
        .total_area = 100.0,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
    defer result.deinit();

    try std.testing.expectEqual(@as(f64, 100.0), result.total_area);
    try std.testing.expectEqual(@as(usize, 5), result.atom_areas.len);
}

test "Config default values" {
    const config = Config{};

    try std.testing.expectEqual(@as(u32, 100), config.n_points);
    try std.testing.expectEqual(@as(f64, 1.4), config.probe_radius);
}

test "Config custom values" {
    const config = Config{
        .n_points = 200,
        .probe_radius = 1.8,
    };

    try std.testing.expectEqual(@as(u32, 200), config.n_points);
    try std.testing.expectEqual(@as(f64, 1.8), config.probe_radius);
}

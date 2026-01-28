const std = @import("std");

/// Precision mode for SASA calculation
pub const Precision = enum {
    f32,
    f64,

    pub fn fromString(s: []const u8) ?Precision {
        if (std.mem.eql(u8, s, "f32") or std.mem.eql(u8, s, "single")) {
            return .f32;
        } else if (std.mem.eql(u8, s, "f64") or std.mem.eql(u8, s, "double")) {
            return .f64;
        }
        return null;
    }
};

/// Generic 3D vector/point representation
pub fn Vec3Gen(comptime T: type) type {
    return struct {
        const Self = @This();

        x: T,
        y: T,
        z: T,

        /// Add two vectors
        pub fn add(self: Self, other: Self) Self {
            return Self{
                .x = self.x + other.x,
                .y = self.y + other.y,
                .z = self.z + other.z,
            };
        }

        /// Subtract two vectors
        pub fn sub(self: Self, other: Self) Self {
            return Self{
                .x = self.x - other.x,
                .y = self.y - other.y,
                .z = self.z - other.z,
            };
        }

        /// Scale vector by scalar
        pub fn scale(self: Self, scalar: T) Self {
            return Self{
                .x = self.x * scalar,
                .y = self.y * scalar,
                .z = self.z * scalar,
            };
        }

        /// Calculate vector length (magnitude)
        pub fn length(self: Self) T {
            return @sqrt(self.x * self.x + self.y * self.y + self.z * self.z);
        }

        /// Calculate dot product
        pub fn dot(self: Self, other: Self) T {
            return self.x * other.x + self.y * other.y + self.z * other.z;
        }

        /// Convert from f64 Vec3
        pub fn fromF64(v: Vec3Gen(f64)) Self {
            if (T == f64) {
                return v;
            }
            return Self{
                .x = @floatCast(v.x),
                .y = @floatCast(v.y),
                .z = @floatCast(v.z),
            };
        }

        /// Convert to f64 Vec3
        pub fn toF64(self: Self) Vec3Gen(f64) {
            if (T == f64) {
                return self;
            }
            return Vec3Gen(f64){
                .x = @floatCast(self.x),
                .y = @floatCast(self.y),
                .z = @floatCast(self.z),
            };
        }
    };
}

/// Default Vec3 type (f64 for backward compatibility)
pub const Vec3 = Vec3Gen(f64);
/// Single precision Vec3
pub const Vec3f32 = Vec3Gen(f32);

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
    /// Element atomic numbers (e.g., 6=C, 7=N, 8=O, 20=Ca) - optional
    /// Used to disambiguate atom names like "CA" (Carbon alpha vs Calcium)
    element: ?[]const u8 = null,
    /// Chain IDs (e.g., "A", "B") - optional, for per-chain analysis
    chain_id: ?[]const []const u8 = null,
    /// Residue sequence numbers - optional, for per-residue analysis
    residue_num: ?[]const i32 = null,
    /// Insertion codes (e.g., "", "A", "B") - optional, for per-residue analysis
    insertion_code: ?[]const []const u8 = null,
    allocator: std.mem.Allocator,

    /// Get number of atoms
    pub fn atomCount(self: AtomInput) usize {
        return self.x.len;
    }

    /// Check if residue/atom names are available for classification
    pub fn hasClassificationInfo(self: AtomInput) bool {
        return self.residue != null and self.atom_name != null;
    }

    /// Check if element info is available
    pub fn hasElementInfo(self: AtomInput) bool {
        return self.element != null;
    }

    /// Check if per-residue analysis info is available
    pub fn hasResidueInfo(self: AtomInput) bool {
        return self.residue != null and self.residue_num != null and
            self.chain_id != null and self.insertion_code != null;
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
        if (self.element) |elem| {
            self.allocator.free(elem);
        }
        if (self.chain_id) |chains| {
            for (chains) |s| {
                self.allocator.free(s);
            }
            self.allocator.free(chains);
        }
        if (self.residue_num) |nums| {
            self.allocator.free(nums);
        }
        if (self.insertion_code) |codes| {
            for (codes) |s| {
                self.allocator.free(s);
            }
            self.allocator.free(codes);
        }
    }
};

/// Generic output data structure for SASA calculation
pub fn SasaResultGen(comptime T: type) type {
    return struct {
        const Self = @This();

        total_area: T,
        atom_areas: []T,
        allocator: std.mem.Allocator,

        /// Free allocated memory
        pub fn deinit(self: *Self) void {
            self.allocator.free(self.atom_areas);
        }

        /// Convert to f64 result (for output compatibility)
        pub fn toF64(self: Self, allocator: std.mem.Allocator) !SasaResultGen(f64) {
            if (T == f64) {
                // For f64, just copy the slice reference
                const areas = try allocator.alloc(f64, self.atom_areas.len);
                @memcpy(areas, self.atom_areas);
                return SasaResultGen(f64){
                    .total_area = self.total_area,
                    .atom_areas = areas,
                    .allocator = allocator,
                };
            }
            const areas = try allocator.alloc(f64, self.atom_areas.len);
            for (self.atom_areas, 0..) |area, i| {
                areas[i] = @floatCast(area);
            }
            return SasaResultGen(f64){
                .total_area = @floatCast(self.total_area),
                .atom_areas = areas,
                .allocator = allocator,
            };
        }
    };
}

/// Default SasaResult type (f64 for backward compatibility)
pub const SasaResult = SasaResultGen(f64);
/// Single precision SasaResult
pub const SasaResultf32 = SasaResultGen(f32);

/// Generic configuration parameters for SASA calculation
pub fn ConfigGen(comptime T: type) type {
    return struct {
        /// Number of test points per atom
        n_points: u32 = 100,
        /// Water probe radius in Angstroms
        probe_radius: T = 1.4,
    };
}

/// Default Config type (f64 for backward compatibility)
pub const Config = ConfigGen(f64);
/// Single precision Config
pub const Configf32 = ConfigGen(f32);

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

// Tests for f32 types
test "Vec3f32 add" {
    const v1 = Vec3f32{ .x = 1.0, .y = 2.0, .z = 3.0 };
    const v2 = Vec3f32{ .x = 4.0, .y = 5.0, .z = 6.0 };
    const result = v1.add(v2);

    try std.testing.expectEqual(@as(f32, 5.0), result.x);
    try std.testing.expectEqual(@as(f32, 7.0), result.y);
    try std.testing.expectEqual(@as(f32, 9.0), result.z);
}

test "Vec3f32 conversion" {
    const v64 = Vec3{ .x = 1.5, .y = 2.5, .z = 3.5 };
    const v32 = Vec3f32.fromF64(v64);

    try std.testing.expectEqual(@as(f32, 1.5), v32.x);
    try std.testing.expectEqual(@as(f32, 2.5), v32.y);
    try std.testing.expectEqual(@as(f32, 3.5), v32.z);

    const back = v32.toF64();
    try std.testing.expectEqual(@as(f64, 1.5), back.x);
    try std.testing.expectEqual(@as(f64, 2.5), back.y);
    try std.testing.expectEqual(@as(f64, 3.5), back.z);
}

test "Precision fromString" {
    try std.testing.expectEqual(Precision.f32, Precision.fromString("f32").?);
    try std.testing.expectEqual(Precision.f32, Precision.fromString("single").?);
    try std.testing.expectEqual(Precision.f64, Precision.fromString("f64").?);
    try std.testing.expectEqual(Precision.f64, Precision.fromString("double").?);
    try std.testing.expectEqual(@as(?Precision, null), Precision.fromString("invalid"));
}

test "SasaResultf32 toF64" {
    const allocator = std.testing.allocator;

    const atom_areas = try allocator.alloc(f32, 3);
    atom_areas[0] = 10.5;
    atom_areas[1] = 20.5;
    atom_areas[2] = 30.5;

    var result32 = SasaResultf32{
        .total_area = 61.5,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
    defer result32.deinit();

    var result64 = try result32.toF64(allocator);
    defer result64.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 61.5), result64.total_area, 0.001);
    try std.testing.expectEqual(@as(usize, 3), result64.atom_areas.len);
    try std.testing.expectApproxEqAbs(@as(f64, 10.5), result64.atom_areas[0], 0.001);
}

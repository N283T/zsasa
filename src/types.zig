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

/// Epsilon values for floating-point comparisons.
/// Different use cases may require different tolerances.
pub fn Epsilon(comptime T: type) type {
    return struct {
        /// Default epsilon for general near-zero comparisons (e.g., distance checks)
        /// f32: 1e-6, f64: 1e-10
        pub const default: T = if (T == f32) 1e-6 else 1e-10;

        /// Stricter epsilon for trigonometric functions (e.g., atan2 near-zero)
        /// f32: 1e-7, f64: 1e-10
        pub const trig: T = if (T == f32) 1e-7 else 1e-10;
    };
}

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

/// Fixed-size string for atom/chain names (max 4 characters)
/// Eliminates per-atom string allocations
pub const FixedString4 = struct {
    data: [4]u8 = .{ 0, 0, 0, 0 },
    len: u8 = 0,

    /// Create from a slice (truncates if > 4 chars)
    pub fn fromSlice(s: []const u8) FixedString4 {
        var result = FixedString4{};
        const copy_len: u8 = @intCast(@min(s.len, 4));
        @memcpy(result.data[0..copy_len], s[0..copy_len]);
        result.len = copy_len;
        return result;
    }

    /// Get as slice
    pub fn slice(self: *const FixedString4) []const u8 {
        return self.data[0..self.len];
    }

    /// Check equality with another FixedString4
    pub fn eql(self: *const FixedString4, other: *const FixedString4) bool {
        return self.len == other.len and std.mem.eql(u8, self.slice(), other.slice());
    }

    /// Check equality with a slice
    pub fn eqlSlice(self: *const FixedString4, s: []const u8) bool {
        return self.len == s.len and std.mem.eql(u8, self.slice(), s);
    }
};

/// Fixed-size string for residue names (max 5 characters)
/// Supports mmCIF 5-character comp_id (extended CCD IDs)
pub const FixedString5 = struct {
    data: [5]u8 = .{ 0, 0, 0, 0, 0 },
    len: u8 = 0,

    /// Create from a slice (truncates if > 5 chars)
    pub fn fromSlice(s: []const u8) FixedString5 {
        var result = FixedString5{};
        const copy_len: u8 = @intCast(@min(s.len, 5));
        @memcpy(result.data[0..copy_len], s[0..copy_len]);
        result.len = copy_len;
        return result;
    }

    /// Get as slice
    pub fn slice(self: *const FixedString5) []const u8 {
        return self.data[0..self.len];
    }

    /// Check equality with another FixedString5
    pub fn eql(self: *const FixedString5, other: *const FixedString5) bool {
        return self.len == other.len and std.mem.eql(u8, self.slice(), other.slice());
    }

    /// Check equality with a slice
    pub fn eqlSlice(self: *const FixedString5, s: []const u8) bool {
        return self.len == s.len and std.mem.eql(u8, self.slice(), s);
    }
};

/// Input data structure for SASA calculation
pub const AtomInput = struct {
    x: []const f64,
    y: []const f64,
    z: []const f64,
    r: []f64,
    /// Residue names (e.g., "ALA", "GLY") - optional, for classifier
    /// Uses FixedString5 for mmCIF 5-character comp_id support
    residue: ?[]const FixedString5 = null,
    /// Atom names (e.g., "CA", "CB") - optional, for classifier
    atom_name: ?[]const FixedString4 = null,
    /// Element atomic numbers (e.g., 6=C, 7=N, 8=O, 20=Ca) - optional
    /// Used to disambiguate atom names like "CA" (Carbon alpha vs Calcium)
    element: ?[]const u8 = null,
    /// Chain IDs (e.g., "A", "B") - optional, for per-chain analysis
    chain_id: ?[]const FixedString4 = null,
    /// Extended chain IDs for formats that can exceed four characters (e.g. mmCIF).
    /// When present, this mirrors `chain_id` length and owns each string.
    /// Chain filtering, residue grouping and the output writers use this
    /// field when it is available; `chain_id` holds only the first four
    /// characters.
    chain_id_full: ?[]const []const u8 = null,
    /// Residue sequence numbers - optional, for per-residue analysis
    residue_num: ?[]const i32 = null,
    /// Insertion codes (e.g., "", "A", "B") - optional, for per-residue analysis
    insertion_code: ?[]const FixedString4 = null,
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

    /// Validate that required coordinate/radius arrays are aligned and finite.
    /// Radii must be non-negative; zero-radius points are accepted.
    pub fn validateFiniteAndRange(self: AtomInput) !void {
        const n_atoms = self.atomCount();
        if (self.y.len != n_atoms or self.z.len != n_atoms or self.r.len != n_atoms) {
            return error.InvalidInput;
        }

        for (0..n_atoms) |i| {
            if (!std.math.isFinite(self.x[i]) or
                !std.math.isFinite(self.y[i]) or
                !std.math.isFinite(self.z[i]) or
                !std.math.isFinite(self.r[i]) or
                self.r[i] < 0.0)
            {
                return error.InvalidInput;
            }
        }
    }

    /// Remove the atoms of the residues `exclusion` names and return how
    /// many were removed. Inputs without residue names are left unchanged.
    /// Removing every atom is `error.NoAtomsFound`.
    pub fn removeResidues(self: *AtomInput, exclusion: ResidueExclusion) !usize {
        const residues = self.residue orelse return 0;
        if (exclusion.isEmpty()) return 0;
        var n_removed: usize = 0;
        for (residues) |res| {
            if (exclusion.excludes(res.slice())) n_removed += 1;
        }
        if (n_removed == 0) return 0;
        const n = residues.len;
        if (n_removed == n) return error.NoAtomsFound;

        const indices = try self.allocator.alloc(usize, n - n_removed);
        defer self.allocator.free(indices);
        var k: usize = 0;
        for (residues, 0..) |res, i| {
            if (exclusion.excludes(res.slice())) continue;
            indices[k] = i;
            k += 1;
        }

        var out = AtomInput{ .x = &.{}, .y = &.{}, .z = &.{}, .r = &.{}, .allocator = self.allocator };
        errdefer out.deinit();
        out.x = try gather(f64, self.allocator, self.x, indices);
        out.y = try gather(f64, self.allocator, self.y, indices);
        out.z = try gather(f64, self.allocator, self.z, indices);
        out.r = try gather(f64, self.allocator, self.r, indices);
        out.residue = try gather(FixedString5, self.allocator, residues, indices);
        if (self.atom_name) |v| out.atom_name = try gather(FixedString4, self.allocator, v, indices);
        if (self.element) |v| out.element = try gather(u8, self.allocator, v, indices);
        if (self.chain_id) |v| out.chain_id = try gather(FixedString4, self.allocator, v, indices);
        if (self.residue_num) |v| out.residue_num = try gather(i32, self.allocator, v, indices);
        if (self.insertion_code) |v| out.insertion_code = try gather(FixedString4, self.allocator, v, indices);
        if (self.chain_id_full) |v| {
            const chains = try self.allocator.alloc([]const u8, indices.len);
            errdefer self.allocator.free(chains);
            for (indices, chains) |i, *dst| dst.* = try self.allocator.dupe(u8, v[i]);
            out.chain_id_full = chains;
        }

        self.deinit();
        self.* = out;
        return n_removed;
    }

    fn gather(comptime T: type, allocator: std.mem.Allocator, src: []const T, indices: []const usize) ![]T {
        const out = try allocator.alloc(T, indices.len);
        for (indices, out) |i, *dst| dst.* = src[i];
        return out;
    }

    /// Free allocated memory
    pub fn deinit(self: *AtomInput) void {
        self.allocator.free(self.x);
        self.allocator.free(self.y);
        self.allocator.free(self.z);
        self.allocator.free(self.r);
        // FixedString4 arrays: single allocation, no per-element free needed
        if (self.residue) |res| {
            self.allocator.free(res);
        }
        if (self.atom_name) |names| {
            self.allocator.free(names);
        }
        if (self.element) |elem| {
            self.allocator.free(elem);
        }
        if (self.chain_id) |chains| {
            self.allocator.free(chains);
        }
        if (self.chain_id_full) |chains| {
            for (chains) |chain| {
                self.allocator.free(chain);
            }
            self.allocator.free(chains);
        }
        if (self.residue_num) |nums| {
            self.allocator.free(nums);
        }
        if (self.insertion_code) |codes| {
            self.allocator.free(codes);
        }
    }
};

/// Residue names of water: the PDB names `HOH`, `DOD` (heavy water), `WAT`,
/// `H2O` and `D2O`, and the names MD force fields give their water models.
const water_residues = std.StaticStringMap(void).initComptime(.{
    .{"HOH"}, .{"DOD"}, .{"WAT"},  .{"H2O"},  .{"D2O"},
    .{"SOL"}, .{"TIP"}, .{"TIP3"}, .{"TIP4"}, .{"TIP5"},
    .{"T3P"}, .{"T4P"}, .{"T5P"},  .{"SPC"},  .{"SPCE"},
});

/// Whether a residue name is water (`water_residues`), ignoring surrounding
/// blanks and case.
pub fn isWaterResidue(name: []const u8) bool {
    const trimmed = std.mem.trim(u8, name, " ");
    if (trimmed.len == 0 or trimmed.len > 5) return false;
    var buf: [5]u8 = undefined;
    const upper = std.ascii.upperString(buf[0..trimmed.len], trimmed);
    return water_residues.has(upper);
}

/// Residues left out of a calculation: water unless `--include-water`, and
/// the names of `--exclude-residues` (comma-separated) or of the workflow
/// `exclude_residues` list. Names compare without surrounding blanks and
/// case.
pub const ResidueExclusion = struct {
    water: bool = true,
    /// Comma-separated residue names, as given to `--exclude-residues`
    names_csv: []const u8 = "",
    names: []const []const u8 = &.{},

    pub fn isEmpty(self: ResidueExclusion) bool {
        return !self.water and std.mem.trim(u8, self.names_csv, " ,").len == 0 and self.names.len == 0;
    }

    pub fn excludes(self: ResidueExclusion, residue: []const u8) bool {
        if (self.water and isWaterResidue(residue)) return true;
        const name = std.mem.trim(u8, residue, " ");
        var it = std.mem.tokenizeScalar(u8, self.names_csv, ',');
        while (it.next()) |listed| {
            if (sameResidueName(listed, name)) return true;
        }
        for (self.names) |listed| {
            if (sameResidueName(listed, name)) return true;
        }
        return false;
    }

    fn sameResidueName(listed: []const u8, name: []const u8) bool {
        const trimmed = std.mem.trim(u8, listed, " ");
        return trimmed.len > 0 and std.ascii.eqlIgnoreCase(trimmed, name);
    }
};

/// Atom centres and radii in the layout the SASA kernels read: one
/// `Vec3Gen(T)` per atom, so the three coordinates of a neighbour are fetched
/// together. A caller that keeps its coordinates in this form passes them
/// without a copy; `gatherPositions` fills one from separate x/y/z arrays.
pub fn AtomsGen(comptime T: type) type {
    return struct {
        const Self = @This();

        positions: []const Vec3Gen(T),
        radii: []const T,

        pub fn atomCount(self: Self) usize {
            return self.positions.len;
        }

        /// The checks of `AtomInput.validateFiniteAndRange`, made on the
        /// values in `T`: a coordinate that does not fit in `T` is rejected
        /// instead of being calculated as infinity.
        pub fn validateFiniteAndRange(self: Self) !void {
            if (self.radii.len != self.positions.len) return error.InvalidInput;
            for (self.positions, self.radii) |p, r| {
                if (!std.math.isFinite(p.x) or
                    !std.math.isFinite(p.y) or
                    !std.math.isFinite(p.z) or
                    !std.math.isFinite(r) or
                    r < 0.0)
                {
                    return error.InvalidInput;
                }
            }
        }
    };
}

/// Write the coordinates `x`, `y`, `z` (slices of f32 or f64) into `dst`,
/// converted to `T`. With `indices`, `dst[i]` takes atom `indices[i]`;
/// without, atom `i`. This is the one copy between coordinates kept as
/// separate arrays and the layout of `AtomsGen(T)`.
pub fn gatherPositions(
    comptime T: type,
    dst: []Vec3Gen(T),
    x: anytype,
    y: anytype,
    z: anytype,
    indices: ?[]const u32,
) !void {
    if (y.len != x.len or z.len != x.len) return error.InvalidInput;
    if (indices) |idx| {
        if (idx.len != dst.len) return error.InvalidInput;
        for (idx, dst) |src, *p| {
            if (src >= x.len) return error.InvalidInput;
            p.* = .{ .x = @floatCast(x[src]), .y = @floatCast(y[src]), .z = @floatCast(z[src]) };
        }
    } else {
        if (x.len != dst.len) return error.InvalidInput;
        for (dst, x, y, z) |*p, xi, yi, zi| {
            p.* = .{ .x = @floatCast(xi), .y = @floatCast(yi), .z = @floatCast(zi) };
        }
    }
}

/// `AtomsGen(T)` with storage of its own, for the entry points that take an
/// `AtomInput`.
pub fn AtomsBufferGen(comptime T: type) type {
    return struct {
        const Self = @This();

        positions: []Vec3Gen(T),
        radii: []T,
        allocator: std.mem.Allocator,

        /// Copy the coordinates and radii of `input`, converted to `T`.
        /// Values are not checked; validate `input` first.
        pub fn fromInput(allocator: std.mem.Allocator, input: AtomInput) !Self {
            const n_atoms = input.atomCount();
            if (input.r.len != n_atoms) return error.InvalidInput;
            const positions = try allocator.alloc(Vec3Gen(T), n_atoms);
            errdefer allocator.free(positions);
            try gatherPositions(T, positions, input.x, input.y, input.z, null);
            const radii = try allocator.alloc(T, n_atoms);
            for (radii, input.r) |*dst, r| dst.* = @floatCast(r);
            return .{ .positions = positions, .radii = radii, .allocator = allocator };
        }

        pub fn atoms(self: Self) AtomsGen(T) {
            return .{ .positions = self.positions, .radii = self.radii };
        }

        pub fn deinit(self: *Self) void {
            self.allocator.free(self.positions);
            self.allocator.free(self.radii);
        }
    };
}

/// Generic output data structure for SASA calculation
pub fn SasaResultGen(comptime T: type) type {
    return struct {
        const Self = @This();

        /// Sum of `atom_areas`, accumulated in f64 whatever `T` is. An f32
        /// accumulator loses about seven digits per addition and made the
        /// total depend on the number of threads.
        total_area: f64,
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
                .total_area = self.total_area,
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

/// Validate a positive, finite probe radius at public calculation boundaries.
pub fn validateProbeRadius(comptime T: type, probe_radius: T) !void {
    if (!std.math.isFinite(probe_radius) or !(probe_radius > 0.0)) {
        return error.InvalidInput;
    }
}

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

test "AtomInput validateFiniteAndRange rejects invalid numeric inputs" {
    const allocator = std.testing.allocator;
    const x = [_]f64{0.0};
    const y = [_]f64{std.math.inf(f64)};
    const z = [_]f64{0.0};
    const r = [_]f64{1.5};
    const input = AtomInput{
        .x = &x,
        .y = &y,
        .z = &z,
        .r = @constCast(&r),
        .allocator = allocator,
    };

    try std.testing.expectError(error.InvalidInput, input.validateFiniteAndRange());
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

test "validateProbeRadius rejects non-positive or non-finite values" {
    try std.testing.expectError(error.InvalidInput, validateProbeRadius(f64, 0.0));
    try std.testing.expectError(error.InvalidInput, validateProbeRadius(f64, std.math.nan(f64)));
    try validateProbeRadius(f32, 1.4);
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

test "FixedString4 fromSlice and slice" {
    const s1 = FixedString4.fromSlice("CA");
    try std.testing.expectEqualStrings("CA", s1.slice());
    try std.testing.expectEqual(@as(u8, 2), s1.len);

    const s2 = FixedString4.fromSlice("ALA");
    try std.testing.expectEqualStrings("ALA", s2.slice());
    try std.testing.expectEqual(@as(u8, 3), s2.len);

    // Truncation for > 4 chars
    const s3 = FixedString4.fromSlice("WATER");
    try std.testing.expectEqualStrings("WATE", s3.slice());
    try std.testing.expectEqual(@as(u8, 4), s3.len);

    // Empty string
    const s4 = FixedString4.fromSlice("");
    try std.testing.expectEqualStrings("", s4.slice());
    try std.testing.expectEqual(@as(u8, 0), s4.len);
}

test "FixedString4 equality" {
    const s1 = FixedString4.fromSlice("CA");
    const s2 = FixedString4.fromSlice("CA");
    const s3 = FixedString4.fromSlice("CB");

    try std.testing.expect(s1.eql(&s2));
    try std.testing.expect(!s1.eql(&s3));
    try std.testing.expect(s1.eqlSlice("CA"));
    try std.testing.expect(!s1.eqlSlice("CB"));
}

test "FixedString5 fromSlice and slice" {
    const s1 = FixedString5.fromSlice("ALA");
    try std.testing.expectEqualStrings("ALA", s1.slice());
    try std.testing.expectEqual(@as(u8, 3), s1.len);

    // 5-character comp_id (mmCIF extended CCD)
    const s2 = FixedString5.fromSlice("A1234");
    try std.testing.expectEqualStrings("A1234", s2.slice());
    try std.testing.expectEqual(@as(u8, 5), s2.len);

    // Truncation for > 5 chars
    const s3 = FixedString5.fromSlice("ABCDEF");
    try std.testing.expectEqualStrings("ABCDE", s3.slice());
    try std.testing.expectEqual(@as(u8, 5), s3.len);

    // Empty string
    const s4 = FixedString5.fromSlice("");
    try std.testing.expectEqualStrings("", s4.slice());
    try std.testing.expectEqual(@as(u8, 0), s4.len);
}

test "FixedString5 equality" {
    const s1 = FixedString5.fromSlice("ALA");
    const s2 = FixedString5.fromSlice("ALA");
    const s3 = FixedString5.fromSlice("GLY");

    try std.testing.expect(s1.eql(&s2));
    try std.testing.expect(!s1.eql(&s3));
    try std.testing.expect(s1.eqlSlice("ALA"));
    try std.testing.expect(!s1.eqlSlice("GLY"));
}

test "isWaterResidue" {
    for ([_][]const u8{ "HOH", "DOD", "WAT", "H2O", "D2O", "SOL", "TIP3", "SPCE", "hoh", " HOH " }) |name| {
        try std.testing.expect(isWaterResidue(name));
    }
    for ([_][]const u8{ "ALA", "HEM", "SO4", "GOL", "", "HOHOHO" }) |name| {
        try std.testing.expect(!isWaterResidue(name));
    }
}

test "AtomInput.removeResidues keeps every field of the other atoms" {
    const allocator = std.testing.allocator;
    var input = AtomInput{
        .x = try allocator.dupe(f64, &.{ 1, 2, 3 }),
        .y = try allocator.dupe(f64, &.{ 4, 5, 6 }),
        .z = try allocator.dupe(f64, &.{ 7, 8, 9 }),
        .r = try allocator.dupe(f64, &.{ 1.5, 1.6, 1.7 }),
        .residue = try allocator.dupe(FixedString5, &.{ FixedString5.fromSlice("ALA"), FixedString5.fromSlice("HOH"), FixedString5.fromSlice("ZN") }),
        .atom_name = try allocator.dupe(FixedString4, &.{ FixedString4.fromSlice("CA"), FixedString4.fromSlice("O"), FixedString4.fromSlice("ZN") }),
        .element = try allocator.dupe(u8, &.{ 6, 8, 30 }),
        .chain_id = try allocator.dupe(FixedString4, &.{ FixedString4.fromSlice("A"), FixedString4.fromSlice("A"), FixedString4.fromSlice("B") }),
        .residue_num = try allocator.dupe(i32, &.{ 1, 101, 201 }),
        .insertion_code = try allocator.dupe(FixedString4, &.{ FixedString4.fromSlice(""), FixedString4.fromSlice(""), FixedString4.fromSlice("") }),
        .allocator = allocator,
    };
    const chains = try allocator.alloc([]const u8, 3);
    chains[0] = try allocator.dupe(u8, "A");
    chains[1] = try allocator.dupe(u8, "A");
    chains[2] = try allocator.dupe(u8, "BLONG");
    input.chain_id_full = chains;
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), try input.removeResidues(.{}));
    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectEqualSlices(f64, &.{ 1, 3 }, input.x);
    try std.testing.expectEqualSlices(f64, &.{ 1.5, 1.7 }, input.r);
    try std.testing.expectEqualStrings("ZN", input.residue.?[1].slice());
    try std.testing.expectEqualSlices(u8, &.{ 6, 30 }, input.element.?);
    try std.testing.expectEqualSlices(i32, &.{ 1, 201 }, input.residue_num.?);
    try std.testing.expectEqualStrings("BLONG", input.chain_id_full.?[1]);

    // Nothing left to remove, unless a name is listed
    try std.testing.expectEqual(@as(usize, 0), try input.removeResidues(.{}));
    try std.testing.expectEqual(@as(usize, 1), try input.removeResidues(.{ .water = false, .names_csv = "SO4, zn" }));
    try std.testing.expectEqualSlices(f64, &.{1}, input.x);
}

test "AtomInput.removeResidues refuses to remove every atom" {
    const allocator = std.testing.allocator;
    var input = AtomInput{
        .x = try allocator.dupe(f64, &.{1}),
        .y = try allocator.dupe(f64, &.{2}),
        .z = try allocator.dupe(f64, &.{3}),
        .r = try allocator.dupe(f64, &.{1.5}),
        .residue = try allocator.dupe(FixedString5, &.{FixedString5.fromSlice("HOH")}),
        .allocator = allocator,
    };
    defer input.deinit();
    try std.testing.expectError(error.NoAtomsFound, input.removeResidues(.{}));
    try std.testing.expectEqual(@as(usize, 0), try input.removeResidues(.{ .water = false }));
    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
}

test "ResidueExclusion matches water by default and listed names" {
    const default: ResidueExclusion = .{};
    try std.testing.expect(default.excludes("HOH"));
    try std.testing.expect(!default.excludes("SO4"));
    try std.testing.expect(!default.isEmpty());

    const keep_water: ResidueExclusion = .{ .water = false };
    try std.testing.expect(!keep_water.excludes("HOH"));
    try std.testing.expect(keep_water.isEmpty());

    const listed: ResidueExclusion = .{ .water = false, .names_csv = "SO4, gol,", .names = &.{"EDO"} };
    try std.testing.expect(listed.excludes("SO4"));
    try std.testing.expect(listed.excludes(" GOL"));
    try std.testing.expect(listed.excludes("edo"));
    try std.testing.expect(!listed.excludes("HOH"));
    try std.testing.expect(!listed.excludes("SO"));
    try std.testing.expect(!(ResidueExclusion{ .water = false, .names_csv = " , " }).excludes(""));
}

test "gatherPositions converts the coordinates and selects atoms" {
    const x = [_]f32{ 1.5, 2.5, 3.5 };
    const y = [_]f32{ -1.0, 0.0, 1.0 };
    const z = [_]f32{ 10.0, 20.0, 30.0 };

    var all: [3]Vec3 = undefined;
    try gatherPositions(f64, &all, &x, &y, &z, null);
    try std.testing.expectEqual(Vec3{ .x = 2.5, .y = 0.0, .z = 20.0 }, all[1]);

    var picked: [2]Vec3f32 = undefined;
    try gatherPositions(f32, &picked, &x, &y, &z, &[_]u32{ 2, 0 });
    try std.testing.expectEqual(Vec3f32{ .x = 3.5, .y = 1.0, .z = 30.0 }, picked[0]);
    try std.testing.expectEqual(Vec3f32{ .x = 1.5, .y = -1.0, .z = 10.0 }, picked[1]);
}

test "gatherPositions rejects mismatched lengths and out-of-range indices" {
    const x = [_]f64{ 1.0, 2.0, 3.0 };
    const short = [_]f64{ 1.0, 2.0 };
    var dst: [2]Vec3 = undefined;

    try std.testing.expectError(error.InvalidInput, gatherPositions(f64, &dst, &x, &x, &x, null));
    try std.testing.expectError(error.InvalidInput, gatherPositions(f64, &dst, &x, &short, &x, &[_]u32{ 0, 1 }));
    try std.testing.expectError(error.InvalidInput, gatherPositions(f64, &dst, &x, &x, &x, &[_]u32{ 0, 3 }));
    try std.testing.expectError(error.InvalidInput, gatherPositions(f64, &dst, &x, &x, &x, &[_]u32{0}));
}

test "AtomsGen.validateFiniteAndRange checks the values in the calculation precision" {
    // 1e39 is finite in f64 but does not fit in f32.
    const big = [_]f64{1e39};
    const zero = [_]f64{0.0};
    var pos: [1]Vec3f32 = undefined;
    try gatherPositions(f32, &pos, &big, &zero, &zero, null);
    try std.testing.expectError(error.InvalidInput, (AtomsGen(f32){ .positions = &pos, .radii = &[_]f32{1.5} }).validateFiniteAndRange());

    const ok = [_]Vec3f32{.{ .x = 0, .y = 0, .z = 0 }};
    try (AtomsGen(f32){ .positions = &ok, .radii = &[_]f32{0.0} }).validateFiniteAndRange();
    try std.testing.expectError(error.InvalidInput, (AtomsGen(f32){ .positions = &ok, .radii = &[_]f32{-1.0} }).validateFiniteAndRange());
    try std.testing.expectError(error.InvalidInput, (AtomsGen(f32){ .positions = &ok, .radii = &[_]f32{ 1.0, 1.0 } }).validateFiniteAndRange());
    try std.testing.expectError(error.InvalidInput, (AtomsGen(f32){ .positions = &.{.{ .x = std.math.nan(f32), .y = 0, .z = 0 }}, .radii = &[_]f32{1.0} }).validateFiniteAndRange());
}

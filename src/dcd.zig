// DCD trajectory file reader
// Binary trajectory format used by NAMD and CHARMM
//
// DCD uses Fortran-style record format:
// - Each data block is bracketed by int32 byte-count markers
// - Endianness auto-detected from first int (should be 84)
// - Coordinates stored as separate X, Y, Z float32 arrays
// - Coordinates in Angstroms (unlike XTC which uses nanometers)
//
const std = @import("std");
const Allocator = std.mem.Allocator;

pub const DcdError = error{
    FileNotFound,
    InvalidMagic,
    EndOfFile,
    ReadError,
    BadFormat,
    OutOfMemory,
    FixedAtomsNotSupported,
};

/// CHARMM extension flags
const DCD_IS_CHARMM: i32 = 0x01;
const DCD_HAS_4DIMS: i32 = 0x02;
const DCD_HAS_EXTRA_BLOCK: i32 = 0x04;

/// DCD frame data
pub const DcdFrame = struct {
    step: i32,
    time: f32,
    coords: []f32, // flat array of x,y,z coordinates (length = natoms * 3)
    unitcell: ?[6]f64, // a, gamma, b, beta, alpha, c (CHARMM convention)

    pub fn deinit(self: *DcdFrame, allocator: Allocator) void {
        allocator.free(self.coords);
    }
};

/// DCD file reader
pub const DcdReader = struct {
    file: std.fs.File,
    allocator: Allocator,
    natoms: i32,
    nsets: i32, // number of frames in file
    istart: i32, // starting timestep
    nsavc: i32, // timesteps between saves
    delta: f32, // timestep value
    charmm: i32, // CHARMM flags
    reverse_endian: bool,
    frame_count: usize = 0,

    const Self = @This();

    pub fn open(allocator: Allocator, path: []const u8) !Self {
        const file = std.fs.cwd().openFile(path, .{}) catch {
            return DcdError.FileNotFound;
        };

        var reader = Self{
            .file = file,
            .allocator = allocator,
            .natoms = 0,
            .nsets = 0,
            .istart = 0,
            .nsavc = 0,
            .delta = 0,
            .charmm = 0,
            .reverse_endian = false,
        };

        reader.readHeader() catch |err| {
            file.close();
            return err;
        };

        return reader;
    }

    pub fn close(self: *Self) void {
        self.file.close();
    }

    pub fn getNumAtoms(self: *const Self) i32 {
        return self.natoms;
    }

    /// Read the DCD file header
    fn readHeader(self: *Self) !void {
        // 1. Read first int32 — should be 84
        var first_int = self.readRawInt32() catch return DcdError.ReadError;
        if (first_int != 84) {
            // Check if byte-swapped
            first_int = @bitCast(@byteSwap(@as(u32, @bitCast(first_int))));
            if (first_int == 84) {
                self.reverse_endian = true;
            } else {
                return DcdError.InvalidMagic;
            }
        }

        // 2. Read 84-byte header block
        var hdrbuf: [84]u8 = undefined;
        const hdr_n = self.file.readAll(&hdrbuf) catch return DcdError.ReadError;
        if (hdr_n < 84) return DcdError.BadFormat;

        // Check "CORD" magic
        if (hdrbuf[0] != 'C' or hdrbuf[1] != 'O' or hdrbuf[2] != 'R' or hdrbuf[3] != 'D') {
            return DcdError.InvalidMagic;
        }

        // Parse header fields (all at 4-byte aligned offsets within hdrbuf)
        // CHARMM detection: last int in header (offset 80) nonzero means CHARMM
        const charmm_ver = self.readIntFromBuf(hdrbuf[80..84]);
        if (charmm_ver != 0) {
            self.charmm = DCD_IS_CHARMM;
            if (self.readIntFromBuf(hdrbuf[44..48]) != 0) {
                self.charmm |= DCD_HAS_EXTRA_BLOCK;
            }
            if (self.readIntFromBuf(hdrbuf[48..52]) == 1) {
                self.charmm |= DCD_HAS_4DIMS;
            }
        }

        self.nsets = self.readIntFromBuf(hdrbuf[4..8]);
        self.istart = self.readIntFromBuf(hdrbuf[8..12]);
        self.nsavc = self.readIntFromBuf(hdrbuf[12..16]);

        const namnf = self.readIntFromBuf(hdrbuf[36..40]);

        // DELTA: stored as float for CHARMM, double for X-PLOR
        if ((self.charmm & DCD_IS_CHARMM) != 0) {
            self.delta = self.readFloatFromBuf(hdrbuf[40..44]);
        } else {
            self.delta = @floatCast(self.readDoubleFromBuf(hdrbuf[40..48]));
        }

        // Read trailing marker of first block (should be 84)
        const trail1 = try self.readInt32();
        if (trail1 != 84) return DcdError.BadFormat;

        // 3. Read title block
        const title_block_size = try self.readInt32();
        if (@mod(title_block_size - 4, 80) != 0) return DcdError.BadFormat;

        const ntitle = try self.readInt32();
        const title_bytes: usize = @intCast(ntitle * 80);
        if (title_bytes > 0) {
            // Skip title strings
            self.file.seekBy(@intCast(title_bytes)) catch return DcdError.ReadError;
        }

        // Read trailing marker of title block
        _ = try self.readRawInt32();

        // 4. Read natoms block
        const natom_marker = try self.readInt32();
        if (natom_marker != 4) return DcdError.BadFormat;

        self.natoms = try self.readInt32();

        const natom_trail = try self.readInt32();
        if (natom_trail != 4) return DcdError.BadFormat;

        // 5. Fixed atoms check
        if (namnf != 0) {
            return DcdError.FixedAtomsNotSupported;
        }
    }

    /// Read next frame
    pub fn readFrame(self: *Self) !DcdFrame {
        const natoms: usize = @intCast(self.natoms);

        // Calculate step and time
        const step: i32 = self.istart + @as(i32, @intCast(self.frame_count)) * self.nsavc;
        const time: f32 = @as(f32, @floatFromInt(self.frame_count)) * self.delta;

        // 1. Read unitcell if CHARMM + HAS_EXTRA_BLOCK
        var unitcell: ?[6]f64 = null;
        if ((self.charmm & DCD_IS_CHARMM) != 0 and (self.charmm & DCD_HAS_EXTRA_BLOCK) != 0) {
            const uc_marker = self.readInt32() catch return DcdError.EndOfFile;
            if (uc_marker == 48) {
                var uc: [6]f64 = undefined;
                for (0..6) |i| {
                    uc[i] = try self.readFloat64();
                }
                unitcell = uc;
            } else {
                // Unknown block, skip
                self.file.seekBy(@intCast(uc_marker)) catch return DcdError.ReadError;
            }
            // Read trailing marker
            _ = try self.readRawInt32();
        }

        // Allocate coordinates (interleaved x,y,z)
        const coords = self.allocator.alloc(f32, natoms * 3) catch return DcdError.OutOfMemory;
        errdefer self.allocator.free(coords);

        // 2. Read X coordinates: marker, float32[N], marker
        try self.readCoordBlock(coords, natoms, 0);

        // 3. Read Y coordinates
        try self.readCoordBlock(coords, natoms, 1);

        // 4. Read Z coordinates
        try self.readCoordBlock(coords, natoms, 2);

        // 5. Skip 4th dimension if present
        if ((self.charmm & DCD_IS_CHARMM) != 0 and (self.charmm & DCD_HAS_4DIMS) != 0) {
            const dim4_marker = try self.readInt32();
            self.file.seekBy(@intCast(dim4_marker)) catch return DcdError.ReadError;
            _ = try self.readRawInt32();
        }

        self.frame_count += 1;

        return DcdFrame{
            .step = step,
            .time = time,
            .coords = coords,
            .unitcell = unitcell,
        };
    }

    /// Read a coordinate block (X, Y, or Z) with Fortran record markers
    /// Coordinates are stored as separate arrays in DCD; we interleave them
    /// into the output coords array as [x0,y0,z0, x1,y1,z1, ...]
    fn readCoordBlock(self: *Self, coords: []f32, natoms: usize, component: usize) !void {
        const expected_size: i32 = @intCast(natoms * 4);
        const lead_marker = try self.readInt32();
        if (lead_marker != expected_size) return DcdError.BadFormat;

        if (self.reverse_endian) {
            // Read and byte-swap each float individually
            var raw: [4]u8 = undefined;
            for (0..natoms) |i| {
                const n = self.file.read(&raw) catch return DcdError.ReadError;
                if (n < 4) return DcdError.EndOfFile;
                const swapped = @byteSwap(std.mem.readInt(u32, &raw, .little));
                coords[i * 3 + component] = @bitCast(swapped);
            }
        } else {
            // Read all floats at once using a temporary buffer
            const byte_count = natoms * 4;
            const buf_bytes: [*]u8 = @ptrCast(coords.ptr);
            // We can't read directly into interleaved positions, use temp buffer
            const tmp = self.allocator.alloc(f32, natoms) catch return DcdError.OutOfMemory;
            defer self.allocator.free(tmp);

            const tmp_bytes: [*]u8 = @ptrCast(tmp.ptr);
            const n = self.file.readAll(tmp_bytes[0..byte_count]) catch return DcdError.ReadError;
            _ = buf_bytes;
            if (n < byte_count) return DcdError.EndOfFile;

            for (0..natoms) |i| {
                coords[i * 3 + component] = tmp[i];
            }
        }

        // Read trailing marker
        const trail_marker = try self.readInt32();
        if (trail_marker != expected_size) return DcdError.BadFormat;
    }

    // ============================================
    // Low-level I/O helpers
    // ============================================

    /// Read a raw int32 without endian swapping
    fn readRawInt32(self: *Self) !i32 {
        var buf: [4]u8 = undefined;
        const n = self.file.read(&buf) catch return DcdError.ReadError;
        if (n < 4) return DcdError.EndOfFile;
        return @bitCast(std.mem.readInt(u32, &buf, .little));
    }

    /// Read int32 with endian handling
    fn readInt32(self: *Self) !i32 {
        var buf: [4]u8 = undefined;
        const n = self.file.read(&buf) catch return DcdError.ReadError;
        if (n < 4) return DcdError.EndOfFile;
        if (self.reverse_endian) {
            return @bitCast(std.mem.readInt(u32, &buf, .big));
        }
        return @bitCast(std.mem.readInt(u32, &buf, .little));
    }

    /// Read float64 with endian handling
    fn readFloat64(self: *Self) !f64 {
        var buf: [8]u8 = undefined;
        const n = self.file.read(&buf) catch return DcdError.ReadError;
        if (n < 8) return DcdError.EndOfFile;
        if (self.reverse_endian) {
            return @bitCast(std.mem.readInt(u64, &buf, .big));
        }
        return @bitCast(std.mem.readInt(u64, &buf, .little));
    }

    /// Read int32 from a byte buffer with endian handling
    fn readIntFromBuf(self: *const Self, buf: *const [4]u8) i32 {
        if (self.reverse_endian) {
            return @bitCast(std.mem.readInt(u32, buf, .big));
        }
        return @bitCast(std.mem.readInt(u32, buf, .little));
    }

    /// Read float32 from a byte buffer with endian handling
    fn readFloatFromBuf(self: *const Self, buf: *const [4]u8) f32 {
        if (self.reverse_endian) {
            return @bitCast(std.mem.readInt(u32, buf, .big));
        }
        return @bitCast(std.mem.readInt(u32, buf, .little));
    }

    /// Read float64 from a byte buffer with endian handling
    fn readDoubleFromBuf(self: *const Self, buf: *const [8]u8) f64 {
        if (self.reverse_endian) {
            return @bitCast(std.mem.readInt(u64, buf, .big));
        }
        return @bitCast(std.mem.readInt(u64, buf, .little));
    }
};

// ============================================
// Tests
// ============================================

test "read 1l2y.dcd header" {
    const allocator = std.testing.allocator;

    var reader = DcdReader.open(allocator, "test_data/1l2y.dcd") catch |err| {
        if (err == DcdError.FileNotFound) return; // Skip if not available
        return err;
    };
    defer reader.close();

    // 1l2y has 304 atoms
    try std.testing.expectEqual(@as(i32, 304), reader.getNumAtoms());

    // Should have 38 frames (same as XTC)
    try std.testing.expectEqual(@as(i32, 38), reader.nsets);
}

test "read 1l2y.dcd first frame" {
    const allocator = std.testing.allocator;

    var reader = DcdReader.open(allocator, "test_data/1l2y.dcd") catch |err| {
        if (err == DcdError.FileNotFound) return;
        return err;
    };
    defer reader.close();

    const natoms = reader.getNumAtoms();
    try std.testing.expectEqual(@as(i32, 304), natoms);

    // Read first frame
    var frame = try reader.readFrame();
    defer frame.deinit(allocator);

    // Verify coordinate count
    try std.testing.expectEqual(@as(usize, 304 * 3), frame.coords.len);

    // DCD coordinates are in Angstroms
    // XTC coordinates for atom[0] are [-0.8901, 0.4127, -0.0555] in nm
    // In Angstroms: [-8.901, 4.127, -0.555]
    const tolerance: f32 = 0.05; // Slightly looser tolerance for format conversion
    try std.testing.expectApproxEqAbs(@as(f32, -8.901), frame.coords[0], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, 4.127), frame.coords[1], tolerance);
    try std.testing.expectApproxEqAbs(@as(f32, -0.555), frame.coords[2], tolerance);
}

test "read 1l2y.dcd all frames" {
    const allocator = std.testing.allocator;

    var reader = DcdReader.open(allocator, "test_data/1l2y.dcd") catch |err| {
        if (err == DcdError.FileNotFound) return;
        return err;
    };
    defer reader.close();

    var frame_count: usize = 0;
    while (true) {
        var frame = reader.readFrame() catch |err| {
            if (err == DcdError.EndOfFile) break;
            return err;
        };
        defer frame.deinit(allocator);
        frame_count += 1;

        // Verify frame has correct number of coordinates
        try std.testing.expectEqual(@as(usize, 304 * 3), frame.coords.len);
    }

    // 1l2y.dcd should have 38 frames (same as XTC)
    try std.testing.expectEqual(@as(usize, 38), frame_count);
}

//! Binary CCD dictionary format (ZSDC = Z-SASA Dictionary Compiled).
//!
//! A compact binary format for storing CCD component data needed for SASA
//! classification. This is a freesasa-zig original format, containing only
//! the fields required for ProtOr-compatible VdW radius assignment.
//!
//! ## Format Layout
//!
//! Header (12 bytes):
//!   [4B magic "ZSDC"] [1B version] [3B reserved] [4B component count LE]
//!
//! Component Record (variable length):
//!   [1B comp_id_len] [N bytes comp_id]
//!   [2B atom_count LE] [12B x atom_count packed atoms]
//!   [2B bond_count LE] [6B x bond_count packed bonds]

const std = @import("std");
const Allocator = std.mem.Allocator;
const ccd_parser = @import("ccd_parser.zig");
const hyb = @import("hybridization.zig");
const CompAtom = hyb.CompAtom;
const CompBond = hyb.CompBond;
const BondOrder = hyb.BondOrder;

pub const MAGIC: [4]u8 = .{ 'Z', 'S', 'D', 'C' };
pub const FORMAT_VERSION: u8 = 1;
const HEADER_SIZE: usize = 12;

pub const ReadError = error{
    InvalidMagic,
    UnsupportedVersion,
    UnexpectedEof,
    OutOfMemory,
    CountTooLarge,
};

// =============================================================================
// Packed binary representations
// =============================================================================

/// Packed atom record (12 bytes, extern struct for stable layout).
pub const PackedAtom = extern struct {
    atom_id: [4]u8,
    atom_id_len: u8,
    type_symbol: [4]u8,
    type_symbol_len: u8,
    flags: u8, // bit0 = leaving, bit1 = aromatic
    _pad: u8,

    comptime {
        std.debug.assert(@sizeOf(PackedAtom) == 12);
    }

    pub fn fromCompAtom(a: CompAtom) PackedAtom {
        var flags: u8 = 0;
        if (a.leaving) flags |= 0x01;
        if (a.aromatic) flags |= 0x02;
        return .{
            .atom_id = a.atom_id,
            .atom_id_len = @intCast(a.atom_id_len),
            .type_symbol = a.type_symbol,
            .type_symbol_len = @intCast(a.type_symbol_len),
            .flags = flags,
            ._pad = 0,
        };
    }

    pub fn toCompAtom(self: PackedAtom) CompAtom {
        const aid_len: u3 = @intCast(self.atom_id_len & 0x07);
        const ts_len: u3 = @intCast(self.type_symbol_len & 0x07);
        return .{
            .atom_id = self.atom_id,
            .atom_id_len = aid_len,
            .type_symbol = self.type_symbol,
            .type_symbol_len = ts_len,
            .leaving = (self.flags & 0x01) != 0,
            .aromatic = (self.flags & 0x02) != 0,
        };
    }
};

/// Packed bond record (6 bytes, extern struct for stable layout).
pub const PackedBond = extern struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: u8,
    flags: u8, // bit0 = aromatic

    comptime {
        std.debug.assert(@sizeOf(PackedBond) == 6);
    }

    pub fn fromCompBond(b: CompBond) PackedBond {
        var flags: u8 = 0;
        if (b.aromatic) flags |= 0x01;
        return .{
            .atom_idx_1 = std.mem.nativeToLittle(u16, b.atom_idx_1),
            .atom_idx_2 = std.mem.nativeToLittle(u16, b.atom_idx_2),
            .order = @intFromEnum(b.order),
            .flags = flags,
        };
    }

    pub fn toCompBond(self: PackedBond) CompBond {
        return .{
            .atom_idx_1 = std.mem.littleToNative(u16, self.atom_idx_1),
            .atom_idx_2 = std.mem.littleToNative(u16, self.atom_idx_2),
            .order = @enumFromInt(self.order),
            .aromatic = (self.flags & 0x01) != 0,
        };
    }
};

// =============================================================================
// Public API
// =============================================================================

/// Check if data starts with ZSDC magic bytes.
pub fn isBinaryDict(data: []const u8) bool {
    if (data.len < 4) return false;
    return std.mem.eql(u8, data[0..4], &MAGIC);
}

fn componentIdLessThan(_: void, a: []const u8, b: []const u8) bool {
    return std.mem.lessThan(u8, a, b);
}

/// Write a ComponentDict to binary ZSDC format.
pub fn writeDict(writer: *std.Io.Writer, dict: *const ccd_parser.ComponentDict) !void {
    // Count components
    const comp_count_usize = dict.components.count();
    const comp_count: u32 = @intCast(comp_count_usize);

    // Write header
    try writer.writeAll(&MAGIC);
    try writer.writeByte(FORMAT_VERSION);
    try writer.writeAll(&[3]u8{ 0, 0, 0 }); // reserved
    try writer.writeAll(&std.mem.toBytes(std.mem.nativeToLittle(u32, comp_count)));

    // Hash-map iteration order is intentionally not stable, but compiled ZSDC
    // bytes are build artifacts and should be reproducible. Write components in
    // component-ID order instead of raw map iteration order.
    const comp_ids = try dict.allocator.alloc([]const u8, comp_count_usize);
    defer dict.allocator.free(comp_ids);

    var idx: usize = 0;
    var it = dict.components.iterator();
    while (it.next()) |entry| : (idx += 1) {
        comp_ids[idx] = entry.key_ptr.*;
    }
    std.mem.sort([]const u8, comp_ids, {}, componentIdLessThan);

    // Write each component
    for (comp_ids) |comp_id| {
        const stored = dict.components.getPtr(comp_id) orelse return error.InvalidDict;

        // Component ID
        const cid_len: u8 = @intCast(@min(comp_id.len, 255));
        try writer.writeByte(cid_len);
        try writer.writeAll(comp_id[0..cid_len]);

        // Atom count and atoms
        const atom_count: u16 = @intCast(stored.atoms.len);
        try writer.writeAll(&std.mem.toBytes(std.mem.nativeToLittle(u16, atom_count)));
        for (stored.atoms) |atom| {
            const pa = PackedAtom.fromCompAtom(atom);
            try writer.writeAll(std.mem.asBytes(&pa));
        }

        // Bond count and bonds
        const bond_count: u16 = @intCast(stored.bonds.len);
        try writer.writeAll(&std.mem.toBytes(std.mem.nativeToLittle(u16, bond_count)));
        for (stored.bonds) |bond| {
            const pb = PackedBond.fromCompBond(bond);
            try writer.writeAll(std.mem.asBytes(&pb));
        }
    }
}

/// Read a ComponentDict from binary ZSDC format.
pub fn readDict(allocator: Allocator, reader: *std.Io.Reader) ReadError!ccd_parser.ComponentDict {
    var dict = ccd_parser.ComponentDict.init(allocator);
    errdefer dict.deinit();

    // Read header
    var header: [HEADER_SIZE]u8 = undefined;
    reader.readSliceAll(&header) catch return error.UnexpectedEof;

    if (!std.mem.eql(u8, header[0..4], &MAGIC)) return error.InvalidMagic;
    if (header[4] != FORMAT_VERSION) return error.UnsupportedVersion;

    const comp_count = std.mem.littleToNative(u32, std.mem.bytesToValue(u32, header[8..12]));
    if (comp_count > 1_000_000) return error.CountTooLarge;

    // Read components
    for (0..comp_count) |_| {
        // Read comp_id
        const cid_len_byte = (reader.takeArray(1) catch return error.UnexpectedEof)[0];
        const cid_len: usize = cid_len_byte;
        if (cid_len == 0 or cid_len > 255) return error.UnexpectedEof;
        var cid_buf: [255]u8 = undefined;
        reader.readSliceAll(cid_buf[0..cid_len]) catch return error.UnexpectedEof;

        // Read atom count
        var atom_count_buf: [2]u8 = undefined;
        reader.readSliceAll(&atom_count_buf) catch return error.UnexpectedEof;
        const atom_count: u16 = std.mem.littleToNative(u16, std.mem.bytesToValue(u16, &atom_count_buf));

        // Read atoms
        const atoms = allocator.alloc(CompAtom, atom_count) catch return error.OutOfMemory;
        var atoms_owned = true;
        defer if (atoms_owned) allocator.free(atoms);
        for (0..atom_count) |i| {
            var packed_buf: [@sizeOf(PackedAtom)]u8 = undefined;
            reader.readSliceAll(&packed_buf) catch return error.UnexpectedEof;
            const pa: PackedAtom = @bitCast(packed_buf);
            atoms[i] = pa.toCompAtom();
        }

        // Read bond count
        var bond_count_buf: [2]u8 = undefined;
        reader.readSliceAll(&bond_count_buf) catch return error.UnexpectedEof;
        const bond_count: u16 = std.mem.littleToNative(u16, std.mem.bytesToValue(u16, &bond_count_buf));

        // Read bonds
        const bonds = allocator.alloc(CompBond, bond_count) catch return error.OutOfMemory;
        var bonds_owned = true;
        defer if (bonds_owned) allocator.free(bonds);
        for (0..bond_count) |i| {
            var packed_buf: [@sizeOf(PackedBond)]u8 = undefined;
            reader.readSliceAll(&packed_buf) catch return error.UnexpectedEof;
            const pb: PackedBond = @bitCast(packed_buf);
            bonds[i] = pb.toCompBond();
        }

        // Build comp_id key
        const key = allocator.dupe(u8, cid_buf[0..cid_len]) catch return error.OutOfMemory;
        var key_owned = true;
        defer if (key_owned) allocator.free(key);

        // Build stored component
        var comp_id_fixed: [5]u8 = .{ 0, 0, 0, 0, 0 };
        const fixed_len: usize = @min(cid_len, 5);
        @memcpy(comp_id_fixed[0..fixed_len], cid_buf[0..fixed_len]);

        const stored = ccd_parser.StoredComponent{
            .comp_id = comp_id_fixed,
            .comp_id_len = @intCast(fixed_len),
            .atoms = atoms,
            .bonds = bonds,
            .allocator = allocator,
        };

        // Transfer ownership to dict
        dict.components.put(allocator, key, stored) catch return error.OutOfMemory;
        dict.owned_keys.append(allocator, key) catch {
            // Key is already in the map, so we should not free it here;
            // dict.deinit() via errdefer will handle cleanup.
            return error.OutOfMemory;
        };
        // Ownership successfully transferred to dict
        atoms_owned = false;
        bonds_owned = false;
        key_owned = false;
    }

    return dict;
}

/// Auto-detect format and load a ComponentDict from raw data.
/// If data starts with ZSDC magic, decode as binary; otherwise parse as CIF text.
pub fn loadDict(allocator: Allocator, data: []const u8) !ccd_parser.ComponentDict {
    if (isBinaryDict(data)) {
        var reader = std.Io.Reader.fixed(data);
        return readDict(allocator, &reader);
    } else {
        return ccd_parser.parseCcdData(allocator, data, null);
    }
}

// =============================================================================
// Tests
// =============================================================================

test "isBinaryDict — valid magic" {
    const data = [_]u8{ 'Z', 'S', 'D', 'C', 1, 0, 0, 0, 0, 0, 0, 0 };
    try std.testing.expect(isBinaryDict(&data));
}

test "isBinaryDict — wrong magic" {
    const data = [_]u8{ 'Z', 'S', 'D', 'X', 1, 0, 0, 0, 0, 0, 0, 0 };
    try std.testing.expect(!isBinaryDict(&data));
}

test "isBinaryDict — empty" {
    const data = [_]u8{};
    try std.testing.expect(!isBinaryDict(&data));
}

test "isBinaryDict — too short" {
    const data = [_]u8{ 'Z', 'S', 'D' };
    try std.testing.expect(!isBinaryDict(&data));
}

test "round-trip: create ComponentDict -> writeDict -> readDict -> verify" {
    const allocator = std.testing.allocator;

    // Build a synthetic ComponentDict
    var dict = ccd_parser.ComponentDict.init(allocator);
    defer dict.deinit();

    // Create atoms for component "ALA"
    const atoms = try allocator.alloc(CompAtom, 2);
    atoms[0] = CompAtom.init("CA", "C");
    atoms[1] = CompAtom.init("N", "N");
    atoms[1].leaving = true;
    atoms[0].aromatic = true;

    // Create bonds
    const bonds = try allocator.alloc(CompBond, 1);
    bonds[0] = .{
        .atom_idx_1 = 0,
        .atom_idx_2 = 1,
        .order = .single,
        .aromatic = false,
    };

    const key = try allocator.dupe(u8, "ALA");
    const stored = ccd_parser.StoredComponent{
        .comp_id = .{ 'A', 'L', 'A', 0, 0 },
        .comp_id_len = 3,
        .atoms = atoms,
        .bonds = bonds,
        .allocator = allocator,
    };
    try dict.components.put(allocator, key, stored);
    try dict.owned_keys.append(allocator, key);

    // Write to buffer
    var buf: [4096]u8 = undefined;
    var w = std.Io.Writer.fixed(&buf);
    try writeDict(&w, &dict);

    // Read back
    const written = buf[0..w.end];
    var read_reader = std.Io.Reader.fixed(written);
    var read_dict = try readDict(allocator, &read_reader);
    defer read_dict.deinit();

    // Verify
    try std.testing.expectEqual(@as(usize, 1), read_dict.components.count());

    const comp = read_dict.get("ALA") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 2), comp.atoms.len);
    try std.testing.expectEqual(@as(usize, 1), comp.bonds.len);

    // Verify atom fields
    try std.testing.expectEqualSlices(u8, "CA", comp.atoms[0].atomIdSlice());
    try std.testing.expectEqualSlices(u8, "C", comp.atoms[0].typeSymbolSlice());
    try std.testing.expect(comp.atoms[0].aromatic);
    try std.testing.expect(!comp.atoms[0].leaving);

    try std.testing.expectEqualSlices(u8, "N", comp.atoms[1].atomIdSlice());
    try std.testing.expectEqualSlices(u8, "N", comp.atoms[1].typeSymbolSlice());
    try std.testing.expect(!comp.atoms[1].aromatic);
    try std.testing.expect(comp.atoms[1].leaving);

    // Verify bond fields
    try std.testing.expectEqual(@as(u16, 0), comp.bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), comp.bonds[0].atom_idx_2);
    try std.testing.expectEqual(BondOrder.single, comp.bonds[0].order);
    try std.testing.expect(!comp.bonds[0].aromatic);
}

test "round-trip: empty dictionary" {
    const allocator = std.testing.allocator;

    var dict = ccd_parser.ComponentDict.init(allocator);
    defer dict.deinit();

    var buf: [256]u8 = undefined;
    var w = std.Io.Writer.fixed(&buf);
    try writeDict(&w, &dict);

    var read_reader = std.Io.Reader.fixed(buf[0..w.end]);
    var read_dict = try readDict(allocator, &read_reader);
    defer read_dict.deinit();

    try std.testing.expectEqual(@as(usize, 0), read_dict.components.count());
}

fn addTestComponent(allocator: Allocator, dict: *ccd_parser.ComponentDict, comp_id: []const u8) !void {
    const atoms = try allocator.alloc(CompAtom, 1);
    atoms[0] = CompAtom.init("C1", "C");

    const bonds = try allocator.alloc(CompBond, 0);
    const key = try allocator.dupe(u8, comp_id);
    errdefer allocator.free(key);

    var stored_id = [_]u8{ 0, 0, 0, 0, 0 };
    const id_len = @min(comp_id.len, stored_id.len);
    @memcpy(stored_id[0..id_len], comp_id[0..id_len]);

    const stored = ccd_parser.StoredComponent{
        .comp_id = stored_id,
        .comp_id_len = @intCast(id_len),
        .atoms = atoms,
        .bonds = bonds,
        .allocator = allocator,
    };
    try dict.components.put(allocator, key, stored);
    try dict.owned_keys.append(allocator, key);
}

fn writeTestDictAlloc(allocator: Allocator, dict: *const ccd_parser.ComponentDict) ![]u8 {
    var out: std.Io.Writer.Allocating = .init(allocator);
    errdefer out.deinit();
    try writeDict(&out.writer, dict);
    return out.toOwnedSlice();
}

test "writeDict emits deterministic bytes independent of insertion order" {
    const allocator = std.testing.allocator;
    const ids = [_][]const u8{ "ZZZ", "AAA", "MSE", "HOH", "ATP", "GLY" };

    var forward = ccd_parser.ComponentDict.init(allocator);
    defer forward.deinit();
    for (ids) |id| try addTestComponent(allocator, &forward, id);

    var reverse = ccd_parser.ComponentDict.init(allocator);
    defer reverse.deinit();
    var i = ids.len;
    while (i > 0) {
        i -= 1;
        try addTestComponent(allocator, &reverse, ids[i]);
    }

    const forward_bytes = try writeTestDictAlloc(allocator, &forward);
    defer allocator.free(forward_bytes);
    const reverse_bytes = try writeTestDictAlloc(allocator, &reverse);
    defer allocator.free(reverse_bytes);

    try std.testing.expectEqualSlices(u8, forward_bytes, reverse_bytes);
}

test "writeDict emits components sorted by component ID" {
    const allocator = std.testing.allocator;
    const ids = [_][]const u8{ "ZZZ", "AAA", "MSE", "HOH", "ATP", "GLY" };

    var dict = ccd_parser.ComponentDict.init(allocator);
    defer dict.deinit();
    for (ids) |id| try addTestComponent(allocator, &dict, id);

    const bytes = try writeTestDictAlloc(allocator, &dict);
    defer allocator.free(bytes);

    var reader = std.Io.Reader.fixed(bytes);
    var loaded = try readDict(allocator, &reader);
    defer loaded.deinit();

    const sorted_ids = [_][]const u8{ "AAA", "ATP", "GLY", "HOH", "MSE", "ZZZ" };
    var offset: usize = HEADER_SIZE;
    for (sorted_ids) |expected_id| {
        const actual_len = bytes[offset];
        offset += 1;
        try std.testing.expectEqual(expected_id.len, actual_len);
        try std.testing.expectEqualStrings(expected_id, bytes[offset .. offset + actual_len]);
        offset += actual_len;
        const atom_count = std.mem.littleToNative(u16, std.mem.bytesToValue(u16, bytes[offset .. offset + 2]));
        offset += 2 + (@as(usize, atom_count) * @sizeOf(PackedAtom));
        const bond_count = std.mem.littleToNative(u16, std.mem.bytesToValue(u16, bytes[offset .. offset + 2]));
        offset += 2 + (@as(usize, bond_count) * @sizeOf(PackedBond));
    }
    try std.testing.expectEqual(bytes.len, offset);
    try std.testing.expectEqual(@as(usize, sorted_ids.len), loaded.components.count());
}

test "writeDict synthetic fixture hash is stable" {
    const allocator = std.testing.allocator;
    const ids = [_][]const u8{ "ZZZ", "AAA", "MSE", "HOH", "ATP", "GLY" };

    var dict = ccd_parser.ComponentDict.init(allocator);
    defer dict.deinit();
    for (ids) |id| try addTestComponent(allocator, &dict, id);

    const bytes = try writeTestDictAlloc(allocator, &dict);
    defer allocator.free(bytes);

    const Sha256 = std.crypto.hash.sha2.Sha256;
    var digest: [Sha256.digest_length]u8 = undefined;
    Sha256.hash(bytes, &digest, .{});
    const hex = std.fmt.bytesToHex(digest, .lower);
    try std.testing.expectEqualStrings("3f759a9f4de717c8d3232694cfaca01bddf48a2ca041591064c6b035e2e10913", &hex);
}

test "loadDict — auto-detect binary" {
    const allocator = std.testing.allocator;

    // Build a minimal binary dict in memory
    var dict = ccd_parser.ComponentDict.init(allocator);
    defer dict.deinit();

    var buf: [256]u8 = undefined;
    var w = std.Io.Writer.fixed(&buf);
    try writeDict(&w, &dict);

    const data = buf[0..w.end];
    var loaded = try loadDict(allocator, data);
    defer loaded.deinit();

    try std.testing.expectEqual(@as(usize, 0), loaded.components.count());
}

test "loadDict — auto-detect CIF text" {
    const allocator = std.testing.allocator;

    // Minimal CIF data with no actual loops (just enough to parse)
    const cif_text = "data_test\n";
    var loaded = try loadDict(allocator, cif_text);
    defer loaded.deinit();

    try std.testing.expectEqual(@as(usize, 0), loaded.components.count());
}

test "PackedAtom size is 12 bytes" {
    try std.testing.expectEqual(@as(usize, 12), @sizeOf(PackedAtom));
}

test "PackedBond size is 6 bytes" {
    try std.testing.expectEqual(@as(usize, 6), @sizeOf(PackedBond));
}

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
//!
//! ## Validation
//!
//! The reader treats every count, length and enum byte in the file as
//! untrusted: lengths must fit the fixed arrays they index (atom IDs and type
//! symbols hold at most 4 bytes), bond orders must name a `BondOrder`, bond
//! atom indices must refer to atoms of the same component, counts must fit in
//! the bytes that remain, and component IDs must be unique. The writer refuses
//! to emit anything the reader would reject instead of truncating a count.

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

/// Maximum number of components accepted by the reader.
const MAX_COMPONENTS: u32 = 1_000_000;
/// Smallest possible component record: 1B id length + 1B id + 2B atom count + 2B bond count.
const MIN_COMPONENT_RECORD_SIZE: usize = 6;
/// Longest atom ID / type symbol a packed atom can hold.
const MAX_NAME_LEN: u8 = 4;

pub const ReadError = error{
    InvalidMagic,
    UnsupportedVersion,
    /// The data ends before a declared count, length or record is complete.
    UnexpectedEof,
    OutOfMemory,
    /// The component count in the header exceeds the supported maximum.
    CountTooLarge,
    /// A component ID has length zero.
    InvalidComponentId,
    /// An atom ID or type symbol length exceeds the 4 bytes it indexes.
    InvalidAtomLength,
    /// A bond order byte does not name a `BondOrder`.
    InvalidBondOrder,
    /// A bond refers to an atom index outside its component.
    InvalidBondIndex,
    /// Two component records share the same ID.
    DuplicateComponent,
    /// Bytes remain after the last declared component.
    TrailingData,
};

pub const WriteError = error{
    /// The dictionary holds more components than the 32-bit header count can store.
    TooManyComponents,
    /// A component has more atoms than the 16-bit atom count can store.
    TooManyAtoms,
    /// A component has more bonds than the 16-bit bond count can store.
    TooManyBonds,
    /// A component ID is empty or longer than the 255 bytes its length byte can store.
    InvalidComponentId,
    /// A bond refers to an atom index outside its component.
    InvalidBondIndex,
    OutOfMemory,
    WriteFailed,
};

/// Names the component that made `writeDictDiag` fail (borrowed from the dictionary).
pub const WriteDiagnostic = struct {
    comp_id: []const u8 = "",
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

    /// Convert to a `CompAtom`. Fails with `InvalidAtomLength` when a length
    /// byte exceeds the 4 bytes of the array it indexes.
    pub fn toCompAtom(self: PackedAtom) error{InvalidAtomLength}!CompAtom {
        if (self.atom_id_len > MAX_NAME_LEN or self.type_symbol_len > MAX_NAME_LEN) {
            return error.InvalidAtomLength;
        }
        const aid_len: u3 = @intCast(self.atom_id_len);
        const ts_len: u3 = @intCast(self.type_symbol_len);
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

    /// Convert to a `CompBond`. Fails with `InvalidBondOrder` when the order
    /// byte does not name a `BondOrder`.
    pub fn toCompBond(self: PackedBond) error{InvalidBondOrder}!CompBond {
        const order = std.enums.fromInt(BondOrder, self.order) orelse return error.InvalidBondOrder;
        return .{
            .atom_idx_1 = std.mem.littleToNative(u16, self.atom_idx_1),
            .atom_idx_2 = std.mem.littleToNative(u16, self.atom_idx_2),
            .order = order,
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
///
/// Fails instead of wrapping a count: see `WriteError`. The whole dictionary
/// is validated before the first byte is written.
pub fn writeDict(writer: *std.Io.Writer, dict: *const ccd_parser.ComponentDict) WriteError!void {
    return writeDictDiag(writer, dict, null);
}

/// Like `writeDict`, but on a per-component error stores the offending
/// component ID in `diag` so callers can report it.
pub fn writeDictDiag(
    writer: *std.Io.Writer,
    dict: *const ccd_parser.ComponentDict,
    diag: ?*WriteDiagnostic,
) WriteError!void {
    const comp_count_usize = dict.components.count();
    const comp_count = std.math.cast(u32, comp_count_usize) orelse return error.TooManyComponents;

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

    // Validate every component before the first byte is written.
    for (comp_ids) |comp_id| {
        const stored = dict.components.getPtr(comp_id) orelse unreachable;
        checkComponentWritable(comp_id, stored) catch |err| {
            if (diag) |d| d.comp_id = comp_id;
            return err;
        };
    }

    try writer.writeAll(&MAGIC);
    try writer.writeByte(FORMAT_VERSION);
    try writer.writeAll(&[3]u8{ 0, 0, 0 }); // reserved
    try writer.writeAll(&std.mem.toBytes(std.mem.nativeToLittle(u32, comp_count)));

    for (comp_ids) |comp_id| {
        const stored = dict.components.getPtr(comp_id) orelse unreachable;

        // Lengths below were checked by checkComponentWritable.
        try writer.writeByte(@intCast(comp_id.len));
        try writer.writeAll(comp_id);

        const atom_count: u16 = @intCast(stored.atoms.len);
        try writer.writeAll(&std.mem.toBytes(std.mem.nativeToLittle(u16, atom_count)));
        for (stored.atoms) |atom| {
            const pa = PackedAtom.fromCompAtom(atom);
            try writer.writeAll(std.mem.asBytes(&pa));
        }

        const bond_count: u16 = @intCast(stored.bonds.len);
        try writer.writeAll(&std.mem.toBytes(std.mem.nativeToLittle(u16, bond_count)));
        for (stored.bonds) |bond| {
            const pb = PackedBond.fromCompBond(bond);
            try writer.writeAll(std.mem.asBytes(&pb));
        }
    }
}

/// Check that a component fits the on-disk fields and that the reader accepts it.
fn checkComponentWritable(comp_id: []const u8, stored: *const ccd_parser.StoredComponent) WriteError!void {
    if (comp_id.len == 0 or comp_id.len > std.math.maxInt(u8)) return error.InvalidComponentId;
    if (stored.atoms.len > std.math.maxInt(u16)) return error.TooManyAtoms;
    if (stored.bonds.len > std.math.maxInt(u16)) return error.TooManyBonds;
    for (stored.bonds) |bond| {
        if (bond.atom_idx_1 >= stored.atoms.len or bond.atom_idx_2 >= stored.atoms.len) {
            return error.InvalidBondIndex;
        }
    }
}

/// Bounds-checked cursor over the bytes of a ZSDC file.
const Cursor = struct {
    data: []const u8,
    pos: usize = 0,

    fn remaining(self: *const Cursor) usize {
        return self.data.len - self.pos;
    }

    /// Take `n` bytes, or fail without consuming anything if fewer remain.
    fn take(self: *Cursor, n: usize) ReadError![]const u8 {
        if (n > self.remaining()) return error.UnexpectedEof;
        const out = self.data[self.pos..][0..n];
        self.pos += n;
        return out;
    }

    fn takeU16(self: *Cursor) ReadError!u16 {
        const b = try self.take(2);
        return std.mem.readInt(u16, b[0..2], .little);
    }
};

/// Read a ComponentDict from binary ZSDC format.
///
/// The reader is drained into memory first so that every count in the file can
/// be checked against the number of bytes that actually remain.
pub fn readDict(allocator: Allocator, reader: *std.Io.Reader) ReadError!ccd_parser.ComponentDict {
    const data = reader.allocRemaining(allocator, .unlimited) catch |err| switch (err) {
        error.OutOfMemory => return error.OutOfMemory,
        else => return error.UnexpectedEof,
    };
    defer allocator.free(data);
    return parseDict(allocator, data);
}

/// Parse a ZSDC image held in memory. Every length, count and enum byte is
/// validated before use; see the module documentation.
pub fn parseDict(allocator: Allocator, data: []const u8) ReadError!ccd_parser.ComponentDict {
    var dict = ccd_parser.ComponentDict.init(allocator);
    errdefer dict.deinit();

    var cur = Cursor{ .data = data };

    const header = try cur.take(HEADER_SIZE);
    if (!std.mem.eql(u8, header[0..4], &MAGIC)) return error.InvalidMagic;
    if (header[4] != FORMAT_VERSION) return error.UnsupportedVersion;

    const comp_count = std.mem.readInt(u32, header[8..12], .little);
    if (comp_count > MAX_COMPONENTS) return error.CountTooLarge;
    // Each record takes at least MIN_COMPONENT_RECORD_SIZE bytes, so a count the
    // file cannot hold is rejected before any per-component work.
    if (@as(usize, comp_count) * MIN_COMPONENT_RECORD_SIZE > cur.remaining()) return error.UnexpectedEof;

    for (0..comp_count) |_| {
        const cid_len = (try cur.take(1))[0];
        if (cid_len == 0) return error.InvalidComponentId;
        const cid = try cur.take(cid_len);

        // Take the whole atom and bond blocks before allocating so an oversized
        // count fails on the remaining size, not on the allocator.
        const atom_count = try cur.takeU16();
        const atom_bytes = try cur.take(@as(usize, atom_count) * @sizeOf(PackedAtom));

        const bond_count = try cur.takeU16();
        const bond_bytes = try cur.take(@as(usize, bond_count) * @sizeOf(PackedBond));

        const key = try allocator.dupe(u8, cid);
        errdefer allocator.free(key);
        const atoms = try allocator.alloc(CompAtom, atom_count);
        errdefer allocator.free(atoms);
        const bonds = try allocator.alloc(CompBond, bond_count);
        errdefer allocator.free(bonds);

        for (atoms, 0..) |*atom, i| {
            const raw = atom_bytes[i * @sizeOf(PackedAtom) ..][0..@sizeOf(PackedAtom)];
            const pa = std.mem.bytesToValue(PackedAtom, raw);
            atom.* = try pa.toCompAtom();
        }
        for (bonds, 0..) |*bond, i| {
            const raw = bond_bytes[i * @sizeOf(PackedBond) ..][0..@sizeOf(PackedBond)];
            const pb = std.mem.bytesToValue(PackedBond, raw);
            bond.* = try pb.toCompBond();
            if (bond.atom_idx_1 >= atom_count or bond.atom_idx_2 >= atom_count) return error.InvalidBondIndex;
        }

        var comp_id_fixed: [5]u8 = .{ 0, 0, 0, 0, 0 };
        const fixed_len: usize = @min(cid.len, comp_id_fixed.len);
        @memcpy(comp_id_fixed[0..fixed_len], cid[0..fixed_len]);

        try dict.owned_keys.ensureUnusedCapacity(allocator, 1);
        const gop = try dict.components.getOrPut(allocator, key);
        if (gop.found_existing) return error.DuplicateComponent;
        gop.value_ptr.* = .{
            .comp_id = comp_id_fixed,
            .comp_id_len = @intCast(fixed_len),
            .atoms = atoms,
            .bonds = bonds,
            .allocator = allocator,
        };
        // Ownership of key, atoms and bonds now belongs to the dictionary.
        dict.owned_keys.appendAssumeCapacity(key);
    }

    if (cur.remaining() != 0) return error.TrailingData;

    return dict;
}

/// Auto-detect format and load a ComponentDict from raw data.
/// If data starts with ZSDC magic, decode as binary; otherwise parse as CIF text.
pub fn loadDict(allocator: Allocator, data: []const u8) !ccd_parser.ComponentDict {
    if (isBinaryDict(data)) {
        return parseDict(allocator, data);
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

// =============================================================================
// Malformed input and limit tests
// =============================================================================

/// Builds raw ZSDC images for tests, including ones the writer would refuse.
const TestImage = struct {
    buf: std.ArrayList(u8) = .empty,

    fn deinit(self: *TestImage, allocator: Allocator) void {
        self.buf.deinit(allocator);
    }

    fn header(self: *TestImage, allocator: Allocator, comp_count: u32) !void {
        try self.buf.appendSlice(allocator, &MAGIC);
        try self.buf.appendSlice(allocator, &[_]u8{ FORMAT_VERSION, 0, 0, 0 });
        try self.buf.appendSlice(allocator, &std.mem.toBytes(std.mem.nativeToLittle(u32, comp_count)));
    }

    fn u16le(self: *TestImage, allocator: Allocator, v: u16) !void {
        try self.buf.appendSlice(allocator, &std.mem.toBytes(std.mem.nativeToLittle(u16, v)));
    }

    /// Append a component record whose declared counts match the slices given.
    fn component(
        self: *TestImage,
        allocator: Allocator,
        id: []const u8,
        atoms: []const PackedAtom,
        bonds: []const PackedBond,
    ) !void {
        try self.buf.append(allocator, @intCast(id.len));
        try self.buf.appendSlice(allocator, id);
        try self.u16le(allocator, @intCast(atoms.len));
        for (atoms) |a| try self.buf.appendSlice(allocator, std.mem.asBytes(&a));
        try self.u16le(allocator, @intCast(bonds.len));
        for (bonds) |b| try self.buf.appendSlice(allocator, std.mem.asBytes(&b));
    }
};

fn testAtom(id: []const u8, symbol: []const u8) PackedAtom {
    return PackedAtom.fromCompAtom(CompAtom.init(id, symbol));
}

fn testBond(i: u16, j: u16, order: u8) PackedBond {
    return .{
        .atom_idx_1 = std.mem.nativeToLittle(u16, i),
        .atom_idx_2 = std.mem.nativeToLittle(u16, j),
        .order = order,
        .flags = 0,
    };
}

/// A valid two-component image used as the base for corruption tests.
fn validTestImage(allocator: Allocator) !TestImage {
    var img = TestImage{};
    errdefer img.deinit(allocator);
    try img.header(allocator, 2);
    try img.component(
        allocator,
        "ALA",
        &.{ testAtom("N", "N"), testAtom("CA", "C"), testAtom("CB", "C") },
        &.{ testBond(0, 1, 0), testBond(1, 2, 1) },
    );
    try img.component(allocator, "HOH", &.{testAtom("O", "O")}, &.{});
    return img;
}

fn expectParseError(expected: ReadError, allocator: Allocator, data: []const u8) !void {
    try std.testing.expectError(expected, parseDict(allocator, data));
}

test "parseDict accepts a valid image built by hand" {
    const allocator = std.testing.allocator;
    var img = try validTestImage(allocator);
    defer img.deinit(allocator);

    var dict = try parseDict(allocator, img.buf.items);
    defer dict.deinit();
    try std.testing.expectEqual(@as(usize, 2), dict.components.count());
    const ala = dict.get("ALA") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 3), ala.atoms.len);
    try std.testing.expectEqual(BondOrder.double, ala.bonds[1].order);
}

test "parseDict rejects atom lengths above the 4-byte arrays" {
    const allocator = std.testing.allocator;
    // 5, 6 and 7 used to pass `len & 0x07`; 8 wrapped to 0.
    for ([_]u8{ 5, 6, 7, 8, 255 }) |bad| {
        var atom = testAtom("CA", "C");
        atom.atom_id_len = bad;
        var img = TestImage{};
        defer img.deinit(allocator);
        try img.header(allocator, 1);
        try img.component(allocator, "ALA", &.{atom}, &.{});
        try expectParseError(error.InvalidAtomLength, allocator, img.buf.items);

        atom = testAtom("CA", "C");
        atom.type_symbol_len = bad;
        var img2 = TestImage{};
        defer img2.deinit(allocator);
        try img2.header(allocator, 1);
        try img2.component(allocator, "ALA", &.{atom}, &.{});
        try expectParseError(error.InvalidAtomLength, allocator, img2.buf.items);
    }
}

test "parseDict rejects bond order bytes that name no BondOrder" {
    const allocator = std.testing.allocator;
    const valid_max: u8 = @typeInfo(BondOrder).@"enum".fields.len - 1;
    for ([_]u8{ valid_max + 1, 7, 200, 255 }) |bad| {
        var img = TestImage{};
        defer img.deinit(allocator);
        try img.header(allocator, 1);
        try img.component(allocator, "ALA", &.{ testAtom("N", "N"), testAtom("CA", "C") }, &.{testBond(0, 1, bad)});
        try expectParseError(error.InvalidBondOrder, allocator, img.buf.items);
    }
    // The largest real order still loads.
    var ok = TestImage{};
    defer ok.deinit(allocator);
    try ok.header(allocator, 1);
    try ok.component(allocator, "ALA", &.{ testAtom("N", "N"), testAtom("CA", "C") }, &.{testBond(0, 1, valid_max)});
    var dict = try parseDict(allocator, ok.buf.items);
    dict.deinit();
}

test "parseDict rejects bond atom indices outside the component" {
    const allocator = std.testing.allocator;
    var img = TestImage{};
    defer img.deinit(allocator);
    try img.header(allocator, 1);
    try img.component(allocator, "ALA", &.{ testAtom("N", "N"), testAtom("CA", "C") }, &.{testBond(0, 2, 0)});
    try expectParseError(error.InvalidBondIndex, allocator, img.buf.items);
}

test "parseDict rejects counts larger than the remaining bytes" {
    const allocator = std.testing.allocator;

    // Component count that cannot fit in the file.
    var many = TestImage{};
    defer many.deinit(allocator);
    try many.header(allocator, 900_000);
    try expectParseError(error.UnexpectedEof, allocator, many.buf.items);

    // Component count over the hard limit.
    var huge = TestImage{};
    defer huge.deinit(allocator);
    try huge.header(allocator, std.math.maxInt(u32));
    try expectParseError(error.CountTooLarge, allocator, huge.buf.items);

    // Atom count of 65535 followed by almost nothing.
    var atoms = TestImage{};
    defer atoms.deinit(allocator);
    try atoms.header(allocator, 1);
    try atoms.buf.appendSlice(allocator, &[_]u8{ 3, 'A', 'L', 'A' });
    try atoms.u16le(allocator, std.math.maxInt(u16));
    try atoms.buf.appendSlice(allocator, &[_]u8{ 0, 0, 0 });
    try expectParseError(error.UnexpectedEof, allocator, atoms.buf.items);

    // Bond count of 65535 followed by almost nothing.
    var bonds = TestImage{};
    defer bonds.deinit(allocator);
    try bonds.header(allocator, 1);
    try bonds.component(allocator, "ALA", &.{testAtom("N", "N")}, &.{});
    bonds.buf.items.len -= 2; // drop the zero bond count
    try bonds.u16le(allocator, std.math.maxInt(u16));
    try bonds.buf.appendSlice(allocator, &[_]u8{ 0, 0, 0 });
    try expectParseError(error.UnexpectedEof, allocator, bonds.buf.items);
}

test "parseDict rejects empty and duplicate component IDs and trailing bytes" {
    const allocator = std.testing.allocator;

    var empty_id = TestImage{};
    defer empty_id.deinit(allocator);
    try empty_id.header(allocator, 1);
    try empty_id.component(allocator, "", &.{testAtom("N", "N")}, &.{});
    try expectParseError(error.InvalidComponentId, allocator, empty_id.buf.items);

    var dup = TestImage{};
    defer dup.deinit(allocator);
    try dup.header(allocator, 2);
    try dup.component(allocator, "ALA", &.{testAtom("N", "N")}, &.{});
    try dup.component(allocator, "ALA", &.{testAtom("CA", "C")}, &.{});
    try expectParseError(error.DuplicateComponent, allocator, dup.buf.items);

    var img = try validTestImage(allocator);
    defer img.deinit(allocator);
    try img.buf.append(allocator, 0);
    try expectParseError(error.TrailingData, allocator, img.buf.items);
}

test "parseDict rejects every truncation of a valid image" {
    const allocator = std.testing.allocator;
    var img = try validTestImage(allocator);
    defer img.deinit(allocator);

    for (0..img.buf.items.len) |len| {
        try std.testing.expectError(error.UnexpectedEof, parseDict(allocator, img.buf.items[0..len]));
    }
}

test "parseDict survives any single corrupted byte without illegal behavior" {
    const allocator = std.testing.allocator;
    var img = try validTestImage(allocator);
    defer img.deinit(allocator);

    // Tests run with safety checks, so any out-of-range index or enum would panic.
    for (0..img.buf.items.len) |i| {
        for ([_]u8{ 0x00, 0x07, 0x80, 0xFF }) |v| {
            const saved = img.buf.items[i];
            img.buf.items[i] = v;
            defer img.buf.items[i] = saved;
            if (parseDict(allocator, img.buf.items)) |loaded| {
                var dict = loaded;
                dict.deinit();
            } else |_| {}
        }
    }
}

test "parseDict has no leaks on allocation failure" {
    var img = try validTestImage(std.testing.allocator);
    defer img.deinit(std.testing.allocator);

    const Check = struct {
        fn run(allocator: Allocator, data: []const u8) !void {
            var dict = try parseDict(allocator, data);
            dict.deinit();
        }
    };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, Check.run, .{img.buf.items});

    // Error paths free everything too: duplicate ID after a successful record.
    var dup = TestImage{};
    defer dup.deinit(std.testing.allocator);
    try dup.header(std.testing.allocator, 2);
    try dup.component(std.testing.allocator, "ALA", &.{testAtom("N", "N")}, &.{});
    try dup.component(std.testing.allocator, "ALA", &.{testAtom("N", "N")}, &.{});
    const CheckErr = struct {
        fn run(allocator: Allocator, data: []const u8) !void {
            var dict = parseDict(allocator, data) catch |err| switch (err) {
                error.DuplicateComponent => return,
                else => |e| return e,
            };
            dict.deinit();
            return error.TestUnexpectedResult;
        }
    };
    try std.testing.checkAllAllocationFailures(std.testing.allocator, CheckErr.run, .{dup.buf.items});
}

test "readDict and loadDict apply the same validation" {
    const allocator = std.testing.allocator;
    var atom = testAtom("CA", "C");
    atom.atom_id_len = 7;
    var img = TestImage{};
    defer img.deinit(allocator);
    try img.header(allocator, 1);
    try img.component(allocator, "ALA", &.{atom}, &.{});

    var reader = std.Io.Reader.fixed(img.buf.items);
    try std.testing.expectError(error.InvalidAtomLength, readDict(allocator, &reader));
    try std.testing.expectError(error.InvalidAtomLength, loadDict(allocator, img.buf.items));
}

fn addSizedTestComponent(
    allocator: Allocator,
    dict: *ccd_parser.ComponentDict,
    comp_id: []const u8,
    atom_count: usize,
    bond_count: usize,
) !void {
    const atoms = try allocator.alloc(CompAtom, atom_count);
    errdefer allocator.free(atoms);
    @memset(atoms, CompAtom.init("C", "C"));
    const bonds = try allocator.alloc(CompBond, bond_count);
    errdefer allocator.free(bonds);
    @memset(bonds, .{ .atom_idx_1 = 0, .atom_idx_2 = 0, .order = .single, .aromatic = false });
    const key = try allocator.dupe(u8, comp_id);
    errdefer allocator.free(key);
    var stored_id = [_]u8{ 0, 0, 0, 0, 0 };
    const id_len = @min(comp_id.len, stored_id.len);
    @memcpy(stored_id[0..id_len], comp_id[0..id_len]);
    try dict.components.put(allocator, key, .{
        .comp_id = stored_id,
        .comp_id_len = @intCast(id_len),
        .atoms = atoms,
        .bonds = bonds,
        .allocator = allocator,
    });
    try dict.owned_keys.append(allocator, key);
}

test "writeDict fails instead of wrapping atom and bond counts" {
    const allocator = std.testing.allocator;
    var sink: [64]u8 = undefined;

    var too_many_atoms = ccd_parser.ComponentDict.init(allocator);
    defer too_many_atoms.deinit();
    try addSizedTestComponent(allocator, &too_many_atoms, "BIG", @as(usize, std.math.maxInt(u16)) + 1, 0);
    var w = std.Io.Writer.fixed(&sink);
    var diag: WriteDiagnostic = .{};
    try std.testing.expectError(error.TooManyAtoms, writeDictDiag(&w, &too_many_atoms, &diag));
    try std.testing.expectEqualStrings("BIG", diag.comp_id);
    try std.testing.expectEqual(@as(usize, 0), w.end); // nothing written

    var too_many_bonds = ccd_parser.ComponentDict.init(allocator);
    defer too_many_bonds.deinit();
    try addSizedTestComponent(allocator, &too_many_bonds, "BND", 1, @as(usize, std.math.maxInt(u16)) + 1);
    w = std.Io.Writer.fixed(&sink);
    try std.testing.expectError(error.TooManyBonds, writeDict(&w, &too_many_bonds));
    try std.testing.expectEqual(@as(usize, 0), w.end);
}

test "writeDict accepts the largest component the format can hold" {
    const allocator = std.testing.allocator;
    var dict = ccd_parser.ComponentDict.init(allocator);
    defer dict.deinit();
    try addSizedTestComponent(allocator, &dict, "MAX", std.math.maxInt(u16), std.math.maxInt(u16));
    // Bonds must stay inside the component: keep the all-zero indices (atom 0).

    const bytes = try writeTestDictAlloc(allocator, &dict);
    defer allocator.free(bytes);
    var loaded = try parseDict(allocator, bytes);
    defer loaded.deinit();
    const comp = loaded.get("MAX") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, std.math.maxInt(u16)), comp.atoms.len);
    try std.testing.expectEqual(@as(usize, std.math.maxInt(u16)), comp.bonds.len);
}

test "writeDict rejects component IDs the length byte cannot hold and dangling bonds" {
    const allocator = std.testing.allocator;
    var sink: [64]u8 = undefined;

    var long_id = ccd_parser.ComponentDict.init(allocator);
    defer long_id.deinit();
    const long_name = [_]u8{'X'} ** 256;
    try addSizedTestComponent(allocator, &long_id, &long_name, 1, 0);
    var w = std.Io.Writer.fixed(&sink);
    try std.testing.expectError(error.InvalidComponentId, writeDict(&w, &long_id));

    var empty_id = ccd_parser.ComponentDict.init(allocator);
    defer empty_id.deinit();
    try addSizedTestComponent(allocator, &empty_id, "", 1, 0);
    w = std.Io.Writer.fixed(&sink);
    try std.testing.expectError(error.InvalidComponentId, writeDict(&w, &empty_id));

    var dangling = ccd_parser.ComponentDict.init(allocator);
    defer dangling.deinit();
    try addSizedTestComponent(allocator, &dangling, "BAD", 1, 1);
    dangling.components.getPtr("BAD").?.bonds[0].atom_idx_2 = 5;
    w = std.Io.Writer.fixed(&sink);
    try std.testing.expectError(error.InvalidBondIndex, writeDict(&w, &dangling));
}

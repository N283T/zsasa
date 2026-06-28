const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const mmap_reader = @import("mmap_reader.zig");
const compressed = @import("compressed.zig");
const types = @import("types.zig");
const ccd_parser = @import("ccd_parser.zig");
const hyb = @import("hybridization.zig");
const AtomInput = types.AtomInput;
pub const ParseError = error{
    InvalidMessagePack,
    UnexpectedEof,
    UnsupportedMessagePackType,
    MissingDataBlocks,
    NoAtomSiteCategory,
    MissingCoordinateField,
    UnsupportedEncoding,
    InvalidColumnData,
    InvalidCoordinate,
    ColumnLengthMismatch,
};

const MsgError = ParseError || Allocator.Error;
const DecodeError = ParseError || Allocator.Error;

const MsgValue = union(enum) {
    nil,
    bool: bool,
    int: i64,
    uint: u64,
    float: f64,
    str: []const u8,
    bin: []const u8,
    array: []MsgValue,
    map: []MsgPair,
};

const MsgPair = struct { key: MsgValue, value: MsgValue };

const MsgReader = struct {
    allocator: Allocator,
    data: []const u8,
    pos: usize = 0,

    fn init(allocator: Allocator, data: []const u8) MsgReader {
        return .{ .allocator = allocator, .data = data };
    }

    fn readValue(self: *MsgReader) MsgError!MsgValue {
        const marker = try self.readByte();

        if (marker <= 0x7f) return .{ .uint = marker };
        if (marker >= 0x80 and marker <= 0x8f) return self.readMap(marker & 0x0f);
        if (marker >= 0x90 and marker <= 0x9f) return self.readArray(marker & 0x0f);
        if (marker >= 0xa0 and marker <= 0xbf) return .{ .str = try self.readBytes(marker & 0x1f) };
        if (marker >= 0xe0) return .{ .int = @as(i8, @bitCast(marker)) };

        return switch (marker) {
            0xc0 => .nil,
            0xc2 => .{ .bool = false },
            0xc3 => .{ .bool = true },
            0xc4 => .{ .bin = try self.readBytes(try self.readU8()) },
            0xc5 => .{ .bin = try self.readBytes(try self.readU16()) },
            0xc6 => .{ .bin = try self.readBytes(try self.readU32AsUsize()) },
            0xca => .{ .float = try self.readF32() },
            0xcb => .{ .float = try self.readF64() },
            0xcc => .{ .uint = try self.readU8() },
            0xcd => .{ .uint = try self.readU16() },
            0xce => .{ .uint = try self.readU32() },
            0xcf => .{ .uint = try self.readU64() },
            0xd0 => .{ .int = try self.readI8() },
            0xd1 => .{ .int = try self.readI16() },
            0xd2 => .{ .int = try self.readI32() },
            0xd3 => .{ .int = try self.readI64() },
            0xd9 => .{ .str = try self.readBytes(try self.readU8()) },
            0xda => .{ .str = try self.readBytes(try self.readU16()) },
            0xdb => .{ .str = try self.readBytes(try self.readU32AsUsize()) },
            0xdc => self.readArray(try self.readU16()),
            0xdd => self.readArray(try self.readU32AsUsize()),
            0xde => self.readMap(try self.readU16()),
            0xdf => self.readMap(try self.readU32AsUsize()),
            else => ParseError.UnsupportedMessagePackType,
        };
    }

    fn readArray(self: *MsgReader, len: usize) MsgError!MsgValue {
        const items = try self.allocator.alloc(MsgValue, len);
        var initialized: usize = 0;
        errdefer {
            for (items[0..initialized]) |item| freeMsgValue(self.allocator, item);
            self.allocator.free(items);
        }

        while (initialized < len) : (initialized += 1) {
            items[initialized] = try self.readValue();
        }

        return .{ .array = items };
    }

    fn readMap(self: *MsgReader, len: usize) MsgError!MsgValue {
        const pairs = try self.allocator.alloc(MsgPair, len);
        var initialized: usize = 0;
        errdefer {
            for (pairs[0..initialized]) |pair| {
                freeMsgValue(self.allocator, pair.key);
                freeMsgValue(self.allocator, pair.value);
            }
            self.allocator.free(pairs);
        }

        while (initialized < len) : (initialized += 1) {
            const key = try self.readValue();
            const value = self.readValue() catch |err| {
                freeMsgValue(self.allocator, key);
                return err;
            };
            pairs[initialized] = .{ .key = key, .value = value };
        }

        return .{ .map = pairs };
    }

    fn readByte(self: *MsgReader) MsgError!u8 {
        if (self.pos >= self.data.len) return ParseError.UnexpectedEof;
        const byte = self.data[self.pos];
        self.pos += 1;
        return byte;
    }

    fn readBytes(self: *MsgReader, len: usize) MsgError![]const u8 {
        if (len > self.data.len -| self.pos) return ParseError.UnexpectedEof;
        const bytes = self.data[self.pos .. self.pos + len];
        self.pos += len;
        return bytes;
    }

    fn readU8(self: *MsgReader) MsgError!u8 {
        return try self.readByte();
    }

    fn readU16(self: *MsgReader) MsgError!u16 {
        const bytes = try self.readBytes(2);
        return std.mem.readInt(u16, bytes[0..2], .big);
    }

    fn readU32(self: *MsgReader) MsgError!u32 {
        const bytes = try self.readBytes(4);
        return std.mem.readInt(u32, bytes[0..4], .big);
    }

    fn readU32AsUsize(self: *MsgReader) MsgError!usize {
        return std.math.cast(usize, try self.readU32()) orelse ParseError.InvalidMessagePack;
    }

    fn readU64(self: *MsgReader) MsgError!u64 {
        const bytes = try self.readBytes(8);
        return std.mem.readInt(u64, bytes[0..8], .big);
    }

    fn readI8(self: *MsgReader) MsgError!i8 {
        return @bitCast(try self.readByte());
    }

    fn readI16(self: *MsgReader) MsgError!i16 {
        const bytes = try self.readBytes(2);
        return std.mem.readInt(i16, bytes[0..2], .big);
    }

    fn readI32(self: *MsgReader) MsgError!i32 {
        const bytes = try self.readBytes(4);
        return std.mem.readInt(i32, bytes[0..4], .big);
    }

    fn readI64(self: *MsgReader) MsgError!i64 {
        const bytes = try self.readBytes(8);
        return std.mem.readInt(i64, bytes[0..8], .big);
    }

    fn readF32(self: *MsgReader) MsgError!f64 {
        const bits = try self.readU32();
        return @as(f64, @floatCast(@as(f32, @bitCast(bits))));
    }

    fn readF64(self: *MsgReader) MsgError!f64 {
        const bits = try self.readU64();
        return @bitCast(bits);
    }
};

fn freeMsgValue(allocator: Allocator, value: MsgValue) void {
    switch (value) {
        .array => |items| {
            for (items) |item| freeMsgValue(allocator, item);
            allocator.free(items);
        },
        .map => |pairs| {
            for (pairs) |pair| {
                freeMsgValue(allocator, pair.key);
                freeMsgValue(allocator, pair.value);
            }
            allocator.free(pairs);
        },
        else => {},
    }
}

fn getMapValue(map: []const MsgPair, key: []const u8) ?MsgValue {
    for (map) |pair| {
        switch (pair.key) {
            .str => |candidate| {
                if (std.mem.eql(u8, candidate, key)) return pair.value;
            },
            else => {},
        }
    }
    return null;
}

const NullKind = enum {
    present,
    not_present,
    unknown,
};

const Scalar = union(enum) {
    int: i64,
    float: f64,
    str: []const u8,
};

const DecodedColumn = struct {
    values: []Scalar,
    nulls: ?[]NullKind = null,

    fn deinit(self: *DecodedColumn, allocator: Allocator) void {
        allocator.free(self.values);
        if (self.nulls) |nulls| allocator.free(nulls);
        self.* = .{ .values = &.{}, .nulls = null };
    }
};

const Encoding = union(enum) {
    byte_array: struct { type_code: u8 },
    integer_packing: struct { byte_count: u8, is_unsigned: bool, src_size: usize },
    run_length: struct { src_size: usize },
    delta: struct { origin: i64 },
    fixed_point: struct { factor: f64 },
    interval_quantization: struct { min: f64, max: f64, num_steps: i64 },
    string_array: struct {
        string_data: []const u8,
        offset_data: []const u8,
        offset_encoding: []const Encoding,
        data_encoding: []const Encoding,
    },
};

const PackingSentinel = struct { positive: i64, negative: i64 };

fn decodeColumn(allocator: Allocator, data: []const u8, encodings: []const Encoding) DecodeError!DecodedColumn {
    return decodeColumnWithNullHints(allocator, data, encodings, null);
}

fn decodeColumnWithNullHints(
    allocator: Allocator,
    data: []const u8,
    encodings: []const Encoding,
    null_hints: ?[]const NullKind,
) DecodeError!DecodedColumn {
    var column = DecodedColumn{ .values = &.{} };
    var initialized = false;
    errdefer if (initialized) column.deinit(allocator);

    var i = encodings.len;
    while (i > 0) {
        i -= 1;
        if (!initialized) {
            column = try decodeInitialEncoding(allocator, data, encodings[i], null_hints);
            initialized = true;
        } else {
            initialized = false;
            column = try applyEncoding(allocator, column, encodings[i]);
            initialized = true;
        }
    }

    if (!initialized) {
        column.values = try allocator.alloc(Scalar, data.len);
        initialized = true;
        for (data, 0..) |byte, idx| column.values[idx] = .{ .int = byte };
    }
    return column;
}

fn decodeColumnWithMask(
    allocator: Allocator,
    data: []const u8,
    encodings: []const Encoding,
    mask_data: []const u8,
    mask_encodings: []const Encoding,
) DecodeError!DecodedColumn {
    var mask = try decodeColumn(allocator, mask_data, mask_encodings);
    defer mask.deinit(allocator);

    const nulls = try allocator.alloc(NullKind, mask.values.len);
    errdefer allocator.free(nulls);
    for (mask.values, nulls) |mask_value, *null_kind| {
        const value = scalarToI64(mask_value) orelse return ParseError.InvalidColumnData;
        null_kind.* = switch (value) {
            0 => .present,
            1 => .not_present,
            2 => .unknown,
            else => return ParseError.InvalidColumnData,
        };
    }

    var column = try decodeColumnWithNullHints(allocator, data, encodings, nulls);
    errdefer column.deinit(allocator);
    if (mask.values.len != column.values.len) return ParseError.ColumnLengthMismatch;

    column.nulls = nulls;
    return column;
}

fn decodeInitialEncoding(allocator: Allocator, data: []const u8, encoding: Encoding, null_hints: ?[]const NullKind) DecodeError!DecodedColumn {
    return switch (encoding) {
        .byte_array => |params| decodeByteArray(allocator, data, params.type_code),
        .string_array => |params| decodeStringArray(allocator, data, params, null_hints),
        else => ParseError.UnsupportedEncoding,
    };
}

fn applyEncoding(allocator: Allocator, column: DecodedColumn, encoding: Encoding) DecodeError!DecodedColumn {
    return switch (encoding) {
        .integer_packing => |params| decodeIntegerPacking(allocator, column, params.byte_count, params.is_unsigned, params.src_size),
        .run_length => |params| decodeRunLength(allocator, column, params.src_size),
        .delta => |params| decodeDelta(allocator, column, params.origin),
        .fixed_point => |params| decodeFixedPoint(allocator, column, params.factor),
        .interval_quantization => |params| decodeIntervalQuantization(allocator, column, params.min, params.max, params.num_steps),
        .byte_array, .string_array => {
            var old = column;
            old.deinit(allocator);
            return ParseError.UnsupportedEncoding;
        },
    };
}

fn decodeByteArray(allocator: Allocator, data: []const u8, type_code: u8) DecodeError!DecodedColumn {
    const item_size = byteArrayItemSize(type_code) orelse return ParseError.UnsupportedEncoding;
    if (data.len % item_size != 0) return ParseError.InvalidColumnData;
    const values = try allocator.alloc(Scalar, data.len / item_size);
    errdefer allocator.free(values);

    for (values, 0..) |*value, idx| {
        const bytes = data[idx * item_size ..][0..item_size];
        value.* = switch (type_code) {
            1 => .{ .int = std.mem.readInt(i8, bytes[0..1], .little) },
            2 => .{ .int = std.mem.readInt(i16, bytes[0..2], .little) },
            3 => .{ .int = std.mem.readInt(i32, bytes[0..4], .little) },
            4 => .{ .int = std.mem.readInt(u8, bytes[0..1], .little) },
            5 => .{ .int = std.mem.readInt(u16, bytes[0..2], .little) },
            6 => .{ .int = std.mem.readInt(u32, bytes[0..4], .little) },
            32 => .{ .float = @as(f64, @floatCast(@as(f32, @bitCast(std.mem.readInt(u32, bytes[0..4], .little))))) },
            33 => .{ .float = @bitCast(std.mem.readInt(u64, bytes[0..8], .little)) },
            else => unreachable,
        };
    }
    return .{ .values = values };
}

fn byteArrayItemSize(type_code: u8) ?usize {
    return switch (type_code) {
        1, 4 => 1,
        2, 5 => 2,
        3, 6, 32 => 4,
        33 => 8,
        else => null,
    };
}

fn decodeIntegerPacking(allocator: Allocator, column: DecodedColumn, byte_count: u8, is_unsigned: bool, src_size: usize) DecodeError!DecodedColumn {
    defer {
        var old = column;
        old.deinit(allocator);
    }

    const sentinel = packingSentinel(byte_count, is_unsigned) orelse return ParseError.UnsupportedEncoding;
    var out = try std.ArrayListUnmanaged(Scalar).initCapacity(allocator, src_size);
    errdefer out.deinit(allocator);

    var acc: i64 = 0;
    for (column.values) |value| {
        const part = scalarToI64(value) orelse return ParseError.InvalidColumnData;
        try appendIntegerPackingPart(allocator, &out, &acc, part, sentinel);
    }
    if (out.items.len != src_size) return ParseError.InvalidColumnData;
    return .{ .values = try out.toOwnedSlice(allocator) };
}

fn appendIntegerPackingPart(
    allocator: Allocator,
    out: *std.ArrayListUnmanaged(Scalar),
    acc: *i64,
    part: i64,
    sentinel: PackingSentinel,
) DecodeError!void {
    acc.* = checkedAddI64(acc.*, part) catch return ParseError.InvalidColumnData;
    if (part != sentinel.positive and part != sentinel.negative) {
        try out.append(allocator, .{ .int = acc.* });
        acc.* = 0;
    }
}

fn packingSentinel(byte_count: u8, is_unsigned: bool) ?PackingSentinel {
    return switch (byte_count) {
        1 => if (is_unsigned)
            .{ .positive = std.math.maxInt(u8), .negative = -1 }
        else
            .{ .positive = std.math.maxInt(i8), .negative = std.math.minInt(i8) },
        2 => if (is_unsigned)
            .{ .positive = std.math.maxInt(u16), .negative = -1 }
        else
            .{ .positive = std.math.maxInt(i16), .negative = std.math.minInt(i16) },
        4 => if (is_unsigned)
            .{ .positive = std.math.maxInt(u32), .negative = -1 }
        else
            .{ .positive = std.math.maxInt(i32), .negative = std.math.minInt(i32) },
        else => null,
    };
}

fn decodeRunLength(allocator: Allocator, column: DecodedColumn, src_size: usize) DecodeError!DecodedColumn {
    defer {
        var old = column;
        old.deinit(allocator);
    }
    if (column.values.len % 2 != 0) return ParseError.InvalidColumnData;

    const values = try allocator.alloc(Scalar, src_size);
    errdefer allocator.free(values);
    var out_idx: usize = 0;
    var in_idx: usize = 0;
    while (in_idx < column.values.len) : (in_idx += 2) {
        const repeat_value = column.values[in_idx];
        const repeat_count = scalarToUsize(column.values[in_idx + 1]) orelse return ParseError.InvalidColumnData;
        if (repeat_count > src_size - out_idx) return ParseError.InvalidColumnData;
        for (values[out_idx .. out_idx + repeat_count]) |*value| value.* = repeat_value;
        out_idx += repeat_count;
    }
    if (out_idx != src_size) return ParseError.InvalidColumnData;
    return .{ .values = values };
}

fn decodeDelta(allocator: Allocator, column: DecodedColumn, origin: i64) DecodeError!DecodedColumn {
    defer {
        var old = column;
        old.deinit(allocator);
    }
    const values = try allocator.alloc(Scalar, column.values.len);
    errdefer allocator.free(values);

    var current = origin;
    for (column.values, values) |value, *out| {
        current = checkedAddI64(current, scalarToI64(value) orelse return ParseError.InvalidColumnData) catch return ParseError.InvalidColumnData;
        out.* = .{ .int = current };
    }
    return .{ .values = values };
}

fn decodeFixedPoint(allocator: Allocator, column: DecodedColumn, factor: f64) DecodeError!DecodedColumn {
    defer {
        var old = column;
        old.deinit(allocator);
    }
    if (factor == 0.0) return ParseError.InvalidColumnData;
    const values = try allocator.alloc(Scalar, column.values.len);
    errdefer allocator.free(values);

    for (column.values, values) |value, *out| {
        const int_value = scalarToI64(value) orelse return ParseError.InvalidColumnData;
        out.* = .{ .float = @as(f64, @floatFromInt(int_value)) / factor };
    }
    return .{ .values = values };
}

fn decodeIntervalQuantization(allocator: Allocator, column: DecodedColumn, min: f64, max: f64, num_steps: i64) DecodeError!DecodedColumn {
    defer {
        var old = column;
        old.deinit(allocator);
    }
    if (num_steps <= 1) return ParseError.InvalidColumnData;
    const values = try allocator.alloc(Scalar, column.values.len);
    errdefer allocator.free(values);

    const delta = (max - min) / @as(f64, @floatFromInt(num_steps - 1));
    for (column.values, values) |value, *out| {
        const int_value = scalarToI64(value) orelse return ParseError.InvalidColumnData;
        out.* = .{ .float = min + delta * @as(f64, @floatFromInt(int_value)) };
    }
    return .{ .values = values };
}

fn decodeStringArray(allocator: Allocator, data: []const u8, params: @FieldType(Encoding, "string_array"), null_hints: ?[]const NullKind) DecodeError!DecodedColumn {
    var indices = try decodeColumn(allocator, data, params.data_encoding);
    defer indices.deinit(allocator);
    var offsets = try decodeColumn(allocator, params.offset_data, params.offset_encoding);
    defer offsets.deinit(allocator);
    if (offsets.values.len == 0) return ParseError.InvalidColumnData;

    const values = try allocator.alloc(Scalar, indices.values.len);
    errdefer allocator.free(values);
    for (indices.values, values, 0..) |index_value, *out, row| {
        const index_int = scalarToI64(index_value) orelse return ParseError.InvalidColumnData;
        if (index_int < 0) {
            const hints = null_hints orelse return ParseError.InvalidColumnData;
            if (row >= hints.len or hints[row] == .present) return ParseError.InvalidColumnData;
            out.* = .{ .str = "" };
            continue;
        }
        const index = std.math.cast(usize, index_int) orelse return ParseError.InvalidColumnData;
        if (index + 1 >= offsets.values.len) return ParseError.InvalidColumnData;
        const start = scalarToUsize(offsets.values[index]) orelse return ParseError.InvalidColumnData;
        var end = scalarToUsize(offsets.values[index + 1]) orelse return ParseError.InvalidColumnData;
        if (index + 1 == offsets.values.len - 1 and end > params.string_data.len) end = params.string_data.len;
        if (end < start or end > params.string_data.len) return ParseError.InvalidColumnData;
        out.* = .{ .str = params.string_data[start..end] };
    }
    return .{ .values = values };
}

fn scalarToI64(value: Scalar) ?i64 {
    return switch (value) {
        .int => |int| int,
        else => null,
    };
}

fn scalarToUsize(value: Scalar) ?usize {
    const int = scalarToI64(value) orelse return null;
    return std.math.cast(usize, int);
}

fn checkedAddI64(a: i64, b: i64) error{Overflow}!i64 {
    return std.math.add(i64, a, b) catch error.Overflow;
}

const AtomSiteColumns = struct {
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    type_symbol: ?usize = null,
    label_atom_id: ?usize = null,
    auth_atom_id: ?usize = null,
    label_comp_id: ?usize = null,
    auth_comp_id: ?usize = null,
    label_asym_id: ?usize = null,
    auth_asym_id: ?usize = null,
    label_seq_id: ?usize = null,
    auth_seq_id: ?usize = null,
    pdbx_pdb_ins_code: ?usize = null,
    group_pdb: ?usize = null,
    label_alt_id: ?usize = null,
    occupancy: ?usize = null,
    pdbx_pdb_model_num: ?usize = null,

    fn hasRequiredFields(self: AtomSiteColumns) bool {
        return self.cartn_x != null and self.cartn_y != null and self.cartn_z != null;
    }

    fn getAtomNameCol(self: AtomSiteColumns) ?usize {
        return self.label_atom_id orelse self.auth_atom_id;
    }

    fn getResNameCol(self: AtomSiteColumns) ?usize {
        return self.label_comp_id orelse self.auth_comp_id;
    }

    fn getChainCol(self: AtomSiteColumns, use_auth: bool) ?usize {
        if (use_auth) return self.auth_asym_id orelse self.label_asym_id;
        return self.label_asym_id orelse self.auth_asym_id;
    }

    fn getResSeqCol(self: AtomSiteColumns) ?usize {
        return self.label_seq_id orelse self.auth_seq_id;
    }

    fn getInsCodeCol(self: AtomSiteColumns) ?usize {
        return self.pdbx_pdb_ins_code;
    }
};

const CcdAtomColumns = struct {
    comp_id: ?usize = null,
    atom_id: ?usize = null,
    type_symbol: ?usize = null,
    pdbx_aromatic_flag: ?usize = null,
    pdbx_leaving_atom_flag: ?usize = null,

    fn hasRequiredFields(self: CcdAtomColumns) bool {
        return self.comp_id != null and self.atom_id != null and self.type_symbol != null;
    }
};

const CcdBondColumns = struct {
    comp_id: ?usize = null,
    atom_id_1: ?usize = null,
    atom_id_2: ?usize = null,
    value_order: ?usize = null,
    pdbx_aromatic_flag: ?usize = null,

    fn hasRequiredFields(self: CcdBondColumns) bool {
        return self.comp_id != null and self.atom_id_1 != null and self.atom_id_2 != null;
    }
};

const BcifCcdBuilder = struct {
    atoms: std.ArrayListUnmanaged(hyb.CompAtom) = .empty,
    bonds: std.ArrayListUnmanaged(hyb.CompBond) = .empty,
    atom_name_map: std.StringHashMapUnmanaged(u16) = .empty,

    fn deinit(self: *BcifCcdBuilder, allocator: Allocator) void {
        self.atoms.deinit(allocator);
        self.bonds.deinit(allocator);
        var it = self.atom_name_map.iterator();
        while (it.next()) |entry| allocator.free(entry.key_ptr.*);
        self.atom_name_map.deinit(allocator);
    }
};

pub const BcifParser = struct {
    allocator: Allocator,
    atom_only: bool = true,
    skip_hydrogens: bool = true,
    first_alt_loc_only: bool = true,
    model_num: ?u32 = null,
    chain_filter: ?[]const []const u8 = null,
    use_auth_chain: bool = false,
    inline_ccd: ?ccd_parser.ComponentDict = null,

    pub fn init(allocator: Allocator) BcifParser {
        return .{ .allocator = allocator };
    }

    /// Return a pointer to inline CCD data decoded from `_chem_comp_atom` and
    /// `_chem_comp_bond` BinaryCIF categories, or null when no components were
    /// present.
    pub fn getInlineCcd(self: *BcifParser) ?*ccd_parser.ComponentDict {
        if (self.inline_ccd != null) return &self.inline_ccd.?;
        return null;
    }

    /// Clean up inline CCD data. Must be called before the parser goes out of
    /// scope when `parse()` has been called and `getInlineCcd()` is non-null.
    pub fn deinitCcd(self: *BcifParser) void {
        if (self.inline_ccd) |*dict| {
            dict.deinit();
            self.inline_ccd = null;
        }
    }

    pub fn parse(self: *BcifParser, source: []const u8) !AtomInput {
        self.deinitCcd();

        var reader = MsgReader.init(self.allocator, source);
        const root = try reader.readValue();
        defer freeMsgValue(self.allocator, root);

        var inline_ccd = try parseInlineCcd(self.allocator, root);
        if (inline_ccd.count() == 0) {
            inline_ccd.deinit();
            self.inline_ccd = null;
        } else {
            self.inline_ccd = inline_ccd;
        }
        errdefer self.deinitCcd();

        const category = try findCategory(root, "atom_site");
        const row_count = try categoryRowCount(category);
        const raw_columns = try categoryColumns(category);

        var columns = AtomSiteColumns{};
        for (raw_columns, 0..) |column, idx| {
            const name = try columnName(column);
            mapAtomSiteColumn(&columns, name, idx);
        }
        if (!columns.hasRequiredFields()) {
            self.deinitCcd();
            return ParseError.MissingCoordinateField;
        }

        var decoded = try self.allocator.alloc(DecodedColumn, raw_columns.len);
        var decoded_count: usize = 0;
        errdefer {
            for (decoded[0..decoded_count]) |*column| column.deinit(self.allocator);
            self.allocator.free(decoded);
        }
        for (raw_columns) |column| {
            decoded[decoded_count] = decodeBcifColumn(self.allocator, column) catch |err| {
                self.deinitCcd();
                return err;
            };
            decoded_count += 1;
            if (decoded[decoded_count - 1].values.len != row_count) {
                self.deinitCcd();
                return ParseError.ColumnLengthMismatch;
            }
        }
        defer {
            for (decoded[0..decoded_count]) |*column| column.deinit(self.allocator);
            self.allocator.free(decoded);
        }

        var x_list = std.ArrayListUnmanaged(f64).empty;
        defer x_list.deinit(self.allocator);
        var y_list = std.ArrayListUnmanaged(f64).empty;
        defer y_list.deinit(self.allocator);
        var z_list = std.ArrayListUnmanaged(f64).empty;
        defer z_list.deinit(self.allocator);
        var r_list = std.ArrayListUnmanaged(f64).empty;
        defer r_list.deinit(self.allocator);
        var residue_list = std.ArrayListUnmanaged(types.FixedString5).empty;
        defer residue_list.deinit(self.allocator);
        var atom_name_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer atom_name_list.deinit(self.allocator);
        var element_list = std.ArrayListUnmanaged(u8).empty;
        defer element_list.deinit(self.allocator);
        var chain_id_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer chain_id_list.deinit(self.allocator);
        var chain_id_full_list = std.ArrayListUnmanaged([]const u8).empty;
        defer chain_id_full_list.deinit(self.allocator);
        errdefer freeStringItems(self.allocator, chain_id_full_list.items);
        var has_extended_chain = false;
        var residue_num_list = std.ArrayListUnmanaged(i32).empty;
        defer residue_num_list.deinit(self.allocator);
        var insertion_code_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer insertion_code_list.deinit(self.allocator);
        var atom_records = std.ArrayListUnmanaged(AtomRecord).empty;
        defer atom_records.deinit(self.allocator);

        var row: usize = 0;
        while (row < row_count) : (row += 1) {
            if (!self.shouldIncludeAtom(decoded, columns, row)) continue;
            try atom_records.append(self.allocator, try self.atomRecordFromRow(decoded, columns, row));
        }

        for (atom_records.items, 0..) |atom, i| {
            if (!self.shouldKeepAltLoc(atom_records.items, i)) continue;

            try x_list.append(self.allocator, atom.x);
            try y_list.append(self.allocator, atom.y);
            try z_list.append(self.allocator, atom.z);
            try r_list.append(self.allocator, atom.radius);
            try element_list.append(self.allocator, atom.element.atomicNumber());
            try residue_list.append(self.allocator, types.FixedString5.fromSlice(atom.residue));
            try atom_name_list.append(self.allocator, types.FixedString4.fromSlice(atom.atom_name));
            try chain_id_list.append(self.allocator, types.FixedString4.fromSlice(atom.chain_id));
            try chain_id_full_list.append(self.allocator, try self.allocator.dupe(u8, atom.chain_id));
            has_extended_chain = has_extended_chain or atom.chain_id.len > 4;
            try residue_num_list.append(self.allocator, atom.residue_num);
            try insertion_code_list.append(self.allocator, types.FixedString4.fromSlice(atom.insertion_code));
        }

        const x = try x_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(x);
        const y = try y_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(y);
        const z = try z_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(z);
        const r = try r_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(r);
        const residue = try residue_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(residue);
        const atom_name = try atom_name_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(atom_name);
        const element = try element_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(element);
        const chain_id = try chain_id_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(chain_id);
        const chain_id_full_owned = try chain_id_full_list.toOwnedSlice(self.allocator);
        const chain_id_full: ?[]const []const u8 = if (has_extended_chain) chain_id_full_owned else blk: {
            freeStringItems(self.allocator, chain_id_full_owned);
            self.allocator.free(chain_id_full_owned);
            break :blk null;
        };
        errdefer if (chain_id_full) |chains| {
            freeStringItems(self.allocator, chains);
            self.allocator.free(chains);
        };
        const residue_num = try residue_num_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(residue_num);
        const insertion_code = try insertion_code_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(insertion_code);

        return AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = r,
            .residue = residue,
            .atom_name = atom_name,
            .element = element,
            .chain_id = chain_id,
            .chain_id_full = chain_id_full,
            .residue_num = residue_num,
            .insertion_code = insertion_code,
            .allocator = self.allocator,
        };
    }

    pub fn parseFile(self: *BcifParser, io: std.Io, path: []const u8) !AtomInput {
        if (compressed.isCompressed(path)) {
            const data = try compressed.read(self.allocator, path);
            defer self.allocator.free(data);
            return self.parse(data);
        }
        const mapped = try mmap_reader.mmapFile(self.allocator, io, path);
        defer mapped.deinit();
        return self.parse(mapped.data);
    }

    const AtomRecord = struct {
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        element: elem.Element,
        atom_name: []const u8,
        residue: []const u8,
        chain_id: []const u8,
        residue_num: i32,
        insertion_code: []const u8,
        alt_loc: u8,
        occupancy: f64,
        model_num: ?u32,
    };

    fn atomRecordFromRow(self: *BcifParser, decoded: []const DecodedColumn, columns: AtomSiteColumns, row: usize) !AtomRecord {
        const x = try scalarFloat(decoded[columns.cartn_x.?], row);
        const y = try scalarFloat(decoded[columns.cartn_y.?], row);
        const z = try scalarFloat(decoded[columns.cartn_z.?], row);

        var element_enum = elem.Element.X;
        if (columns.type_symbol) |type_col| {
            if (scalarString(decoded[type_col], row)) |symbol| {
                element_enum = elem.fromSymbol(symbol);
            }
        }

        const atom_name = if (columns.getAtomNameCol()) |atom_col| scalarString(decoded[atom_col], row) orelse "X" else "X";
        const residue = if (columns.getResNameCol()) |res_col| scalarString(decoded[res_col], row) orelse "UNK" else "UNK";
        const chain_id = if (columns.getChainCol(self.use_auth_chain)) |chain_col| scalarString(decoded[chain_col], row) orelse "" else "";
        const residue_num = if (columns.getResSeqCol()) |seq_col| blk: {
            const seq = scalarInt(decoded[seq_col], row) orelse 0;
            break :blk std.math.cast(i32, seq) orelse 0;
        } else 0;
        const insertion_code = if (columns.getInsCodeCol()) |ins_col| scalarString(decoded[ins_col], row) orelse "" else "";
        const alt_loc: u8 = if (columns.label_alt_id) |alt_col| blk: {
            const alt_id = scalarString(decoded[alt_col], row) orelse "";
            break :blk if (alt_id.len == 0) ' ' else alt_id[0];
        } else ' ';
        const occupancy = if (columns.occupancy) |occ_col| scalarFloat(decoded[occ_col], row) catch 0.0 else 0.0;
        const model_num = if (columns.pdbx_pdb_model_num) |model_col| blk: {
            const model = scalarInt(decoded[model_col], row) orelse break :blk null;
            break :blk std.math.cast(u32, model);
        } else null;

        return .{
            .x = x,
            .y = y,
            .z = z,
            .radius = element_enum.vdwRadius(),
            .element = element_enum,
            .atom_name = atom_name,
            .residue = residue,
            .chain_id = chain_id,
            .residue_num = residue_num,
            .insertion_code = insertion_code,
            .alt_loc = alt_loc,
            .occupancy = occupancy,
            .model_num = model_num,
        };
    }

    fn sameAltLocSite(a: AtomRecord, b: AtomRecord) bool {
        return a.model_num == b.model_num and
            a.residue_num == b.residue_num and
            std.mem.eql(u8, a.chain_id, b.chain_id) and
            std.mem.eql(u8, a.residue, b.residue) and
            std.mem.eql(u8, a.insertion_code, b.insertion_code) and
            std.mem.eql(u8, a.atom_name, b.atom_name);
    }

    fn shouldKeepAltLoc(self: *BcifParser, atoms: []const AtomRecord, index: usize) bool {
        if (!self.first_alt_loc_only) return true;

        const atom = atoms[index];
        if (atom.alt_loc == ' ') return true;

        var best_non_preferred: ?usize = null;
        for (atoms, 0..) |other, other_index| {
            if (!sameAltLocSite(atom, other)) continue;
            if (other.alt_loc == ' ') return false;
            if (other.alt_loc == 'A') return atom.alt_loc == 'A';
            if (best_non_preferred) |best_index| {
                if (other.occupancy > atoms[best_index].occupancy) {
                    best_non_preferred = other_index;
                }
            } else {
                best_non_preferred = other_index;
            }
        }
        return best_non_preferred == index;
    }

    fn shouldIncludeAtom(self: *BcifParser, decoded: []const DecodedColumn, columns: AtomSiteColumns, row: usize) bool {
        if (self.atom_only) {
            if (columns.group_pdb) |col| {
                const group = scalarString(decoded[col], row) orelse return false;
                if (!std.mem.eql(u8, group, "ATOM")) return false;
            }
        }

        if (self.skip_hydrogens) {
            if (columns.type_symbol) |col| {
                if (scalarString(decoded[col], row)) |symbol| {
                    if (std.mem.eql(u8, symbol, "H") or std.mem.eql(u8, symbol, "D")) return false;
                }
            }
        }

        if (self.model_num) |target_model| {
            if (columns.pdbx_pdb_model_num) |col| {
                if (scalarInt(decoded[col], row)) |model| {
                    if (model != target_model) return false;
                }
            }
        }

        if (self.chain_filter) |chains| {
            if (columns.getChainCol(self.use_auth_chain)) |col| {
                const chain = scalarString(decoded[col], row) orelse return false;
                for (chains) |target_chain| {
                    if (std.mem.eql(u8, chain, target_chain)) return true;
                }
                return false;
            }
        }

        return true;
    }
};

fn freeStringItems(allocator: Allocator, items: []const []const u8) void {
    for (items) |item| allocator.free(item);
}

fn parseEncoding(allocator: Allocator, value: MsgValue) DecodeError!Encoding {
    const map = switch (value) {
        .map => |map| map,
        else => return ParseError.InvalidColumnData,
    };
    const kind = try msgString(getMapValue(map, "kind") orelse return ParseError.InvalidColumnData);
    if (std.mem.eql(u8, kind, "ByteArray")) {
        const type_code = std.math.cast(u8, try msgInt(getMapValue(map, "type") orelse return ParseError.InvalidColumnData)) orelse return ParseError.InvalidColumnData;
        return .{ .byte_array = .{ .type_code = type_code } };
    }
    if (std.mem.eql(u8, kind, "FixedPoint")) {
        return .{ .fixed_point = .{ .factor = try msgFloat(getMapValue(map, "factor") orelse return ParseError.InvalidColumnData) } };
    }
    if (std.mem.eql(u8, kind, "IntervalQuantization")) {
        return .{ .interval_quantization = .{
            .min = try msgFloat(getMapValue(map, "min") orelse return ParseError.InvalidColumnData),
            .max = try msgFloat(getMapValue(map, "max") orelse return ParseError.InvalidColumnData),
            .num_steps = try msgInt(getMapValue(map, "numSteps") orelse return ParseError.InvalidColumnData),
        } };
    }
    if (std.mem.eql(u8, kind, "RunLength")) {
        const src_size = std.math.cast(usize, try msgInt(getMapValue(map, "srcSize") orelse return ParseError.InvalidColumnData)) orelse return ParseError.InvalidColumnData;
        return .{ .run_length = .{ .src_size = src_size } };
    }
    if (std.mem.eql(u8, kind, "Delta")) {
        return .{ .delta = .{ .origin = try msgInt(getMapValue(map, "origin") orelse return ParseError.InvalidColumnData) } };
    }
    if (std.mem.eql(u8, kind, "IntegerPacking")) {
        const byte_count = std.math.cast(u8, try msgInt(getMapValue(map, "byteCount") orelse return ParseError.InvalidColumnData)) orelse return ParseError.InvalidColumnData;
        const src_size = std.math.cast(usize, try msgInt(getMapValue(map, "srcSize") orelse return ParseError.InvalidColumnData)) orelse return ParseError.InvalidColumnData;
        return .{ .integer_packing = .{
            .byte_count = byte_count,
            .is_unsigned = try msgBool(getMapValue(map, "isUnsigned") orelse return ParseError.InvalidColumnData),
            .src_size = src_size,
        } };
    }
    if (std.mem.eql(u8, kind, "StringArray")) {
        const data_encoding = try parseEncodingList(allocator, getMapValue(map, "dataEncoding") orelse return ParseError.InvalidColumnData);
        errdefer freeEncodingList(allocator, data_encoding);
        const offset_encoding = try parseEncodingList(allocator, getMapValue(map, "offsetEncoding") orelse return ParseError.InvalidColumnData);
        errdefer freeEncodingList(allocator, offset_encoding);
        return .{ .string_array = .{
            .string_data = try msgString(getMapValue(map, "stringData") orelse return ParseError.InvalidColumnData),
            .offset_data = try msgBytes(getMapValue(map, "offsets") orelse return ParseError.InvalidColumnData),
            .offset_encoding = offset_encoding,
            .data_encoding = data_encoding,
        } };
    }
    return ParseError.UnsupportedEncoding;
}

fn parseEncodingList(allocator: Allocator, value: MsgValue) DecodeError![]Encoding {
    const array = switch (value) {
        .array => |array| array,
        else => return ParseError.InvalidColumnData,
    };
    const encodings = try allocator.alloc(Encoding, array.len);
    var initialized: usize = 0;
    errdefer {
        for (encodings[0..initialized]) |encoding| freeEncoding(allocator, encoding);
        allocator.free(encodings);
    }
    for (array) |item| {
        encodings[initialized] = try parseEncoding(allocator, item);
        initialized += 1;
    }
    return encodings;
}

fn freeEncodingList(allocator: Allocator, encodings: []const Encoding) void {
    for (encodings) |encoding| freeEncoding(allocator, encoding);
    allocator.free(encodings);
}

fn freeEncoding(allocator: Allocator, encoding: Encoding) void {
    switch (encoding) {
        .string_array => |params| {
            freeEncodingList(allocator, params.offset_encoding);
            freeEncodingList(allocator, params.data_encoding);
        },
        else => {},
    }
}

fn findCategory(root: MsgValue, name: []const u8) !MsgValue {
    return (try findCategoryOptional(root, name)) orelse ParseError.NoAtomSiteCategory;
}

fn findCategoryOptional(root: MsgValue, name: []const u8) !?MsgValue {
    const root_map = switch (root) {
        .map => |map| map,
        else => return ParseError.InvalidMessagePack,
    };
    const blocks = switch (getMapValue(root_map, "dataBlocks") orelse return ParseError.MissingDataBlocks) {
        .array => |array| array,
        else => return ParseError.InvalidMessagePack,
    };
    for (blocks) |block| {
        const block_map = switch (block) {
            .map => |map| map,
            else => continue,
        };
        const categories = switch (getMapValue(block_map, "categories") orelse continue) {
            .array => |array| array,
            else => continue,
        };
        for (categories) |category| {
            const candidate = categoryRowName(category) catch continue;
            if (categoryNameEql(candidate, name)) return category;
        }
    }
    return null;
}

fn categoryRowCount(category: MsgValue) !usize {
    const map = switch (category) {
        .map => |map| map,
        else => return ParseError.InvalidMessagePack,
    };
    const row_count = try msgInt(getMapValue(map, "rowCount") orelse return ParseError.InvalidColumnData);
    return std.math.cast(usize, row_count) orelse ParseError.InvalidColumnData;
}

fn categoryColumns(category: MsgValue) ![]const MsgValue {
    const map = switch (category) {
        .map => |map| map,
        else => return ParseError.InvalidMessagePack,
    };
    return switch (getMapValue(map, "columns") orelse return ParseError.InvalidColumnData) {
        .array => |array| array,
        else => ParseError.InvalidColumnData,
    };
}

fn columnName(column: MsgValue) ![]const u8 {
    const map = switch (column) {
        .map => |map| map,
        else => return ParseError.InvalidColumnData,
    };
    return msgString(getMapValue(map, "name") orelse return ParseError.InvalidColumnData);
}

fn columnData(column: MsgValue) !MsgValue {
    const map = switch (column) {
        .map => |map| map,
        else => return ParseError.InvalidColumnData,
    };
    return getMapValue(map, "data") orelse return ParseError.InvalidColumnData;
}

fn columnMask(column: MsgValue) ?MsgValue {
    const map = switch (column) {
        .map => |map| map,
        else => return null,
    };
    const mask = getMapValue(map, "mask") orelse return null;
    return switch (mask) {
        .nil => null,
        else => mask,
    };
}

fn decodeBcifColumn(allocator: Allocator, column: MsgValue) !DecodedColumn {
    const data_obj = try columnData(column);
    const data_map = switch (data_obj) {
        .map => |map| map,
        else => return ParseError.InvalidColumnData,
    };
    const data = try msgBytes(getMapValue(data_map, "data") orelse return ParseError.InvalidColumnData);
    const encodings = try parseEncodingList(allocator, getMapValue(data_map, "encoding") orelse return ParseError.InvalidColumnData);
    defer freeEncodingList(allocator, encodings);

    if (columnMask(column)) |mask_obj| {
        const mask_map = switch (mask_obj) {
            .map => |map| map,
            else => return ParseError.InvalidColumnData,
        };
        const mask_data = try msgBytes(getMapValue(mask_map, "data") orelse return ParseError.InvalidColumnData);
        const mask_encodings = try parseEncodingList(allocator, getMapValue(mask_map, "encoding") orelse return ParseError.InvalidColumnData);
        defer freeEncodingList(allocator, mask_encodings);
        return decodeColumnWithMask(allocator, data, encodings, mask_data, mask_encodings);
    }

    return decodeColumn(allocator, data, encodings);
}

fn parseInlineCcd(allocator: Allocator, root: MsgValue) !ccd_parser.ComponentDict {
    var dict = ccd_parser.ComponentDict.init(allocator);
    errdefer dict.deinit();

    var builders: std.StringHashMapUnmanaged(BcifCcdBuilder) = .empty;
    defer {
        var it = builders.iterator();
        while (it.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            entry.value_ptr.deinit(allocator);
        }
        builders.deinit(allocator);
    }

    if (try findCategoryOptional(root, "chem_comp_atom")) |atom_category| {
        try parseInlineCcdAtoms(allocator, atom_category, &builders);
    }

    if (builders.count() == 0) return dict;

    if (try findCategoryOptional(root, "chem_comp_bond")) |bond_category| {
        try parseInlineCcdBonds(allocator, bond_category, &builders);
    }

    var it = builders.iterator();
    while (it.next()) |entry| {
        const comp_id = entry.key_ptr.*;
        const builder = entry.value_ptr;

        const atoms = try builder.atoms.toOwnedSlice(allocator);
        var atoms_owned = true;
        errdefer if (atoms_owned) allocator.free(atoms);
        const bonds = try builder.bonds.toOwnedSlice(allocator);
        var bonds_owned = true;
        errdefer if (bonds_owned) allocator.free(bonds);

        var stored = ccd_parser.StoredComponent{
            .comp_id = .{ 0, 0, 0, 0, 0 },
            .comp_id_len = @intCast(@min(comp_id.len, 5)),
            .atoms = atoms,
            .bonds = bonds,
            .allocator = allocator,
        };
        atoms_owned = false;
        bonds_owned = false;
        var stored_owned = true;
        errdefer if (stored_owned) stored.deinit();

        for (comp_id[0..stored.comp_id_len], 0..) |c, idx| {
            stored.comp_id[idx] = c;
        }

        const dict_key = try allocator.dupe(u8, comp_id);
        var key_owned_by_dict = false;
        errdefer if (!key_owned_by_dict) allocator.free(dict_key);
        try dict.owned_keys.append(allocator, dict_key);
        key_owned_by_dict = true;

        try dict.components.put(allocator, dict_key, stored);
        stored_owned = false;
    }

    return dict;
}

fn parseInlineCcdAtoms(
    allocator: Allocator,
    category: MsgValue,
    builders: *std.StringHashMapUnmanaged(BcifCcdBuilder),
) !void {
    const row_count = try categoryRowCount(category);
    const raw_columns = try categoryColumns(category);

    var cols = CcdAtomColumns{};
    for (raw_columns, 0..) |column, idx| {
        mapCcdAtomColumn(&cols, try columnName(column), idx);
    }
    if (!cols.hasRequiredFields()) return;

    var decoded = try allocator.alloc(DecodedColumn, raw_columns.len);
    var decoded_count: usize = 0;
    errdefer {
        for (decoded[0..decoded_count]) |*column| column.deinit(allocator);
        allocator.free(decoded);
    }
    for (raw_columns) |column| {
        decoded[decoded_count] = try decodeBcifColumn(allocator, column);
        decoded_count += 1;
        if (decoded[decoded_count - 1].values.len != row_count) return ParseError.ColumnLengthMismatch;
    }
    defer {
        for (decoded[0..decoded_count]) |*column| column.deinit(allocator);
        allocator.free(decoded);
    }

    var row: usize = 0;
    while (row < row_count) : (row += 1) {
        const comp_id = scalarString(decoded[cols.comp_id.?], row) orelse continue;
        const atom_id = scalarString(decoded[cols.atom_id.?], row) orelse continue;
        const type_symbol = scalarString(decoded[cols.type_symbol.?], row) orelse continue;

        const builder = try getOrCreateCcdBuilder(allocator, builders, comp_id);
        var atom = hyb.CompAtom.init(atom_id, type_symbol);
        if (cols.pdbx_aromatic_flag) |col| atom.aromatic = isYesFlag(scalarString(decoded[col], row));
        if (cols.pdbx_leaving_atom_flag) |col| atom.leaving = isYesFlag(scalarString(decoded[col], row));

        const atom_idx = std.math.cast(u16, builder.atoms.items.len) orelse return ParseError.InvalidColumnData;
        try builder.atoms.append(allocator, atom);

        const name_key = try allocator.dupe(u8, atom_id);
        builder.atom_name_map.put(allocator, name_key, atom_idx) catch |err| {
            allocator.free(name_key);
            return err;
        };
    }
}

fn parseInlineCcdBonds(
    allocator: Allocator,
    category: MsgValue,
    builders: *std.StringHashMapUnmanaged(BcifCcdBuilder),
) !void {
    const row_count = try categoryRowCount(category);
    const raw_columns = try categoryColumns(category);

    var cols = CcdBondColumns{};
    for (raw_columns, 0..) |column, idx| {
        mapCcdBondColumn(&cols, try columnName(column), idx);
    }
    if (!cols.hasRequiredFields()) return;

    var decoded = try allocator.alloc(DecodedColumn, raw_columns.len);
    var decoded_count: usize = 0;
    errdefer {
        for (decoded[0..decoded_count]) |*column| column.deinit(allocator);
        allocator.free(decoded);
    }
    for (raw_columns) |column| {
        decoded[decoded_count] = try decodeBcifColumn(allocator, column);
        decoded_count += 1;
        if (decoded[decoded_count - 1].values.len != row_count) return ParseError.ColumnLengthMismatch;
    }
    defer {
        for (decoded[0..decoded_count]) |*column| column.deinit(allocator);
        allocator.free(decoded);
    }

    var row: usize = 0;
    while (row < row_count) : (row += 1) {
        const comp_id = scalarString(decoded[cols.comp_id.?], row) orelse continue;
        const atom_id_1 = scalarString(decoded[cols.atom_id_1.?], row) orelse continue;
        const atom_id_2 = scalarString(decoded[cols.atom_id_2.?], row) orelse continue;

        const builder = builders.getPtr(comp_id) orelse continue;
        const idx1 = builder.atom_name_map.get(atom_id_1) orelse continue;
        const idx2 = builder.atom_name_map.get(atom_id_2) orelse continue;
        const order = if (cols.value_order) |col|
            if (scalarString(decoded[col], row)) |value| hyb.BondOrder.fromString(value) else hyb.BondOrder.unknown
        else
            hyb.BondOrder.unknown;

        try builder.bonds.append(allocator, .{
            .atom_idx_1 = idx1,
            .atom_idx_2 = idx2,
            .order = order,
            .aromatic = if (cols.pdbx_aromatic_flag) |col| isYesFlag(scalarString(decoded[col], row)) else false,
        });
    }
}

fn getOrCreateCcdBuilder(
    allocator: Allocator,
    builders: *std.StringHashMapUnmanaged(BcifCcdBuilder),
    comp_id: []const u8,
) !*BcifCcdBuilder {
    const gop = try builders.getOrPut(allocator, comp_id);
    if (!gop.found_existing) {
        const key = try allocator.dupe(u8, comp_id);
        gop.key_ptr.* = key;
        gop.value_ptr.* = .{};
    }
    return gop.value_ptr;
}

fn mapCcdAtomColumn(columns: *CcdAtomColumns, name: []const u8, idx: usize) void {
    if (eqlIgnoreCase(name, "comp_id")) {
        columns.comp_id = idx;
    } else if (eqlIgnoreCase(name, "atom_id")) {
        columns.atom_id = idx;
    } else if (eqlIgnoreCase(name, "type_symbol")) {
        columns.type_symbol = idx;
    } else if (eqlIgnoreCase(name, "pdbx_aromatic_flag")) {
        columns.pdbx_aromatic_flag = idx;
    } else if (eqlIgnoreCase(name, "pdbx_leaving_atom_flag")) {
        columns.pdbx_leaving_atom_flag = idx;
    }
}

fn mapCcdBondColumn(columns: *CcdBondColumns, name: []const u8, idx: usize) void {
    if (eqlIgnoreCase(name, "comp_id")) {
        columns.comp_id = idx;
    } else if (eqlIgnoreCase(name, "atom_id_1")) {
        columns.atom_id_1 = idx;
    } else if (eqlIgnoreCase(name, "atom_id_2")) {
        columns.atom_id_2 = idx;
    } else if (eqlIgnoreCase(name, "value_order")) {
        columns.value_order = idx;
    } else if (eqlIgnoreCase(name, "pdbx_aromatic_flag")) {
        columns.pdbx_aromatic_flag = idx;
    }
}

fn isYesFlag(value: ?[]const u8) bool {
    const flag = value orelse return false;
    return flag.len == 1 and (flag[0] == 'Y' or flag[0] == 'y');
}

fn scalarString(column: DecodedColumn, row: usize) ?[]const u8 {
    if (isNull(column, row)) return null;
    if (row >= column.values.len) return null;
    return switch (column.values[row]) {
        .str => |value| if (isCifNullString(value)) null else value,
        else => null,
    };
}

fn scalarFloat(column: DecodedColumn, row: usize) !f64 {
    if (isNull(column, row) or row >= column.values.len) return ParseError.InvalidCoordinate;
    return switch (column.values[row]) {
        .float => |value| value,
        .int => |value| @floatFromInt(value),
        .str => |value| if (isCifNullString(value)) ParseError.InvalidCoordinate else std.fmt.parseFloat(f64, value) catch ParseError.InvalidCoordinate,
    };
}

fn scalarInt(column: DecodedColumn, row: usize) ?i64 {
    if (isNull(column, row) or row >= column.values.len) return null;
    return switch (column.values[row]) {
        .int => |value| value,
        .float => null,
        .str => |value| if (isCifNullString(value)) null else std.fmt.parseInt(i64, value, 10) catch null,
    };
}

fn isNull(column: DecodedColumn, row: usize) bool {
    const nulls = column.nulls orelse return false;
    if (row >= nulls.len) return true;
    return nulls[row] != .present;
}

fn mapAtomSiteColumn(columns: *AtomSiteColumns, name: []const u8, idx: usize) void {
    if (eqlIgnoreCase(name, "Cartn_x")) {
        columns.cartn_x = idx;
    } else if (eqlIgnoreCase(name, "Cartn_y")) {
        columns.cartn_y = idx;
    } else if (eqlIgnoreCase(name, "Cartn_z")) {
        columns.cartn_z = idx;
    } else if (eqlIgnoreCase(name, "type_symbol")) {
        columns.type_symbol = idx;
    } else if (eqlIgnoreCase(name, "label_atom_id")) {
        columns.label_atom_id = idx;
    } else if (eqlIgnoreCase(name, "auth_atom_id")) {
        columns.auth_atom_id = idx;
    } else if (eqlIgnoreCase(name, "label_comp_id")) {
        columns.label_comp_id = idx;
    } else if (eqlIgnoreCase(name, "auth_comp_id")) {
        columns.auth_comp_id = idx;
    } else if (eqlIgnoreCase(name, "label_asym_id")) {
        columns.label_asym_id = idx;
    } else if (eqlIgnoreCase(name, "auth_asym_id")) {
        columns.auth_asym_id = idx;
    } else if (eqlIgnoreCase(name, "label_seq_id")) {
        columns.label_seq_id = idx;
    } else if (eqlIgnoreCase(name, "auth_seq_id")) {
        columns.auth_seq_id = idx;
    } else if (eqlIgnoreCase(name, "pdbx_PDB_ins_code")) {
        columns.pdbx_pdb_ins_code = idx;
    } else if (eqlIgnoreCase(name, "group_PDB")) {
        columns.group_pdb = idx;
    } else if (eqlIgnoreCase(name, "label_alt_id")) {
        columns.label_alt_id = idx;
    } else if (eqlIgnoreCase(name, "occupancy")) {
        columns.occupancy = idx;
    } else if (eqlIgnoreCase(name, "pdbx_PDB_model_num")) {
        columns.pdbx_pdb_model_num = idx;
    }
}

fn categoryRowName(category: MsgValue) ![]const u8 {
    const map = switch (category) {
        .map => |map| map,
        else => return ParseError.InvalidMessagePack,
    };
    return msgString(getMapValue(map, "name") orelse return ParseError.InvalidColumnData);
}

fn categoryNameEql(a: []const u8, b: []const u8) bool {
    return eqlIgnoreCase(stripLeadingUnderscore(a), stripLeadingUnderscore(b));
}

fn stripLeadingUnderscore(value: []const u8) []const u8 {
    if (value.len > 0 and value[0] == '_') return value[1..];
    return value;
}

fn msgString(value: MsgValue) ![]const u8 {
    return switch (value) {
        .str => |str| str,
        else => ParseError.InvalidColumnData,
    };
}

fn msgBytes(value: MsgValue) ![]const u8 {
    return switch (value) {
        .bin => |bin| bin,
        .str => |str| str,
        else => ParseError.InvalidColumnData,
    };
}

fn msgInt(value: MsgValue) !i64 {
    return switch (value) {
        .int => |int| int,
        .uint => |uint| std.math.cast(i64, uint) orelse ParseError.InvalidColumnData,
        else => ParseError.InvalidColumnData,
    };
}

fn msgFloat(value: MsgValue) !f64 {
    return switch (value) {
        .float => |float| float,
        .int => |int| @floatFromInt(int),
        .uint => |uint| @floatFromInt(uint),
        else => ParseError.InvalidColumnData,
    };
}

fn msgBool(value: MsgValue) !bool {
    return switch (value) {
        .bool => |boolean| boolean,
        else => ParseError.InvalidColumnData,
    };
}

fn isCifNullString(value: []const u8) bool {
    return value.len == 1 and (value[0] == '.' or value[0] == '?');
}

fn eqlIgnoreCase(a: []const u8, b: []const u8) bool {
    if (a.len != b.len) return false;
    for (a, b) |ca, cb| {
        if (std.ascii.toLower(ca) != std.ascii.toLower(cb)) return false;
    }
    return true;
}

const BuildBcifOptions = struct {
    include_hydrogen: bool = false,
    include_hetatm: bool = false,
    include_second_model: bool = false,
    omit_z: bool = false,
    omit_model: bool = false,
    omit_chain: bool = false,
    null_first_model: bool = false,
    include_long_chain: bool = false,
    include_very_long_chain: bool = false,
    include_alt_later_b_only: bool = false,
    include_alt_across_models: bool = false,
};

fn buildMinimalBcif(options: BuildBcifOptions) ![]u8 {
    const allocator = std.testing.allocator;
    var rows = std.ArrayListUnmanaged(TestAtomRow).empty;
    defer rows.deinit(allocator);

    const first_chain = if (options.include_very_long_chain)
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg"
    else if (options.include_long_chain)
        "ABCDE"
    else
        "A";
    const second_chain = if (options.include_very_long_chain)
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef"
    else if (options.include_long_chain)
        "ABCD"
    else
        "A";
    if (options.include_alt_across_models) {
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = "A", .seq = 1, .alt = "A", .x = 10.0, .y = 20.0, .z = 30.0, .model = 1 });
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = "A", .seq = 1, .alt = "B", .x = 12.0, .y = 22.0, .z = 32.0, .model = 1 });
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = "A", .seq = 1, .alt = "B", .x = 14.0, .y = 24.0, .z = 34.0, .model = 2 });
    } else if (options.include_alt_later_b_only) {
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = "A", .seq = 1, .alt = "A", .x = 10.0, .y = 20.0, .z = 30.0, .model = 1 });
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = "A", .seq = 1, .alt = "B", .x = 12.0, .y = 22.0, .z = 32.0, .model = 1 });
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "GLY", .chain = "A", .seq = 2, .alt = "B", .x = 14.0, .y = 24.0, .z = 34.0, .model = 1 });
    } else {
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = first_chain, .seq = 1, .x = 10.0, .y = 20.0, .z = 30.0, .model = 1 });
        try rows.append(allocator, .{ .group = "ATOM", .element = "N", .atom = "N", .residue = "ALA", .chain = second_chain, .seq = 1, .x = 11.0, .y = 21.0, .z = 32.0, .model = 1 });
    }
    if (options.include_hydrogen) {
        try rows.append(allocator, .{ .group = "ATOM", .element = "H", .atom = "H", .residue = "ALA", .chain = "A", .seq = 1, .x = 12.0, .y = 22.0, .z = 33.0, .model = 1 });
    }
    if (options.include_hetatm) {
        try rows.append(allocator, .{ .group = "HETATM", .element = "O", .atom = "O", .residue = "HOH", .chain = "A", .seq = 2, .x = 13.0, .y = 23.0, .z = 34.0, .model = 1 });
    }
    if (options.include_second_model) {
        try rows.append(allocator, .{ .group = "ATOM", .element = "C", .atom = "CA", .residue = "ALA", .chain = "B", .seq = 1, .x = 14.0, .y = 24.0, .z = 35.0, .model = 2 });
    }

    var bytes = std.ArrayListUnmanaged(u8).empty;
    errdefer bytes.deinit(allocator);

    try packMapHeader(allocator, &bytes, 3);
    try packStr(allocator, &bytes, "version");
    try packStr(allocator, &bytes, "0.3.0");
    try packStr(allocator, &bytes, "encoder");
    try packStr(allocator, &bytes, "zsasa-test");
    try packStr(allocator, &bytes, "dataBlocks");
    try packArrayHeader(allocator, &bytes, 1);
    try packMapHeader(allocator, &bytes, 2);
    try packStr(allocator, &bytes, "header");
    try packStr(allocator, &bytes, "TEST");
    try packStr(allocator, &bytes, "categories");
    try packArrayHeader(allocator, &bytes, 1);
    try packMapHeader(allocator, &bytes, 3);
    try packStr(allocator, &bytes, "name");
    try packStr(allocator, &bytes, "_atom_site");
    try packStr(allocator, &bytes, "rowCount");
    try packInt(allocator, &bytes, @intCast(rows.items.len));
    try packStr(allocator, &bytes, "columns");

    var column_count: usize = 9;
    if (!options.omit_z) column_count += 1;
    if (!options.omit_model) column_count += 1;
    if (!options.omit_chain) column_count += 2;
    try packArrayHeader(allocator, &bytes, column_count);
    try packStringColumn(allocator, &bytes, "group_PDB", rows.items, TestAtomRow.groupValue);
    try packStringColumn(allocator, &bytes, "type_symbol", rows.items, TestAtomRow.elementValue);
    try packStringColumn(allocator, &bytes, "label_atom_id", rows.items, TestAtomRow.atomValue);
    try packStringColumn(allocator, &bytes, "label_comp_id", rows.items, TestAtomRow.residueValue);
    if (!options.omit_chain) try packStringColumn(allocator, &bytes, "label_asym_id", rows.items, TestAtomRow.chainValue);
    try packIntColumn(allocator, &bytes, "label_seq_id", rows.items, TestAtomRow.seqValue, false);
    try packStringColumn(allocator, &bytes, "pdbx_PDB_ins_code", rows.items, TestAtomRow.emptyValue);
    try packStringColumn(allocator, &bytes, "label_alt_id", rows.items, TestAtomRow.altValue);
    if (!options.omit_model) try packIntColumn(allocator, &bytes, "pdbx_PDB_model_num", rows.items, TestAtomRow.modelValue, options.null_first_model);
    try packFloatColumn(allocator, &bytes, "Cartn_x", rows.items, TestAtomRow.xValue);
    try packFloatColumn(allocator, &bytes, "Cartn_y", rows.items, TestAtomRow.yValue);
    if (!options.omit_z) try packFloatColumn(allocator, &bytes, "Cartn_z", rows.items, TestAtomRow.zValue);
    if (!options.omit_chain) try packStringColumn(allocator, &bytes, "auth_asym_id", rows.items, TestAtomRow.chainValue);

    return bytes.toOwnedSlice(allocator);
}

fn buildMinimalBcifWithInlineCcd() ![]u8 {
    const allocator = std.testing.allocator;

    const atom_comp_ids = [_][]const u8{ "LIG", "LIG", "LIG" };
    const atom_ids = [_][]const u8{ "C1", "O1", "H1" };
    const atom_types = [_][]const u8{ "C", "O", "H" };
    const atom_aromatic = [_][]const u8{ "N", "N", "N" };
    const atom_leaving = [_][]const u8{ "N", "N", "N" };

    const bond_comp_ids = [_][]const u8{ "LIG", "LIG" };
    const bond_atom_1 = [_][]const u8{ "C1", "O1" };
    const bond_atom_2 = [_][]const u8{ "O1", "H1" };
    const bond_orders = [_][]const u8{ "DOUB", "SING" };
    const bond_aromatic = [_][]const u8{ "N", "N" };

    var bytes = std.ArrayListUnmanaged(u8).empty;
    errdefer bytes.deinit(allocator);

    try packMapHeader(allocator, &bytes, 3);
    try packStr(allocator, &bytes, "version");
    try packStr(allocator, &bytes, "0.3.0");
    try packStr(allocator, &bytes, "encoder");
    try packStr(allocator, &bytes, "zsasa-test");
    try packStr(allocator, &bytes, "dataBlocks");
    try packArrayHeader(allocator, &bytes, 1);
    try packMapHeader(allocator, &bytes, 2);
    try packStr(allocator, &bytes, "header");
    try packStr(allocator, &bytes, "TEST");
    try packStr(allocator, &bytes, "categories");
    try packArrayHeader(allocator, &bytes, 3);

    try packStringOnlyCategory(allocator, &bytes, "_chem_comp_atom", atom_comp_ids.len, &.{
        .{ .name = "comp_id", .values = atom_comp_ids[0..] },
        .{ .name = "atom_id", .values = atom_ids[0..] },
        .{ .name = "type_symbol", .values = atom_types[0..] },
        .{ .name = "pdbx_aromatic_flag", .values = atom_aromatic[0..] },
        .{ .name = "pdbx_leaving_atom_flag", .values = atom_leaving[0..] },
    });
    try packStringOnlyCategory(allocator, &bytes, "_chem_comp_bond", bond_comp_ids.len, &.{
        .{ .name = "comp_id", .values = bond_comp_ids[0..] },
        .{ .name = "atom_id_1", .values = bond_atom_1[0..] },
        .{ .name = "atom_id_2", .values = bond_atom_2[0..] },
        .{ .name = "value_order", .values = bond_orders[0..] },
        .{ .name = "pdbx_aromatic_flag", .values = bond_aromatic[0..] },
    });

    const atom_site = try buildMinimalBcif(.{ .include_hetatm = true });
    defer allocator.free(atom_site);
    try appendAtomSiteCategoryFromMinimalBcif(allocator, &bytes, atom_site);

    return bytes.toOwnedSlice(allocator);
}

const TestStringColumn = struct {
    name: []const u8,
    values: []const []const u8,
};

fn packStringOnlyCategory(
    allocator: Allocator,
    bytes: *std.ArrayListUnmanaged(u8),
    category_name: []const u8,
    row_count: usize,
    columns: []const TestStringColumn,
) !void {
    try packMapHeader(allocator, bytes, 3);
    try packStr(allocator, bytes, "name");
    try packStr(allocator, bytes, category_name);
    try packStr(allocator, bytes, "rowCount");
    try packInt(allocator, bytes, @intCast(row_count));
    try packStr(allocator, bytes, "columns");
    try packArrayHeader(allocator, bytes, columns.len);
    for (columns) |column| {
        try packStringValuesColumn(allocator, bytes, column.name, column.values);
    }
}

fn packStringValuesColumn(
    allocator: Allocator,
    bytes: *std.ArrayListUnmanaged(u8),
    name: []const u8,
    values: []const []const u8,
) !void {
    var string_data = std.ArrayListUnmanaged(u8).empty;
    defer string_data.deinit(allocator);
    var offsets = std.ArrayListUnmanaged(u8).empty;
    defer offsets.deinit(allocator);
    var indices = std.ArrayListUnmanaged(u8).empty;
    defer indices.deinit(allocator);

    try offsets.append(allocator, 0);
    for (values, 0..) |value, idx| {
        try string_data.appendSlice(allocator, value);
        try offsets.append(allocator, @intCast(string_data.items.len));
        try indices.append(allocator, @intCast(idx));
    }

    try packColumnHeader(allocator, bytes, name);
    try packEncodedDataMapHeader(allocator, bytes);
    try packBin(allocator, bytes, indices.items);
    try packStr(allocator, bytes, "encoding");
    try packStringArrayEncoding(allocator, bytes, string_data.items, offsets.items);
}

fn appendAtomSiteCategoryFromMinimalBcif(
    allocator: Allocator,
    bytes: *std.ArrayListUnmanaged(u8),
    source: []const u8,
) !void {
    var reader = MsgReader.init(allocator, source);
    const root = try reader.readValue();
    defer freeMsgValue(allocator, root);
    const atom_site = try findCategory(root, "atom_site");
    try packMsgValue(allocator, bytes, atom_site);
}

fn packMsgValue(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), value: MsgValue) !void {
    switch (value) {
        .nil => try bytes.append(allocator, 0xc0),
        .bool => |b| try bytes.append(allocator, if (b) 0xc3 else 0xc2),
        .int => |i| try packInt(allocator, bytes, i),
        .uint => |u| try packInt(allocator, bytes, std.math.cast(i64, u) orelse return error.Overflow),
        .float => |f| {
            try bytes.append(allocator, 0xcb);
            var buf: [8]u8 = undefined;
            std.mem.writeInt(u64, &buf, @bitCast(f), .big);
            try bytes.appendSlice(allocator, &buf);
        },
        .str => |s| try packStr(allocator, bytes, s),
        .bin => |b| try packBin(allocator, bytes, b),
        .array => |items| {
            try packArrayHeader(allocator, bytes, items.len);
            for (items) |item| try packMsgValue(allocator, bytes, item);
        },
        .map => |pairs| {
            try packMapHeader(allocator, bytes, pairs.len);
            for (pairs) |pair| {
                try packMsgValue(allocator, bytes, pair.key);
                try packMsgValue(allocator, bytes, pair.value);
            }
        },
    }
}

const TestAtomRow = struct {
    group: []const u8,
    element: []const u8,
    atom: []const u8,
    residue: []const u8,
    chain: []const u8,
    seq: i32,
    x: f32,
    y: f32,
    z: f32,
    model: i32,
    alt: []const u8 = "",

    fn groupValue(row: TestAtomRow) []const u8 {
        return row.group;
    }
    fn elementValue(row: TestAtomRow) []const u8 {
        return row.element;
    }
    fn atomValue(row: TestAtomRow) []const u8 {
        return row.atom;
    }
    fn residueValue(row: TestAtomRow) []const u8 {
        return row.residue;
    }
    fn chainValue(row: TestAtomRow) []const u8 {
        return row.chain;
    }
    fn emptyValue(_: TestAtomRow) []const u8 {
        return "";
    }
    fn altValue(row: TestAtomRow) []const u8 {
        return row.alt;
    }
    fn seqValue(row: TestAtomRow) i32 {
        return row.seq;
    }
    fn modelValue(row: TestAtomRow) i32 {
        return row.model;
    }
    fn xValue(row: TestAtomRow) f32 {
        return row.x;
    }
    fn yValue(row: TestAtomRow) f32 {
        return row.y;
    }
    fn zValue(row: TestAtomRow) f32 {
        return row.z;
    }
};

fn packStringColumn(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), name: []const u8, rows: []const TestAtomRow, comptime get: fn (TestAtomRow) []const u8) !void {
    var string_data = std.ArrayListUnmanaged(u8).empty;
    defer string_data.deinit(allocator);
    var offsets = std.ArrayListUnmanaged(u8).empty;
    defer offsets.deinit(allocator);
    var indices = std.ArrayListUnmanaged(u8).empty;
    defer indices.deinit(allocator);

    try offsets.append(allocator, 0);
    for (rows, 0..) |row, idx| {
        const value = get(row);
        try string_data.appendSlice(allocator, value);
        try offsets.append(allocator, @intCast(string_data.items.len));
        try indices.append(allocator, @intCast(idx));
    }

    try packColumnHeader(allocator, bytes, name);
    try packEncodedDataMapHeader(allocator, bytes);
    try packBin(allocator, bytes, indices.items);
    try packStr(allocator, bytes, "encoding");
    try packStringArrayEncoding(allocator, bytes, string_data.items, offsets.items);
}

fn packIntColumn(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), name: []const u8, rows: []const TestAtomRow, comptime get: fn (TestAtomRow) i32, null_first: bool) !void {
    var data = std.ArrayListUnmanaged(u8).empty;
    defer data.deinit(allocator);
    for (rows) |row| {
        const value = get(row);
        var buf: [4]u8 = undefined;
        std.mem.writeInt(i32, &buf, value, .little);
        try data.appendSlice(allocator, &buf);
    }
    try packColumnHeaderWithFieldCount(allocator, bytes, name, if (null_first) 3 else 2);
    try packEncodedDataMapHeader(allocator, bytes);
    try packBin(allocator, bytes, data.items);
    try packStr(allocator, bytes, "encoding");
    try packArrayHeader(allocator, bytes, 1);
    try packByteArrayEncoding(allocator, bytes, 3);
    if (null_first) try packNullFirstMask(allocator, bytes, rows.len);
}

fn packFloatColumn(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), name: []const u8, rows: []const TestAtomRow, comptime get: fn (TestAtomRow) f32) !void {
    var data = std.ArrayListUnmanaged(u8).empty;
    defer data.deinit(allocator);
    for (rows) |row| {
        const bits: u32 = @bitCast(get(row));
        var buf: [4]u8 = undefined;
        std.mem.writeInt(u32, &buf, bits, .little);
        try data.appendSlice(allocator, &buf);
    }
    try packColumnHeader(allocator, bytes, name);
    try packEncodedDataMapHeader(allocator, bytes);
    try packBin(allocator, bytes, data.items);
    try packStr(allocator, bytes, "encoding");
    try packArrayHeader(allocator, bytes, 1);
    try packByteArrayEncoding(allocator, bytes, 32);
}

fn packColumnHeader(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), name: []const u8) !void {
    try packColumnHeaderWithFieldCount(allocator, bytes, name, 2);
}

fn packColumnHeaderWithFieldCount(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), name: []const u8, field_count: usize) !void {
    try packMapHeader(allocator, bytes, field_count);
    try packStr(allocator, bytes, "name");
    try packStr(allocator, bytes, name);
    try packStr(allocator, bytes, "data");
}

fn packNullFirstMask(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), row_count: usize) !void {
    try packStr(allocator, bytes, "mask");
    try packEncodedDataMapHeader(allocator, bytes);
    var mask = std.ArrayListUnmanaged(u8).empty;
    defer mask.deinit(allocator);
    var row: usize = 0;
    while (row < row_count) : (row += 1) {
        try mask.append(allocator, if (row == 0) 1 else 0);
    }
    try packBin(allocator, bytes, mask.items);
    try packStr(allocator, bytes, "encoding");
    try packArrayHeader(allocator, bytes, 1);
    try packByteArrayEncoding(allocator, bytes, 4);
}

fn packEncodedDataMapHeader(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8)) !void {
    try packMapHeader(allocator, bytes, 2);
    try packStr(allocator, bytes, "data");
}

fn packByteArrayEncoding(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), type_code: u8) !void {
    try packMapHeader(allocator, bytes, 2);
    try packStr(allocator, bytes, "kind");
    try packStr(allocator, bytes, "ByteArray");
    try packStr(allocator, bytes, "type");
    try packInt(allocator, bytes, type_code);
}

fn packStringArrayEncoding(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), string_data: []const u8, offsets: []const u8) !void {
    try packArrayHeader(allocator, bytes, 1);
    try packMapHeader(allocator, bytes, 5);
    try packStr(allocator, bytes, "kind");
    try packStr(allocator, bytes, "StringArray");
    try packStr(allocator, bytes, "stringData");
    try packStr(allocator, bytes, string_data);
    try packStr(allocator, bytes, "offsets");
    try packBin(allocator, bytes, offsets);
    try packStr(allocator, bytes, "offsetEncoding");
    try packArrayHeader(allocator, bytes, 1);
    try packByteArrayEncoding(allocator, bytes, 4);
    try packStr(allocator, bytes, "dataEncoding");
    try packArrayHeader(allocator, bytes, 1);
    try packByteArrayEncoding(allocator, bytes, 4);
}

fn packMapHeader(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), len: usize) !void {
    if (len <= 15) return bytes.append(allocator, 0x80 | @as(u8, @intCast(len)));
    try bytes.append(allocator, 0xde);
    var buf: [2]u8 = undefined;
    std.mem.writeInt(u16, &buf, @intCast(len), .big);
    try bytes.appendSlice(allocator, &buf);
}

fn packArrayHeader(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), len: usize) !void {
    if (len <= 15) return bytes.append(allocator, 0x90 | @as(u8, @intCast(len)));
    try bytes.append(allocator, 0xdc);
    var buf: [2]u8 = undefined;
    std.mem.writeInt(u16, &buf, @intCast(len), .big);
    try bytes.appendSlice(allocator, &buf);
}

fn packStr(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), value: []const u8) !void {
    if (value.len <= 31) {
        try bytes.append(allocator, 0xa0 | @as(u8, @intCast(value.len)));
    } else {
        try bytes.append(allocator, 0xd9);
        try bytes.append(allocator, @intCast(value.len));
    }
    try bytes.appendSlice(allocator, value);
}

fn packBin(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), value: []const u8) !void {
    try bytes.append(allocator, 0xc4);
    try bytes.append(allocator, @intCast(value.len));
    try bytes.appendSlice(allocator, value);
}

fn packInt(allocator: Allocator, bytes: *std.ArrayListUnmanaged(u8), value: i64) !void {
    if (value >= 0 and value <= 127) {
        try bytes.append(allocator, @intCast(value));
    } else {
        try bytes.append(allocator, 0xd2);
        var buf: [4]u8 = undefined;
        std.mem.writeInt(i32, &buf, @intCast(value), .big);
        try bytes.appendSlice(allocator, &buf);
    }
}

test "parse minimal BinaryCIF atom_site" {
    const source = try buildMinimalBcif(.{});
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 21.0), input.y[1], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 32.0), input.z[1], 0.001);
    try std.testing.expectEqual(@as(u8, 6), input.element.?[0]);
    try std.testing.expectEqualStrings("ALA", input.residue.?[0].slice());
    try std.testing.expectEqualStrings("CA", input.atom_name.?[0].slice());
    try std.testing.expectEqualStrings("A", input.chain_id.?[0].slice());
    try std.testing.expectEqual(@as(i32, 1), input.residue_num.?[0]);
}

test "parse BinaryCIF keeps extended chain IDs lossless beyond fixed buffers" {
    const source = try buildMinimalBcif(.{ .include_very_long_chain = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expect(input.chain_id_full != null);
    try std.testing.expectEqualStrings("ABCD", input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg", input.chain_id_full.?[0]);
    try std.testing.expectEqualStrings("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef", input.chain_id_full.?[1]);
}

test "parse BinaryCIF with inline CCD data" {
    const source = try buildMinimalBcifWithInlineCcd();
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.atom_only = false;
    var input = try parser.parse(source);
    defer input.deinit();
    defer parser.deinitCcd();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());

    const ccd = parser.getInlineCcd() orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 1), ccd.count());

    const comp = ccd.get("LIG") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("LIG", comp.compIdSlice());
    try std.testing.expectEqual(@as(usize, 3), comp.atoms.len);
    try std.testing.expectEqual(@as(usize, 2), comp.bonds.len);

    try std.testing.expectEqualStrings("C1", comp.atoms[0].atomIdSlice());
    try std.testing.expectEqualStrings("C", comp.atoms[0].typeSymbolSlice());
    try std.testing.expectEqualStrings("O1", comp.atoms[1].atomIdSlice());
    try std.testing.expectEqualStrings("O", comp.atoms[1].typeSymbolSlice());

    try std.testing.expectEqual(@as(u16, 0), comp.bonds[0].atom_idx_1);
    try std.testing.expectEqual(@as(u16, 1), comp.bonds[0].atom_idx_2);
    try std.testing.expectEqual(.double, comp.bonds[0].order);
}

test "parse BinaryCIF default model selection includes all models" {
    const source = try buildMinimalBcif(.{ .include_second_model = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
    try std.testing.expectEqualStrings("B", input.chain_id.?[2].slice());
}

test "parse BinaryCIF explicit model selection filters requested model" {
    const source = try buildMinimalBcif(.{ .include_second_model = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.model_num = 2;
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectEqualStrings("B", input.chain_id.?[0].slice());
}

test "parse BinaryCIF altLoc selection is per atom site and keeps later B-only sites" {
    const source = try buildMinimalBcif(.{ .include_alt_later_b_only = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 14.0), input.x[1], 0.001);
    try std.testing.expectEqual(@as(i32, 2), input.residue_num.?[1]);
}

test "parse BinaryCIF altLoc selection is scoped by model" {
    const source = try buildMinimalBcif(.{ .include_alt_across_models = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 14.0), input.x[1], 0.001);
}

test "parse BinaryCIF applies filters" {
    const source = try buildMinimalBcif(.{ .include_hydrogen = true, .include_hetatm = true, .include_second_model = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.model_num = 1;
    parser.chain_filter = &[_][]const u8{"A"};
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    for (input.chain_id.?) |chain| {
        try std.testing.expectEqualStrings("A", chain.slice());
    }
}

test "parse BinaryCIF reports missing coordinate field" {
    const source = try buildMinimalBcif(.{ .omit_z = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    try std.testing.expectError(ParseError.MissingCoordinateField, parser.parse(source));
}

test "parse BinaryCIF model filter ignores absent model column" {
    const source = try buildMinimalBcif(.{ .include_second_model = true, .omit_model = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.model_num = 2;
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
}

test "parse BinaryCIF model filter ignores null model value" {
    const source = try buildMinimalBcif(.{ .null_first_model = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.model_num = 2;
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectEqualStrings("CA", input.atom_name.?[0].slice());
}

test "parse BinaryCIF chain filter ignores absent chain column" {
    const source = try buildMinimalBcif(.{ .omit_chain = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.chain_filter = &[_][]const u8{"A"};
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectEqualStrings("", input.chain_id.?[0].slice());
}

fn expectMinimalBcifFileParses(path: []const u8) !void {
    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parseFile(std.testing.io, path);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 21.0), input.y[1], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 32.0), input.z[1], 0.001);
}

test "parseFile reads plain BinaryCIF fixture" {
    try expectMinimalBcifFileParses("test_data/bcif/minimal.bcif");
}

test "parseFile reads gzip BinaryCIF fixture" {
    try expectMinimalBcifFileParses("test_data/bcif/minimal.bcif.gz");
}

test "parseFile reads zstd BinaryCIF fixture" {
    try expectMinimalBcifFileParses("test_data/bcif/minimal.bcif.zst");
}

test "bcif scalarInt rejects out-of-range float integer field" {
    var values = [_]Scalar{.{ .float = 1.0e100 }};
    const column = DecodedColumn{ .values = &values };
    try std.testing.expectEqual(@as(?i64, null), scalarInt(column, 0));
}

test "msgpack reader decodes primitives" {
    const bytes = [_]u8{
        0x88,
        0xa3,
        'n',
        'i',
        'l',
        0xc0,
        0xa4,
        't',
        'r',
        'u',
        'e',
        0xc3,
        0xa5,
        'f',
        'a',
        'l',
        's',
        'e',
        0xc2,
        0xa3,
        'i',
        'n',
        't',
        0xd0,
        0xfe,
        0xa4,
        'u',
        'i',
        'n',
        't',
        0xcc,
        200,
        0xa3,
        'f',
        '3',
        '2',
        0xca,
        0x3f,
        0x80,
        0x00,
        0x00,
        0xa3,
        's',
        't',
        'r',
        0xa3,
        'a',
        'b',
        'c',
        0xa3,
        'b',
        'i',
        'n',
        0xc4,
        0x03,
        1,
        2,
        3,
    };
    var reader = MsgReader.init(std.testing.allocator, &bytes);
    const value = try reader.readValue();
    defer freeMsgValue(std.testing.allocator, value);

    const map = value.map;
    try std.testing.expectEqual(@as(usize, 8), map.len);
    try std.testing.expectEqual(MsgValue.nil, getMapValue(map, "nil").?);
    try std.testing.expectEqual(true, getMapValue(map, "true").?.bool);
    try std.testing.expectEqual(false, getMapValue(map, "false").?.bool);
    try std.testing.expectEqual(@as(i64, -2), getMapValue(map, "int").?.int);
    try std.testing.expectEqual(@as(u64, 200), getMapValue(map, "uint").?.uint);
    try std.testing.expectApproxEqAbs(@as(f64, 1.0), getMapValue(map, "f32").?.float, 0.000001);
    try std.testing.expectEqualStrings("abc", getMapValue(map, "str").?.str);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 1, 2, 3 }, getMapValue(map, "bin").?.bin);
}

test "msgpack reader decodes arrays and 16 bit lengths" {
    var bytes = std.ArrayListUnmanaged(u8).empty;
    defer bytes.deinit(std.testing.allocator);
    try bytes.appendSlice(std.testing.allocator, &[_]u8{ 0xdc, 0x00, 0x02 });
    try bytes.appendSlice(std.testing.allocator, &[_]u8{ 0xda, 0x00, 0x05, 'h', 'e', 'l', 'l', 'o' });
    try bytes.appendSlice(std.testing.allocator, &[_]u8{ 0xc5, 0x00, 0x04, 9, 8, 7, 6 });

    var reader = MsgReader.init(std.testing.allocator, bytes.items);
    const value = try reader.readValue();
    defer freeMsgValue(std.testing.allocator, value);

    try std.testing.expectEqual(@as(usize, 2), value.array.len);
    try std.testing.expectEqualStrings("hello", value.array[0].str);
    try std.testing.expectEqualSlices(u8, &[_]u8{ 9, 8, 7, 6 }, value.array[1].bin);
}

test "bcif decoders decode numeric encoding chains" {
    const allocator = std.testing.allocator;
    const byte_array = [_]u8{ 120, 0, 123, 0, 12, 0 };
    const encs = [_]Encoding{
        .{ .fixed_point = .{ .factor = 100.0 } },
        .{ .byte_array = .{ .type_code = 2 } },
    };
    var column = try decodeColumn(allocator, &byte_array, &encs);
    defer column.deinit(allocator);

    try std.testing.expectEqual(@as(usize, 3), column.values.len);
    try std.testing.expectApproxEqAbs(@as(f64, 1.20), column.values[0].float, 0.000001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.23), column.values[1].float, 0.000001);
    try std.testing.expectApproxEqAbs(@as(f64, 0.12), column.values[2].float, 0.000001);
}

test "bcif decoders decode integer packing run length delta" {
    const allocator = std.testing.allocator;
    const packed_bytes = [_]u8{ 1, 3, 2, 1 };
    const encs = [_]Encoding{
        .{ .delta = .{ .origin = 10 } },
        .{ .run_length = .{ .src_size = 4 } },
        .{ .integer_packing = .{ .byte_count = 1, .is_unsigned = false, .src_size = 4 } },
        .{ .byte_array = .{ .type_code = 1 } },
    };
    var column = try decodeColumn(allocator, &packed_bytes, &encs);
    defer column.deinit(allocator);

    try std.testing.expectEqual(@as(i64, 11), column.values[0].int);
    try std.testing.expectEqual(@as(i64, 12), column.values[1].int);
    try std.testing.expectEqual(@as(i64, 13), column.values[2].int);
    try std.testing.expectEqual(@as(i64, 15), column.values[3].int);
}

test "bcif decoders decode string arrays" {
    const allocator = std.testing.allocator;
    const index_bytes = [_]u8{ 0, 1, 0, 2 };
    const offset_bytes = [_]u8{ 0, 1, 4, 7 };
    const encoding = Encoding{ .string_array = .{
        .string_data = "ACAGLY",
        .offset_data = &offset_bytes,
        .offset_encoding = &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
        .data_encoding = &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
    } };
    var column = try decodeColumn(allocator, &index_bytes, &[_]Encoding{encoding});
    defer column.deinit(allocator);

    try std.testing.expectEqualStrings("A", column.values[0].str);
    try std.testing.expectEqualStrings("CAG", column.values[1].str);
    try std.testing.expectEqualStrings("A", column.values[2].str);
    try std.testing.expectEqualStrings("LY", column.values[3].str);
}

test "bcif mask maps dot and question to null states" {
    const allocator = std.testing.allocator;
    const data = [_]u8{ 1, 2, 3 };
    const mask_data = [_]u8{ 0, 1, 2 };
    var column = try decodeColumnWithMask(
        allocator,
        &data,
        &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
        &mask_data,
        &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
    );
    defer column.deinit(allocator);

    try std.testing.expect(column.values[0] == .int);
    try std.testing.expectEqual(NullKind.not_present, column.nulls.?[1]);
    try std.testing.expectEqual(NullKind.unknown, column.nulls.?[2]);
}

test "bcif column treats nil mask as absent mask" {
    const allocator = std.testing.allocator;
    const data = [_]u8{ 42, 7 };

    var encoding_pairs = [_]MsgPair{
        .{ .key = .{ .str = "kind" }, .value = .{ .str = "ByteArray" } },
        .{ .key = .{ .str = "type" }, .value = .{ .int = 4 } },
    };
    var encodings = [_]MsgValue{.{ .map = &encoding_pairs }};
    var data_pairs = [_]MsgPair{
        .{ .key = .{ .str = "data" }, .value = .{ .bin = &data } },
        .{ .key = .{ .str = "encoding" }, .value = .{ .array = &encodings } },
    };
    var column_pairs = [_]MsgPair{
        .{ .key = .{ .str = "name" }, .value = .{ .str = "test" } },
        .{ .key = .{ .str = "data" }, .value = .{ .map = &data_pairs } },
        .{ .key = .{ .str = "mask" }, .value = .nil },
    };

    var column = try decodeBcifColumn(allocator, .{ .map = &column_pairs });
    defer column.deinit(allocator);

    try std.testing.expect(column.nulls == null);
    try std.testing.expectEqual(@as(i64, 42), column.values[0].int);
    try std.testing.expectEqual(@as(i64, 7), column.values[1].int);
}

test "bcif masked string arrays tolerate negative null sentinels" {
    const allocator = std.testing.allocator;
    const index_bytes = [_]u8{ 0xff, 1, 0xff };
    const offset_bytes = [_]u8{ 0, 3, 6 };
    const mask_data = [_]u8{ 1, 0, 2 };
    const encoding = Encoding{ .string_array = .{
        .string_data = "ABCDEF",
        .offset_data = &offset_bytes,
        .offset_encoding = &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
        .data_encoding = &[_]Encoding{.{ .byte_array = .{ .type_code = 1 } }},
    } };

    var column = try decodeColumnWithMask(
        allocator,
        &index_bytes,
        &[_]Encoding{encoding},
        &mask_data,
        &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
    );
    defer column.deinit(allocator);

    try std.testing.expectEqual(NullKind.not_present, column.nulls.?[0]);
    try std.testing.expectEqualStrings("DEF", column.values[1].str);
    try std.testing.expectEqual(NullKind.unknown, column.nulls.?[2]);
}

test "bcif decoder chain frees intermediate column once on transform error" {
    const allocator = std.testing.allocator;
    const bytes = [_]u8{ 1, 0 };
    try std.testing.expectError(
        ParseError.InvalidColumnData,
        decodeColumn(allocator, &bytes, &[_]Encoding{
            .{ .fixed_point = .{ .factor = 0.0 } },
            .{ .byte_array = .{ .type_code = 2 } },
        }),
    );
}

test "bcif integer packing accumulation rejects overflow" {
    const allocator = std.testing.allocator;
    var out = std.ArrayListUnmanaged(Scalar).empty;
    defer out.deinit(allocator);
    var acc: i64 = std.math.maxInt(i64);

    try std.testing.expectError(
        ParseError.InvalidColumnData,
        appendIntegerPackingPart(
            allocator,
            &out,
            &acc,
            1,
            .{ .positive = 127, .negative = -128 },
        ),
    );
}

test "bcif delta rejects accumulation overflow" {
    const allocator = std.testing.allocator;
    const values = try allocator.dupe(Scalar, &[_]Scalar{.{ .int = 1 }});
    const column = DecodedColumn{ .values = values };

    try std.testing.expectError(
        ParseError.InvalidColumnData,
        decodeDelta(allocator, column, std.math.maxInt(i64)),
    );
}

test "bcif run length rejects malformed pairs" {
    const allocator = std.testing.allocator;
    const bytes = [_]u8{ 7, 2, 9 };

    try std.testing.expectError(
        ParseError.InvalidColumnData,
        decodeColumn(allocator, &bytes, &[_]Encoding{
            .{ .run_length = .{ .src_size = 2 } },
            .{ .byte_array = .{ .type_code = 4 } },
        }),
    );
}

test "bcif mask rejects mismatched lengths" {
    const allocator = std.testing.allocator;
    const data = [_]u8{ 1, 2 };
    const mask_data = [_]u8{0};

    try std.testing.expectError(
        ParseError.ColumnLengthMismatch,
        decodeColumnWithMask(
            allocator,
            &data,
            &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
            &mask_data,
            &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
        ),
    );
}

test "bcif mask rejects non-integer masks" {
    const allocator = std.testing.allocator;
    const data = [_]u8{1};
    const mask_data = [_]u8{ 0, 0, 0, 0 };

    try std.testing.expectError(
        ParseError.InvalidColumnData,
        decodeColumnWithMask(
            allocator,
            &data,
            &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
            &mask_data,
            &[_]Encoding{.{ .byte_array = .{ .type_code = 32 } }},
        ),
    );
}

test "bcif string arrays reject bad offsets" {
    const allocator = std.testing.allocator;
    const index_bytes = [_]u8{0};
    const offset_bytes = [_]u8{ 2, 1 };
    const encoding = Encoding{ .string_array = .{
        .string_data = "ABC",
        .offset_data = &offset_bytes,
        .offset_encoding = &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
        .data_encoding = &[_]Encoding{.{ .byte_array = .{ .type_code = 4 } }},
    } };

    try std.testing.expectError(
        ParseError.InvalidColumnData,
        decodeColumn(allocator, &index_bytes, &[_]Encoding{encoding}),
    );
}

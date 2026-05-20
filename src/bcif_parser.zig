const std = @import("std");
const Allocator = std.mem.Allocator;
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

fn decodeColumn(allocator: Allocator, data: []const u8, encodings: []const Encoding) DecodeError!DecodedColumn {
    var column = DecodedColumn{ .values = &.{} };
    var initialized = false;
    errdefer if (initialized) column.deinit(allocator);

    var i = encodings.len;
    while (i > 0) {
        i -= 1;
        if (!initialized) {
            column = try decodeInitialEncoding(allocator, data, encodings[i]);
            initialized = true;
        } else {
            column = try applyEncoding(allocator, column, encodings[i]);
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
    var column = try decodeColumn(allocator, data, encodings);
    errdefer column.deinit(allocator);

    var mask = try decodeColumn(allocator, mask_data, mask_encodings);
    defer mask.deinit(allocator);
    if (mask.values.len != column.values.len) return ParseError.ColumnLengthMismatch;

    const nulls = try allocator.alloc(NullKind, column.values.len);
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
    column.nulls = nulls;
    return column;
}

fn decodeInitialEncoding(allocator: Allocator, data: []const u8, encoding: Encoding) DecodeError!DecodedColumn {
    return switch (encoding) {
        .byte_array => |params| decodeByteArray(allocator, data, params.type_code),
        .string_array => |params| decodeStringArray(allocator, data, params),
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
        .byte_array, .string_array => ParseError.UnsupportedEncoding,
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
        acc += part;
        if (part != sentinel.positive and part != sentinel.negative) {
            try out.append(allocator, .{ .int = acc });
            acc = 0;
        }
    }
    if (out.items.len != src_size) return ParseError.InvalidColumnData;
    return .{ .values = try out.toOwnedSlice(allocator) };
}

fn packingSentinel(byte_count: u8, is_unsigned: bool) ?struct { positive: i64, negative: i64 } {
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
        current += scalarToI64(value) orelse return ParseError.InvalidColumnData;
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

fn decodeStringArray(allocator: Allocator, data: []const u8, params: @FieldType(Encoding, "string_array")) DecodeError!DecodedColumn {
    var indices = try decodeColumn(allocator, data, params.data_encoding);
    defer indices.deinit(allocator);
    var offsets = try decodeColumn(allocator, params.offset_data, params.offset_encoding);
    defer offsets.deinit(allocator);
    if (offsets.values.len == 0) return ParseError.InvalidColumnData;

    const values = try allocator.alloc(Scalar, indices.values.len);
    errdefer allocator.free(values);
    for (indices.values, values) |index_value, *out| {
        const index = scalarToUsize(index_value) orelse return ParseError.InvalidColumnData;
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

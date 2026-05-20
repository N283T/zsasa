const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const mmap_reader = @import("mmap_reader.zig");
const compressed = @import("compressed.zig");
const types = @import("types.zig");
const AtomInput = types.AtomInput;

comptime {
    _ = elem;
    _ = mmap_reader;
    _ = compressed;
    _ = AtomInput;
}

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

# BinaryCIF Input Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add native `.bcif`, `.bcif.gz`, and `.bcif.zst` input support for `zsasa calc` and `zsasa batch` by decoding BinaryCIF `_atom_site` into the existing `AtomInput` pipeline.

**Architecture:** Add a focused `bcif_parser` module that contains a small MessagePack reader, BinaryCIF column decoders, and `_atom_site` row conversion with the same filters as `MmcifParser`. Wire the new parser into format detection, `calc`, `batch`, and public module exports; keep `traj` out of scope except for explicit unsupported handling if compilation requires it.

**Tech Stack:** Zig 0.16.0+, native `std` only, existing `compressed.zig`, `mmap_reader.zig`, `types.zig`, `element.zig`, CLI docs in Markdown.

---

## File structure

- Create `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`: owns BinaryCIF parsing, MessagePack subset reading, BinaryCIF column decoding, `_atom_site` filtering, and unit tests for these pieces.
- Modify `/Users/nagaet/freesasa-zig/src/root.zig`: exports `bcif_parser` and includes it in compile-reference tests.
- Modify `/Users/nagaet/freesasa-zig/src/format_detect.zig`: adds `.bcif` input kind and supported filename extensions.
- Modify `/Users/nagaet/freesasa-zig/src/calc.zig`: routes `InputFormat.bcif` to `BcifParser` with existing `calc` filters.
- Modify `/Users/nagaet/freesasa-zig/src/batch.zig`: routes `InputFormat.bcif` to `BcifParser` with existing batch filters and includes files in discovery through `format_detect`.
- Modify `/Users/nagaet/freesasa-zig/src/traj.zig` only if adding `InputFormat.bcif` makes existing switches non-exhaustive; return a clear unsupported topology-format error instead of parsing BCIF topology.
- Add small fixtures under `/Users/nagaet/freesasa-zig/test_data/bcif/` only if synthetic in-test byte construction becomes too large for readable tests.
- Modify docs: `/Users/nagaet/freesasa-zig/README.md`, `/Users/nagaet/freesasa-zig/website/docs/cli/input.md`, `/Users/nagaet/freesasa-zig/website/docs/python-api/core.md`, `/Users/nagaet/freesasa-zig/website/docs/changelog.md`.

## Interfaces to preserve

Use this public parser shape in `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`:

```zig
pub const BcifParser = struct {
    allocator: Allocator,
    atom_only: bool = true,
    skip_hydrogens: bool = true,
    first_alt_loc_only: bool = true,
    model_num: ?u32 = null,
    chain_filter: ?[]const []const u8 = null,
    use_auth_chain: bool = false,

    pub fn init(allocator: Allocator) BcifParser {
        return .{ .allocator = allocator };
    }

    pub fn parse(self: *BcifParser, source: []const u8) !AtomInput {
        // Body is introduced in Task 3 after the reader and decoders exist.
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
};
```

Reuse the field preferences and filter semantics from `MmcifParser`:

- atom name: `label_atom_id` then `auth_atom_id`
- residue name: `label_comp_id` then `auth_comp_id`
- chain: `auth_asym_id` then `label_asym_id` when `use_auth_chain`; otherwise `label_asym_id` then `auth_asym_id`
- residue number: `label_seq_id` then `auth_seq_id`
- insertion code: `pdbx_PDB_ins_code`

---

### Task 1: Add MessagePack subset reader

**Files:**
- Create: `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`

- [ ] **Step 1: Write failing MessagePack reader tests**

Add `/Users/nagaet/freesasa-zig/src/bcif_parser.zig` with imports, error set, a `MsgValue` union skeleton, a `MsgReader` skeleton, and these tests. Use these exact test names so later steps can run focused commands:

```zig
const std = @import("std");
const Allocator = std.mem.Allocator;
const elem = @import("element.zig");
const mmap_reader = @import("mmap_reader.zig");
const compressed = @import("compressed.zig");
const types = @import("types.zig");
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

    fn readValue(self: *MsgReader) !MsgValue {
        _ = self;
        return ParseError.InvalidMessagePack;
    }
};

test "msgpack reader decodes primitives" {
    const bytes = [_]u8{
        0x88,
        0xa3, 'n', 'i', 'l', 0xc0,
        0xa4, 't', 'r', 'u', 'e', 0xc3,
        0xa5, 'f', 'a', 'l', 's', 'e', 0xc2,
        0xa3, 'i', 'n', 't', 0xd0, 0xfe,
        0xa4, 'u', 'i', 'n', 't', 0xcc, 200,
        0xa3, 'f', '3', '2', 0xca, 0x3f, 0x80, 0x00, 0x00,
        0xa3, 's', 't', 'r', 0xa3, 'a', 'b', 'c',
        0xa3, 'b', 'i', 'n', 0xc4, 0x03, 1, 2, 3,
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
```

Include helper signatures below the skeleton so the tests compile after implementation:

```zig
fn freeMsgValue(allocator: Allocator, value: MsgValue) void {
    _ = allocator;
    _ = value;
}

fn getMapValue(map: []const MsgPair, key: []const u8) ?MsgValue {
    _ = map;
    _ = key;
    return null;
}
```

- [ ] **Step 2: Run the focused tests to verify they fail for missing implementation**

Run:

```bash
zig build test --summary all 2>&1 | rg "msgpack reader|InvalidMessagePack|FAIL"
```

Expected: tests with names `msgpack reader decodes primitives` and `msgpack reader decodes arrays and 16 bit lengths` fail because `MsgReader.readValue()` returns `InvalidMessagePack`.

- [ ] **Step 3: Implement the MessagePack reader**

Replace the skeletons in `/Users/nagaet/freesasa-zig/src/bcif_parser.zig` with a real reader that supports the BinaryCIF subset:

```zig
fn readByte(self: *MsgReader) !u8 {
    if (self.pos >= self.data.len) return ParseError.UnexpectedEof;
    const b = self.data[self.pos];
    self.pos += 1;
    return b;
}

fn readBytes(self: *MsgReader, len: usize) ![]const u8 {
    if (self.pos + len > self.data.len) return ParseError.UnexpectedEof;
    const out = self.data[self.pos .. self.pos + len];
    self.pos += len;
    return out;
}

fn readIntBig(comptime T: type, self: *MsgReader) !T {
    var buf: [@sizeOf(T)]u8 = undefined;
    @memcpy(&buf, try self.readBytes(@sizeOf(T)));
    return std.mem.readInt(T, &buf, .big);
}

fn readValue(self: *MsgReader) !MsgValue {
    const tag = try self.readByte();
    if (tag <= 0x7f) return .{ .uint = tag };
    if (tag >= 0xe0) return .{ .int = @as(i8, @bitCast(tag)) };
    if (tag >= 0xa0 and tag <= 0xbf) return .{ .str = try self.readBytes(tag & 0x1f) };
    if (tag >= 0x90 and tag <= 0x9f) return self.readArray(tag & 0x0f);
    if (tag >= 0x80 and tag <= 0x8f) return self.readMap(tag & 0x0f);

    return switch (tag) {
        0xc0 => .nil,
        0xc2 => .{ .bool = false },
        0xc3 => .{ .bool = true },
        0xc4 => .{ .bin = try self.readBytes(try self.readIntBig(u8)) },
        0xc5 => .{ .bin = try self.readBytes(try self.readIntBig(u16)) },
        0xc6 => .{ .bin = try self.readBytes(try self.readIntBig(u32)) },
        0xca => blk: {
            const bits = try self.readIntBig(u32);
            break :blk .{ .float = @as(f64, @floatCast(@as(f32, @bitCast(bits)))) };
        },
        0xcb => blk: {
            const bits = try self.readIntBig(u64);
            break :blk .{ .float = @as(f64, @bitCast(bits)) };
        },
        0xcc => .{ .uint = try self.readIntBig(u8) },
        0xcd => .{ .uint = try self.readIntBig(u16) },
        0xce => .{ .uint = try self.readIntBig(u32) },
        0xcf => .{ .uint = try self.readIntBig(u64) },
        0xd0 => .{ .int = try self.readIntBig(i8) },
        0xd1 => .{ .int = try self.readIntBig(i16) },
        0xd2 => .{ .int = try self.readIntBig(i32) },
        0xd3 => .{ .int = try self.readIntBig(i64) },
        0xd9 => .{ .str = try self.readBytes(try self.readIntBig(u8)) },
        0xda => .{ .str = try self.readBytes(try self.readIntBig(u16)) },
        0xdb => .{ .str = try self.readBytes(try self.readIntBig(u32)) },
        0xdc => self.readArray(try self.readIntBig(u16)),
        0xdd => self.readArray(try self.readIntBig(u32)),
        0xde => self.readMap(try self.readIntBig(u16)),
        0xdf => self.readMap(try self.readIntBig(u32)),
        else => ParseError.UnsupportedMessagePackType,
    };
}

fn readArray(self: *MsgReader, len: usize) !MsgValue {
    const items = try self.allocator.alloc(MsgValue, len);
    errdefer self.allocator.free(items);
    for (items) |*item| item.* = try self.readValue();
    return .{ .array = items };
}

fn readMap(self: *MsgReader, len: usize) !MsgValue {
    const pairs = try self.allocator.alloc(MsgPair, len);
    errdefer self.allocator.free(pairs);
    for (pairs) |*pair| pair.* = .{ .key = try self.readValue(), .value = try self.readValue() };
    return .{ .map = pairs };
}
```

Implement `freeMsgValue()` and `getMapValue()`:

```zig
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
        if (pair.key == .str and std.mem.eql(u8, pair.key.str, key)) return pair.value;
    }
    return null;
}
```

- [ ] **Step 4: Run focused tests to verify they pass**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bcif-task1.log
rg "msgpack reader" /tmp/zsasa-bcif-task1.log
```

Expected: `zig build test` exits 0. The log includes the MessagePack test names as passed or the Zig summary reports all tests passed.

- [ ] **Step 5: Commit Task 1**

```bash
git add /Users/nagaet/freesasa-zig/src/bcif_parser.zig
git commit -m "feat: add binary cif messagepack reader"
```

---

### Task 2: Add BinaryCIF column decoders

**Files:**
- Modify: `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`

- [ ] **Step 1: Add failing decoder tests**

Append these tests to `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`:

```zig
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
    const packed = [_]u8{ 1, 3, 2, 1 };
    const encs = [_]Encoding{
        .{ .delta = .{ .origin = 10 } },
        .{ .run_length = .{ .src_size = 4 } },
        .{ .integer_packing = .{ .byte_count = 1, .is_unsigned = false, .src_size = 4 } },
        .{ .byte_array = .{ .type_code = 1 } },
    };
    var column = try decodeColumn(allocator, &packed, &encs);
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
```

- [ ] **Step 2: Run tests to verify decoder APIs are missing**

Run:

```bash
zig build test --summary all 2>&1 | rg "decodeColumn|Encoding|NullKind|FAIL|error:"
```

Expected: compilation fails because `Encoding`, `decodeColumn()`, `decodeColumnWithMask()`, and `NullKind` are not defined.

- [ ] **Step 3: Implement decoder data types**

Add these types above the tests:

```zig
const NullKind = enum { present, not_present, unknown };

const Scalar = union(enum) {
    int: i64,
    uint: u64,
    float: f64,
    str: []const u8,
};

const DecodedColumn = struct {
    values: []Scalar,
    nulls: ?[]NullKind = null,

    fn deinit(self: *DecodedColumn, allocator: Allocator) void {
        allocator.free(self.values);
        if (self.nulls) |n| allocator.free(n);
    }
};

const Encoding = union(enum) {
    byte_array: struct { type_code: u8 },
    fixed_point: struct { factor: f64 },
    interval_quantization: struct { min: f64, max: f64, num_steps: u64 },
    run_length: struct { src_size: usize },
    delta: struct { origin: i64 },
    integer_packing: struct { byte_count: u8, src_size: usize, is_unsigned: bool },
    string_array: struct {
        string_data: []const u8,
        offset_data: []const u8,
        offset_encoding: []const Encoding,
        data_encoding: []const Encoding,
    },
};
```

- [ ] **Step 4: Implement decoder functions**

Add decoder functions that operate on `[]const u8` and `[]const Encoding`. The implementation must apply encodings from last to first. Use explicit intermediate arrays owned by the allocator; free replaced intermediates before returning.

Key function signatures:

```zig
fn decodeColumn(allocator: Allocator, data: []const u8, encodings: []const Encoding) !DecodedColumn;
fn decodeColumnWithMask(allocator: Allocator, data: []const u8, encodings: []const Encoding, mask_data: []const u8, mask_encodings: []const Encoding) !DecodedColumn;
fn decodeByteArray(allocator: Allocator, data: []const u8, type_code: u8) ![]Scalar;
fn decodeIntegerPacking(allocator: Allocator, input: []const Scalar, byte_count: u8, is_unsigned: bool) ![]Scalar;
fn decodeRunLength(allocator: Allocator, input: []const Scalar, src_size: usize) ![]Scalar;
fn decodeDelta(allocator: Allocator, input: []const Scalar, origin: i64) ![]Scalar;
fn decodeFixedPoint(allocator: Allocator, input: []const Scalar, factor: f64) ![]Scalar;
fn decodeIntervalQuantization(allocator: Allocator, input: []const Scalar, min: f64, max: f64, num_steps: u64) ![]Scalar;
fn decodeStringArray(allocator: Allocator, data: []const u8, spec: anytype) ![]Scalar;
```

Use these BinaryCIF type codes in `decodeByteArray()`:

```zig
// 1 Int8, 2 Int16, 3 Int32, 4 Uint8, 5 Uint16, 6 Uint32, 32 Float32, 33 Float64
```

Use little-endian reads:

```zig
const value = std.mem.readInt(i16, data[i .. i + 2][0..2], .little);
```

For integer packing, reproduce the BinaryCIF rule:

```zig
const limit: i64 = if (is_unsigned)
    (if (byte_count == 1) 255 else 65535)
else
    (if (byte_count == 1) 127 else 32767);
const lower: i64 = if (byte_count == 1) -128 else -32768;
```

Continue accumulating while a packed value equals the positive limit, or for signed input equals the negative limit. Emit the accumulated value plus the first non-limit value.

- [ ] **Step 5: Run focused tests to verify decoders pass**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bcif-task2.log
rg "bcif decoders|bcif mask" /tmp/zsasa-bcif-task2.log
```

Expected: `zig build test` exits 0.

- [ ] **Step 6: Commit Task 2**

```bash
git add /Users/nagaet/freesasa-zig/src/bcif_parser.zig
git commit -m "feat: decode binary cif columns"
```

---

### Task 3: Parse BinaryCIF container and `_atom_site`

**Files:**
- Modify: `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`

- [ ] **Step 1: Add failing parser tests with synthetic BinaryCIF bytes**

Add a test helper that builds a minimal BinaryCIF MessagePack object using `ByteArray` and `StringArray` encodings. Place helpers near tests so they are not exported.

Add these tests:

```zig
test "parse minimal BinaryCIF atom_site" {
    const source = try buildMinimalBcif(std.testing.allocator, .{});
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.len());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.000001);
    try std.testing.expectApproxEqAbs(@as(f64, 21.0), input.y[1], 0.000001);
    try std.testing.expectApproxEqAbs(@as(f64, 32.0), input.z[1], 0.000001);
    try std.testing.expectEqual(@as(u8, 6), input.element[0]);
    try std.testing.expectEqualStrings("ALA", input.residue[0].slice());
    try std.testing.expectEqualStrings("CA", input.atom_name[0].slice());
    try std.testing.expectEqualStrings("A", input.chain_id[0].slice());
    try std.testing.expectEqual(@as(i32, 1), input.residue_num[0]);
}

test "parse BinaryCIF applies filters" {
    const source = try buildMinimalBcif(std.testing.allocator, .{ .include_hydrogen = true, .include_hetatm = true, .include_second_model = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    parser.model_num = 1;
    parser.chain_filter = &[_][]const u8{"A"};
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.len());
    for (input.chain_id) |chain| try std.testing.expectEqualStrings("A", chain.slice());
}

test "parse BinaryCIF reports missing coordinate field" {
    const source = try buildMinimalBcif(std.testing.allocator, .{ .omit_z = true });
    defer std.testing.allocator.free(source);

    var parser = BcifParser.init(std.testing.allocator);
    try std.testing.expectError(ParseError.MissingCoordinateField, parser.parse(source));
}
```

Define `BuildBcifOptions` for the helper:

```zig
const BuildBcifOptions = struct {
    include_hydrogen: bool = false,
    include_hetatm: bool = false,
    include_second_model: bool = false,
    omit_z: bool = false,
};
```

- [ ] **Step 2: Run tests to verify parser APIs are missing**

Run:

```bash
zig build test --summary all 2>&1 | rg "buildMinimalBcif|BcifParser|MissingCoordinateField|FAIL|error:"
```

Expected: compilation fails because `BcifParser.parse()` and `buildMinimalBcif()` are incomplete.

- [ ] **Step 3: Implement `BcifParser` fields and `parseFile()`**

Add the `BcifParser` type with the exact public fields from the spec. Use the `parseFile()` implementation from the interface section above.

- [ ] **Step 4: Implement BinaryCIF structure extraction**

Add helper functions with these signatures:

```zig
fn parseEncoding(allocator: Allocator, value: MsgValue) !Encoding;
fn parseEncodingList(allocator: Allocator, value: MsgValue) ![]Encoding;
fn findCategory(root: MsgValue, name: []const u8) !MsgValue;
fn categoryRowCount(category: MsgValue) !usize;
fn categoryColumns(category: MsgValue) ![]const MsgValue;
fn columnName(column: MsgValue) ![]const u8;
fn columnData(column: MsgValue) !MsgValue;
fn columnMask(column: MsgValue) ?MsgValue;
fn decodeBcifColumn(allocator: Allocator, column: MsgValue) !DecodedColumn;
```

`parseEncoding()` must handle maps with `kind` values `ByteArray`, `FixedPoint`, `IntervalQuantization`, `RunLength`, `Delta`, `IntegerPacking`, and `StringArray`. For `StringArray`, decode nested `dataEncoding`, `offsetEncoding`, `offsets`, and `stringData` fields into the `Encoding.string_array` variant.

- [ ] **Step 5: Implement `_atom_site` column mapping and row conversion**

Add an `AtomSiteColumns` equivalent for decoded columns:

```zig
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
    pdbx_pdb_model_num: ?usize = null,

    fn hasRequiredFields(self: AtomSiteColumns) bool {
        return self.cartn_x != null and self.cartn_y != null and self.cartn_z != null;
    }
};
```

Implement helper accessors mirroring `MmcifParser`:

```zig
fn getAtomNameCol(columns: AtomSiteColumns) ?usize { return columns.label_atom_id orelse columns.auth_atom_id; }
fn getResNameCol(columns: AtomSiteColumns) ?usize { return columns.label_comp_id orelse columns.auth_comp_id; }
fn getChainCol(columns: AtomSiteColumns, use_auth: bool) ?usize {
    return if (use_auth) columns.auth_asym_id orelse columns.label_asym_id else columns.label_asym_id orelse columns.auth_asym_id;
}
fn getResSeqCol(columns: AtomSiteColumns) ?usize { return columns.label_seq_id orelse columns.auth_seq_id; }
```

Implement row accessors:

```zig
fn scalarString(column: DecodedColumn, row: usize) ?[]const u8;
fn scalarFloat(column: DecodedColumn, row: usize) !f64;
fn scalarInt(column: DecodedColumn, row: usize) ?i64;
fn isNull(column: DecodedColumn, row: usize) bool;
```

`scalarString()` returns `null` for mask values `not_present` and `unknown`. `scalarFloat()` returns `ParseError.InvalidCoordinate` for masked coordinates.

- [ ] **Step 6: Implement filter logic**

Port `MmcifParser.shouldIncludeAtom()` into `BcifParser` using decoded columns. Keep behavior identical:

```zig
fn shouldIncludeAtom(self: *BcifParser, decoded: []const DecodedColumn, columns: AtomSiteColumns, row: usize, first_alt_loc: *?u8) bool;
```

Rules:

- If `atom_only` and `group_PDB` exists, include only exact string `ATOM`.
- If `skip_hydrogens` and `type_symbol` exists, exclude exact strings `H` and `D`.
- If `first_alt_loc_only` and `label_alt_id` exists, remember the first non-null alt ID byte and exclude other alt IDs.
- If `model_num` is set and `pdbx_PDB_model_num` exists, compare parsed integer model number, defaulting invalid model values to `1` to match mmCIF behavior.
- If `chain_filter` is set and a chain column exists, include only exact chain matches.

- [ ] **Step 7: Implement `parse()` to build `AtomInput`**

`parse()` flow:

```zig
pub fn parse(self: *BcifParser, source: []const u8) !AtomInput {
    var reader = MsgReader.init(self.allocator, source);
    const root = try reader.readValue();
    defer freeMsgValue(self.allocator, root);

    const atom_site = try findCategory(root, "atom_site");
    const row_count = try categoryRowCount(atom_site);
    const columns_value = try categoryColumns(atom_site);

    // map names, decode selected columns, validate required coordinates,
    // append included rows to ArrayListUnmanaged values,
    // return AtomInput with owned slices.
}
```

Use the same fallback metadata defaults as `MmcifParser`:

- missing residue name -> `UNK`
- missing atom name -> `X`
- missing chain -> empty string
- missing residue number -> `0`
- missing insertion code -> empty string
- missing element -> `Element.X`

- [ ] **Step 8: Implement synthetic BCIF builder for tests**

Implement `buildMinimalBcif()` using local MessagePack packing helpers in tests. It should create a top-level map equivalent to:

```text
{
  "version": "0.3.0",
  "encoder": "zsasa-test",
  "dataBlocks": [
    {
      "header": "TEST",
      "categories": [
        {
          "name": "_atom_site",
          "rowCount": N,
          "columns": [ ... selected columns ... ]
        }
      ]
    }
  ]
}
```

Use `ByteArray` for numeric columns and `StringArray` for string columns. Add helper packers:

```zig
fn packMapHeader(list: *std.ArrayListUnmanaged(u8), allocator: Allocator, len: usize) !void;
fn packArrayHeader(list: *std.ArrayListUnmanaged(u8), allocator: Allocator, len: usize) !void;
fn packStr(list: *std.ArrayListUnmanaged(u8), allocator: Allocator, value: []const u8) !void;
fn packBin(list: *std.ArrayListUnmanaged(u8), allocator: Allocator, value: []const u8) !void;
fn packUint(list: *std.ArrayListUnmanaged(u8), allocator: Allocator, value: u64) !void;
fn packBool(list: *std.ArrayListUnmanaged(u8), allocator: Allocator, value: bool) !void;
```

- [ ] **Step 9: Run parser tests**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bcif-task3.log
rg "parse .*BinaryCIF" /tmp/zsasa-bcif-task3.log
```

Expected: `zig build test` exits 0.

- [ ] **Step 10: Commit Task 3**

```bash
git add /Users/nagaet/freesasa-zig/src/bcif_parser.zig
git commit -m "feat: parse binary cif atom site"
```

---

### Task 4: Wire BinaryCIF into format detection, root exports, calc, and batch

**Files:**
- Modify: `/Users/nagaet/freesasa-zig/src/format_detect.zig`
- Modify: `/Users/nagaet/freesasa-zig/src/root.zig`
- Modify: `/Users/nagaet/freesasa-zig/src/calc.zig`
- Modify: `/Users/nagaet/freesasa-zig/src/batch.zig`
- Modify: `/Users/nagaet/freesasa-zig/src/traj.zig` if switch exhaustiveness requires it

- [ ] **Step 1: Add failing format detection tests**

In `/Users/nagaet/freesasa-zig/src/format_detect.zig`, update tests before implementation:

```zig
test "isSupportedFile accepts BinaryCIF files" {
    try std.testing.expect(isSupportedFile("structure.bcif"));
    try std.testing.expect(isSupportedFile("structure.bcif.gz"));
    try std.testing.expect(isSupportedFile("structure.bcif.zst"));
    try std.testing.expect(isSupportedFile("structure.BCIF"));
}

test "detectInputFormat handles BinaryCIF files" {
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.bcif"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.bcif.gz"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.bcif.zst"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.BCIF"));
}
```

- [ ] **Step 2: Run tests to verify `InputFormat.bcif` is missing**

Run:

```bash
zig build test --summary all 2>&1 | rg "bcif|InputFormat|error:"
```

Expected: compilation fails because `InputFormat.bcif` does not exist.

- [ ] **Step 3: Implement `InputFormat.bcif`**

Change enum and extension list in `/Users/nagaet/freesasa-zig/src/format_detect.zig`:

```zig
pub const InputFormat = enum {
    json,
    mmcif,
    bcif,
    pdb,
    sdf,
};
```

Add to `supported_extensions`:

```zig
".bcif",
".bcif.gz",
".bcif.zst",
```

Add uppercase support in `isSupportedFile()`:

```zig
if (std.mem.endsWith(u8, name, ".BCIF")) return true;
```

Add detection before mmCIF detection:

```zig
if (std.mem.endsWith(u8, base, ".bcif")) return .bcif;
if (std.mem.endsWith(u8, base, ".BCIF")) return .bcif;
```

- [ ] **Step 4: Export the parser from root**

In `/Users/nagaet/freesasa-zig/src/root.zig`, add:

```zig
pub const bcif_parser = @import("bcif_parser.zig");
```

In the compile-reference test block, add:

```zig
_ = bcif_parser;
```

- [ ] **Step 5: Wire `calc`**

In `/Users/nagaet/freesasa-zig/src/calc.zig`, add import near other parsers:

```zig
const bcif_parser = @import("bcif_parser.zig");
```

In `readInputFile()`, add a `.bcif` switch case next to `.mmcif`:

```zig
.bcif => blk: {
    var parser = bcif_parser.BcifParser.init(allocator);
    parser.model_num = args.model_num;
    parser.use_auth_chain = args.use_auth_chain;
    parser.skip_hydrogens = !args.include_hydrogens;
    parser.atom_only = !args.include_hetatm;

    var chain_filter_slice: ?[]const []const u8 = null;
    if (args.chain_filter) |filter_str| {
        chain_filter_slice = try parseChainFilter(allocator, filter_str);
        parser.chain_filter = chain_filter_slice;
    }
    defer if (chain_filter_slice) |s| allocator.free(s);

    break :blk .{ .input = try parser.parseFile(io, path) };
},
```

- [ ] **Step 6: Wire `batch`**

In `/Users/nagaet/freesasa-zig/src/batch.zig`, add import near other parsers:

```zig
const bcif_parser = @import("bcif_parser.zig");
```

In `readInputFile()`, add:

```zig
.bcif => blk: {
    var parser = bcif_parser.BcifParser.init(allocator);
    parser.skip_hydrogens = !config.include_hydrogens;
    parser.atom_only = !config.include_hetatm;
    parser.chain_filter = config.chain_filter;
    parser.use_auth_chain = config.use_auth_chain;
    break :blk parser.parseFile(io, path);
},
```

- [ ] **Step 7: Handle `traj` switch exhaustiveness if needed**

Run `zig build test` after adding `.bcif`. If `src/traj.zig` has a non-exhaustive switch over `InputFormat`, add an explicit case that returns the same kind of unsupported-format error used for JSON/SDF topology. If there is no switch error, leave `traj.zig` unchanged.

- [ ] **Step 8: Run tests**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bcif-task4.log
```

Expected: `zig build test` exits 0 and the new format detection tests pass.

- [ ] **Step 9: Commit Task 4**

```bash
git add /Users/nagaet/freesasa-zig/src/format_detect.zig /Users/nagaet/freesasa-zig/src/root.zig /Users/nagaet/freesasa-zig/src/calc.zig /Users/nagaet/freesasa-zig/src/batch.zig /Users/nagaet/freesasa-zig/src/traj.zig
git commit -m "feat: route binary cif inputs"
```

---

### Task 5: Add compressed-file and CLI smoke coverage

**Files:**
- Modify: `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`
- Create if needed: `/Users/nagaet/freesasa-zig/test_data/bcif/minimal.bcif`
- Create if needed: `/Users/nagaet/freesasa-zig/test_data/bcif/minimal.bcif.gz`
- Create if needed: `/Users/nagaet/freesasa-zig/test_data/bcif/minimal.bcif.zst`

- [ ] **Step 1: Add parseFile tests for plain, gzip, and zstd**

If the synthetic builder can write files during tests, keep fixtures out of git and add this test to `/Users/nagaet/freesasa-zig/src/bcif_parser.zig`:

```zig
test "parseFile reads plain gzip and zstd BinaryCIF" {
    const source = try buildMinimalBcif(std.testing.allocator, .{});
    defer std.testing.allocator.free(source);

    var tmp = std.testing.tmpDir(.{});
    defer tmp.cleanup();

    try tmp.dir.writeFile(.{ .sub_path = "minimal.bcif", .data = source });
    // If project compression helpers expose writers, use them here. Otherwise create committed fixtures in test_data/bcif and parse those paths.

    var parser = BcifParser.init(std.testing.allocator);
    var input = try parser.parseFile(std.testing.io, "test_data/bcif/minimal.bcif");
    defer input.deinit();
    try std.testing.expectEqual(@as(usize, 2), input.len());
}
```

If writing gzip/zstd in Zig tests is not practical with existing helpers, create small committed fixtures with Python using `py-mmcif` or by writing the synthetic builder output from a short local script. Keep each fixture under 50 KB.

- [ ] **Step 2: Build the CLI**

Run:

```bash
zig build -Doptimize=ReleaseFast
```

Expected: command exits 0 and produces `/Users/nagaet/freesasa-zig/zig-out/bin/zsasa`.

- [ ] **Step 3: Run CLI smoke commands**

Run these commands after creating or committing fixtures:

```bash
./zig-out/bin/zsasa calc test_data/bcif/minimal.bcif /tmp/zsasa-bcif.json
./zig-out/bin/zsasa calc test_data/bcif/minimal.bcif.gz /tmp/zsasa-bcif-gz.json
./zig-out/bin/zsasa calc test_data/bcif/minimal.bcif.zst /tmp/zsasa-bcif-zst.json
```

Expected: each command exits 0 and writes a JSON output file.

- [ ] **Step 4: Run batch discovery smoke command**

Run:

```bash
rm -rf /tmp/zsasa-bcif-batch
mkdir -p /tmp/zsasa-bcif-batch
./zig-out/bin/zsasa batch test_data/bcif /tmp/zsasa-bcif-batch --quiet
find /tmp/zsasa-bcif-batch -type f | sort
```

Expected: outputs exist for the `.bcif`, `.bcif.gz`, and `.bcif.zst` fixtures.

- [ ] **Step 5: Run full focused verification**

Run:

```bash
zig fmt --check src/
zig build test
zig build -Doptimize=ReleaseFast
```

Expected: all commands exit 0.

- [ ] **Step 6: Commit Task 5**

```bash
git add /Users/nagaet/freesasa-zig/src/bcif_parser.zig /Users/nagaet/freesasa-zig/test_data/bcif
git commit -m "test: add binary cif smoke coverage"
```

---

### Task 6: Update documentation

**Files:**
- Modify: `/Users/nagaet/freesasa-zig/README.md`
- Modify: `/Users/nagaet/freesasa-zig/website/docs/cli/input.md`
- Modify: `/Users/nagaet/freesasa-zig/website/docs/python-api/core.md`
- Modify: `/Users/nagaet/freesasa-zig/website/docs/changelog.md`

- [ ] **Step 1: Update README supported format line**

In `/Users/nagaet/freesasa-zig/README.md`, change the multiple-input-formats bullet to include BinaryCIF:

```markdown
- **Multiple input formats**: mmCIF, BinaryCIF (`.bcif`), PDB, SDF/MOL, JSON, XTC, DCD
```

- [ ] **Step 2: Update CLI input table**

In `/Users/nagaet/freesasa-zig/website/docs/cli/input.md`, add a BinaryCIF row near mmCIF:

```markdown
| `.bcif`, `.bcif.gz`, `.bcif.zst`, `.BCIF` | BinaryCIF (`_atom_site` support) |
```

Add a short note below the table:

```markdown
BinaryCIF input currently decodes `_atom_site` for SASA calculation. Inline CCD data embedded in BinaryCIF is not used yet; provide external CCD or SDF topology when needed for non-standard compounds.
```

- [ ] **Step 3: Update Python batch format docs**

In `/Users/nagaet/freesasa-zig/website/docs/python-api/core.md`, update the supported directory formats sentence so it includes `.bcif`, `.bcif.gz`, and `.bcif.zst`:

```markdown
Process all supported structure files (.pdb, .cif, .mmcif, .bcif, .ent, .json, .sdf, .mol, and .gz/.zst variants) in a directory.
```

- [ ] **Step 4: Update changelog**

In `/Users/nagaet/freesasa-zig/website/docs/changelog.md`, add an unreleased bullet near other input-format entries:

```markdown
- **BinaryCIF input support**: `calc` and `batch` now accept `.bcif`, `.bcif.gz`, and `.bcif.zst` files by decoding `_atom_site` directly in Zig. Inline CCD extraction from BinaryCIF remains out of scope for this release.
```

- [ ] **Step 5: Run documentation-adjacent checks**

Run:

```bash
rg -n "BinaryCIF|bcif" README.md website/docs/cli/input.md website/docs/python-api/core.md website/docs/changelog.md
```

Expected: all four files contain BinaryCIF or `.bcif` references.

- [ ] **Step 6: Commit Task 6**

```bash
git add /Users/nagaet/freesasa-zig/README.md /Users/nagaet/freesasa-zig/website/docs/cli/input.md /Users/nagaet/freesasa-zig/website/docs/python-api/core.md /Users/nagaet/freesasa-zig/website/docs/changelog.md
git commit -m "docs: document binary cif input support"
```

---

### Task 7: Final verification

**Files:**
- No planned source edits unless verification finds a defect.

- [ ] **Step 1: Run formatting check**

```bash
zig fmt --check src/
```

Expected: exits 0.

- [ ] **Step 2: Run Zig tests**

```bash
zig build test
```

Expected: exits 0.

- [ ] **Step 3: Run release build**

```bash
zig build -Doptimize=ReleaseFast
```

Expected: exits 0.

- [ ] **Step 4: Run CLI smoke checks**

```bash
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa --version
./zig-out/bin/zsasa calc test_data/bcif/minimal.bcif /tmp/zsasa-bcif-final.json
```

Expected: help and version print successfully; BinaryCIF calc exits 0.

- [ ] **Step 5: Inspect git status and recent commits**

```bash
git status --short --branch
git log --oneline -8
```

Expected: branch is `feat/binary-cif-support`. Working tree is clean after final fixes are committed.

## Self-review notes

- Spec coverage: Tasks 1-3 cover native decoder and `_atom_site`; Task 4 covers `.bcif` format detection plus `calc` and `batch`; Task 5 covers gzip/zstd and smoke checks; Task 6 covers docs.
- Scope control: The plan does not add BinaryCIF writing, full CIF API, inline CCD extraction, or `traj` topology support.
- Type consistency: The public parser is consistently `BcifParser`; the format enum is consistently `InputFormat.bcif`; decoded values use `DecodedColumn`, `Scalar`, and `NullKind` throughout.

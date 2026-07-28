const std = @import("std");

const Allocator = std.mem.Allocator;

pub const AsymIdType = enum {
    label,
    auth,

    fn parse(value: []const u8) !AsymIdType {
        if (std.ascii.eqlIgnoreCase(value, "label")) return .label;
        if (std.ascii.eqlIgnoreCase(value, "auth")) return .auth;
        return error.InvalidAsymIdType;
    }
};

pub const Entry = struct {
    filename: []const u8,
    id: ?[]const u8,
    chains: []const []const u8,
    asym_id_type: AsymIdType,

    fn deinit(self: Entry, allocator: Allocator) void {
        allocator.free(self.filename);
        if (self.id) |id| allocator.free(id);
        for (self.chains) |chain| allocator.free(chain);
        allocator.free(self.chains);
    }
};

pub const ChainMap = struct {
    allocator: Allocator,
    entries: []Entry,
    index: std.StringHashMapUnmanaged(std.ArrayListUnmanaged(usize)),

    pub fn deinit(self: *ChainMap) void {
        var values = self.index.valueIterator();
        while (values.next()) |indices| indices.deinit(self.allocator);
        self.index.deinit(self.allocator);
        for (self.entries) |entry| entry.deinit(self.allocator);
        self.allocator.free(self.entries);
        self.* = undefined;
    }

    pub fn getIndices(self: *const ChainMap, filename: []const u8) ?[]const usize {
        const indices = self.index.get(filename) orelse return null;
        return indices.items;
    }

    /// Return the first selection for compatibility with one-row-per-file maps.
    pub fn get(self: *const ChainMap, filename: []const u8) ?Entry {
        const indices = self.getIndices(filename) orelse return null;
        return self.entries[indices[0]];
    }

    pub fn hasMultipleSelections(self: *const ChainMap) bool {
        var values = self.index.valueIterator();
        while (values.next()) |indices| {
            if (indices.items.len > 1) return true;
        }
        return false;
    }
};

pub const InterfaceEntry = struct {
    filename: []const u8,
    id: ?[]const u8,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
    asym_id_type: AsymIdType,

    fn deinit(self: InterfaceEntry, allocator: Allocator) void {
        allocator.free(self.filename);
        if (self.id) |id| allocator.free(id);
        freeChains(allocator, self.partner_a);
        freeChains(allocator, self.partner_b);
    }
};

pub const InterfaceMap = struct {
    allocator: Allocator,
    entries: []InterfaceEntry,
    index: std.StringHashMapUnmanaged(std.ArrayListUnmanaged(usize)),

    pub fn deinit(self: *InterfaceMap) void {
        var values = self.index.valueIterator();
        while (values.next()) |indices| indices.deinit(self.allocator);
        self.index.deinit(self.allocator);
        for (self.entries) |entry| entry.deinit(self.allocator);
        self.allocator.free(self.entries);
        self.* = undefined;
    }

    pub fn getIndices(self: *const InterfaceMap, filename: []const u8) ?[]const usize {
        const indices = self.index.get(filename) orelse return null;
        return indices.items;
    }

    /// Return the first interface for compatibility with one-row-per-file maps.
    pub fn get(self: *const InterfaceMap, filename: []const u8) ?InterfaceEntry {
        const indices = self.getIndices(filename) orelse return null;
        return self.entries[indices[0]];
    }
};

pub fn loadFile(allocator: Allocator, io: std.Io, path: []const u8) !ChainMap {
    const content = try readFile(allocator, io, path);
    defer allocator.free(content);

    if (endsWithIgnoreCase(path, ".csv")) return parseCsv(allocator, content);
    if (endsWithIgnoreCase(path, ".json")) return parseJson(allocator, content);
    return error.UnsupportedChainMapFormat;
}

pub fn loadInterfaceFile(allocator: Allocator, io: std.Io, path: []const u8) !InterfaceMap {
    const content = try readFile(allocator, io, path);
    defer allocator.free(content);

    if (endsWithIgnoreCase(path, ".csv")) return parseInterfaceCsv(allocator, content);
    if (endsWithIgnoreCase(path, ".json")) return parseInterfaceJson(allocator, content);
    return error.UnsupportedChainMapFormat;
}

pub fn parseCsv(allocator: Allocator, content: []const u8) !ChainMap {
    var entries = std.ArrayListUnmanaged(Entry).empty;
    errdefer deinitEntryList(allocator, &entries);

    var lines = std.mem.splitScalar(u8, content, '\n');
    var header: ?[]const []const u8 = null;
    defer if (header) |fields| freeFields(allocator, fields);

    var filename_col: ?usize = null;
    var id_col: ?usize = null;
    var chains_col: ?usize = null;
    var asym_id_type_col: ?usize = null;

    while (lines.next()) |raw_line| {
        const line = std.mem.trimEnd(u8, raw_line, "\r");
        if (std.mem.trim(u8, line, " \t").len == 0) continue;

        const fields = try parseCsvRecord(allocator, line);
        if (header == null) {
            header = fields;
            for (fields, 0..) |field, i| {
                const normalized = if (i == 0) std.mem.trimStart(u8, field, "\xEF\xBB\xBF") else field;
                if (std.ascii.eqlIgnoreCase(normalized, "filename")) filename_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "id")) id_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "chains")) chains_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "asym_id_type")) asym_id_type_col = i;
            }
            if (filename_col == null or chains_col == null) return error.MissingChainMapColumn;
            continue;
        }
        defer freeFields(allocator, fields);

        if (fields.len != header.?.len) return error.InvalidCsv;
        const required_max = @max(filename_col.?, chains_col.?);
        if (fields.len <= required_max) return error.MissingChainMapField;
        if (asym_id_type_col) |col| {
            if (fields.len <= col) return error.MissingChainMapField;
        }

        const asym_id_type = if (asym_id_type_col) |col|
            try parseAsymIdTypeOrDefault(fields[col])
        else
            AsymIdType.label;
        try appendEntry(
            allocator,
            &entries,
            fields[filename_col.?],
            if (id_col) |col| fields[col] else null,
            fields[chains_col.?],
            asym_id_type,
        );
    }

    if (header == null) return error.MissingChainMapHeader;
    return finish(allocator, &entries);
}

const JsonEntry = struct {
    filename: []const u8,
    id: ?[]const u8 = null,
    chains: []const []const u8,
    asym_id_type: ?[]const u8 = null,
};

pub fn parseJson(allocator: Allocator, content: []const u8) !ChainMap {
    const parsed = try std.json.parseFromSlice([]const JsonEntry, allocator, content, .{});
    defer parsed.deinit();

    var entries = std.ArrayListUnmanaged(Entry).empty;
    errdefer deinitEntryList(allocator, &entries);

    for (parsed.value) |json_entry| {
        const asym_id_type = try parseAsymIdTypeOrDefault(json_entry.asym_id_type orelse "");
        try appendEntryFromChains(allocator, &entries, json_entry.filename, json_entry.id, json_entry.chains, asym_id_type);
    }
    return finish(allocator, &entries);
}

pub fn parseInterfaceCsv(allocator: Allocator, content: []const u8) !InterfaceMap {
    var entries = std.ArrayListUnmanaged(InterfaceEntry).empty;
    errdefer deinitInterfaceEntryList(allocator, &entries);

    var lines = std.mem.splitScalar(u8, content, '\n');
    var header: ?[]const []const u8 = null;
    defer if (header) |fields| freeFields(allocator, fields);

    var filename_col: ?usize = null;
    var id_col: ?usize = null;
    var partner_a_col: ?usize = null;
    var partner_b_col: ?usize = null;
    var asym_id_type_col: ?usize = null;

    while (lines.next()) |raw_line| {
        const line = std.mem.trimEnd(u8, raw_line, "\r");
        if (std.mem.trim(u8, line, " \t").len == 0) continue;

        const fields = try parseCsvRecord(allocator, line);
        if (header == null) {
            header = fields;
            for (fields, 0..) |field, i| {
                const normalized = if (i == 0) std.mem.trimStart(u8, field, "\xEF\xBB\xBF") else field;
                if (std.ascii.eqlIgnoreCase(normalized, "filename")) filename_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "id")) id_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "partner_a")) partner_a_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "partner_b")) partner_b_col = i;
                if (std.ascii.eqlIgnoreCase(normalized, "asym_id_type")) asym_id_type_col = i;
            }
            if (filename_col == null or partner_a_col == null or partner_b_col == null) {
                return error.MissingChainMapColumn;
            }
            continue;
        }
        defer freeFields(allocator, fields);

        if (fields.len != header.?.len) return error.InvalidCsv;
        const required_max = @max(filename_col.?, @max(partner_a_col.?, partner_b_col.?));
        if (fields.len <= required_max) return error.MissingChainMapField;
        if (asym_id_type_col) |col| {
            if (fields.len <= col) return error.MissingChainMapField;
        }

        const asym_id_type = if (asym_id_type_col) |col|
            try parseAsymIdTypeOrDefault(fields[col])
        else
            AsymIdType.label;
        try appendInterfaceEntryCsv(
            allocator,
            &entries,
            fields[filename_col.?],
            if (id_col) |col| fields[col] else null,
            fields[partner_a_col.?],
            fields[partner_b_col.?],
            asym_id_type,
        );
    }

    if (header == null) return error.MissingChainMapHeader;
    return finishInterface(allocator, &entries);
}

const JsonInterfaceEntry = struct {
    filename: []const u8,
    id: ?[]const u8 = null,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
    asym_id_type: ?[]const u8 = null,
};

pub fn parseInterfaceJson(allocator: Allocator, content: []const u8) !InterfaceMap {
    const parsed = try std.json.parseFromSlice([]const JsonInterfaceEntry, allocator, content, .{});
    defer parsed.deinit();

    var entries = std.ArrayListUnmanaged(InterfaceEntry).empty;
    errdefer deinitInterfaceEntryList(allocator, &entries);

    for (parsed.value) |json_entry| {
        const asym_id_type = try parseAsymIdTypeOrDefault(json_entry.asym_id_type orelse "");
        try appendInterfaceEntryJson(
            allocator,
            &entries,
            json_entry.filename,
            json_entry.id,
            json_entry.partner_a,
            json_entry.partner_b,
            asym_id_type,
        );
    }
    return finishInterface(allocator, &entries);
}

fn parseAsymIdTypeOrDefault(value: []const u8) !AsymIdType {
    const trimmed = std.mem.trim(u8, value, " \t\r");
    if (trimmed.len == 0) return .label;
    return AsymIdType.parse(trimmed);
}

fn appendEntry(
    allocator: Allocator,
    entries: *std.ArrayListUnmanaged(Entry),
    filename_value: []const u8,
    id_value: ?[]const u8,
    chains_value: []const u8,
    asym_id_type: AsymIdType,
) !void {
    const owned_chains = try parseCommaSeparatedChains(allocator, chains_value);
    errdefer freeChains(allocator, owned_chains);
    try appendOwnedEntry(allocator, entries, filename_value, id_value, owned_chains, asym_id_type);
}

fn appendEntryFromChains(
    allocator: Allocator,
    entries: *std.ArrayListUnmanaged(Entry),
    filename_value: []const u8,
    id_value: ?[]const u8,
    chain_values: []const []const u8,
    asym_id_type: AsymIdType,
) !void {
    const chains = try dupeChains(allocator, chain_values);
    errdefer freeChains(allocator, chains);
    try appendOwnedEntry(allocator, entries, filename_value, id_value, chains, asym_id_type);
}

fn appendInterfaceEntryCsv(
    allocator: Allocator,
    entries: *std.ArrayListUnmanaged(InterfaceEntry),
    filename_value: []const u8,
    id_value: ?[]const u8,
    partner_a_value: []const u8,
    partner_b_value: []const u8,
    asym_id_type: AsymIdType,
) !void {
    const partner_a = try parseCommaSeparatedChains(allocator, partner_a_value);
    errdefer freeChains(allocator, partner_a);
    const partner_b = try parseCommaSeparatedChains(allocator, partner_b_value);
    errdefer freeChains(allocator, partner_b);
    try appendOwnedInterfaceEntry(allocator, entries, filename_value, id_value, partner_a, partner_b, asym_id_type);
}

fn appendInterfaceEntryJson(
    allocator: Allocator,
    entries: *std.ArrayListUnmanaged(InterfaceEntry),
    filename_value: []const u8,
    id_value: ?[]const u8,
    partner_a_values: []const []const u8,
    partner_b_values: []const []const u8,
    asym_id_type: AsymIdType,
) !void {
    const partner_a = try dupeChains(allocator, partner_a_values);
    errdefer freeChains(allocator, partner_a);
    const partner_b = try dupeChains(allocator, partner_b_values);
    errdefer freeChains(allocator, partner_b);
    try appendOwnedInterfaceEntry(allocator, entries, filename_value, id_value, partner_a, partner_b, asym_id_type);
}

fn appendOwnedEntry(
    allocator: Allocator,
    entries: *std.ArrayListUnmanaged(Entry),
    filename_value: []const u8,
    id_value: ?[]const u8,
    chains: []const []const u8,
    asym_id_type: AsymIdType,
) !void {
    const filename = try validateFilename(filename_value);
    const owned_filename = try allocator.dupe(u8, filename);
    errdefer allocator.free(owned_filename);
    const owned_id = if (id_value) |value| blk: {
        const id = std.mem.trim(u8, value, " \t\r");
        break :blk if (id.len == 0) null else try allocator.dupe(u8, try validateSelectionId(id));
    } else null;
    errdefer if (owned_id) |id| allocator.free(id);
    try entries.append(allocator, .{
        .filename = owned_filename,
        .id = owned_id,
        .chains = chains,
        .asym_id_type = asym_id_type,
    });
}

fn appendOwnedInterfaceEntry(
    allocator: Allocator,
    entries: *std.ArrayListUnmanaged(InterfaceEntry),
    filename_value: []const u8,
    id_value: ?[]const u8,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
    asym_id_type: AsymIdType,
) !void {
    const filename = try validateFilename(filename_value);
    const owned_filename = try allocator.dupe(u8, filename);
    errdefer allocator.free(owned_filename);
    const owned_id = if (id_value) |value| blk: {
        const id = std.mem.trim(u8, value, " \t\r");
        break :blk if (id.len == 0) null else try allocator.dupe(u8, try validateInterfaceId(id));
    } else null;
    errdefer if (owned_id) |id| allocator.free(id);
    try entries.append(allocator, .{
        .filename = owned_filename,
        .id = owned_id,
        .partner_a = partner_a,
        .partner_b = partner_b,
        .asym_id_type = asym_id_type,
    });
}

fn finish(allocator: Allocator, entries: *std.ArrayListUnmanaged(Entry)) !ChainMap {
    const owned_entries = try entries.toOwnedSlice(allocator);
    errdefer {
        for (owned_entries) |entry| entry.deinit(allocator);
        allocator.free(owned_entries);
    }

    var index = std.StringHashMapUnmanaged(std.ArrayListUnmanaged(usize)){};
    errdefer {
        var values = index.valueIterator();
        while (values.next()) |indices| indices.deinit(allocator);
        index.deinit(allocator);
    }
    for (owned_entries, 0..) |entry, i| {
        const result = try index.getOrPut(allocator, entry.filename);
        if (!result.found_existing) result.value_ptr.* = .empty;
        try result.value_ptr.append(allocator, i);
    }
    var values = index.valueIterator();
    while (values.next()) |indices| {
        if (indices.items.len <= 1) continue;
        const first_type = owned_entries[indices.items[0]].asym_id_type;
        for (indices.items) |entry_index| {
            const entry = owned_entries[entry_index];
            if (entry.id == null) return error.MissingSelectionId;
            if (entry.asym_id_type != first_type) return error.MixedSelectionAsymIdType;
        }
    }
    var ids = std.StringHashMapUnmanaged(void){};
    defer ids.deinit(allocator);
    for (owned_entries) |entry| {
        const id = entry.id orelse entry.filename;
        if (ids.contains(id)) return error.DuplicateSelectionId;
        try ids.put(allocator, id, {});
    }

    return .{
        .allocator = allocator,
        .entries = owned_entries,
        .index = index,
    };
}

fn finishInterface(allocator: Allocator, entries: *std.ArrayListUnmanaged(InterfaceEntry)) !InterfaceMap {
    const owned_entries = try entries.toOwnedSlice(allocator);
    errdefer {
        for (owned_entries) |entry| entry.deinit(allocator);
        allocator.free(owned_entries);
    }

    var index = std.StringHashMapUnmanaged(std.ArrayListUnmanaged(usize)){};
    errdefer {
        var values = index.valueIterator();
        while (values.next()) |indices| indices.deinit(allocator);
        index.deinit(allocator);
    }
    for (owned_entries, 0..) |entry, i| {
        const result = try index.getOrPut(allocator, entry.filename);
        if (!result.found_existing) result.value_ptr.* = .empty;
        try result.value_ptr.append(allocator, i);
    }
    var values = index.valueIterator();
    while (values.next()) |indices| {
        if (indices.items.len <= 1) continue;
        for (indices.items) |entry_index| {
            if (owned_entries[entry_index].id == null) return error.MissingInterfaceId;
        }
    }
    var ids = std.StringHashMapUnmanaged(void){};
    defer ids.deinit(allocator);
    for (owned_entries) |entry| {
        const id = entry.id orelse entry.filename;
        if (ids.contains(id)) return error.DuplicateInterfaceId;
        try ids.put(allocator, id, {});
    }

    return .{
        .allocator = allocator,
        .entries = owned_entries,
        .index = index,
    };
}

fn readFile(allocator: Allocator, io: std.Io, path: []const u8) ![]u8 {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var read_buf: [65536]u8 = undefined;
    var reader = file.reader(io, &read_buf);
    return reader.interface.allocRemaining(allocator, .unlimited);
}

fn validateFilename(filename_value: []const u8) ![]const u8 {
    const filename = std.mem.trim(u8, filename_value, " \t\r");
    if (filename.len == 0) return error.EmptyFilename;
    if (std.mem.findAny(u8, filename, "/\\") != null or std.mem.indexOf(u8, filename, "..") != null) {
        return error.UnsafeFilename;
    }
    return filename;
}

fn validateInterfaceId(id_value: []const u8) ![]const u8 {
    const id = std.mem.trim(u8, id_value, " \t\r");
    if (id.len == 0) return error.EmptyInterfaceId;
    return id;
}

fn validateSelectionId(id_value: []const u8) ![]const u8 {
    const id = std.mem.trim(u8, id_value, " \t\r");
    if (id.len == 0) return error.EmptySelectionId;
    return id;
}

fn parseCommaSeparatedChains(allocator: Allocator, value: []const u8) ![]const []const u8 {
    var chains = std.ArrayListUnmanaged([]const u8).empty;
    errdefer {
        for (chains.items) |chain| allocator.free(chain);
        chains.deinit(allocator);
    }

    var iter = std.mem.splitScalar(u8, value, ',');
    while (iter.next()) |raw_chain| {
        const chain = std.mem.trim(u8, raw_chain, " \t\r");
        if (chain.len == 0) return error.EmptyChainId;
        const owned_chain = try allocator.dupe(u8, chain);
        chains.append(allocator, owned_chain) catch |err| {
            allocator.free(owned_chain);
            return err;
        };
    }
    if (chains.items.len == 0) return error.EmptyChains;
    return chains.toOwnedSlice(allocator);
}

fn dupeChains(allocator: Allocator, values: []const []const u8) ![]const []const u8 {
    if (values.len == 0) return error.EmptyChains;
    const chains = try allocator.alloc([]const u8, values.len);
    var initialized: usize = 0;
    errdefer {
        for (chains[0..initialized]) |chain| allocator.free(chain);
        allocator.free(chains);
    }
    for (values, 0..) |value, i| {
        const chain = std.mem.trim(u8, value, " \t\r");
        if (chain.len == 0) return error.EmptyChainId;
        chains[i] = try allocator.dupe(u8, chain);
        initialized += 1;
    }
    return chains;
}

fn freeChains(allocator: Allocator, chains: []const []const u8) void {
    for (chains) |chain| allocator.free(chain);
    allocator.free(chains);
}

fn parseCsvRecord(allocator: Allocator, line: []const u8) ![]const []const u8 {
    var fields = std.ArrayListUnmanaged([]const u8).empty;
    errdefer freeFieldList(allocator, &fields);

    var i: usize = 0;
    while (true) {
        while (i < line.len and (line[i] == ' ' or line[i] == '\t')) : (i += 1) {}

        var value = std.ArrayListUnmanaged(u8).empty;
        errdefer value.deinit(allocator);
        if (i < line.len and line[i] == '"') {
            i += 1;
            var closed = false;
            while (i < line.len) {
                if (line[i] == '"') {
                    if (i + 1 < line.len and line[i + 1] == '"') {
                        try value.append(allocator, '"');
                        i += 2;
                        continue;
                    }
                    i += 1;
                    closed = true;
                    break;
                }
                try value.append(allocator, line[i]);
                i += 1;
            }
            if (!closed) return error.InvalidCsv;
            while (i < line.len and (line[i] == ' ' or line[i] == '\t')) : (i += 1) {}
            if (i < line.len and line[i] != ',') return error.InvalidCsv;
        } else {
            const start = i;
            while (i < line.len and line[i] != ',') : (i += 1) {}
            try value.appendSlice(allocator, std.mem.trim(u8, line[start..i], " \t"));
        }

        const owned_value = try value.toOwnedSlice(allocator);
        fields.append(allocator, owned_value) catch |err| {
            allocator.free(owned_value);
            return err;
        };
        if (i >= line.len) break;
        i += 1;
        if (i == line.len) {
            const empty = try allocator.alloc(u8, 0);
            fields.append(allocator, empty) catch |err| {
                allocator.free(empty);
                return err;
            };
            break;
        }
    }
    return fields.toOwnedSlice(allocator);
}

fn freeFields(allocator: Allocator, fields: []const []const u8) void {
    for (fields) |field| allocator.free(field);
    allocator.free(fields);
}

fn freeFieldList(allocator: Allocator, fields: *std.ArrayListUnmanaged([]const u8)) void {
    for (fields.items) |field| allocator.free(field);
    fields.deinit(allocator);
}

fn deinitEntryList(allocator: Allocator, entries: *std.ArrayListUnmanaged(Entry)) void {
    for (entries.items) |entry| entry.deinit(allocator);
    entries.deinit(allocator);
}

fn deinitInterfaceEntryList(allocator: Allocator, entries: *std.ArrayListUnmanaged(InterfaceEntry)) void {
    for (entries.items) |entry| entry.deinit(allocator);
    entries.deinit(allocator);
}

fn endsWithIgnoreCase(value: []const u8, suffix: []const u8) bool {
    if (value.len < suffix.len) return false;
    return std.ascii.eqlIgnoreCase(value[value.len - suffix.len ..], suffix);
}

test "parse CSV chain map with multi-chain auth selection" {
    var map = try parseCsv(
        std.testing.allocator,
        "filename,chains,asym_id_type\r\n" ++
            "1abc.cif,\"A,C\",auth\r\n" ++
            "2xyz.pdb,B,label\r\n",
    );
    defer map.deinit();

    const first = map.get("1abc.cif").?;
    try std.testing.expectEqual(@as(usize, 2), first.chains.len);
    try std.testing.expectEqualStrings("A", first.chains[0]);
    try std.testing.expectEqualStrings("C", first.chains[1]);
    try std.testing.expectEqual(AsymIdType.auth, first.asym_id_type);
    try std.testing.expectEqual(AsymIdType.label, map.get("2xyz.pdb").?.asym_id_type);
}

test "parse JSON chain map defaults to label IDs" {
    var map = try parseJson(std.testing.allocator,
        \\[
        \\  {"filename":"1abc.cif","chains":["A","C"]},
        \\  {"filename":"2xyz.cif","chains":["X"],"asym_id_type":"auth"}
        \\]
    );
    defer map.deinit();

    try std.testing.expectEqual(AsymIdType.label, map.get("1abc.cif").?.asym_id_type);
    try std.testing.expectEqual(AsymIdType.auth, map.get("2xyz.cif").?.asym_id_type);
}

test "parse multiple CSV selections for one filename with globally unique IDs" {
    var map = try parseCsv(
        std.testing.allocator,
        "filename,id,chains,asym_id_type\n" ++
            "1abc.cif,a,A,label\n" ++
            "1abc.cif,bc,\"B,C\",label\n" ++
            "2xyz.cif,other,X,auth\n",
    );
    defer map.deinit();

    const indices = map.getIndices("1abc.cif").?;
    try std.testing.expectEqual(@as(usize, 2), indices.len);
    try std.testing.expectEqualStrings("a", map.entries[indices[0]].id.?);
    try std.testing.expectEqualStrings("bc", map.entries[indices[1]].id.?);
    try std.testing.expectEqualStrings("B", map.entries[indices[1]].chains[0]);
    try std.testing.expectEqualStrings("C", map.entries[indices[1]].chains[1]);
}

test "parse multiple JSON selections for one filename with globally unique IDs" {
    var map = try parseJson(std.testing.allocator,
        \\[
        \\  {"filename":"1abc.cif","id":"a","chains":["A"]},
        \\  {"filename":"1abc.cif","id":"bc","chains":["B","C"]}
        \\]
    );
    defer map.deinit();

    const indices = map.getIndices("1abc.cif").?;
    try std.testing.expectEqual(@as(usize, 2), indices.len);
    try std.testing.expectEqualStrings("a", map.entries[indices[0]].id.?);
    try std.testing.expectEqualStrings("bc", map.entries[indices[1]].id.?);
}

test "require IDs when a chain map filename appears more than once" {
    try std.testing.expectError(
        error.MissingSelectionId,
        parseCsv(
            std.testing.allocator,
            "filename,id,chains\n" ++
                "1abc.cif,a,A\n" ++
                "1abc.cif,,B\n",
        ),
    );
}

test "reject duplicate selection IDs across filenames" {
    try std.testing.expectError(
        error.DuplicateSelectionId,
        parseCsv(
            std.testing.allocator,
            "filename,id,chains\n" ++
                "1abc.cif,same,A\n" ++
                "2xyz.cif,same,B\n",
        ),
    );
}

test "reject mixed asym ID types within one selection-map filename" {
    try std.testing.expectError(
        error.MixedSelectionAsymIdType,
        parseCsv(
            std.testing.allocator,
            "filename,id,chains,asym_id_type\n" ++
                "1abc.cif,a,A,label\n" ++
                "1abc.cif,b,B,auth\n",
        ),
    );
}

test "reject unquoted multi-chain CSV field" {
    try std.testing.expectError(
        error.InvalidCsv,
        parseCsv(
            std.testing.allocator,
            "filename,chains\n" ++
                "1abc.cif,A,C\n",
        ),
    );
}

test "parse CSV interface map with multi-chain partners" {
    var map = try parseInterfaceCsv(
        std.testing.allocator,
        "filename,partner_a,partner_b,asym_id_type\n" ++
            "1abc.cif,\"A,B\",\"C,D\",auth\n",
    );
    defer map.deinit();

    const entry = map.entries[map.getIndices("1abc.cif").?[0]];
    try std.testing.expectEqual(@as(usize, 2), entry.partner_a.len);
    try std.testing.expectEqualStrings("A", entry.partner_a[0]);
    try std.testing.expectEqualStrings("B", entry.partner_a[1]);
    try std.testing.expectEqual(@as(usize, 2), entry.partner_b.len);
    try std.testing.expectEqualStrings("C", entry.partner_b[0]);
    try std.testing.expectEqualStrings("D", entry.partner_b[1]);
    try std.testing.expectEqual(AsymIdType.auth, entry.asym_id_type);
}

test "parse JSON interface map" {
    var map = try parseInterfaceJson(std.testing.allocator,
        \\[
        \\  {
        \\    "filename":"1abc.cif",
        \\    "partner_a":["A","B"],
        \\    "partner_b":["C","D"],
        \\    "asym_id_type":"label"
        \\  }
        \\]
    );
    defer map.deinit();

    const entry = map.entries[map.getIndices("1abc.cif").?[0]];
    try std.testing.expectEqualStrings("B", entry.partner_a[1]);
    try std.testing.expectEqualStrings("D", entry.partner_b[1]);
    try std.testing.expectEqual(AsymIdType.label, entry.asym_id_type);
}

test "parse multiple CSV interfaces for one filename with stable IDs" {
    var map = try parseInterfaceCsv(
        std.testing.allocator,
        "filename,id,partner_a,partner_b,asym_id_type\n" ++
            "1abc.cif,interaction-001,A,\"B,C\",label\n" ++
            "1abc.cif,interaction-002,D,\"B,C\",label\n",
    );
    defer map.deinit();

    const indices = map.getIndices("1abc.cif").?;
    try std.testing.expectEqual(@as(usize, 2), indices.len);
    try std.testing.expectEqualStrings("interaction-001", map.entries[indices[0]].id.?);
    try std.testing.expectEqualStrings("interaction-002", map.entries[indices[1]].id.?);
}

test "parse multiple JSON interfaces for one filename with stable IDs" {
    var map = try parseInterfaceJson(std.testing.allocator,
        \\[
        \\  {"filename":"1abc.cif","id":"one","partner_a":["A"],"partner_b":["B"]},
        \\  {"filename":"1abc.cif","id":"two","partner_a":["C"],"partner_b":["D"]}
        \\]
    );
    defer map.deinit();

    const indices = map.getIndices("1abc.cif").?;
    try std.testing.expectEqual(@as(usize, 2), indices.len);
    try std.testing.expectEqualStrings("one", map.entries[indices[0]].id.?);
    try std.testing.expectEqualStrings("two", map.entries[indices[1]].id.?);
}

test "require IDs when a filename has multiple interfaces" {
    try std.testing.expectError(
        error.MissingInterfaceId,
        parseInterfaceCsv(
            std.testing.allocator,
            "filename,id,partner_a,partner_b\n" ++
                "1abc.cif,one,A,B\n" ++
                "1abc.cif,,C,D\n",
        ),
    );
}

test "reject duplicate interface IDs" {
    try std.testing.expectError(
        error.DuplicateInterfaceId,
        parseInterfaceCsv(
            std.testing.allocator,
            "filename,id,partner_a,partner_b\n" ++
                "1abc.cif,same,A,B\n" ++
                "2xyz.cif,same,C,D\n",
        ),
    );
}

test "reject explicit interface ID that collides with legacy filename ID" {
    try std.testing.expectError(
        error.DuplicateInterfaceId,
        parseInterfaceCsv(
            std.testing.allocator,
            "filename,id,partner_a,partner_b\n" ++
                "1abc.cif,,A,B\n" ++
                "2xyz.cif,1abc.cif,C,D\n",
        ),
    );
}

const std = @import("std");
const Allocator = std.mem.Allocator;
const toml_parser = @import("toml_parser.zig");

pub const ManifestError = error{
    UnsupportedVersion,
    MissingJobName,
    DuplicateJobName,
    UnsafeJobName,
    InvalidFieldType,
    NoJobs,
};

pub const Error = ManifestError || Allocator.Error || toml_parser.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError;

pub const Globals = struct {
    input_dir: ?[]const u8 = null,
    output_dir: ?[]const u8 = null,
    algorithm: ?[]const u8 = null,
    classifier: ?[]const u8 = null,
    ccd: ?[]const u8 = null,
    sdf: ?[]const []const u8 = null,
    threads: ?usize = null,
    probe_radius: ?f64 = null,
    n_points: ?u32 = null,
    n_slices: ?u32 = null,
    precision: ?[]const u8 = null,
    format: ?[]const u8 = null,
    include_hydrogens: ?bool = null,
    include_hetatm: ?bool = null,
    use_bitmask: ?bool = null,
    timing: ?bool = null,
    quiet: ?bool = null,
    auth_chain: ?bool = null,
};

pub const Job = struct {
    name: []const u8,
    chains: ?[]const []const u8 = null,
    auth_chain: ?bool = null,
};

pub const Manifest = struct {
    allocator: Allocator,
    content: []const u8,
    globals: Globals,
    jobs: []Job,

    pub fn deinit(self: *Manifest) void {
        if (self.globals.sdf) |items| self.allocator.free(items);
        for (self.jobs) |job| {
            if (job.chains) |chains| self.allocator.free(chains);
        }
        self.allocator.free(self.jobs);
        self.allocator.free(self.content);
        self.* = undefined;
    }
};

pub fn parse(allocator: Allocator, content: []const u8) Error!Manifest {
    const owned_content = try allocator.dupe(u8, content);
    errdefer allocator.free(owned_content);
    return parseOwned(allocator, owned_content);
}

pub fn parseFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Manifest {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var read_buf: [65536]u8 = undefined;
    var reader = file.reader(io, &read_buf);
    const content = try reader.interface.allocRemaining(allocator, .unlimited);
    errdefer allocator.free(content);
    return parseOwned(allocator, content);
}

fn parseOwned(allocator: Allocator, owned_content: []const u8) Error!Manifest {
    var doc = try toml_parser.parse(allocator, owned_content);
    defer doc.deinit();

    const root = doc.getTable("") orelse return error.UnsupportedVersion;
    try validateVersion(root);

    const globals = try parseGlobals(allocator, root);
    errdefer if (globals.sdf) |items| allocator.free(items);

    var jobs = std.ArrayListUnmanaged(Job).empty;
    errdefer {
        for (jobs.items) |job| {
            if (job.chains) |chains| allocator.free(chains);
        }
        jobs.deinit(allocator);
    }

    for (doc.array_tables) |table| {
        if (!std.mem.eql(u8, table.name, "jobs")) continue;
        {
            const job = try parseJob(allocator, table.entries, jobs.items);
            errdefer if (job.chains) |chains| allocator.free(chains);
            try jobs.append(allocator, job);
        }
    }

    if (jobs.items.len == 0) return error.NoJobs;

    return .{
        .allocator = allocator,
        .content = owned_content,
        .globals = globals,
        .jobs = try jobs.toOwnedSlice(allocator),
    };
}

fn validateVersion(root: toml_parser.Table) ManifestError!void {
    const value = findValue(root.entries, "version") orelse return error.UnsupportedVersion;
    switch (value) {
        .integer => |version| if (version == 1) return else return error.UnsupportedVersion,
        else => return error.UnsupportedVersion,
    }
}

fn parseGlobals(allocator: Allocator, root: toml_parser.Table) Error!Globals {
    var globals = Globals{};
    errdefer if (globals.sdf) |items| allocator.free(items);

    globals.input_dir = try optionalString(root.entries, "input_dir");
    globals.output_dir = try optionalString(root.entries, "output_dir");
    globals.algorithm = try optionalString(root.entries, "algorithm");
    globals.classifier = try optionalString(root.entries, "classifier");
    globals.ccd = try optionalString(root.entries, "ccd");
    globals.sdf = try optionalStringOrStringArray(allocator, root.entries, "sdf");
    globals.threads = try optionalUsize(root.entries, "threads");
    globals.probe_radius = try optionalFloat(root.entries, "probe_radius");
    globals.n_points = try optionalU32(root.entries, "n_points");
    globals.n_slices = try optionalU32(root.entries, "n_slices");
    globals.precision = try optionalString(root.entries, "precision");
    globals.format = try optionalString(root.entries, "format");
    globals.include_hydrogens = try optionalBool(root.entries, "include_hydrogens");
    globals.include_hetatm = try optionalBool(root.entries, "include_hetatm");
    globals.use_bitmask = try optionalBool(root.entries, "use_bitmask");
    globals.timing = try optionalBool(root.entries, "timing");
    globals.quiet = try optionalBool(root.entries, "quiet");
    globals.auth_chain = try optionalBool(root.entries, "auth_chain");
    return globals;
}

fn parseJob(allocator: Allocator, entries: []const toml_parser.Value.Entry, existing_jobs: []const Job) Error!Job {
    const name = try optionalString(entries, "name") orelse return error.MissingJobName;
    if (name.len == 0) return error.MissingJobName;
    if (!isSafeJobName(name)) return error.UnsafeJobName;
    for (existing_jobs) |job| {
        if (std.mem.eql(u8, job.name, name)) return error.DuplicateJobName;
    }

    const chains = try optionalStringArray(allocator, entries, "chains");
    errdefer if (chains) |items| allocator.free(items);
    const auth_chain = try optionalBool(entries, "auth_chain");

    return .{
        .name = name,
        .chains = chains,
        .auth_chain = auth_chain,
    };
}

fn isSafeJobName(name: []const u8) bool {
    return std.mem.indexOfScalar(u8, name, '/') == null and
        std.mem.indexOfScalar(u8, name, '\\') == null and
        std.mem.indexOf(u8, name, "..") == null;
}

fn findValue(entries: []const toml_parser.Value.Entry, key: []const u8) ?toml_parser.Value {
    for (entries) |entry| {
        if (std.mem.eql(u8, entry.key, key)) return entry.value;
    }
    return null;
}

fn optionalString(entries: []const toml_parser.Value.Entry, key: []const u8) ManifestError!?[]const u8 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .string => |s| s,
        else => error.InvalidFieldType,
    };
}

fn optionalBool(entries: []const toml_parser.Value.Entry, key: []const u8) ManifestError!?bool {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .boolean => |b| b,
        else => error.InvalidFieldType,
    };
}

fn optionalFloat(entries: []const toml_parser.Value.Entry, key: []const u8) ManifestError!?f64 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .float => |f| f,
        .integer => |i| @floatFromInt(i),
        else => error.InvalidFieldType,
    };
}

fn optionalUsize(entries: []const toml_parser.Value.Entry, key: []const u8) ManifestError!?usize {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .integer => |i| if (i >= 0 and i <= std.math.maxInt(usize)) @intCast(i) else error.InvalidFieldType,
        else => error.InvalidFieldType,
    };
}

fn optionalU32(entries: []const toml_parser.Value.Entry, key: []const u8) ManifestError!?u32 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .integer => |i| if (i >= 0 and i <= std.math.maxInt(u32)) @intCast(i) else error.InvalidFieldType,
        else => error.InvalidFieldType,
    };
}

fn optionalStringArray(allocator: Allocator, entries: []const toml_parser.Value.Entry, key: []const u8) Error!?[]const []const u8 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .string_array => |items| try allocator.dupe([]const u8, items),
        else => error.InvalidFieldType,
    };
}

fn optionalStringOrStringArray(allocator: Allocator, entries: []const toml_parser.Value.Entry, key: []const u8) Error!?[]const []const u8 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .string => |s| blk: {
            const items = try allocator.alloc([]const u8, 1);
            items[0] = s;
            break :blk items;
        },
        .string_array => |items| try allocator.dupe([]const u8, items),
        else => error.InvalidFieldType,
    };
}

test "parse A B AB manifest" {
    const input =
        \\version = 1
        \\
        \\input_dir = "structures"
        \\output_dir = "results"
        \\
        \\format = "jsonl"
        \\use_bitmask = true
        \\n_points = 128
        \\classifier = "ccd"
        \\
        \\[[jobs]]
        \\name = "chain_A"
        \\chains = ["A"]
        \\
        \\[[jobs]]
        \\name = "chain_B"
        \\chains = ["B"]
        \\
        \\[[jobs]]
        \\name = "complex_AB"
        \\chains = ["A", "B"]
    ;

    var manifest = try parse(std.testing.allocator, input);
    defer manifest.deinit();

    try std.testing.expectEqualStrings("structures", manifest.globals.input_dir.?);
    try std.testing.expectEqualStrings("results", manifest.globals.output_dir.?);
    try std.testing.expectEqualStrings("jsonl", manifest.globals.format.?);
    try std.testing.expectEqual(true, manifest.globals.use_bitmask.?);
    try std.testing.expectEqual(@as(u32, 128), manifest.globals.n_points.?);
    try std.testing.expectEqualStrings("ccd", manifest.globals.classifier.?);
    try std.testing.expectEqual(@as(usize, 3), manifest.jobs.len);
    try std.testing.expectEqualStrings("chain_A", manifest.jobs[0].name);
    try std.testing.expectEqualStrings("A", manifest.jobs[0].chains.?[0]);
    try std.testing.expectEqualStrings("chain_B", manifest.jobs[1].name);
    try std.testing.expectEqualStrings("B", manifest.jobs[1].chains.?[0]);
    try std.testing.expectEqualStrings("complex_AB", manifest.jobs[2].name);
    try std.testing.expectEqualStrings("A", manifest.jobs[2].chains.?[0]);
    try std.testing.expectEqualStrings("B", manifest.jobs[2].chains.?[1]);
}

test "reject duplicate and unsafe job names" {
    const duplicate =
        \\version = 1
        \\[[jobs]]
        \\name = "chain_A"
        \\[[jobs]]
        \\name = "chain_A"
    ;
    try std.testing.expectError(error.DuplicateJobName, parse(std.testing.allocator, duplicate));

    const unsafe_slash =
        \\version = 1
        \\[[jobs]]
        \\name = "bad/name"
    ;
    try std.testing.expectError(error.UnsafeJobName, parse(std.testing.allocator, unsafe_slash));

    const unsafe_parent =
        \\version = 1
        \\[[jobs]]
        \\name = ".."
    ;
    try std.testing.expectError(error.UnsafeJobName, parse(std.testing.allocator, unsafe_parent));
}

test "parse sdf string and array globals" {
    const single =
        \\version = 1
        \\sdf = "lig.sdf"
        \\[[jobs]]
        \\name = "all"
    ;
    var single_manifest = try parse(std.testing.allocator, single);
    defer single_manifest.deinit();
    try std.testing.expectEqual(@as(usize, 1), single_manifest.globals.sdf.?.len);
    try std.testing.expectEqualStrings("lig.sdf", single_manifest.globals.sdf.?[0]);

    const array =
        \\version = 1
        \\sdf = ["a.sdf", "b.sdf"]
        \\[[jobs]]
        \\name = "all"
    ;
    var array_manifest = try parse(std.testing.allocator, array);
    defer array_manifest.deinit();
    try std.testing.expectEqual(@as(usize, 2), array_manifest.globals.sdf.?.len);
    try std.testing.expectEqualStrings("a.sdf", array_manifest.globals.sdf.?[0]);
    try std.testing.expectEqualStrings("b.sdf", array_manifest.globals.sdf.?[1]);
}

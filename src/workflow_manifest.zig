const std = @import("std");
const Allocator = std.mem.Allocator;
const toml_parser = @import("toml_parser.zig");

pub const WorkflowError = error{
    UnsupportedVersion,
    InvalidKind,
    MissingJobName,
    DuplicateJobName,
    UnsafeJobName,
    InvalidFieldType,
    InvalidClassifierConfig,
    UnknownField,
    NoJobs,
};

pub const Error = WorkflowError || Allocator.Error || toml_parser.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError;

pub const Input = struct {
    path: ?[]const u8 = null,
    dir: ?[]const u8 = null,
    chain: ?[]const u8 = null,
    model: ?u32 = null,
    mol: ?[]const u8 = null,
};

pub const Output = struct {
    path: ?[]const u8 = null,
    dir: ?[]const u8 = null,
    format: ?[]const u8 = null,
};

pub const Calculation = struct {
    algorithm: ?[]const u8 = null,
    threads: ?usize = null,
    probe_radius: ?f64 = null,
    n_points: ?u32 = null,
    n_slices: ?u32 = null,
    precision: ?[]const u8 = null,
    include_hydrogens: ?bool = null,
    include_hetatm: ?bool = null,
    use_bitmask: ?bool = null,
    timing: ?bool = null,
    quiet: ?bool = null,
    auth_chain: ?bool = null,
    residue_map: ?bool = null,
    per_residue: ?bool = null,
    rsa: ?bool = null,
    polar: ?bool = null,
    validate_only: ?bool = null,
};

pub const ClassifierConfig = struct {
    type: ?[]const u8 = null,
    config: ?[]const u8 = null,
    ccd: ?[]const u8 = null,
    sdf: ?[]const []const u8 = null,
};

pub const Job = struct {
    name: []const u8,
    chains: ?[]const []const u8 = null,
    auth_chain: ?bool = null,
};

pub const Workflow = struct {
    allocator: Allocator,
    content: []const u8,
    input: Input = .{},
    output: Output = .{},
    calculation: Calculation = .{},
    classifier: ClassifierConfig = .{},
    jobs: []Job = &.{},
    is_legacy_batch_manifest: bool = false,

    pub fn deinit(self: *Workflow) void {
        if (self.classifier.sdf) |items| self.allocator.free(items);
        for (self.jobs) |job| {
            if (job.chains) |chains| self.allocator.free(chains);
        }
        self.allocator.free(self.jobs);
        self.allocator.free(self.content);
        self.* = undefined;
    }
};

pub fn parse(allocator: Allocator, content: []const u8) Error!Workflow {
    const owned_content = try allocator.dupe(u8, content);
    return parseOwned(allocator, owned_content);
}

pub fn parseFile(allocator: Allocator, io: std.Io, path: []const u8) Error!Workflow {
    const file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var read_buf: [65536]u8 = undefined;
    var reader = file.reader(io, &read_buf);
    const content = try reader.interface.allocRemaining(allocator, .unlimited);
    return parseOwned(allocator, content);
}

fn parseOwned(allocator: Allocator, owned_content: []const u8) Error!Workflow {
    var workflow_initialized = false;
    errdefer if (!workflow_initialized) allocator.free(owned_content);

    var doc = try toml_parser.parse(allocator, owned_content);
    defer doc.deinit();

    const root = doc.getTable("") orelse return error.UnsupportedVersion;
    try validateVersion(root);
    try validateKind(root);

    const is_legacy = isLegacyBatchManifest(root);
    try validateKnownDocumentShape(doc, is_legacy);

    var workflow = Workflow{
        .allocator = allocator,
        .content = owned_content,
        .is_legacy_batch_manifest = is_legacy,
    };
    workflow_initialized = true;
    errdefer workflow.deinit();

    if (is_legacy) {
        try parseLegacyRootIntoWorkflow(allocator, root, &workflow);
    } else {
        if (doc.getTable("input")) |table| workflow.input = try parseInput(table);
        if (doc.getTable("output")) |table| workflow.output = try parseOutput(table);
        if (doc.getTable("calculation")) |table| workflow.calculation = try parseCalculation(table);
        if (doc.getTable("classifier")) |table| workflow.classifier = try parseClassifier(allocator, table);
    }

    workflow.jobs = try parseJobs(allocator, doc.array_tables);
    try validateClassifier(workflow.classifier);
    return workflow;
}

fn validateVersion(root: toml_parser.Table) WorkflowError!void {
    const value = findValue(root.entries, "version") orelse return error.UnsupportedVersion;
    switch (value) {
        .integer => |version| if (version == 1) return else return error.UnsupportedVersion,
        else => return error.UnsupportedVersion,
    }
}

fn validateKind(root: toml_parser.Table) WorkflowError!void {
    const value = findValue(root.entries, "kind") orelse return;
    switch (value) {
        .string => |kind| if (std.mem.eql(u8, kind, "workflow")) return else return error.InvalidKind,
        else => return error.InvalidKind,
    }
}

fn isLegacyBatchManifest(root: toml_parser.Table) bool {
    const legacy_fields = [_][]const u8{
        "input_dir",   "output_dir", "classifier",        "format",
        "algorithm",   "threads",    "probe_radius",      "n_points",
        "n_slices",    "precision",  "include_hydrogens", "include_hetatm",
        "use_bitmask", "timing",     "quiet",             "auth_chain",
        "residue_map", "ccd",        "sdf",
    };
    for (legacy_fields) |field| {
        if (findValue(root.entries, field) != null) return true;
    }
    return false;
}

fn parseLegacyRootIntoWorkflow(allocator: Allocator, root: toml_parser.Table, workflow: *Workflow) Error!void {
    try rejectUnknownFields(root.entries, &.{
        "version",     "kind",         "input_dir", "output_dir", "algorithm",   "classifier", "ccd",               "sdf",
        "threads",     "probe_radius", "n_points",  "n_slices",   "precision",   "format",     "include_hydrogens", "include_hetatm",
        "use_bitmask", "timing",       "quiet",     "auth_chain", "residue_map",
    });

    workflow.input = .{ .dir = try optionalString(root.entries, "input_dir") };
    workflow.output = .{
        .dir = try optionalString(root.entries, "output_dir"),
        .format = try optionalString(root.entries, "format"),
    };
    workflow.calculation = .{
        .algorithm = try optionalString(root.entries, "algorithm"),
        .threads = try optionalUsize(root.entries, "threads"),
        .probe_radius = try optionalFloat(root.entries, "probe_radius"),
        .n_points = try optionalU32(root.entries, "n_points"),
        .n_slices = try optionalU32(root.entries, "n_slices"),
        .precision = try optionalString(root.entries, "precision"),
        .include_hydrogens = try optionalBool(root.entries, "include_hydrogens"),
        .include_hetatm = try optionalBool(root.entries, "include_hetatm"),
        .use_bitmask = try optionalBool(root.entries, "use_bitmask"),
        .timing = try optionalBool(root.entries, "timing"),
        .quiet = try optionalBool(root.entries, "quiet"),
        .auth_chain = try optionalBool(root.entries, "auth_chain"),
        .residue_map = try optionalBool(root.entries, "residue_map"),
    };
    workflow.classifier = .{
        .type = try optionalString(root.entries, "classifier"),
        .ccd = try optionalString(root.entries, "ccd"),
        .sdf = try optionalStringOrStringArray(allocator, root.entries, "sdf"),
    };
}

fn parseInput(table: toml_parser.Table) WorkflowError!Input {
    try rejectUnknownFields(table.entries, &.{ "path", "dir", "chain", "model", "mol" });
    return .{
        .path = try optionalString(table.entries, "path"),
        .dir = try optionalString(table.entries, "dir"),
        .chain = try optionalString(table.entries, "chain"),
        .model = try optionalU32(table.entries, "model"),
        .mol = try optionalString(table.entries, "mol"),
    };
}

fn parseOutput(table: toml_parser.Table) WorkflowError!Output {
    try rejectUnknownFields(table.entries, &.{ "path", "dir", "format" });
    return .{
        .path = try optionalString(table.entries, "path"),
        .dir = try optionalString(table.entries, "dir"),
        .format = try optionalString(table.entries, "format"),
    };
}

fn parseCalculation(table: toml_parser.Table) WorkflowError!Calculation {
    try rejectUnknownFields(table.entries, &.{
        "algorithm",         "threads",        "probe_radius", "n_points", "n_slices",      "precision",
        "include_hydrogens", "include_hetatm", "use_bitmask",  "timing",   "quiet",         "auth_chain",
        "residue_map",       "per_residue",    "rsa",          "polar",    "validate_only",
    });
    return .{
        .algorithm = try optionalString(table.entries, "algorithm"),
        .threads = try optionalUsize(table.entries, "threads"),
        .probe_radius = try optionalFloat(table.entries, "probe_radius"),
        .n_points = try optionalU32(table.entries, "n_points"),
        .n_slices = try optionalU32(table.entries, "n_slices"),
        .precision = try optionalString(table.entries, "precision"),
        .include_hydrogens = try optionalBool(table.entries, "include_hydrogens"),
        .include_hetatm = try optionalBool(table.entries, "include_hetatm"),
        .use_bitmask = try optionalBool(table.entries, "use_bitmask"),
        .timing = try optionalBool(table.entries, "timing"),
        .quiet = try optionalBool(table.entries, "quiet"),
        .auth_chain = try optionalBool(table.entries, "auth_chain"),
        .residue_map = try optionalBool(table.entries, "residue_map"),
        .per_residue = try optionalBool(table.entries, "per_residue"),
        .rsa = try optionalBool(table.entries, "rsa"),
        .polar = try optionalBool(table.entries, "polar"),
        .validate_only = try optionalBool(table.entries, "validate_only"),
    };
}

fn parseClassifier(allocator: Allocator, table: toml_parser.Table) Error!ClassifierConfig {
    try rejectUnknownFields(table.entries, &.{ "type", "config", "ccd", "sdf" });
    return .{
        .type = try optionalString(table.entries, "type"),
        .config = try optionalString(table.entries, "config"),
        .ccd = try optionalString(table.entries, "ccd"),
        .sdf = try optionalStringOrStringArray(allocator, table.entries, "sdf"),
    };
}

fn validateClassifier(config: ClassifierConfig) WorkflowError!void {
    const classifier_type = config.type orelse return;
    if (std.mem.eql(u8, classifier_type, "custom")) {
        if (config.config == null) return error.InvalidClassifierConfig;
        if (config.ccd != null or config.sdf != null) return error.InvalidClassifierConfig;
        return;
    }
    if (config.config != null) return error.InvalidClassifierConfig;
}

fn parseJobs(allocator: Allocator, array_tables: []const toml_parser.Document.ArrayTable) Error![]Job {
    var jobs = std.ArrayListUnmanaged(Job).empty;
    errdefer {
        for (jobs.items) |job| {
            if (job.chains) |chains| allocator.free(chains);
        }
        jobs.deinit(allocator);
    }

    for (array_tables) |table| {
        if (!std.mem.eql(u8, table.name, "jobs")) return error.UnknownField;
        {
            const job = try parseJob(allocator, table.entries, jobs.items);
            errdefer if (job.chains) |chains| allocator.free(chains);
            try jobs.append(allocator, job);
        }
    }

    return jobs.toOwnedSlice(allocator);
}

fn parseJob(allocator: Allocator, entries: []const toml_parser.Value.Entry, existing_jobs: []const Job) Error!Job {
    try rejectUnknownFields(entries, &.{ "name", "chains", "auth_chain" });

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

fn validateKnownDocumentShape(doc: toml_parser.Document, is_legacy: bool) WorkflowError!void {
    for (doc.tables) |table| {
        if (table.name.len == 0) continue;
        if (is_legacy) return error.UnknownField;
        if (std.mem.eql(u8, table.name, "input") or
            std.mem.eql(u8, table.name, "output") or
            std.mem.eql(u8, table.name, "calculation") or
            std.mem.eql(u8, table.name, "classifier"))
        {
            continue;
        }
        return error.UnknownField;
    }

    if (!is_legacy) {
        const root = doc.getTable("") orelse return error.UnsupportedVersion;
        try rejectUnknownFields(root.entries, &.{ "version", "kind" });
    }
}

fn rejectUnknownFields(entries: []const toml_parser.Value.Entry, allowed: []const []const u8) WorkflowError!void {
    for (entries) |entry| {
        if (!isAllowedField(entry.key, allowed)) return error.UnknownField;
    }
}

fn isAllowedField(key: []const u8, allowed: []const []const u8) bool {
    for (allowed) |allowed_key| {
        if (std.mem.eql(u8, key, allowed_key)) return true;
    }
    return false;
}

fn findValue(entries: []const toml_parser.Value.Entry, key: []const u8) ?toml_parser.Value {
    for (entries) |entry| {
        if (std.mem.eql(u8, entry.key, key)) return entry.value;
    }
    return null;
}

fn optionalString(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?[]const u8 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .string => |s| s,
        else => error.InvalidFieldType,
    };
}

fn optionalBool(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?bool {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .boolean => |b| b,
        else => error.InvalidFieldType,
    };
}

fn optionalFloat(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?f64 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .float => |f| f,
        .integer => |i| @floatFromInt(i),
        else => error.InvalidFieldType,
    };
}

fn optionalUsize(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?usize {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .integer => |i| if (i >= 0 and i <= std.math.maxInt(usize)) @intCast(i) else error.InvalidFieldType,
        else => error.InvalidFieldType,
    };
}

fn optionalU32(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?u32 {
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

fn isSafeJobName(name: []const u8) bool {
    return std.mem.indexOfScalar(u8, name, '/') == null and
        std.mem.indexOfScalar(u8, name, '\\') == null and
        std.mem.indexOf(u8, name, "..") == null;
}

test "parse sectioned calc workflow" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\path = "structure.cif"
        \\chain = "A,B"
        \\model = 1
        \\mol = "LIG"
        \\
        \\[output]
        \\path = "result.json"
        \\format = "json"
        \\
        \\[calculation]
        \\algorithm = "sr"
        \\probe_radius = 1.4
        \\n_points = 128
        \\n_slices = 20
        \\precision = "f64"
        \\use_bitmask = true
        \\include_hydrogens = false
        \\include_hetatm = true
        \\auth_chain = true
        \\per_residue = true
        \\rsa = true
        \\polar = true
        \\validate_only = false
        \\
        \\[classifier]
        \\type = "custom"
        \\config = "my_radii.toml"
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqual(false, workflow.is_legacy_batch_manifest);
    try std.testing.expectEqualStrings("structure.cif", workflow.input.path.?);
    try std.testing.expectEqualStrings("A,B", workflow.input.chain.?);
    try std.testing.expectEqual(@as(u32, 1), workflow.input.model.?);
    try std.testing.expectEqualStrings("LIG", workflow.input.mol.?);
    try std.testing.expectEqualStrings("result.json", workflow.output.path.?);
    try std.testing.expectEqualStrings("json", workflow.output.format.?);
    try std.testing.expectEqualStrings("sr", workflow.calculation.algorithm.?);
    try std.testing.expectEqual(@as(f64, 1.4), workflow.calculation.probe_radius.?);
    try std.testing.expectEqual(@as(u32, 128), workflow.calculation.n_points.?);
    try std.testing.expectEqual(@as(u32, 20), workflow.calculation.n_slices.?);
    try std.testing.expectEqualStrings("f64", workflow.calculation.precision.?);
    try std.testing.expectEqual(true, workflow.calculation.use_bitmask.?);
    try std.testing.expectEqual(false, workflow.calculation.include_hydrogens.?);
    try std.testing.expectEqual(true, workflow.calculation.include_hetatm.?);
    try std.testing.expectEqual(true, workflow.calculation.auth_chain.?);
    try std.testing.expectEqual(true, workflow.calculation.per_residue.?);
    try std.testing.expectEqual(true, workflow.calculation.rsa.?);
    try std.testing.expectEqual(true, workflow.calculation.polar.?);
    try std.testing.expectEqual(false, workflow.calculation.validate_only.?);
    try std.testing.expectEqualStrings("custom", workflow.classifier.type.?);
    try std.testing.expectEqualStrings("my_radii.toml", workflow.classifier.config.?);
    try std.testing.expectEqual(@as(usize, 0), workflow.jobs.len);
}

test "parse sectioned batch workflow with jobs" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "structures"
        \\
        \\[output]
        \\dir = "results"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\threads = 8
        \\n_points = 128
        \\residue_map = true
        \\
        \\[classifier]
        \\type = "ccd"
        \\ccd = "components.zsdc"
        \\sdf = ["ligand.sdf", "cofactor.sdf"]
        \\
        \\[[jobs]]
        \\name = "chain_A"
        \\chains = ["A"]
        \\
        \\[[jobs]]
        \\name = "complex_AB"
        \\chains = ["A", "B"]
        \\auth_chain = true
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqualStrings("structures", workflow.input.dir.?);
    try std.testing.expectEqualStrings("results", workflow.output.dir.?);
    try std.testing.expectEqualStrings("jsonl", workflow.output.format.?);
    try std.testing.expectEqual(@as(usize, 8), workflow.calculation.threads.?);
    try std.testing.expectEqual(true, workflow.calculation.residue_map.?);
    try std.testing.expectEqualStrings("ccd", workflow.classifier.type.?);
    try std.testing.expectEqualStrings("components.zsdc", workflow.classifier.ccd.?);
    try std.testing.expectEqual(@as(usize, 2), workflow.classifier.sdf.?.len);
    try std.testing.expectEqualStrings("ligand.sdf", workflow.classifier.sdf.?[0]);
    try std.testing.expectEqual(@as(usize, 2), workflow.jobs.len);
    try std.testing.expectEqualStrings("chain_A", workflow.jobs[0].name);
    try std.testing.expectEqualStrings("A", workflow.jobs[0].chains.?[0]);
    try std.testing.expectEqualStrings("complex_AB", workflow.jobs[1].name);
    try std.testing.expectEqualStrings("B", workflow.jobs[1].chains.?[1]);
    try std.testing.expectEqual(true, workflow.jobs[1].auth_chain.?);
}

test "parse legacy flat batch manifest" {
    const input =
        \\version = 1
        \\input_dir = "structures"
        \\output_dir = "results"
        \\format = "jsonl"
        \\classifier = "ccd"
        \\ccd = "components.zsdc"
        \\sdf = "ligand.sdf"
        \\threads = 4
        \\n_points = 128
        \\use_bitmask = true
        \\
        \\[[jobs]]
        \\name = "chain_A"
        \\chains = ["A"]
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqual(true, workflow.is_legacy_batch_manifest);
    try std.testing.expectEqualStrings("structures", workflow.input.dir.?);
    try std.testing.expectEqualStrings("results", workflow.output.dir.?);
    try std.testing.expectEqualStrings("jsonl", workflow.output.format.?);
    try std.testing.expectEqualStrings("ccd", workflow.classifier.type.?);
    try std.testing.expectEqualStrings("components.zsdc", workflow.classifier.ccd.?);
    try std.testing.expectEqual(@as(usize, 1), workflow.classifier.sdf.?.len);
    try std.testing.expectEqualStrings("ligand.sdf", workflow.classifier.sdf.?[0]);
    try std.testing.expectEqual(@as(usize, 4), workflow.calculation.threads.?);
    try std.testing.expectEqual(@as(u32, 128), workflow.calculation.n_points.?);
    try std.testing.expectEqual(true, workflow.calculation.use_bitmask.?);
    try std.testing.expectEqual(@as(usize, 1), workflow.jobs.len);
    try std.testing.expectEqualStrings("chain_A", workflow.jobs[0].name);
    try std.testing.expectEqualStrings("A", workflow.jobs[0].chains.?[0]);
}

test "reject invalid classifier combinations" {
    const custom_with_ccd =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "custom"
        \\config = "my_radii.toml"
        \\ccd = "components.zsdc"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, custom_with_ccd));

    const custom_without_config =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "custom"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, custom_without_config));

    const ccd_with_config =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "ccd"
        \\config = "my_radii.toml"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, ccd_with_config));
}

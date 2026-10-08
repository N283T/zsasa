const std = @import("std");
const builtin = @import("builtin");
const Allocator = std.mem.Allocator;
const toml_parser = @import("toml_parser.zig");
const altloc = @import("altloc.zig");

pub const WorkflowError = error{
    UnsupportedVersion,
    InvalidKind,
    MissingJobName,
    DuplicateJobName,
    UnsafeJobName,
    InvalidFieldType,
    InvalidClassifierConfig,
    InvalidAnalysisConfig,
    UnknownField,
    NoJobs,
    EmptyJobChains,
    AnalysisWithJobs,
};

/// Explanation of a workflow error whose name does not say what to change,
/// or null. Printed by the commands below the error name.
pub fn errorHint(err: anyerror) ?[]const u8 {
    return switch (err) {
        error.EmptyJobChains => "a [[jobs]] entry has an empty chains array; " ++
            "list at least one chain ID, or remove the key to select every chain",
        error.AnalysisWithJobs => "[analysis] and [[jobs]] cannot be combined: an [analysis] " ++
            "workflow runs one BSA analysis and has no jobs; keep one of them and move " ++
            "the other into a separate workflow file",
        else => null,
    };
}

pub const Error = WorkflowError || Allocator.Error || toml_parser.Error || std.Io.File.OpenError || std.Io.File.ReadStreamingError || error{ ReadFailed, StreamTooLong };

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
    jsonl: JsonlOutput = .{},
};

pub const JsonlOutput = struct {
    atom_areas: ?bool = null,
    atom_identity: ?bool = null,
    total_area: ?bool = null,
    decimals: ?u8 = null,
    metadata: ?[]const u8 = null,
};

pub const Calculation = struct {
    algorithm: ?[]const u8 = null,
    threads: ?usize = null,
    probe_radius: ?f64 = null,
    n_points: ?u32 = null,
    n_slices: ?u32 = null,
    /// Lee-Richards arc angles: "exact" or "fast"
    lr_trig: ?[]const u8 = null,
    precision: ?[]const u8 = null,
    include_hydrogens: ?bool = null,
    include_hetatm: ?bool = null,
    use_bitmask: ?bool = null,
    timing: ?bool = null,
    quiet: ?bool = null,
    auth_chain: ?bool = null,
    /// `altloc`: "auto", "none", "all", "highest-occupancy" or one altLoc ID,
    /// as for the `--altloc` flag.
    altloc: ?altloc.AltLocSetting = null,
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

pub const Analysis = struct {
    type: ?[]const u8 = null,
    name: ?[]const u8 = null,
    partner_a: ?[]const []const u8 = null,
    partner_b: ?[]const []const u8 = null,
    chain_map: ?[]const u8 = null,
    level: ?[]const u8 = null,
    atom_output: ?bool = null,
};

pub const Job = struct {
    name: []const u8,
    chains: ?[]const []const u8 = null,
    chain_map: ?[]const u8 = null,
    auth_chain: ?bool = null,
};

pub const Workflow = struct {
    allocator: Allocator,
    content: []const u8,
    input: Input = .{},
    output: Output = .{},
    calculation: Calculation = .{},
    classifier: ClassifierConfig = .{},
    analysis: ?Analysis = null,
    jobs: []Job = &.{},
    is_legacy_batch_workflow: bool = false,

    pub fn deinit(self: *Workflow) void {
        if (self.classifier.sdf) |items| self.allocator.free(items);
        if (self.analysis) |analysis| {
            if (analysis.partner_a) |items| self.allocator.free(items);
            if (analysis.partner_b) |items| self.allocator.free(items);
        }
        for (self.jobs) |job| {
            if (job.chains) |chains| self.allocator.free(chains);
        }
        self.allocator.free(self.jobs);
        self.allocator.free(self.content);
        self.* = undefined;
    }

    /// Make `[input] chain` the chain selection of every job that has none of
    /// its own (neither `chains` nor `chain_map`), as batch workflows read it.
    /// The chain IDs are comma separated, as for `--chain`. `checkKeys`
    /// rejects a blank value and a value no job would use.
    pub fn applyInputChainToJobs(self: *Workflow) Error!void {
        const chain = self.input.chain orelse return;
        for (self.jobs) |*job| {
            if (job.chains != null or job.chain_map != null) continue;
            var ids = std.ArrayListUnmanaged([]const u8).empty;
            errdefer ids.deinit(self.allocator);
            var parts = std.mem.splitScalar(u8, chain, ',');
            while (parts.next()) |part| {
                const id = std.mem.trim(u8, part, " ");
                if (id.len > 0) try ids.append(self.allocator, id);
            }
            if (ids.items.len == 0) return error.InvalidFieldType;
            job.chains = try ids.toOwnedSlice(self.allocator);
        }
    }
};

/// The command a manifest is run by and, for batch, the kind of workflow.
pub const Mode = enum {
    /// `zsasa calc --workflow`: one structure.
    calc,
    /// `zsasa batch --workflow` with `[[jobs]]`.
    batch_jobs,
    /// `zsasa batch --workflow` with `[analysis]`.
    batch_analysis,
};

pub const Severity = enum {
    /// The key would change the results, or the manifest does not fit the
    /// command: the command stops.
    err,
    /// The key only affects reporting or where output goes: the command goes on.
    warning,
};

pub const Finding = struct {
    severity: Severity,
    /// Names the key and says what to do about it.
    message: []const u8,
};

/// Keys a command does not read. Nothing a user writes is dropped silently:
/// each of them is an error or a warning.
pub const Findings = struct {
    items: [max_findings]Finding = undefined,
    len: usize = 0,

    /// More than any manifest can produce (`checkKeys` has fewer checks).
    pub const max_findings = 32;

    pub fn add(self: *Findings, severity: Severity, message: []const u8) void {
        std.debug.assert(self.len < max_findings);
        self.items[self.len] = .{ .severity = severity, .message = message };
        self.len += 1;
    }

    pub fn slice(self: *const Findings) []const Finding {
        return self.items[0..self.len];
    }

    pub fn errorCount(self: *const Findings) usize {
        var count: usize = 0;
        for (self.slice()) |finding| {
            if (finding.severity == .err) count += 1;
        }
        return count;
    }

    /// Print every finding to stderr (also in quiet mode: they are not
    /// progress output) and fail if any is an error. Prints nothing in tests.
    pub fn report(self: *const Findings) error{InvalidArgument}!void {
        for (self.slice()) |finding| {
            if (builtin.is_test) continue;
            const label = switch (finding.severity) {
                .err => "Error",
                .warning => "Warning",
            };
            std.debug.print("{s}: {s}\n", .{ label, finding.message });
        }
        if (self.errorCount() > 0) return error.InvalidArgument;
    }
};

const batch_hint = "run it with 'zsasa batch --workflow'";
const calc_hint = "run it on one structure with 'zsasa calc --workflow'";

/// The keys of `workflow` that the command for `mode` does not read: errors
/// where honoring them would change the results, warnings where they only
/// affect reporting or the output location. A boolean key set to false asks
/// for what the command does anyway and is not reported.
pub fn checkKeys(workflow: Workflow, mode: Mode) Findings {
    var findings = Findings{};
    switch (mode) {
        .calc => checkCalcKeys(workflow, &findings),
        .batch_jobs, .batch_analysis => checkBatchKeys(workflow, mode, &findings),
    }
    return findings;
}

fn checkCalcKeys(workflow: Workflow, findings: *Findings) void {
    const legacy = workflow.is_legacy_batch_workflow;
    if (workflow.analysis != null) {
        findings.add(.err, "[analysis] is read only by 'zsasa batch --workflow' and calc would ignore it: " ++
            "this is a batch manifest; " ++ batch_hint);
    }
    if (workflow.jobs.len > 0) {
        findings.add(.err, "[[jobs]] is read only by 'zsasa batch --workflow' and calc would ignore it: " ++
            "this is a batch manifest; " ++ batch_hint);
    }
    if (workflow.input.dir != null) {
        findings.add(.err, if (legacy)
            "input_dir is read only by 'zsasa batch --workflow' and calc would ignore it: " ++
                "this is a batch manifest; " ++ batch_hint ++ ", or name one structure with [input] path"
        else
            "[input] dir is read only by 'zsasa batch --workflow' and calc would ignore it: " ++
                "this is a batch manifest; " ++ batch_hint ++ ", or name one structure with [input] path");
    }
    if (workflow.calculation.residue_map orelse false) {
        findings.add(.err, "[calculation] residue_map = true is read only by 'zsasa batch --workflow' " ++
            "(calc writes no residue map): remove it, or " ++ batch_hint);
    }
    if (workflow.output.dir != null) {
        findings.add(.warning, if (legacy)
            "output_dir is read only by 'zsasa batch --workflow' and calc ignores it: " ++
                "the result goes to [output] path or the output argument"
        else
            "[output] dir is read only by 'zsasa batch --workflow' and calc ignores it: " ++
                "the result goes to [output] path or the output argument");
    }
    const jsonl = workflow.output.jsonl;
    if (jsonl.atom_areas != null or jsonl.atom_identity != null or jsonl.total_area != null or
        jsonl.decimals != null or jsonl.metadata != null)
    {
        findings.add(.warning, "[output.jsonl] is read only by 'zsasa batch --workflow' and calc ignores it " ++
            "(calc writes no JSONL): remove the table");
    }
}

fn checkBatchKeys(workflow: Workflow, mode: Mode, findings: *Findings) void {
    if (workflow.input.path != null) {
        findings.add(.err, "[input] path is not read by 'zsasa batch --workflow' (batch processes a directory): " ++
            "set [input] dir, or " ++ calc_hint);
    }
    if (workflow.input.model != null) {
        findings.add(.err, "[input] model is not supported by batch (it cannot select a model): " ++
            "remove it, or " ++ calc_hint);
    }
    if (workflow.input.mol != null) {
        findings.add(.err, "[input] mol is not supported by batch (every molecule of an SDF file is processed): " ++
            "remove it, or " ++ calc_hint);
    }
    const calculation = workflow.calculation;
    if (calculation.rsa orelse false) {
        findings.add(.err, "[calculation] rsa = true is not supported by batch (only calc writes RSA tables): " ++
            "remove it, or " ++ calc_hint);
    }
    if (calculation.per_residue orelse false) {
        findings.add(.err, "[calculation] per_residue = true is not supported by batch " ++
            "(only calc writes per-residue output; batch has residue_map with JSONL output): " ++
            "remove it, or " ++ calc_hint);
    }
    if (calculation.polar orelse false) {
        findings.add(.err, "[calculation] polar = true is not supported by batch " ++
            "(only calc writes the polar/apolar split): remove it, or " ++ calc_hint);
    }
    if (calculation.validate_only orelse false) {
        findings.add(.err, "[calculation] validate_only = true is not supported by batch " ++
            "(batch would run the calculation): remove it, or " ++ calc_hint);
    }
    if (workflow.output.path != null) {
        findings.add(.warning, "[output] path is not read by 'zsasa batch --workflow' (batch writes into a directory): " ++
            "the output goes to [output] dir or the output argument");
    }

    if (mode == .batch_analysis) {
        if (workflow.input.chain != null) {
            findings.add(.err, "[input] chain is not supported by an [analysis] workflow: " ++
                "select the chains of the interface with partner_a and partner_b");
        }
        if (calculation.residue_map orelse false) {
            findings.add(.err, "[calculation] residue_map = true is not supported by an [analysis] workflow " ++
                "(its rows carry no residue map): remove it");
        }
        if (workflow.output.jsonl.atom_identity orelse false) {
            findings.add(.err, "[output.jsonl] atom_identity = true is not supported by an [analysis] workflow " ++
                "(its rows carry no atom list; atom_identity belongs to chain_map jobs): remove it");
        }
    } else if (workflow.input.chain) |chain| {
        checkInputChainOfJobs(workflow, chain, findings);
    }
}

/// `[input] chain` is the default chain selection of the jobs of a batch
/// workflow (see `Workflow.applyInputChainToJobs`).
fn checkInputChainOfJobs(workflow: Workflow, chain: []const u8, findings: *Findings) void {
    var has_id = false;
    var parts = std.mem.splitScalar(u8, chain, ',');
    while (parts.next()) |part| {
        if (std.mem.trim(u8, part, " ").len > 0) has_id = true;
    }
    if (!has_id) {
        findings.add(.err, "[input] chain needs at least one chain ID (for example \"A\" or \"A,B\"): " ++
            "list the chains or remove the key");
        return;
    }
    var uses_default = false;
    for (workflow.jobs) |job| {
        if (job.chain_map != null) {
            findings.add(.err, "[input] chain cannot be combined with a job that sets chain_map " ++
                "(the map selects the chains per file): remove [input] chain, or give the other jobs their own chains");
            return;
        }
        if (job.chains == null) uses_default = true;
    }
    if (!uses_default and workflow.jobs.len > 0) {
        findings.add(.err, "[input] chain is used by no job because every [[jobs]] entry sets its own chains: " ++
            "remove [input] chain");
    }
}

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

    try rejectDuplicateTableHeaders(owned_content);

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
        .is_legacy_batch_workflow = is_legacy,
    };
    workflow_initialized = true;
    errdefer workflow.deinit();

    if (is_legacy) {
        try parseLegacyRootIntoWorkflow(allocator, root, &workflow);
    } else {
        if (doc.getTable("input")) |table| workflow.input = try parseInput(table);
        if (doc.getTable("output")) |table| workflow.output = try parseOutput(table);
        if (doc.getTable("output.jsonl")) |table| workflow.output.jsonl = try parseOutputJsonl(table);
        if (doc.getTable("calculation")) |table| workflow.calculation = try parseCalculation(table);
        if (doc.getTable("classifier")) |table| workflow.classifier = try parseClassifier(allocator, table);
        if (doc.getTable("analysis")) |table| workflow.analysis = try parseAnalysis(allocator, table);
    }

    workflow.jobs = try parseJobs(allocator, doc.array_tables);
    // The TOML parser drops an array table without keys, so a bare `[[jobs]]`
    // is only visible in the text.
    if (try countJobHeaders(owned_content) != workflow.jobs.len) return error.MissingJobName;
    if (workflow.analysis != null and workflow.jobs.len > 0) return error.AnalysisWithJobs;
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

fn parseOutputJsonl(table: toml_parser.Table) WorkflowError!JsonlOutput {
    try rejectUnknownFields(table.entries, &.{ "atom_areas", "atom_identity", "total_area", "decimals", "metadata" });
    const output = JsonlOutput{
        .atom_areas = try optionalBool(table.entries, "atom_areas"),
        .atom_identity = try optionalBool(table.entries, "atom_identity"),
        .total_area = try optionalBool(table.entries, "total_area"),
        .decimals = try optionalU8(table.entries, "decimals"),
        .metadata = try optionalString(table.entries, "metadata"),
    };
    try validateOutputJsonl(output);
    return output;
}

fn validateOutputJsonl(output: JsonlOutput) WorkflowError!void {
    if (output.decimals) |decimals| {
        if (decimals > 15) return error.InvalidFieldType;
    }
    if (output.metadata) |metadata| {
        if (!std.mem.eql(u8, metadata, "none") and !std.mem.eql(u8, metadata, "sidecar")) {
            return error.InvalidFieldType;
        }
    }
}

fn parseCalculation(table: toml_parser.Table) WorkflowError!Calculation {
    try rejectUnknownFields(table.entries, &.{
        "algorithm",         "threads",        "probe_radius", "n_points", "n_slices",      "precision",
        "include_hydrogens", "include_hetatm", "use_bitmask",  "timing",   "quiet",         "auth_chain",
        "residue_map",       "per_residue",    "rsa",          "polar",    "validate_only", "altloc",
        "lr_trig",
    });
    const lr_trig = try optionalString(table.entries, "lr_trig");
    if (lr_trig) |value| {
        if (!std.mem.eql(u8, value, "exact") and !std.mem.eql(u8, value, "fast")) {
            return error.InvalidFieldType;
        }
    }
    return .{
        .algorithm = try optionalString(table.entries, "algorithm"),
        .threads = try optionalUsize(table.entries, "threads"),
        .probe_radius = try optionalFloat(table.entries, "probe_radius"),
        .n_points = try optionalU32(table.entries, "n_points"),
        .n_slices = try optionalU32(table.entries, "n_slices"),
        .lr_trig = lr_trig,
        .precision = try optionalString(table.entries, "precision"),
        .include_hydrogens = try optionalBool(table.entries, "include_hydrogens"),
        .include_hetatm = try optionalBool(table.entries, "include_hetatm"),
        .use_bitmask = try optionalBool(table.entries, "use_bitmask"),
        .timing = try optionalBool(table.entries, "timing"),
        .quiet = try optionalBool(table.entries, "quiet"),
        .auth_chain = try optionalBool(table.entries, "auth_chain"),
        .altloc = try optionalAltLoc(table.entries, "altloc"),
        .residue_map = try optionalBool(table.entries, "residue_map"),
        .per_residue = try optionalBool(table.entries, "per_residue"),
        .rsa = try optionalBool(table.entries, "rsa"),
        .polar = try optionalBool(table.entries, "polar"),
        .validate_only = try optionalBool(table.entries, "validate_only"),
    };
}

/// An altLoc setting written as for the `--altloc` flag. Any other string
/// is an error, so that a typo does not fall back to the default.
fn optionalAltLoc(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?altloc.AltLocSetting {
    const value = try optionalString(entries, key) orelse return null;
    return altloc.parseSetting(value) orelse error.InvalidFieldType;
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

fn parseAnalysis(allocator: Allocator, table: toml_parser.Table) Error!Analysis {
    try rejectUnknownFields(table.entries, &.{ "type", "name", "partner_a", "partner_b", "chain_map", "level", "atom_output" });
    const analysis = Analysis{
        .type = try optionalString(table.entries, "type"),
        .name = try optionalString(table.entries, "name"),
        .partner_a = try optionalStringArray(allocator, table.entries, "partner_a"),
        .partner_b = try optionalStringArray(allocator, table.entries, "partner_b"),
        .chain_map = try optionalString(table.entries, "chain_map"),
        .level = try optionalString(table.entries, "level"),
        .atom_output = try optionalBool(table.entries, "atom_output"),
    };
    errdefer {
        if (analysis.partner_a) |items| allocator.free(items);
        if (analysis.partner_b) |items| allocator.free(items);
    }
    try validateAnalysis(analysis);
    return analysis;
}

fn validateClassifier(config: ClassifierConfig) WorkflowError!void {
    const classifier_type = config.type orelse return;
    if (std.mem.eql(u8, classifier_type, "custom")) {
        if (config.config == null) return error.InvalidClassifierConfig;
        if (config.ccd != null or config.sdf != null) return error.InvalidClassifierConfig;
        return;
    }
    if (config.config != null) return error.InvalidClassifierConfig;
    if (!classifierAllowsCcdResources(classifier_type) and (config.ccd != null or config.sdf != null)) {
        return error.InvalidClassifierConfig;
    }
}

fn classifierAllowsCcdResources(classifier_type: []const u8) bool {
    return std.mem.eql(u8, classifier_type, "ccd");
}

fn validateAnalysis(analysis: Analysis) WorkflowError!void {
    const analysis_type = analysis.type orelse return error.InvalidAnalysisConfig;
    if (!std.mem.eql(u8, analysis_type, "bsa")) return error.InvalidAnalysisConfig;
    if (analysis.name) |name| {
        if (name.len == 0 or !isSafeJobName(name)) return error.InvalidAnalysisConfig;
    }
    if (analysis.chain_map) |path| {
        if (path.len == 0 or analysis.partner_a != null or analysis.partner_b != null) {
            return error.InvalidAnalysisConfig;
        }
    } else if (analysis.partner_a == null or analysis.partner_b == null) {
        return error.InvalidAnalysisConfig;
    }
    if (analysis.chain_map == null) {
        const partner_a = analysis.partner_a.?;
        const partner_b = analysis.partner_b.?;
        if (partner_a.len == 0 or partner_b.len == 0) return error.InvalidAnalysisConfig;
        for (partner_a) |chain| {
            if (chain.len == 0) return error.InvalidAnalysisConfig;
        }
        for (partner_b) |chain| {
            if (chain.len == 0) return error.InvalidAnalysisConfig;
        }
    }
    if (analysis.level) |level| {
        if (!std.mem.eql(u8, level, "total") and !std.mem.eql(u8, level, "residue")) {
            return error.InvalidAnalysisConfig;
        }
    }
    if ((analysis.atom_output orelse false) and !std.mem.eql(u8, analysis.level orelse "total", "residue")) {
        return error.InvalidAnalysisConfig;
    }
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
    try rejectUnknownFields(entries, &.{ "name", "chains", "chain_map", "auth_chain" });

    const name = try optionalString(entries, "name") orelse return error.MissingJobName;
    if (name.len == 0) return error.MissingJobName;
    if (!isSafeJobName(name)) return error.UnsafeJobName;
    for (existing_jobs) |job| {
        if (std.mem.eql(u8, job.name, name)) return error.DuplicateJobName;
    }

    const chains = try optionalStringArray(allocator, entries, "chains");
    errdefer if (chains) |items| allocator.free(items);
    // An empty list is neither "every chain" (that is the absent key) nor a
    // selection: the runners would disagree on what it means.
    if (chains) |items| {
        if (items.len == 0) return error.EmptyJobChains;
    }
    const chain_map = try optionalString(entries, "chain_map");
    if (chains != null and chain_map != null) return error.InvalidFieldType;
    const auth_chain = try optionalBool(entries, "auth_chain");
    if (chain_map != null and auth_chain != null) return error.InvalidFieldType;

    return .{
        .name = name,
        .chains = chains,
        .chain_map = chain_map,
        .auth_chain = auth_chain,
    };
}

/// The number of `[[jobs]]` headers in `content`. Any other array-of-tables
/// header is an unknown field.
fn countJobHeaders(content: []const u8) WorkflowError!usize {
    var count: usize = 0;
    var lines = std.mem.splitScalar(u8, content, '\n');
    while (lines.next()) |raw_line| {
        const line = std.mem.trim(u8, toml_parser.stripComment(raw_line), " \t\r");
        if (!std.mem.startsWith(u8, line, "[[")) continue;
        const end = std.mem.find(u8, line[2..], "]]") orelse continue;
        if (!std.mem.eql(u8, std.mem.trim(u8, line[2 .. 2 + end], " \t"), "jobs")) return error.UnknownField;
        count += 1;
    }
    return count;
}

fn rejectDuplicateTableHeaders(content: []const u8) WorkflowError!void {
    var start: usize = 0;
    while (start <= content.len) {
        const end = std.mem.indexOfScalarPos(u8, content, start, '\n') orelse content.len;
        const raw_line = content[start..end];
        if (tableHeaderName(raw_line)) |name| {
            if (hasEarlierTableHeader(content[0..start], name)) return error.UnknownField;
        }
        if (end == content.len) break;
        start = end + 1;
    }
}

fn hasEarlierTableHeader(content: []const u8, name: []const u8) bool {
    var start: usize = 0;
    while (start <= content.len) {
        const end = std.mem.indexOfScalarPos(u8, content, start, '\n') orelse content.len;
        const raw_line = content[start..end];
        if (tableHeaderName(raw_line)) |earlier_name| {
            if (std.mem.eql(u8, earlier_name, name)) return true;
        }
        if (end == content.len) break;
        start = end + 1;
    }
    return false;
}

fn tableHeaderName(raw_line: []const u8) ?[]const u8 {
    const line_without_comment = if (std.mem.indexOfScalar(u8, raw_line, '#')) |comment_start|
        raw_line[0..comment_start]
    else
        raw_line;
    const line = std.mem.trim(u8, line_without_comment, " \t\r");
    if (line.len == 0 or !std.mem.startsWith(u8, line, "[") or std.mem.startsWith(u8, line, "[[")) return null;
    const end = std.mem.indexOfScalar(u8, line, ']') orelse return null;
    const after = std.mem.trim(u8, line[end + 1 ..], " \t\r");
    if (after.len != 0) return null;
    return std.mem.trim(u8, line[1..end], " \t");
}

fn validateKnownDocumentShape(doc: toml_parser.Document, is_legacy: bool) WorkflowError!void {
    for (doc.tables, 0..) |table, index| {
        if (hasEarlierTableName(doc.tables[0..index], table.name)) return error.UnknownField;
        if (table.name.len == 0) continue;
        if (is_legacy) return error.UnknownField;
        if (std.mem.eql(u8, table.name, "input") or
            std.mem.eql(u8, table.name, "output") or
            std.mem.eql(u8, table.name, "output.jsonl") or
            std.mem.eql(u8, table.name, "calculation") or
            std.mem.eql(u8, table.name, "classifier") or
            std.mem.eql(u8, table.name, "analysis"))
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
    for (entries, 0..) |entry, index| {
        if (!isAllowedField(entry.key, allowed)) return error.UnknownField;
        if (hasEarlierKey(entries[0..index], entry.key)) return error.UnknownField;
    }
}

fn hasEarlierTableName(tables: []const toml_parser.Table, name: []const u8) bool {
    for (tables) |table| {
        if (std.mem.eql(u8, table.name, name)) return true;
    }
    return false;
}

fn hasEarlierKey(entries: []const toml_parser.Value.Entry, key: []const u8) bool {
    for (entries) |entry| {
        if (std.mem.eql(u8, entry.key, key)) return true;
    }
    return false;
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

fn optionalU8(entries: []const toml_parser.Value.Entry, key: []const u8) WorkflowError!?u8 {
    const value = findValue(entries, key) orelse return null;
    return switch (value) {
        .integer => |i| if (i >= 0 and i <= std.math.maxInt(u8)) @intCast(i) else error.InvalidFieldType,
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
        \\lr_trig = "fast"
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

    try std.testing.expectEqual(false, workflow.is_legacy_batch_workflow);
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
    try std.testing.expectEqualStrings("fast", workflow.calculation.lr_trig.?);
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

test "parse workflow job with per-file chain map" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "structures"
        \\
        \\[[jobs]]
        \\name = "selected_complexes"
        \\chain_map = "chains.csv"
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqual(@as(usize, 1), workflow.jobs.len);
    try std.testing.expectEqualStrings("chains.csv", workflow.jobs[0].chain_map.?);
    try std.testing.expect(workflow.jobs[0].chains == null);
}

test "reject workflow job with an empty chains array" {
    const header =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[[jobs]]
        \\name = "nothing"
        \\
    ;
    try std.testing.expectError(error.EmptyJobChains, parse(std.testing.allocator, header ++ "chains = []\n"));
    try std.testing.expectError(error.EmptyJobChains, parse(std.testing.allocator, header ++ "chains = [ ]\nauth_chain = true\n"));
    // An earlier job with chains must be released when a later job is rejected.
    try std.testing.expectError(error.EmptyJobChains, parse(
        std.testing.allocator,
        header ++ "chains = [\"A\"]\n\n[[jobs]]\nname = \"second\"\nchains = []\n",
    ));
    try std.testing.expect(errorHint(error.EmptyJobChains) != null);
    try std.testing.expect(errorHint(error.MissingJobName) == null);

    // Without the key the job selects every chain.
    var workflow = try parse(std.testing.allocator, header);
    defer workflow.deinit();
    try std.testing.expect(workflow.jobs[0].chains == null);
}

test "reject workflow job with both chains and chain map" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[[jobs]]
        \\name = "ambiguous"
        \\chains = ["A"]
        \\chain_map = "chains.csv"
    ;
    try std.testing.expectError(error.InvalidFieldType, parse(std.testing.allocator, input));
}

test "parse output jsonl workflow options" {
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
        \\[output.jsonl]
        \\atom_areas = false
        \\atom_identity = true
        \\total_area = true
        \\decimals = 3
        \\metadata = "sidecar"
        \\
        \\[[jobs]]
        \\name = "all"
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqual(false, workflow.output.jsonl.atom_areas.?);
    try std.testing.expectEqual(true, workflow.output.jsonl.atom_identity.?);
    try std.testing.expectEqual(true, workflow.output.jsonl.total_area.?);
    try std.testing.expectEqual(@as(u8, 3), workflow.output.jsonl.decimals.?);
    try std.testing.expectEqualStrings("sidecar", workflow.output.jsonl.metadata.?);
}

test "reject unknown output jsonl workflow option" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[output.jsonl]
        \\unexpected = true
    ;

    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, input));
}

test "parse BSA analysis workflow" {
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
        \\[analysis]
        \\type = "bsa"
        \\name = "interface_ab"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\level = "residue"
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqualStrings("bsa", workflow.analysis.?.type.?);
    try std.testing.expectEqualStrings("interface_ab", workflow.analysis.?.name.?);
    try std.testing.expectEqual(@as(usize, 1), workflow.analysis.?.partner_a.?.len);
    try std.testing.expectEqualStrings("A", workflow.analysis.?.partner_a.?[0]);
    try std.testing.expectEqual(@as(usize, 1), workflow.analysis.?.partner_b.?.len);
    try std.testing.expectEqualStrings("B", workflow.analysis.?.partner_b.?[0]);
    try std.testing.expectEqualStrings("residue", workflow.analysis.?.level.?);
}

test "parse BSA analysis workflow with per-file chain map" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "interfaces"
        \\chain_map = "interfaces.csv"
        \\level = "total"
    ;
    var workflow = try parse(std.testing.allocator, input);
    defer workflow.deinit();

    try std.testing.expectEqualStrings("interfaces.csv", workflow.analysis.?.chain_map.?);
    try std.testing.expect(workflow.analysis.?.partner_a == null);
    try std.testing.expect(workflow.analysis.?.partner_b == null);
}

test "reject BSA analysis with fixed partners and chain map" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\chain_map = "interfaces.csv"
    ;
    try std.testing.expectError(error.InvalidAnalysisConfig, parse(std.testing.allocator, input));
}

test "reject BSA analysis without both partners" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
    ;
    try std.testing.expectError(error.InvalidAnalysisConfig, parse(std.testing.allocator, input));
}

test "reject BSA analysis invalid level" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\level = "chain"
    ;
    try std.testing.expectError(error.InvalidAnalysisConfig, parse(std.testing.allocator, input));
}

test "parse BSA analysis with opt-in atom output" {
    const content =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\level = "residue"
        \\atom_output = true
    ;
    var workflow = try parse(std.testing.allocator, content);
    defer workflow.deinit();

    try std.testing.expectEqual(true, workflow.analysis.?.atom_output.?);
}

test "reject BSA atom output without residue detail" {
    const content =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\atom_output = true
    ;
    try std.testing.expectError(error.InvalidAnalysisConfig, parse(std.testing.allocator, content));
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

    try std.testing.expectEqual(true, workflow.is_legacy_batch_workflow);
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

test "calculation lr_trig accepts exact and fast, is optional, and rejects anything else" {
    const allocator = std.testing.allocator;
    const header =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[calculation]
        \\algorithm = "lr"
        \\
    ;

    {
        var workflow = try parse(allocator, header ++ "lr_trig = \"exact\"\n");
        defer workflow.deinit();
        try std.testing.expectEqualStrings("exact", workflow.calculation.lr_trig.?);
    }
    {
        var workflow = try parse(allocator, header ++ "lr_trig = \"fast\"\n");
        defer workflow.deinit();
        try std.testing.expectEqualStrings("fast", workflow.calculation.lr_trig.?);
    }
    {
        var workflow = try parse(allocator, header);
        defer workflow.deinit();
        try std.testing.expect(workflow.calculation.lr_trig == null);
    }

    try std.testing.expectError(error.InvalidFieldType, parse(allocator, header ++ "lr_trig = \"approximate\"\n"));
    try std.testing.expectError(error.InvalidFieldType, parse(allocator, header ++ "lr_trig = \"\"\n"));
    try std.testing.expectError(error.InvalidFieldType, parse(allocator, header ++ "lr_trig = true\n"));
}

test "lr_trig is a [calculation] key only: the legacy root form rejects it" {
    const input =
        \\version = 1
        \\input_dir = "structures"
        \\algorithm = "lr"
        \\lr_trig = "fast"
    ;
    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, input));
}

test "reject duplicate section names" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[calculation]
        \\n_points = 128
        \\
        \\[calculation]
        \\n_slices = 20
    ;

    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, input));
}

test "reject trailing empty duplicate section" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\path = "a.cif"
        \\
        \\[input]
    ;

    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, input));
}

test "reject leading empty duplicate section" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\
        \\[input]
        \\path = "a.cif"
    ;

    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, input));
}

test "reject duplicate keys in section" {
    const input =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\path = "first.cif"
        \\path = "second.cif"
    ;

    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, input));
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

    const naccess_with_ccd =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "naccess"
        \\ccd = "components.zsdc"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, naccess_with_ccd));

    const naccess_with_sdf =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "naccess"
        \\sdf = "ligand.sdf"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, naccess_with_sdf));

    const protor_with_ccd =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "protor"
        \\ccd = "components.zsdc"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, protor_with_ccd));

    const protor_with_sdf =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[classifier]
        \\type = "protor"
        \\sdf = "ligand.sdf"
    ;
    try std.testing.expectError(error.InvalidClassifierConfig, parse(std.testing.allocator, protor_with_sdf));
}

test "parse workflow calculation altloc" {
    const allocator = std.testing.allocator;
    const Case = struct { value: []const u8, mode: altloc.AltLocMode, id: u8 = 'A' };
    const cases = [_]Case{
        .{ .value = "auto", .mode = .auto },
        .{ .value = "none", .mode = .none },
        .{ .value = "all", .mode = .all },
        .{ .value = "highest-occupancy", .mode = .highest_occupancy },
        .{ .value = "B", .mode = .selected, .id = 'B' },
    };
    for (cases) |case| {
        const content = try std.fmt.allocPrint(allocator,
            \\version = 1
            \\
            \\[calculation]
            \\altloc = "{s}"
            \\
        , .{case.value});
        defer allocator.free(content);

        var workflow = try parse(allocator, content);
        defer workflow.deinit();
        try std.testing.expectEqual(case.mode, workflow.calculation.altloc.?.mode);
        try std.testing.expectEqual(case.id, workflow.calculation.altloc.?.id);
    }

    var without = try parse(allocator,
        \\version = 1
        \\
        \\[calculation]
        \\n_points = 8
        \\
    );
    defer without.deinit();
    try std.testing.expect(without.calculation.altloc == null);
}

test "parse workflow rejects an invalid calculation altloc" {
    const allocator = std.testing.allocator;
    // A typo must not fall back to the default
    try std.testing.expectError(error.InvalidFieldType, parse(allocator,
        \\version = 1
        \\
        \\[calculation]
        \\altloc = "first"
        \\
    ));
    try std.testing.expectError(error.InvalidFieldType, parse(allocator,
        \\version = 1
        \\
        \\[calculation]
        \\altloc = true
        \\
    ));
}

test "parse rejects an [analysis] workflow that also has [[jobs]]" {
    const allocator = std.testing.allocator;
    const analysis =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\
    ;
    // Either order, with a chains array that has to be released
    try std.testing.expectError(error.AnalysisWithJobs, parse(allocator, analysis ++ "\n[[jobs]]\nname = \"j\"\nchains = [\"A\"]\n"));
    try std.testing.expectError(error.AnalysisWithJobs, parse(
        allocator,
        "version = 1\n\n[[jobs]]\nname = \"j\"\nchains = [\"A\"]\n\n[analysis]\ntype = \"bsa\"\npartner_a = [\"A\"]\npartner_b = [\"B\"]\n",
    ));
    try std.testing.expect(errorHint(error.AnalysisWithJobs) != null);

    // Each alone is fine
    var only_analysis = try parse(allocator, analysis);
    defer only_analysis.deinit();
    try std.testing.expect(only_analysis.analysis != null);
    var only_jobs = try parse(allocator, "version = 1\n\n[[jobs]]\nname = \"j\"\n");
    defer only_jobs.deinit();
    try std.testing.expectEqual(@as(usize, 1), only_jobs.jobs.len);
}

test "parse rejects a [[jobs]] table without keys wherever it is" {
    const allocator = std.testing.allocator;
    const named = "[[jobs]]\nname = \"named\"\nchains = [\"A\"]\n";
    // First, in the middle and last: the TOML parser drops the empty table
    try std.testing.expectError(error.MissingJobName, parse(allocator, "version = 1\n[[jobs]]\n" ++ named));
    try std.testing.expectError(error.MissingJobName, parse(allocator, "version = 1\n" ++ named ++ "[[jobs]]\n" ++ "[[jobs]]\nname = \"other\"\n"));
    try std.testing.expectError(error.MissingJobName, parse(allocator, "version = 1\n" ++ named ++ "[[jobs]]\n"));
    try std.testing.expectError(error.MissingJobName, parse(allocator, "version = 1\n[[ jobs ]] # no name\n"));
    // A table with keys but no name was always rejected
    try std.testing.expectError(error.MissingJobName, parse(allocator, "version = 1\n[[jobs]]\nchains = [\"A\"]\n[[jobs]]\nname = \"x\"\n"));

    // A header in a comment is not a table
    var commented = try parse(allocator, "version = 1\n# [[jobs]]\n" ++ named);
    defer commented.deinit();
    try std.testing.expectEqual(@as(usize, 1), commented.jobs.len);
}

test "parse rejects an empty array table that is not [[jobs]]" {
    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, "version = 1\n[[job]]\n"));
    try std.testing.expectError(error.UnknownField, parse(std.testing.allocator, "version = 1\n[[jobs]]\nname = \"a\"\n[[analysis]]\n"));
}

fn expectFinding(findings: Findings, severity: Severity, needle: []const u8) !void {
    for (findings.slice()) |finding| {
        if (finding.severity == severity and std.mem.find(u8, finding.message, needle) != null) return;
    }
    std.debug.print("no {s} containing '{s}' in:\n", .{ @tagName(severity), needle });
    for (findings.slice()) |finding| std.debug.print("  {s}: {s}\n", .{ @tagName(finding.severity), finding.message });
    return error.TestExpectedFinding;
}

fn expectFindingCount(content: []const u8, mode: Mode, errors: usize, warnings: usize) !void {
    var workflow = try parse(std.testing.allocator, content);
    defer workflow.deinit();
    const findings = checkKeys(workflow, mode);
    try std.testing.expectEqual(errors, findings.errorCount());
    try std.testing.expectEqual(warnings, findings.len - findings.errorCount());
}

const test_legacy_manifest =
    \\version = 1
    \\input_dir = "examples"
    \\output_dir = "out"
    \\format = "jsonl"
    \\n_points = 32
    \\
    \\[[jobs]]
    \\name = "all"
    \\
;

test "checkKeys reports nothing for manifests that fit their command" {
    const calc =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\path = "examples/1crn.pdb"
        \\chain = "A"
        \\model = 1
        \\
        \\[output]
        \\path = "out.json"
        \\format = "json"
        \\
        \\[calculation]
        \\n_points = 32
        \\timing = true
        \\rsa = true
        \\per_residue = true
        \\polar = true
        \\validate_only = false
        \\
        \\[classifier]
        \\type = "ccd"
        \\
    ;
    const batch =
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "examples"
        \\
        \\[output]
        \\dir = "out"
        \\format = "jsonl"
        \\
        \\[output.jsonl]
        \\decimals = 3
        \\
        \\[calculation]
        \\n_points = 32
        \\timing = false
        \\residue_map = true
        \\rsa = false
        \\per_residue = false
        \\polar = false
        \\validate_only = false
        \\
        \\[[jobs]]
        \\name = "all"
        \\
    ;
    const analysis =
        \\version = 1
        \\
        \\[input]
        \\dir = "examples"
        \\
        \\[output]
        \\dir = "out"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\
    ;
    try expectFindingCount(calc, .calc, 0, 0);
    try expectFindingCount(batch, .batch_jobs, 0, 0);
    try expectFindingCount(analysis, .batch_analysis, 0, 0);
    try expectFindingCount(test_legacy_manifest, .batch_jobs, 0, 0);
}

test "checkKeys: calc rejects batch-only structure and warns about batch-only output" {
    const allocator = std.testing.allocator;

    var analysis = try parse(allocator,
        \\version = 1
        \\
        \\[input]
        \\path = "a.pdb"
        \\
        \\[analysis]
        \\type = "bsa"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\
    );
    defer analysis.deinit();
    const analysis_findings = checkKeys(analysis, .calc);
    try std.testing.expectEqual(@as(usize, 1), analysis_findings.len);
    try expectFinding(analysis_findings, .err, "[analysis]");
    try expectFinding(analysis_findings, .err, "zsasa batch --workflow");

    var jobs = try parse(allocator, "version = 1\n[input]\npath = \"a.pdb\"\n[[jobs]]\nname = \"j\"\n");
    defer jobs.deinit();
    const jobs_findings = checkKeys(jobs, .calc);
    try std.testing.expectEqual(@as(usize, 1), jobs_findings.len);
    try expectFinding(jobs_findings, .err, "[[jobs]]");

    var input_dir = try parse(allocator, "version = 1\n[input]\npath = \"a.pdb\"\ndir = \"d\"\n");
    defer input_dir.deinit();
    const dir_findings = checkKeys(input_dir, .calc);
    try std.testing.expectEqual(@as(usize, 1), dir_findings.len);
    try expectFinding(dir_findings, .err, "[input] dir");
    try expectFinding(dir_findings, .err, "zsasa batch --workflow");

    // The flat batch manifest names the key as it is written there
    var legacy = try parse(allocator, test_legacy_manifest);
    defer legacy.deinit();
    const legacy_findings = checkKeys(legacy, .calc);
    try expectFinding(legacy_findings, .err, "input_dir");
    try expectFinding(legacy_findings, .err, "[[jobs]]");
    try expectFinding(legacy_findings, .warning, "output_dir");

    // calc writes no residue map; false asks for what calc does
    try expectFindingCount("version = 1\n[calculation]\nresidue_map = true\n", .calc, 1, 0);
    try expectFindingCount("version = 1\n[calculation]\nresidue_map = false\n", .calc, 0, 0);

    // Output location and JSONL options only warn
    try expectFindingCount("version = 1\n[output]\ndir = \"out\"\n", .calc, 0, 1);
    try expectFindingCount("version = 1\n[output.jsonl]\natom_areas = false\n", .calc, 0, 1);
    try expectFindingCount("version = 1\n[output.jsonl]\nmetadata = \"none\"\ndecimals = 2\n", .calc, 0, 1);
}

test "checkKeys: batch rejects keys it cannot honor and says what to do" {
    const allocator = std.testing.allocator;
    const header = "version = 1\n[input]\ndir = \"d\"\n";
    const job = "\n[[jobs]]\nname = \"j\"\n";

    inline for (.{
        .{ "path = \"a.pdb\"\n", "[input] path", "set [input] dir" },
        .{ "model = 2\n", "[input] model", "zsasa calc --workflow" },
        .{ "mol = \"1\"\n", "[input] mol", "zsasa calc --workflow" },
    }) |case| {
        var workflow = try parse(allocator, header ++ case[0] ++ job);
        defer workflow.deinit();
        const findings = checkKeys(workflow, .batch_jobs);
        try std.testing.expectEqual(@as(usize, 1), findings.len);
        try expectFinding(findings, .err, case[1]);
        try expectFinding(findings, .err, case[2]);
    }

    inline for (.{ "rsa", "per_residue", "polar", "validate_only" }) |key| {
        inline for (.{ Mode.batch_jobs, Mode.batch_analysis }) |mode| {
            const tail = if (mode == .batch_jobs) job else "\n[analysis]\ntype = \"bsa\"\npartner_a = [\"A\"]\npartner_b = [\"B\"]\n";
            var workflow = try parse(allocator, header ++ "[calculation]\n" ++ key ++ " = true\n" ++ tail);
            defer workflow.deinit();
            const findings = checkKeys(workflow, mode);
            try std.testing.expectEqual(@as(usize, 1), findings.len);
            try expectFinding(findings, .err, "[calculation] " ++ key ++ " = true");
            try expectFinding(findings, .err, "zsasa calc --workflow");
        }
    }

    // Where the output goes is a warning
    var output_path = try parse(allocator, header ++ "[output]\npath = \"o.json\"\n" ++ job);
    defer output_path.deinit();
    const findings = checkKeys(output_path, .batch_jobs);
    try std.testing.expectEqual(@as(usize, 1), findings.len);
    try expectFinding(findings, .warning, "[output] path");
    try expectFinding(findings, .warning, "[output] dir");
}

test "checkKeys: an [analysis] workflow rejects chain, residue_map and atom_identity" {
    const base = "version = 1\n[analysis]\ntype = \"bsa\"\npartner_a = [\"A\"]\npartner_b = [\"B\"]\n";
    try expectFindingCount("version = 1\n[input]\nchain = \"A\"\n[analysis]\ntype = \"bsa\"\npartner_a = [\"A\"]\npartner_b = [\"B\"]\n", .batch_analysis, 1, 0);
    try expectFindingCount(base ++ "[calculation]\nresidue_map = true\n", .batch_analysis, 1, 0);
    try expectFindingCount(base ++ "[output.jsonl]\natom_identity = true\n", .batch_analysis, 1, 0);
    try expectFindingCount(base ++ "[output.jsonl]\natom_identity = false\ndecimals = 2\n", .batch_analysis, 0, 0);
}

test "checkKeys: [input] chain is the default selection of the jobs that have none" {
    const allocator = std.testing.allocator;
    const header = "version = 1\n[input]\ndir = \"d\"\nchain = \"A, B\"\n";

    // A job without a selection takes it; a job with its own keeps it
    var workflow = try parse(allocator, header ++ "[[jobs]]\nname = \"default\"\n[[jobs]]\nname = \"own\"\nchains = [\"C\"]\n");
    defer workflow.deinit();
    try std.testing.expectEqual(@as(usize, 0), checkKeys(workflow, .batch_jobs).len);
    try workflow.applyInputChainToJobs();
    try std.testing.expectEqual(@as(usize, 2), workflow.jobs[0].chains.?.len);
    try std.testing.expectEqualStrings("A", workflow.jobs[0].chains.?[0]);
    try std.testing.expectEqualStrings("B", workflow.jobs[0].chains.?[1]);
    try std.testing.expectEqual(@as(usize, 1), workflow.jobs[1].chains.?.len);
    try std.testing.expectEqualStrings("C", workflow.jobs[1].chains.?[0]);

    // Without the key nothing changes
    var plain = try parse(allocator, "version = 1\n[[jobs]]\nname = \"default\"\n");
    defer plain.deinit();
    try plain.applyInputChainToJobs();
    try std.testing.expect(plain.jobs[0].chains == null);

    // A chain_map job selects per file, so the default cannot apply
    var mapped = try parse(allocator, header ++ "[[jobs]]\nname = \"m\"\nchain_map = \"c.csv\"\n[[jobs]]\nname = \"default\"\n");
    defer mapped.deinit();
    try expectFinding(checkKeys(mapped, .batch_jobs), .err, "chain_map");

    // A value no job would use
    var overridden = try parse(allocator, header ++ "[[jobs]]\nname = \"own\"\nchains = [\"C\"]\n");
    defer overridden.deinit();
    try expectFinding(checkKeys(overridden, .batch_jobs), .err, "used by no job");

    // A value without a chain ID
    inline for (.{ "\"\"", "\" , \"" }) |blank| {
        var blank_chain = try parse(allocator, "version = 1\n[input]\nchain = " ++ blank ++ "\n[[jobs]]\nname = \"default\"\n");
        defer blank_chain.deinit();
        try expectFinding(checkKeys(blank_chain, .batch_jobs), .err, "at least one chain ID");
    }
}

test "Findings.report fails on errors and passes warnings" {
    var warnings = Findings{};
    warnings.add(.warning, "only a warning");
    try warnings.report();

    var errors = Findings{};
    errors.add(.warning, "a warning");
    errors.add(.err, "an error");
    try std.testing.expectError(error.InvalidArgument, errors.report());
    try std.testing.expectEqual(@as(usize, 1), errors.errorCount());
}

test "checkKeys never needs more room than Findings has" {
    var workflow = try parse(std.testing.allocator,
        \\version = 1
        \\
        \\[input]
        \\path = "a.pdb"
        \\dir = "d"
        \\chain = "A"
        \\model = 1
        \\mol = "1"
        \\
        \\[output]
        \\path = "o.json"
        \\dir = "o"
        \\
        \\[output.jsonl]
        \\atom_identity = true
        \\
        \\[calculation]
        \\residue_map = true
        \\rsa = true
        \\per_residue = true
        \\polar = true
        \\validate_only = true
        \\
        \\[[jobs]]
        \\name = "j"
        \\
    );
    defer workflow.deinit();
    for ([_]Mode{ .calc, .batch_jobs, .batch_analysis }) |mode| {
        try std.testing.expect(checkKeys(workflow, mode).len < Findings.max_findings);
    }
}

# BSA Analysis Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a BSA/ΔSASA analysis mode to batch workflow files with analysis-specific JSONL output.

**Architecture:** Extend workflow parsing with an optional `[analysis]` table, add a JSONL serializer for BSA analysis rows, and add a file-first batch workflow analysis runner that computes partner A, partner B, and AB complex SASA using existing selection and SASA calculation helpers. Keep normal workflow jobs unchanged when `[analysis]` is absent.

**Tech Stack:** Zig 0.16, existing TOML parser, existing batch workflow file-first parsing/classification path, existing JSON stringify helpers, Docusaurus Markdown docs.

---

### Task 1: Parse `[analysis]` in workflow manifests

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/workflow_manifest.zig`

- [ ] **Step 1: Write failing parser tests**

Add tests in `src/workflow_manifest.zig` near existing workflow parser tests:

```zig
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
```

- [ ] **Step 2: Run parser tests to verify RED**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bsa-parser-red.log
```

Expected: FAIL because `Workflow` has no `analysis` field and `InvalidAnalysisConfig` is undefined.

- [ ] **Step 3: Implement manifest parsing**

Add `InvalidAnalysisConfig` to `WorkflowError`; add:

```zig
pub const Analysis = struct {
    type: ?[]const u8 = null,
    name: ?[]const u8 = null,
    partner_a: ?[]const []const u8 = null,
    partner_b: ?[]const []const u8 = null,
    level: ?[]const u8 = null,
};
```

Add `analysis: ?Analysis = null` to `Workflow`, free `partner_a` and `partner_b` in `Workflow.deinit`, allow table name `analysis` in `validateKnownDocumentShape`, parse it in `parseOwned`, and implement:

```zig
fn parseAnalysis(allocator: Allocator, table: toml_parser.Table) Error!Analysis {
    try rejectUnknownFields(table.entries, &.{ "type", "name", "partner_a", "partner_b", "level" });
    const analysis = Analysis{
        .type = try optionalString(table.entries, "type"),
        .name = try optionalString(table.entries, "name"),
        .partner_a = try optionalStringArray(allocator, table.entries, "partner_a"),
        .partner_b = try optionalStringArray(allocator, table.entries, "partner_b"),
        .level = try optionalString(table.entries, "level"),
    };
    try validateAnalysis(analysis);
    return analysis;
}
```

`validateAnalysis` accepts only `type = "bsa"`, requires non-empty partner arrays, checks `name` with `isSafeJobName` when present, and accepts only absent/`"total"`/`"residue"` level.

- [ ] **Step 4: Run parser tests to verify GREEN**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bsa-parser-green.log
```

Expected: PASS.

### Task 2: Add BSA analysis JSONL serialization

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/json_writer.zig`

- [ ] **Step 1: Write failing serializer test**

Add a test near existing JSONL tests:

```zig
test "BSA analysis JSONL includes total and residue delta fields" {
    const allocator = std.testing.allocator;
    const partner_a = [_][]const u8{"A"};
    const partner_b = [_][]const u8{"B"};
    const residue_chain = [_][]const u8{ "A", "B" };
    const residue_name = [_][]const u8{ "GLY", "ALA" };
    const residue_number = [_]i32{ 1, 2 };
    const residue_insertion_code = [_][]const u8{ "", "" };
    const residue_delta_sasa = [_]f64{ 3.0, 5.0 };

    const line = try bsaAnalysisToJsonlLine(allocator, .{
        .filename = "tiny.pdb",
        .name = "interface_ab",
        .partner_a = partner_a[0..],
        .partner_b = partner_b[0..],
        .sasa_partner_a = 10.0,
        .sasa_partner_b = 20.0,
        .sasa_complex = 14.0,
        .delta_sasa_total = 16.0,
        .bsa = 8.0,
        .delta_sasa_level = "residue",
        .residue_chain = residue_chain[0..],
        .residue_name = residue_name[0..],
        .residue_number = residue_number[0..],
        .residue_insertion_code = residue_insertion_code[0..],
        .residue_delta_sasa = residue_delta_sasa[0..],
    });
    defer allocator.free(line);

    try std.testing.expect(std.mem.indexOf(u8, line, "\"analysis\":\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"delta_sasa_total\":16") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"bsa\":8") != null);
    try std.testing.expect(std.mem.indexOf(u8, line, "\"residue_delta_sasa\":[3,5]") != null);
}
```

- [ ] **Step 2: Run serializer test to verify RED**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bsa-json-red.log
```

Expected: FAIL because `bsaAnalysisToJsonlLine` is undefined.

- [ ] **Step 3: Implement serializer**

Add public input struct and function in `src/json_writer.zig`:

```zig
pub const BsaAnalysisJsonl = struct {
    filename: []const u8,
    name: []const u8,
    partner_a: []const []const u8,
    partner_b: []const []const u8,
    sasa_partner_a: f64,
    sasa_partner_b: f64,
    sasa_complex: f64,
    delta_sasa_total: f64,
    bsa: f64,
    delta_sasa_level: []const u8,
    residue_chain: []const []const u8 = &.{},
    residue_name: []const []const u8 = &.{},
    residue_number: []const i32 = &.{},
    residue_insertion_code: []const []const u8 = &.{},
    residue_delta_sasa: []const f64 = &.{},
};
```

Implement `bsaAnalysisToJsonlLine` with separate total-only and residue entry structs so residue arrays are omitted when `delta_sasa_level` is `"total"`.

- [ ] **Step 4: Run serializer test to verify GREEN**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bsa-json-green.log
```

Expected: PASS.

### Task 3: Implement BSA analysis workflow runner

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/src/batch.zig`

- [ ] **Step 1: Write failing workflow integration test**

Add a test near existing workflow integration tests in `src/batch.zig` that creates a temporary input directory with a two-chain PDB, writes a workflow with `[analysis]`, runs `runWorkflow`, and checks `<output>/interface_ab.jsonl` contains BSA fields:

```zig
test "workflow BSA analysis writes analysis JSONL" {
    const allocator = std.testing.allocator;
    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();

    var root_buf: [std.fs.max_path_bytes]u8 = undefined;
    const root_len = try tmp_dir.dir.realPath(std.testing.io, &root_buf);
    const root = root_buf[0..root_len];

    const input_dir = try std.fs.path.join(allocator, &.{ root, "input" });
    defer allocator.free(input_dir);
    const output_dir = try std.fs.path.join(allocator, &.{ root, "output" });
    defer allocator.free(output_dir);
    try std.Io.Dir.cwd().createDirPath(std.testing.io, input_dir);

    const pdb_path = try std.fs.path.join(allocator, &.{ input_dir, "tiny_ab.pdb" });
    defer allocator.free(pdb_path);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = pdb_path, .data =
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      2  CA  GLY A   1       1.500   0.000   0.000  1.00 20.00           C  \n" ++
        "ATOM      3  N   ALA B   2       3.000   0.000   0.000  1.00 20.00           N  \n" ++
        "ATOM      4  CA  ALA B   2       4.500   0.000   0.000  1.00 20.00           C  \n" ++
        "END\n",
    });

    const workflow_path = try std.fs.path.join(allocator, &.{ root, "bsa.toml" });
    defer allocator.free(workflow_path);
    const workflow = try std.fmt.allocPrint(allocator,
        \\version = 1
        \\kind = "workflow"
        \\
        \\[input]
        \\dir = "{s}"
        \\
        \\[output]
        \\dir = "{s}"
        \\format = "jsonl"
        \\
        \\[calculation]
        \\n_points = 16
        \\quiet = true
        \\
        \\[analysis]
        \\type = "bsa"
        \\name = "interface_ab"
        \\partner_a = ["A"]
        \\partner_b = ["B"]
        \\level = "residue"
    , .{ input_dir, output_dir });
    defer allocator.free(workflow);
    try std.Io.Dir.cwd().writeFile(std.testing.io, .{ .sub_path = workflow_path, .data = workflow });

    try runWorkflow(allocator, std.testing.io, .{ .workflow_path = workflow_path });

    const jsonl_path = try std.fs.path.join(allocator, &.{ output_dir, "interface_ab.jsonl" });
    defer allocator.free(jsonl_path);
    const content = try std.Io.Dir.cwd().readFileAlloc(std.testing.io, jsonl_path, allocator, .limited(8192));
    defer allocator.free(content);

    try std.testing.expect(std.mem.indexOf(u8, content, "\"analysis\":\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"delta_sasa_total\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"bsa\"") != null);
    try std.testing.expect(std.mem.indexOf(u8, content, "\"residue_delta_sasa\"") != null);
}
```

- [ ] **Step 2: Run integration test to verify RED**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bsa-workflow-red.log
```

Expected: FAIL because workflow analysis runner is not implemented.

- [ ] **Step 3: Implement analysis runner**

In `src/batch.zig`, add helper structs/functions for BSA analysis:

- `analysisOutputName(analysis)` returns `analysis.name orelse "bsa"`.
- `analysisOutputPath(allocator, output_dir, name)` returns `<output_dir>/<name>.jsonl`.
- `appendStringArrays(allocator, a, b)` creates complex chain selection.
- `calculateSelectedSasaForAnalysis` wraps `copySelectedAtomInput` plus `calculatePreparedInputResult` and returns selected input and result while preserving atom areas for ΔSASA.
- `buildResidueDeltaArrays` maps selected partner atom deltas to residue arrays using the selected inputs' residue metadata.
- `runWorkflowBsaAnalysis` mirrors the sequential file-first path: load resources, scan files, parse/classify once, calculate A/B/AB SASA, write one BSA JSONL row per successful file, and print `Workflow complete` counts.

At the start of `runWorkflow`, after parsing the workflow, dispatch to `runWorkflowBsaAnalysis` when `workflow.analysis != null`.

- [ ] **Step 4: Run integration test to verify GREEN**

Run:

```bash
zig build test --summary all 2>&1 | tee /tmp/zsasa-bsa-workflow-green.log
```

Expected: PASS.

### Task 4: Document BSA analysis workflows

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/guide/workflows.md`

- [ ] **Step 1: Add workflow docs**

Add a "BSA / ΔSASA Analysis" section after the batch workflow example. Include the TOML example from the design, the formulas `delta_sasa_total = sasa_partner_a + sasa_partner_b - sasa_complex` and `bsa = delta_sasa_total / 2`, and mention that analysis JSONL does not use normal SASA `total_area`/`atom_areas` fields.

- [ ] **Step 2: Run docs-relevant grep check**

Run:

```bash
rg -n "BSA|ΔSASA|delta_sasa_total|partner_a" website/docs/guide/workflows.md
```

Expected: Shows the new section and fields.

### Task 5: Format and verify

**Files:**
- Verify all modified files.

- [ ] **Step 1: Format Zig files**

Run:

```bash
zig fmt src/workflow_manifest.zig src/json_writer.zig src/batch.zig
```

- [ ] **Step 2: Run focused checks**

Run:

```bash
zig fmt --check src/
zig build test
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa --version
mkdir -p /tmp/zsasa-check
./zig-out/bin/zsasa calc examples/1ubq.pdb /tmp/zsasa-check/output.json
```

Expected: all commands exit 0.

- [ ] **Step 3: Review diff**

Run:

```bash
git diff -- src/workflow_manifest.zig src/json_writer.zig src/batch.zig website/docs/guide/workflows.md docs/superpowers/specs/2026-06-28-bsa-analysis-workflow-design.md docs/superpowers/plans/2026-06-28-bsa-analysis-workflow.md
```

Expected: Changes match the design and contain no debug prints, placeholders, or unrelated edits.

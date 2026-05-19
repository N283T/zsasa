# Workflow Manifest Design

## Summary

Promote the current batch-only TOML manifest into a general zsasa workflow file that can configure `calc`, `batch`, and future trajectory workflows through one consistent schema. The workflow file should own reproducible calculation settings, input/output routing, classifier selection, and batch job definitions.

Custom classifiers will be referenced from the workflow file instead of embedded inline. Custom classifier configuration will be TOML-only; the legacy FreeSASA-style custom classifier format will be removed in this change and documented with a migration example.

## Goals

- Support a shared workflow schema for `calc` and `batch`, with room for future trajectory workflows.
- Keep `zsasa batch --manifest existing.toml` working for existing batch manifest files where practical.
- Add `--workflow=PATH` as the preferred CLI spelling for workflow files.
- Keep `--manifest=PATH` as a batch compatibility alias.
- Configure custom classifiers from workflow files via a referenced TOML classifier config.
- Remove FreeSASA-format custom classifier parsing and auto-detection.
- Update docs and help text with TOML-only custom classifier migration guidance.

## Non-goals

- Do not embed custom classifier `[types]` and `[[atoms]]` definitions directly inside workflow files in the first implementation.
- Do not add YAML or JSON workflow support.
- Do not redesign output result schemas beyond the existing manifest-related routing.
- Do not add a full workflow engine or dependency graph. A workflow file describes one calculation mode at a time.

## User-facing CLI

Preferred new usage:

```bash
zsasa calc --workflow sasa.toml
zsasa batch --workflow bsa.toml
```

Compatibility usage:

```bash
zsasa batch --manifest bsa.toml
```

`--manifest` should remain accepted by `batch` as an alias for `--workflow`. New docs should prefer `--workflow` except where explaining backwards compatibility.

For `calc`, add `--workflow=PATH` / `--workflow PATH`. Do not add `--manifest` to `calc`; the old name was batch-specific and should not spread to new surfaces.

CLI precedence remains:

```text
built-in defaults
  < workflow globals
  < workflow job-specific settings, for batch jobs
  < explicit CLI options
```

Explicit CLI options should override workflow values. Existing direct CLI usage without a workflow file must keep current behavior.

## Workflow schema

Use TOML with `version = 1`. `kind = "workflow"` is optional in v1 so legacy batch manifests with only `version = 1` can continue to parse.

### Single calculation example

```toml
version = 1
kind = "workflow"

[input]
path = "structure.cif"

[output]
path = "result.json"
format = "json"

[calculation]
algorithm = "sr"
probe_radius = 1.4
n_points = 128
n_slices = 20
precision = "f64"
use_bitmask = true
include_hydrogens = false
include_hetatm = false
auth_chain = false
per_residue = true
rsa = true
polar = true

[classifier]
type = "ccd"
ccd = "components.zsdc"
sdf = ["ligand.sdf"]
```

### Batch example

```toml
version = 1
kind = "workflow"

[input]
dir = "structures"

[output]
dir = "results"
format = "jsonl"

[calculation]
algorithm = "sr"
n_points = 128
use_bitmask = true
residue_map = true
threads = 8

[classifier]
type = "custom"
config = "my_radii.toml"

[[jobs]]
name = "chain_A"
chains = ["A"]

[[jobs]]
name = "chain_B"
chains = ["B"]

[[jobs]]
name = "complex_AB"
chains = ["A", "B"]
```

### Sections

#### `[input]`

- `path`: single input file for `calc`.
- `dir`: input directory for `batch`.
- `chain`: optional chain filter for `calc`, as a string such as `"A"` or `"A,B"`.
- `model`: optional model number for NMR structures in `calc`.
- `mol`: optional SDF molecule selector for `calc`.

A workflow used with `calc` must resolve to `input.path`. A workflow used with `batch` must resolve to `input.dir`.

#### `[output]`

- `path`: output file for `calc`.
- `dir`: output directory for `batch`.
- `format`: `json`, `compact`, `csv`, or `jsonl`.

For batch JSONL output, each job continues to write `<output.dir>/<job.name>.jsonl`. For per-file batch formats, each job continues to write under `<output.dir>/<job.name>/`.

#### `[calculation]`

Shared options:

- `algorithm`: `sr` or `lr`.
- `probe_radius`: finite positive radius in Å.
- `n_points`: SR point count.
- `n_slices`: LR slice count.
- `precision`: `f32` or `f64`.
- `use_bitmask`: boolean.
- `include_hydrogens`: boolean.
- `include_hetatm`: boolean.
- `auth_chain`: boolean.
- `timing`: boolean.
- `quiet`: boolean.

Calc-specific options:

- `per_residue`: boolean.
- `rsa`: boolean.
- `polar`: boolean.
- `validate_only`: boolean.

Batch-specific options:

- `threads`: file-level worker count.
- `residue_map`: boolean, valid only with `format = "jsonl"`.

Unknown fields should be rejected with a clear error so misspellings do not silently change calculations.

#### `[classifier]`

Built-in classifiers:

```toml
[classifier]
type = "ccd"
ccd = "components.zsdc"
sdf = ["ligand.sdf"]
```

- `type`: `ccd`, `protor`, `naccess`, `oons`, or `custom`.
- `ccd`: optional external CCD dictionary path, valid with `ccd` / `protor`.
- `sdf`: optional string or string array of SDF topology files, valid with `ccd` / `protor`.

Custom classifier:

```toml
[classifier]
type = "custom"
config = "my_radii.toml"
```

`type = "custom"` requires `config`. The referenced config must be TOML and parsed by the existing TOML classifier parser.

If both `type = "custom"` and built-in-only fields such as `ccd` are present, reject the workflow. If `config` is present for a built-in classifier, reject the workflow.

#### `[[jobs]]`

Batch workflows can define one or more jobs:

- `name`: required, non-empty, path-safe.
- `chains`: optional string array. Omitted means all chains.
- `auth_chain`: optional boolean override for that job.

The first implementation should keep job-specific overrides limited to chain selection and `auth_chain`, matching current batch manifest behavior.

## Legacy batch manifest compatibility

Existing flat batch manifests should continue to parse as legacy v1 batch workflows:

```toml
version = 1
input_dir = "structures"
output_dir = "results"
format = "jsonl"
classifier = "ccd"
use_bitmask = true
n_points = 128

[[jobs]]
name = "chain_A"
chains = ["A"]
```

Map legacy fields to the new internal representation:

- `input_dir` -> `[input].dir`
- `output_dir` -> `[output].dir`
- `format` -> `[output].format`
- calculation fields -> `[calculation]`
- `classifier`, `ccd`, `sdf` -> `[classifier]`

Docs should show the new sectioned form. Legacy flat form should be documented only as accepted compatibility syntax, not as the recommended style.

## FreeSASA custom classifier format removal

Remove support for legacy FreeSASA-style custom classifier config files such as:

```text
name: my-classifier

types:
C_ALI 1.87 apolar

atoms:
ANY CA C_ALI
```

TOML becomes the only supported custom classifier format:

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"
```

Behavior changes:

- `--config=FILE` only accepts TOML classifier files.
- Workflow `[classifier].config` only accepts TOML classifier files.
- Extension-based auto-detection is removed.
- If a non-TOML-looking config is supplied, return an error that says custom classifier configs are TOML-only and points to the migration docs.

Because CCD is the recommended classifier for most use cases, docs should position custom classifiers as an advanced feature.

## Internal architecture

### Module structure

Introduce a shared workflow manifest module, for example `src/workflow_manifest.zig`, and keep batch-specific compatibility contained there.

Suggested public data model:

```zig
pub const Workflow = struct {
    allocator: Allocator,
    content: []const u8,
    input: Input,
    output: Output,
    calculation: Calculation,
    classifier: ClassifierConfig,
    jobs: []Job,
    is_legacy_batch_manifest: bool,
};

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
```

`src/batch_manifest.zig` can either be removed or temporarily become a small compatibility wrapper around `workflow_manifest.zig`. Prefer updating imports to the new module and deleting the old module if the patch remains manageable.

### Resolution layer

Keep parsing separate from command-specific resolution:

- `calc` resolves `Workflow` + `CalcArgs` into `CalcArgs` or a new resolved calc config.
- `batch` resolves `Workflow` + `BatchArgs` into per-job `BatchConfig` values.

This prevents workflow parsing from depending on execution modules and keeps parser tests cheap.

### Classifier loading

Create or reuse a helper that loads custom classifier TOML by path. Both `calc` direct `--config` and workflow `[classifier].config` should use the same TOML-only path.

The existing `toml_classifier_parser.zig` should remain the parser for TOML classifier configs. `classifier_parser.zig` should be removed or reduced to a TOML-only compatibility facade if that minimizes churn.

## Error handling

Return concise errors for:

- Unsupported workflow `version`.
- Invalid or unsupported `kind`.
- Missing `input.path` for `calc` after CLI/workflow merging.
- Missing `input.dir` for `batch` after CLI/workflow merging.
- Both `[input].path` and `[input].dir` used with the wrong command.
- Multiple batch jobs without `output.dir` when output would be ambiguous.
- Duplicate or unsafe batch job names.
- `--chain` combined with batch jobs, unless the CLI chain is explicitly treated as a global override. Preserve current safer behavior: reject it and direct users to `[[jobs]].chains`.
- `residue_map = true` with a non-JSONL output format.
- Invalid classifier combinations such as `type = "custom"` without `config`.
- Legacy FreeSASA-format config paths, with a TOML migration message.

## Documentation updates

Update at least:

- CLI command docs for `calc --workflow` and `batch --workflow`.
- Batch manifest docs, renamed or reframed as workflow docs.
- Classifier guide to remove FreeSASA format as a supported option.
- CLI input docs for TOML-only custom classifier configs.
- Changelog with a breaking-change entry.

Include a migration example:

```text
# Old FreeSASA format
C_ALI 1.87 apolar
ANY CA C_ALI
```

becomes:

```toml
[types]
C_ALI = { radius = 1.87, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"
```

## Testing

Add focused tests for:

- Parsing new sectioned calc workflow.
- Parsing new sectioned batch workflow with jobs.
- Parsing existing flat batch manifest as legacy compatibility syntax.
- Rejecting invalid classifier combinations.
- Resolving custom classifier config paths for `calc`.
- Resolving custom classifier config paths for `batch`.
- Rejecting FreeSASA-format custom classifier files with a migration-oriented error.
- CLI parsing for `calc --workflow`.
- CLI parsing for `batch --workflow` and `batch --manifest` alias.
- Existing batch `--manifest` behavior with a legacy flat fixture.

Run at least:

```bash
zig build test
zig build run -- calc --workflow <fixture.toml>
zig build run -- batch --workflow <fixture.toml>
zig build run -- batch --manifest <legacy-fixture.toml>
```

Use the smallest fixtures that verify config merge and classifier loading without making tests slow.

## Implementation notes

Keep the first implementation small by:

- Supporting referenced custom classifier files only.
- Keeping batch job overrides limited to `chains` and `auth_chain`.
- Preserving current batch execution flow.
- Avoiding trajectory support until a future workflow consumer is designed.
- Updating docs in the same change so the breaking classifier format change is clear.

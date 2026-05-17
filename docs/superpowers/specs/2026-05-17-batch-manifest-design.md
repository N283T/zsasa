# Batch TOML Manifest Design

## Summary

Add a TOML manifest mode for `zsasa batch` so one command can run multiple named batch jobs over the same structure directory. The motivating workflow is homo-dimer/interface analysis: calculate atom-level SASA for chain A, chain B, and the AB complex without writing an external loop.

The manifest is self-contained when it defines input and output paths, but CLI positional arguments and explicit options can override manifest values. This keeps simple scripted use possible while supporting reusable workflow files.

## Goals

- Support `zsasa batch --manifest bsa.toml` as a self-contained batch workflow.
- Support `zsasa batch structures/ results/ --manifest bsa.toml` as CLI path override/fallback.
- Support multiple named jobs, each with its own chain selection.
- Keep output separated by job name.
- Preserve existing `zsasa batch` behavior when no manifest is provided.
- Reuse existing project patterns and avoid adding a YAML parser dependency.

## Non-goals

- Do not add YAML support in the first implementation.
- Do not compute buried SASA/BSA directly in this feature; the manifest produces the required A, B, and AB atom-level SASA outputs for downstream analysis.
- Do not redesign the batch output schema beyond job-level output placement.
- Do not add Python API manifest support in the first implementation.

## Manifest format

The first implementation uses TOML and a required `version = 1` field.

Example:

```toml
version = 1

input_dir = "structures"
output_dir = "results"

format = "jsonl"
use_bitmask = true
n_points = 128
classifier = "ccd"

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

### Global fields

Manifest globals map to existing batch settings:

- `input_dir`: input directory, optional if CLI positional input is supplied.
- `output_dir`: output directory, optional; current no-file-output behavior remains available when neither manifest nor CLI supplies output.
- `algorithm`: `sr` or `lr`.
- `classifier`: `ccd`, `protor`, `naccess`, or `oons`.
- `ccd`: external CCD dictionary path.
- `sdf`: string or array of strings for SDF topology files.
- `threads`: number of file-level worker threads.
- `probe_radius`: probe radius in Å.
- `n_points`: SR point count.
- `n_slices`: LR slice count.
- `precision`: `f32` or `f64`.
- `format`: `json`, `compact`, `csv`, or `jsonl`.
- `include_hydrogens`: boolean.
- `include_hetatm`: boolean.
- `use_bitmask`: boolean.
- `timing`: boolean.
- `quiet`: boolean.
- `auth_chain`: boolean, using `auth_asym_id` instead of `label_asym_id` for mmCIF chain matching.

### Job fields

Each `[[jobs]]` table defines one run over the input directory:

- `name`: required, used for output placement and user-facing status.
- `chains`: optional array of chain IDs. If omitted, the job uses all chains.
- Any safe global calculation option may be overridden per job, but the first implementation can limit job-specific overrides to `chains` and `auth_chain` if that keeps the patch small.

Job names must be non-empty and path-safe. The parser should reject names containing `/`, `\\`, or `..` path traversal components.

## Precedence rules

Configuration is resolved in this order:

```text
built-in defaults
  < manifest global defaults
  < manifest job-specific settings
  < explicit CLI options
```

Path handling:

- `zsasa batch --manifest bsa.toml` uses `input_dir` and `output_dir` from the manifest.
- `zsasa batch structures/ results/ --manifest bsa.toml` overrides manifest `input_dir` and `output_dir`.
- If no input directory is available after merging, the command errors with a clear message.

Option handling:

- Existing CLI options continue to work and override manifest globals when explicitly supplied.
- To distinguish explicit CLI options from defaults, `BatchArgs` needs presence flags for options that can come from the manifest.
- `--chain` is useful for non-manifest single-job batch mode, but in manifest mode it should not silently override all job chain selections. The initial behavior should reject `--manifest` combined with `--chain` and explain that chain selections belong in `[[jobs]].chains`.

## Output layout

Each job writes to a separate output target under the resolved output directory.

For JSONL:

```text
results/
  chain_A.jsonl
  chain_B.jsonl
  complex_AB.jsonl
```

For per-file formats (`json`, `compact`, `csv`):

```text
results/
  chain_A/
    input1.json
    input2.json
  chain_B/
    input1.json
    input2.json
  complex_AB/
    input1.json
    input2.json
```

If no output directory is configured and `format = "jsonl"`, existing stdout behavior can remain for a single job. Multiple manifest jobs should require an output directory to avoid interleaving unrelated JSONL streams on stdout.

## Architecture

### Parser module

Add a focused manifest parser module, for example `src/batch_manifest.zig`, rather than expanding `src/batch.zig` further. It should:

1. Read a TOML file.
2. Validate `version = 1`.
3. Parse globals into a manifest config structure.
4. Parse `[[jobs]]` into named job definitions.
5. Return structured validation errors that can be rendered as concise CLI messages.

The project already has a minimal TOML parser for classifier configuration. Reuse it if it supports arrays of tables and the needed value types; otherwise extend it narrowly for this manifest shape or use `std.json` only for tests/fixtures if TOML support is insufficient. Do not introduce YAML support.

### Batch configuration merge

Add an intermediate resolved-job structure:

```zig
const ResolvedBatchJob = struct {
    name: []const u8,
    input_dir: []const u8,
    output_dir: ?[]const u8,
    config: BatchConfig,
};
```

Resolution should be separate from execution so it can be unit-tested without running SASA.

### Chain filtering

Extend `BatchConfig` with:

- `chain_filter: ?[]const []const u8`
- `use_auth_chain: bool`

When reading PDB/mmCIF files in batch mode, pass these settings to `PdbParser` and `MmcifParser`. JSON input is not chain-filtered in the first implementation unless a safe existing helper is available.

### Execution flow

Without a manifest, `zsasa batch` follows the current path.

With a manifest:

1. Parse CLI args including `--manifest`.
2. Parse the TOML manifest.
3. Resolve one `ResolvedBatchJob` per manifest job.
4. For each resolved job, call the existing `runBatch` path with that job's `BatchConfig` and output target.
5. Aggregate and print a concise manifest summary after all jobs complete.

The first implementation can run jobs sequentially at the manifest-job level. File-level parallelism within each job remains controlled by `threads`.

## Error handling

- Missing manifest file: clear file-open error including the path.
- Unsupported manifest version: explain supported version is `1`.
- Missing input directory after merge: explain how to fix with manifest `input_dir` or CLI positional input.
- Duplicate job names: reject before running any job.
- Unsafe job names: reject before running any job.
- Multiple jobs without output directory: reject if output would otherwise go to stdout/no files in a confusing way.
- `--manifest` with `--chain`: reject and direct users to `[[jobs]].chains`.

## Documentation

Update CLI docs and help text with:

- `--manifest=PATH` / `--manifest PATH`.
- A homo-dimer A/B/AB example.
- Precedence rules.
- Output layout for JSONL and per-file formats.

## Testing

Add focused tests for:

- Manifest parsing of the A/B/AB example.
- CLI positional input/output overriding manifest paths.
- CLI options overriding manifest globals.
- Duplicate and unsafe job names.
- `--manifest` plus `--chain` rejection.
- Batch parser accepting `--chain` in non-manifest mode.
- PDB/mmCIF batch chain filtering smoke tests using a small multi-chain fixture.
- JSONL output path generation per job.

Run at least:

```bash
zig build test
zig build run -- batch --manifest <fixture.toml>
```

Use the narrowest smoke fixture that verifies chain A, chain B, and AB produce distinct atom counts or chain IDs in output.

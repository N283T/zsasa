---
sidebar_position: 2
---

# Batch Processing

Use batch mode when you need to calculate SASA for many structure files with the same settings.

## Basic Directory Batch

```bash
zsasa batch structures/ results/
```

This scans `structures/` for supported input files and writes per-file outputs under `results/`.

Common options:

```bash
zsasa batch structures/ results/ --threads=8 --format=json
zsasa batch structures/ results/ --format=jsonl --output=results.jsonl
zsasa batch structures/ results/ --classifier=ccd --ccd=components.zsdc
```

## JSONL for Large Runs

For large datasets, prefer JSONL because each structure result is written as one line and can be streamed by downstream tools:

```bash
zsasa batch structures/ results/ --format=jsonl --output=results.jsonl
```

JSONL is especially useful when you want to concatenate, filter, or process results incrementally.

Successful JSONL rows include `status: "ok"` plus the result fields:

```json
{"status":"ok","filename":"1ubq.pdb","total_area":5656.6511885225,"atom_areas":[0,6.759702080876085]}
```

Failed structures are emitted as `status: "err"` rows instead of being available only in the batch summary:

```json
{"status":"err","filename":"bad.pdb","error":"read/parse failed: InvalidFormat"}
```

Parallel JSONL is streamed in completion order for throughput; input-order output is not currently guaranteed.

Use `--jsonl-decimals=N` to round JSONL floating-point values and reduce output size:

```bash
zsasa batch structures/ --format=jsonl --output=results.jsonl --jsonl-decimals=3
```

Rounding applies to JSONL floating-point fields such as `total_area`, `atom_areas`, residue SASA arrays, and BSA analysis values. Omitting the option preserves full precision.

## Thread Count for Large File Sets

By default, `zsasa batch` uses the detected CPU count. For large directories with many small structure files, the workload can spend significant time waiting on file open/read/parse operations. In those I/O-bound cases, you can explicitly set `--threads=N` above the CPU count to keep more files in flight:

```bash
zsasa batch structures/ results.jsonl \
  --format=jsonl \
  --threads=40 \
  --precision=f32 \
  --use-bitmask
```

This can improve throughput on fast local SSDs, but the best value is machine- and dataset-dependent. Higher thread counts increase memory use and may stop helping once storage or CPU scheduling becomes saturated.

## Experimental Adaptive Bitmask SR

For large SR batch jobs that already use bitmask mode, `--adaptive-sr` runs a coarse bitmask pass for every atom and recomputes only intermediate-exposure atoms with a fine point count:

```bash
zsasa batch structures/ results/ \
  --use-bitmask \
  --adaptive-sr \
  --coarse-points=64 \
  --fine-points=256 \
  --adaptive-low=0.10 \
  --adaptive-high=0.90
```

Adaptive mode is currently available for `zsasa batch` only. It requires `--use-bitmask` and `--algorithm=sr`. The output schema is unchanged; compare against fixed fine-point bitmask runs when validating a new dataset.

## Residue Maps in JSONL

Add `--residue-map` to include compact residue-level arrays in each JSONL row:

```bash
zsasa batch structures/ results/ \
  --format=jsonl \
  --output=results.jsonl \
  --residue-map
```

This adds these arrays to each JSONL result:

- `residue_chain`
- `residue_name`
- `residue_number`
- `residue_insertion_code`
- `residue_atom_start`
- `residue_atom_count`
- `residue_sasa`

Without `--residue-map`, result rows include only `status`, `filename`, `total_area`, and `atom_areas`.

## Chain Filters

For a single non-workflow batch job, use `--chain` or `--auth-chain`:

```bash
zsasa batch structures/ results/ --chain=A
zsasa batch structures/ results/ --auth-chain --chain=A
```

Use `--chain` for label/asym chain IDs. Add `--auth-chain` as a boolean modifier when the `--chain` value should match author-provided chain IDs in mmCIF inputs.

For named multi-job runs such as chain A, chain B, and AB complex calculations, use [Workflow Files](workflows.md) instead of repeating shell commands.

## When to Use Workflow Files

Use workflow files when you need:

- Reproducible settings checked into a project.
- Multiple named jobs in one batch run.
- Per-job chain filters.
- Shared classifier, output, and calculation settings.
- Custom classifier TOML configs for batch jobs.

See [Workflow Files](workflows.md) for TOML examples.

## Reference

- [CLI Commands](../cli/commands.md)
- [Input Formats](../cli/input.md)
- [Output Formats](../cli/output.md)

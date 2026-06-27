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

## Experimental Large-Directory Chunking

For very large directories where per-file scheduling or JSONL write overhead may dominate, `zsasa batch` exposes experimental chunk options for benchmarking:

```bash
# Simple chunk scheduling baseline: output path and schema stay the same
zsasa batch structures/ results.jsonl \
  --format=jsonl \
  --threads=10 \
  --chunk-size=256

# Chunked JSONL writer: serialize rows normally, but write/flush once per chunk
zsasa batch structures/ results.jsonl \
  --format=jsonl \
  --threads=10 \
  --chunk-size=256 \
  --chunked-jsonl
```

`--chunk-size=N` changes only how parallel workers claim work ranges. `--chunked-jsonl` additionally buffers JSONL output per chunk and requires `--format=jsonl`. These options are experimental and intended for comparing SwissProt-scale runs before choosing a stable large-batch interface.

To test macro-sharding separately, `--shard-size=N` splits a JSONL batch into large output shards:

```bash
zsasa batch structures/ results.jsonl \
  --format=jsonl \
  --threads=10 \
  --shard-size=50000
```

This writes `results.part-0.jsonl`, `results.part-1.jsonl`, and so on. The JSONL row schema is unchanged; concatenate the shard files if a single stream is needed.

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

The default JSONL schema is unchanged unless `--residue-map` is enabled.

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

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

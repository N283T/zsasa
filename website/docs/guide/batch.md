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

Each output file is named after the input stem: `1ubq.pdb` and `1ubq.cif.gz` would both be written to `results/1ubq.json`. If two inputs would write the same output file, batch mode lists them and exits with an error before processing anything, instead of letting one result overwrite another:

```text
Error: 1 output name is shared by more than one input:
  1ubq.json <- 1ubq.cif.gz, 1ubq.pdb
```

Output names are compared without regard to ASCII case, because names that differ only in case are one file on the default macOS and Windows filesystems. `PROT.pdb` and `prot.cif.gz` are rejected on every platform, also where both files could be written, and the message says that case is the only difference:

```text
Error: 1 output name is shared by more than one input:
  PROT.json <- PROT.pdb, prot.cif.gz (the output names differ only in case)
```

Split those inputs into separate directories, or use JSONL output, which keeps one record per input.

### SDF and MOL Output Names

An SDF or MOL file is expanded into one result per molecule. Each molecule gets a name made of the file stem and the molecule title. This name is the `filename` of the molecule in JSONL rows and in `process_directory()` results, and the per-file output is named after it.

| Molecule | Name |
|----------|------|
| Has a title | `stem_title` |
| Title is blank | `stem_N`, where `N` is the position of the molecule in the file, starting at 1 |
| Shares that name with other molecules of the file | `stem_title_N` for each of them |

`dup.sdf` with two molecules titled `ethanol` gives `dup_ethanol_1` and `dup_ethanol_2`. If a name with the position appended is still the name of another molecule of the file, `_N` is appended again: the titles `x`, `x` and `x_2` in `c.sdf` give `c_x_1`, `c_x_2_2` and `c_x_2`. A molecule whose name is not shared always keeps the plain `stem_title` or `stem_N`.

The stem is the file name without `.sdf` or `.mol` (and `.gz` or `.zst`). The output file is the molecule name with the output extension appended, so dots in the stem or the title are kept: the molecule `v1.5` in `lig.v2.sdf` is written to `lig.v2_v1.5.json`.

A title is used in a file name, so characters that cannot be part of one are replaced by `_` in the output file name: `/`, `\` and control characters on every platform, and `< > : " | ? *` on Windows. The molecule `a/b` in `lig.sdf` is written to `lig_a_b.json`, and no title can write outside the output directory. The `filename` in JSONL rows and API results keeps the title as written (`lig_a/b`). Two names count as shared when they would be written to the same file, so titles that differ only in ASCII case or only in a replaced character also get their position appended (`Ethanol` and `ethanol` in `c.sdf` give `c_Ethanol_1` and `c_ethanol_2`).

Molecule outputs take part in the collision check with their real names. `lig.sdf` and `lig.mol` can be processed together as long as their molecules have different names, while an unnamed first molecule in `lig.sdf` collides with `lig_1.pdb`:

```text
Error: 1 output name is shared by more than one input:
  lig_1.json <- lig.sdf (molecule 1), lig_1.pdb
```

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
{"status":"ok","filename":"1ubq.pdb","total_area":4834.716264864688,"atom_areas":[17.420005600449258,16.223284994102602]}
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

For workflow-driven batch runs, use `[output.jsonl]` to also omit bulky fields
or write metadata sidecars:

```toml
[output]
format = "jsonl"

[output.jsonl]
atom_areas = false
total_area = true
decimals = 3
metadata = "sidecar"
```

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

## Experimental AlphaFold mmCIF Fast Parser

For AlphaFold-like protein mmCIF batches, `--af-model-fast` enables an
experimental parser for simple, consecutive `_atom_site` `ATOM` rows:

```bash
zsasa batch afdb/ results.jsonl \
  --format=jsonl \
  --af-model-fast
```

This path preserves chain, residue, atom, element, and residue-number metadata,
including multiple chains. It falls back to the generic mmCIF parser for
unsupported layouts, alternate locations, multiple models, hydrogens, and
extended chain IDs. Malformed coordinates and I/O errors are reported instead
of being hidden by fallback. Chain filters, author-chain matching,
alternate-location overrides, and explicit hydrogen inclusion use the generic
parser directly.

The input strategy can be selected independently with
`--input-io=auto|mmap|read`. `auto` uses whole-file reads for the AF fast path
and keeps each generic parser's existing default. Explicit `mmap` or `read`
values are useful when benchmarking storage and operating-system behavior.

For parser and output profiling, combine `--timing` with `--profile-stages`.
The additional machine-readable lines report aggregate read/parse, classifier,
and JSONL write times across the batch.

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

Use `--chain` for label/asym chain IDs. Add `--auth-chain` as a boolean modifier when the `--chain` value should match author-provided chain IDs in mmCIF inputs. With `--auth-chain`, reported residue numbers are author-provided too (`auth_seq_id`).

For named multi-job runs such as chain A, chain B, and AB complex calculations, use [Workflow Files](workflows.md) instead of repeating shell commands.

When the desired chain set differs by input file, use a workflow
[`chain_map`](workflows.md#per-file-chain-maps). Chain maps accept CSV or JSON
and can choose `label` or `auth` chain IDs independently for each mmCIF or
BinaryCIF file. Multiple stable-ID rows for one filename are calculated from
one parsed/classified structure, with optional source-indexed atom identity for
joining primitive SASA selections downstream.

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

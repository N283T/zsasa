---
sidebar_position: 3
---

# Workflow Files

Workflow files keep input paths, output paths, calculation settings, and classifier settings in a reproducible TOML file.

Use workflows when command lines become long, when you need named batch jobs, or when you want settings committed alongside an analysis.

## Commands

```bash
zsasa calc --workflow sasa.toml
zsasa batch --workflow bsa.toml
```

For batch mode, `--manifest` is still accepted as a compatibility alias for `--workflow`:

```bash
zsasa batch --manifest bsa.toml
```

Prefer `--workflow` in new scripts and docs.

## Calc Workflow Example

```toml
version = 1
kind = "workflow"

[input]
path = "structure.cif"

[output]
path = "sasa.json"
format = "json"

[calculation]
algorithm = "sr"
rsa = true
n_points = 100
probe_radius = 1.4
threads = 0

[classifier]
type = "ccd"
ccd = "components.zsdc"
```

Run it with:

```bash
zsasa calc --workflow sasa.toml
```

## Batch Workflow Example

Batch workflows contain one or more named `[[jobs]]` entries:

```toml
version = 1
kind = "workflow"

[input]
dir = "structures"

[output]
dir = "results"
format = "jsonl"

[output.jsonl]
decimals = 3
atom_areas = false
total_area = true
metadata = "sidecar"

[calculation]
algorithm = "sr"
residue_map = true
n_points = 100
threads = 8

[classifier]
type = "custom"
config = "my_classifier.toml"

[[jobs]]
name = "chain_a"
chains = ["A"]

[[jobs]]
name = "chain_b"
chains = ["B"]

[[jobs]]
name = "complex_ab"
chains = ["A", "B"]
```

Run it with:

```bash
zsasa batch --workflow bsa.toml
```

For ordinary PDB, JSON, and unfiltered mmCIF/BinaryCIF workflow batch runs with compatible chain-ID settings, zsasa reuses each parsed input structure across jobs internally. For named chain analyses such as chain A, chain B, and complex AB, list only the jobs you want; eligible runs parse each input structure once and then compute each requested chain selection independently. Some inputs or settings, such as SDF files, per-job `auth_chain` changes, or mmCIF/BinaryCIF workflows with chain filters, use the compatibility job-first path instead so full chain-ID selection matches parser behavior.

## BSA / ΔSASA Analysis

Batch workflows can also write an analysis JSONL file for a two-partner interface:

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
n_points = 100
threads = 8

[classifier]
type = "ccd"

[analysis]
type = "bsa"
name = "interface_ab"
partner_a = ["A"]
partner_b = ["B"]
level = "residue"
```

Run it with:

```bash
zsasa batch --workflow bsa.toml
```

This writes `results/interface_ab.jsonl`, with one row per input structure. The initial implementation computes partner A, partner B, and the AB complex internally, then reports:

```text
delta_sasa_total = sasa_partner_a + sasa_partner_b - sasa_complex
bsa = delta_sasa_total / 2
```

`ΔSASA` and `BSA` are deliberately separate fields. `delta_sasa_total` and `residue_delta_sasa` are not halved; `bsa` is the two-partner buried surface area after the `1/2` factor.

BSA analysis JSONL uses analysis-specific fields such as `sasa_partner_a`, `sasa_partner_b`, `sasa_complex`, `delta_sasa_total`, `bsa`, and `residue_delta_sasa`. It does not use the normal SASA JSONL `total_area` and `atom_areas` fields, because those names are ambiguous for interface analysis.

## Override Precedence

When the same setting appears in multiple places, zsasa applies this order:

```text
built-in defaults < workflow settings < job settings < explicit CLI options
```

For example, this command uses the workflow but overrides the thread count:

```bash
zsasa batch structures/ results/ --workflow bsa.toml --threads=16
```

## Workflow Jobs vs CLI Chain Filters

Use CLI flags for a single ad hoc chain-filtered batch run:

```bash
zsasa batch structures/ results/ --chain=A
```

Use workflow jobs for named, repeatable multi-chain analyses:

```toml
[[jobs]]
name = "chain_a"
chains = ["A"]

[[jobs]]
name = "complex_ab"
chains = ["A", "B"]
```

## Custom Classifier Configs

Custom classifier configs are TOML-only. In batch workflows, set them in the workflow classifier section:

```toml
[classifier]
type = "custom"
config = "my_classifier.toml"
```

For single `calc` commands, you can also use the CLI option:

```bash
zsasa calc --config=my_classifier.toml structure.cif output.json
```

## Residue Maps

Set `residue_map = true` under `[calculation]` to add compact residue arrays to JSONL rows:

```toml
[output]
format = "jsonl"

[calculation]
residue_map = true
```

This is equivalent to passing `--residue-map` with `--format=jsonl` in non-workflow batch mode.

## JSONL Output Options

Workflow files can tune batch JSONL output under `[output.jsonl]`:

```toml
[output]
dir = "results"
format = "jsonl"

[output.jsonl]
atom_areas = false    # omit per-atom SASA arrays
total_area = true     # keep per-structure totals
decimals = 3          # round JSONL floating-point values
metadata = "sidecar"  # none | sidecar
```

Defaults preserve the CLI JSONL schema: `atom_areas = true`, `total_area = true`,
full precision, and no metadata sidecar. When `metadata = "sidecar"` is set,
workflow batch jobs write a `<job>.meta.json` file next to `<job>.jsonl` with
the effective JSONL and calculation settings.

## Reference

- [Batch Processing](batch.md)
- [CLI Commands](../cli/commands.md)
- [Input Formats](../cli/input.md)
- [Output Formats](../cli/output.md)

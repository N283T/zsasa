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

## Per-file Chain Maps

Use a chain map when each structure needs a different chain selection. Set
`chain_map` on a workflow job instead of `chains`:

```toml
version = 1
kind = "workflow"

[input]
dir = "structures"

[output]
dir = "results"
format = "jsonl"

[[jobs]]
name = "selected_complexes"
chain_map = "chains.csv"
```

The CSV must contain `filename` and `chains` columns. `asym_id_type` is optional
and defaults to `label`:

```csv
filename,chains,asym_id_type
1abc.pdb,"A,C",label
2xyz.cif,"A,C",auth
3def.cif,"B,D",label
```

`chains` is a comma-separated set. A row containing `A,C` calculates the SASA
of the A+C complex: atoms in both chains are included and occlude each other.
It does not calculate A and C independently. A map has one selection per
filename; to calculate multiple selections for each file, add multiple named
workflow jobs with separate chain maps.

JSON chain maps are also accepted:

```json
[
  {"filename": "1abc.pdb", "chains": ["A", "C"], "asym_id_type": "label"},
  {"filename": "2xyz.cif", "chains": ["A", "C"], "asym_id_type": "auth"}
]
```

For mmCIF and BinaryCIF, `label` matches `_atom_site.label_asym_id` and `auth`
matches `_atom_site.auth_asym_id`. PDB has one chain-ID field, so `label` and
`auth` select the same value for PDB input.

Filenames must be basenames that exactly match files in the input directory,
including their structure and compression extensions. Every discovered
structure must have one map entry; a missing entry is written as a failed batch
result. Keep a JSON chain map outside the input directory because `.json` is
also a supported structure-input extension. Duplicate filenames, empty chain
lists, and jobs that specify both
`chains` and `chain_map` are rejected. Chain maps are supported for PDB,
mmCIF, and BinaryCIF inputs, not JSON atom arrays or SDF molecule batches.

## BSA / ΔSASA Analysis {#bsa-analysis}

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

When the two partner groups differ by structure, replace `partner_a` and
`partner_b` with an analysis `chain_map`:

```toml
[analysis]
type = "bsa"
name = "interfaces"
chain_map = "interfaces.csv"
level = "residue"
```

CSV interface maps use one comma-separated chain set for each partner:

```csv
filename,id,partner_a,partner_b,asym_id_type
1abc.cif,interaction-001,"A,B","C,D",auth
1abc.cif,interaction-002,E,"C,D",auth
2xyz.cif,interaction-003,A,"B,C",label
```

The equivalent JSON is:

```json
[
  {
    "filename": "1abc.cif",
    "id": "interaction-001",
    "partner_a": ["A", "B"],
    "partner_b": ["C", "D"],
    "asym_id_type": "auth"
  },
  {
    "filename": "1abc.cif",
    "id": "interaction-002",
    "partner_a": ["E"],
    "partner_b": ["C", "D"],
    "asym_id_type": "auth"
  }
]
```

Multiple rows may name the same input file. `id` is required for every row
when a filename occurs more than once and must be unique across the map. A
legacy one-row-per-file map may omit `id`; its stable output ID defaults to the
filename. Fixed `partner_a`/`partner_b` workflows also use the filename as the
interface ID. All rows for one file must use the same `asym_id_type`, allowing
zsasa to read, decompress, parse, and classify that structure once before
processing its interfaces.

For the first row, zsasa calculates isolated A+B, isolated C+D, and the
A+B+C+D complex. The reported values are therefore:

```text
delta_sasa_total = sasa(A+B) + sasa(C+D) - sasa(A+B+C+D)
bsa = delta_sasa_total / 2
```

`asym_id_type` defaults to `label` and may be set independently for each
mmCIF or BinaryCIF file. Fixed `partner_a`/`partner_b` settings and
`analysis.chain_map` are mutually exclusive.

Run it with:

```bash
zsasa batch --workflow bsa.toml
```

This writes `results/interface_ab.jsonl`, with one row per requested
interface. Each success row has `status = "ok"` and its stable `id`. Invalid
chain selections and read, parse, classification, or calculation failures are
written as `status = "err"` rows with the same `filename` and `id`. A map row
whose source file is missing also produces an error row. A discovered file
without a map entry uses its filename as the error-row ID.

The workflow computes partner A, partner B, and the AB complex internally,
then reports:

```text
delta_sasa_total = sasa_partner_a + sasa_partner_b - sasa_complex
bsa = delta_sasa_total / 2
```

`ΔSASA` and `BSA` are deliberately separate fields. `delta_sasa_total`,
`residue_delta_sasa`, and `atom_delta_sasa` are not halved; `bsa` is the
two-partner buried surface area after the `1/2` factor.

BSA analysis JSONL uses analysis-specific fields such as `sasa_partner_a`,
`sasa_partner_b`, `sasa_complex`, `delta_sasa_total`, and `bsa`. With
`level = "residue"`, parallel residue arrays include `residue_partner`,
residue identity, `residue_sasa_isolated`, `residue_sasa_complex`, and
`residue_delta_sasa`.

Atom detail is disabled by default because it can greatly increase JSONL
volume. Enable it only with residue detail:

```toml
[analysis]
type = "bsa"
chain_map = "interfaces.csv"
level = "residue"
atom_output = true
```

Atom-detail rows include `atom_index` (zero-based in the selected interface
complex), `atom_partner`, chain, residue and atom identity, element,
`atom_sasa_isolated`, `atom_sasa_complex`, and
`atom_delta_sasa`. Atom and residue ΔSASA arrays reconcile to the unhalved
`delta_sasa_total` within floating-point tolerance before optional JSONL
decimal rounding. Polar/apolar decomposition
is not emitted; the atom metadata is available for downstream classification.
The schema does not use normal SASA JSONL `total_area` and `atom_areas` fields,
because those names are ambiguous for interface analysis.

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

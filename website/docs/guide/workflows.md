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

The CSV must contain `filename` and `chains` columns. Add a globally unique
`id` for each requested selection when a filename appears more than once.
`asym_id_type` is optional and defaults to `label`:

```csv
filename,id,chains,asym_id_type
1111.cif,a,A,label
1111.cif,bc,"B,C",label
1111.cif,abc,"A,B,C",label
2222.cif,a_2222,X,label
2222.cif,b_2222,Y,label
2222.cif,complex_2222,"X,Y",label
```

`chains` is a comma-separated set. A row containing `A,C` calculates the SASA
of the A+C complex: atoms in both chains are included and occlude each other.
It does not calculate A and C independently. The first three rows above emit
the reusable primitive results SASA(A), SASA(B+C), and SASA(A+B+C). This path
does not calculate ΔSASA or BSA; derive those quantities downstream.

JSON chain maps are also accepted:

```json
[
  {"filename": "1111.cif", "id": "a", "chains": ["A"], "asym_id_type": "label"},
  {"filename": "1111.cif", "id": "bc", "chains": ["B", "C"], "asym_id_type": "label"},
  {"filename": "1111.cif", "id": "abc", "chains": ["A", "B", "C"], "asym_id_type": "label"}
]
```

For a repeated filename, every row must have an `id`, IDs must be unique across
the whole map, and all rows must use the same `asym_id_type`. Legacy maps with
one row per filename may omit `id`; the output ID then defaults to the
filename, preserving one result per structure.

zsasa groups map rows by filename. It reads/decompresses, parses, and classifies
each discovered structure once, then calculates its selections from that
prepared input. Chain sets are canonicalized as sets, so requests such as
`["B", "C"]` and `["C", "B"]` share one SASA calculation while still emitting
separate rows with their requested IDs and chain order.

Workflow `threads` controls concurrent structure-file workers. With more than
one file worker, every individual SASA calculation uses one internal thread to
avoid nested oversubscription. With one file, the configured threads are
available to its SASA calculations. Progress advances once per discovered
structure after all of its selections finish.

Multi-selection maps require JSONL output. Each requested selection produces a
complete, non-interleaved success or error row:

```json
{"status":"ok","filename":"1111.cif","id":"bc","chains":["B","C"],"total_area":1234.5,"atom_areas":[12.3,0]}
{"status":"err","filename":"1111.cif","id":"missing","chains":["Z"],"error":"selected chain not found: Z"}
```

Success rows always include `filename`, `id`, `chains`, and `total_area`.
`atom_areas` follows `[output.jsonl].atom_areas`, and residue arrays follow
`[calculation].residue_map`. Missing map rows, missing input structures,
missing selected chains, and per-selection calculation failures are emitted as
error rows without discarding other selections for the same structure.

For stable atom joins across selections, enable identity metadata together with
atom areas:

```toml
[output.jsonl]
atom_areas = true
atom_identity = true
```

This adds parallel `source_atom_index`, `atom_chain`, `atom_residue_name`,
`atom_residue_number`, `atom_insertion_code`, `atom_name`, and `atom_element`
arrays. `source_atom_index` is zero-based in the once-parsed source atom array
after configured model, alternate-location, hydrogen, and HETATM filtering. It
therefore aligns atoms across A, B+C, and A+B+C rows from the same source and
settings. `atom_element` supports downstream polar/apolar classification
without adding this metadata when atom output is disabled. `atom_identity`
currently applies only to `chain_map` jobs and requires `atom_areas = true`.

For mmCIF and BinaryCIF, `label` matches `_atom_site.label_asym_id` and `auth`
matches `_atom_site.auth_asym_id`. PDB has one chain-ID field, so `label` and
`auth` select the same value for PDB input.

Filenames must be basenames that exactly match files in the input directory,
including their structure and compression extensions. Every discovered
structure must have one map entry; a missing entry is written as a failed batch
result. Keep a JSON chain map outside the input directory because `.json` is
also a supported structure-input extension. Missing/duplicate IDs, mixed
`asym_id_type` values for one filename, empty chain lists, and jobs that specify
both `chains` and `chain_map` are rejected. Chain maps are supported for PDB,
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

For BSA workflows, `threads` and `--threads=N` control the number of structure
files processed concurrently. Explicit values may exceed the detected CPU
count, which is useful for large I/O-bound datasets. In the multi-file worker
path, each partner A, partner B, and complex SASA calculation uses one internal
thread to avoid nested oversubscription; all interfaces for a structure remain
serial and reuse the same parsed and classified input. When only one file or
one worker is available, the configured thread count is instead available to
the SASA calculations within that file.

BSA workflows show standard terminal progress by completed structure file, not
by interface row. Each discovered file advances progress once after all of its
interfaces finish, including files that produce parse, classification,
selection, or calculation errors. Set `quiet = true` or pass `--quiet` to
suppress progress. The standard progress display requires stderr to be a
terminal: redirecting stderr, including with `2>&1 | tee`, disables the
interactive bar, while piping stdout alone leaves the bar on the terminal but
does not capture it in `tee`.

BSA JSONL rows from different files are streamed in completion order, so global
input order is not guaranteed. Rows are written atomically, and interface order
within each file follows the interface-map order.

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
atom_identity = false # chain-map-only stable atom identity arrays
total_area = true     # keep per-structure totals
decimals = 3          # round JSONL floating-point values
metadata = "sidecar"  # none | sidecar
```

Defaults preserve the CLI JSONL schema: `atom_areas = true`,
`atom_identity = false`, `total_area = true`, full precision, and no metadata
sidecar. `atom_identity` is available only for chain-map jobs and requires atom
areas. Selection-map JSONL also requires `total_area = true`. When
`metadata = "sidecar"` is set, workflow batch jobs write a `<job>.meta.json`
file next to `<job>.jsonl` with the effective JSONL and calculation settings.

## Reference

- [Batch Processing](batch.md)
- [CLI Commands](../cli/commands.md)
- [Input Formats](../cli/input.md)
- [Output Formats](../cli/output.md)

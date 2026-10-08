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

Prefer `--workflow` in new scripts and docs. `--manifest` also reads the older [flat-root manifest layout](#legacy-flat-root-manifests).

Every key a workflow file accepts is listed in the [Key Reference](#key-reference).

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

## Alternate Locations {#alternate-locations}

Set `altloc` under `[calculation]` to choose which alternate locations are used. It takes the values of the `--altloc` option: `"auto"` (the default), `"none"`, `"all"`, `"highest-occupancy"`, or one altloc ID such as `"B"`.

```toml
[calculation]
altloc = "highest-occupancy"
```

The key applies to `calc --workflow` and to every batch workflow: jobs, chain maps and BSA analysis. `--altloc` on the command line takes precedence over the key, and a value that is not one of the above is rejected when the workflow is read. The modes are described in [Alternate Locations](../cli/input.md#alternate-locations).

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

For multi-selection maps, zsasa schedules discovered structures by descending
selection count before file workers claim them. Filename order breaks equal-cost
ties deterministically. This longest-processing-time-first estimate starts
selection-heavy structures early to reduce the file-level tail without reading
or parsing inputs during scheduling. Generic batch scans, BSA workflows, and
legacy one-row-per-file maps keep their existing order. Selection calculations
within a structure remain serial, so exceptionally heavy final structures can
still leave residual tail latency; bounded parse-once work stealing is tracked
in [GitHub issue #413](https://github.com/N283T/zsasa/issues/413).

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
`auth` select the same value for PDB input. The choice also sets the residue
numbers in residue-level and atom-level output: `label_seq_id` with `label`
(`auth_seq_id` for waters and ligands, which have no `label_seq_id`) and
`auth_seq_id` with `auth`.

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

`--chain` belongs to this non-workflow form: it cannot be combined with `--workflow` (the command stops with `--workflow cannot be combined with --chain`), so put chain selections in the jobs instead.

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

These keys choose the fields of each row, not where the rows go: with
`atom_areas = false` a job still writes one row per input to `<job>.jsonl`, or
to standard output when the workflow has one job and no output directory.

## Key Reference

A workflow file starts with `version = 1` (required) and an optional `kind = "workflow"`, followed by the sections below. The parser rejects the whole file when it meets an unknown section or key, a repeated section or key, or a value of the wrong type or out of range, and prints `Error reading workflow file '<path>': <name>` with `UnknownField`, `InvalidFieldType`, `UnsupportedVersion` or `InvalidKind`.

The **Used by** column says which command reads the key. `calc` and `batch` accept the keys of the other command and ignore them, so one file can serve both. A value given on the command line overrides the same setting in the file; see [Override Precedence](#override-precedence).

### Top level

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `version` | integer | required | calc, batch | Must be `1` |
| `kind` | string | none | calc, batch | If present, must be `"workflow"` |

### `[input]`

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `path` | string | none | calc | Input structure file. A positional input argument overrides it |
| `dir` | string | none | batch | Input directory. A positional `input_dir` overrides it |
| `chain` | string | all chains | calc | Chain filter such as `"A"` or `"A,B"` (`--chain`). In batch workflows chains are chosen per job with `chains` or `chain_map` |
| `model` | integer ≥ 1 | all models | calc | Model number (`--model`) |
| `mol` | string | first molecule | calc | Molecule title or 1-based index in a multi-molecule SDF (`--mol`) |

### `[output]` and `[output.jsonl]`

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `path` | string | `output.json` | calc | Output file. A positional output argument or `-o` overrides it |
| `dir` | string | none | batch | Output directory; each job writes into its own `<job name>` subdirectory, or to `<job name>.jsonl` for JSONL output. A positional `output_dir` overrides it. A workflow with several jobs requires it; with one job and no directory, JSONL rows go to standard output and `json`/`csv` write no files |
| `format` | string | `"json"` | calc, batch | `calc`: `json`, `compact`, `csv`, `freesasa` or `rsa`. `batch`: `json`, `compact`, `csv` or `jsonl` |

The `[output.jsonl]` keys apply to batch JSONL output only. See [JSONL Output Options](#jsonl-output-options).

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `atom_areas` | boolean | `true` | batch | Write the per-atom `atom_areas` array |
| `atom_identity` | boolean | `false` | batch | Write stable atom identity arrays; `chain_map` jobs only, requires `atom_areas = true` |
| `total_area` | boolean | `true` | batch | Write `total_area`; selection-map output requires it |
| `decimals` | integer 0-15 | full precision | batch | Round floating-point values (`--jsonl-decimals`) |
| `metadata` | string | `"none"` | batch | `"none"` or `"sidecar"` (write `<job>.meta.json`) |

### `[calculation]`

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `algorithm` | string | `"sr"` | calc, batch | `"sr"` or `"lr"` |
| `threads` | integer ≥ 0 | `0` (auto) | calc, batch | Worker threads; concurrent file workers in batch |
| `probe_radius` | number | `1.4` | calc, batch | Probe radius in Å, above 0 and at most 10 |
| `n_points` | integer | `100` | calc, batch | Test points per atom for SR, 1-10000 |
| `n_slices` | integer | `20` | calc, batch | Slices per atom diameter for LR, 1-1000 |
| `lr_trig` | string | `"exact"` | calc, batch | [Arc angles](algorithms.mdx#lee-richards-arc-angles) of LR: `"exact"`, or `"fast"` for the approximation of zsasa 0.9.1 and earlier (`--lr-trig`) |
| `precision` | string | `"f64"` | calc, batch | `"f32"` or `"f64"` |
| `include_hydrogens` | boolean | `false` | calc, batch | Include hydrogen atoms |
| `include_hetatm` | boolean | `false` | calc, batch | Include HETATM records |
| `use_bitmask` | boolean | `false` | calc, batch | Bitmask LUT optimization (SR only, `n_points` 1-1024) |
| `timing` | boolean | `false` | calc, batch | Print the timing breakdown |
| `quiet` | boolean | `false` | calc, batch | Suppress progress output |
| `auth_chain` | boolean | `false` | calc, batch | Match chains and number residues by `auth_asym_id` / `auth_seq_id` (mmCIF/BinaryCIF). A job can override it with its own `auth_chain` |
| `altloc` | string | `"auto"` | calc, batch | Alternate-location handling: `"auto"`, `"none"`, `"all"`, `"highest-occupancy"` or one altloc ID such as `"A"`; see [Alternate Locations](../cli/input.md#alternate-locations). `--altloc` on the command line takes precedence |
| `residue_map` | boolean | `false` | batch | Add residue map arrays to JSONL rows (`--residue-map`) |
| `per_residue` | boolean | `false` | calc | Per-residue aggregation |
| `rsa` | boolean | `false` | calc | Relative solvent accessibility; implies `per_residue` |
| `polar` | boolean | `false` | calc | Polar/nonpolar summary; implies `per_residue` |
| `validate_only` | boolean | `false` | calc | Validate the input without calculating (`--validate`) |

### `[classifier]`

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `type` | string | as the command line | calc, batch | `"ccd"`, `"protor"`, `"naccess"`, `"oons"` or `"custom"`. Without it the command default applies (`ccd` for structure files) |
| `config` | string | none | calc, batch | Path of a custom classifier TOML file. Required with `type = "custom"` and not allowed with any other `type` |
| `ccd` | string | none | calc, batch | External CCD dictionary (CIF, `.gz`/`.zst` or ZSDC). Only with `type = "ccd"` |
| `sdf` | string or array of strings | none | calc, batch | SDF file or files with bond topology. Only with `type = "ccd"` |

Set `type = "ccd"` whenever `ccd` or `sdf` is given. A file that combines them with another classifier type is rejected with `InvalidClassifierConfig`.

### `[analysis]`

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `type` | string | required | batch | Must be `"bsa"` |
| `name` | string | `"bsa"` | batch | Output name: results go to `<name>.jsonl`. Must not contain `/`, `\` or `..` |
| `partner_a`, `partner_b` | arrays of strings | none | batch | Chain IDs of the two partners. Required unless `chain_map` is set, and not allowed together with it |
| `chain_map` | string | none | batch | CSV or JSON file with one interface per row, instead of `partner_a`/`partner_b` |
| `level` | string | `"total"` | batch | `"total"` or `"residue"` |
| `atom_output` | boolean | `false` | batch | Add atom-level ΔSASA; requires `level = "residue"` |

See [BSA / ΔSASA Analysis](#bsa-analysis) for the output.

### `[[jobs]]`

A batch workflow needs at least one job unless it has an `[analysis]` section (`NoJobs` otherwise). `calc` ignores jobs.

| Key | Type | Default | Used by | Description |
|-----|------|---------|---------|-------------|
| `name` | string | required | batch | Job name, unique within the file; used as the output subdirectory or file name, so it must not contain `/`, `\` or `..` |
| `chains` | array of strings | all chains | batch | Chain IDs to calculate together as one complex. Leave the key out to select every chain: an empty array (`chains = []`) is rejected with `EmptyJobChains` |
| `chain_map` | string | none | batch | Per-file chain map; see [Per-file Chain Maps](#per-file-chain-maps). Not allowed together with `chains` or `auth_chain` |
| `auth_chain` | boolean | from `[calculation]` | batch | Use author chain IDs for this job |

## Legacy Flat-root Manifests

Batch also reads the older manifest layout, in which the settings sit at the root of the file instead of in sections. It is selected automatically when the root contains any of the flat keys below, and it is read by `--manifest` as well as `--workflow`:

```toml
version = 1
input_dir = "examples"
output_dir = "zig-out/workflow-smoke/legacy"
format = "jsonl"
classifier = "ccd"
n_points = 32
quiet = true

[[jobs]]
name = "all"
```

(`test_data/legacy-batch-manifest.toml` in the repository is a working copy of this file.)

| Flat key | Same as |
|----------|---------|
| `input_dir` | `[input] dir` |
| `output_dir` | `[output] dir` |
| `format` | `[output] format` |
| `classifier` | `[classifier] type` (a built-in name) |
| `ccd`, `sdf` | `[classifier] ccd`, `sdf` |
| `algorithm`, `threads`, `probe_radius`, `n_points`, `n_slices`, `precision`, `include_hydrogens`, `include_hetatm`, `use_bitmask`, `timing`, `quiet`, `auth_chain`, `residue_map` | the key of the same name in `[calculation]` |

`version` is required, `kind` is optional, and `[[jobs]]` entries work as above. A legacy manifest cannot also contain sections such as `[input]` or `[calculation]`, and it cannot use the keys that only exist in the sectioned layout (`custom` classifiers, `[output.jsonl]`, `[analysis]`, `lr_trig`, `per_residue`, `rsa`, `polar`, `validate_only`): the file is rejected (`UnknownField`, or `InvalidClassifierConfig` for `classifier = "custom"`). New files should use the sectioned layout.

## Reference

- [Batch Processing](batch.md)
- [CLI Commands](../cli/commands.md)
- [Input Formats](../cli/input.md)
- [Output Formats](../cli/output.md)

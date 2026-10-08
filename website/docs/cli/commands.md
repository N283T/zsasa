---
sidebar_position: 1
---

# Commands & Options

## Synopsis

```bash
zsasa calc <input> [output] [OPTIONS]
zsasa calc --workflow <workflow.toml>
zsasa batch <input_dir> [output_dir] [OPTIONS]
zsasa batch --workflow <workflow.toml>
zsasa traj <trajectory> <topology> [OPTIONS]
zsasa compile-dict <input.cif[.gz|.zst]> -o <output.zsdc>
```

## Workflow Files

```bash
zsasa calc --workflow <workflow.toml>
zsasa batch --workflow <workflow.toml>
```

Workflow files are TOML files that keep input, output, calculation, classifier, and batch job settings together. In batch mode, `--manifest` remains a compatibility alias for `--workflow`; it also accepts the older flat-root manifest layout.

See [Workflow Files](../guide/workflows.md) for full examples, the [key reference](../guide/workflows.md#key-reference), precedence rules, custom classifier config, and residue-map usage.

## Subcommands

### `calc` - Single File SASA Calculation

Calculate SASA for a single structure file.

```bash
zsasa calc structure.cif output.json
zsasa calc --workflow sasa.toml
```

### `batch` - Directory Batch Processing

Process all structure files in a directory.

```bash
zsasa batch structures/ results/
zsasa batch --workflow bsa.toml
```

Positional paths are required for ordinary batch mode, but can be supplied by the workflow file when using `--workflow`.

Batch mode uses file-level parallelism: multiple files are processed simultaneously, one thread per file. Use `--threads` to control the number of concurrent files.

For task-oriented examples, see [Batch Processing](../guide/batch.md) and [Workflow Files](../guide/workflows.md).

### `traj` - Trajectory Analysis

Calculate SASA for each frame in a trajectory file (XTC, TRR, DCD, or AMBER NetCDF).

```bash
zsasa traj trajectory.xtc topology.pdb --output=sasa.csv
```

See [Trajectory Options](#trajectory-options) for traj-specific options and details.

### `compile-dict` - Compile CCD Dictionary

Convert a CCD dictionary from CIF text to compact binary ZSDC format for faster loading.

```bash
zsasa compile-dict components.cif.gz -o components.zsdc
```

The compiled ZSDC file can then be used with `--ccd=components.zsdc` for faster dictionary loading compared to parsing CIF text at runtime.

## Basic Usage

```bash
# Basic SASA calculation
./zig-out/bin/zsasa calc structure.cif output.json

# With algorithm selection
./zig-out/bin/zsasa calc --algorithm=lr structure.cif output.json

# Multi-threaded
./zig-out/bin/zsasa calc --threads=4 structure.cif output.json

# With analysis features
./zig-out/bin/zsasa calc --rsa --polar structure.cif output.json

# With workflow file
./zig-out/bin/zsasa calc --workflow sasa.toml
```

## Calc-only Options

### Custom classifier config

`--config=FILE` is a `calc`-only option for advanced custom classifiers. Files must use TOML format and the `.toml` extension. If both `--classifier` and `--config` are specified for `calc`, `--config` takes precedence and a warning is printed.

```bash
zsasa calc --config=my_classifier.toml structure.cif output.json
```

For batch custom classifier configs, use a workflow `[classifier]` section. See [Workflow Files](../guide/workflows.md#custom-classifier-configs).

## Options

The three commands share most of their options but not all of them. The **Commands** column of every table below lists the commands that accept the option; a command rejects any other option with `Error: Unknown option: <option>`. `calc`, `batch` and `traj` all accept `-h, --help`.

### Algorithm Options

| Option | Description | Default | Commands |
|--------|-------------|---------|----------|
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) or `lr` (Lee-Richards) | `sr` | calc, batch, traj |
| `--precision=P` | Floating-point precision: `f32` or `f64` | `f64` (`traj`: `f32`) | calc, batch, traj |
| `--probe-radius=R` | Probe radius in Å (0 < R ≤ 10) | `1.4` | calc, batch, traj |
| `--n-points=N` | Test points per atom (SR only, 1-10000) | `100` | calc, batch, traj |
| `--n-slices=N` | Slices per atom diameter (LR only, 1-1000) | `20` | calc, batch, traj |
| `--lr-trig=MODE` | [Arc angles](../guide/algorithms.mdx#lee-richards-arc-angles) of LR: `exact` (`acos`/`atan2`) or `fast` (the polynomial approximation of zsasa 0.9.1 and earlier, totals a few tenths of a percent too high). Accepted and unused with `--algorithm=sr`, like `--n-slices` | `exact` | calc, batch, traj |
| `--use-bitmask` | Use [bitmask LUT optimization](../guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points 1-1024) | off | calc, batch, traj |
| `--bitmask-lut-mode=MODE` | Bitmask LUT reuse mode: `single`, `per-frame`, or `cycle`; non-default modes require `--use-bitmask` | `single` | traj |
| `--bitmask-correction` | Experimental exposed-fraction correction for bitmask quantization bias; requires `--use-bitmask` | off | calc, batch, traj |
| `--bitmask-correction-coeff=V` | Override the experimental correction coefficient | `0.020` | calc, batch, traj |
| `--adaptive-sr` | Experimental two-stage bitmask SR; requires `--use-bitmask` and `--algorithm=sr` | off | batch |
| `--coarse-points=N` | Coarse adaptive SR test points (bitmask range 1-1024) | `64` | batch |
| `--fine-points=N` | Fine adaptive SR test points; defaults to explicit `--n-points`, otherwise 256 | `256` | batch |
| `--adaptive-low=X` | Accept coarse result when exposed fraction is ≤ X | `0.10` | batch |
| `--adaptive-high=X` | Accept coarse result when exposed fraction is ≥ X | `0.90` | batch |
| `--threads=N` | Number of threads (0 = auto-detect). `calc` and `traj` use them for the calculation; `batch` runs that many files at once, and an explicit value may exceed the CPU count for I/O-bound file sets | `0` | calc, batch, traj |

### Classifier Options

| Option | Description | Default | Commands |
|--------|-------------|---------|----------|
| `--classifier=TYPE` | Built-in classifier: `ccd`, `protor`, `naccess`, or `oons` | calc/batch: `ccd` for PDB/mmCIF/BinaryCIF/SDF/MOL, none for JSON; traj: `naccess` | calc, batch, traj |
| `--ccd=FILE` | External CCD dictionary (CIF text, optionally `.gz`/`.zst`, or ZSDC binary); used with the `ccd` classifier | none | calc, batch, traj |
| `--sdf=PATH` | SDF file with bond topology for the `ccd` classifier; repeat the option for several ligands | none | calc, batch, traj |
| `--mol=NAME\|N` | Select one molecule of a multi-molecule SDF input by title or 1-based index | first molecule | calc |
| `--config=FILE` | Custom classifier in TOML format; see [Custom classifier config](#custom-classifier-config) | none | calc |

When `--classifier` is used, atom radii are assigned based on residue and atom names. For PDB/mmCIF/BinaryCIF/SDF/MOL input, `ccd` is used by default. HETATM records are excluded unless `--include-hetatm` is given, whichever classifier is used.

`--ccd` and `--sdf` only matter for the `ccd` classifier; with another classifier they are ignored.

See [Classifiers](../guide/classifiers.mdx) for detailed classifier documentation.

### Structure Filtering

| Option | Description | Default | Commands |
|--------|-------------|---------|----------|
| `--chain=ID` | Filter by chain ID (e.g., `A` or `A,B,C`). In `batch` it cannot be combined with `--workflow`; use the `chains` of a workflow job. `batch` rejects a value without any chain ID, such as `--chain=,` | all chains | calc, batch |
| `--model=N` | Model number for NMR structures (≥1) | all models | calc |
| `--auth-chain` | Use auth_asym_id for chain IDs and auth_seq_id for residue numbers (mmCIF/BinaryCIF) | label_asym_id, label_seq_id | calc, batch |
| `--altloc=MODE` | [Alternate-location handling](input.md#alternate-locations) for PDB, mmCIF and BinaryCIF input: `auto`, `none`, `all`, `highest-occupancy`, or one ID such as `A` | `auto` | calc, batch, traj |
| `--include-hydrogens` | Include hydrogen atoms | calc/batch: excluded; traj: included | calc, batch, traj |
| `--no-hydrogens` | Exclude hydrogen atoms (`--exclude-hydrogens` is a synonym) | — | traj |
| `--include-hetatm` | Include HETATM records | excluded | calc, batch |

`traj` always reads every `ATOM` and `HETATM` record of the first model of its topology; see [Trajectory Options](#trajectory-options).

### Analysis Options

| Option | Description | Commands |
|--------|-------------|----------|
| `--per-residue` | Show per-residue SASA aggregation | calc |
| `--rsa` | Calculate Relative Solvent Accessibility (enables `--per-residue`) | calc |
| `--polar` | Show polar/nonpolar SASA summary (enables `--per-residue`) | calc |

See [Output & Analysis](output.md#analysis-features) for detailed output descriptions.

### Workflow Options

| Option | Description | Commands |
|--------|-------------|----------|
| `--workflow=PATH` | Read settings from a TOML [workflow file](../guide/workflows.md); explicit command-line options override it | calc, batch |
| `--manifest=PATH` | Compatibility alias for `--workflow` that also accepts legacy flat-root manifests | batch |

### Batch Input Options

| Option | Description | Default | Commands |
|--------|-------------|---------|----------|
| `--input-io=MODE` | File input strategy where supported: `auto`, `mmap`, or `read` (`auto` reads whole files for the AlphaFold fast path and keeps each parser's own default otherwise) | `auto` | batch |
| `--af-model-fast` | Experimental AlphaFold-model mmCIF fast parser with multi-chain metadata and a safe fallback to the generic parser; see [Batch Processing](../guide/batch.md#experimental-alphafold-mmcif-fast-parser) | off | batch |

### Output Options

| Option | Description | Default | Commands |
|--------|-------------|---------|----------|
| `-o, --output=FILE` | Output location, as an alternative to the positional argument: `calc` writes this file; `batch` writes per-file results into this directory, or all rows to this file with `--format=jsonl`; `traj` writes this CSV file. The option wins over a positional path. `-o FILE` and `--output FILE` also work, `-o=FILE` only for `traj` | calc: `output.json`; batch: no per-file output; traj: `traj_sasa.csv` | calc, batch, traj |
| `--format=FMT` | Output format. `calc`: `json`, `compact`, `csv`, `freesasa`, `rsa`; `batch`: `json`, `compact`, `csv`, `jsonl` | `json` | calc, batch |
| `--residue-map` | Add compact residue map arrays to batch JSONL output (`--format=jsonl` only) | off | batch |
| `--jsonl-decimals=N` | Round JSONL floating-point values to `N` decimal places (`0..15`) | full precision | batch |
| `--timing` | Show timing breakdown (for benchmarking) | off | calc, batch |
| `--profile-stages` | With `--timing`, add read/parse, classifier, and JSONL write timings to the batch summary | off | batch |
| `-q, --quiet` | Suppress progress output, including standard progress bars shown by `batch` and `traj`, and the `batch` summary. Errors are not suppressed: `batch` still lists the [inputs that failed](../guide/batch.md#failed-inputs) on standard error | off | calc, batch, traj |
| `--validate` | Validate input only, do not calculate | off | calc |

`batch --format=jsonl` without an output path writes the rows to standard output.

Prefer `json` or `jsonl` for machine-readable pipelines. The `freesasa` and
`rsa` formats are legacy compatibility text formats; `rsa` warns if values do
not fit legacy NACCESS fixed-width columns.

### Information Options

| Option | Description | Commands |
|--------|-------------|----------|
| `-h, --help` | Show help message (`zsasa <command> --help` for a command) | all commands |
| `-V, --version` | Show version (`zsasa --version`) | top level only |

---

## Trajectory Options

The `traj` subcommand has additional options specific to trajectory processing.

### Arguments

| Argument | Description |
|----------|-------------|
| `<trajectory>` | Trajectory file (`.xtc`/`.trr` for GROMACS, `.dcd` for NAMD/CHARMM, `.nc`/`.ncdf` for AMBER NetCDF) |
| `<topology>` | Topology file (PDB or mmCIF) for atom names and radii. It must list the atoms of the trajectory in the same order; see [Notes](#notes) |

### Options

Every option marked `traj` in the [algorithm](#algorithm-options), [classifier](#classifier-options), [structure filtering](#structure-filtering) and [output](#output-options) tables applies, plus the trajectory-specific options below. For trajectory calculations, `--classifier` defaults to `naccess` (not `ccd`); `traj` has no custom config CLI option, so use a built-in classifier. The structure filters `--chain`, `--model`, `--auth-chain` and `--include-hetatm`, the analysis options, `--format`, `--timing` and `--validate` are not available (`Error: Unknown option`): `traj` writes the total SASA of all topology atoms per frame, and the output file is always given with `-o`/`--output` (there is no positional output argument).

| Option | Description | Default |
|--------|-------------|---------|
| `--precision=P` | Floating-point precision: `f32` or `f64` | `f32` (note: different from calc/batch) |
| `--no-hydrogens`, `--exclude-hydrogens` | Exclude hydrogen atoms from the calculation. They stay in the topology and trajectory files and are skipped in every frame | included |
| `--include-hydrogens` | Include hydrogen atoms (default, for backward compat) | included |
| `--altloc=MODE` | [Alternate-location handling](input.md#alternate-locations) for the topology (PDB or mmCIF): `auto`, `none`, `all`, `highest-occupancy`, or one ID such as `A` | `auto` |
| `--stride=N` | Process every Nth frame (N ≥ 1) | `1` |
| `--start=N` | Start from frame N | `0` |
| `--end=N` | End at frame N (inclusive) | all |
| `--batch-size=N` | Frames per batch for parallel processing (omit for auto, `threads * 2`) | auto |
| `--bitmask-lut-mode=MODE` | Bitmask LUT reuse mode (`single`, `per-frame`, `cycle`); see [Algorithms](../guide/algorithms.mdx#trajectory-lut-reuse-modes) | `single` |
| `-o FILE`, `--output=FILE` | Output CSV file (`-o=FILE` and `--output FILE` are also accepted) | `traj_sasa.csv` |

### Notes

- Supported trajectory formats: **XTC** and **TRR** (GROMACS), **DCD** (NAMD/CHARMM), and **AMBER NetCDF** (`.nc`, `.ncdf`), auto-detected from extension
- Coordinates are normalized to Å internally by the ztraj readers before SASA calculation.
- Hydrogen atoms are **included** by default in trajectory mode; use `--no-hydrogens` to exclude them. The option removes hydrogens from the topology and from every frame, so the trajectory itself must still contain them
- Topology file provides atom names for radius classification. Every `ATOM` and `HETATM` record of its **first model** is read, in file order; later models of an NMR ensemble are ignored. Solvent, ions and ligands in the topology are therefore part of the calculation
- The number of atoms in the trajectory must match the topology (before hydrogens are removed)
- Default precision is `f32` (faster for trajectory processing). `--algorithm` and `--precision` are independent: all four combinations are available
- Options are validated before the output file is created, so a rejected command leaves an existing output file untouched
- If reading or calculating a frame fails, the command exits with an error after writing the frames completed before it

---

## Examples

### Basic Calculations

```bash
# mmCIF input
./zig-out/bin/zsasa calc structure.cif output.json

# PDB input
./zig-out/bin/zsasa calc structure.pdb output.json

# JSON input
./zig-out/bin/zsasa calc atoms.json output.json
```

### Algorithm Selection

```bash
# Lee-Richards with 50 slices
./zig-out/bin/zsasa calc --algorithm=lr --n-slices=50 structure.cif output.json

# Lee-Richards with the approximate arc angles of zsasa 0.9.1 and earlier
./zig-out/bin/zsasa calc --algorithm=lr --lr-trig=fast structure.cif output.json

# Shrake-Rupley with 200 test points
./zig-out/bin/zsasa calc --algorithm=sr --n-points=200 structure.cif output.json
```

### Performance Tuning

```bash
# Fast mode: f32 precision
./zig-out/bin/zsasa calc --precision=f32 structure.cif output.json

# Multi-threaded
./zig-out/bin/zsasa calc --threads=4 structure.cif output.json

# Show timing breakdown
./zig-out/bin/zsasa calc --timing structure.cif output.json
```

### Classifier Usage

```bash
# NACCESS classifier
./zig-out/bin/zsasa calc --classifier=naccess structure.cif output.json

# Custom config
./zig-out/bin/zsasa calc --config=custom.toml structure.cif output.json
```

### Small Molecules

```bash
# Second molecule of a multi-molecule SDF file
./zig-out/bin/zsasa calc --mol=2 molecules.sdf output.json
```

### Chain/Model Filtering

```bash
# Single chain
./zig-out/bin/zsasa calc --chain=A structure.cif output.json

# Multiple chains
./zig-out/bin/zsasa calc --chain=A,B,C structure.cif output.json

# Specific model (NMR)
./zig-out/bin/zsasa calc --model=1 nmr_structure.cif output.json

# Use auth chain IDs
./zig-out/bin/zsasa calc --auth-chain --chain=A structure.cif output.json
```

### Atom Filtering

By default, hydrogen atoms and HETATM records (ligands, ions, waters) are excluded, whichever classifier is used.

```bash
# Include hydrogen atoms
./zig-out/bin/zsasa calc --include-hydrogens structure.pdb output.json

# Include HETATM records (water, ligands, etc.)
./zig-out/bin/zsasa calc --include-hetatm structure.pdb output.json

# Include both
./zig-out/bin/zsasa calc --include-hydrogens --include-hetatm structure.pdb output.json
```

### Batch Processing

```bash
# Basic batch processing
./zig-out/bin/zsasa batch input_dir/ output_dir/

# Multi-threaded (file-level parallelism)
./zig-out/bin/zsasa batch --threads=8 input_dir/ output_dir/

# One JSONL file for all structures (without -o the rows go to standard output)
./zig-out/bin/zsasa batch --format=jsonl -o results.jsonl input_dir/

# With workflow file
./zig-out/bin/zsasa batch --workflow bsa.toml
```

BSA workflow interface maps may contain multiple stable-ID rows for one input
filename. The command writes one JSONL success or error row per requested
interface; residue components and opt-in atom detail are configured in the
workflow. See [Workflow Files](../guide/workflows.md#bsa-analysis).

Ordinary SASA workflow `chain_map` files may also contain multiple stable-ID
rows per filename. zsasa parses and classifies each source once, reuses
identical chain-set calculations, and emits primitive SASA JSONL rows for
downstream analysis. See
[Per-file Chain Maps](../guide/workflows.md#per-file-chain-maps).

### Trajectory Analysis

```bash
# Basic trajectory analysis
zsasa traj trajectory.xtc topology.pdb

# With NACCESS classifier
zsasa traj trajectory.xtc topology.pdb --classifier=naccess

# Every 10th frame
zsasa traj trajectory.xtc topology.pdb --stride=10

# Frames 100-200 only
zsasa traj trajectory.xtc topology.pdb --start=100 --end=200

# Write the CSV to a chosen file
zsasa traj trajectory.xtc topology.pdb -o sasa.csv

# Lee-Richards algorithm with f64 precision
zsasa traj trajectory.xtc topology.pdb --algorithm=lr --precision=f64
```

### Validation Only

```bash
# Validate without calculation
./zig-out/bin/zsasa calc --validate structure.cif
```

---

## Error Messages

Errors are written to standard error and the command exits with status 1. `<...>` marks a value taken from the command line or the input.

### Command line

| Message | Command | Description |
|---------|---------|-------------|
| `Error: Missing input file` | calc | No input file specified (the usage line follows) |
| `Error: Missing input directory` | batch | No input directory specified |
| `Error: Missing trajectory file`, `Error: Missing topology file` | traj | A required argument is missing |
| `Error: Unknown option: <option>` | all | The command does not accept the option; see the **Commands** columns above |
| `Error: Probe radius must be between 0 and 10 Angstroms: <value>` | all | Probe radius out of range (0, 10] |
| `Error: n-points must be between 1 and 10000: <value>` | all | Test points out of range [1, 10000] |
| `Error: n-slices must be between 1 and 1000: <value>` | all | Slices out of range [1, 1000] |
| `Error: Invalid probe radius: <value>`, `Error: Invalid n-points: <value>`, `Error: Invalid thread count: <value>` | all | The value is not a number |
| `Error: Invalid format: <value>` | calc, batch | Unknown output format; the valid formats for the command follow |
| `Error: Invalid algorithm: <value>` | all | Unknown algorithm name (`calc` and `batch` also list the valid names) |
| `Error: Invalid classifier: <value>` | all | Unknown classifier name |
| `Error: Invalid precision: <value>` | all | Precision is not `f32` or `f64` |
| `Error: Invalid lr-trig: <value>` | all | `--lr-trig` is not `exact` or `fast` |
| `Error: Model number must be >= 1` | calc | `--model` is 0 or negative |
| `Error: Stride must be >= 1: <value>` | traj | `--stride=0` |
| `Error: --use-bitmask requires --n-points to be 1..1024 (got <value>)` | calc | Bitmask mode supports 1 to 1024 test points. `batch` reports `Error: UnsupportedNPoints` and `traj` `Error: --use-bitmask requires --n-points=1..1024` |
| `Error: --use-bitmask is only supported with the sr (shrake-rupley) algorithm` | calc | Bitmask mode with `--algorithm=lr`. `batch` reports `Error: BitmaskRequiresSR` and `traj` `Error: --use-bitmask requires --algorithm=sr` |
| `Error: --residue-map is only supported with --format=jsonl` | batch | `--residue-map` without JSONL output |
| `Error: --mol=<value> out of range (SDF has <n> molecules, use 1-based index)` | calc | `--mol` index past the end of the SDF file |
| `Error loading config file '<path>': custom classifier configs are TOML-only; ...` | calc | `--config` file without the `.toml` extension |
| `Error: --chain needs at least one chain ID (for example --chain=A or --chain=A,B), got '<value>'` | batch | `--chain` with an empty list, such as `--chain=,` |
| `Error reading workflow file '<path>': <name>` | calc, batch | The workflow file is missing or invalid; `<name>` is `FileNotFound`, `UnknownField`, `UnsupportedVersion`, `InvalidKind`, `InvalidFieldType`, `InvalidClassifierConfig`, `InvalidAnalysisConfig`, `MissingJobName`, `DuplicateJobName`, `UnsafeJobName`, `NoJobs` or `EmptyJobChains` (a job with `chains = []`; an explanation follows on the next line) |

### Input files

| Message | Command | Description |
|---------|---------|-------------|
| `Error reading input file '<path>': FileNotFound` | calc | Input file not found or not readable. Other names appear for other problems, for example `SyntaxError` (not a valid structure file), `ArrayLengthMismatch` (JSON arrays have different lengths) or `EmptyInput` (no atoms) |
| `Input validation failed with <n> error(s):` | calc | Followed by one line per problem, such as `Atom 1: Radius must be positive (value: 0)`, `Atom 0: Radius too large (max 100 Angstroms) (value: 101)` or `Atom 1: Invalid x coordinate (NaN or Inf) (value: inf)` |
| `Error: Classifier requires 'residue' and 'atom_name' fields in input` | calc | `--classifier` with JSON input that lacks classification info |
| `Error: --format=rsa requires chain, residue name, residue number, and insertion code metadata; ...` | calc | `--format=rsa` with JSON input |
| `Error: Unknown trajectory format. Supported: .xtc, .trr, .dcd, .nc, .ncdf` | traj | The trajectory extension is not recognized |
| `Error: Atom count mismatch - trajectory has <n> atoms, topology has <m>` | traj | The trajectory and the topology (first model, `ATOM` and `HETATM` records) have different numbers of atoms |
| `Error: 1 output name is shared by more than one input:` | batch | Several inputs would write the same output file (`<n> output names are ...` for several names); see [Batch Processing](../guide/batch.md#basic-directory-batch) |
| `Error running workflow job '<job>': <cause>` | batch | A workflow job could not run: `cannot read input directory '<dir>': <name>`, `its inputs share output names (listed above)`, `cannot create output directory '<dir>': <name>` or `cannot create JSONL output '<path>': <name>`. The other jobs still run, `<n> of <m> jobs failed: <jobs>` follows the totals, and the command ends with `Error: WorkflowJobFailed`; see [Workflow Files](../guide/workflows.md#failures) |
| `<n> of <m> inputs failed:` | batch | Not an error exit: inputs that could not be processed, one per line with the reason, printed at the end of the run also with `--quiet`. A workflow prints `Job '<job>': <n> of <m> inputs failed:` for each job; see [Batch Processing](../guide/batch.md#failed-inputs) |

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success. For `batch` this means that the run got through its input directory: inputs that failed are [listed on standard error](../guide/batch.md#failed-inputs) (also with `--quiet`) and have `status: "err"` rows in JSONL output, and they do not change the exit status |
| 1 | Error (invalid input, file not found, etc.). For `batch` also a run that could not start or write its output, and a [workflow job that could not run](../guide/workflows.md#failures) |

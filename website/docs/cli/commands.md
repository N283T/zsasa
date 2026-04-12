---
sidebar_position: 1
---

# Commands & Options

## Synopsis

```
zsasa calc <input> [output] [OPTIONS]
zsasa batch <input_dir> [output_dir] [OPTIONS]
zsasa traj <trajectory> <topology> [output] [OPTIONS]
zsasa compile-dict <input.cif[.gz]> -o <output.zsdc>
```

## Subcommands

### `calc` - Single File SASA Calculation

Calculate SASA for a single structure file.

```bash
zsasa calc structure.cif output.json [OPTIONS]
```

### `batch` - Directory Batch Processing

Process all structure files in a directory.

```bash
zsasa batch input_dir/ output_dir/ [OPTIONS]
```

Batch mode uses file-level parallelism: multiple files are processed simultaneously, one thread per file. Use `--threads` to control the number of concurrent files.

### `traj` - Trajectory Analysis

Calculate SASA for each frame in a trajectory file (XTC or DCD).

```bash
zsasa traj trajectory.xtc topology.pdb [OPTIONS]
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
```

## Common Options

### Algorithm Options

| Option | Description | Default |
|--------|-------------|---------|
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) or `lr` (Lee-Richards) | `sr` |
| `--precision=P` | Floating-point precision: `f32` or `f64` | `f64` |
| `--probe-radius=R` | Probe radius in Å (0 < R ≤ 10) | `1.4` |
| `--n-points=N` | Test points per atom (SR only, 1-10000) | `100` |
| `--n-slices=N` | Slices per atom diameter (LR only, 1-1000) | `20` |
| `--use-bitmask` | Use [bitmask LUT optimization](../guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points 1-1024) | off |
| `--threads=N` | Number of threads (0 = auto-detect) | `0` |

### Classifier Options

| Option | Description | Default |
|--------|-------------|---------|
| `--classifier=TYPE` | Built-in classifier: `naccess`, `protor`, `oons`, or `ccd` | `protor` for PDB/mmCIF, none for JSON |
| `--config=FILE` | Custom classifier config file (TOML or FreeSASA format, auto-detected by extension) | none |
| `--ccd=FILE` | External CCD dictionary (CIF text or ZSDC binary) | none |

When `--classifier` is used, atom radii are assigned based on residue and atom names. For PDB/mmCIF input, `protor` is used by default (matching FreeSASA/RustSASA defaults). If both `--classifier` and `--config` are specified, `--config` takes precedence. When `--classifier=ccd` is used, HETATM records are included automatically without needing `--include-hetatm`.

See [Classifiers](../guide/classifiers.mdx) for detailed classifier documentation.

### Structure Filtering

| Option | Description | Default |
|--------|-------------|---------|
| `--chain=ID` | Filter by chain ID (e.g., `A` or `A,B,C`) | all chains |
| `--model=N` | Model number for NMR structures (≥1) | all models |
| `--auth-chain` | Use auth_asym_id instead of label_asym_id | label_asym_id |
| `--include-hydrogens` | Include hydrogen atoms (calc/batch default: excluded) | excluded |
| `--no-hydrogens` | Exclude hydrogen atoms (traj default: included) | — |
| `--include-hetatm` | Include HETATM records | excluded |

### Analysis Options

| Option | Description |
|--------|-------------|
| `--per-residue` | Show per-residue SASA aggregation |
| `--rsa` | Calculate Relative Solvent Accessibility (enables `--per-residue`) |
| `--polar` | Show polar/nonpolar SASA summary (enables `--per-residue`) |

See [Output & Analysis](output.md#analysis-features) for detailed output descriptions.

### Output Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output=FILE` | Output file path | `output.json` |
| `--format=FMT` | Output format: `json`, `compact`, `csv` | `json` |
| `--timing` | Show timing breakdown (for benchmarking) | off |
| `-q, --quiet` | Suppress progress output | off |
| `--validate` | Validate input only, do not calculate | off |

### Information Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Show help message |
| `-V, --version` | Show version |

---

## Trajectory Options

The `traj` subcommand has additional options specific to trajectory processing.

### Arguments

| Argument | Description |
|----------|-------------|
| `<trajectory>` | Trajectory file (`.xtc` for GROMACS, `.dcd` for NAMD/CHARMM) |
| `<topology>` | Topology file (PDB or mmCIF) for atom names and radii |

### Options

All [common options](#common-options) apply, plus:

| Option | Description | Default |
|--------|-------------|---------|
| `--precision=P` | Floating-point precision: `f32` or `f64` | `f32` (note: different from calc/batch) |
| `--no-hydrogens` | Exclude hydrogen atoms | included |
| `--include-hydrogens` | Include hydrogen atoms (default, for backward compat) | included |
| `--stride=N` | Process every Nth frame | `1` |
| `--start=N` | Start from frame N | `0` |
| `--end=N` | End at frame N | all |
| `--batch-size=N` | Frames per batch for parallel processing (omit for auto) | auto |
| `-o, --output=FILE` | Output CSV file | `traj_sasa.csv` |

### Notes

- Supported trajectory formats: **XTC** (GROMACS) and **DCD** (NAMD/CHARMM), auto-detected from extension
- XTC coordinates are in nm; automatically converted to Å. DCD coordinates are already in Å.
- Hydrogen atoms are **included** by default in trajectory mode; use `--no-hydrogens` to exclude them
- Topology file provides atom names for radius classification
- The number of atoms in XTC must match the topology
- Default precision is `f32` (faster for trajectory processing)

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
./zig-out/bin/zsasa calc --config=custom.config structure.cif output.json
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

By default, hydrogen atoms and HETATM records are excluded (matching FreeSASA/RustSASA defaults).

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
```

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

| Error | Description |
|-------|-------------|
| `Missing input file` | No input file specified |
| `Cannot access '<path>'` | Input file not found or not readable |
| `Invalid probe radius` | Probe radius out of range (0, 10] |
| `Invalid n-points` | Test points out of range [1, 10000] |
| `Invalid n-slices` | Slices out of range [1, 1000] |
| `Invalid format` | Unknown output format |
| `Invalid algorithm` | Unknown algorithm name |
| `Invalid classifier` | Unknown classifier name |
| `Invalid model number` | Model number must be ≥ 1 |
| `Array lengths do not match` | JSON arrays have different lengths |
| `Radius must be positive` | Radius ≤ 0 in input |
| `Coordinate is not finite` | NaN or Inf in coordinates |
| `Classifier requires 'residue' and 'atom_name'` | Missing classification info |

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid input, file not found, etc.) |

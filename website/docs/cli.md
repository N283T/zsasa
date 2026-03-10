# CLI Reference

This document provides complete documentation for the zsasa command-line interface.

## Synopsis

```
zsasa calc <input> [output] [OPTIONS]
zsasa batch <input_dir> [output_dir] [OPTIONS]
zsasa traj <trajectory> <topology> [output] [OPTIONS]
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

### `traj` - Trajectory Analysis

Calculate SASA for each frame in a trajectory file (XTC or DCD).

```bash
zsasa traj trajectory.xtc topology.pdb [OPTIONS]
```

See [Trajectory Mode](#trajectory-mode) section below for details.

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

## Options

### Algorithm Options

| Option | Description | Default |
|--------|-------------|---------|
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) or `lr` (Lee-Richards) | `sr` |
| `--precision=P` | Floating-point precision: `f32` or `f64` | `f64` |
| `--probe-radius=R` | Probe radius in Å (0 < R ≤ 10) | `1.4` |
| `--n-points=N` | Test points per atom (SR only, 1-10000) | `100` |
| `--n-slices=N` | Slices per atom diameter (LR only, 1-1000) | `20` |
| `--use-bitmask` | Use [bitmask LUT optimization](guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points 1-1024) | off |
| `--threads=N` | Number of threads (0 = auto-detect) | `0` |

### Classifier Options

| Option | Description |
|--------|-------------|
| `--classifier=TYPE` | Built-in classifier: `naccess`, `protor`, or `oons` | `protor` for PDB/mmCIF, none for JSON |
| `--config=FILE` | Custom classifier config file (TOML or FreeSASA format, auto-detected by extension) | none |

When `--classifier` is used, atom radii are assigned based on residue and atom names. For PDB/mmCIF input, `protor` is used by default (matching FreeSASA/RustSASA defaults). If both `--classifier` and `--config` are specified, `--config` takes precedence.

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

## Input Formats

The input format is auto-detected from the file extension:

| Extension | Format |
|-----------|--------|
| `.json`, `.json.gz` | JSON |
| `.cif`, `.mmcif`, `.CIF`, `.mmCIF` | mmCIF |
| `.pdb`, `.PDB`, `.ent`, `.ENT` | PDB |

### JSON Format

Minimal JSON input with coordinates and radii:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [4.0, 5.0, 6.0],
  "z": [7.0, 8.0, 9.0],
  "r": [1.7, 1.55, 1.52]
}
```

Extended JSON with classification info (required for `--classifier`):

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [4.0, 5.0, 6.0],
  "z": [7.0, 8.0, 9.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["N", "CA", "C"],
  "element": [7, 6, 6]
}
```

**Fields:**

| Field | Required | Description |
|-------|----------|-------------|
| `x`, `y`, `z` | Yes | Atom coordinates in Å |
| `r` | Yes | Van der Waals radii in Å |
| `residue` | For classifier | 3-letter residue code (e.g., "ALA") |
| `atom_name` | For classifier | Atom name (e.g., "CA", "N") |
| `element` | Optional | Atomic numbers (e.g., 6=C, 7=N, 8=O) |

**Validation rules:**

- All arrays must have the same length
- Arrays cannot be empty
- Coordinates must be finite (no NaN or Inf)
- Radii must be positive and ≤ 100 Å

### mmCIF Format

Standard mmCIF files are supported. The parser extracts:

- `_atom_site.Cartn_x/y/z` - Coordinates
- `_atom_site.type_symbol` - Element (for VdW radius)
- `_atom_site.label_atom_id` / `auth_atom_id` - Atom name
- `_atom_site.label_comp_id` / `auth_comp_id` - Residue name
- `_atom_site.label_asym_id` / `auth_asym_id` - Chain ID
- `_atom_site.label_seq_id` / `auth_seq_id` - Residue number
- `_atom_site.pdbx_PDB_ins_code` - Insertion code
- `_atom_site.pdbx_PDB_model_num` - Model number
- `_atom_site.label_alt_id` - Alternate location (first kept by default)

### PDB Format

Standard PDB format files are supported with ATOM and HETATM records.

---

## Output Formats

### JSON (default)

Pretty-printed JSON with 2-space indentation:

```json
{
  "total_area": 18923.28,
  "atom_areas": [
    32.47,
    0.25,
    15.82
  ]
}
```

### Compact JSON

Single-line JSON without whitespace:

```json
{"total_area":18923.28,"atom_areas":[32.47,0.25,15.82]}
```

### CSV

Basic CSV with atom index and area:

```csv
atom_index,area
0,32.470000
1,0.250000
2,15.820000
total,18923.280000
```

When input has structural info (mmCIF/PDB), rich CSV is generated:

```csv
chain,residue,resnum,atom_name,x,y,z,radius,area
A,ALA,1,N,1.000,3.000,5.000,1.500,32.470000
A,ALA,1,CA,2.000,4.000,6.000,1.700,0.250000
,,,,,,,,18923.280000
```

---

## Analysis Features

### Per-Residue Aggregation (`--per-residue`)

Groups atom SASA by residue (chain + residue number + insertion code):

```
Per-residue SASA:
Chain  Res    Num       SASA  Atoms
----- ---- ------ ---------- ------
    A  MET      1     198.52     19
    A  LYS      2     142.31     22
    A  ALA      3      45.67      5
```

### RSA Calculation (`--rsa`)

Calculates Relative Solvent Accessibility (RSA = SASA / MaxSASA).

MaxSASA reference values from Tien et al. (2013):

| Residue | MaxSASA (Å²) | Residue | MaxSASA (Å²) |
|---------|-------------|---------|-------------|
| ALA | 129.0 | LEU | 201.0 |
| ARG | 274.0 | LYS | 236.0 |
| ASN | 195.0 | MET | 224.0 |
| ASP | 193.0 | PHE | 240.0 |
| CYS | 167.0 | PRO | 159.0 |
| GLN | 225.0 | SER | 155.0 |
| GLU | 223.0 | THR | 172.0 |
| GLY | 104.0 | TRP | 285.0 |
| HIS | 224.0 | TYR | 263.0 |
| ILE | 197.0 | VAL | 174.0 |

Output with RSA values (can exceed 1.0 for exposed terminal residues):

```
Per-residue SASA with RSA:
Chain  Res    Num       SASA    RSA  Atoms
----- ---- ------ ---------- ------ ------
    A  MET      1     198.52   0.89     19
    A  LYS      2     142.31   0.60     22
    A  ALA      3      45.67   0.35      5
```

### Polar/Nonpolar Summary (`--polar`)

Classifies residues and shows SASA breakdown. Automatically enables `--per-residue`.

- **Polar**: ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR, TYR
- **Nonpolar**: ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, TRP, VAL
- **Unknown**: Non-standard residues (ligands, modified residues, etc.)

```
Polar/Nonpolar SASA:
  Polar:       2345.67 Å² ( 45.2%) - 42 residues
  Nonpolar:    2845.23 Å² ( 54.8%) - 58 residues
  Unknown:        0.00 Å² (  0.0%) -  0 residues
```

---

## Batch Mode

Process all structure files in a directory using `zsasa batch`:

```bash
# Basic batch processing
./zig-out/bin/zsasa batch input_dir/ output_dir/

# Multi-threaded (file-level parallelism: N files in parallel)
./zig-out/bin/zsasa batch --threads=8 input_dir/ output_dir/
```

Batch mode uses file-level parallelism: multiple files are processed simultaneously, one thread per file. Use `--threads` to control the number of concurrent files.

---

## Classifiers

Built-in classifiers assign atom radii based on residue and atom names.

### NACCESS (`--classifier=naccess`)

NACCESS-compatible radii. Default classifier used by many tools.

### ProtOr (`--classifier=protor`)

ProtOr radii from Tsai et al. (1999).

### OONS (`--classifier=oons`)

OONS radii from Ooi et al.

### Custom Config (`--config=FILE`)

Load a custom classifier from a config file. The format is auto-detected by extension:

- `.toml` files use TOML format
- All other extensions use FreeSASA format

#### TOML Format

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }
C_CAR = { radius = 1.76, class = "apolar" }
N     = { radius = 1.65, class = "polar" }
O     = { radius = 1.40, class = "polar" }
S     = { radius = 1.85, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"

[[atoms]]
residue = "ALA"
atom = "CB"
type = "C_ALI"
```

- `name` - Classifier name (optional, default: "custom")
- `[types]` - Define atom types with radius (angstrom) and class (`"polar"` or `"apolar"`)
- `[[atoms]]` - Map (residue, atom) pairs to defined types. Use `"ANY"` for fallback entries.

#### FreeSASA Format

```
name: my-classifier

types:
C_ALI 1.87 apolar
C_CAR 1.76 apolar
O     1.40 polar

atoms:
ANY CA  C_ALI
ALA CB  C_ALI
```

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

### Analysis Features

```bash
# Per-residue output
./zig-out/bin/zsasa calc --per-residue structure.cif output.json

# RSA calculation
./zig-out/bin/zsasa calc --rsa structure.cif output.json

# Polar/nonpolar summary
./zig-out/bin/zsasa calc --polar structure.cif output.json

# Combined analysis
./zig-out/bin/zsasa calc --rsa --polar structure.cif output.json
```

### Output Formats

```bash
# CSV output
./zig-out/bin/zsasa calc --format=csv structure.cif output.csv

# Compact JSON
./zig-out/bin/zsasa calc --format=compact structure.cif output.json
```

### Validation Only

```bash
# Validate without calculation
./zig-out/bin/zsasa calc --validate structure.cif
```

---

## Trajectory Mode

Calculate SASA for each frame in a trajectory file.

### Usage

```bash
zsasa traj <trajectory> <topology> [OPTIONS]
```

### Arguments

| Argument | Description |
|----------|-------------|
| `<trajectory>` | Trajectory file (`.xtc` for GROMACS, `.dcd` for NAMD/CHARMM) |
| `<topology>` | Topology file (PDB or mmCIF) for atom names and radii |

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) or `lr` (Lee-Richards) | `sr` |
| `--classifier=TYPE` | Built-in classifier: `naccess`, `protor`, `oons` | none |
| `--threads=N` | Number of threads (0 = auto-detect) | `0` |
| `--probe-radius=R` | Probe radius in Å | `1.4` |
| `--n-points=N` | Test points per atom (SR only) | `100` |
| `--n-slices=N` | Slices per atom diameter (LR only) | `20` |
| `--precision=P` | Floating-point precision: `f32` or `f64` | `f32` |
| `--no-hydrogens` | Exclude hydrogen atoms | included |
| `--include-hydrogens` | Include hydrogen atoms (default, for backward compat) | included |
| `--stride=N` | Process every Nth frame | `1` |
| `--start=N` | Start from frame N | `0` |
| `--end=N` | End at frame N | all |
| `--batch-size=N` | Frames per batch for parallel processing (omit for auto) | auto |
| `-o, --output=FILE` | Output CSV file | `traj_sasa.csv` |
| `-q, --quiet` | Suppress progress output | off |
| `-h, --help` | Show help message | |

### Output Format (CSV)

```csv
frame,step,time,total_sasa
0,1,1.000,1866.44
1,2,2.000,1977.96
2,3,3.000,1884.93
...
```

### Examples

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

### Notes

- Supported trajectory formats: **XTC** (GROMACS) and **DCD** (NAMD/CHARMM), auto-detected from extension
- XTC coordinates are in nm; automatically converted to Å. DCD coordinates are already in Å.
- Hydrogen atoms are **included** by default in trajectory mode; use `--no-hydrogens` to exclude them
- Topology file provides atom names for radius classification
- The number of atoms in XTC must match the topology
- Default precision is `f32` (faster for trajectory processing)

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

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid input, file not found, etc.) |

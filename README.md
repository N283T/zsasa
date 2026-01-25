# freesasa-zig

[日本語](README.ja.md) | English

Solvent Accessible Surface Area (SASA) calculator implemented in Zig, ported from [FreeSASA](https://github.com/mittinatten/freesasa).

## Overview

SASA (Solvent Accessible Surface Area) measures the surface area of a biomolecule that is accessible to solvent molecules. This implementation provides two algorithms:

- **Shrake-Rupley (SR)**: Test point method with Golden Section Spiral
- **Lee-Richards (LR)**: Slice-based method with exact arc integration

## Features

- **Two SASA algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Direct structure input**: mmCIF file format supported
- **Chain/Model selection**: Filter by chain ID, model number, or auth chain ID
- **Atom radius classifier**: NACCESS/ProtOr/OONS classifiers with CLI and library support
- **Custom config files**: FreeSASA-compatible configuration format
- **Analysis features**:
  - Per-residue SASA aggregation
  - RSA (Relative Solvent Accessibility) calculation
  - Polar/Nonpolar surface classification
- JSON input/output format with multiple output options
- Configurable parameters (test points, slices, probe radius)
- Input validation with detailed error messages
- Memory-safe implementation with explicit allocators
- **Multi-threading support** for parallel atom processing
- **SIMD optimization** using `@Vector(8, f64)` for 8-wide batch calculations
- **Fast trigonometry** using polynomial approximations for LR algorithm
- **Neighbor list optimization** for O(N) instead of O(N²) neighbor checking

## Building

Requires Zig 0.15.2 or later.

**Supported platforms**: Linux, macOS. Windows users should use [WSL](https://learn.microsoft.com/en-us/windows/wsl/) (Windows Subsystem for Linux) with the Linux build.

```bash
zig build
```

Run tests:

```bash
zig build test
```

## Usage

```bash
freesasa_zig [OPTIONS] <input> [output.json]
```

Supported input formats: JSON, mmCIF (.cif, .cif.gz)

### Examples

```bash
# Basic usage - Shrake-Rupley (default)
./zig-out/bin/freesasa_zig input.json output.json

# Lee-Richards algorithm
./zig-out/bin/freesasa_zig --algorithm=lr input.json output.json

# Lee-Richards with custom slice count
./zig-out/bin/freesasa_zig --algorithm=lr --n-slices=50 input.json output.json

# Multi-threaded (both algorithms support parallel processing)
./zig-out/bin/freesasa_zig --threads=4 input.json output.json
./zig-out/bin/freesasa_zig --algorithm=lr --threads=4 input.json output.json

# Custom Shrake-Rupley parameters
./zig-out/bin/freesasa_zig --probe-radius=1.5 --n-points=200 input.json output.json

# CSV output format
./zig-out/bin/freesasa_zig --format=csv input.json output.csv

# Compact JSON (single line)
./zig-out/bin/freesasa_zig --format=compact input.json output.json

# Validate input only
./zig-out/bin/freesasa_zig --validate input.json

# Quiet mode
./zig-out/bin/freesasa_zig --quiet input.json output.json

# Use NACCESS classifier for automatic radius assignment
./zig-out/bin/freesasa_zig --classifier=naccess input.json output.json

# Use custom config file
./zig-out/bin/freesasa_zig --config=custom.config input.json output.json

# Direct mmCIF input (auto-detects format)
./zig-out/bin/freesasa_zig structure.cif output.json
./zig-out/bin/freesasa_zig structure.cif.gz output.json

# Chain/Model selection
./zig-out/bin/freesasa_zig --chain=A structure.cif output.json
./zig-out/bin/freesasa_zig --model=1 structure.cif output.json
./zig-out/bin/freesasa_zig --auth-chain --chain=A structure.cif output.json

# Per-residue analysis
./zig-out/bin/freesasa_zig --per-residue structure.cif output.json

# RSA (Relative Solvent Accessibility)
./zig-out/bin/freesasa_zig --rsa structure.cif output.json

# Polar/Nonpolar surface analysis
./zig-out/bin/freesasa_zig --polar structure.cif output.json
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--algorithm=ALGO` | Algorithm: `sr` (Shrake-Rupley), `lr` (Lee-Richards) | sr |
| `--classifier=TYPE` | Built-in classifier: `naccess`, `protor`, `oons` | - |
| `--config=FILE` | Custom classifier config file (FreeSASA format) | - |
| `--threads=N` | Number of threads | auto-detect |
| `--probe-radius=R` | Probe radius in Angstroms (0 < R ≤ 10) | 1.4 |
| `--n-points=N` | Test points per atom (SR only, 1-10000) | 100 |
| `--n-slices=N` | Slices per atom diameter (LR only, 1-1000) | 20 |
| `--format=FORMAT` | Output format: `json`, `compact`, `csv` | json |
| `--chain=ID` | Filter by label chain ID (e.g., "A", "B") | - |
| `--model=N` | Select model number (for multi-model files) | all |
| `--auth-chain` | Use auth chain ID instead of label chain ID | - |
| `--per-residue` | Output per-residue SASA aggregation | - |
| `--rsa` | Calculate RSA (implies --per-residue) | - |
| `--polar` | Show polar/nonpolar summary (implies --per-residue) | - |
| `--validate` | Validate input only, do not calculate SASA | - |
| `-q, --quiet` | Suppress progress output | - |
| `-h, --help` | Show help message | - |
| `-V, --version` | Show version | - |

### Input Format

JSON file with atom coordinates and van der Waals radii:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52]
}
```

**Extended format** (for classifier support):

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["CA", "CB", "C"],
  "element": [6, 6, 6]
}
```

The `residue`, `atom_name`, and `element` fields are optional:
- **`residue` + `atom_name`**: Required for `--classifier` or `--config` options
- **`element`**: Atomic numbers (6=C, 7=N, 8=O, etc.) for unambiguous element identification

**mmCIF input**:

Files with `.cif` or `.cif.gz` extensions are automatically recognized:

```bash
# Direct structure file input
./zig-out/bin/freesasa_zig 1CRN.cif output.json
./zig-out/bin/freesasa_zig 1CRN.cif.gz output.json
```

When using structure files, a classifier is automatically applied (default: NACCESS).

**Validation rules:**
- All arrays must have the same length
- Input must contain at least one atom
- Coordinates must be finite (not NaN or Inf)
- Radii must be positive and ≤ 100 Å

### Output Formats

**JSON (default)** - Pretty-printed with indentation:
```json
{
  "total_area": 1234.56,
  "atom_areas": [
    10.5,
    20.3
  ]
}
```

**Compact** - Single-line JSON:
```json
{"total_area":1234.56,"atom_areas":[10.5,20.3]}
```

**CSV** - Comma-separated values:
```csv
atom_index,area
0,10.500000
1,20.300000
total,1234.560000
```

## Preparing Input Data

Use the provided Python scripts to convert structure files:

```bash
# Convert mmCIF/PDB to input JSON
./scripts/data/cif_to_json.py structure.cif output.json

# Generate reference SASA using FreeSASA (for validation)
./scripts/data/calc_reference.py structure.cif reference.json

# Run benchmark
./scripts/benchmark.py --runs 5 --threads 4
```

Requirements: Python 3.11+, gemmi, freesasa (installed automatically via PEP 723)

## Algorithms

This implementation provides two algorithms. See [docs/algorithm.md](docs/algorithm.md) for details.

### Shrake-Rupley (SR)

Test point method - fast and simple.

1. Generate test points on unit sphere (Golden Section Spiral)
2. For each atom, count exposed points not buried by neighbors
3. SASA = 4π × r² × (exposed / total)

### Lee-Richards (LR)

Slice-based method - mathematically precise.

1. Slice each atom sphere along Z-axis
2. Calculate exposed arc length on each slice
3. Integrate arc lengths to get surface area

### Algorithm Comparison

| Aspect | Shrake-Rupley | Lee-Richards |
|--------|---------------|--------------|
| Method | Test points | Slice integration |
| Precision control | `--n-points` | `--n-slices` |
| Speed (1A0Q, 4 threads) | 2.6ms | 9.9ms |
| vs FreeSASA C | 1.2x-2.3x faster | 1.1x-1.7x faster |
| Best for | Large structures, quick analysis | High precision requirements |

### Parameters

| Parameter | Default | Algorithm | Description |
|-----------|---------|-----------|-------------|
| `n_points` | 100 | SR | Test points per atom |
| `n_slices` | 20 | LR | Slices per atom diameter |
| `probe_radius` | 1.4 Å | Both | Water molecule radius |

## Validation

Tested against FreeSASA reference implementation using ProtOr classifier:

| Structure | Atoms | FreeSASA (Å²) | Zig (Å²) | Difference |
|-----------|------:|-------------:|--------:|----------:|
| 1CRN | 327 | 3,001.13 | 3,001.13 | 0.000% |
| 1UBQ | 602 | 4,834.72 | 4,834.72 | 0.000% |
| 1A0Q | 3,183 | 18,908.90 | 18,908.90 | 0.000% |
| 3HHB | 4,384 | 25,527.36 | 25,527.36 | 0.000% |
| 1AON | 58,674 | 316,879.14 | 316,879.14 | 0.000% |
| 4V6X | 237,685 | 1,325,369.25 | 1,325,369.25 | 0.000% |

Run validation: `./scripts/validate.py`

## Performance

Benchmark comparing Zig (ReleaseFast) vs FreeSASA C (native binary), SASA calculation time only (4 threads):

| Structure | Atoms | SR Zig (ms) | SR FS-C (ms) | SR Speedup | LR Zig (ms) | LR FS-C (ms) | LR Speedup |
|-----------|------:|------------:|-------------:|-----------:|------------:|-------------:|-----------:|
| 1CRN | 327 | 0.44 | 0.53 | 1.20x | 1.45 | 1.64 | 1.13x |
| 1UBQ | 602 | 0.63 | 0.88 | 1.40x | 2.19 | 2.99 | 1.37x |
| 1A0Q | 3,183 | 2.64 | 4.42 | 1.67x | 9.90 | 16.09 | 1.62x |
| 3HHB | 4,384 | 3.54 | 5.95 | 1.68x | 14.15 | 23.44 | 1.66x |
| 1AON | 58,674 | 45.61 | 98.28 | 2.15x | 182.83 | 317.41 | 1.74x |
| 4V6X | 237,685 | 189.06 | 424.53 | 2.25x | 741.86 | 1293.48 | 1.74x |

**Summary**: Both algorithms are faster than FreeSASA C. Shrake-Rupley is **1.2x-2.3x faster**, Lee-Richards is **1.1x-1.7x faster**. Speedup increases with structure size.

### Comparison with RustSASA

[RustSASA](https://github.com/maxall41/RustSASA) is a Rust implementation that claims to be "5x faster than FreeSASA". However, this claim is for **batch processing of the entire E. coli proteome** (4,391 structures), not for single-protein calculations.

According to RustSASA's own benchmark paper, single-protein performance is essentially identical to FreeSASA:

| Implementation | Single Protein (ms) | E. coli Proteome (s) |
|----------------|--------------------:|---------------------:|
| FreeSASA | 4.0 | 28.0 |
| RustSASA | 4.0 | 5.2 |
| **Speedup** | 1.0x | 5.4x |

The "5x" speedup comes from RustSASA's CLI processing multiple files in parallel, not from faster SASA algorithm.

**SASA-only comparison** (single-threaded, n_points=100, Shrake-Rupley):

| Structure | Atoms | Zig (ms) | RustSASA (ms) | FreeSASA C (ms) | vs Rust | vs FS-C |
|-----------|------:|---------:|--------------:|----------------:|--------:|--------:|
| 1CRN | 327 | 0.40 | 0.69 | 0.79 | 1.7x | 2.0x |
| 4V6X | 237,685 | 158.58 | 500.11 | 747.95 | 3.2x | 4.7x |

**Summary**: For single-protein SASA calculation, freesasa-zig is **1.7x-3.2x faster than RustSASA** and **2.0x-4.7x faster than FreeSASA C**. RustSASA only supports Shrake-Rupley algorithm.

Run benchmark: `./scripts/benchmark.py`

### Optimization Techniques

1. **Neighbor List**: O(N) neighbor lookup instead of O(N²)
2. **8-wide SIMD**: Process 8 calculations in parallel using `@Vector(8, f64)`
3. **Fast Trigonometry**: Polynomial approximations for acos/atan2 in LR algorithm (~37% faster)
4. **Multi-threading**: Parallel atom processing with work-stealing thread pool

## Python Bindings

Use freesasa-zig from Python with NumPy arrays:

```python
import numpy as np
from freesasa_zig import calculate_sasa

# Atom coordinates (N, 3) and radii (N,)
coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

# Calculate SASA
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")

# Use Lee-Richards algorithm
result_lr = calculate_sasa(coords, radii, algorithm="lr")
```

### Installation

```bash
# Build shared library
zig build -Doptimize=ReleaseFast

# Install Python package
cd python
pip install -e .
```

See [python/README.md](python/README.md) for full API documentation.

## Atom Classifier

The classifier module assigns van der Waals radii and polarity classes to atoms based on residue and atom names. Three built-in classifiers are available, plus support for custom config files.

### Available Classifiers

| Classifier | Description | C aliphatic | ANY fallback |
|------------|-------------|-------------|--------------|
| **NACCESS** | Default, NACCESS-compatible | 1.87 Å | Yes |
| **ProtOr** | Hybridization-based (Tsai et al. 1999) | 1.88 Å | No |
| **OONS** | Older FreeSASA default | 2.00 Å | Yes |

### CLI Usage

```bash
# Use built-in classifier (requires residue/atom_name in input JSON)
./zig-out/bin/freesasa_zig --classifier=naccess input.json output.json
./zig-out/bin/freesasa_zig --classifier=protor input.json output.json
./zig-out/bin/freesasa_zig --classifier=oons input.json output.json

# Use custom config file (FreeSASA format)
./zig-out/bin/freesasa_zig --config=my_radii.config input.json output.json
```

When a classifier is specified:
1. The `residue` and `atom_name` fields in the input JSON are used to look up radii
2. If lookup fails, falls back to element-based guessing
3. If `--config` and `--classifier` are both specified, `--config` takes precedence

### Library Usage

```zig
const classifier = @import("classifier.zig");
const naccess = @import("classifier_naccess.zig");
const protor = @import("classifier_protor.zig");
const oons = @import("classifier_oons.zig");

// NACCESS classifier (default)
const radius = naccess.getRadius("ALA", "CA");  // 1.87

// ProtOr classifier (hybridization-based)
const radius = protor.getRadius("ALA", "CA");   // 1.88

// OONS classifier (larger aliphatic carbon)
const radius = oons.getRadius("ALA", "CA");     // 2.00

// Element-based fallback (for unknown atoms)
const r = classifier.guessRadiusFromAtomName(" CA "); // 1.70 (carbon)
```

### Lookup Order

1. **Residue-specific**: Exact (residue, atom) match
2. **ANY fallback**: Generic atom definition (NACCESS/OONS only)
3. **Element guess**: van der Waals radius from element symbol

See [docs/classifier.md](docs/classifier.md) for details.

## Project Structure

```
freesasa-zig/
├── src/
│   ├── main.zig              # CLI entry point
│   ├── root.zig              # Library root module
│   ├── types.zig             # Data structures (Vec3, AtomInput, etc.)
│   ├── json_parser.zig       # JSON input parsing and validation
│   ├── json_writer.zig       # Output writing (JSON, CSV)
│   ├── mmcif_parser.zig      # mmCIF file parser
│   ├── analysis.zig          # Analysis features (per-residue, RSA, polar)
│   ├── classifier.zig        # Atom classifier core (types, element guessing)
│   ├── classifier_naccess.zig # NACCESS built-in classifier
│   ├── classifier_protor.zig  # ProtOr built-in classifier
│   ├── classifier_oons.zig    # OONS built-in classifier
│   ├── classifier_parser.zig  # FreeSASA config file parser
│   ├── test_points.zig       # Golden Section Spiral generation
│   ├── neighbor_list.zig     # Spatial neighbor list (O(N) lookup)
│   ├── simd.zig              # SIMD batch operations
│   ├── thread_pool.zig       # Generic thread pool for parallelization
│   ├── shrake_rupley.zig     # Shrake-Rupley algorithm
│   └── lee_richards.zig      # Lee-Richards algorithm
├── scripts/
│   ├── benchmark.py               # Unified benchmark across all structures
│   ├── validate.py                # Validate against FreeSASA references
│   └── data/                      # Data preparation scripts
│       ├── cif_to_json.py         # Structure to JSON converter
│       ├── calc_reference.py      # Reference SASA calculator
│       ├── generate.py            # Download structures and generate references
│       ├── generate_protor.py     # Generate inputs with ProtOr radii
│       └── compare_classifiers.py # Compare classifier radii
├── benchmarks/
│   ├── structures/            # Downloaded PDB structures (.cif.gz)
│   ├── inputs/                # Generated input JSONs (element-based radii)
│   ├── inputs_protor/         # Generated input JSONs (ProtOr radii)
│   └── references/            # FreeSASA reference SASA values
├── examples/
│   ├── 1A0Q.cif.gz        # Original structure file (PDB 1A0Q)
│   ├── input_1a0q.json    # Example input (converted from cif)
│   └── 1A0Q_sasa.json     # Reference SASA from FreeSASA
├── python/
│   ├── freesasa_zig/      # Python bindings package
│   │   ├── __init__.py    # Package exports
│   │   └── core.py        # ctypes bindings with NumPy
│   └── tests/             # Python tests
├── docs/                  # Technical documentation (Japanese)
│   ├── architecture.md    # Architecture overview
│   ├── algorithm.md       # Algorithm details
│   ├── optimizations.md   # Optimization techniques
│   ├── cli-io.md          # CLI and I/O details
│   └── classifier.md      # Atom classifier details
└── plans/
    └── *.md               # Implementation plans
```

## Roadmap

- [x] Phase 1: Basic Shrake-Rupley implementation
- [x] Phase 2: Neighbor list optimization (O(N²) → O(N))
- [x] Phase 3: SIMD optimization
- [x] Phase 4: Multi-threading
- [x] Phase 5: Production features (CLI, output formats, validation)
- [x] Phase 6: CI/CD pipeline
- [x] Phase 7: Fair benchmark (timing breakdown, SASA-only measurement)
- [x] Phase 8: Benchmark dataset (6 structures from tiny to xlarge)
- [x] Phase 9: Radius classifier (auto-detect atom radii)
  - [x] Core data structures and element-based guessing
  - [x] NACCESS/ProtOr/OONS built-in classifiers
  - [x] Custom config file parser (FreeSASA format)
  - [x] CLI integration (`--classifier`, `--config`)
- [x] Phase 10: Direct structure file input (mmCIF)
- [x] Phase 11: Lee-Richards algorithm (with multi-threading & SIMD)
- [x] Phase 13: Python bindings (NumPy integration via C API)
- [x] Phase 15: Chain/Model selection (`--chain`, `--model`, `--auth-chain`)
- [x] Phase 16: Analysis features
  - [x] Per-residue SASA aggregation (`--per-residue`)
  - [x] RSA calculation (`--rsa`)
  - [x] Polar/Nonpolar classification (`--polar`)

## License

MIT

## References

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - Original C implementation
- [RustSASA](https://github.com/maxall41/RustSASA) - Rust implementation (batch processing optimized)

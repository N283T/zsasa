# freesasa-zig

[日本語](README.ja.md) | English

Solvent Accessible Surface Area (SASA) calculator implemented in Zig, ported from [FreeSASA](https://github.com/mittinatten/freesasa).

## Overview

SASA (Solvent Accessible Surface Area) measures the surface area of a biomolecule that is accessible to solvent molecules. This implementation provides two algorithms:

- **Shrake-Rupley (SR)**: Test point method with Golden Section Spiral
- **Lee-Richards (LR)**: Slice-based method with exact arc integration

## Features

- **Two SASA algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Atom radius classifier**: NACCESS/ProtOr/OONS classifiers with CLI and library support
- **Custom config files**: FreeSASA-compatible configuration format
- JSON input/output format with multiple output options
- Configurable parameters (test points, slices, probe radius)
- Input validation with detailed error messages
- Memory-safe implementation with explicit allocators
- **Multi-threading support** for parallel atom processing
- **SIMD optimization** using `@Vector(4, f64)` for batch calculations
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
freesasa_zig [OPTIONS] <input.json> [output.json]
```

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
./scripts/cif_to_input_json.py structure.cif output.json

# Generate reference SASA using FreeSASA (for validation)
./scripts/calc_reference_sasa.py structure.cif reference.json

# Run benchmark
./scripts/benchmark.py examples/1A0Q.cif.gz --runs 5 --threads 4
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
| Speed (1A0Q, 4 threads) | 13ms | 26ms |
| Best for | Large structures, quick analysis | High precision requirements |

### Parameters

| Parameter | Default | Algorithm | Description |
|-----------|---------|-----------|-------------|
| `n_points` | 100 | SR | Test points per atom |
| `n_slices` | 20 | LR | Slices per atom diameter |
| `probe_radius` | 1.4 Å | Both | Water molecule radius |

## Validation

Tested against FreeSASA reference implementation using ProtOr classifier (default for FreeSASA Python):

| Structure | Atoms | FreeSASA (Å²) | Zig (Å²) | Difference |
|-----------|------:|-------------:|--------:|----------:|
| 1CRN | 327 | 2,407.47 | 2,407.47 | 0.000% |
| 1UBQ | 602 | 4,803.47 | 4,803.47 | 0.000% |
| 1A0Q | 3,183 | 18,923.28 | 18,923.28 | 0.000% |
| 3HHB | 4,384 | 26,047.03 | 26,047.03 | 0.000% |
| 1AON | 58,674 | 291,728.28 | 291,728.28 | 0.000% |
| 4V6X | 237,685 | 1,141,803.28 | 1,141,803.28 | 0.000% |

Run validation: `./scripts/validate_accuracy.py`

## Performance

Benchmark comparing Zig (ReleaseFast) vs FreeSASA Python, SASA calculation time only:

| Structure | Atoms | SR Zig (ms) | SR FS (ms) | SR Speedup | LR Zig (ms) | LR FS (ms) | LR Speedup |
|-----------|------:|-----------:|-----------:|----------:|-----------:|-----------:|----------:|
| 1CRN | 327 | 0.21 | 0.40 | 1.90x | 0.44 | 0.46 | 1.05x |
| 1UBQ | 602 | 0.28 | 0.74 | 2.69x | 0.87 | 0.90 | 1.04x |
| 1A0Q | 3,183 | 1.28 | 4.71 | 3.67x | 4.86 | 5.37 | 1.10x |
| 3HHB | 4,384 | 2.02 | 7.20 | 3.57x | 4.64 | 7.72 | 1.66x |
| 1AON | 58,674 | 32.82 | 140.50 | 4.28x | 92.11 | 102.23 | 1.11x |
| 4V6X | 237,685 | 161.93 | 707.48 | 4.37x | 495.74 | 554.38 | 1.12x |

**Summary**: Shrake-Rupley is **1.9x-4.4x faster** than FreeSASA. Speedup increases with structure size.

Run benchmark: `./scripts/benchmark_all.py`

### Optimization Techniques

1. **Neighbor List**: O(N) neighbor lookup instead of O(N²)
2. **SIMD**: Process 4 calculations in parallel using `@Vector(4, f64)`
3. **Multi-threading**: Parallel atom processing with work-stealing thread pool

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
│   ├── cif_to_input_json.py       # Structure to JSON converter
│   ├── calc_reference_sasa.py     # Reference SASA calculator
│   ├── benchmark.py               # Single-structure benchmark
│   ├── benchmark_all.py           # Unified benchmark across all structures
│   ├── validate_accuracy.py       # Validate against FreeSASA references
│   ├── generate_benchmark_data.py # Download structures and generate references
│   └── compare_classifiers.py     # Compare classifier radii
├── benchmarks/
│   ├── structures/            # Downloaded PDB structures (.cif.gz)
│   ├── inputs/                # Generated input JSONs
│   └── references/            # FreeSASA reference SASA values
├── examples/
│   ├── 1A0Q.cif.gz        # Original structure file (PDB 1A0Q)
│   ├── input_1a0q.json    # Example input (converted from cif)
│   └── 1A0Q_sasa.json     # Reference SASA from FreeSASA
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
- [x] Phase 11: Lee-Richards algorithm (with multi-threading & SIMD)
- [ ] Phase 10: Direct mmCIF input support

## License

MIT

## References

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - Original C implementation

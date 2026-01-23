# freesasa-zig

[日本語](README.ja.md) | English

Solvent Accessible Surface Area (SASA) calculator implemented in Zig, ported from [FreeSASA](https://github.com/mittinatten/freesasa).

## Overview

SASA (Solvent Accessible Surface Area) measures the surface area of a biomolecule that is accessible to solvent molecules. This implementation provides two algorithms:

- **Shrake-Rupley (SR)**: Test point method with Golden Section Spiral
- **Lee-Richards (LR)**: Slice-based method with exact arc integration

## Features

- **Two SASA algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
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
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--algorithm=ALGO` | Algorithm: `sr` (Shrake-Rupley), `lr` (Lee-Richards) | sr |
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

**Extended format** (for future classifier support):

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["CA", "CB", "C"]
}
```

The `residue` and `atom_name` fields are optional. When provided, they enable automatic radius assignment via classifiers (planned for Phase 9).

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

Tested against FreeSASA reference implementation:

| Structure | Atoms | Reference (Å²) | Zig (Å²) | Difference |
|-----------|-------|----------------|----------|------------|
| 1A0Q | 3183 | 18923.28 | 19211.19 | 1.52% |

## Performance

Benchmark on PDB 1A0Q (3,183 atoms), ReleaseFast build:

| Algorithm | 1 thread | 4 threads | vs FreeSASA |
|-----------|----------|-----------|-------------|
| Shrake-Rupley | 17ms | 13ms | **3.9x faster** |
| Lee-Richards | 53ms | 26ms | 2.0x faster |
| FreeSASA (Python) | ~51ms | - | 1.0x |

### Optimization Techniques

1. **Neighbor List**: O(N) neighbor lookup instead of O(N²)
2. **SIMD**: Process 4 calculations in parallel using `@Vector(4, f64)`
3. **Multi-threading**: Parallel atom processing with work-stealing thread pool

## Project Structure

```
freesasa-zig/
├── src/
│   ├── main.zig           # CLI entry point
│   ├── root.zig           # Library root module
│   ├── types.zig          # Data structures (Vec3, AtomInput, etc.)
│   ├── json_parser.zig    # JSON input parsing and validation
│   ├── json_writer.zig    # Output writing (JSON, CSV)
│   ├── test_points.zig    # Golden Section Spiral generation
│   ├── neighbor_list.zig  # Spatial neighbor list (O(N) lookup)
│   ├── simd.zig           # SIMD batch operations
│   ├── thread_pool.zig    # Generic thread pool for parallelization
│   ├── shrake_rupley.zig  # Shrake-Rupley algorithm
│   └── lee_richards.zig   # Lee-Richards algorithm
├── scripts/
│   ├── cif_to_input_json.py   # Structure to JSON converter
│   ├── calc_reference_sasa.py # Reference SASA calculator
│   └── benchmark.py           # Performance benchmark
├── examples/
│   ├── 1A0Q.cif.gz        # Original structure file (PDB 1A0Q)
│   ├── input_1a0q.json    # Example input (converted from cif)
│   └── 1A0Q_sasa.json     # Reference SASA from FreeSASA
├── docs/                  # Technical documentation (Japanese)
│   ├── architecture.md    # Architecture overview
│   ├── algorithm.md       # Algorithm details
│   ├── optimizations.md   # Optimization techniques
│   └── cli-io.md          # CLI and I/O details
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
- [x] Phase 11: Lee-Richards algorithm (with multi-threading & SIMD)
- [ ] Phase 9: Radius classifier (auto-detect atom radii)
- [ ] Phase 10: Direct mmCIF input support

## License

MIT

## References

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - Original C implementation

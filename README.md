# freesasa-zig

Solvent Accessible Surface Area (SASA) calculator implemented in Zig, ported from [FreeSASA](https://github.com/mittinatten/freesasa).

## Overview

SASA (Solvent Accessible Surface Area) measures the surface area of a biomolecule that is accessible to solvent molecules. This implementation uses the Shrake-Rupley algorithm with Golden Section Spiral test points.

## Features

- Shrake-Rupley algorithm implementation
- JSON input/output format
- Configurable test points and probe radius
- Memory-safe implementation with explicit allocators
- **Multi-threading support** for parallel atom processing
- **SIMD optimization** using `@Vector(4, f64)` for batch distance calculations
- **Neighbor list optimization** for O(N) instead of O(N²) neighbor checking

## Building

Requires Zig 0.14.0 or later.

```bash
zig build
```

Run tests:

```bash
zig build test
```

## Usage

```bash
# Basic usage (auto-detect thread count)
./zig-out/bin/freesasa_zig <input.json> [output.json]

# Specify thread count
./zig-out/bin/freesasa_zig --threads=4 <input.json> [output.json]

# Single-threaded mode
./zig-out/bin/freesasa_zig --threads=1 <input.json> [output.json]

# Example
./zig-out/bin/freesasa_zig examples/input_1a0q.json output.json
```

### Options

| Option | Description |
|--------|-------------|
| `--threads=N` | Number of threads (default: auto-detect CPU count) |

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

### Output Format

```json
{
  "total_area": 1234.56,
  "atom_areas": [10.5, 20.3, ...]
}
```

## Preparing Input Data

Use the provided Python scripts to convert structure files:

```bash
# Convert mmCIF/PDB to input JSON
./scripts/cif_to_input_json.py structure.cif output.json

# Generate reference SASA using FreeSASA (for validation)
./scripts/calc_reference_sasa.py structure.cif reference.json
```

Requirements: Python 3.11+, gemmi, freesasa (installed automatically via PEP 723)

## Algorithm

### Shrake-Rupley Method

1. Generate uniformly distributed test points on a unit sphere using Golden Section Spiral
2. For each atom:
   - Scale test points by (van der Waals radius + probe radius)
   - Translate to atom position
   - Count points not buried inside neighboring atoms
   - SASA = 4π × r² × (exposed points / total points)

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `n_points` | 100 | Test points per atom |
| `probe_radius` | 1.4 Å | Water molecule radius |

## Validation

Tested against FreeSASA reference implementation:

| Structure | Atoms | Reference (Å²) | Zig (Å²) | Difference |
|-----------|-------|----------------|----------|------------|
| 1A0Q | 3183 | 18923.28 | 19211.19 | 1.52% |

## Performance

Benchmark on PDB 1A0Q (3,183 atoms):

| Configuration | Time | Speedup vs FreeSASA |
|---------------|------|---------------------|
| FreeSASA (Python) | ~51ms | 1.0x |
| Zig (single-threaded) | ~13ms | 3.9x |
| Zig (multi-threaded) | ~8ms | **6.4x** |

### Optimization Techniques

1. **Neighbor List**: O(N) neighbor lookup instead of O(N²)
2. **SIMD**: Process 4 distance calculations in parallel using `@Vector(4, f64)`
3. **Multi-threading**: Parallel atom processing with work-stealing thread pool

## Project Structure

```
freesasa-zig/
├── src/
│   ├── main.zig           # CLI entry point
│   ├── types.zig          # Data structures (Vec3, AtomInput, etc.)
│   ├── json_parser.zig    # JSON input parsing
│   ├── json_writer.zig    # JSON output writing
│   ├── test_points.zig    # Golden Section Spiral generation
│   ├── neighbor_list.zig  # Spatial neighbor list (O(N) lookup)
│   ├── simd.zig           # SIMD batch operations
│   ├── thread_pool.zig    # Generic thread pool for parallelization
│   └── shrake_rupley.zig  # Core SASA algorithm
├── scripts/
│   ├── cif_to_input_json.py   # Structure to JSON converter
│   ├── calc_reference_sasa.py # Reference SASA calculator
│   └── benchmark.py           # Performance benchmark
├── examples/
│   └── input_1a0q.json    # Example input (PDB 1A0Q)
└── plans/
    ├── phase-1-shrake-rupley.md
    ├── phase-3-simd-optimization.md
    └── phase-4-multithreading.md
```

## Roadmap

- [x] Phase 1: Basic Shrake-Rupley implementation
- [x] Phase 2: Neighbor list optimization (O(N²) → O(N))
- [x] Phase 3: SIMD optimization (4.5x faster than FreeSASA)
- [x] Phase 4: Multi-threading (6.4x faster than FreeSASA)
- [ ] Phase 5: Production features (more CLI options, configuration file)

## License

MIT

## References

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- [FreeSASA](https://github.com/mittinatten/freesasa) - Original C implementation

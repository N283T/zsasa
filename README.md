# freesasa-zig

[日本語](README.ja.md) | English

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.

## Features

- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Multiple input formats**: mmCIF, PDB, JSON
- **Analysis features**: Per-residue aggregation, RSA, polar/nonpolar classification
- **High performance**: SIMD optimization, multi-threading, neighbor list O(N)
- **Python bindings**: NumPy integration with BioPython/Biotite/Gemmi support

## Quick Start

**Requirements**: Zig 0.15.2+ ([download](https://ziglang.org/download/))

```bash
# Build
zig build -Doptimize=ReleaseFast

# Run
./zig-out/bin/freesasa_zig structure.cif output.json
```

## Installation

### CLI

```bash
git clone https://github.com/N283T/freesasa-zig.git
cd freesasa-zig
zig build -Doptimize=ReleaseFast
```

### Python

```bash
cd python
pip install -e .
```

Requires Zig 0.15.2+ for building the native library.

## Usage

### CLI Examples

```bash
# Basic SASA calculation
freesasa_zig structure.cif output.json

# Lee-Richards algorithm
freesasa_zig --algorithm=lr structure.cif output.json

# Multi-threaded
freesasa_zig --threads=4 structure.cif output.json

# Per-residue analysis with RSA
freesasa_zig --rsa structure.cif output.json

# CSV output
freesasa_zig --format=csv structure.cif output.csv
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) or `lr` (Lee-Richards) | sr |
| `--threads=N` | Number of threads (0 = auto) | auto |
| `--classifier=TYPE` | Atom classifier: `naccess`, `protor`, `oons` | - |
| `--chain=ID` | Filter by chain ID (e.g., `A` or `A,B,C`) | all |
| `--model=N` | Model number for NMR structures | all |
| `--per-residue` | Output per-residue SASA aggregation | - |
| `--rsa` | Calculate Relative Solvent Accessibility | - |
| `--polar` | Show polar/nonpolar summary | - |
| `--format=FMT` | Output: `json`, `compact`, `csv` | json |
| `--probe-radius=R` | Probe radius in Å | 1.4 |
| `--n-points=N` | Test points per atom (SR) | 100 |
| `--n-slices=N` | Slices per atom (LR) | 20 |
| `--precision=P` | `f32` (fast) or `f64` (precise) | f64 |

See [CLI Reference](docs/cli-io.md) for full options and input/output formats.

### Python Examples

```python
import numpy as np
from freesasa_zig import calculate_sasa

# Calculate SASA
coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])
result = calculate_sasa(coords, radii)
print(f"Total: {result.total_area:.2f} Å²")
```

With structure files (requires `pip install freesasa-zig[gemmi]`):

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
print(f"Total: {result.total_area:.1f} Å²")
```

See [Python README](python/README.md) for full documentation.

## Documentation

| Document | Description |
|----------|-------------|
| [CLI & I/O](docs/cli-io.md) | Command-line options, input/output formats |
| [Python Bindings](python/README.md) | Python API and integrations |
| [Algorithms](docs/algorithm.md) | Shrake-Rupley and Lee-Richards details |
| [Classifiers](docs/classifier.md) | NACCESS, ProtOr, OONS atom classifiers |
| [Optimizations](docs/optimizations.md) | SIMD, threading, performance techniques |
| [Benchmarks](docs/benchmark/) | Methodology and results |

## Performance

Validated against FreeSASA C on ~100k structures (mean error: 0.0004%).

**Speedup** (vs FreeSASA C, single-threaded):
- Shrake-Rupley: 1.2x - 2.3x faster
- Lee-Richards: ~1.7x faster

Multi-threading provides significant gains for structures with 500+ atoms.

See [benchmark results](docs/benchmark/results.md) for details.

## Project Structure

```
freesasa-zig/
├── src/           # Zig source (algorithms, parsers, CLI)
├── python/        # Python bindings
├── docs/          # Documentation
└── benchmarks/    # Benchmark tools and results
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.

## License

MIT

## References

- Shrake & Rupley (1973). Environment and exposure to solvent of protein atoms. *J. Mol. Biol.* 79(2), 351-371.
- Lee & Richards (1971). The interpretation of protein structures. *J. Mol. Biol.* 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - Original C implementation

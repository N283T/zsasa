# freesasa-zig

[日本語](README.ja.md) | English

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.

## Features

- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Multiple input formats**: mmCIF, PDB, JSON
- **Analysis features**: Per-residue aggregation, RSA, polar/nonpolar classification
- **High performance**: SIMD optimization, multi-threading, neighbor list O(N)
- **Python bindings**: NumPy integration with BioPython/Biotite/Gemmi support

## Benchmark Highlights

**Up to 3x faster** than FreeSASA C while maintaining **f64 precision** (mean error: <0.001%).

| Speedup (threads=10) | Thread Scaling (100k+ atoms) |
|:--------------------:|:----------------------------:|
| ![Speedup](benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](benchmarks/results/plots/large/speedup_by_threads.png) |

**Key Results (100k+ atoms, threads=10):**
- **2.3x** median speedup vs FreeSASA and RustSASA
- Speedup increases with thread count (superior parallel efficiency)

See [benchmark results](docs/benchmark/results.md) for detailed analysis.

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
./zig-out/bin/freesasa_zig structure.cif output.json

# Lee-Richards algorithm
./zig-out/bin/freesasa_zig --algorithm=lr structure.cif output.json

# Multi-threaded
./zig-out/bin/freesasa_zig --threads=4 structure.cif output.json

# Per-residue analysis with RSA
./zig-out/bin/freesasa_zig --rsa structure.cif output.json

# CSV output
./zig-out/bin/freesasa_zig --format=csv structure.cif output.csv
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

Run `./zig-out/bin/freesasa_zig --help` for all options. See [CLI Reference](docs/cli.md) for detailed documentation.

### Python Examples

```python
import numpy as np
from freesasa_zig import calculate_sasa

# Calculate SASA
coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])
result = calculate_sasa(coords, radii)
print(f"Total: {result.total_area:.2f} Å²")
print(f"Per-atom: {result.atom_areas}")
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
| [CLI Reference](docs/cli.md) | Command-line options, input/output formats |
| [Python Bindings](python/README.md) | Python API and integrations |
| [Algorithms](docs/algorithm.md) | Shrake-Rupley and Lee-Richards details |
| [Classifiers](docs/classifier.md) | NACCESS, ProtOr, OONS atom classifiers |
| [Optimizations](docs/optimizations.md) | SIMD, threading, performance techniques |
| [Benchmarks](docs/benchmark/) | Methodology and results |

## Performance

See [Benchmark Highlights](#benchmark-highlights) above and [detailed results](docs/benchmark/results.md).

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

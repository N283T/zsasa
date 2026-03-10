# zsasa

[![CI](https://github.com/N283T/zsasa/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/zsasa/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/zsasa?color=blue)](https://pypi.org/project/zsasa/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Zig](https://img.shields.io/badge/Zig-0.15.2+-f7a41d?logo=zig&logoColor=white)](https://ziglang.org/)
[![Python](https://img.shields.io/badge/Python-3.11+-3776ab?logo=python&logoColor=white)](https://www.python.org/)

[日本語](README.ja.md) | English

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.
**Up to 3x faster** than FreeSASA C with f64 precision.

**[Documentation](https://n283t.github.io/zsasa/)**

## Features

- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise), with bitmask LUT optimization
- **Multiple input formats**: mmCIF, PDB, JSON
- **MD trajectory analysis**: Native XTC and DCD readers, MDTraj and MDAnalysis integration
- **Batch processing**: Native directory processing for proteome-scale datasets
- **Analysis features**: Per-residue aggregation, RSA, polar/nonpolar classification
- **High performance**: SIMD optimization, multi-threading, f64/f32 selectable precision
- **Cross-platform**: Linux, macOS, and Windows (pre-built wheels via `pip install zsasa`)
- **Python bindings**: NumPy integration with Gemmi/BioPython/Biotite support

## Quick Start

### Python

```bash
pip install zsasa
# or
uv add zsasa
```

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
```

### CLI

Requires [Zig 0.15.2+](https://ziglang.org/download/).

```bash
git clone https://github.com/N283T/zsasa.git
cd zsasa && zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa calc structure.cif output.json
```

## Benchmarks

### Single-File (2,013 structures, n_points=100, threads=10)

| Metric | vs FreeSASA | vs RustSASA |
|--------|-------------|-------------|
| Median speedup | **1.88x** | **1.84x** |
| Thread scaling (t=1→10) | 2.71x | 1.39x |

### Batch (E. coli K-12 proteome, 4,370 structures)

| Tool | Time | RSS |
|------|------|-----|
| zsasa bitmask (f32) | **1.42s** | 43 MB |
| Lahuta bitmask | 2.01s | 291 MB |
| RustSASA | 5.24s | 169 MB |

### MD Trajectory

**82–96x less memory** than mdsasa-bolt (RustSASA). On 10K-frame datasets, mdsasa-bolt timed out (>2 hours) where zsasa completed in 38 seconds.

See [full benchmarks](https://n283t.github.io/zsasa/docs/benchmarks) for methodology and results.

## Documentation

Full documentation is available at **[n283t.github.io/zsasa](https://n283t.github.io/zsasa/)**.

| Section | Description |
|---------|-------------|
| [Getting Started](https://n283t.github.io/zsasa/docs/getting-started) | Installation and first calculation |
| [Comparison](https://n283t.github.io/zsasa/docs/comparison) | How zsasa compares to FreeSASA, RustSASA, and Lahuta |
| [CLI Reference](https://n283t.github.io/zsasa/docs/cli/commands) | Full CLI options and output formats |
| [Python API](https://n283t.github.io/zsasa/docs/python-api) | Core API, integrations, MD trajectory |
| [Benchmarks](https://n283t.github.io/zsasa/docs/benchmarks) | Performance and accuracy results |

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.

## License

[MIT](LICENSE)

## References

- Shrake, A.; Rupley, J. A. Environment and Exposure to Solvent of Protein Atoms. *J. Mol. Biol.* 1973, 79(2), 351-371. [doi:10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
- Lee, B.; Richards, F. M. The Interpretation of Protein Structures: Estimation of Static Accessibility. *J. Mol. Biol.* 1971, 55(3), 379-400. [doi:10.1016/0022-2836(71)90324-x](https://doi.org/10.1016/0022-2836(71)90324-x)
- Mitternacht, S. FreeSASA: An Open Source C Library for Solvent Accessible Surface Area Calculations. *F1000Res.* 2016, 5, 189. [doi:10.12688/f1000research.7931.1](https://doi.org/10.12688/f1000research.7931.1)
- Campbell, M. J. RustSASA: A Rust Crate for Accelerated Solvent Accessible Surface Area Calculations. *J. Open Source Softw.* 2026, 11(117), 9537. [doi:10.21105/joss.09537](https://doi.org/10.21105/joss.09537)

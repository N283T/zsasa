# zsasa

[![CI](https://github.com/N283T/zsasa/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/zsasa/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/zsasa?color=blue)](https://pypi.org/project/zsasa/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Zig](https://img.shields.io/badge/Zig-0.15.2+-f7a41d?logo=zig&logoColor=white)](https://ziglang.org/)
[![Python](https://img.shields.io/badge/Python-3.11+-3776ab?logo=python&logoColor=white)](https://www.python.org/)
[![Nix](https://img.shields.io/badge/Nix-flake-5277C3?logo=nixos&logoColor=white)](https://nixos.org/)

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.
**Up to 3x faster** than FreeSASA C with f64 precision.

**[Documentation](https://n283t.github.io/zsasa/)** · **[Benchmarks](https://n283t.github.io/zsasa/docs/benchmarks)** · **[Comparison](https://n283t.github.io/zsasa/docs/comparison)**

## Features

- **Two algorithms**: Shrake-Rupley and Lee-Richards, with bitmask LUT optimization
- **Multiple input formats**: mmCIF, PDB, JSON, XTC, DCD
- **Batch & trajectory**: Proteome-scale directory processing, MD trajectory analysis
- **Python bindings**: NumPy, Gemmi, BioPython, Biotite, MDTraj, MDAnalysis
- **High performance**: SIMD, multi-threading, f64/f32 selectable, zero dependencies
- **Cross-platform**: Linux, macOS, Windows (pre-built wheels on PyPI)

## Quick Start

### Python

```bash
pip install zsasa
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

```bash
# One-line install (downloads pre-built binary)
curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | sh

# Or with custom install directory
curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | INSTALL_DIR=/usr/local/bin sh

# Or with Nix
nix run github:N283T/zsasa -- calc structure.cif output.json

# Or build from source (requires Zig 0.15.2)
git clone https://github.com/N283T/zsasa.git
cd zsasa && zig build -Doptimize=ReleaseFast
```

```bash
zsasa calc structure.cif output.json
```

## Documentation

| | |
|---|---|
| [Getting Started](https://n283t.github.io/zsasa/docs/getting-started) | Installation and first calculation |
| [CLI Reference](https://n283t.github.io/zsasa/docs/cli/commands) | Commands, input formats, output options |
| [Python API](https://n283t.github.io/zsasa/docs/python-api) | Core API, integrations, trajectory |
| [Benchmarks](https://n283t.github.io/zsasa/docs/benchmarks) | Performance and accuracy results |
| [Comparison](https://n283t.github.io/zsasa/docs/comparison) | vs FreeSASA, RustSASA, Lahuta |

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and guidelines.

## License

[MIT](LICENSE)

## Acknowledgments

zsasa builds on ideas and methods from several excellent SASA implementations:

- **[FreeSASA](https://github.com/mittinatten/freesasa)** by Simon Mitternacht — the foundational C library that established open-source SASA calculation. zsasa uses FreeSASA as the accuracy reference for validation.
- **[RustSASA](https://github.com/maxall41/RustSASA)** by Max Campbell — a modern Rust implementation that demonstrated SIMD-accelerated SASA calculation.
- **[Lahuta](https://github.com/bisejdiu/lahuta)** by Besian I. Sejdiu — a C++ toolkit that pioneered bitmask LUT optimization for SASA.

We also thank the developers of [MDTraj](https://github.com/mdtraj/mdtraj), [MDAnalysis](https://github.com/MDAnalysis/mdanalysis), [Gemmi](https://github.com/project-gemmi/gemmi), [BioPython](https://github.com/biopython/biopython), and [Biotite](https://github.com/biotite-tools/biotite) for their libraries that make structural biology analysis accessible.

## References

- Shrake, A.; Rupley, J. A. *J. Mol. Biol.* 1973, 79(2), 351-371. [doi:10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
- Lee, B.; Richards, F. M. *J. Mol. Biol.* 1971, 55(3), 379-400. [doi:10.1016/0022-2836(71)90324-x](https://doi.org/10.1016/0022-2836(71)90324-x)
- Mitternacht, S. *F1000Res.* 2016, 5, 189. [doi:10.12688/f1000research.7931.1](https://doi.org/10.12688/f1000research.7931.1)
- Campbell, M. J. *J. Open Source Softw.* 2026, 11(117), 9537. [doi:10.21105/joss.09537](https://doi.org/10.21105/joss.09537)

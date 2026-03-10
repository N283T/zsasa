---
sidebar_position: 1
---

# zsasa

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.
**Up to 3x faster** than FreeSASA C with f64 precision.

## Features

- **Two algorithms**: Shrake-Rupley (fast, recommended) and Lee-Richards (precise)
- **Multiple input formats**: mmCIF, PDB, JSON
- **Analysis features**: Per-residue aggregation, RSA, polar/nonpolar classification
- **High performance**: SIMD optimization, multi-threading
- **Cross-platform**: Linux, macOS, and Windows
- **Python bindings**: NumPy integration with Gemmi/BioPython/Biotite support
- **MD trajectory analysis**: Native XTC reader, MDTraj and MDAnalysis integration

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

```bash
git clone https://github.com/N283T/zsasa.git
cd zsasa && zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa calc structure.cif output.json
```

## Documentation

| Section | Description |
|---------|-------------|
| [Getting Started](getting-started.md) | Installation and first calculation |
| [Algorithms](guide/algorithms.mdx) | SR vs LR algorithm comparison |
| [Classifiers](guide/classifiers.mdx) | Atom radius assignment |
| [Trajectory Analysis](guide/trajectory.md) | MD trajectory SASA |
| [CLI Reference](cli/commands.md) | Full CLI options |
| [Python API](python-api/) | Python bindings documentation |

---
sidebar_position: 1
---

# zsasa

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.
**Up to 3x faster** than FreeSASA C with f64 precision.

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
# One-line install
curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | sh
```

```bash
zsasa calc structure.cif output.json
```

## Documentation

| Section | Description |
|---------|-------------|
| [Getting Started](getting-started.md) | Installation and first calculation |
| [Comparison](comparison.md) | How zsasa compares to FreeSASA, RustSASA, and Lahuta |
| [Algorithms](guide/algorithms.mdx) | SR vs LR algorithm comparison |
| [Classifiers](guide/classifiers.mdx) | Atom radius assignment |
| [Trajectory Analysis](guide/trajectory.md) | MD trajectory SASA |
| [CLI Reference](cli/commands.md) | Full CLI options |
| [Python API](python-api/) | Python bindings documentation |
| [Benchmarks](benchmarks/single-file.md) | Performance and accuracy benchmarks |

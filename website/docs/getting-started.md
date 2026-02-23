---
sidebar_position: 2
---

# Getting Started

## Installation

### Python (Recommended)

```bash
pip install zsasa
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, arm64), and Windows (x86_64).
Python 3.11-3.13 supported.

For structure file support (mmCIF/PDB), install with an integration:

```bash
pip install zsasa[gemmi]      # Gemmi (fast mmCIF/PDB parsing)
pip install zsasa[biopython]  # BioPython
pip install zsasa[biotite]    # Biotite
```

### CLI

Requires [Zig 0.15.2+](https://ziglang.org/download/).

```bash
git clone https://github.com/N283T/zsasa.git
cd zsasa
zig build -Doptimize=ReleaseFast
```

The binary is at `./zig-out/bin/zsasa`.

## First Calculation

### Python

```python
import numpy as np
from zsasa import calculate_sasa

# Two atoms with coordinates and van der Waals radii
coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
print(f"Per-atom: {result.atom_areas}")
```

With a structure file (mmCIF/PDB):

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
print(f"Total SASA: {result.total_area:.1f} Å²")

# Per-residue analysis
for res in result.residue_areas:
    print(f"  {res.chain}:{res.name}{res.number} = {res.area:.1f} Å²")
```

### CLI

```bash
# Calculate SASA from a structure file
./zig-out/bin/zsasa structure.cif output.json

# With atom classifier and per-residue output
./zig-out/bin/zsasa --classifier=naccess --rsa structure.cif output.json

# CSV output
./zig-out/bin/zsasa --format=csv structure.cif output.csv
```

## What's Next?

- [User Guide: Algorithms](guide/algorithms.md) - Choosing between Shrake-Rupley and Lee-Richards
- [User Guide: Classifiers](guide/classifiers.md) - Atom radius assignment and its impact on results
- [CLI Reference](cli.md) - Full CLI options
- [Python API](python-api/) - Complete Python API documentation

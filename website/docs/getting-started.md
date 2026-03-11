---
sidebar_position: 2
---

# Getting Started

## Installation

### Python (Recommended)

```bash
pip install zsasa
# or
uv add zsasa
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, arm64), and Windows (x86_64).
Python 3.11-3.13 supported.

For structure file support (mmCIF/PDB), install with an integration:

```bash
pip install zsasa[gemmi]      # Gemmi (fast mmCIF/PDB parsing)
pip install zsasa[biopython]  # BioPython
pip install zsasa[biotite]    # Biotite
# or with uv
uv add zsasa[gemmi]
```

### CLI

```bash
# One-line install (downloads pre-built binary)
curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | sh

# Or with custom install directory
curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | INSTALL_DIR=/usr/local/bin sh

# Or build from source (requires Zig 0.15.2)
git clone https://github.com/N283T/zsasa.git
cd zsasa
zig build -Doptimize=ReleaseFast
```

The binary is at `./zig-out/bin/zsasa` (source build) or `~/.local/bin/zsasa` (install.sh).

### Nix

```bash
# Run directly
nix run github:N283T/zsasa -- calc structure.cif output.json

# Install to profile
nix profile install github:N283T/zsasa
```

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
from zsasa import aggregate_from_result

for res in aggregate_from_result(result):
    print(f"  {res.chain_id}:{res.residue_name}{res.residue_id} = {res.total_area:.1f} Å²")
```

### CLI

```bash
# Calculate SASA from a structure file
./zig-out/bin/zsasa calc structure.cif output.json

# With atom classifier and per-residue output
./zig-out/bin/zsasa calc --classifier=naccess --rsa structure.cif output.json

# CSV output
./zig-out/bin/zsasa calc --format=csv structure.cif output.csv
```

## What's Next?

- [User Guide: Algorithms](guide/algorithms.mdx) - Choosing between Shrake-Rupley and Lee-Richards
- [User Guide: Classifiers](guide/classifiers.mdx) - Atom radius assignment and its impact on results
- [CLI Reference](cli/commands.md) - Full CLI options
- [Python API](python-api/) - Complete Python API documentation

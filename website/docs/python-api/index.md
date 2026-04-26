---
sidebar_position: 1
sidebar_label: Overview
---

# Python API Reference

Comprehensive documentation for the zsasa Python bindings.

> Looking for the auto-generated API reference? See the **[Python Autodoc (pdoc)](pathname:///zsasa/python-autodoc/zsasa.html)**.

## Overview

The Python bindings provide:

- **NumPy integration**: Pass coordinates and radii as NumPy arrays
- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Multi-threading**: Automatic parallelization across CPU cores
- **Atom classification**: NACCESS, PROTOR, and OONS classifiers
- **RSA calculation**: Relative Solvent Accessibility
- **Per-residue aggregation**: Aggregate atom SASA to residue level
- **Directory batch processing**: Process entire directories of structure files
- **Library integrations**: gemmi, BioPython, Biotite, MDTraj, MDAnalysis (see [Integrations](../integrations/))

## Contents

| Document | Description |
|----------|-------------|
| [Core API](core.md) | `calculate_sasa`, batch API, directory processing |
| [Classifier](classifier.md) | Atom classification, RSA calculation |
| [Analysis](analysis.md) | Per-residue aggregation, examples |
| [Native XTC Reader](xtc.md) | Standalone XTC reading, no dependencies |

For structure file parsing and MD trajectory integrations, see [Integrations](../integrations/).

## Installation

### From PyPI (Recommended)

```bash
pip install zsasa
# or
uv add zsasa
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, arm64), and Windows (x86_64).
Python 3.11-3.13 supported.

### From Source (Development)

Requires Zig 0.16.0+.

```bash
cd zsasa/python
pip install -e .
# or
uv pip install -e .
```

### Optional Dependencies

```bash
# For structure file support
pip install zsasa[gemmi]     # gemmi (fast mmCIF/PDB)
pip install zsasa[biopython] # BioPython
pip install zsasa[biotite]   # Biotite (also works with AtomWorks)

# For MD trajectory analysis
pip install mdtraj                   # MDTraj
pip install MDAnalysis               # MDAnalysis

# All integrations
pip install zsasa[all]

# With uv
uv add zsasa[gemmi]          # or any extra above
```

---

## Quick Start

### Basic SASA Calculation

```python
import numpy as np
from zsasa import calculate_sasa

# Define atom coordinates (N, 3) and radii (N,)
coords = np.array([
    [0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
])
radii = np.array([1.5, 1.5])

# Calculate SASA
result = calculate_sasa(coords, radii)
print(f"Total: {result.total_area:.2f} Å²")
print(f"Per-atom: {result.atom_areas}")
```

### From Structure Files (Recommended)

```python
# With gemmi (pip install zsasa[gemmi])
from zsasa.integrations.gemmi import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
print(f"Total: {result.total_area:.1f} Å²")
print(f"Polar: {result.polar_area:.1f} Å²")
print(f"Apolar: {result.apolar_area:.1f} Å²")
```

### Per-Residue Analysis with RSA

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure
from zsasa import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

for res in residues:
    rsa_str = f"{res.rsa:.1%}" if res.rsa is not None else "N/A"
    print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: RSA={rsa_str}")
```

### Directory Batch Processing

```python
from zsasa import process_directory

result = process_directory("path/to/structures/")
print(f"Processed {result.successful}/{result.total_files} files")

for i in range(result.total_files):
    if result.status[i] == 1:
        print(f"  {result.filenames[i]}: SASA={result.total_sasa[i]:.1f} Å²")
```

### MD Trajectory Analysis

```python
import MDAnalysis as mda
from zsasa.mdanalysis import SASAAnalysis

u = mda.Universe("topology.pdb", "trajectory.xtc")
sasa = SASAAnalysis(u, select="protein")
sasa.run()

print(f"Mean SASA: {sasa.results.mean_total_area:.2f} Å²")
```

---

## Error Handling

### Common Exceptions

| Exception | Cause |
|-----------|-------|
| `ValueError` | Invalid input (shape mismatch, negative radii) |
| `FileNotFoundError` | Structure file not found |
| `IndexError` | Invalid model index |
| `ImportError` | Missing optional dependency (gemmi, biopython, biotite) |
| `MemoryError` | Out of memory during calculation |
| `RuntimeError` | Calculation error |

---

## Utility Functions

### get_version

```python
def get_version() -> str
```

Returns the library version string.

```python
from zsasa import get_version
print(get_version())  # e.g., "0.1.2"
```

---

## Generate Autodoc Locally

You can generate the full API reference locally using [pdoc](https://pdoc.dev/):

```bash
cd python
uv run pdoc zsasa --output-dir /tmp/python-autodoc
# Open /tmp/python-autodoc/zsasa.html in browser
```

Or serve with live reload:

```bash
cd python
uv run pdoc zsasa
# Opens http://localhost:8080 automatically
```

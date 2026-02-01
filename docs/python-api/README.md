# Python API Reference

Comprehensive documentation for the zsasa Python bindings.

## Overview

The Python bindings provide:

- **NumPy integration**: Pass coordinates and radii as NumPy arrays
- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Multi-threading**: Automatic parallelization across CPU cores
- **Atom classification**: NACCESS, PROTOR, and OONS classifiers
- **RSA calculation**: Relative Solvent Accessibility
- **Per-residue aggregation**: Aggregate atom SASA to residue level
- **Library integrations**: gemmi, BioPython, and Biotite support
- **MD trajectory analysis**: MDTraj and MDAnalysis integration (3.4x faster than mdsasa-bolt)

## Contents

| Document | Description |
|----------|-------------|
| [Core API](core.md) | `calculate_sasa`, batch API, precision |
| [Classifier](classifier.md) | Atom classification, RSA calculation |
| [Analysis](analysis.md) | Per-residue aggregation, examples |
| **Integrations** | |
| [Common Interface](integrations/README.md) | Shared API for all integrations |
| [Gemmi](integrations/gemmi.md) | Fast mmCIF/PDB parsing |
| [BioPython](integrations/biopython.md) | BioPython Structure support |
| [Biotite](integrations/biotite.md) | Biotite AtomArray support |
| [MDTraj](integrations/mdtraj.md) | MD trajectory analysis |
| [MDAnalysis](integrations/mdanalysis.md) | MD trajectory analysis |

## Installation

### Requirements

- Python 3.11+
- NumPy 1.24+
- Zig 0.15.2+ (for building the native library)

### From Source

```bash
# 1. Build the Zig shared library
cd zsasa
zig build -Doptimize=ReleaseFast

# 2. Install Python package
cd python
pip install -e .
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
```

### Library Location

The bindings look for the native library in these locations:

1. `ZSASA_LIB` environment variable (if set)
2. Bundled in package (wheel installation)
3. `../zig-out/lib/` (relative to package, for development)
4. `/usr/local/lib/`
5. `/usr/lib/`
6. Current working directory
7. `./zig-out/lib/` (current directory's zig-out)

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

## Performance

Comparison with FreeSASA Python bindings (single-threaded):

| Structure | Atoms | Zig SR | FS SR | Speedup |
|-----------|------:|-------:|------:|--------:|
| 1CRN | 327 | 0.5ms | 0.7ms | 1.4x |
| 1UBQ | 602 | 0.6ms | 1.3ms | 2.2x |
| 1A0Q | 3,183 | 2.4ms | 7.6ms | 3.2x |
| 3HHB | 4,384 | 3.4ms | 11ms | 3.2x |
| 1AON | 58,674 | 44ms | 163ms | 3.7x |

- **SR algorithm**: 1.4-3.7x faster (speedup increases with size)
- **LR algorithm**: 2.9-5.5x faster
- **Accuracy**: Results match FreeSASA (< 0.01% difference)

### MD Trajectory Analysis

Comparison with mdsasa-bolt (RustSASA) for MD trajectory SASA:

| Implementation | Time (20k atoms × 1k frames) | Notes |
|----------------|------------------------------|-------|
| zsasa   | 8.8 s                        | f64, controllable threads |
| mdsasa-bolt    | 30.3 s                       | f32, rayon global pool |
| **Speedup**    | **3.4x**                     | |

*Benchmark: MD ATLAS 6qfk_A trajectory (20,391 atoms, 1,001 frames, n_points=100, 8 threads)*

**Key advantages:**

- **Controllable parallelism**: Set exact thread count with `n_threads`, unlike rayon's global pool
- **Higher precision**: f64 by default (f32 available for comparison)
- **Efficient scaling**: 5.2x speedup from 1→8 threads

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

### Example

```python
from zsasa import calculate_sasa
import numpy as np

try:
    result = calculate_sasa(
        np.array([[0, 0, 0]]),
        np.array([-1.0])  # Invalid: negative radius
    )
except ValueError as e:
    print(f"Error: {e}")
```

---

## Utility Functions

### get_version

```python
def get_version() -> str
```

Returns the library version string (e.g., "0.1.0").

```python
from zsasa import get_version
print(get_version())  # "0.1.0"
```

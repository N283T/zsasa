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
- **MD trajectory analysis**: MDTraj and MDAnalysis integration (4x faster than mdsasa-bolt)

## Contents

| Document | Description |
|----------|-------------|
| [Core API](core.md) | `calculate_sasa`, batch API, precision |
| [Classifier](classifier.md) | Atom classification, RSA calculation |
| [Analysis](analysis.md) | Per-residue aggregation, examples |
| [Native XTC Reader](xtc.md) | Standalone XTC reading, no dependencies |
| **Integrations** | |
| [Common Interface](integrations/README.md) | Shared API for all integrations |
| [Gemmi](integrations/gemmi.md) | Fast mmCIF/PDB parsing |
| [BioPython](integrations/biopython.md) | BioPython Structure support |
| [Biotite](integrations/biotite.md) | Biotite AtomArray support |
| [MDTraj](integrations/mdtraj.md) | MD trajectory analysis |
| [MDAnalysis](integrations/mdanalysis.md) | MD trajectory analysis |

## Installation

### From PyPI (Recommended)

```bash
pip install zsasa
```

Pre-built wheels are available for Linux (x86_64, aarch64), macOS (x86_64, arm64), and Windows (x86_64).
Python 3.11-3.13 supported.

### From Source (Development)

Requires Zig 0.15.2+.

```bash
cd zsasa/python
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

### Native XTC Reading (No Dependencies)

```python
import numpy as np
from zsasa.xtc import XtcReader, compute_sasa_trajectory

# Low-level reader
with XtcReader("trajectory.xtc") as reader:
    for frame in reader:
        print(f"Step {frame.step}: {frame.natoms} atoms")

# High-level SASA calculation
radii = np.full(304, 1.7)  # Atomic radii in Angstroms
result = compute_sasa_trajectory("trajectory.xtc", radii)
print(f"Total SASA: {result.total_areas}")
```

---

## Performance

The Python bindings call the same high-performance Zig engine as the CLI. See [benchmark results](../benchmark/) for detailed methodology.

### Single-File Performance

| Speedup (threads=10) | Thread Scaling (100k+ atoms) |
|:--------------------:|:----------------------------:|
| ![Speedup](../../benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](../../benchmarks/results/plots/thread_scaling/individual/sr.png) |

**Key Results (100k+ atoms, threads=10):**
- **2.3x** median speedup vs FreeSASA and RustSASA
- Speedup increases with thread count (superior parallel efficiency)
- Accuracy: Results match FreeSASA (< 0.01% difference)

> **Note**: Zig/FreeSASA use f64, RustSASA uses f32.

See [single-file benchmark results](../benchmark/single-file.md) for detailed analysis.

### MD Trajectory Performance

**4.3x faster** than mdsasa-bolt (RustSASA) on real MD trajectory data.

![MD Trajectory Benchmark](../../benchmarks/results/md/6sup_A_analysis/plots/bar.png)

*33,377 atoms, 1,001 frames, n_points=100*

**Key advantages:**

- **Controllable parallelism**: Set exact thread count with `n_threads`, unlike rayon's global pool
- **Higher precision**: f64 by default (f32 available for comparison)
- **Thread scaling**: 6.0x speedup from 1→10 threads

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

Returns the library version string (e.g., "0.1.1").

```python
from zsasa import get_version
print(get_version())  # "0.1.1"
```

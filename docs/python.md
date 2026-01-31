# Python API Reference

This document provides comprehensive documentation for the freesasa-zig Python bindings.

## Overview

The Python bindings provide:

- **NumPy integration**: Pass coordinates and radii as NumPy arrays
- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Multi-threading**: Automatic parallelization across CPU cores
- **Atom classification**: NACCESS, PROTOR, and OONS classifiers
- **RSA calculation**: Relative Solvent Accessibility
- **Per-residue aggregation**: Aggregate atom SASA to residue level
- **Library integrations**: gemmi, BioPython, and Biotite support

## Installation

### Requirements

- Python 3.11+
- NumPy 1.24+
- Zig 0.15.2+ (for building the native library)

### From Source

```bash
# 1. Build the Zig shared library
cd freesasa-zig
zig build -Doptimize=ReleaseFast

# 2. Install Python package
cd python
pip install -e .
```

### Optional Dependencies

```bash
# For structure file support
pip install freesasa-zig[gemmi]     # gemmi (fast mmCIF/PDB)
pip install freesasa-zig[biopython] # BioPython
pip install freesasa-zig[biotite]   # Biotite (also works with AtomWorks)

# All integrations
pip install freesasa-zig[all]
```

### Library Location

The bindings look for the native library in these locations:

1. `FREESASA_ZIG_LIB` environment variable (if set)
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
from freesasa_zig import calculate_sasa

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
# With gemmi (pip install freesasa-zig[gemmi])
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
print(f"Total: {result.total_area:.1f} Å²")
print(f"Polar: {result.polar_area:.1f} Å²")
print(f"Apolar: {result.apolar_area:.1f} Å²")
```

### Per-Residue Analysis with RSA

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
from freesasa_zig import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

for res in residues:
    rsa_str = f"{res.rsa:.1%}" if res.rsa is not None else "N/A"
    print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: RSA={rsa_str}")
```

---

## Core API

### calculate_sasa

```python
def calculate_sasa(
    coords: NDArray[np.float64],
    radii: NDArray[np.float64],
    *,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
) -> SasaResult
```

Calculate Solvent Accessible Surface Area.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `coords` | `NDArray[float64]` | required | Atom coordinates as (N, 3) array |
| `radii` | `NDArray[float64]` | required | Atom radii as (N,) array |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Shrake-Rupley or Lee-Richards |
| `n_points` | `int` | `100` | Test points per atom (SR only) |
| `n_slices` | `int` | `20` | Slices per atom (LR only) |
| `probe_radius` | `float` | `1.4` | Water probe radius in Å |
| `n_threads` | `int` | `0` | Number of threads (0 = auto) |

**Returns:** `SasaResult`

**Raises:** `ValueError` for invalid input

### SasaResult

```python
@dataclass
class SasaResult:
    total_area: float           # Total SASA in Å²
    atom_areas: NDArray[float64] # Per-atom SASA values
```

---

## Classifier API

### ClassifierType

```python
class ClassifierType(IntEnum):
    NACCESS = 0  # NACCESS-compatible radii (default)
    PROTOR = 1   # ProtOr radii
    OONS = 2     # OONS radii
```

### AtomClass

```python
class AtomClass(IntEnum):
    POLAR = 0    # Polar atoms (N, O, etc.)
    APOLAR = 1   # Apolar atoms (C, etc.)
    UNKNOWN = 2  # Unknown classification
```

### classify_atoms

```python
def classify_atoms(
    residues: list[str],
    atoms: list[str],
    classifier_type: ClassifierType = ClassifierType.NACCESS,
    *,
    include_classes: bool = True,
) -> ClassificationResult
```

Classify multiple atoms at once (batch operation).

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `residues` | `list[str]` | required | Residue names (e.g., ["ALA", "GLY"]) |
| `atoms` | `list[str]` | required | Atom names (e.g., ["CA", "N"]) |
| `classifier_type` | `ClassifierType` | `NACCESS` | Classifier to use |
| `include_classes` | `bool` | `True` | Whether to compute atom classes |

**Returns:** `ClassificationResult`

**Example:**

```python
from freesasa_zig import classify_atoms, ClassifierType

result = classify_atoms(
    ["ALA", "ALA", "GLY"],
    ["CA", "O", "N"],
    ClassifierType.NACCESS
)
print(result.radii)   # [1.87, 1.4, 1.65]
print(result.classes) # [1, 0, 0] (APOLAR, POLAR, POLAR)
```

### ClassificationResult

```python
@dataclass
class ClassificationResult:
    radii: NDArray[float64]  # Per-atom radii (NaN for unknown)
    classes: NDArray[int32]  # Per-atom polarity classes
```

### Utility Functions

| Function | Description |
|----------|-------------|
| `get_radius(residue, atom, classifier)` | Get radius for a specific atom |
| `get_atom_class(residue, atom, classifier)` | Get polarity class for an atom |
| `guess_radius(element)` | Guess radius from element symbol |
| `guess_radius_from_atom_name(atom_name)` | Guess radius from PDB atom name |

---

## RSA API

### MAX_SASA

Reference maximum SASA values from Tien et al. (2013):

```python
MAX_SASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 225.0, "GLU": 223.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
```

### calculate_rsa

```python
def calculate_rsa(sasa: float, residue_name: str) -> float | None
```

Calculate Relative Solvent Accessibility (RSA = SASA / MaxSASA).

**Returns:** RSA value (0.0-1.0+), or `None` for non-standard amino acids.

**Example:**

```python
from freesasa_zig import calculate_rsa

rsa = calculate_rsa(64.5, "ALA")  # 64.5 / 129.0 = 0.5
```

### calculate_rsa_batch

```python
def calculate_rsa_batch(
    sasas: NDArray[float64] | list[float],
    residue_names: list[str],
) -> NDArray[float64]
```

Batch RSA calculation. Returns NaN for non-standard amino acids.

### get_max_sasa

```python
def get_max_sasa(residue_name: str) -> float | None
```

Get maximum SASA reference value for a residue.

---

## Analysis API

### aggregate_by_residue

```python
def aggregate_by_residue(
    atom_areas: NDArray[float64],
    chain_ids: list[str],
    residue_ids: list[int],
    residue_names: list[str],
    atom_classes: NDArray[int32] | None = None,
) -> list[ResidueResult]
```

Aggregate per-atom SASA values to per-residue.

### aggregate_from_result

```python
def aggregate_from_result(result: SasaResultWithAtoms) -> list[ResidueResult]
```

Convenience wrapper for `SasaResultWithAtoms` from integration modules.

### ResidueResult

```python
@dataclass
class ResidueResult:
    chain_id: str       # Chain identifier
    residue_id: int     # Residue sequence number
    residue_name: str   # 3-letter residue code
    total_area: float   # Total SASA in Å²
    polar_area: float   # Polar SASA in Å²
    apolar_area: float  # Apolar SASA in Å²
    rsa: float | None   # Relative Solvent Accessibility
    n_atoms: int        # Number of atoms
```

---

## Integration Modules

All integration modules provide the same interface:

### Common Functions

| Function | Description |
|----------|-------------|
| `extract_atoms_from_model(model)` | Extract atom data (gemmi, BioPython) |
| `extract_atoms_from_atom_array(atom_array)` | Extract atom data (Biotite) |
| `calculate_sasa_from_model(model)` | Calculate SASA from model (gemmi, BioPython) |
| `calculate_sasa_from_atom_array(atom_array)` | Calculate SASA from AtomArray (Biotite) |
| `calculate_sasa_from_structure(source)` | Calculate SASA from file or structure (all) |

### Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `classifier` | `ClassifierType` | `NACCESS` | Atom radius classifier |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | SASA algorithm |
| `n_points` | `int` | `100` | Test points (SR) |
| `n_slices` | `int` | `20` | Slices (LR) |
| `probe_radius` | `float` | `1.4` | Probe radius in Å |
| `n_threads` | `int` | `0` | Threads (0 = auto) |
| `include_hetatm` | `bool` | `True` | Include HETATM records |
| `include_hydrogens` | `bool` | `False` | Include hydrogen atoms |
| `model_index` | `int` | `0` | Model index (NMR) |

### SasaResultWithAtoms

```python
@dataclass
class SasaResultWithAtoms(SasaResult):
    atom_classes: NDArray[int32]  # Per-atom polarity classes
    atom_data: AtomData           # Atom metadata
    polar_area: float             # Total polar SASA
    apolar_area: float            # Total apolar SASA
```

### AtomData

```python
@dataclass
class AtomData:
    coords: NDArray[float64]   # (N, 3) coordinates
    residue_names: list[str]   # Residue names
    atom_names: list[str]      # Atom names
    chain_ids: list[str]       # Chain IDs
    residue_ids: list[int]     # Residue numbers
    elements: list[str]        # Element symbols
```

---

## Gemmi Integration

```python
from freesasa_zig.integrations.gemmi import (
    calculate_sasa_from_structure,
    calculate_sasa_from_model,
    extract_atoms_from_model,
)
```

**Requires:** `pip install freesasa-zig[gemmi]`

**Supported formats:** mmCIF, PDB

**Example:**

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure

# From file
result = calculate_sasa_from_structure("protein.cif")

# From gemmi Structure
import gemmi
structure = gemmi.read_structure("protein.pdb")
result = calculate_sasa_from_structure(structure)
```

---

## BioPython Integration

```python
from freesasa_zig.integrations.biopython import (
    calculate_sasa_from_structure,
    calculate_sasa_from_model,
    extract_atoms_from_model,
)
```

**Requires:** `pip install freesasa-zig[biopython]`

**Supported formats:** PDB, mmCIF

**Example:**

```python
from freesasa_zig.integrations.biopython import calculate_sasa_from_structure

# From file
result = calculate_sasa_from_structure("protein.pdb")

# From BioPython Structure
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "protein.pdb")
result = calculate_sasa_from_structure(structure)
```

---

## Biotite Integration

```python
from freesasa_zig.integrations.biotite import (
    calculate_sasa_from_structure,
    calculate_sasa_from_atom_array,
    extract_atoms_from_atom_array,
)
```

**Requires:** `pip install freesasa-zig[biotite]`

**Supported formats:** PDB, mmCIF, BinaryCIF

**Also works with AtomWorks** (built on Biotite):

```python
from atomworks.io.utils.io_utils import load_any
from freesasa_zig.integrations.biotite import calculate_sasa_from_atom_array

atom_array = load_any("protein.cif.gz")
result = calculate_sasa_from_atom_array(atom_array)
```

---

## Examples

### Algorithm Comparison

```python
import numpy as np
from freesasa_zig import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

# Shrake-Rupley (faster)
result_sr = calculate_sasa(coords, radii, algorithm="sr", n_points=100)
print(f"SR: {result_sr.total_area:.2f} Å²")

# Lee-Richards (more precise)
result_lr = calculate_sasa(coords, radii, algorithm="lr", n_slices=20)
print(f"LR: {result_lr.total_area:.2f} Å²")
```

### Multi-threading

```python
from freesasa_zig import calculate_sasa

# Auto-detect CPU cores (default)
result = calculate_sasa(coords, radii, n_threads=0)

# Specific thread count
result = calculate_sasa(coords, radii, n_threads=4)

# Single-threaded
result = calculate_sasa(coords, radii, n_threads=1)
```

### Finding Buried Residues

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
from freesasa_zig import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

# Find buried residues (RSA < 20%)
buried = [r for r in residues if r.rsa is not None and r.rsa < 0.2]

print(f"Buried residues ({len(buried)}):")
for r in buried:
    print(f"  {r.chain_id}:{r.residue_name}{r.residue_id} - RSA: {r.rsa:.1%}")
```

### Custom Radii

```python
import numpy as np
from freesasa_zig import calculate_sasa, classify_atoms

# Get radii from classifier
residues = ["ALA", "ALA", "ALA"]
atoms = ["N", "CA", "C"]
classification = classify_atoms(residues, atoms)

# Or use custom radii
custom_radii = np.array([1.65, 1.87, 1.76])

# Calculate with either
coords = np.array([[0, 0, 0], [1.5, 0, 0], [2.3, 1.2, 0]])
result = calculate_sasa(coords, classification.radii)
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

### Example

```python
from freesasa_zig import calculate_sasa
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

---

## Utility Functions

### get_version

```python
def get_version() -> str
```

Returns the library version string (e.g., "0.1.0").

```python
from freesasa_zig import get_version
print(get_version())  # "0.1.0"
```

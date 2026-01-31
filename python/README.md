# freesasa-zig Python Bindings

Python bindings for the freesasa-zig library - a high-performance implementation of Solvent Accessible Surface Area (SASA) calculation.

## Features

- **NumPy integration**: Pass coordinates and radii as NumPy arrays
- **Two algorithms**: Shrake-Rupley (point-based) and Lee-Richards (slice-based)
- **Multi-threading**: Automatic parallelization across CPU cores
- **Atom classification**: NACCESS, PROTOR, and OONS classifiers with polar/apolar assignment
- **RSA calculation**: Relative Solvent Accessibility using standard reference values
- **Per-residue aggregation**: Aggregate atom-level SASA to residue-level with RSA
- **Library integrations**: Built-in support for gemmi, BioPython, and Biotite
- **Type-safe**: Full type hints and runtime validation

## Installation

### From Source (Development)

```bash
# 1. Build the Zig shared library
cd /path/to/freesasa-zig
zig build -Doptimize=ReleaseFast

# 2. Install Python package
cd python
pip install -e .
```

### Optional Dependencies

```bash
# For gemmi integration
pip install freesasa-zig[gemmi]

# For BioPython integration
pip install freesasa-zig[biopython]

# For Biotite integration (also works with AtomWorks)
pip install freesasa-zig[biotite]

# All integrations
pip install freesasa-zig[all]
```

### Library Location

The Python bindings look for `libfreesasa_zig.dylib` (macOS), `libfreesasa_zig.so` (Linux), or `freesasa_zig.dll` (Windows) in these locations:

1. `FREESASA_ZIG_LIB` environment variable (if set)
2. `../zig-out/lib/` (relative to package)
3. `/usr/local/lib/`
4. `/usr/lib/`
5. Current working directory

## Quick Start

### Basic SASA Calculation

```python
import numpy as np
from freesasa_zig import calculate_sasa, get_version

print(f"Library version: {get_version()}")

# Define atom coordinates (N, 3) and radii (N,)
coords = np.array([
    [0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [6.0, 0.0, 0.0],
])
radii = np.array([1.5, 1.5, 1.5])

# Calculate SASA
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
print(f"Per-atom areas: {result.atom_areas}")
```

### Using Library Integrations (Recommended)

The easiest way to calculate SASA from structure files:

```python
# With gemmi (fast mmCIF/PDB parser)
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
result = calculate_sasa_from_structure("protein.cif")

# With BioPython
from freesasa_zig.integrations.biopython import calculate_sasa_from_structure
result = calculate_sasa_from_structure("protein.pdb")

# With Biotite (also works with AtomWorks)
from freesasa_zig.integrations.biotite import calculate_sasa_from_structure
result = calculate_sasa_from_structure("protein.pdb")

print(f"Total: {result.total_area:.1f} Å²")
print(f"Polar: {result.polar_area:.1f} Å²")
print(f"Apolar: {result.apolar_area:.1f} Å²")
```

### Per-Residue Analysis with RSA

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
from freesasa_zig.analysis import aggregate_from_result

# Calculate SASA
result = calculate_sasa_from_structure("protein.cif")

# Aggregate to per-residue
residues = aggregate_from_result(result)

for res in residues:
    rsa_str = f"{res.rsa:.1%}" if res.rsa is not None else "N/A"
    print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: "
          f"{res.total_area:.1f} Å² (RSA: {rsa_str})")
```

## API Reference

### Core Functions

#### `calculate_sasa(coords, radii, *, algorithm="sr", n_points=100, n_slices=20, probe_radius=1.4, n_threads=0)`

Calculate Solvent Accessible Surface Area.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `coords` | `NDArray[float64]` | required | Atom coordinates as (N, 3) array |
| `radii` | `NDArray[float64]` | required | Atom radii as (N,) array |
| `algorithm` | `"sr" \| "lr"` | `"sr"` | Shrake-Rupley or Lee-Richards |
| `n_points` | `int` | `100` | Test points per atom (SR only) |
| `n_slices` | `int` | `20` | Slices per atom (LR only) |
| `probe_radius` | `float` | `1.4` | Water probe radius in Å |
| `n_threads` | `int` | `0` | Number of threads (0 = auto-detect) |

**Returns:** `SasaResult`

#### `SasaResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_area` | `float` | Total SASA in Å² |
| `atom_areas` | `NDArray[float64]` | Per-atom SASA values |

#### `classify_atoms(residue_names, atom_names, classifier=NACCESS)`

Classify atoms and get radii.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `residue_names` | `list[str]` | required | Residue names (e.g., ["ALA", "GLY"]) |
| `atom_names` | `list[str]` | required | Atom names (e.g., ["CA", "N"]) |
| `classifier` | `ClassifierType` | `NACCESS` | NACCESS, PROTOR, or OONS |
| `include_classes` | `bool` | `True` | Whether to include atom classes |

**Returns:** `ClassificationResult` with `radii` and `classes` arrays

#### `get_version()`

Returns the library version string (e.g., "0.1.0").

#### `calculate_rsa(sasa, residue_name)`

Calculate Relative Solvent Accessibility.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `sasa` | `float` | SASA value in Å² |
| `residue_name` | `str` | 3-letter residue code |

**Returns:** `float | None` (None for non-standard amino acids)

### Utility Functions

These functions are available for advanced use cases:

| Function | Description |
|----------|-------------|
| `get_radius(residue, atom, classifier)` | Get radius for a specific atom |
| `get_atom_class(residue, atom, classifier)` | Get polarity class for an atom |
| `guess_radius(element)` | Guess radius from element symbol |
| `guess_radius_from_atom_name(atom_name)` | Guess radius from PDB atom name |
| `get_max_sasa(residue_name)` | Get max SASA reference value for RSA |
| `calculate_rsa_batch(sasa_values, residue_names)` | Batch RSA calculation |

### Analysis Functions

#### `aggregate_from_result(result)`

Aggregate per-atom SASA to per-residue.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `result` | `SasaResultWithAtoms` | Result from integration module |

**Returns:** `list[ResidueResult]`

#### `ResidueResult`

| Attribute | Type | Description |
|-----------|------|-------------|
| `chain_id` | `str` | Chain identifier |
| `residue_id` | `int` | Residue sequence number |
| `residue_name` | `str` | 3-letter residue code |
| `total_area` | `float` | Total SASA in Å² |
| `polar_area` | `float` | Polar SASA in Å² |
| `apolar_area` | `float` | Apolar SASA in Å² |
| `rsa` | `float \| None` | Relative Solvent Accessibility |
| `n_atoms` | `int` | Number of atoms |

### Integration Modules

All integration modules provide the same API:

- `extract_atoms_from_model(model)` - Extract atom data
- `calculate_sasa_from_model(model)` - Calculate SASA from model object
- `calculate_sasa_from_structure(source)` - Calculate SASA from file or structure

#### `SasaResultWithAtoms`

Extended result with atom metadata:

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_area` | `float` | Total SASA in Å² |
| `atom_areas` | `NDArray[float64]` | Per-atom SASA values |
| `atom_classes` | `NDArray[int32]` | Per-atom polarity classes |
| `atom_data` | `AtomData` | Atom metadata |
| `polar_area` | `float` | Total polar SASA in Å² |
| `apolar_area` | `float` | Total apolar SASA in Å² |

### Enums

#### `ClassifierType`

| Value | Description |
|-------|-------------|
| `NACCESS` | NACCESS classifier (default) |
| `PROTOR` | ProtOr classifier |
| `OONS` | OONS classifier |

#### `AtomClass`

| Value | Description |
|-------|-------------|
| `POLAR` | Polar atom (N, O, etc.) |
| `APOLAR` | Apolar atom (C, S, etc.) |
| `UNKNOWN` | Unknown classification |

## Examples

### Algorithm Comparison

```python
import numpy as np
from freesasa_zig import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

# Shrake-Rupley (default, faster)
result_sr = calculate_sasa(coords, radii, algorithm="sr", n_points=100)
print(f"SR: {result_sr.total_area:.2f} Å²")

# Lee-Richards (more precise)
result_lr = calculate_sasa(coords, radii, algorithm="lr", n_slices=20)
print(f"LR: {result_lr.total_area:.2f} Å²")
```

### Multi-threading

```python
# Auto-detect CPU cores (default)
result = calculate_sasa(coords, radii, n_threads=0)

# Use specific number of threads
result = calculate_sasa(coords, radii, n_threads=4)

# Single-threaded
result = calculate_sasa(coords, radii, n_threads=1)
```

### Atom Classification

```python
from freesasa_zig import classify_atoms, ClassifierType, AtomClass

residue_names = ["ALA", "ALA", "ALA"]
atom_names = ["N", "CA", "O"]

# With default classifier (NACCESS)
result = classify_atoms(residue_names, atom_names)

# Or with explicit classifier
result = classify_atoms(residue_names, atom_names, ClassifierType.PROTOR)

print(f"Radii: {result.radii}")
print(f"Classes: {result.classes}")

# Count polar/apolar
polar_count = sum(1 for c in result.classes if c == AtomClass.POLAR)
apolar_count = sum(1 for c in result.classes if c == AtomClass.APOLAR)
print(f"Polar: {polar_count}, Apolar: {apolar_count}")
```

### Finding Buried Residues

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
from freesasa_zig.analysis import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

# Find buried residues (RSA < 20%)
buried = [r for r in residues if r.rsa is not None and r.rsa < 0.2]

print(f"Buried residues ({len(buried)}):")
for r in buried:
    print(f"  {r.chain_id}:{r.residue_name}{r.residue_id} - RSA: {r.rsa:.1%}")
```

### Working with AtomWorks

```python
# AtomWorks is built on Biotite, so use the biotite integration
from atomworks.io.utils.io_utils import load_any
from freesasa_zig.integrations.biotite import calculate_sasa_from_atom_array
from freesasa_zig.analysis import aggregate_from_result

# Load with AtomWorks
atom_array = load_any("protein.cif.gz")

# Calculate SASA with freesasa-zig
result = calculate_sasa_from_atom_array(atom_array)
residues = aggregate_from_result(result)

print(f"Total SASA: {result.total_area:.1f} Å²")
```

## Performance

Library-to-library comparison (Python bindings vs FreeSASA Python):

| Structure | Atoms | Zig SR | FS SR | SR Speedup | Zig LR | FS LR | LR Speedup |
|-----------|------:|-------:|------:|----------:|-------:|------:|----------:|
| 1CRN | 327 | 0.5ms | 0.7ms | 1.4x | 1.5ms | 4.4ms | 2.9x |
| 1UBQ | 602 | 0.6ms | 1.3ms | 2.2x | 2.0ms | 8.4ms | 4.2x |
| 1A0Q | 3,183 | 2.4ms | 7.6ms | 3.2x | 8.9ms | 48ms | 5.5x |
| 3HHB | 4,384 | 3.4ms | 11ms | 3.2x | 12ms | 69ms | 5.5x |
| 1AON | 58,674 | 44ms | 163ms | 3.7x | 171ms | 931ms | 5.5x |

- **SR algorithm**: Zig is 1.4-3.7x faster (speedup increases with size)
- **LR algorithm**: Zig is 2.9-5.5x faster
- **Accuracy**: Results match FreeSASA (< 0.01% difference)

Run benchmark: `./benchmarks/scripts/run.py --tool zig --algorithm sr`

## Development

### Running Tests

```bash
cd python
uv run --with pytest pytest tests/ -v
```

### Type Checking

```bash
cd python
uv run --with ty ty check freesasa_zig/
```

### Linting

```bash
ruff format .
ruff check --fix .
```

## References

### Integration Libraries

- **Gemmi**: Wojdyr, M. GEMMI: A Library for Structural Biology. *J. Open Source Softw.* 2022, 7(73), 4200. [doi:10.21105/joss.04200](https://doi.org/10.21105/joss.04200)
- **BioPython**: Cock, P. J. A. et al. Biopython: Freely Available Python Tools for Computational Molecular Biology and Bioinformatics. *Bioinformatics* 2009, 25(11), 1422–1423. [doi:10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)
- **Biotite**: Kunzmann, P.; Hamacher, K. Biotite: A Unifying Open Source Computational Biology Framework in Python. *BMC Bioinformatics* 2018, 19, 346. [doi:10.1186/s12859-018-2367-z](https://doi.org/10.1186/s12859-018-2367-z)

## License

MIT

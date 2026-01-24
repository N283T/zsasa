# freesasa-zig Python Bindings

Python bindings for the freesasa-zig library - a high-performance implementation of Solvent Accessible Surface Area (SASA) calculation.

## Features

- **NumPy integration**: Pass coordinates and radii as NumPy arrays
- **Two algorithms**: Shrake-Rupley (fast) and Lee-Richards (precise)
- **Multi-threading**: Automatic parallelization across CPU cores
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

### Library Location

The Python bindings look for `libfreesasa_zig.dylib` (macOS), `libfreesasa_zig.so` (Linux), or `freesasa_zig.dll` (Windows) in these locations:

1. `FREESASA_ZIG_LIB` environment variable (if set)
2. `../zig-out/lib/` (relative to package)
3. `/usr/local/lib/`
4. `/usr/lib/`
5. Current working directory

## Quick Start

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

## API Reference

### `calculate_sasa(coords, radii, *, algorithm="sr", n_points=100, n_slices=20, probe_radius=1.4, n_threads=0)`

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

**Raises:**
- `ValueError`: Invalid input (wrong shape, negative radii, invalid parameters)
- `RuntimeError`: Calculation error
- `MemoryError`: Out of memory

### `SasaResult`

Dataclass containing calculation results.

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_area` | `float` | Total SASA in Å² |
| `atom_areas` | `NDArray[float64]` | Per-atom SASA values in Å² |

### `get_version()`

Returns the library version string.

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

### Working with Protein Structures

```python
import json
import numpy as np
from freesasa_zig import calculate_sasa

# Load input JSON (from cif_to_input_json.py)
with open("input.json") as f:
    data = json.load(f)

coords = np.column_stack([data["x"], data["y"], data["z"]])
radii = np.array(data["r"])

result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
```

## Performance

Library-to-library comparison (Python bindings vs FreeSASA Python):

| Structure | Atoms | Zig SR | FreeSASA SR | LR Speedup |
|-----------|------:|-------:|------------:|----------:|
| 1CRN | 327 | 1.9ms | 0.7ms | 0.8x |
| 1UBQ | 602 | 2.3ms | 1.2ms | 1.2x |
| 1A0Q | 3,183 | 10.6ms | 7.5ms | 1.4x |
| 3HHB | 4,384 | 14.4ms | 11.0ms | 1.4x |
| 1AON | 58,674 | 201ms | 163ms | 1.4x |

- **SR algorithm**: FreeSASA's C extension is slightly faster
- **LR algorithm**: Zig is 1.2-1.4x faster
- **Accuracy**: Results are identical (0.00 Å² difference)

Run benchmark: `./scripts/benchmark_python.py`

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

## License

MIT

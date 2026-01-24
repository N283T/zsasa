# freesasa-zig

Python bindings for the freesasa-zig library - a high-performance implementation of Solvent Accessible Surface Area (SASA) calculation.

## Installation

```bash
pip install freesasa-zig
```

**Note:** You need to have the `libfreesasa_zig` shared library built and available. Set the `FREESASA_ZIG_LIB` environment variable to the library path if needed.

## Usage

```python
import numpy as np
from freesasa_zig import calculate_sasa

# Define atom coordinates and radii
coords = np.array([
    [0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
])
radii = np.array([1.5, 1.5])

# Calculate SASA using Shrake-Rupley algorithm
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
print(f"Per-atom areas: {result.atom_areas}")

# Use Lee-Richards algorithm
result_lr = calculate_sasa(coords, radii, algorithm="lr")
print(f"Total SASA (LR): {result_lr.total_area:.2f} Å²")
```

## API

### `calculate_sasa(coords, radii, *, algorithm="sr", n_points=100, n_slices=20, probe_radius=1.4, n_threads=0)`

Calculate Solvent Accessible Surface Area.

**Parameters:**
- `coords`: Atom coordinates as (N, 3) NumPy array
- `radii`: Atom radii as (N,) NumPy array
- `algorithm`: "sr" (Shrake-Rupley) or "lr" (Lee-Richards)
- `n_points`: Test points per atom (SR algorithm, default: 100)
- `n_slices`: Slices per atom (LR algorithm, default: 20)
- `probe_radius`: Water probe radius in Å (default: 1.4)
- `n_threads`: Number of threads (0 = auto-detect)

**Returns:**
- `SasaResult` with `total_area` and `atom_areas` attributes

### `get_version()`

Get the library version string.

## Building from Source

```bash
# Build the Zig library
zig build -Doptimize=ReleaseFast

# Install Python package
cd python
pip install -e .
```

## License

MIT

# Core API

Basic SASA calculation functions.

## calculate_sasa

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
    use_bitmask: bool = False,
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
| `use_bitmask` | `bool` | `False` | Use bitmask LUT optimization (SR only, n_points must be 1..1024) |

**Returns:** `SasaResult`

**Raises:** `ValueError` for invalid input

## SasaResult

```python
@dataclass
class SasaResult:
    total_area: float           # Total SASA in Å²
    atom_areas: NDArray[float64] # Per-atom SASA values
```

---

## Batch API

For trajectory analysis with multiple frames.

### calculate_sasa_batch

```python
def calculate_sasa_batch(
    coordinates: NDArray[np.floating],
    radii: NDArray[np.floating],
    *,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    precision: Literal["f64", "f32"] = "f64",
    use_bitmask: bool = False,
) -> BatchSasaResult
```

Calculate SASA for multiple frames in batch.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `coordinates` | `NDArray[floating]` | required | Coordinates as (n_frames, n_atoms, 3) array |
| `radii` | `NDArray[floating]` | required | Atom radii as (n_atoms,) array |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Shrake-Rupley or Lee-Richards |
| `n_points` | `int` | `100` | Test points per atom (SR only) |
| `n_slices` | `int` | `20` | Slices per atom (LR only) |
| `probe_radius` | `float` | `1.4` | Water probe radius in Å |
| `n_threads` | `int` | `0` | Number of threads (0 = auto) |
| `precision` | `"f64"` or `"f32"` | `"f64"` | Floating-point precision |
| `use_bitmask` | `bool` | `False` | Use bitmask LUT optimization (SR only, n_points must be 1..1024) |

**Returns:** `BatchSasaResult`

### BatchSasaResult

```python
@dataclass
class BatchSasaResult:
    atom_areas: NDArray[float32]  # Per-atom SASA, shape (n_frames, n_atoms)

    # Properties
    n_frames: int                  # Number of frames
    n_atoms: int                   # Number of atoms
    total_areas: NDArray[float32]  # Total SASA per frame, shape (n_frames,)
```

### Precision Parameter

The `precision` parameter controls the floating-point precision used internally:

| Precision | Description | Use Case |
|-----------|-------------|----------|
| `"f64"` (default) | 64-bit double precision | Higher accuracy, recommended |
| `"f32"` | 32-bit single precision | Comparison with RustSASA/mdsasa-bolt |

**Note:** The difference between f64 and f32 is typically < 0.01% for SASA calculations. f64 is recommended unless you need exact comparison with f32-based implementations.

**Example:**

```python
import numpy as np
from zsasa import calculate_sasa_batch

# Multiple frames
n_frames = 100
n_atoms = 1000
coords = np.random.rand(n_frames, n_atoms, 3).astype(np.float32) * 50
radii = np.full(n_atoms, 1.5, dtype=np.float32)

# f64 precision (default)
result = calculate_sasa_batch(coords, radii, n_threads=4)
print(f"Shape: {result.atom_areas.shape}")  # (100, 1000)

# f32 precision (for comparison with RustSASA)
result_f32 = calculate_sasa_batch(coords, radii, precision="f32")
```

---

## Directory Batch Processing

Process all supported structure files in a directory at once.

### process_directory

```python
def process_directory(
    input_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
    classifier: ClassifierType | None = ClassifierType.PROTOR,
    include_hydrogens: bool = False,
    include_hetatm: bool = False,
) -> BatchDirResult
```

Process all supported structure files (.pdb, .cif, .mmcif, .ent, .json, and .gz variants) in a directory.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `input_dir` | `str \| Path` | required | Directory containing structure files |
| `output_dir` | `str \| Path \| None` | `None` | Output directory for per-file results |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Shrake-Rupley or Lee-Richards |
| `n_points` | `int` | `100` | Test points per atom (SR only) |
| `n_slices` | `int` | `20` | Slices per atom (LR only) |
| `probe_radius` | `float` | `1.4` | Water probe radius in Å |
| `n_threads` | `int` | `0` | Number of threads (0 = auto) |
| `classifier` | `ClassifierType \| None` | `PROTOR` | Classifier for radii. `None` = use input radii |
| `include_hydrogens` | `bool` | `False` | Include hydrogen atoms |
| `include_hetatm` | `bool` | `False` | Include HETATM records |

**Returns:** `BatchDirResult`

**Raises:**

| Exception | Cause |
|-----------|-------|
| `ValueError` | Invalid parameters (algorithm, probe_radius, n_points, etc.) |
| `FileNotFoundError` | Directory does not exist or is not readable |
| `MemoryError` | Out of memory |
| `RuntimeError` | Calculation or other processing error |

### BatchDirResult

```python
@dataclass
class BatchDirResult:
    total_files: int        # Number of supported files found
    successful: int         # Number successfully processed
    failed: int             # Number that failed
    filenames: list[str]    # Per-file basenames
    n_atoms: list[int]      # Per-file atom counts
    total_sasa: list[float] # Per-file SASA in Å² (NaN if failed)
    status: list[int]       # Per-file status (1=ok, 0=failed)
```

**Example:**

```python
from zsasa import process_directory

# Process all structure files in a directory
result = process_directory("path/to/structures/")
print(f"Processed {result.successful}/{result.total_files} files")

# Per-file results
for i in range(result.total_files):
    if result.status[i] == 1:
        print(f"  {result.filenames[i]}: {result.n_atoms[i]} atoms, "
              f"SASA={result.total_sasa[i]:.1f} Å²")

# With output files and LR algorithm
result = process_directory(
    "input/",
    output_dir="output/",
    algorithm="lr",
    n_threads=4,
)
```

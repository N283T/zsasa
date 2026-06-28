# Legacy Native XTC Reader

The `zsasa.xtc` module provides a standalone XTC trajectory reader kept for compatibility with existing zsasa users. For new Python workflows that read trajectory files directly, prefer [pyztraj](https://github.com/N283T/ztraj), which centralizes trajectory I/O and trajectory-native analysis across XTC, TRR, DCD, and AMBER NetCDF.

## When to Use

| Use Case | Recommended Module |
|----------|-------------------|
| Direct trajectory-file I/O or trajectory-native analysis | `pyztraj` |
| Existing code already using `zsasa.xtc` | `zsasa.xtc` |
| Need MDTraj/MDAnalysis ecosystem objects | `zsasa.mdtraj` or `zsasa.mdanalysis` |

## Compatibility Status

`zsasa.xtc` remains supported for existing users, but new trajectory formats are not planned for zsasa Python modules. Use pyztraj for new direct trajectory-file workflows.

## XtcReader

Low-level XTC file reader with iterator support.

### Constructor

```python
XtcReader(path: str | Path)
```

Opens an XTC file for reading.

**Parameters:**
- `path`: Path to XTC trajectory file

**Raises:**
- `FileNotFoundError`: If file doesn't exist
- `RuntimeError`: If file is invalid or corrupted

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `natoms` | `int` | Number of atoms in trajectory |

### Methods

#### read_frame

```python
def read_frame(self) -> XtcFrame | None
```

Read the next frame from the trajectory.

**Returns:**
- `XtcFrame` object, or `None` if end of file reached

**Raises:**
- `RuntimeError`: If reader is closed or read error occurs

#### close

```python
def close(self) -> None
```

Close the file and release resources. Safe to call multiple times.

### Example: Basic Reading

```python
from zsasa.xtc import XtcReader

# Using context manager (recommended)
with XtcReader("trajectory.xtc") as reader:
    print(f"Trajectory has {reader.natoms} atoms")

    for frame in reader:
        print(f"Step {frame.step}, time {frame.time} ps")
        print(f"First atom: {frame.coords[0]}")
```

### Example: Manual Control

```python
from zsasa.xtc import XtcReader

reader = XtcReader("trajectory.xtc")
try:
    frame = reader.read_frame()
    while frame is not None:
        # Process frame...
        frame = reader.read_frame()
finally:
    reader.close()
```

---

## XtcFrame

Represents a single frame from an XTC trajectory.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `step` | `int` | Simulation step number |
| `time` | `float` | Simulation time in picoseconds |
| `coords` | `NDArray[float32]` | Coordinates as (n_atoms, 3) array in **nanometers** |
| `box` | `NDArray[float32]` | Box matrix as 3x3 array in **nanometers** |
| `precision` | `float` | XTC compression precision |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `natoms` | `int` | Number of atoms |

---

## compute_sasa_trajectory

High-level function for SASA calculation on XTC trajectories.

```python
def compute_sasa_trajectory(
    xtc_path: str | Path,
    radii: NDArray[np.floating] | list[float],
    *,
    probe_radius: float = 1.4,
    n_points: int = 100,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    use_bitmask: bool = False,
    bitmask_correction: bool = False,
    bitmask_correction_coeff: float | None = None,
) -> TrajectorySasaResult
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `xtc_path` | `str \| Path` | Required | Path to XTC file |
| `radii` | `array-like` | Required | Atomic radii in Angstroms (n_atoms,) |
| `probe_radius` | `float` | `1.4` | Water probe radius in Angstroms |
| `n_points` | `int` | `100` | Test points per atom (SR algorithm) |
| `algorithm` | `str` | `"sr"` | `"sr"` (Shrake-Rupley) or `"lr"` (Lee-Richards) |
| `n_slices` | `int` | `20` | Slices per atom (LR algorithm) |
| `n_threads` | `int` | `0` | Thread count (0 = auto-detect) |
| `start` | `int` | `0` | First frame to process |
| `stop` | `int \| None` | `None` | Stop before this frame (None = all) |
| `step` | `int` | `1` | Process every Nth frame |
| `use_bitmask` | `bool` | `False` | Use [bitmask LUT optimization](../guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points must be 1..1024) |
| `bitmask_correction` | `bool` | `False` | Experimental exposed-fraction correction for bitmask quantization bias; requires `use_bitmask=True` |
| `bitmask_correction_coeff` | `float \| None` | `None` | Override the experimental correction coefficient (`None` uses library default) |

**Returns:** `TrajectorySasaResult`

**Raises:**
- `FileNotFoundError`: If XTC file doesn't exist
- `ValueError`: If radii length doesn't match trajectory atoms

### Unit Conversion

XTC coordinates are in **nanometers** (GROMACS convention). This function automatically converts to **Angstroms** for SASA calculation. Output SASA values are in **Angstroms²**.

### Example: Basic Usage

```python
import numpy as np
from zsasa.xtc import compute_sasa_trajectory

# Define radii (must match trajectory atom count)
# Typically you'd get these from a topology file
radii = np.full(304, 1.7)  # 1.7 Å for all atoms (simplified)

result = compute_sasa_trajectory("trajectory.xtc", radii)

print(f"Frames: {result.n_frames}")
print(f"Total SASA per frame: {result.total_areas}")
```

### Example: With Frame Selection

```python
# Process every 10th frame, starting from frame 100
result = compute_sasa_trajectory(
    "trajectory.xtc",
    radii,
    start=100,
    stop=1000,
    step=10,
)
```

### Example: With Topology Radii

```python
import numpy as np
from zsasa import classify_atoms
from zsasa.xtc import compute_sasa_trajectory

# Get radii from topology (e.g., from a PDB file)
# This example assumes you have residue and atom names
residues = ["ALA", "ALA", "ALA", ...]
atoms = ["N", "CA", "C", ...]

classification = classify_atoms(residues, atoms)
radii = classification.radii

# Handle unknown atoms
radii = np.where(np.isnan(radii), 1.7, radii)  # Default to 1.7 Å

result = compute_sasa_trajectory("trajectory.xtc", radii)
```

---

## TrajectorySasaResult

Result container for trajectory SASA calculation.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `atom_areas` | `NDArray[float32]` | Per-atom SASA, shape (n_frames, n_atoms) in Å² |
| `steps` | `NDArray[int32]` | Step numbers for each frame |
| `times` | `NDArray[float32]` | Time values in picoseconds |

### Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_frames` | `int` | Number of frames |
| `n_atoms` | `int` | Number of atoms |
| `total_areas` | `NDArray[float32]` | Total SASA per frame, shape (n_frames,) |

---

## Comparison with MDTraj/MDAnalysis Integration

| Feature | pyztraj | zsasa.xtc | zsasa.mdtraj | zsasa.mdanalysis |
|---------|---------|-----------|--------------|------------------|
| Dependencies | pyztraj | None (only NumPy) | mdtraj | MDAnalysis |
| Trajectory formats | XTC, TRR, DCD, AMBER NetCDF | XTC only | Many (XTC, TRR, DCD, ...) | Many |
| Topology support | From ztraj loaders | Manual radii | From topology | From topology |
| Atom selection | Yes | No | Yes | Yes |
| Status | Preferred for direct trajectory files | Legacy compatibility | Ecosystem integration | Ecosystem integration |

### When to Use pyztraj

- You need direct trajectory-file I/O in Python
- You work with TRR, DCD, AMBER NetCDF, or multiple trajectory formats
- You want trajectory-native analysis such as SASA near the I/O layer

### When to Use zsasa.xtc

- You are maintaining existing code that already imports `zsasa.xtc`
- You only need the legacy XTC compatibility API

### When to Use MDTraj/MDAnalysis

- You need MDTraj or MDAnalysis topology/selection objects
- You already use those ecosystems for upstream trajectory handling

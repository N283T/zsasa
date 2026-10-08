# Legacy Native XTC Reader

The `zsasa.xtc` module provides a standalone XTC trajectory reader kept for compatibility with existing zsasa users; the `zsasa.dcd` module offers the same API for DCD files (see [DCD](#dcd-zsasadcd)). For new Python workflows that read trajectory files directly, prefer [pyztraj](https://github.com/N283T/ztraj), which centralizes trajectory I/O and trajectory-native analysis across XTC, TRR, DCD, and AMBER NetCDF.

## When to Use

| Use Case | Recommended Module |
|----------|-------------------|
| Direct trajectory-file I/O or trajectory-native analysis | `pyztraj` |
| Existing code already using `zsasa.xtc` or `zsasa.dcd` | `zsasa.xtc` / `zsasa.dcd` |
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
- `FileNotFoundError`: If the file doesn't exist
- `IsADirectoryError`: If `path` is a directory
- `PermissionError`: If the file cannot be read
- `ValueError`: If the file exists but is not a valid XTC file: empty, truncated, or in another format (for example a DCD file)
- `MemoryError`: If memory cannot be allocated

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
- `RuntimeError`: If the reader is closed, or the frame is corrupt or truncated
- `MemoryError`: If memory for the frame cannot be allocated

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
| `start` | `int` | `0` | First frame to process (non-negative) |
| `stop` | `int \| None` | `None` | Stop before this frame (non-negative; None = all) |
| `step` | `int` | `1` | Process every Nth frame (positive) |
| `use_bitmask` | `bool` | `False` | Use [bitmask LUT optimization](../guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points must be 1..1024) |
| `bitmask_correction` | `bool` | `False` | Experimental exposed-fraction correction for bitmask quantization bias; requires `use_bitmask=True` |
| `bitmask_correction_coeff` | `float \| None` | `None` | Override the experimental correction coefficient (`None` uses library default) |

**Returns:** `TrajectorySasaResult`

**Raises:**
- `ValueError`: If `step` is less than 1 or `start` or `stop` is negative (checked before the file is opened), if the file is not a valid XTC file, if radii length doesn't match trajectory atoms, or if no frame is selected
- `FileNotFoundError`: If the XTC file doesn't exist

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

## compute_sasa_trajectory_summary

Memory-efficient variant of `compute_sasa_trajectory`: frames are read and processed in chunks, and the per-atom SASA arrays are discarded after the per-frame totals (and optional per-residue sums) have been accumulated. Use it for long trajectories when you do not need per-atom values.

```python
def compute_sasa_trajectory_summary(
    xtc_path: str | Path,
    radii: NDArray[np.floating] | list[float],
    *,
    atom_to_residue: NDArray[np.integer] | list[int] | None = None,
    probe_radius: float = 1.4,
    n_points: int = 100,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    chunk_size: int = 16,
    use_bitmask: bool = False,
    bitmask_correction: bool = False,
    bitmask_correction_coeff: float | None = None,
) -> TrajectorySasaSummaryResult
```

It takes the parameters of [`compute_sasa_trajectory`](#compute_sasa_trajectory) and two more:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `atom_to_residue` | `array-like of int \| None` | `None` | Residue index of each atom, shape `(n_atoms,)`, non-negative. When given, per-residue SASA is returned too. The indices are remapped to consecutive columns in ascending order of the original index |
| `chunk_size` | `int` | `16` | Frames processed per native batch (must be positive). Only one chunk of per-atom areas is held in memory at a time |

**Returns:** `TrajectorySasaSummaryResult`

**Raises:**
- `ValueError`: If `radii` does not match the trajectory atoms, `chunk_size` is not positive, `step` is less than 1 or `start` or `stop` is negative, `atom_to_residue` has the wrong shape or negative values, or no frame is selected

### TrajectorySasaSummaryResult

| Attribute | Type | Description |
|-----------|------|-------------|
| `total_areas` | `NDArray[float32]` | Total SASA per frame, shape (n_frames,) in Å² |
| `steps` | `NDArray[int32]` | Step numbers for each frame |
| `times` | `NDArray[float32]` | Time values in picoseconds |
| `residue_areas` | `NDArray[float32] \| None` | Per-residue SASA, shape (n_frames, n_residues) in Å²; `None` without `atom_to_residue` |

Properties: `n_frames` and `n_residues` (`0` when `residue_areas` is `None`).

```python
import numpy as np
from zsasa.xtc import compute_sasa_trajectory_summary

radii = np.full(304, 1.7)
atom_to_residue = np.repeat(np.arange(19), 16)[:304]  # one residue index per atom

result = compute_sasa_trajectory_summary(
    "trajectory.xtc", radii, atom_to_residue=atom_to_residue, chunk_size=8
)
print(result.total_areas.shape)     # (n_frames,)
print(result.residue_areas.shape)   # (n_frames, n_residues)
```

---

## DCD: zsasa.dcd

`zsasa.dcd` provides the same reader and SASA functions for DCD files (NAMD/CHARMM). Unlike XTC, DCD coordinates are already in **Angstroms**, so no unit conversion is applied.

```python
from zsasa.dcd import (
    DcdFrame,
    DcdReader,
    compute_sasa_trajectory,
    compute_sasa_trajectory_summary,
)
```

| Name | Description |
|------|-------------|
| `DcdReader(path)` | Frame reader with the same interface as `XtcReader`: `natoms`, `read_frame()`, `close()`, iteration and context-manager use. Raises the same exceptions as `XtcReader` (`FileNotFoundError` for a missing file, `ValueError` for an empty, truncated or non-DCD file) |
| `DcdFrame` | One frame: `step`, `time`, `coords` (`NDArray[float32]`, (n_atoms, 3), Å), `unitcell` (six doubles, or `None` when the file has none) and the `natoms` property |
| `compute_sasa_trajectory(dcd_path, radii, ...)` | Same parameters and `TrajectorySasaResult` as the XTC function; `dcd_path` replaces `xtc_path` |
| `compute_sasa_trajectory_summary(dcd_path, radii, ...)` | Same parameters and `TrajectorySasaSummaryResult` as the XTC function, including `atom_to_residue` and `chunk_size` |

`TrajectorySasaResult` and `TrajectorySasaSummaryResult` are the classes of `zsasa.xtc`.

```python
import numpy as np
from zsasa.dcd import DcdReader, compute_sasa_trajectory

with DcdReader("trajectory.dcd") as reader:
    print(reader.natoms)
    frame = reader.read_frame()
    print(frame.step, frame.coords.shape)

radii = np.full(304, 1.7)
result = compute_sasa_trajectory("trajectory.dcd", radii, step=2)
print(result.total_areas)
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
| Trajectory formats | XTC, TRR, DCD, AMBER NetCDF | XTC only (`zsasa.dcd` for DCD) | Many (XTC, TRR, DCD, ...) | Many |
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

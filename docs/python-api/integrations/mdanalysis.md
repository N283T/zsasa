# MDAnalysis Integration

High-performance SASA analysis compatible with MDAnalysis' `AnalysisBase` pattern.

## Installation

```bash
pip install MDAnalysis
```

## Import

```python
from freesasa_zig.mdanalysis import SASAAnalysis, compute_sasa
```

## Overview

**3.4x faster** than mdsasa-bolt on real MD data with controllable thread count.

## SASAAnalysis

Class-based API following MDAnalysis conventions.

```python
class SASAAnalysis:
    """SASA analysis for MDAnalysis trajectories."""

    def __init__(
        self,
        universe_or_atomgroup: Universe | AtomGroup,
        select: str = "all",
    ) -> None
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `universe_or_atomgroup` | `Universe` or `AtomGroup` | required | MDAnalysis object to analyze |
| `select` | `str` | `"all"` | Atom selection string |

**Attributes:**

| Attribute | Type | Description |
|-----------|------|-------------|
| `atomgroup` | `AtomGroup` | The atoms being analyzed |
| `results` | `Results` | Results object (after `run()`) |
| `n_frames` | `int` | Number of frames analyzed |
| `times` | `NDArray[float64]` | Frame times |
| `frames` | `NDArray[int64]` | Frame indices |

### run()

```python
def run(
    self,
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    *,
    probe_radius: float = 1.4,
    n_points: int = 960,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
) -> SASAAnalysis
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `start` | `int` | `0` | First frame to analyze |
| `stop` | `int` | `None` | Last frame (None = last) |
| `step` | `int` | `1` | Step between frames |
| `probe_radius` | `float` | `1.4` | Probe radius in Å |
| `n_points` | `int` | `960` | Test points per atom (SR) |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Algorithm to use |
| `n_slices` | `int` | `20` | Slices per atom (LR) |
| `n_threads` | `int` | `0` | Threads (0 = auto) |

**Returns:** `self` for method chaining.

### Results

After calling `run()`, results are available in the `results` attribute:

| Attribute | Type | Description |
|-----------|------|-------------|
| `atom_area` | `NDArray[float32]` | Per-atom SASA, shape `(n_frames, n_atoms)` |
| `residue_area` | `NDArray[float32]` | Per-residue SASA, shape `(n_frames, n_residues)` |
| `total_area` | `NDArray[float32]` | Total SASA, shape `(n_frames,)` |
| `mean_total_area` | `float` | Mean total SASA across all frames |

**Units:** All SASA values are in Å² (matching MDAnalysis conventions).

## Example

```python
import MDAnalysis as mda
from freesasa_zig.mdanalysis import SASAAnalysis

# Load trajectory
u = mda.Universe("topology.pdb", "trajectory.xtc")

# Analyze protein only
sasa = SASAAnalysis(u, select="protein")
sasa.run(start=0, stop=100, step=10)

# Access results
print(f"Mean SASA: {sasa.results.mean_total_area:.2f} Å²")
print(f"Per-frame: {sasa.results.total_area}")
print(f"Per-residue shape: {sasa.results.residue_area.shape}")
print(f"Per-atom shape: {sasa.results.atom_area.shape}")
```

## compute_sasa (function)

Convenience function for simple use cases.

```python
def compute_sasa(
    universe_or_atomgroup: Universe | AtomGroup,
    *,
    select: str = "all",
    start: int = 0,
    stop: int | None = None,
    step: int = 1,
    probe_radius: float = 1.4,
    n_points: int = 960,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    mode: Literal["atom", "residue", "total"] = "atom",
) -> NDArray[np.float32]
```

**Returns:** SASA values in Å².

| Mode | Shape |
|------|-------|
| `"atom"` | `(n_frames, n_atoms)` |
| `"residue"` | `(n_frames, n_residues)` |
| `"total"` | `(n_frames,)` |

**Example:**

```python
from freesasa_zig.mdanalysis import compute_sasa
import MDAnalysis as mda

u = mda.Universe("topology.pdb", "trajectory.xtc")

# Simple one-liner
total_sasa = compute_sasa(u, select="protein", mode="total")
```

## Key Advantages

- **Controllable parallelism**: Set exact thread count with `n_threads`
- **Higher precision**: f64 by default
- **Efficient scaling**: 5.2x speedup from 1→8 threads
- **MDAnalysis compatible**: Works with selections, AtomGroups

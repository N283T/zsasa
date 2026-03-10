# MDTraj Integration

High-performance MD trajectory SASA analysis.

## Installation

```bash
pip install mdtraj
# or
uv add mdtraj
```

## Import

```python
from zsasa.mdtraj import compute_sasa, shrake_rupley
```

## Overview

A drop-in replacement for `mdtraj.shrake_rupley()` with better performance through SIMD optimization and frame-level parallelization.

**3.4x faster** than mdsasa-bolt on real MD data.

## compute_sasa

```python
def compute_sasa(
    traj: md.Trajectory,
    *,
    probe_radius: float = 1.4,
    n_points: int = 960,
    algorithm: Literal["sr", "lr"] = "sr",
    n_slices: int = 20,
    n_threads: int = 0,
    mode: Literal["atom", "residue", "total"] = "atom",
    use_bitmask: bool = False,
) -> NDArray[np.float32]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `traj` | `md.Trajectory` | required | MDTraj Trajectory object |
| `probe_radius` | `float` | `1.4` | Water probe radius in Å (MDTraj uses 0.14 nm) |
| `n_points` | `int` | `960` | Test points per atom (SR) |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Algorithm to use |
| `n_slices` | `int` | `20` | Slices per atom (LR) |
| `n_threads` | `int` | `0` | Threads (0 = auto) |
| `mode` | `"atom"`, `"residue"`, `"total"` | `"atom"` | Output mode |
| `use_bitmask` | `bool` | `False` | Use [bitmask LUT optimization](../guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points must be 1..1024) |

**Returns:** SASA values in nm² (matching MDTraj's output units).

| Mode | Shape |
|------|-------|
| `"atom"` | `(n_frames, n_atoms)` |
| `"residue"` | `(n_frames, n_residues)` |
| `"total"` | `(n_frames,)` |

**Note:** Coordinates are automatically converted from nm (MDTraj) to Angstrom internally, and output is converted back to nm².

## Example

```python
import mdtraj as md
from zsasa.mdtraj import compute_sasa

# Load trajectory
traj = md.load('trajectory.xtc', top='topology.pdb')

# Per-atom SASA (default)
sasa_atom = compute_sasa(traj, mode='atom')
print(f"Shape: {sasa_atom.shape}")  # (n_frames, n_atoms)

# Per-residue SASA
sasa_residue = compute_sasa(traj, mode='residue')

# Total SASA per frame
sasa_total = compute_sasa(traj, mode='total')

# Use Lee-Richards algorithm
sasa_lr = compute_sasa(traj, algorithm='lr', n_slices=40)

# Control thread count
sasa = compute_sasa(traj, n_threads=8)
```

## shrake_rupley (alias)

`shrake_rupley` is an alias for `compute_sasa` for compatibility with MDTraj's naming convention.

```python
from zsasa.mdtraj import shrake_rupley

# Same as compute_sasa
sasa = shrake_rupley(traj)
```

## Migration from MDTraj

```python
# Before (MDTraj native)
import mdtraj as md
sasa = md.shrake_rupley(traj, probe_radius=0.14)  # nm

# After (zsasa)
from zsasa.mdtraj import compute_sasa
sasa = compute_sasa(traj, probe_radius=1.4)  # Å (same value, different units)
```

Output units are the same (nm²).

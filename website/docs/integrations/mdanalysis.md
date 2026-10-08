# MDAnalysis Integration

High-performance SASA analysis compatible with MDAnalysis' `AnalysisBase` pattern.

## Installation

```bash
pip install MDAnalysis
# or
uv add MDAnalysis
```

## Import

```python
from zsasa.mdanalysis import SASAAnalysis, compute_sasa
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
    chunk_size: int | None = None,
    store_atom_areas: bool = True,
    use_bitmask: bool = False,
    bitmask_correction: bool = False,
    bitmask_correction_coeff: float | None = None,
) -> SASAAnalysis
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `start` | `int` | `0` | First frame to analyze (not negative) |
| `stop` | `int \| None` | `None` | Stop before this frame (not negative; None = run through the last frame) |
| `step` | `int` | `1` | Step between frames (at least 1; `0` or a negative value raises `ValueError` before any frame is read) |
| `probe_radius` | `float` | `1.4` | Probe radius in Å |
| `n_points` | `int` | `960` | Test points per atom (SR) |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Algorithm to use |
| `n_slices` | `int` | `20` | Slices per atom (LR) |
| `n_threads` | `int` | `0` | Threads (0 = auto) |
| `chunk_size` | `int \| None` | `None` | Frames per native batch (must be positive). `None` processes all selected frames in one batch; a smaller value lowers peak memory |
| `store_atom_areas` | `bool` | `True` | Keep per-atom SASA in `results.atom_area`. With `False`, only the totals and per-residue sums are kept and `results.atom_area` is `None` |
| `use_bitmask` | `bool` | `False` | Use [bitmask LUT optimization](../guide/algorithms.mdx#bitmask-lut-optimization) (SR only, n_points must be 1..1024) |
| `bitmask_correction` | `bool` | `False` | Experimental exposed-fraction correction for bitmask quantization bias; requires `use_bitmask=True` |
| `bitmask_correction_coeff` | `float \| None` | `None` | Override the experimental correction coefficient (`None` uses library default) |

**Returns:** `self` for method chaining.

### Results

After calling `run()`, results are available in the `results` attribute:

| Attribute | Type | Description |
|-----------|------|-------------|
| `atom_area` | `NDArray[float32] \| None` | Per-atom SASA, shape `(n_frames, n_atoms)`; `None` when `run()` was called with `store_atom_areas=False` |
| `residue_area` | `NDArray[float32]` | Per-residue SASA, shape `(n_frames, n_residues)` |
| `total_area` | `NDArray[float32]` | Total SASA, shape `(n_frames,)` |
| `mean_total_area` | `float` | Mean total SASA across all frames |

**Units:** All SASA values are in Å² (matching MDAnalysis conventions).

### Radii

The radius of each atom comes from the MDAnalysis van der Waals table (`MDAnalysis.guesser.tables.vdwradii`; `MDAnalysis.topology.tables.vdwradii` in MDAnalysis releases before 2.8, which zsasa falls back to) and the atom's element. An element missing from that table (for example Fe, which MDAnalysis does not list) gets 2.0 Å.

The element is the `element` attribute when the topology has one (for example a PDB file with an element column). Without it (for example GRO files, or PDB files without an element column), it is inferred from the atom type and name together with the residue name:

- A name or type that starts with `CL`, `BR`, `FE`, `ZN`, `MG`, `MN`, `NI`, `CU`, `LI` or `AL` is that element (`CL1` is chlorine, `ZN` is zinc).
- A two-letter symbol that is also an ordinary atom name (`CA`, `CD`, `NA`, `HG`, `SE`, `CO`, ...) is that element only when the residue is the ion itself (atom `CA` in residue `CA`, atom `NA` in residue `NA+`). The alpha carbon `CA` of `ALA`, the `CD` of `GLN` and the `NA` nitrogen of `HEM` stay carbon, carbon and nitrogen, and `HG` of `SER` stays hydrogen.
- CHARMM ion names (`SOD`, `POT`, `CLA`, `CAL`, `CES`) give sodium, potassium, chlorine, calcium and caesium.
- The type and the name are both read, and a two-letter reading from either wins: a parser may guess the type of an ion from its name and get it wrong (MDAnalysis types the calcium ion `CA` as carbon).

The atom mass is not used. MDAnalysis guesses masses from the same types when a topology has none, so a guessed mass cannot confirm a guessed type.

Anything else takes the first letter of the type, or the name when there is no type, and carbon when there is neither.

## Example

```python
import MDAnalysis as mda
from zsasa.mdanalysis import SASAAnalysis

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
    chunk_size: int | None = None,
    mode: Literal["atom", "residue", "total"] = "atom",
    use_bitmask: bool = False,
    bitmask_correction: bool = False,
    bitmask_correction_coeff: float | None = None,
) -> NDArray[np.float32]
```

`chunk_size` has the same meaning as in `run()`. Per-atom areas are only kept for `mode="atom"`, so `"residue"` and `"total"` need less memory.

**Returns:** SASA values in Å².

| Mode | Shape |
|------|-------|
| `"atom"` | `(n_frames, n_atoms)` |
| `"residue"` | `(n_frames, n_residues)` |
| `"total"` | `(n_frames,)` |

**Example:**

```python
from zsasa.mdanalysis import compute_sasa
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

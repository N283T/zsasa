---
sidebar_position: 4
---

# Trajectory Analysis

zsasa supports SASA calculation over MD trajectory frames using the CLI or Python bindings.

Supported trajectory formats:
- **XTC** (GROMACS) — compressed, coordinates in nm (auto-converted to Å)
- **DCD** (NAMD/CHARMM) — uncompressed, coordinates in Å

Format is auto-detected from file extension.

## CLI: `traj` Subcommand

```bash
zsasa traj <trajectory> <topology> [OPTIONS]
```

The topology file (PDB or mmCIF) provides atom names and radii.

### Example

```bash
# XTC trajectory
zsasa traj trajectory.xtc topology.pdb

# DCD trajectory (NAMD/CHARMM)
zsasa traj trajectory.dcd topology.pdb

# With classifier and frame selection
zsasa traj trajectory.xtc topology.pdb \
    --classifier=naccess \
    --stride=10 \
    --start=100 --end=500

# Exclude hydrogens (included by default)
zsasa traj trajectory.xtc topology.pdb --no-hydrogens

# Output to specific file
zsasa traj trajectory.xtc topology.pdb -o sasa_results.csv
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--stride=N` | Process every Nth frame | `1` |
| `--start=N` | Start from frame N | `0` |
| `--end=N` | End at frame N | all |
| `--classifier=TYPE` | `naccess`, `protor`, `oons` | none |
| `--threads=N` | Thread count (0 = auto) | `0` |
| `--precision=P` | `f32` (fast) or `f64` (precise) | `f32` |
| `--no-hydrogens` | Exclude hydrogen atoms | included |
| `--batch-size=N` | Frames per batch (omit for auto) | auto |
| `-o, --output=FILE` | Output CSV file | `traj_sasa.csv` |

### Output Format

CSV with per-frame total SASA:

```csv
frame,step,time,total_sasa
0,1,1.000,1866.44
1,2,2.000,1977.96
2,3,3.000,1884.93
```

## Python: MDAnalysis Integration

```python
import MDAnalysis as mda
from zsasa.mdanalysis import SASAAnalysis

u = mda.Universe("topology.pdb", "trajectory.xtc")
sasa = SASAAnalysis(u, select="protein")
sasa.run()

print(f"Mean SASA: {sasa.results.mean_total_area:.2f} Å²")
print(f"Per-frame: {sasa.results.total_area}")
```

See [MDAnalysis Integration](../integrations/mdanalysis.md) for full API details.

## Python: MDTraj Integration

```python
import mdtraj as md
from zsasa.mdtraj import compute_sasa

traj = md.load("trajectory.xtc", top="topology.pdb")
sasa = compute_sasa(traj)  # returns (n_frames, n_atoms) array
```

A drop-in replacement for `mdtraj.shrake_rupley()`.

See [MDTraj Integration](../integrations/mdtraj.md) for full API details.

## Python: Native XTC Reader

For simple XTC reading without MDTraj or MDAnalysis dependencies:

```python
from zsasa.xtc import XtcReader

reader = XtcReader("trajectory.xtc")
for frame in reader:
    print(f"Step {frame.step}, Time {frame.time:.1f} ps")
    coords = frame.coords  # numpy array (n_atoms, 3) in nm
```

See [Native XTC Reader](../python-api/xtc.md) for full API details.

## Python: Native DCD Reader

For DCD trajectories without external dependencies:

```python
from zsasa.dcd import DcdReader

reader = DcdReader("trajectory.dcd")
for frame in reader:
    print(f"Step {frame.step}, Time {frame.time:.1f}")
    coords = frame.coords  # numpy array (n_atoms, 3) in Å
```

DCD coordinates are already in Angstroms (no unit conversion needed, unlike XTC).

## Choosing an Approach

| Approach | Best For | Formats | Dependencies |
|----------|----------|---------|-------------|
| CLI `traj` | Quick analysis, scripting | XTC, DCD | None (Zig binary) |
| MDAnalysis | Complex selections, multi-format | XTC, DCD, + many more | `MDAnalysis` |
| MDTraj | Drop-in replacement for `mdtraj.shrake_rupley` | XTC, DCD, + many more | `mdtraj` |
| Native XTC | Simple XTC reading, no extra deps | XTC | None |
| Native DCD | Simple DCD reading, no extra deps | DCD | None |

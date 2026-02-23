# Trajectory Analysis

zsasa supports SASA calculation over MD trajectory frames using the CLI or Python bindings.

## CLI: `traj` Subcommand

```bash
zsasa traj <xtc_file> <topology_file> [OPTIONS]
```

The topology file (PDB or mmCIF) provides atom names and radii. XTC coordinates are automatically converted from nm to Å.

### Example

```bash
# Basic trajectory SASA
zsasa traj trajectory.xtc topology.pdb

# With classifier and frame selection
zsasa traj trajectory.xtc topology.pdb \
    --classifier=naccess \
    --stride=10 \
    --start=100 --end=500

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

See [MDAnalysis Integration](../python-api/integrations/mdanalysis.md) for full API details.

## Python: MDTraj Integration

```python
import mdtraj as md
from zsasa.mdtraj import compute_sasa

traj = md.load("trajectory.xtc", top="topology.pdb")
sasa = compute_sasa(traj)  # returns (n_frames, n_atoms) array
```

A drop-in replacement for `mdtraj.shrake_rupley()`.

See [MDTraj Integration](../python-api/integrations/mdtraj.md) for full API details.

## Python: Native XTC Reader

For simple XTC reading without MDTraj or MDAnalysis dependencies:

```python
from zsasa.xtc import XtcReader

reader = XtcReader("trajectory.xtc")
for frame in reader:
    print(f"Step {frame.step}, Time {frame.time:.1f} ps")
    coords = frame.coords  # numpy array (n_atoms, 3)
```

See [Native XTC Reader](../python-api/xtc.md) for full API details.

## Choosing an Approach

| Approach | Best For | Dependencies |
|----------|----------|-------------|
| CLI `traj` | Quick analysis, scripting | None (Zig binary) |
| MDAnalysis | Complex selections, multi-format | `MDAnalysis` |
| MDTraj | Drop-in replacement for `mdtraj.shrake_rupley` | `mdtraj` |
| Native XTC | Simple reading, no extra deps | None |

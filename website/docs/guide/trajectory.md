---
sidebar_position: 4
---

# Trajectory Analysis

zsasa supports SASA calculation over MD trajectory frames using the CLI or Python bindings.

Supported trajectory formats:
- **XTC** (GROMACS) — compressed trajectory
- **TRR** (GROMACS) — full-precision trajectory with coordinates
- **DCD** (NAMD/CHARMM) — uncompressed trajectory
- **AMBER NetCDF** (`.nc`, `.ncdf`) — AMBER convention NetCDF trajectory

Coordinates are normalized to Å internally by the ztraj readers before SASA calculation.

Format is auto-detected from file extension.

## CLI: `traj` Subcommand

```bash
zsasa traj <trajectory> <topology> [OPTIONS]
```

The topology file (PDB or mmCIF) provides atom names and radii. See [Topology Requirements](#topology-requirements) for how it has to match the trajectory.

### Example

```bash
# XTC trajectory
zsasa traj trajectory.xtc topology.pdb

# DCD trajectory (NAMD/CHARMM)
zsasa traj trajectory.dcd topology.pdb

# TRR trajectory (GROMACS)
zsasa traj trajectory.trr topology.pdb

# AMBER NetCDF trajectory
zsasa traj trajectory.nc topology.pdb

# With classifier and frame selection
zsasa traj trajectory.xtc topology.pdb \
    --classifier=naccess \
    --stride=10 \
    --start=100 --end=500

# Exclude hydrogens (included by default)
zsasa traj trajectory.xtc topology.pdb --no-hydrogens

# Lee-Richards instead of Shrake-Rupley
zsasa traj trajectory.xtc topology.pdb --algorithm=lr --n-slices=50

# Output to specific file
zsasa traj trajectory.xtc topology.pdb -o sasa_results.csv
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--stride=N` | Process every Nth frame (N ≥ 1) | `1` |
| `--start=N` | Start from frame N | `0` |
| `--end=N` | End at frame N (inclusive) | all |
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) or `lr` (Lee-Richards) | `sr` |
| `--n-points=N` | Test points per atom (SR, 1-10000) | `100` |
| `--n-slices=N` | Slices per atom diameter (LR, 1-1000) | `20` |
| `--lr-trig=MODE` | LR arc angles: `exact` or `fast` (approximation of zsasa 0.9.1 and earlier) | `exact` |
| `--probe-radius=R` | Probe radius in Å (0 < R ≤ 10) | `1.4` |
| `--classifier=TYPE` | `ccd`, `naccess`, `protor`, `oons` | `naccess` |
| `--threads=N` | Thread count (0 = auto) | `0` |
| `--precision=P` | `f32` (fast) or `f64` (precise) | `f32` |
| `--no-hydrogens` | Exclude hydrogen atoms from the calculation (`--exclude-hydrogens` is a synonym) | included |
| `--include-hydrogens` | Include hydrogen atoms (the default) | included |
| `--ccd=PATH`, `--sdf=PATH` | External CCD dictionary and SDF bond topology for the `ccd` classifier | none |
| `--use-bitmask` | [Bitmask LUT optimization](algorithms.mdx#bitmask-lut-optimization) (SR only, `--n-points` 1-1024) | off |
| `--bitmask-lut-mode=MODE`, `--bitmask-correction`, `--bitmask-correction-coeff=V` | Bitmask variants; see [Algorithms](algorithms.mdx#experimental-bias-correction) | `single`, off, `0.020` |
| `--batch-size=N` | Frames per batch (omit for auto) | auto |
| `-q, --quiet` | Suppress progress output | off |
| `-o FILE`, `--output=FILE` | Output CSV file | `traj_sasa.csv` |

`--algorithm` and `--precision` are independent: Lee-Richards runs at `f32` by default and at `f64` with `--precision=f64`. See [Commands & Options](../cli/commands.md#trajectory-options) for the full list.

### Output Format

CSV with per-frame total SASA:

```csv
frame,step,time,total_sasa
0,1,1.000,1840.88
1,2,2.000,1944.47
2,3,3.000,1848.46
```

These are the first frames of `test_data/1l2y.xtc` with `test_data/1l2y.pdb` as topology, at the defaults (NACCESS classifier, 100 test points, `f32`).

### Topology Requirements

The topology must describe exactly the atoms stored in the trajectory, in the same order. zsasa reads every `ATOM` and `HETATM` record of the **first model** of the topology file and compares the atom count with the trajectory before any frame is processed.

- A multi-model file (for example an NMR ensemble) can be used directly; only its first model is the topology.
- Solvent, ions and ligands that are part of the trajectory must be listed in the topology, and they are part of the reported SASA. To restrict the calculation to a subset, write a trajectory and topology that contain only those atoms, or use the [MDAnalysis integration](../integrations/mdanalysis.md) with a selection.
- For residues the classifier does not know, radii fall back to generic atom-name or element-based values, as in `calc`.
- `--altloc` selects which [alternate locations](../cli/input.md#alternate-locations) of a PDB or mmCIF topology are kept (default `auto`, one per atom). The trajectory must hold the atoms that remain.

### Hydrogens

Hydrogens are included by default, because MD trajectories carry explicit hydrogens and the default `naccess` classifier has radii for them. `--no-hydrogens` removes the hydrogen atoms of the topology and skips their coordinates in every frame, so it works on a normal all-atom trajectory: the trajectory and the topology still have to contain the same atoms, hydrogens included. If the trajectory has already been stripped of hydrogens, pass a hydrogen-free topology and leave the option out.

### Failures

All options are checked before the output file is created, so a rejected command never overwrites an existing results file. If reading or calculating a frame fails, zsasa writes the frames completed before it, reports the error and exits with a non-zero status.

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

## Python: Direct Trajectory I/O with pyztraj

For Python workflows that read trajectory files directly, use [pyztraj](https://github.com/N283T/ztraj). pyztraj centralizes trajectory I/O and trajectory-native analysis for XTC, TRR, DCD, and AMBER NetCDF. zsasa keeps `zsasa.xtc` and `zsasa.dcd` for compatibility, but new trajectory formats are handled in pyztraj.

```python
import pyztraj

structure = pyztraj.load_pdb("topology.pdb")
with pyztraj.open_trr("trajectory.trr", structure.n_atoms) as reader:
    for frame in reader:
        sasa = pyztraj.compute_sasa(structure, frame.coords)
        print(frame.step, sasa.total_area)
```

See [Legacy Native XTC Reader](../python-api/xtc.md) for the compatibility XTC API.

## Choosing an Approach

| Approach | Best For | Formats | Dependencies |
|----------|----------|---------|-------------|
| CLI `traj` | Quick analysis, scripting | XTC, TRR, DCD, AMBER NetCDF | None (Zig binary) |
| MDAnalysis | Complex selections, multi-format | XTC, DCD, + many more | `MDAnalysis` |
| MDTraj | Drop-in replacement for `mdtraj.shrake_rupley` | XTC, DCD, + many more | `mdtraj` |
| pyztraj | Direct trajectory-file I/O and trajectory-native analysis | XTC, TRR, DCD, AMBER NetCDF | `pyztraj` |
| Legacy `zsasa.xtc`/`zsasa.dcd` | Existing compatibility code | XTC, DCD | None |

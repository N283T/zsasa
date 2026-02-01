# zsasa Python Examples

This directory contains example scripts demonstrating various features of the zsasa Python package.

## Quick Start

```bash
# Install zsasa with all integrations
pip install zsasa[gemmi,biopython,biotite]

# Or minimal install
pip install zsasa

# Run an example
python basic_sasa.py
```

## Examples Overview

| Example | Description | Requirements |
|---------|-------------|--------------|
| `basic_sasa.py` | Core SASA calculation with numpy | `zsasa` |
| `classifier.py` | Atom radii and classification | `zsasa` |
| `gemmi_example.py` | SASA from structure files (mmCIF/PDB) | `zsasa[gemmi]` |
| `biopython_example.py` | BioPython structure integration | `zsasa[biopython]` |
| `biotite_example.py` | Biotite/AtomWorks integration | `zsasa[biotite]` |
| `mdtraj_example.py` | MD trajectory analysis with MDTraj | `zsasa`, `mdtraj` |
| `mdanalysis_example.py` | MD trajectory analysis with MDAnalysis | `zsasa`, `MDAnalysis` |

## Example Files

The examples use structure files from the `examples/` directory at the project root:

- `1ubq.pdb` - Ubiquitin (76 residues)
- `1ubq.cif` - Ubiquitin (mmCIF format)
- `1crn.pdb` - Crambin (46 residues)
- `3hhb.cif.gz` - Hemoglobin (compressed mmCIF)

## Core Examples

### basic_sasa.py

Demonstrates fundamental SASA calculation:

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} A^2")
```

Features covered:
- Single atom SASA
- Multi-atom systems
- Algorithm comparison (SR vs LR)
- Multi-threaded calculation

### classifier.py

Demonstrates atom classification and radius lookup:

```python
from zsasa import classify_atoms, ClassifierType

residues = ["ALA", "ALA", "GLY"]
atoms = ["CA", "O", "N"]

result = classify_atoms(residues, atoms, ClassifierType.NACCESS)
print(result.radii)    # Per-atom radii
print(result.classes)  # Polar/apolar classification
```

Features covered:
- Single atom radius lookup
- Batch classification
- Polar/apolar classification
- Element-based radius guessing
- RSA calculation

## Structure File Integrations

### gemmi_example.py (Recommended)

The fastest and most feature-complete structure file integration:

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
print(f"Total: {result.total_area:.1f} A^2")
print(f"Polar: {result.polar_area:.1f} A^2")
print(f"Apolar: {result.apolar_area:.1f} A^2")
```

### biopython_example.py

For users already using BioPython:

```python
from zsasa.integrations.biopython import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.pdb")
```

### biotite_example.py

For users using Biotite or AtomWorks:

```python
from zsasa.integrations.biotite import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
```

## MD Trajectory Analysis

### mdtraj_example.py

Drop-in replacement for MDTraj's `shrake_rupley`:

```python
import mdtraj as md
from zsasa.mdtraj import compute_sasa

traj = md.load('trajectory.xtc', top='topology.pdb')
sasa = compute_sasa(traj)  # Shape: (n_frames, n_atoms), units: nm^2
```

### mdanalysis_example.py

MDAnalysis-style analysis class:

```python
import MDAnalysis as mda
from zsasa.mdanalysis import SASAAnalysis

u = mda.Universe('topology.pdb', 'trajectory.xtc')
sasa = SASAAnalysis(u, select='protein')
sasa.run()

print(sasa.results.total_area)  # Per-frame total SASA in A^2
```

## Performance Tips

1. **Use multi-threading**: Set `n_threads=0` for auto-detection
2. **Use Lee-Richards for large systems**: `algorithm="lr"` is faster for >10k atoms
3. **Use f32 precision for speed**: `precision="f32"` is ~2x faster with minimal accuracy loss
4. **Batch trajectories**: zsasa automatically parallelizes across frames

## Common Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `algorithm` | `"sr"` | `"sr"` (Shrake-Rupley) or `"lr"` (Lee-Richards) |
| `n_points` | `100` | Test points per atom (SR algorithm) |
| `n_slices` | `20` | Slices per atom (LR algorithm) |
| `probe_radius` | `1.4` | Water probe radius in Angstroms |
| `n_threads` | `0` | Number of threads (0 = auto-detect) |
| `classifier` | `NACCESS` | Atom radius classifier |

## See Also

- [Python API Documentation](../../docs/python-api/)
- [CLI Documentation](../../docs/cli.md)
- [Algorithm Details](../../docs/algorithm.md)

---
sidebar_position: 1
sidebar_label: Choosing CLI vs Python
---

# Choosing CLI vs Python

zsasa provides two main interfaces. The right choice depends on your input data and workflow.

## Decision Guide

```
Is your structure pre-processed?
(single model, no ligands, no altloc, standard residues only)
│
├── Yes → CLI is the simplest option
│         zsasa structure.cif output.json
│
└── No  → Python integrations handle pre-processing
          calculate_sasa_from_structure("structure.cif")
```

## CLI: Pre-processed Structures

The CLI works best with **clean, simple structures** that need no pre-processing:

- AlphaFold/ESMFold predicted structures (single model, protein only)
- Pre-cleaned PDB/mmCIF files
- Batch processing of many files

```bash
# Single structure
zsasa structure.cif output.json

# With classifier
zsasa --classifier=naccess structure.cif output.json

# Batch: process all files in a directory
for f in structures/*.cif; do
    zsasa "$f" "results/$(basename "$f" .cif).json"
done
```

The CLI does **not** handle:
- Alternate conformations (altloc) — uses the first encountered
- HETATM filtering — includes all atoms
- Model selection for NMR ensembles — uses the first model
- Ligand removal or chain selection

## Python Integrations: Structure Processing Required

Use Python integrations when your structures need pre-processing or you're combining zsasa with other analysis tools.

### Structure pre-processing (Gemmi, BioPython, Biotite)

These integrations parse structure files through their respective libraries, which handle:

- **Altloc selection** — choose specific conformers
- **HETATM filtering** — include/exclude ligands, waters, ions
- **Model selection** — pick specific models from NMR ensembles
- **Chain/residue selection** — analyze specific parts of a structure

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure

# Basic usage (handles altloc, HETATM automatically)
result = calculate_sasa_from_structure("complex.cif")

# Exclude ligands and hydrogens
result = calculate_sasa_from_structure(
    "complex.cif",
    include_hetatm=False,
    include_hydrogens=False,
)
```

### MD trajectory analysis (MDAnalysis, MDTraj)

For molecular dynamics trajectories with atom selections:

```python
import MDAnalysis as mda
from zsasa.mdanalysis import SASAAnalysis

u = mda.Universe("topology.pdb", "trajectory.xtc")
sasa = SASAAnalysis(u, select="protein and not resname HOH")
sasa.run()
```

### Combining with other tools

When zsasa is part of a larger analysis pipeline:

```python
import gemmi
from zsasa.integrations.gemmi import calculate_sasa_from_model

# Load and manipulate structure with Gemmi
st = gemmi.read_structure("complex.cif")
st.remove_waters()
st.remove_ligands_and_waters()

# Then calculate SASA
result = calculate_sasa_from_model(st[0])
```

## Quick Reference

| Scenario | Recommended | Why |
|----------|------------|-----|
| AlphaFold structures | CLI | Pre-processed, no cleanup needed |
| Batch processing (clean files) | CLI | Fast, scriptable |
| PDB with ligands/waters | Python (Gemmi) | Filter HETATM, select chains |
| NMR ensemble (specific model) | Python (Gemmi/BioPython) | Model selection |
| MD trajectory | CLI `traj` or Python | CLI for quick analysis, Python for selections |
| Integration with BioPython pipeline | Python (BioPython) | Native object passing |
| Biotite/AtomWorks workflow | Python (Biotite) | Native AtomArray support |

## Available Integrations

| Integration | Library | Formats | Best For |
|-------------|---------|---------|----------|
| [Gemmi](../integrations/gemmi.md) | gemmi | mmCIF, PDB | Fast parsing, large files |
| [BioPython](../integrations/biopython.md) | BioPython | PDB, mmCIF | Existing BioPython workflows |
| [Biotite](../integrations/biotite.md) | Biotite | PDB, mmCIF, BinaryCIF | AtomWorks compatibility |
| [MDAnalysis](../integrations/mdanalysis.md) | MDAnalysis | XTC, TRR, DCD, etc. | MD trajectory + selections |
| [MDTraj](../integrations/mdtraj.md) | MDTraj | XTC, TRR, DCD, etc. | Drop-in for `mdtraj.shrake_rupley` |

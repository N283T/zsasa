---
sidebar_label: Overview
sidebar_position: 1
---

# Integration Modules

All integration modules provide a consistent interface for calculating SASA from structure files.

## Available Integrations

| Module | Library | Formats | Best For |
|--------|---------|---------|----------|
| [gemmi](gemmi.md) | gemmi | mmCIF, PDB | Fast parsing, large files |
| [biopython](biopython.md) | BioPython | PDB, mmCIF | Existing BioPython workflows |
| [biotite](biotite.md) | Biotite | PDB, mmCIF, BinaryCIF | AtomWorks compatibility |
| [mdtraj](mdtraj.md) | MDTraj | XTC, TRR, DCD, etc. | MD trajectory analysis |
| [mdanalysis](mdanalysis.md) | MDAnalysis | XTC, TRR, DCD, etc. | MD trajectory analysis |

## Common Functions

All structure file integrations (gemmi, biopython, biotite) provide:

| Function | Description |
|----------|-------------|
| `extract_atoms_from_model(model)` | Extract atom data from model |
| `calculate_sasa_from_model(model)` | Calculate SASA from model |
| `calculate_sasa_from_structure(source)` | Calculate SASA from file or structure |

Biotite uses `atom_array` instead of `model`:

| Function | Description |
|----------|-------------|
| `extract_atoms_from_atom_array(atom_array)` | Extract atom data from AtomArray |
| `calculate_sasa_from_atom_array(atom_array)` | Calculate SASA from AtomArray |

## Common Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `classifier` | `ClassifierType` | `NACCESS` | Atom radius classifier |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | SASA algorithm |
| `n_points` | `int` | `100` | Test points (SR) |
| `n_slices` | `int` | `20` | Slices (LR) |
| `probe_radius` | `float` | `1.4` | Probe radius in Å |
| `n_threads` | `int` | `0` | Threads (0 = auto) |
| `include_hetatm` | `bool` | `True` | Include HETATM records |
| `include_hydrogens` | `bool` | `False` | Include hydrogen atoms |
| `model_index` | `int` | `0` | Model index (NMR) |

## SasaResultWithAtoms

Extended result class with atom metadata:

```python
@dataclass
class SasaResultWithAtoms(SasaResult):
    atom_classes: NDArray[int32]  # Per-atom polarity classes
    atom_data: AtomData           # Atom metadata
    polar_area: float             # Total polar SASA
    apolar_area: float            # Total apolar SASA
```

## AtomData

```python
@dataclass
class AtomData:
    coords: NDArray[float64]   # (N, 3) coordinates
    residue_names: list[str]   # Residue names
    atom_names: list[str]      # Atom names
    chain_ids: list[str]       # Chain IDs
    residue_ids: list[int]     # Residue numbers
    elements: list[str]        # Element symbols
```

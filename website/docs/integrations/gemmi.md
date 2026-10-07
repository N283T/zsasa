# Gemmi Integration

Fast mmCIF and PDB parsing with gemmi.

## Installation

```bash
pip install zsasa[gemmi]
# or
uv add zsasa[gemmi]
```

## Import

```python
from zsasa.integrations.gemmi import (
    calculate_sasa_from_structure,
    calculate_sasa_from_model,
    extract_atoms_from_model,
)
```

## Supported Formats

- mmCIF (.cif, .cif.gz)
- PDB (.pdb, .pdb.gz)

## Example

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure

# From file
result = calculate_sasa_from_structure("protein.cif")

# From gemmi Structure
import gemmi
structure = gemmi.read_structure("protein.pdb")
result = calculate_sasa_from_structure(structure)

# Access results
print(f"Total: {result.total_area:.1f} Å²")
print(f"Polar: {result.polar_area:.1f} Å²")
print(f"Apolar: {result.apolar_area:.1f} Å²")
```

## Alternate Locations

A gemmi model holds every alternate conformer. The integration keeps one conformer per site by the rules of [`--altloc=auto`](../cli/input.md#alternate-locations) in the CLI (an atom without an altloc ID, then altloc `A`, then the highest occupancy, and one whole residue where the alternates of a position are different residues), so it counts the same atoms as `zsasa calc`. To use other conformers, edit the model first:

```python
import gemmi
from zsasa.integrations.gemmi import calculate_sasa_from_model

structure = gemmi.read_structure("protein.cif")
for chain in structure[0]:
    for residue in chain:
        # Keep the atoms without an altloc ID and those of conformer B
        for i in reversed(range(len(residue))):
            if residue[i].altloc not in ("\0", "B"):
                del residue[i]

result = calculate_sasa_from_model(structure[0])
```

## Why Gemmi?

- **Fast**: C++ backend, optimized for large files
- **Memory efficient**: Streaming parser
- **Comprehensive**: Full mmCIF dictionary support

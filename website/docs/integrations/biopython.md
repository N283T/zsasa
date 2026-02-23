# BioPython Integration

Support for BioPython Structure objects.

## Installation

```bash
pip install zsasa[biopython]
# or
uv add zsasa[biopython]
```

## Import

```python
from zsasa.integrations.biopython import (
    calculate_sasa_from_structure,
    calculate_sasa_from_model,
    extract_atoms_from_model,
)
```

## Supported Formats

- PDB (.pdb)
- mmCIF (.cif)

## Example

```python
from zsasa.integrations.biopython import calculate_sasa_from_structure

# From file
result = calculate_sasa_from_structure("protein.pdb")

# From BioPython Structure
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "protein.pdb")
result = calculate_sasa_from_structure(structure)

# Access results
print(f"Total: {result.total_area:.1f} Å²")
```

## When to Use

- Existing BioPython workflows
- Need BioPython's other features (sequence alignment, etc.)

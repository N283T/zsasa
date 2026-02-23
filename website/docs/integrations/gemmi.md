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

## Why Gemmi?

- **Fast**: C++ backend, optimized for large files
- **Memory efficient**: Streaming parser
- **Comprehensive**: Full mmCIF dictionary support

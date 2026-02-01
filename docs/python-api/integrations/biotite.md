# Biotite Integration

Support for Biotite AtomArray objects.

## Installation

```bash
pip install freesasa-zig[biotite]
```

## Import

```python
from freesasa_zig.integrations.biotite import (
    calculate_sasa_from_structure,
    calculate_sasa_from_atom_array,
    extract_atoms_from_atom_array,
)
```

## Supported Formats

- PDB (.pdb)
- mmCIF (.cif)
- BinaryCIF (.bcif)

## Example

```python
from freesasa_zig.integrations.biotite import calculate_sasa_from_structure

# From file
result = calculate_sasa_from_structure("protein.cif")

# From Biotite AtomArray
import biotite.structure.io as strucio
atom_array = strucio.load_structure("protein.pdb")
from freesasa_zig.integrations.biotite import calculate_sasa_from_atom_array
result = calculate_sasa_from_atom_array(atom_array)

# Access results
print(f"Total: {result.total_area:.1f} Å²")
```

## AtomWorks Compatibility

Also works with AtomWorks (built on Biotite):

```python
from atomworks.io.utils.io_utils import load_any
from freesasa_zig.integrations.biotite import calculate_sasa_from_atom_array

atom_array = load_any("protein.cif.gz")
result = calculate_sasa_from_atom_array(atom_array)
```

# Analysis API

Per-residue aggregation and analysis utilities.

## aggregate_by_residue

```python
def aggregate_by_residue(
    atom_areas: NDArray[float64],
    chain_ids: list[str],
    residue_ids: list[int],
    residue_names: list[str],
    atom_classes: NDArray[int32] | None = None,
    insertion_codes: list[str] | None = None,
) -> list[ResidueResult]
```

Aggregate per-atom SASA values to per-residue. Atoms are grouped by chain ID,
residue number and insertion code, and the residues are returned in the order
in which they first appear.

Pass `insertion_codes` (one entry per atom, `""` for a residue without an
insertion code) for structures that use them, such as antibodies: residues
`10`, `10A` and `10B` have the same `residue_ids` entry, and without their
insertion codes they are merged into one result. `aggregate_from_result` does
this for you.

The residue name is not part of the grouping key, and atoms of one residue need
not be contiguous. In both points this function differs from the CLI, whose
[residue identity](../cli/output.md#residue-identity) is a run of consecutive
atoms with the same chain ID, residue number, insertion code and residue name.

## aggregate_from_result

```python
def aggregate_from_result(result: SasaResultWithAtoms) -> list[ResidueResult]
```

Convenience wrapper for `SasaResultWithAtoms` from integration modules. It
passes the insertion codes that the gemmi, BioPython and Biotite integrations
record in `result.atom_data.insertion_codes`.

## ResidueResult

```python
@dataclass
class ResidueResult:
    chain_id: str       # Chain identifier
    residue_id: int     # Residue sequence number
    residue_name: str   # 3-letter residue code
    total_area: float   # Total SASA in Å²
    polar_area: float   # Polar SASA in Å²
    apolar_area: float  # Apolar SASA in Å²
    rsa: float | None   # Relative Solvent Accessibility
    n_atoms: int        # Number of atoms
    insertion_code: str = ""  # Insertion code, "" for a residue without one
```

`residue_id` and `insertion_code` together give the residue number as it is
written in the structure file: residue `10A` has `residue_id == 10` and
`insertion_code == "A"`.

---

## Examples

### Algorithm Comparison

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

# Shrake-Rupley (faster)
result_sr = calculate_sasa(coords, radii, algorithm="sr", n_points=100)
print(f"SR: {result_sr.total_area:.2f} Å²")

# Lee-Richards (more precise)
result_lr = calculate_sasa(coords, radii, algorithm="lr", n_slices=20)
print(f"LR: {result_lr.total_area:.2f} Å²")
```

### Multi-threading

```python
from zsasa import calculate_sasa

# Auto-detect CPU cores (default)
result = calculate_sasa(coords, radii, n_threads=0)

# Specific thread count
result = calculate_sasa(coords, radii, n_threads=4)

# Single-threaded
result = calculate_sasa(coords, radii, n_threads=1)
```

### Finding Buried Residues

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure
from zsasa import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

# Find buried residues (RSA < 20%)
buried = [r for r in residues if r.rsa is not None and r.rsa < 0.2]

print(f"Buried residues ({len(buried)}):")
for r in buried:
    print(f"  {r.chain_id}:{r.residue_name}{r.residue_id} - RSA: {r.rsa:.1%}")
```

### Custom Radii

```python
import numpy as np
from zsasa import calculate_sasa, classify_atoms

# Get radii from classifier
residues = ["ALA", "ALA", "ALA"]
atoms = ["N", "CA", "C"]
classification = classify_atoms(residues, atoms)

# Or use custom radii
custom_radii = np.array([1.65, 1.87, 1.76])

# Calculate with either
coords = np.array([[0, 0, 0], [1.5, 0, 0], [2.3, 1.2, 0]])
result = calculate_sasa(coords, classification.radii)
```

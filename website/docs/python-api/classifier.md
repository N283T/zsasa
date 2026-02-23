# Classifier API

Atom classification and RSA calculation.

## ClassifierType

```python
class ClassifierType(IntEnum):
    NACCESS = 0  # NACCESS-compatible radii (default)
    PROTOR = 1   # ProtOr radii
    OONS = 2     # OONS radii
```

## AtomClass

```python
class AtomClass(IntEnum):
    POLAR = 0    # Polar atoms (N, O, etc.)
    APOLAR = 1   # Apolar atoms (C, etc.)
    UNKNOWN = 2  # Unknown classification
```

## classify_atoms

```python
def classify_atoms(
    residues: list[str],
    atoms: list[str],
    classifier_type: ClassifierType = ClassifierType.NACCESS,
    *,
    include_classes: bool = True,
) -> ClassificationResult
```

Classify multiple atoms at once (batch operation).

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `residues` | `list[str]` | required | Residue names (e.g., ["ALA", "GLY"]) |
| `atoms` | `list[str]` | required | Atom names (e.g., ["CA", "N"]) |
| `classifier_type` | `ClassifierType` | `NACCESS` | Classifier to use |
| `include_classes` | `bool` | `True` | Whether to compute atom classes |

**Returns:** `ClassificationResult`

**Example:**

```python
from zsasa import classify_atoms, ClassifierType

result = classify_atoms(
    ["ALA", "ALA", "GLY"],
    ["CA", "O", "N"],
    ClassifierType.NACCESS
)
print(result.radii)   # [1.87, 1.4, 1.65]
print(result.classes) # [1, 0, 0] (APOLAR, POLAR, POLAR)
```

## ClassificationResult

```python
@dataclass
class ClassificationResult:
    radii: NDArray[float64]  # Per-atom radii (NaN for unknown)
    classes: NDArray[int32]  # Per-atom polarity classes
```

## Utility Functions

| Function | Description |
|----------|-------------|
| `get_radius(residue, atom, classifier)` | Get radius for a specific atom |
| `get_atom_class(residue, atom, classifier)` | Get polarity class for an atom |
| `guess_radius(element)` | Guess radius from element symbol |
| `guess_radius_from_atom_name(atom_name)` | Guess radius from PDB atom name |

---

# RSA API

Relative Solvent Accessibility calculation.

## MAX_SASA

Reference maximum SASA values from Tien et al. (2013):

```python
MAX_SASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 225.0, "GLU": 223.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
```

## calculate_rsa

```python
def calculate_rsa(sasa: float, residue_name: str) -> float | None
```

Calculate Relative Solvent Accessibility (RSA = SASA / MaxSASA).

**Returns:** RSA value (0.0-1.0+), or `None` for non-standard amino acids.

**Example:**

```python
from zsasa import calculate_rsa

rsa = calculate_rsa(64.5, "ALA")  # 64.5 / 129.0 = 0.5
```

## calculate_rsa_batch

```python
def calculate_rsa_batch(
    sasas: NDArray[float64] | list[float],
    residue_names: list[str],
) -> NDArray[float64]
```

Batch RSA calculation. Returns NaN for non-standard amino acids.

## get_max_sasa

```python
def get_max_sasa(residue_name: str) -> float | None
```

Get maximum SASA reference value for a residue.

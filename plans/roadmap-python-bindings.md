# Roadmap: Python Bindings Enhancement

## Goal

Python bindingsを強化し、主要な構造生物学ライブラリとの統合を提供する。

## Current State

```python
from freesasa_zig import calculate_sasa, classify_atoms

# 現在: numpy配列を直接渡す
coords = np.array([[x1,y1,z1], [x2,y2,z2], ...])
radii = np.array([r1, r2, ...])
result = calculate_sasa(coords, radii)
```

## Requirements

### Phase 1: Core API Enhancement

現在未実装または改善が必要な機能:

| Feature | Status | Priority |
|---------|--------|----------|
| `calculate_sasa()` | ✅ | - |
| `classify_atoms()` | ✅ | - |
| Per-residue aggregation | ❌ | High |
| RSA calculation | ❌ | High |
| Selection support | ❌ | Medium |
| Batch processing | ❌ | Medium |

### Phase 2: BioPython Integration

```python
from freesasa_zig.integrations import biopython as fsz_bp

# BioPython Structure → SASA
from Bio.PDB import PDBParser
parser = PDBParser()
structure = parser.get_structure("1crn", "1crn.pdb")

result = fsz_bp.calculate_sasa(structure)
# or
result = fsz_bp.from_structure(structure).calculate()
```

**Tasks:**
- [ ] `freesasa_zig/integrations/biopython.py`
- [ ] `from_structure()` / `from_model()` / `from_chain()`
- [ ] Atom coordinate extraction
- [ ] Residue/atom name mapping for classification
- [ ] Tests with various PDB structures

### Phase 3: Biotite Integration

```python
from freesasa_zig.integrations import biotite as fsz_bt
import biotite.structure.io as strucio

structure = strucio.load_structure("1crn.cif")
result = fsz_bt.calculate_sasa(structure)
```

**Tasks:**
- [ ] `freesasa_zig/integrations/biotite.py`
- [ ] AtomArray → coords/radii conversion
- [ ] Filter support (protein only, etc.)
- [ ] Tests

### Phase 4: Gemmi Integration

```python
from freesasa_zig.integrations import gemmi as fsz_gm
import gemmi

structure = gemmi.read_structure("1crn.cif")
result = fsz_gm.calculate_sasa(structure[0])  # Model
```

**Tasks:**
- [ ] `freesasa_zig/integrations/gemmi.py`
- [ ] Model/Chain/Residue level support
- [ ] Tests

### Phase 5: AtomWorks Integration

```python
from freesasa_zig.integrations import atomworks as fsz_aw

# AtomWorks DataFrame → SASA
result = fsz_aw.calculate_sasa(atom_df)
```

**Tasks:**
- [ ] `freesasa_zig/integrations/atomworks.py`
- [ ] DataFrame column mapping
- [ ] Tests

## Package Structure

```
freesasa_zig/
├── __init__.py           # Core API
├── core.py               # calculate_sasa, classify_atoms
├── analysis.py           # RSA, per-residue, etc.
└── integrations/
    ├── __init__.py
    ├── biopython.py      # Optional: Bio.PDB
    ├── biotite.py        # Optional: biotite
    ├── gemmi.py          # Optional: gemmi
    └── atomworks.py      # Optional: atomworks
```

## Optional Dependencies

```toml
# pyproject.toml
[project.optional-dependencies]
biopython = ["biopython>=1.80"]
biotite = ["biotite>=0.36"]
gemmi = ["gemmi>=0.6"]
atomworks = ["atomworks>=0.1"]
all = ["biopython", "biotite", "gemmi", "atomworks"]
```

```bash
pip install freesasa-zig[biopython]
pip install freesasa-zig[all]
```

## Implementation Order

1. [ ] Core API enhancement (RSA, per-residue)
2. [ ] BioPython integration (most widely used)
3. [ ] Biotite integration (modern, fast)
4. [ ] Gemmi integration (mmCIF専門)
5. [ ] AtomWorks integration

## Success Criteria

```python
# User can do this with any library they prefer
from Bio.PDB import PDBParser
from freesasa_zig.integrations.biopython import calculate_sasa

structure = PDBParser().get_structure("1crn", "1crn.pdb")
result = calculate_sasa(structure)

print(f"Total SASA: {result.total_area:.2f} Å²")
for res in result.residues:
    print(f"{res.name}{res.number}: {res.area:.2f} Å² (RSA: {res.rsa:.1%})")
```

---
- [ ] **DONE** - Phase complete

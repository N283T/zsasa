# Plan: Per-Residue Aggregation for Python Bindings

## Goal

Python bindingsにper-residue aggregation機能を追加。原子レベルのSASA結果を残基レベルに集計し、RSAも計算する。

## Current State

- `SasaResultWithAtoms` (gemmi.py): 原子ごとのSASA + メタデータ (chain, residue_id, residue_name)
- RSA functions (core.py): `calculate_rsa()`, `MAX_SASA` dict は既に実装済み
- 残基ごとの集計機能: 未実装

## Design

### New Module: `analysis.py`

```python
@dataclass
class ResidueResult:
    chain_id: str
    residue_id: int
    residue_name: str
    total_area: float
    polar_area: float
    apolar_area: float
    rsa: float | None  # None for non-standard amino acids
    n_atoms: int

def aggregate_by_residue(result: SasaResultWithAtoms) -> list[ResidueResult]:
    """Aggregate per-atom SASA to per-residue."""
    ...
```

### Why `analysis.py`?

- `core.py`: 低レベルのSASA計算
- `integrations/gemmi.py`: Gemmi固有の処理
- **`analysis.py`**: 高レベルの解析機能 (残基集計、RSAなど)

Gemmiに依存せず、メタデータがあれば任意のソースから使える設計。

## Files to Modify/Create

| File | Action | Description |
|------|--------|-------------|
| `python/freesasa_zig/analysis.py` | Create | ResidueResult, aggregate_by_residue() |
| `python/freesasa_zig/__init__.py` | Edit | Export analysis functions |
| `python/tests/test_analysis.py` | Create | Unit tests for aggregation |

## Implementation Steps

### 1. Create `analysis.py`

```python
# python/freesasa_zig/analysis.py
from dataclasses import dataclass
from collections import defaultdict
import numpy as np
from .core import AtomClass, MAX_SASA

@dataclass
class ResidueResult:
    """Per-residue SASA result."""
    chain_id: str
    residue_id: int
    residue_name: str
    total_area: float
    polar_area: float
    apolar_area: float
    rsa: float | None
    n_atoms: int

def aggregate_by_residue(
    atom_areas: np.ndarray,
    chain_ids: list[str],
    residue_ids: list[int],
    residue_names: list[str],
    atom_classes: np.ndarray | None = None,
) -> list[ResidueResult]:
    """Aggregate per-atom SASA values to per-residue."""
    ...
```

### 2. Helper for SasaResultWithAtoms

```python
def aggregate_from_result(result: "SasaResultWithAtoms") -> list[ResidueResult]:
    """Convenience wrapper for SasaResultWithAtoms."""
    return aggregate_by_residue(
        atom_areas=result.atom_areas,
        chain_ids=result.atom_data.chain_ids,
        residue_ids=result.atom_data.residue_ids,
        residue_names=result.atom_data.residue_names,
        atom_classes=result.atom_classes,
    )
```

### 3. Update `__init__.py`

```python
from .analysis import ResidueResult, aggregate_by_residue, aggregate_from_result
```

### 4. Tests (`test_analysis.py`)

Test cases:
- Single residue (all atoms in one group)
- Multiple residues (proper grouping)
- Multi-chain (same residue_id in different chains)
- Polar/apolar breakdown
- RSA calculation (standard vs non-standard amino acids)
- Empty input handling
- Real structure (1CRN from test fixtures)

## API Usage Example

```python
from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
from freesasa_zig.analysis import aggregate_from_result

result = calculate_sasa_from_structure("1crn.cif")
residues = aggregate_from_result(result)

for res in residues:
    rsa_str = f"{res.rsa:.1%}" if res.rsa else "N/A"
    print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: "
          f"{res.total_area:.1f} Å² (RSA: {rsa_str})")
```

## Verification

```bash
# Run tests
cd python && pixi run pytest tests/test_analysis.py -v

# Type check
cd python && pixi run ty check

# Lint
cd python && pixi run ruff check freesasa_zig/analysis.py
```

---
- [x] **DONE** - Phase complete (PR #93 merged)

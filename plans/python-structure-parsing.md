# Python Bindings: Structure Parsing Strategy

## Background

- CLI: mmCIF only (PDB非対応は設計思想)
- Python bindings: 柔軟性を持たせる余地あり

## Options

### Option A: Internal mmCIF Parser を C API で公開

**Pros:**
- 外部依存なし
- CLI と一貫性
- freesasa-zig 単体で完結

**Cons:**
- C API 追加の実装コスト
- mmCIF のみ (PDB 非対応)
- gemmi ほど機能豊富ではない

**実装内容:**
- `freesasa_parse_mmcif(path)` → atom配列
- Python側で座標・residue/atom名を受け取り `classify_atoms()` + `calculate_sasa()`

### Option B: gemmi に任せる (推奨?)

**Pros:**
- gemmi は高速 (C++)、機能豊富、広く使われている
- mmCIF + PDB 両対応
- Python bindings 側の実装が薄くなる
- ユーザーは gemmi の使い方を既に知っている可能性高い

**Cons:**
- 外部依存 (gemmi)
- freesasa-zig の Python package が gemmi を require する or optional dependency

**実装内容:**
- gemmi で構造を読む (ユーザー側 or ヘルパー関数)
- freesasa-zig は座標 + radii を受け取るだけ (現状のまま)
- Optional: `freesasa_zig.from_gemmi(structure)` ヘルパー

### Option C: Hybrid (A + B)

- Internal mmCIF parser を公開 (Option A)
- gemmi integration もオプションで提供 (Option B)

## Recommendation

**Option B (gemmi に任せる)** を推奨:

1. freesasa-zig の強みは **SASA計算の速度** であり、パーサーではない
2. gemmi は構造生物学の de facto standard になりつつある
3. 実装コストが低い (現状の `calculate_sasa()` で十分)
4. ユーザーは既に gemmi/BioPython で構造を扱っている

### 具体的な方針

```python
# ユーザーコード例
import gemmi
from freesasa_zig import calculate_sasa, classify_atoms
import numpy as np

# gemmi で読み込み
structure = gemmi.read_structure("protein.cif")
model = structure[0]

# 座標と atom 情報を抽出
coords = []
residues = []
atoms = []
for chain in model:
    for residue in chain:
        for atom in residue:
            coords.append(atom.pos.tolist())
            residues.append(residue.name)
            atoms.append(atom.name)

coords = np.array(coords)

# freesasa-zig で分類 + 計算
result = classify_atoms(residues, atoms)
sasa = calculate_sasa(coords, result.radii)
```

### Optional: ヘルパー関数の提供

```python
# freesasa_zig/integrations/gemmi.py (optional dependency)
def from_gemmi_model(model, classifier=ClassifierType.NACCESS):
    """Convert gemmi Model to freesasa-zig input."""
    ...
    return SasaInput(coords=coords, radii=radii, ...)
```

## Decision

**Option B 採用** (2025-01-27)

- gemmi を optional dependency として受け入れる
- Internal parser の C API 公開は不要
- BioPython integration は将来検討 (gemmi の方が高速)

## Next Steps

1. [ ] `freesasa_zig/integrations/gemmi.py` モジュール作成
2. [ ] `from_gemmi_model()` / `from_gemmi_structure()` ヘルパー実装
3. [ ] pyproject.toml に optional dependency 追加 (`pip install freesasa-zig[gemmi]`)
4. [ ] Tests 追加
5. [ ] Examples / ドキュメント追加

# Python API リファレンス

zsasa Python バインディングの包括的なドキュメントです。

## 概要

Python バインディングは以下を提供します:

- **NumPy統合**: 座標と半径をNumPy配列として渡す
- **2つのアルゴリズム**: Shrake-Rupley (高速) と Lee-Richards (高精度)
- **マルチスレッド**: CPUコア間での自動並列化
- **原子分類**: NACCESS、PROTOR、OONS分類器
- **RSA計算**: 相対溶媒接触性
- **残基ごとの集計**: 原子SASAを残基レベルに集計
- **ライブラリ連携**: gemmi、BioPython、Biotiteサポート

## インストール

### 必要条件

- Python 3.11+
- NumPy 1.24+
- Zig 0.15.2+ (ネイティブライブラリのビルド用)

### ソースから

```bash
# 1. Zig共有ライブラリをビルド
cd zsasa
zig build -Doptimize=ReleaseFast

# 2. Pythonパッケージをインストール
cd python
pip install -e .
```

### オプション依存関係

```bash
# 構造ファイルサポート用
pip install zsasa[gemmi]     # gemmi (高速mmCIF/PDB)
pip install zsasa[biopython] # BioPython
pip install zsasa[biotite]   # Biotite (AtomWorksでも動作)

# すべての連携
pip install zsasa[all]
```

### ライブラリの場所

バインディングは以下の場所でネイティブライブラリを検索します:

1. `ZSASA_LIB` 環境変数 (設定されている場合)
2. パッケージにバンドル (wheelインストール)
3. `../zig-out/lib/` (パッケージからの相対パス、開発用)
4. `/usr/local/lib/`
5. `/usr/lib/`
6. カレントディレクトリ
7. `./zig-out/lib/` (カレントディレクトリのzig-out)

---

## クイックスタート

### 基本的なSASA計算

```python
import numpy as np
from zsasa import calculate_sasa

# 原子座標 (N, 3) と半径 (N,) を定義
coords = np.array([
    [0.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
])
radii = np.array([1.5, 1.5])

# SASA計算
result = calculate_sasa(coords, radii)
print(f"合計: {result.total_area:.2f} Å²")
print(f"原子ごと: {result.atom_areas}")
```

### 構造ファイルから (推奨)

```python
# gemmi使用 (pip install zsasa[gemmi])
from zsasa.integrations.gemmi import calculate_sasa_from_structure

result = calculate_sasa_from_structure("protein.cif")
print(f"合計: {result.total_area:.1f} Å²")
print(f"極性: {result.polar_area:.1f} Å²")
print(f"非極性: {result.apolar_area:.1f} Å²")
```

### RSA付き残基ごとの解析

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure
from zsasa import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

for res in residues:
    rsa_str = f"{res.rsa:.1%}" if res.rsa is not None else "N/A"
    print(f"{res.chain_id}:{res.residue_name}{res.residue_id}: RSA={rsa_str}")
```

---

## コアAPI

### calculate_sasa

```python
def calculate_sasa(
    coords: NDArray[np.float64],
    radii: NDArray[np.float64],
    *,
    algorithm: Literal["sr", "lr"] = "sr",
    n_points: int = 100,
    n_slices: int = 20,
    probe_radius: float = 1.4,
    n_threads: int = 0,
) -> SasaResult
```

溶媒接触表面積を計算します。

**パラメータ:**

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `coords` | `NDArray[float64]` | 必須 | 原子座標 (N, 3) 配列 |
| `radii` | `NDArray[float64]` | 必須 | 原子半径 (N,) 配列 |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | Shrake-Rupley または Lee-Richards |
| `n_points` | `int` | `100` | 原子あたりのテスト点数 (SRのみ) |
| `n_slices` | `int` | `20` | 原子あたりのスライス数 (LRのみ) |
| `probe_radius` | `float` | `1.4` | 水プローブ半径 (Å) |
| `n_threads` | `int` | `0` | スレッド数 (0 = 自動) |

**戻り値:** `SasaResult`

**例外:** 無効な入力の場合 `ValueError`

### SasaResult

```python
@dataclass
class SasaResult:
    total_area: float           # 合計SASA (Å²)
    atom_areas: NDArray[float64] # 原子ごとのSASA値
```

---

## 分類器API

### ClassifierType

```python
class ClassifierType(IntEnum):
    NACCESS = 0  # NACCESS互換半径 (デフォルト)
    PROTOR = 1   # ProtOr半径
    OONS = 2     # OONS半径
```

### AtomClass

```python
class AtomClass(IntEnum):
    POLAR = 0    # 極性原子 (N, Oなど)
    APOLAR = 1   # 非極性原子 (Cなど)
    UNKNOWN = 2  # 不明な分類
```

### classify_atoms

```python
def classify_atoms(
    residues: list[str],
    atoms: list[str],
    classifier_type: ClassifierType = ClassifierType.NACCESS,
    *,
    include_classes: bool = True,
) -> ClassificationResult
```

複数の原子を一度に分類します (バッチ操作)。

**パラメータ:**

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `residues` | `list[str]` | 必須 | 残基名 (例: ["ALA", "GLY"]) |
| `atoms` | `list[str]` | 必須 | 原子名 (例: ["CA", "N"]) |
| `classifier_type` | `ClassifierType` | `NACCESS` | 使用する分類器 |
| `include_classes` | `bool` | `True` | 原子クラスを計算するか |

**戻り値:** `ClassificationResult`

**例:**

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

### ClassificationResult

```python
@dataclass
class ClassificationResult:
    radii: NDArray[float64]  # 原子ごとの半径 (不明な場合はNaN)
    classes: NDArray[int32]  # 原子ごとの極性クラス
```

### ユーティリティ関数

| 関数 | 説明 |
|----------|-------------|
| `get_radius(residue, atom, classifier)` | 特定の原子の半径を取得 |
| `get_atom_class(residue, atom, classifier)` | 原子の極性クラスを取得 |
| `guess_radius(element)` | 元素記号から半径を推測 |
| `guess_radius_from_atom_name(atom_name)` | PDB原子名から半径を推測 |

---

## RSA API

### MAX_SASA

Tien et al. (2013) からの最大SASA参照値:

```python
MAX_SASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 225.0, "GLU": 223.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
```

### calculate_rsa

```python
def calculate_rsa(sasa: float, residue_name: str) -> float | None
```

相対溶媒接触性 (RSA = SASA / MaxSASA) を計算します。

**戻り値:** RSA値 (0.0-1.0+)、または非標準アミノ酸の場合は `None`。

**例:**

```python
from zsasa import calculate_rsa

rsa = calculate_rsa(64.5, "ALA")  # 64.5 / 129.0 = 0.5
```

### calculate_rsa_batch

```python
def calculate_rsa_batch(
    sasas: NDArray[float64] | list[float],
    residue_names: list[str],
) -> NDArray[float64]
```

バッチRSA計算。非標準アミノ酸にはNaNを返します。

### get_max_sasa

```python
def get_max_sasa(residue_name: str) -> float | None
```

残基の最大SASA参照値を取得します。

---

## 解析API

### aggregate_by_residue

```python
def aggregate_by_residue(
    atom_areas: NDArray[float64],
    chain_ids: list[str],
    residue_ids: list[int],
    residue_names: list[str],
    atom_classes: NDArray[int32] | None = None,
) -> list[ResidueResult]
```

原子ごとのSASA値を残基ごとに集計します。

### aggregate_from_result

```python
def aggregate_from_result(result: SasaResultWithAtoms) -> list[ResidueResult]
```

連携モジュールの `SasaResultWithAtoms` 用の便利なラッパー。

### ResidueResult

```python
@dataclass
class ResidueResult:
    chain_id: str       # チェーン識別子
    residue_id: int     # 残基配列番号
    residue_name: str   # 3文字残基コード
    total_area: float   # 合計SASA (Å²)
    polar_area: float   # 極性SASA (Å²)
    apolar_area: float  # 非極性SASA (Å²)
    rsa: float | None   # 相対溶媒接触性
    n_atoms: int        # 原子数
```

---

## 連携モジュール

すべての連携モジュールは同じインターフェースを提供します:

### 共通関数

| 関数 | 説明 |
|----------|-------------|
| `extract_atoms_from_model(model)` | 原子データを抽出 (gemmi, BioPython) |
| `extract_atoms_from_atom_array(atom_array)` | 原子データを抽出 (Biotite) |
| `calculate_sasa_from_model(model)` | モデルからSASA計算 (gemmi, BioPython) |
| `calculate_sasa_from_atom_array(atom_array)` | AtomArrayからSASA計算 (Biotite) |
| `calculate_sasa_from_structure(source)` | ファイルまたは構造からSASA計算 (すべて) |

### 共通パラメータ

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `classifier` | `ClassifierType` | `NACCESS` | 原子半径分類器 |
| `algorithm` | `"sr"` or `"lr"` | `"sr"` | SASAアルゴリズム |
| `n_points` | `int` | `100` | テスト点数 (SR) |
| `n_slices` | `int` | `20` | スライス数 (LR) |
| `probe_radius` | `float` | `1.4` | プローブ半径 (Å) |
| `n_threads` | `int` | `0` | スレッド数 (0 = 自動) |
| `include_hetatm` | `bool` | `True` | HETATMレコードを含める |
| `include_hydrogens` | `bool` | `False` | 水素原子を含める |
| `model_index` | `int` | `0` | モデルインデックス (NMR) |

### SasaResultWithAtoms

```python
@dataclass
class SasaResultWithAtoms(SasaResult):
    atom_classes: NDArray[int32]  # 原子ごとの極性クラス
    atom_data: AtomData           # 原子メタデータ
    polar_area: float             # 合計極性SASA
    apolar_area: float            # 合計非極性SASA
```

### AtomData

```python
@dataclass
class AtomData:
    coords: NDArray[float64]   # (N, 3) 座標
    residue_names: list[str]   # 残基名
    atom_names: list[str]      # 原子名
    chain_ids: list[str]       # チェーンID
    residue_ids: list[int]     # 残基番号
    elements: list[str]        # 元素記号
```

---

## Gemmi連携

```python
from zsasa.integrations.gemmi import (
    calculate_sasa_from_structure,
    calculate_sasa_from_model,
    extract_atoms_from_model,
)
```

**必要:** `pip install zsasa[gemmi]`

**サポート形式:** mmCIF, PDB

**例:**

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure

# ファイルから
result = calculate_sasa_from_structure("protein.cif")

# gemmi Structureから
import gemmi
structure = gemmi.read_structure("protein.pdb")
result = calculate_sasa_from_structure(structure)
```

---

## BioPython連携

```python
from zsasa.integrations.biopython import (
    calculate_sasa_from_structure,
    calculate_sasa_from_model,
    extract_atoms_from_model,
)
```

**必要:** `pip install zsasa[biopython]`

**サポート形式:** PDB, mmCIF

**例:**

```python
from zsasa.integrations.biopython import calculate_sasa_from_structure

# ファイルから
result = calculate_sasa_from_structure("protein.pdb")

# BioPython Structureから
from Bio.PDB import PDBParser
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", "protein.pdb")
result = calculate_sasa_from_structure(structure)
```

---

## Biotite連携

```python
from zsasa.integrations.biotite import (
    calculate_sasa_from_structure,
    calculate_sasa_from_atom_array,
    extract_atoms_from_atom_array,
)
```

**必要:** `pip install zsasa[biotite]`

**サポート形式:** PDB, mmCIF, BinaryCIF

**AtomWorksでも動作** (Biotiteベース):

```python
from atomworks.io.utils.io_utils import load_any
from zsasa.integrations.biotite import calculate_sasa_from_atom_array

atom_array = load_any("protein.cif.gz")
result = calculate_sasa_from_atom_array(atom_array)
```

---

## 使用例

### アルゴリズム比較

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

# Shrake-Rupley (高速)
result_sr = calculate_sasa(coords, radii, algorithm="sr", n_points=100)
print(f"SR: {result_sr.total_area:.2f} Å²")

# Lee-Richards (高精度)
result_lr = calculate_sasa(coords, radii, algorithm="lr", n_slices=20)
print(f"LR: {result_lr.total_area:.2f} Å²")
```

### マルチスレッド

```python
from zsasa import calculate_sasa

# CPUコアを自動検出 (デフォルト)
result = calculate_sasa(coords, radii, n_threads=0)

# 特定のスレッド数
result = calculate_sasa(coords, radii, n_threads=4)

# シングルスレッド
result = calculate_sasa(coords, radii, n_threads=1)
```

### 埋没残基の検索

```python
from zsasa.integrations.gemmi import calculate_sasa_from_structure
from zsasa import aggregate_from_result

result = calculate_sasa_from_structure("protein.cif")
residues = aggregate_from_result(result)

# 埋没残基を検索 (RSA < 20%)
buried = [r for r in residues if r.rsa is not None and r.rsa < 0.2]

print(f"埋没残基 ({len(buried)}):")
for r in buried:
    print(f"  {r.chain_id}:{r.residue_name}{r.residue_id} - RSA: {r.rsa:.1%}")
```

### カスタム半径

```python
import numpy as np
from zsasa import calculate_sasa, classify_atoms

# 分類器から半径を取得
residues = ["ALA", "ALA", "ALA"]
atoms = ["N", "CA", "C"]
classification = classify_atoms(residues, atoms)

# またはカスタム半径を使用
custom_radii = np.array([1.65, 1.87, 1.76])

# どちらかで計算
coords = np.array([[0, 0, 0], [1.5, 0, 0], [2.3, 1.2, 0]])
result = calculate_sasa(coords, classification.radii)
```

---

## エラーハンドリング

### 一般的な例外

| 例外 | 原因 |
|-----------|-------|
| `ValueError` | 無効な入力 (形状の不一致、負の半径) |
| `FileNotFoundError` | 構造ファイルが見つからない |
| `IndexError` | 無効なモデルインデックス |
| `ImportError` | オプション依存関係が不足 (gemmi, biopython, biotite) |
| `MemoryError` | 計算中にメモリ不足 |
| `RuntimeError` | 計算エラー |

### 例

```python
from zsasa import calculate_sasa
import numpy as np

try:
    result = calculate_sasa(
        np.array([[0, 0, 0]]),
        np.array([-1.0])  # 無効: 負の半径
    )
except ValueError as e:
    print(f"エラー: {e}")
```

---

## パフォーマンス

FreeSASA Pythonバインディングとの比較 (シングルスレッド):

| 構造 | 原子数 | Zig SR | FS SR | 高速化 |
|-----------|------:|-------:|------:|--------:|
| 1CRN | 327 | 0.5ms | 0.7ms | 1.4x |
| 1UBQ | 602 | 0.6ms | 1.3ms | 2.2x |
| 1A0Q | 3,183 | 2.4ms | 7.6ms | 3.2x |
| 3HHB | 4,384 | 3.4ms | 11ms | 3.2x |
| 1AON | 58,674 | 44ms | 163ms | 3.7x |

- **SRアルゴリズム**: 1.4-3.7倍高速 (サイズが大きいほど高速化)
- **LRアルゴリズム**: 2.9-5.5倍高速
- **精度**: FreeSASAと一致 (0.01%未満の差)

---

## ユーティリティ関数

### get_version

```python
def get_version() -> str
```

ライブラリバージョン文字列を返します (例: "0.1.0")。

```python
from zsasa import get_version
print(get_version())  # "0.1.0"
```

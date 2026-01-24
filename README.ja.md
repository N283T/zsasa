# freesasa-zig

[English](README.md) | 日本語

[FreeSASA](https://github.com/mittinatten/freesasa)からZigに移植した溶媒接触可能表面積（SASA）計算ツール。

## 概要

SASA（Solvent Accessible Surface Area：溶媒接触可能表面積）は、生体分子の表面のうち溶媒分子がアクセス可能な領域の面積を測定する指標。本実装では2つのアルゴリズムを提供：

- **Shrake-Rupley (SR)**: テストポイント法（Golden Section Spiral）
- **Lee-Richards (LR)**: スライスベース法（円弧積分）

## 特徴

- **2つのSASAアルゴリズム**: Shrake-Rupley（高速）とLee-Richards（高精度）
- **構造ファイル直接入力**: PDB/mmCIF形式対応
- **チェーン/モデル選択**: チェーンID、モデル番号、authチェーンIDでフィルタリング
- **原子半径分類器**: NACCESS/ProtOr/OONS分類器（CLI・ライブラリ対応）
- **カスタム設定ファイル**: FreeSASA互換形式
- **解析機能**:
  - 残基単位SASA集計
  - RSA（相対溶媒接触可能性）計算
  - 極性/非極性表面分類
  - 複合体界面SASA
- JSON入出力（複数フォーマット対応）
- 各種パラメータ設定可能（テストポイント数、スライス数、プローブ半径）
- 詳細なエラーメッセージ付き入力バリデーション
- 明示的アロケータによるメモリ安全な実装
- **マルチスレッド対応** - 並列原子処理
- **SIMD最適化** - `@Vector(4, f64)`によるバッチ計算
- **近傍リスト最適化** - O(N²)からO(N)への計算量削減

## ビルド

Zig 0.15.2以上が必要。

**対応プラットフォーム**: Linux、macOS。Windowsユーザーは[WSL](https://learn.microsoft.com/ja-jp/windows/wsl/)（Windows Subsystem for Linux）でLinuxビルドを使用してください。

```bash
zig build
```

テスト実行:

```bash
zig build test
```

## 使い方

```bash
freesasa_zig [OPTIONS] <input> [output.json]
```

対応入力形式: JSON, PDB (.pdb), mmCIF (.cif, .cif.gz)

### 使用例

```bash
# 基本的な使用法 - Shrake-Rupley（デフォルト）
./zig-out/bin/freesasa_zig input.json output.json

# Lee-Richardsアルゴリズム
./zig-out/bin/freesasa_zig --algorithm=lr input.json output.json

# Lee-Richardsでスライス数を指定
./zig-out/bin/freesasa_zig --algorithm=lr --n-slices=50 input.json output.json

# マルチスレッド（両アルゴリズムで対応）
./zig-out/bin/freesasa_zig --threads=4 input.json output.json
./zig-out/bin/freesasa_zig --algorithm=lr --threads=4 input.json output.json

# Shrake-Rupleyパラメータのカスタマイズ
./zig-out/bin/freesasa_zig --probe-radius=1.5 --n-points=200 input.json output.json

# CSV出力フォーマット
./zig-out/bin/freesasa_zig --format=csv input.json output.csv

# コンパクトJSON（1行）
./zig-out/bin/freesasa_zig --format=compact input.json output.json

# 入力バリデーションのみ
./zig-out/bin/freesasa_zig --validate input.json

# 静寂モード
./zig-out/bin/freesasa_zig --quiet input.json output.json

# PDB/mmCIF直接入力（形式自動判別）
./zig-out/bin/freesasa_zig structure.pdb output.json
./zig-out/bin/freesasa_zig structure.cif output.json
./zig-out/bin/freesasa_zig structure.cif.gz output.json

# チェーン/モデル選択
./zig-out/bin/freesasa_zig --chain=A structure.cif output.json
./zig-out/bin/freesasa_zig --model=1 structure.pdb output.json
./zig-out/bin/freesasa_zig --auth-chain=A structure.cif output.json

# 残基単位解析
./zig-out/bin/freesasa_zig --per-residue structure.cif output.json

# RSA（相対溶媒接触可能性）
./zig-out/bin/freesasa_zig --rsa structure.cif output.json

# 極性/非極性表面解析
./zig-out/bin/freesasa_zig --polar structure.cif output.json

# 界面SASA（複合体）
./zig-out/bin/freesasa_zig --interface structure.cif output.json
```

### オプション

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `--algorithm=ALGO` | アルゴリズム: `sr` (Shrake-Rupley), `lr` (Lee-Richards) | sr |
| `--classifier=TYPE` | 組み込み分類器: `naccess`, `protor`, `oons` | - |
| `--config=FILE` | カスタム分類器設定ファイル（FreeSASA形式） | - |
| `--threads=N` | スレッド数 | 自動検出 |
| `--probe-radius=R` | プローブ半径（Å単位、0 < R ≤ 10） | 1.4 |
| `--n-points=N` | テストポイント数（SR専用、1-10000） | 100 |
| `--n-slices=N` | スライス数（LR専用、1-1000） | 20 |
| `--format=FORMAT` | 出力形式: `json`, `compact`, `csv` | json |
| `--chain=ID` | ラベルチェーンIDでフィルタ（例: "A", "B"） | - |
| `--model=N` | モデル番号選択（マルチモデルファイル用） | 1 |
| `--auth-chain=ID` | authチェーンIDでフィルタ（PDBチェーンID） | - |
| `--per-residue` | 残基単位SASA集計を出力 | - |
| `--rsa` | RSA計算（--per-residueを含む） | - |
| `--polar` | 極性/非極性サマリー表示（--per-residueを含む） | - |
| `--interface` | 複合体界面SASAを計算 | - |
| `--validate` | 入力バリデーションのみ実行 | - |
| `-q, --quiet` | 進捗出力を抑制 | - |
| `-h, --help` | ヘルプメッセージを表示 | - |
| `-V, --version` | バージョンを表示 | - |

### 入力形式

原子座標とvan der Waals半径を含むJSONファイル:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52]
}
```

**拡張形式**（分類器対応）:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["CA", "CB", "C"],
  "element": [6, 6, 6]
}
```

`residue`、`atom_name`、`element`フィールドはオプション:
- **`residue` + `atom_name`**: `--classifier`または`--config`オプション使用時に必要
- **`element`**: 元素の原子番号（6=C, 7=N, 8=O等）

**PDB/mmCIF入力**:

`.pdb`、`.cif`、`.cif.gz`拡張子のファイルは自動認識:

```bash
# 構造ファイル直接入力
./zig-out/bin/freesasa_zig 1CRN.pdb output.json
./zig-out/bin/freesasa_zig 1CRN.cif output.json
./zig-out/bin/freesasa_zig 1CRN.cif.gz output.json
```

構造ファイル使用時は分類器が自動適用（デフォルト: NACCESS）。

**バリデーションルール:**
- 全ての配列は同じ長さでなければならない
- 1原子以上のデータが必要
- 座標は有限値でなければならない（NaN/Inf不可）
- 半径は正の値かつ100 Å以下

### 出力形式

**JSON（デフォルト）** - インデント付き整形:
```json
{
  "total_area": 1234.56,
  "atom_areas": [
    10.5,
    20.3
  ]
}
```

**Compact** - 1行JSON:
```json
{"total_area":1234.56,"atom_areas":[10.5,20.3]}
```

**CSV** - カンマ区切り:
```csv
atom_index,area
0,10.500000
1,20.300000
total,1234.560000
```

## 入力データの準備

付属のPythonスクリプトで構造ファイルを変換:

```bash
# mmCIF/PDBを入力JSONに変換
./scripts/data/cif_to_json.py structure.cif output.json

# FreeSASAで参照SASA値を生成（検証用）
./scripts/data/calc_reference.py structure.cif reference.json

# ベンチマーク実行
./scripts/benchmark.py --runs 5 --threads 4
```

必要環境: Python 3.11+, gemmi, freesasa（PEP 723により自動インストール）

## アルゴリズム

2つのアルゴリズムを提供。詳細は[docs/algorithm.md](docs/algorithm.md)を参照。

### Shrake-Rupley (SR)

テストポイント法 - 高速でシンプル。

1. 単位球上にテストポイントを生成（Golden Section Spiral）
2. 各原子について、近傍原子に埋没していないポイントをカウント
3. SASA = 4π × r² × (露出数 / 総数)

### Lee-Richards (LR)

スライスベース法 - 数学的に厳密。

1. 各原子球をZ軸方向にスライス
2. 各スライスで露出円弧の長さを計算
3. 円弧長を積分して表面積を算出

### アルゴリズム比較

| 項目 | Shrake-Rupley | Lee-Richards |
|------|---------------|--------------|
| 手法 | テストポイント | スライス積分 |
| 精度制御 | `--n-points` | `--n-slices` |
| 速度（1A0Q, 4スレッド） | 13ms | 26ms |
| 推奨用途 | 大規模構造、迅速な解析 | 高精度が必要な場合 |

### パラメータ

| パラメータ | デフォルト | アルゴリズム | 説明 |
|-----------|-----------|-------------|------|
| `n_points` | 100 | SR | テストポイント数 |
| `n_slices` | 20 | LR | スライス数 |
| `probe_radius` | 1.4 Å | 両方 | 水分子の半径 |

## 検証

FreeSASA参照実装との比較（ProtOr分類器使用）:

| 構造 | 原子数 | FreeSASA (Å²) | Zig (Å²) | 差分 |
|------|--------|---------------|----------|------|
| 1CRN | 327 | 3,001.13 | 3,001.13 | 0.000% |
| 1UBQ | 602 | 4,834.72 | 4,834.72 | 0.000% |
| 1A0Q | 3,183 | 18,908.90 | 18,908.90 | 0.000% |
| 3HHB | 4,384 | 25,527.36 | 25,527.36 | 0.000% |
| 1AON | 58,674 | 316,879.14 | 316,879.14 | 0.000% |
| 4V6X | 237,685 | 1,325,369.25 | 1,325,369.25 | 0.000% |

検証実行: `./scripts/validate.py`

## 性能

Zig（ReleaseFast）とFreeSASA Pythonの比較、SASA計算時間のみ:

| 構造 | 原子数 | SR Zig (ms) | SR FS (ms) | SR 高速化 | LR Zig (ms) | LR FS (ms) | LR 高速化 |
|------|--------|-------------|------------|-----------|-------------|------------|-----------|
| 1CRN | 327 | 0.53 | 0.77 | 1.45x | 4.60 | 4.85 | 1.05x |
| 1UBQ | 602 | 0.64 | 1.40 | 2.19x | 8.45 | 9.19 | 1.09x |
| 1A0Q | 3,183 | 2.46 | 9.16 | 3.72x | 49.17 | 53.34 | 1.08x |
| 3HHB | 4,384 | 3.56 | 12.14 | 3.41x | 68.56 | 74.89 | 1.09x |
| 1AON | 58,674 | 43.76 | 179.52 | 4.10x | 930.22 | 1030.58 | 1.11x |
| 4V6X | 237,685 | 178.47 | 744.13 | 4.17x | 3740.60 | 4132.06 | 1.10x |

**概要**: Shrake-RupleyはFreeSASAより**1.5x-4.2x高速**。構造サイズが大きいほど高速化率が向上。

ベンチマーク実行: `./scripts/benchmark.py`

### 最適化技術

1. **近傍リスト**: O(N²)からO(N)への近傍探索
2. **SIMD**: `@Vector(4, f64)`によるバッチ計算
3. **マルチスレッド**: Work-stealingスレッドプールによる並列原子処理

## Pythonバインディング

NumPy配列を使ってPythonからfreesasa-zigを利用可能:

```python
import numpy as np
from freesasa_zig import calculate_sasa

# 原子座標 (N, 3) と半径 (N,)
coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])

# SASA計算
result = calculate_sasa(coords, radii)
print(f"総SASA: {result.total_area:.2f} Å²")

# Lee-Richardsアルゴリズムを使用
result_lr = calculate_sasa(coords, radii, algorithm="lr")
```

### インストール

```bash
# 共有ライブラリをビルド
zig build -Doptimize=ReleaseFast

# Pythonパッケージをインストール
cd python
pip install -e .
```

詳細なAPIドキュメントは[python/README.md](python/README.md)を参照。

## 原子分類器

分類器モジュールは残基名と原子名に基づいてvan der Waals半径と極性クラスを割り当てます。3種類の組み込み分類器と、カスタム設定ファイルに対応。

### 利用可能な分類器

| 分類器 | 説明 | C脂肪族 | ANYフォールバック |
|--------|------|---------|-------------------|
| **NACCESS** | デフォルト、NACCESS互換 | 1.87 Å | あり |
| **ProtOr** | ハイブリダイゼーションベース（Tsai et al. 1999） | 1.88 Å | なし |
| **OONS** | 旧FreeSASAデフォルト | 2.00 Å | あり |

### CLI使用方法

```bash
# 組み込み分類器を使用（入力JSONにresidue/atom_nameが必要）
./zig-out/bin/freesasa_zig --classifier=naccess input.json output.json
./zig-out/bin/freesasa_zig --classifier=protor input.json output.json
./zig-out/bin/freesasa_zig --classifier=oons input.json output.json

# カスタム設定ファイルを使用（FreeSASA形式）
./zig-out/bin/freesasa_zig --config=my_radii.config input.json output.json
```

分類器を指定すると:
1. 入力JSONの`residue`と`atom_name`フィールドで半径を検索
2. 見つからない場合は元素ベースの推定にフォールバック
3. `--config`と`--classifier`を両方指定した場合は`--config`が優先

### ライブラリ使用方法

```zig
const classifier = @import("classifier.zig");
const naccess = @import("classifier_naccess.zig");
const protor = @import("classifier_protor.zig");
const oons = @import("classifier_oons.zig");

// NACCESS分類器（デフォルト）
const radius = naccess.getRadius("ALA", "CA");  // 1.87

// ProtOr分類器（ハイブリダイゼーションベース）
const radius = protor.getRadius("ALA", "CA");   // 1.88

// OONS分類器（脂肪族炭素が大きい）
const radius = oons.getRadius("ALA", "CA");     // 2.00

// 元素ベースのフォールバック（未知原子用）
const r = classifier.guessRadiusFromAtomName(" CA "); // 1.70（炭素）
```

### ルックアップ順序

1. **残基固有**: 正確な（残基, 原子）マッチ
2. **ANYフォールバック**: 全残基共通の原子定義（NACCESS/OONSのみ）
3. **元素推定**: 元素記号からvan der Waals半径

詳細は[docs/classifier.md](docs/classifier.md)を参照。

## プロジェクト構成

```
freesasa-zig/
├── src/
│   ├── main.zig              # CLIエントリーポイント
│   ├── root.zig              # ライブラリルートモジュール
│   ├── types.zig             # データ構造（Vec3, AtomInput等）
│   ├── json_parser.zig       # JSON入力パース・バリデーション
│   ├── json_writer.zig       # 出力（JSON, CSV）
│   ├── pdb_parser.zig        # PDBファイルパーサー
│   ├── mmcif_parser.zig      # mmCIFファイルパーサー
│   ├── analysis.zig          # 解析機能（残基単位、RSA、極性、界面）
│   ├── classifier.zig        # 原子分類器コア（型、元素推定）
│   ├── classifier_naccess.zig # NACCESS組み込み分類器
│   ├── classifier_protor.zig  # ProtOr組み込み分類器
│   ├── classifier_oons.zig    # OONS組み込み分類器
│   ├── classifier_parser.zig  # FreeSASA設定ファイルパーサー
│   ├── test_points.zig       # Golden Section Spiral生成
│   ├── neighbor_list.zig     # 空間近傍リスト（O(N)探索）
│   ├── simd.zig              # SIMDバッチ処理
│   ├── thread_pool.zig       # 汎用スレッドプール
│   ├── shrake_rupley.zig     # Shrake-Rupleyアルゴリズム
│   └── lee_richards.zig      # Lee-Richardsアルゴリズム
├── scripts/
│   ├── benchmark.py               # 統合ベンチマーク
│   ├── validate.py                # FreeSASA参照値との検証
│   └── data/                      # データ準備スクリプト
│       ├── cif_to_json.py         # 構造→JSON変換
│       ├── calc_reference.py      # 参照SASA計算
│       ├── generate.py            # 構造DL・参照値生成
│       ├── generate_protor.py     # ProtOr半径入力生成
│       └── compare_classifiers.py # 分類器比較
├── benchmarks/
│   ├── structures/            # DLしたPDB構造 (.cif.gz)
│   ├── inputs/                # 生成した入力JSON（元素ベース半径）
│   ├── inputs_protor/         # 生成した入力JSON（ProtOr半径）
│   └── references/            # FreeSASA参照SASA値
├── examples/
│   ├── 1A0Q.cif.gz        # 元構造ファイル（PDB 1A0Q）
│   ├── input_1a0q.json    # 入力例（cifから変換）
│   └── 1A0Q_sasa.json     # FreeSASAによる参照SASA
├── python/
│   ├── freesasa_zig/      # Pythonバインディング
│   │   ├── __init__.py    # パッケージエクスポート
│   │   └── core.py        # ctypesバインディング
│   └── tests/             # Pythonテスト
├── docs/
│   ├── architecture.md    # アーキテクチャ詳解
│   ├── algorithm.md       # アルゴリズム詳解
│   ├── optimizations.md   # 最適化技術詳解
│   ├── cli-io.md          # CLI・入出力詳解
│   └── classifier.md      # 原子分類器詳解
└── plans/
    └── *.md               # 実装計画
```

## ロードマップ

- [x] Phase 1: 基本Shrake-Rupley実装
- [x] Phase 2: 近傍リスト最適化（O(N²) → O(N)）
- [x] Phase 3: SIMD最適化
- [x] Phase 4: マルチスレッド
- [x] Phase 5: 本番機能（CLI、出力形式、バリデーション）
- [x] Phase 6: CI/CDパイプライン
- [x] Phase 7: 公平なベンチマーク（タイミング分解、SASA計算のみ測定）
- [x] Phase 8: ベンチマークデータセット（6構造、小〜超大）
- [x] Phase 9: 半径分類器（原子半径の自動判定）
  - [x] コアデータ構造と元素ベース推定
  - [x] NACCESS/ProtOr/OONS組み込み分類器
  - [x] カスタム設定ファイルパーサー（FreeSASA形式）
  - [x] CLI統合（`--classifier`, `--config`）
- [x] Phase 10: 構造ファイル直接入力（PDB/mmCIF）
- [x] Phase 11: Lee-Richardsアルゴリズム（マルチスレッド・SIMD対応）
- [x] Phase 13: Pythonバインディング（C API経由のNumPy統合）
- [x] Phase 15: チェーン/モデル選択（`--chain`, `--model`, `--auth-chain`）
- [x] Phase 16: 解析機能
  - [x] 残基単位SASA集計（`--per-residue`）
  - [x] RSA計算（`--rsa`）
  - [x] 極性/非極性分類（`--polar`）
  - [x] 界面SASA（`--interface`）

## ライセンス

MIT

## 参考文献

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - オリジナルC実装

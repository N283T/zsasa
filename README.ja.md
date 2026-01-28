# freesasa-zig

[English](README.md) | 日本語

[FreeSASA](https://github.com/mittinatten/freesasa)からZigに移植した溶媒接触可能表面積（SASA）計算ツール。

## 概要

SASA（Solvent Accessible Surface Area：溶媒接触可能表面積）は、生体分子の表面のうち溶媒分子がアクセス可能な領域の面積を測定する指標。本実装では2つのアルゴリズムを提供：

- **Shrake-Rupley (SR)**: テストポイント法（Golden Section Spiral）
- **Lee-Richards (LR)**: スライスベース法（円弧積分）

## 特徴

- **2つのSASAアルゴリズム**: Shrake-Rupley（高速）とLee-Richards（高精度）
- **構造ファイル直接入力**: mmCIF形式対応
- **チェーン/モデル選択**: チェーンID、モデル番号、authチェーンIDでフィルタリング
- **原子半径分類器**: NACCESS/ProtOr/OONS分類器（CLI・ライブラリ対応）
- **カスタム設定ファイル**: FreeSASA互換形式
- **解析機能**:
  - 残基単位SASA集計
  - RSA（相対溶媒接触可能性）計算
  - 極性/非極性表面分類
- JSON入出力（複数フォーマット対応）
- 各種パラメータ設定可能（テストポイント数、スライス数、プローブ半径）
- 詳細なエラーメッセージ付き入力バリデーション
- 明示的アロケータによるメモリ安全な実装
- **マルチスレッド対応** - 並列原子処理
- **SIMD最適化** - `@Vector`による8並列バッチ計算（f32/f64選択可能）
- **高速三角関数** - LRアルゴリズム向け多項式近似
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

対応入力形式: JSON, mmCIF (.cif, .cif.gz)

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

# f32精度（高速、RustSASAとの公平比較向け）
./zig-out/bin/freesasa_zig --precision=f32 input.json output.json

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

# mmCIF直接入力（形式自動判別）
./zig-out/bin/freesasa_zig structure.cif output.json
./zig-out/bin/freesasa_zig structure.cif.gz output.json

# チェーン/モデル選択
./zig-out/bin/freesasa_zig --chain=A structure.cif output.json
./zig-out/bin/freesasa_zig --model=1 structure.cif output.json
./zig-out/bin/freesasa_zig --auth-chain --chain=A structure.cif output.json

# 残基単位解析
./zig-out/bin/freesasa_zig --per-residue structure.cif output.json

# RSA（相対溶媒接触可能性）
./zig-out/bin/freesasa_zig --rsa structure.cif output.json

# 極性/非極性表面解析
./zig-out/bin/freesasa_zig --polar structure.cif output.json

# バッチモード（ディレクトリ入力）
./zig-out/bin/freesasa_zig input_dir/ output_dir/ --threads=8
./zig-out/bin/freesasa_zig input_dir/ output_dir/ --timing  # ベンチマーク用
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
| `--precision=P` | 浮動小数点精度: `f32`, `f64` | f64 |
| `--format=FORMAT` | 出力形式: `json`, `compact`, `csv` | json |
| `--chain=ID` | ラベルチェーンIDでフィルタ（例: "A", "B"） | - |
| `--model=N` | モデル番号選択（マルチモデルファイル用） | 全て |
| `--auth-chain` | ラベルチェーンIDの代わりにauthチェーンIDを使用 | - |
| `--per-residue` | 残基単位SASA集計を出力 | - |
| `--rsa` | RSA計算（--per-residueを含む） | - |
| `--polar` | 極性/非極性サマリー表示（--per-residueを含む） | - |
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

**mmCIF入力**:

`.cif`、`.cif.gz`拡張子のファイルは自動認識:

```bash
# 構造ファイル直接入力
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
# CIFをベンチマークJSONに変換（ProtOr半径付き）
./benchmarks/scripts/generate_json.py --file structure.cif output.json.gz

# ベンチマーク実行
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10

# 結果解析
./benchmarks/scripts/analyze.py all
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
| 対FreeSASA C | 1.2x-2.3x高速 | 1.7x高速 |
| 推奨用途 | 大規模構造、迅速な解析 | 高精度が必要な場合 |

### パラメータ

| パラメータ | デフォルト | アルゴリズム | 説明 |
|-----------|-----------|-------------|------|
| `n_points` | 100 | SR | テストポイント数 |
| `n_slices` | 20 | LR | スライス数 |
| `probe_radius` | 1.4 Å | 両方 | 水分子の半径 |

## 検証

約10万構造でFreeSASAと検証（R² = 1.000、平均誤差 = 0.0004%）：

![Validation](benchmarks/results/plots/validation/sr.png)

## 性能

約10万PDB構造でベンチマーク。詳細は[測定手法](docs/benchmark/methodology.md)と[結果](docs/benchmark/results.md)を参照。

![構造サイズ・スレッド数別高速化率](benchmarks/results/plots/speedup_by_bin/grid.png)

**主な結果**（SASA計算時間のみ測定、単ファイル: Zig/FreeSASA f64、RustSASA f32）：
- **シングルスレッド**：大規模構造でFreeSASA Cより1.5x以上、RustSASAより1.1x高速
- **マルチスレッド**：効率大幅向上—500原子以上で1.0x超、20k原子以上で**1.9x-2.3x**高速
- **バッチ処理**：高精度（f64）でRustSASAと同等速度、f32モードならスループット優位

### 最適化技術

1. **近傍リスト**: O(N²)からO(N)への近傍探索
2. **8並列SIMD**: `@Vector`による8並列バッチ計算（f32/f64選択可能）
3. **高速三角関数**: LRアルゴリズム向けacos/atan2多項式近似（約37%高速化）
4. **マルチスレッド**: Work-stealingスレッドプールによる並列原子処理
5. **精度選択**: f32（高速）/ f64（高精度）切り替え

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
├── src/                   # Zigソースコード
│   ├── main.zig           # CLIエントリーポイント
│   ├── root.zig           # ライブラリルートモジュール
│   ├── types.zig          # データ構造（Vec3, AtomInput等）
│   ├── json_parser.zig    # JSON入力パース・バリデーション
│   ├── json_writer.zig    # 出力（JSON, CSV）
│   ├── mmcif_parser.zig   # mmCIFファイルパーサー
│   ├── analysis.zig       # 解析機能（残基単位、RSA、極性）
│   ├── classifier*.zig    # 原子分類器（NACCESS, ProtOr, OONS）
│   ├── shrake_rupley.zig  # Shrake-Rupleyアルゴリズム
│   └── lee_richards.zig   # Lee-Richardsアルゴリズム
├── benchmarks/
│   ├── scripts/           # ベンチマークツール
│   │   ├── run.py         # ベンチマーク実行（単一ツール/アルゴリズム）
│   │   ├── analyze.py     # 結果解析・グラフ生成
│   │   └── generate_json.py # CIF→JSON変換
│   ├── dataset/           # 標準6構造（JSON、git管理）
│   ├── cif/               # ソースCIFファイル（圧縮）
│   ├── external/          # 比較ツール（FreeSASA, RustSASAフォーク）
│   ├── inputs/            # 生成データ（gitignore）
│   └── results/           # ベンチマーク結果
├── python/
│   ├── freesasa_zig/      # Pythonバインディング
│   └── tests/             # Pythonテスト
├── docs/                  # 技術ドキュメント
└── plans/                 # 実装計画
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
- [x] Phase 10: 構造ファイル直接入力（mmCIF）
- [x] Phase 11: Lee-Richardsアルゴリズム（マルチスレッド・SIMD対応）
- [x] Phase 13: Pythonバインディング（C API経由のNumPy統合）
- [x] Phase 15: チェーン/モデル選択（`--chain`, `--model`, `--auth-chain`）
- [x] Phase 16: 解析機能
  - [x] 残基単位SASA集計（`--per-residue`）
  - [x] RSA計算（`--rsa`）
  - [x] 極性/非極性分類（`--polar`）

## ライセンス

MIT

## 参考文献

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - オリジナルC実装
- [RustSASA](https://github.com/maxall41/RustSASA) - Rust実装（バッチ処理最適化）

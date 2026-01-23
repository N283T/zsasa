# freesasa-zig

[English](README.md) | 日本語

[FreeSASA](https://github.com/mittinatten/freesasa)からZigに移植した溶媒接触可能表面積（SASA）計算ツール。

## 概要

SASA（Solvent Accessible Surface Area：溶媒接触可能表面積）は、生体分子の表面のうち溶媒分子がアクセス可能な領域の面積を測定する指標。本実装はShrake-Rupleyアルゴリズムを使用し、Golden Section Spiralによるテストポイント生成を行う。

## 特徴

- Shrake-Rupleyアルゴリズム実装
- JSON入出力（複数フォーマット対応）
- テストポイント数・プローブ半径の設定可能
- 詳細なエラーメッセージ付き入力バリデーション
- 明示的アロケータによるメモリ安全な実装
- **マルチスレッド対応** - 並列原子処理
- **SIMD最適化** - `@Vector(4, f64)`による4並列距離計算
- **近傍リスト最適化** - O(N²)からO(N)への計算量削減

## ビルド

Zig 0.15.2以上が必要。

```bash
zig build
```

テスト実行:

```bash
zig build test
```

## 使い方

```bash
freesasa_zig [OPTIONS] <input.json> [output.json]
```

### 使用例

```bash
# 基本的な使用法（スレッド数自動検出、整形済みJSON出力）
./zig-out/bin/freesasa_zig input.json output.json

# スレッド数指定
./zig-out/bin/freesasa_zig --threads=4 input.json output.json

# シングルスレッドモード
./zig-out/bin/freesasa_zig --threads=1 input.json output.json

# アルゴリズムパラメータのカスタマイズ
./zig-out/bin/freesasa_zig --probe-radius=1.5 --n-points=200 input.json output.json

# CSV出力フォーマット
./zig-out/bin/freesasa_zig --format=csv input.json output.csv

# コンパクトJSON（1行、空白なし）
./zig-out/bin/freesasa_zig --format=compact input.json output.json

# 入力バリデーションのみ（計算なし）
./zig-out/bin/freesasa_zig --validate input.json

# 静寂モード（進捗出力を抑制）
./zig-out/bin/freesasa_zig --quiet input.json output.json
```

### オプション

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `--threads=N` | スレッド数 | 自動検出 |
| `--probe-radius=R` | プローブ半径（Å単位、0 < R ≤ 10） | 1.4 |
| `--n-points=N` | 原子あたりのテストポイント数（1-10000） | 100 |
| `--format=FORMAT` | 出力形式: `json`, `compact`, `csv` | json |
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
./scripts/cif_to_input_json.py structure.cif output.json

# FreeSASAで参照SASA値を生成（検証用）
./scripts/calc_reference_sasa.py structure.cif reference.json

# ベンチマーク実行
./scripts/benchmark.py examples/1A0Q.cif.gz --runs 5 --threads 4
```

必要環境: Python 3.11+, gemmi, freesasa（PEP 723により自動インストール）

## アルゴリズム

### Shrake-Rupley法

1. Golden Section Spiralを用いて単位球上に均一分布のテストポイントを生成
2. 各原子について:
   - テストポイントを（van der Waals半径 + プローブ半径）でスケーリング
   - 原子位置に平行移動
   - 隣接原子内部に埋没していないポイントをカウント
   - SASA = 4π × r² × (露出ポイント数 / 総ポイント数)

### パラメータ

| パラメータ | デフォルト | 説明 |
|-----------|-----------|------|
| `n_points` | 100 | 原子あたりのテストポイント数 |
| `probe_radius` | 1.4 Å | 水分子の半径 |

## 検証

FreeSASA参照実装との比較:

| 構造 | 原子数 | 参照値 (Å²) | Zig (Å²) | 差分 |
|------|--------|------------|----------|------|
| 1A0Q | 3183 | 18923.28 | 19211.19 | 1.52% |

## 性能

PDB 1A0Q（3,183原子）でのベンチマーク:

| 構成 | 時間 | FreeSASA比 |
|------|------|------------|
| FreeSASA (Python) | 約51ms | 1.0x |
| Zig (シングルスレッド) | 約13ms | 3.9x |
| Zig (マルチスレッド) | 約8ms | **6.4x** |

### 最適化技術

1. **近傍リスト**: O(N²)からO(N)への近傍探索
2. **SIMD**: `@Vector(4, f64)`による4並列距離計算
3. **マルチスレッド**: Work-stealingスレッドプールによる並列原子処理

## プロジェクト構成

```
freesasa-zig/
├── src/
│   ├── main.zig           # CLIエントリーポイント
│   ├── root.zig           # ライブラリルートモジュール
│   ├── types.zig          # データ構造（Vec3, AtomInput等）
│   ├── json_parser.zig    # JSON入力パース・バリデーション
│   ├── json_writer.zig    # 出力（JSON, CSV）
│   ├── test_points.zig    # Golden Section Spiral生成
│   ├── neighbor_list.zig  # 空間近傍リスト（O(N)探索）
│   ├── simd.zig           # SIMDバッチ処理
│   ├── thread_pool.zig    # 汎用スレッドプール
│   └── shrake_rupley.zig  # コアSASAアルゴリズム
├── scripts/
│   ├── cif_to_input_json.py   # 構造→JSON変換
│   ├── calc_reference_sasa.py # 参照SASA計算
│   └── benchmark.py           # 性能ベンチマーク
├── examples/
│   ├── 1A0Q.cif.gz        # 元構造ファイル（PDB 1A0Q）
│   ├── input_1a0q.json    # 入力例（cifから変換）
│   └── 1A0Q_sasa.json     # FreeSASAによる参照SASA
├── docs/
│   ├── architecture.md    # アーキテクチャ詳解
│   ├── algorithm.md       # アルゴリズム詳解
│   ├── optimizations.md   # 最適化技術詳解
│   └── cli-io.md          # CLI・入出力詳解
└── plans/
    └── *.md               # 実装計画
```

## ロードマップ

- [x] Phase 1: 基本Shrake-Rupley実装
- [x] Phase 2: 近傍リスト最適化（O(N²) → O(N)）
- [x] Phase 3: SIMD最適化（FreeSASA比4.5倍高速化）
- [x] Phase 4: マルチスレッド（FreeSASA比6.4倍高速化）
- [x] Phase 5: 本番機能（CLIオプション、出力形式、バリデーション）

## ライセンス

MIT

## 参考文献

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- [FreeSASA](https://github.com/mittinatten/freesasa) - オリジナルC実装

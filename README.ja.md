# freesasa-zig

[English](README.md) | 日本語

[FreeSASA](https://github.com/mittinatten/freesasa)からZigに移植した溶媒接触可能表面積（SASA）計算ツール。

## 概要

SASA（Solvent Accessible Surface Area：溶媒接触可能表面積）は、生体分子の表面のうち溶媒分子がアクセス可能な領域の面積を測定する指標。本実装では2つのアルゴリズムを提供：

- **Shrake-Rupley (SR)**: テストポイント法（Golden Section Spiral）
- **Lee-Richards (LR)**: スライスベース法（円弧積分）

## 特徴

- **2つのSASAアルゴリズム**: Shrake-Rupley（高速）とLee-Richards（高精度）
- JSON入出力（複数フォーマット対応）
- 各種パラメータ設定可能（テストポイント数、スライス数、プローブ半径）
- 詳細なエラーメッセージ付き入力バリデーション
- 明示的アロケータによるメモリ安全な実装
- **マルチスレッド対応** - 並列原子処理
- **SIMD最適化** - `@Vector(4, f64)`によるバッチ計算
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
```

### オプション

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `--algorithm=ALGO` | アルゴリズム: `sr` (Shrake-Rupley), `lr` (Lee-Richards) | sr |
| `--threads=N` | スレッド数 | 自動検出 |
| `--probe-radius=R` | プローブ半径（Å単位、0 < R ≤ 10） | 1.4 |
| `--n-points=N` | テストポイント数（SR専用、1-10000） | 100 |
| `--n-slices=N` | スライス数（LR専用、1-1000） | 20 |
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

FreeSASA参照実装との比較:

| 構造 | 原子数 | 参照値 (Å²) | Zig (Å²) | 差分 |
|------|--------|------------|----------|------|
| 1A0Q | 3183 | 18923.28 | 19211.19 | 1.52% |

## 性能

PDB 1A0Q（3,183原子）でのベンチマーク、ReleaseFastビルド:

| アルゴリズム | 1スレッド | 4スレッド | FreeSASA比 |
|-------------|-----------|-----------|------------|
| Shrake-Rupley | 17ms | 13ms | **3.9x高速** |
| Lee-Richards | 53ms | 26ms | 2.0x高速 |
| FreeSASA (Python) | 約51ms | - | 1.0x |

### 最適化技術

1. **近傍リスト**: O(N²)からO(N)への近傍探索
2. **SIMD**: `@Vector(4, f64)`によるバッチ計算
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
│   ├── shrake_rupley.zig  # Shrake-Rupleyアルゴリズム
│   └── lee_richards.zig   # Lee-Richardsアルゴリズム
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
- [x] Phase 3: SIMD最適化
- [x] Phase 4: マルチスレッド
- [x] Phase 5: 本番機能（CLI、出力形式、バリデーション）
- [x] Phase 6: CI/CDパイプライン
- [x] Phase 11: Lee-Richardsアルゴリズム（マルチスレッド・SIMD対応）
- [ ] Phase 9: 半径分類器（原子半径の自動判定）
- [ ] Phase 10: mmCIF直接入力対応

## ライセンス

MIT

## 参考文献

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.
- [FreeSASA](https://github.com/mittinatten/freesasa) - オリジナルC実装

# CLI リファレンス

zsasa コマンドラインインターフェースの完全なドキュメントです。

## 概要

```
zsasa [OPTIONS] <input> [output.json]
zsasa [OPTIONS] <input_dir/> <output_dir/>
```

## 基本的な使い方

```bash
# 基本的なSASA計算
./zig-out/bin/zsasa structure.cif output.json

# アルゴリズム選択
./zig-out/bin/zsasa --algorithm=lr structure.cif output.json

# マルチスレッド
./zig-out/bin/zsasa --threads=4 structure.cif output.json

# 解析機能付き
./zig-out/bin/zsasa --rsa --polar structure.cif output.json
```

## オプション

### アルゴリズムオプション

| オプション | 説明 | デフォルト |
|--------|-------------|---------|
| `--algorithm=ALGO` | `sr` (Shrake-Rupley) または `lr` (Lee-Richards) | `sr` |
| `--precision=P` | 浮動小数点精度: `f32` または `f64` | `f64` |
| `--probe-radius=R` | プローブ半径 (Å単位, 0 < R ≤ 10) | `1.4` |
| `--n-points=N` | 原子あたりのテスト点数 (SRのみ, 1-10000) | `100` |
| `--n-slices=N` | 原子直径あたりのスライス数 (LRのみ, 1-1000) | `20` |
| `--threads=N` | スレッド数 (0 = 自動検出) | `0` |

### 分類器オプション

| オプション | 説明 |
|--------|-------------|
| `--classifier=TYPE` | 組み込み分類器: `naccess`, `protor`, `oons` |
| `--config=FILE` | カスタム分類器設定ファイル (FreeSASA形式) |

`--classifier` を使用すると、残基名と原子名に基づいて原子半径が割り当てられます。`--classifier` と `--config` の両方が指定された場合は `--config` が優先されます。

### 構造フィルタリング

| オプション | 説明 | デフォルト |
|--------|-------------|---------|
| `--chain=ID` | チェーンIDでフィルタ (例: `A` または `A,B,C`) | 全チェーン |
| `--model=N` | NMR構造のモデル番号 (≥1) | 全モデル |
| `--auth-chain` | label_asym_idの代わりにauth_asym_idを使用 | label_asym_id |

### 解析オプション

| オプション | 説明 |
|--------|-------------|
| `--per-residue` | 残基ごとのSASA集計を表示 |
| `--rsa` | 相対溶媒接触性 (RSA) を計算 (`--per-residue`を有効化) |
| `--polar` | 極性/非極性SASAサマリーを表示 (`--per-residue`を有効化) |

### 出力オプション

| オプション | 説明 | デフォルト |
|--------|-------------|---------|
| `-o, --output=FILE` | 出力ファイルパス | `output.json` |
| `--format=FMT` | 出力形式: `json`, `compact`, `csv` | `json` |
| `--timing` | タイミング内訳を表示 (ベンチマーク用) | オフ |
| `-q, --quiet` | 進捗出力を抑制 | オフ |
| `--validate` | 入力検証のみ、計算しない | オフ |

### 情報オプション

| オプション | 説明 |
|--------|-------------|
| `-h, --help` | ヘルプメッセージを表示 |
| `-V, --version` | バージョンを表示 |

### バッチモードオプション

| オプション | 説明 | デフォルト |
|--------|-------------|---------|
| `--parallelism=MODE` | 並列化戦略: `file`, `atom`, `pipeline` | `file` |

---

## 入力形式

入力形式はファイル拡張子から自動検出されます:

| 拡張子 | 形式 |
|-----------|--------|
| `.json`, `.json.gz`, `.json.zst` | JSON |
| `.cif`, `.mmcif`, `.CIF`, `.mmCIF` | mmCIF |
| `.pdb`, `.PDB`, `.ent`, `.ENT` | PDB |

### JSON形式

座標と半径のみの最小限のJSON入力:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [4.0, 5.0, 6.0],
  "z": [7.0, 8.0, 9.0],
  "r": [1.7, 1.55, 1.52]
}
```

分類情報を含む拡張JSON (`--classifier`使用時に必要):

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [4.0, 5.0, 6.0],
  "z": [7.0, 8.0, 9.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["N", "CA", "C"],
  "element": [7, 6, 6]
}
```

**フィールド:**

| フィールド | 必須 | 説明 |
|-------|----------|-------------|
| `x`, `y`, `z` | 必須 | 原子座標 (Å単位) |
| `r` | 必須 | ファンデルワールス半径 (Å単位) |
| `residue` | 分類器使用時 | 3文字残基コード (例: "ALA") |
| `atom_name` | 分類器使用時 | 原子名 (例: "CA", "N") |
| `element` | オプション | 原子番号 (例: 6=C, 7=N, 8=O) |

**検証ルール:**

- すべての配列は同じ長さでなければならない
- 配列は空であってはならない
- 座標は有限値でなければならない (NaNやInfは不可)
- 半径は正で100 Å以下でなければならない

### mmCIF形式

標準のmmCIFファイルがサポートされています。パーサーは以下を抽出します:

- `_atom_site.Cartn_x/y/z` - 座標
- `_atom_site.type_symbol` - 元素 (VdW半径用)
- `_atom_site.label_atom_id` / `auth_atom_id` - 原子名
- `_atom_site.label_comp_id` / `auth_comp_id` - 残基名
- `_atom_site.label_asym_id` / `auth_asym_id` - チェーンID
- `_atom_site.label_seq_id` / `auth_seq_id` - 残基番号
- `_atom_site.pdbx_PDB_ins_code` - 挿入コード
- `_atom_site.pdbx_PDB_model_num` - モデル番号
- `_atom_site.label_alt_id` - 代替位置 (デフォルトで最初のものを保持)

### PDB形式

ATOMおよびHETATMレコードを含む標準PDB形式ファイルがサポートされています。

---

## 出力形式

### JSON (デフォルト)

2スペースインデントの整形済みJSON:

```json
{
  "total_area": 18923.28,
  "atom_areas": [
    32.47,
    0.25,
    15.82
  ]
}
```

### コンパクトJSON

空白なしの1行JSON:

```json
{"total_area":18923.28,"atom_areas":[32.47,0.25,15.82]}
```

### CSV

原子インデックスと面積を含む基本CSV:

```csv
atom_index,area
0,32.470000
1,0.250000
2,15.820000
total,18923.280000
```

入力に構造情報がある場合 (mmCIF/PDB)、詳細なCSVが生成されます:

```csv
chain,residue,resnum,atom_name,x,y,z,radius,area
A,ALA,1,N,1.000,3.000,5.000,1.500,32.470000
A,ALA,1,CA,2.000,4.000,6.000,1.700,0.250000
,,,,,,,,18923.280000
```

---

## 解析機能

### 残基ごとの集計 (`--per-residue`)

原子SASAを残基 (チェーン + 残基番号 + 挿入コード) ごとにグループ化:

```
Per-residue SASA:
Chain  Res    Num       SASA  Atoms
----- ---- ------ ---------- ------
    A  MET      1     198.52     19
    A  LYS      2     142.31     22
    A  ALA      3      45.67      5
```

### RSA計算 (`--rsa`)

相対溶媒接触性 (RSA = SASA / MaxSASA) を計算します。

MaxSASA参照値は Tien et al. (2013) より:

| 残基 | MaxSASA (Å²) | 残基 | MaxSASA (Å²) |
|---------|-------------|---------|-------------|
| ALA | 129.0 | LEU | 201.0 |
| ARG | 274.0 | LYS | 236.0 |
| ASN | 195.0 | MET | 224.0 |
| ASP | 193.0 | PHE | 240.0 |
| CYS | 167.0 | PRO | 159.0 |
| GLN | 225.0 | SER | 155.0 |
| GLU | 223.0 | THR | 172.0 |
| GLY | 104.0 | TRP | 285.0 |
| HIS | 224.0 | TYR | 263.0 |
| ILE | 197.0 | VAL | 174.0 |

RSA値付きの出力 (露出した末端残基では1.0を超える場合があります):

```
Per-residue SASA with RSA:
Chain  Res    Num       SASA    RSA  Atoms
----- ---- ------ ---------- ------ ------
    A  MET      1     198.52   0.89     19
    A  LYS      2     142.31   0.60     22
    A  ALA      3      45.67   0.35      5
```

### 極性/非極性サマリー (`--polar`)

残基を分類し、SASAの内訳を表示します。自動的に `--per-residue` を有効にします。

- **極性**: ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR, TYR
- **非極性**: ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, TRP, VAL
- **不明**: 非標準残基 (リガンド、修飾残基など)

```
Polar/Nonpolar SASA:
  Polar:       2345.67 Å² ( 45.2%) - 42 residues
  Nonpolar:    2845.23 Å² ( 54.8%) - 58 residues
  Unknown:        0.00 Å² (  0.0%) -  0 residues
```

---

## バッチモード

ディレクトリ内のすべての構造ファイルを処理:

```bash
# 基本的なバッチ処理
./zig-out/bin/zsasa input_dir/ output_dir/

# ファイルレベル並列処理 (N個のファイルを並列、各1スレッド)
./zig-out/bin/zsasa --threads=8 input_dir/ output_dir/

# 原子レベル並列処理 (1ファイルずつ、SASAにNスレッド)
./zig-out/bin/zsasa --parallelism=atom --threads=8 input_dir/

# パイプラインモード (I/Oプリフェッチ + 原子レベルSASA)
./zig-out/bin/zsasa --parallelism=pipeline --threads=8 input_dir/
```

**並列化戦略:**

| 戦略 | 説明 | 最適な用途 |
|----------|-------------|----------|
| `file` | N個のファイルを並列、ファイルごとに1スレッド | 多数の小さいファイル |
| `atom` | 1ファイルずつ、SASAにNスレッド | 少数の大きいファイル |
| `pipeline` | I/Oプリフェッチ + 原子レベルSASA | バランスの取れたワークロード |

---

## 分類器

組み込み分類器は残基名と原子名に基づいて原子半径を割り当てます。

### NACCESS (`--classifier=naccess`)

NACCESS互換の半径。多くのツールで使用されるデフォルト分類器です。

### ProtOr (`--classifier=protor`)

Tsai et al. (1999) によるProtOr半径。

### OONS (`--classifier=oons`)

Ooi et al. によるOONS半径。

### カスタム設定 (`--config=FILE`)

FreeSASA設定形式のカスタム分類器を読み込み:

```
name custom_classifier
atom ALA CA 1.87
atom ALA CB 1.87
atom ALA C  1.76
atom ALA N  1.65
atom ALA O  1.40
# ... 他の原子定義
```

---

## 使用例

### 基本的な計算

```bash
# mmCIF入力
./zig-out/bin/zsasa structure.cif output.json

# PDB入力
./zig-out/bin/zsasa structure.pdb output.json

# JSON入力
./zig-out/bin/zsasa atoms.json output.json
```

### アルゴリズム選択

```bash
# Lee-Richards (50スライス)
./zig-out/bin/zsasa --algorithm=lr --n-slices=50 structure.cif output.json

# Shrake-Rupley (200テスト点)
./zig-out/bin/zsasa --algorithm=sr --n-points=200 structure.cif output.json
```

### パフォーマンスチューニング

```bash
# 高速モード: f32精度
./zig-out/bin/zsasa --precision=f32 structure.cif output.json

# マルチスレッド
./zig-out/bin/zsasa --threads=4 structure.cif output.json

# タイミング内訳を表示
./zig-out/bin/zsasa --timing structure.cif output.json
```

### 分類器の使用

```bash
# NACCESS分類器
./zig-out/bin/zsasa --classifier=naccess structure.cif output.json

# カスタム設定
./zig-out/bin/zsasa --config=custom.config structure.cif output.json
```

### チェーン/モデルフィルタリング

```bash
# 単一チェーン
./zig-out/bin/zsasa --chain=A structure.cif output.json

# 複数チェーン
./zig-out/bin/zsasa --chain=A,B,C structure.cif output.json

# 特定のモデル (NMR)
./zig-out/bin/zsasa --model=1 nmr_structure.cif output.json

# authチェーンIDを使用
./zig-out/bin/zsasa --auth-chain --chain=A structure.cif output.json
```

### 解析機能

```bash
# 残基ごとの出力
./zig-out/bin/zsasa --per-residue structure.cif output.json

# RSA計算
./zig-out/bin/zsasa --rsa structure.cif output.json

# 極性/非極性サマリー
./zig-out/bin/zsasa --polar structure.cif output.json

# 組み合わせ解析
./zig-out/bin/zsasa --rsa --polar structure.cif output.json
```

### 出力形式

```bash
# CSV出力
./zig-out/bin/zsasa --format=csv structure.cif output.csv

# コンパクトJSON
./zig-out/bin/zsasa --format=compact structure.cif output.json
```

### 検証のみ

```bash
# 計算なしで検証
./zig-out/bin/zsasa --validate structure.cif
```

---

## エラーメッセージ

| エラー | 説明 |
|-------|-------------|
| `Missing input file` | 入力ファイルが指定されていない |
| `Cannot access '<path>'` | 入力ファイルが見つからないか読み取れない |
| `Invalid probe radius` | プローブ半径が範囲外 (0, 10] |
| `Invalid n-points` | テスト点数が範囲外 [1, 10000] |
| `Invalid n-slices` | スライス数が範囲外 [1, 1000] |
| `Invalid format` | 不明な出力形式 |
| `Invalid algorithm` | 不明なアルゴリズム名 |
| `Invalid classifier` | 不明な分類器名 |
| `Invalid model number` | モデル番号は1以上でなければならない |
| `Array lengths do not match` | JSON配列の長さが一致しない |
| `Radius must be positive` | 入力の半径が0以下 |
| `Coordinate is not finite` | 座標にNaNまたはInfが含まれている |
| `Classifier requires 'residue' and 'atom_name'` | 分類情報が不足 |

---

## 終了コード

| コード | 意味 |
|------|---------|
| 0 | 成功 |
| 1 | エラー (無効な入力、ファイルが見つからないなど) |

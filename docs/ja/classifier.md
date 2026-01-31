# 原子分類器

残基名と原子名に基づいてファンデルワールス半径と極性クラスを割り当てるモジュールです。

## 概要

SASA計算には各原子のファンデルワールス半径が必要です。分類器はFreeSASAの設定ファイル形式に従い、以下の順序で半径を決定します:

1. **残基特異的マッチ**: (残基名, 原子名) の完全一致
2. **ANYフォールバック**: すべての残基に共通の原子定義
3. **元素推定**: 原子名から元素記号を抽出してファンデルワールス半径を推定

## モジュール構成

| ファイル | 説明 |
|------|-------------|
| `classifier.zig` | コアデータ構造、元素ベースの推定、ClassifierType |
| `classifier_naccess.zig` | NACCESS互換の組み込み分類器 |
| `classifier_protor.zig` | ProtOr分類器（混成軌道ベース） |
| `classifier_oons.zig` | OONS分類器（旧FreeSASAデフォルト） |
| `classifier_parser.zig` | FreeSASA互換の設定ファイルパーサー |

## CLI使用方法

### 組み込み分類器

```bash
# NACCESS分類器（推奨デフォルト）
./zig-out/bin/freesasa_zig --classifier=naccess input.json output.json

# ProtOr分類器（混成軌道ベース）
./zig-out/bin/freesasa_zig --classifier=protor input.json output.json

# OONS分類器（旧FreeSASAデフォルト）
./zig-out/bin/freesasa_zig --classifier=oons input.json output.json
```

### カスタム設定ファイル

```bash
# FreeSASA形式の設定ファイルを使用
./zig-out/bin/freesasa_zig --config=my_radii.config input.json output.json
```

### 入力要件

分類器を使用するには、入力JSONに`residue`と`atom_name`フィールドが必要です:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["CA", "CB", "C"]
}
```

分類器は元の`r`値を上書きし、残基/原子名に基づいて半径を割り当てます。

### 優先順位

`--config`と`--classifier`の両方が指定された場合、`--config`が優先されます。

## API

### 分類器の選択

```zig
const classifier = @import("classifier.zig");

// 分類器タイプ
const ClassifierType = classifier.ClassifierType;
// .naccess - NACCESS互換（デフォルト）
// .protor  - ProtOr（混成軌道ベース）
// .oons    - OONS（旧FreeSASAデフォルト）

// 文字列から変換
const ct = ClassifierType.fromString("naccess");  // .naccess
const name = ClassifierType.naccess.name();       // "NACCESS"
```

### NACCESS分類器

```zig
const naccess = @import("classifier_naccess.zig");

// 半径を取得（残基特異的 → ANY → 元素推定）
const radius: ?f64 = naccess.getRadius("ALA", "CA");  // 1.87

// 極性クラスを取得（残基特異的 → ANY → unknown）
const class: AtomClass = naccess.getClass("ALA", "O");  // .polar

// 両方を取得
const props: ?AtomProperties = naccess.getProperties("CYS", "SG");
// props.radius = 1.85, props.class = .apolar
```

### ProtOr分類器

混成状態に基づく半径（Tsai et al. 1999）。

```zig
const protor = @import("classifier_protor.zig");

// 半径を取得（残基特異的 → 元素推定、ANYフォールバックなし）
const radius: ?f64 = protor.getRadius("ALA", "CA");  // 1.88

// ProtOrでは硫黄は極性
const class: AtomClass = protor.getClass("CYS", "SG");  // .polar
```

### OONS分類器

OONS半径（旧FreeSASAデフォルト）。脂肪族炭素の半径が大きい。

```zig
const oons = @import("classifier_oons.zig");

// 半径を取得（残基特異的 → ANY → 元素推定）
const radius: ?f64 = oons.getRadius("ALA", "CA");  // 2.00（NACASSは1.87）

// OONSではC_CARとSは極性
const class: AtomClass = oons.getClass("ALA", "C");  // .polar
```

### 元素ベースの推定

```zig
const classifier = @import("classifier.zig");

// 原子番号から半径を推定（CA/Caの曖昧さを解決）
const radius = classifier.guessRadiusFromAtomicNumber(6);   // 1.70（炭素）
const radius = classifier.guessRadiusFromAtomicNumber(20);  // 2.31（カルシウム）

// 元素記号から半径を推定
const radius = classifier.guessRadius("C");   // 1.70
const radius = classifier.guessRadius("FE");  // 1.26

// 原子名から元素を抽出して半径を推定
const radius = classifier.guessRadiusFromAtomName(" CA "); // 1.70 (C)
const radius = classifier.guessRadiusFromAtomName("FE  "); // 1.26 (FE)

// 原子名から元素記号を抽出
const elem = classifier.extractElement(" CA "); // "C"（先頭スペース → 1文字元素）
const elem = classifier.extractElement("FE  "); // "FE"（先頭スペースなし → 2文字チェック）
```

### 原子番号による明示的な元素特定

入力JSONに`element`フィールド（原子番号配列）を含めると、原子名の曖昧さを解決できます:

```json
{
  "x": [1.0, 2.0],
  "y": [3.0, 4.0],
  "z": [5.0, 6.0],
  "r": [1.7, 2.31],
  "atom_name": ["CA", "CA"],
  "element": [6, 20]
}
```

上記の例では:
- 1番目のCA: 原子番号6 = 炭素 (Cα)
- 2番目のCA: 原子番号20 = カルシウム（金属イオン）

### データ型

```zig
// 極性クラス
pub const AtomClass = enum {
    polar,   // 親水性
    apolar,  // 疎水性
    unknown, // 不明
};

// 原子プロパティ
pub const AtomProperties = struct {
    radius: f64,
    class: AtomClass,
};
```

## 分類器比較

| 特性 | NACCESS | ProtOr | OONS |
|----------|---------|--------|------|
| 脂肪族C | 1.87 Å | 1.88 Å | 2.00 Å |
| 芳香族C | 1.76 Å | 1.76 Å | 1.75 Å |
| 窒素 | 1.65 Å | 1.64 Å | 1.55 Å |
| 酸素 | 1.40 Å | 1.42-1.46 Å | 1.40 Å |
| 硫黄 | 1.85 Å (非極性) | 1.77 Å (極性) | 2.00 Å (極性) |
| ANYフォールバック | あり | なし | あり |
| 分類基準 | 原子タイプ | 混成軌道 | 原子タイプ |
| 参考文献 | João Rodrigues | Tsai et al. 1999 | Ooi et al. |

## NACCESS原子タイプ

FreeSASAのnaccess.configに基づく原子タイプ定義:

| タイプ | 半径 (Å) | クラス | 説明 | 例 |
|------|------------|-------|-------------|----------|
| C_ALI | 1.87 | 非極性 | 脂肪族炭素 | CA, CB, CG (脂肪族) |
| C_CAR | 1.76 | 非極性 | カルボニル/芳香族炭素 | C (主鎖), CG (芳香族) |
| C_NUC | 1.80 | 非極性 | 核酸炭素 | C1', C2', 塩基炭素 |
| N_AMD | 1.65 | 極性 | アミド窒素 | N (主鎖), NE, NH1/2 |
| N_AMN | 1.50 | 極性 | アミノ窒素 | NZ (LYS) |
| N_NUC | 1.60 | 極性 | 核酸窒素 | N1, N3, N7, N9 |
| O | 1.40 | 極性 | 酸素 | O, OG, OH |
| S | 1.85 | 非極性 | 硫黄 | SG (CYS), SD (MET) |
| SE | 1.80 | 非極性 | セレン | SE (SEC, MSE) |
| P | 1.90 | 非極性 | リン | P (核酸主鎖) |

## ProtOr原子タイプ

Tsai et al. 1999に基づく混成軌道ベースの原子タイプ:

| タイプ | 半径 (Å) | クラス | 説明 |
|------|------------|-------|-------------|
| C3H0 | 1.61 | 非極性 | sp2炭素、水素なし |
| C3H1 | 1.76 | 非極性 | sp2炭素、水素1個 |
| C4H1 | 1.88 | 非極性 | sp3炭素、水素1個 |
| C4H2 | 1.88 | 非極性 | sp3炭素、水素2個 |
| C4H3 | 1.88 | 非極性 | sp3炭素、水素3個（メチル） |
| N3H0 | 1.64 | 極性 | sp2窒素、水素なし |
| N3H1 | 1.64 | 極性 | sp2窒素、水素1個 |
| N3H2 | 1.64 | 極性 | sp2窒素、水素2個（アミド） |
| N4H3 | 1.64 | 極性 | sp3窒素、水素3個（アミノ） |
| O1H0 | 1.42 | 極性 | カルボニル酸素 |
| O2H1 | 1.46 | 極性 | ヒドロキシル酸素 |
| S2H0 | 1.77 | 極性 | チオエーテル硫黄 |
| S2H1 | 1.77 | 極性 | チオール硫黄 |
| SE2H0 | 1.90 | 極性 | セレノエーテル (MSE) |
| SE2H1 | 1.90 | 極性 | セレノール (SEC) |
| P4H0 | 1.80 | 極性 | リン |

## OONS原子タイプ

Ooi, Oobatake, Nemethy, Scheragaに基づく原子タイプ:

| タイプ | 半径 (Å) | クラス | 説明 |
|------|------------|-------|-------------|
| C_ALI | 2.00 | 非極性 | 脂肪族炭素 |
| C_ARO | 1.75 | 非極性 | 芳香族炭素 |
| C_CAR | 1.55 | 極性 | カルボニル炭素 |
| N | 1.55 | 極性 | 窒素 |
| O | 1.40 | 極性 | 酸素 |
| S | 2.00 | 極性 | 硫黄 |
| P | 1.80 | 極性 | リン |
| SE | 1.90 | 極性 | セレン |
| U_POL | 1.50 | 極性 | 不明極性（ASX, GLX用） |

## 対応残基

### アミノ酸（標準20種 + 非標準2種）

| 残基 | 説明 | 備考 |
|---------|-------------|-------|
| ALA | アラニン | |
| ARG | アルギニン | NH1/NH2はN_AMD |
| ASN | アスパラギン | |
| ASP | アスパラギン酸 | |
| CYS | システイン | SGはS |
| GLN | グルタミン | |
| GLU | グルタミン酸 | |
| GLY | グリシン | |
| HIS | ヒスチジン | 環はC_CAR/N_AMD |
| ILE | イソロイシン | |
| LEU | ロイシン | |
| LYS | リシン | NZはN_AMN |
| MET | メチオニン | SDはS |
| PHE | フェニルアラニン | 環はC_CAR |
| PRO | プロリン | |
| SER | セリン | OGはO |
| THR | トレオニン | OG1はO |
| TRP | トリプトファン | NE1はN_AMD |
| TYR | チロシン | OHはO |
| VAL | バリン | |
| SEC | セレノシステイン | SEはSE |
| MSE | セレノメチオニン | SEはSE |

### 核酸

| 残基 | 説明 |
|---------|-------------|
| A, DA | アデニン (RNA/DNA) |
| C, DC | シトシン (RNA/DNA) |
| G, DG | グアニン (RNA/DNA) |
| I, DI | イノシン (RNA/DNA) |
| T, DT | チミン (RNA/DNA) |
| U, DU | ウラシル (RNA/DNA) |

## PDB原子名規則

PDB形式の原子名は4文字で、位置に意味があります:

```
列13-16: 原子名
- 先頭スペース + 1文字元素: " CA ", " N  ", " O  "
- 2文字元素（先頭スペースなし）: "FE  ", "ZN  ", "CA  " (カルシウム)
```

`extractElement`関数はこの規則に従います:
- 先頭スペース → 1文字元素 (" CA " → "C")
- 先頭スペースなし → 2文字元素をチェック ("CA  " → "CA", "CD1 " → "C")

## 元素半径表

Mantina et al. 2009に基づくファンデルワールス半径:

| 元素 | 半径 (Å) | 元素 | 半径 (Å) |
|---------|------------|---------|------------|
| H | 1.10 | FE | 1.26 |
| C | 1.70 | ZN | 1.39 |
| N | 1.55 | CU | 1.40 |
| O | 1.52 | MG | 1.73 |
| P | 1.80 | CA | 2.31 |
| S | 1.80 | NA | 2.27 |
| SE | 1.90 | K | 2.75 |

## 実装詳細

### O(1)検索

`std.StaticStringMap`を使用してコンパイル時ハッシュマップを構築:

```zig
const residue_atoms = std.StaticStringMap(AtomType).initComptime(.{
    .{ "ALA :CB  ", atom_types.C_ALI },
    .{ "ARG :CG  ", atom_types.C_ALI },
    // ...
});
```

キー形式: `"RES :ATOM"` (9文字、スペースパディング)

### メモリ使用量

- すべてのデータはコンパイル時に埋め込み
- ランタイムでのアロケーション不要
- 約150エントリ（アミノ酸 + 核酸）

## 設定ファイル形式

FreeSASA互換の設定ファイル形式をサポート。カスタム原子半径を定義できます。

### 基本形式

```
# コメント
name: MyClassifier

types:
C_ALI 1.87 apolar    # タイプ名、半径、クラス (polar/apolar)
C_CAR 1.76 apolar
O     1.40 polar

atoms:
ANY C   C_CAR        # 残基、原子名、タイプ参照
ANY O   O
ANY CA  C_ALI
ALA CB  C_ALI        # 残基特異的定義
```

### セクション

| セクション | 説明 |
|---------|-------------|
| `name:` | 分類器名（オプション） |
| `types:` | 原子タイプ定義（名前、半径、クラス） |
| `atoms:` | (残基, 原子名) → タイプマッピング |

### ライブラリ使用

```zig
const parser = @import("classifier_parser.zig");

// ファイルから読み込み
var classifier = try parser.parseConfigFile(allocator, "my_radii.config");
defer classifier.deinit();

// 使用
const radius = classifier.getRadius("ALA", "CA");
```

## 実装済み機能

- [x] ProtOr分類器
- [x] OONS分類器
- [x] 設定ファイルパーサー（FreeSASA互換）
- [x] CLI統合（`--classifier`/`--config`オプション）

## 参考文献

- [FreeSASA naccess.config](https://github.com/mittinatten/freesasa/blob/master/share/naccess.config)
- [FreeSASA protor.config](https://github.com/mittinatten/freesasa/blob/master/share/protor.config)
- [FreeSASA oons.config](https://github.com/mittinatten/freesasa/blob/master/share/oons.config)
- Tsai, J., Taylor, R., Chothia, C., & Gerstein, M. (1999). The packing density in proteins: standard radii and volumes. J. Mol. Biol. 290:253-266.
- Mantina et al. (2009) Consistent van der Waals Radii for the Whole Main Group

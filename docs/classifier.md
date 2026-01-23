# 原子分類器（Atom Classifier）

原子の残基名と原子名に基づいてvan der Waals半径と極性クラスを割り当てるモジュール。

## 概要

SASA計算には各原子のvan der Waals半径が必要。分類器はFreeSASAの設定ファイル形式に準拠し、以下の順序で半径を決定する：

1. **残基固有マッチ**: (残基名, 原子名) の完全一致
2. **ANYフォールバック**: 全残基共通の原子定義
3. **元素推定**: 原子名から元素記号を抽出し、van der Waals半径を推定

## モジュール構成

| ファイル | 説明 |
|---------|------|
| `classifier.zig` | コアデータ構造、元素ベース推定、ClassifierType |
| `classifier_naccess.zig` | NACCESS互換組み込み分類器 |
| `classifier_protor.zig` | ProtOr分類器（ハイブリダイゼーションベース） |
| `classifier_oons.zig` | OONS分類器（旧FreeSASAデフォルト） |

## API

### 分類器の選択

```zig
const classifier = @import("classifier.zig");

// 分類器タイプ
const ClassifierType = classifier.ClassifierType;
// .naccess - NACCESS互換（デフォルト）
// .protor  - ProtOr（ハイブリダイゼーションベース）
// .oons    - OONS（旧FreeSASAデフォルト）

// 文字列から変換
const ct = ClassifierType.fromString("naccess");  // .naccess
const name = ClassifierType.naccess.name();       // "NACCESS"
```

### NACCESS分類器

```zig
const naccess = @import("classifier_naccess.zig");

// 半径取得（残基固有 → ANY → 元素推定）
const radius: ?f64 = naccess.getRadius("ALA", "CA");  // 1.87

// 極性クラス取得（残基固有 → ANY → unknown）
const class: AtomClass = naccess.getClass("ALA", "O");  // .polar

// 両方取得
const props: ?AtomProperties = naccess.getProperties("CYS", "SG");
// props.radius = 1.85, props.class = .apolar
```

### ProtOr分類器

ハイブリダイゼーション状態に基づく半径（Tsai et al. 1999）。

```zig
const protor = @import("classifier_protor.zig");

// 半径取得（残基固有 → 元素推定、ANYフォールバックなし）
const radius: ?f64 = protor.getRadius("ALA", "CA");  // 1.88

// ProtOrでは硫黄は極性
const class: AtomClass = protor.getClass("CYS", "SG");  // .polar
```

### OONS分類器

OONS半径（旧FreeSASAデフォルト）。脂肪族炭素が大きめ。

```zig
const oons = @import("classifier_oons.zig");

// 半径取得（残基固有 → ANY → 元素推定）
const radius: ?f64 = oons.getRadius("ALA", "CA");  // 2.00（NACCESSは1.87）

// OONSではC_CARとSは極性
const class: AtomClass = oons.getClass("ALA", "C");  // .polar
```

### 元素ベース推定

```zig
const classifier = @import("classifier.zig");

// 原子番号から半径推定（CA/Ca問題を解決）
const radius = classifier.guessRadiusFromAtomicNumber(6);   // 1.70 (Carbon)
const radius = classifier.guessRadiusFromAtomicNumber(20);  // 2.31 (Calcium)

// 元素記号から半径推定
const radius = classifier.guessRadius("C");   // 1.70
const radius = classifier.guessRadius("FE");  // 1.26

// 原子名から元素抽出して半径推定
const radius = classifier.guessRadiusFromAtomName(" CA "); // 1.70 (C)
const radius = classifier.guessRadiusFromAtomName("FE  "); // 1.26 (FE)

// 原子名から元素記号抽出
const elem = classifier.extractElement(" CA "); // "C"（先頭スペース→1文字元素）
const elem = classifier.extractElement("FE  "); // "FE"（先頭スペースなし→2文字チェック）
```

### 原子番号による明確な元素識別

入力JSONに`element`フィールド（原子番号配列）を含めることで、原子名の曖昧さを解消できる：

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

上記の例では：
- 1番目のCA: 原子番号6 = 炭素（Cα）
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

## 分類器の比較

| 項目 | NACCESS | ProtOr | OONS |
|------|---------|--------|------|
| C脂肪族 | 1.87 Å | 1.88 Å | 2.00 Å |
| C芳香族 | 1.76 Å | 1.76 Å | 1.75 Å |
| 窒素 | 1.65 Å | 1.64 Å | 1.55 Å |
| 酸素 | 1.40 Å | 1.42-1.46 Å | 1.40 Å |
| 硫黄 | 1.85 Å (apolar) | 1.77 Å (polar) | 2.00 Å (polar) |
| ANYフォールバック | あり | なし | あり |
| 分類方式 | 原子タイプ | ハイブリダイゼーション | 原子タイプ |
| 参照 | João Rodrigues | Tsai et al. 1999 | Ooi et al. |

## NACCESS原子タイプ

FreeSASAのnaccess.configに基づく原子タイプ定義：

| タイプ | 半径 (Å) | クラス | 説明 | 使用例 |
|--------|----------|--------|------|--------|
| C_ALI | 1.87 | apolar | 脂肪族炭素 | CA, CB, CG (脂肪族) |
| C_CAR | 1.76 | apolar | カルボニル/芳香族炭素 | C (主鎖), CG (芳香族) |
| C_NUC | 1.80 | apolar | 核酸炭素 | C1', C2', N塩基炭素 |
| N_AMD | 1.65 | polar | アミド窒素 | N (主鎖), NE, NH1/2 |
| N_AMN | 1.50 | polar | アミノ窒素 | NZ (LYS) |
| N_NUC | 1.60 | polar | 核酸窒素 | N1, N3, N7, N9 |
| O | 1.40 | polar | 酸素 | O, OG, OH |
| S | 1.85 | apolar | 硫黄 | SG (CYS), SD (MET) |
| SE | 1.80 | apolar | セレン | SE (SEC, MSE) |
| P | 1.90 | apolar | リン | P (核酸バックボーン) |

## ProtOr原子タイプ

Tsai et al. 1999に基づくハイブリダイゼーションベースの原子タイプ：

| タイプ | 半径 (Å) | クラス | 説明 |
|--------|----------|--------|------|
| C3H0 | 1.61 | apolar | sp2炭素、水素なし |
| C3H1 | 1.76 | apolar | sp2炭素、水素1個 |
| C4H1 | 1.88 | apolar | sp3炭素、水素1個 |
| C4H2 | 1.88 | apolar | sp3炭素、水素2個 |
| C4H3 | 1.88 | apolar | sp3炭素、水素3個（メチル） |
| N3H0 | 1.64 | polar | sp2窒素、水素なし |
| N3H1 | 1.64 | polar | sp2窒素、水素1個 |
| N3H2 | 1.64 | polar | sp2窒素、水素2個（アミド） |
| N4H3 | 1.64 | polar | sp3窒素、水素3個（アミノ） |
| O1H0 | 1.42 | polar | カルボニル酸素 |
| O2H1 | 1.46 | polar | ヒドロキシル酸素 |
| S2H0 | 1.77 | polar | チオエーテル硫黄 |
| S2H1 | 1.77 | polar | チオール硫黄 |
| SE2H0 | 1.90 | polar | セレノエーテル（MSE） |
| SE2H1 | 1.90 | polar | セレノール（SEC） |
| P4H0 | 1.80 | polar | リン |

## OONS原子タイプ

Ooi, Oobatake, Nemethy, Scheragaに基づく原子タイプ：

| タイプ | 半径 (Å) | クラス | 説明 |
|--------|----------|--------|------|
| C_ALI | 2.00 | apolar | 脂肪族炭素 |
| C_ARO | 1.75 | apolar | 芳香族炭素 |
| C_CAR | 1.55 | polar | カルボニル炭素 |
| N | 1.55 | polar | 窒素 |
| O | 1.40 | polar | 酸素 |
| S | 2.00 | polar | 硫黄 |
| P | 1.80 | polar | リン |
| SE | 1.90 | polar | セレン |
| U_POL | 1.50 | polar | 不明極性（ASX, GLX用） |

## 対応残基

### アミノ酸（20種 + 非標準2種）

| 残基 | 説明 | 特記事項 |
|------|------|----------|
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
| THR | スレオニン | OG1はO |
| TRP | トリプトファン | NE1はN_AMD |
| TYR | チロシン | OHはO |
| VAL | バリン | |
| SEC | セレノシステイン | SEはSE |
| MSE | セレノメチオニン | SEはSE |

### 核酸

| 残基 | 説明 |
|------|------|
| A, DA | アデニン（RNA/DNA） |
| C, DC | シトシン（RNA/DNA） |
| G, DG | グアニン（RNA/DNA） |
| I, DI | イノシン（RNA/DNA） |
| T, DT | チミン（RNA/DNA） |
| U, DU | ウラシル（RNA/DNA） |

## PDB原子名規則

PDB形式の原子名は4文字で、位置に意味がある：

```
Column 13-16: 原子名
- 先頭スペース + 1文字元素: " CA ", " N  ", " O  "
- 2文字元素（スペースなし）: "FE  ", "ZN  ", "CA  "（カルシウム）
```

`extractElement`関数はこの規則に従い：
- 先頭スペースあり → 1文字元素（" CA " → "C"）
- 先頭スペースなし → 2文字元素をチェック（"CA  " → "CA"、"CD1 " → "C"）

## 元素半径テーブル

Mantina et al. 2009に基づくvan der Waals半径：

| 元素 | 半径 (Å) | 元素 | 半径 (Å) |
|------|----------|------|----------|
| H | 1.10 | FE | 1.26 |
| C | 1.70 | ZN | 1.39 |
| N | 1.55 | CU | 1.40 |
| O | 1.52 | MG | 1.73 |
| P | 1.80 | CA | 2.31 |
| S | 1.80 | NA | 2.27 |
| SE | 1.90 | K | 2.75 |

## 実装詳細

### O(1)ルックアップ

`std.StaticStringMap`を使用してコンパイル時にハッシュマップを構築：

```zig
const residue_atoms = std.StaticStringMap(AtomType).initComptime(.{
    .{ "ALA :CB  ", atom_types.C_ALI },
    .{ "ARG :CG  ", atom_types.C_ALI },
    // ...
});
```

キー形式: `"RES :ATOM"` (9文字、スペースパディング)

### メモリ使用

- 全データはコンパイル時に埋め込み
- 実行時アロケーション不要
- 約150エントリ（アミノ酸 + 核酸）

## 今後の予定

- [x] ProtOr分類器
- [x] OONS分類器
- [ ] 設定ファイルパーサー（FreeSASA互換）
- [ ] CLI統合（`--classifier=naccess`オプション）

## 参考

- [FreeSASA naccess.config](https://github.com/mittinatten/freesasa/blob/master/share/naccess.config)
- [FreeSASA protor.config](https://github.com/mittinatten/freesasa/blob/master/share/protor.config)
- [FreeSASA oons.config](https://github.com/mittinatten/freesasa/blob/master/share/oons.config)
- Tsai, J., Taylor, R., Chothia, C., & Gerstein, M. (1999). The packing density in proteins: standard radii and volumes. J. Mol. Biol. 290:253-266.
- Mantina et al. (2009) Consistent van der Waals Radii for the Whole Main Group

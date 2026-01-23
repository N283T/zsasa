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
| `classifier.zig` | コアデータ構造、元素ベース推定 |
| `classifier_naccess.zig` | NACCESS互換組み込み分類器 |

## API

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

### 元素ベース推定

```zig
const classifier = @import("classifier.zig");

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

- [ ] ProtOr分類器
- [ ] OONS分類器
- [ ] 設定ファイルパーサー（FreeSASA互換）
- [ ] CLI統合（`--classifier=naccess`オプション）

## 参考

- [FreeSASA naccess.config](https://github.com/mittinatten/freesasa/blob/master/share/naccess.config)
- Mantina et al. (2009) Consistent van der Waals Radii for the Whole Main Group

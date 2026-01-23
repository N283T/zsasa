# Phase 9: Radius Classifier

## Goal

FreeSASAの原子半径分類アルゴリズムを移植し、PDB/mmCIFから直接SASA計算できるようにする。

## 優先度: 高

## 参考: FreeSASA実装

FreeSASAの半径分類器（`src/classifier*.c`）:
- `classifier.c` - 基本分類ロジック、設定ファイルパーサー
- `classifier_naccess.c` - NACCESS互換半径（自動生成）
- `classifier_oons.c` - OONS半径セット（自動生成）
- `classifier_protor.c` - ProtOr半径セット（自動生成）

設定ファイル（`share/`）:
- `naccess.config` - NACCESS互換（João Rodrigues提供）
- `oons.config` - OONS半径
- `protor.config` - ProtOr半径（Tsai et al. 1999）
- `dssp.config` - DSSP互換

### FreeSASA設定ファイル形式

```
name: NACCESS

types:
C_ALI 1.87 apolar    # Type名 半径 クラス(polar/apolar)
C_CAR 1.76 apolar
N_AMN 1.50 polar
...

atoms:
ANY C   C_CAR        # ANY = 全残基共通
ANY O   O
ANY CA  C_ALI
ALA CB  C_ALI        # 残基固有
ARG CG  C_ALI
...
```

### 半径ルックアップの流れ

1. 残基名 + 原子名で検索
2. 見つからない → `ANY` + 原子名で検索
3. 見つからない → 元素記号から推定（`freesasa_guess_radius`）

---

## Phase 9.0: Input Format Extension ✅

**目標**: 入力JSONに残基名・原子名フィールドを追加

### Tasks

- [x] `cif_to_input_json.py` を更新し`residue`/`atom_name`を出力
- [x] `types.zig` の`AtomInput`に`residue`/`atom_name`フィールド追加
- [x] `json_parser.zig` でオプショナルフィールドをパース
- [x] テスト追加
- [x] README更新

### 変更後の入力形式

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

`residue`と`atom_name`は後方互換性のためオプション。

---

## Phase 9.1: Core Data Structures ✅

**目標**: 分類器の基本データ構造を実装

### Tasks

- [x] `AtomClass` enum定義（polar, apolar, unknown）
- [x] `AtomProperties` struct定義（radius, class）
- [x] `AtomKey` struct定義（residue + atom_name のハッシュマップキー）
- [x] `Classifier` struct定義
- [x] 基本テスト

### API設計

```zig
pub const AtomClass = enum {
    polar,
    apolar,
    unknown,
};

pub const AtomType = struct {
    name: []const u8,
    radius: f64,
    class: AtomClass,
};

pub const Classifier = struct {
    name: []const u8,
    residues: std.StringHashMap(ResidueConfig),
    any_residue: ?*ResidueConfig,  // ANYフォールバック

    pub fn getRadius(self: *const Classifier, res_name: []const u8, atom_name: []const u8) ?f64;
    pub fn getClass(self: *const Classifier, res_name: []const u8, atom_name: []const u8) AtomClass;
    pub fn deinit(self: *Classifier, allocator: Allocator) void;
};
```

### Files

| File | Action |
|------|--------|
| `src/classifier.zig` | CREATE |

### Success Criteria

- [x] データ構造が定義されている
- [x] 単純なルックアップテストが通る

---

## Phase 9.2: Element-Based Radius Guessing

**目標**: 元素記号からvan der Waals半径を推定（最終フォールバック）

### Tasks

- [ ] `guessRadius(element: []const u8) ?f64` 実装
- [ ] 主要元素のvdW半径テーブル（Mantina et al. 2009）
- [ ] PDB原子名から元素記号抽出ロジック
- [ ] テスト

### 元素半径テーブル（一部）

| Element | Radius (Å) | Source |
|---------|------------|--------|
| H | 1.10 | Mantina et al. |
| C | 1.70 | Mantina et al. |
| N | 1.55 | Mantina et al. |
| O | 1.52 | Mantina et al. |
| P | 1.80 | Mantina et al. |
| S | 1.80 | Mantina et al. |
| SE | 1.90 | gemmi |
| F | 1.47 | Mantina et al. |
| CL | 1.75 | Mantina et al. |
| ... | ... | ... |

### API

```zig
/// 元素記号からvdW半径を推定
pub fn guessRadius(element: []const u8) ?f64;

/// PDB原子名から元素記号を抽出
/// e.g., " CA " -> "C", " NE2" -> "N", "SE  " -> "SE"
pub fn extractElement(atom_name: []const u8) []const u8;
```

### Files

| File | Action |
|------|--------|
| `src/classifier.zig` | MODIFY |

### Success Criteria

- [ ] 主要元素（H, C, N, O, S, P）の半径が正しい
- [ ] PDB形式の原子名から元素を正しく抽出

---

## Phase 9.3: Built-in NACCESS Classifier

**目標**: NACCESS互換の組み込み分類器を実装

### Tasks

- [ ] NACCESS半径データを`comptime`で埋め込み
- [ ] `naccess_classifier` を`pub const`で公開
- [ ] ANYフォールバックの実装
- [ ] FreeSASAと同じ半径が返ることを検証

### 実装方針

```zig
// comptime埋め込みデータ
const naccess_types = [_]AtomType{
    .{ .name = "C_ALI", .radius = 1.87, .class = .apolar },
    .{ .name = "C_CAR", .radius = 1.76, .class = .apolar },
    // ...
};

const naccess_atoms = [_]struct { res: []const u8, atom: []const u8, type_idx: usize }{
    .{ .res = "ANY", .atom = "C", .type_idx = 1 },
    .{ .res = "ANY", .atom = "O", .type_idx = 6 },
    // ...
};

pub const naccess_classifier = Classifier.initBuiltin(&naccess_types, &naccess_atoms);
```

### Files

| File | Action |
|------|--------|
| `src/classifier_naccess.zig` | CREATE |

### Success Criteria

- [ ] NACCESS分類器が利用可能
- [ ] 標準アミノ酸の全原子に半径が割り当てられる
- [ ] ANYフォールバックが動作

---

## Phase 9.4: Built-in ProtOr & OONS Classifiers

**目標**: ProtOrとOONSの組み込み分類器を追加

### Tasks

- [ ] ProtOr半径データ埋め込み
- [ ] OONS半径データ埋め込み
- [ ] `ClassifierType` enumと選択API

### API

```zig
pub const ClassifierType = enum {
    naccess,
    protor,
    oons,
};

pub fn getBuiltinClassifier(classifier_type: ClassifierType) *const Classifier;
```

### Files

| File | Action |
|------|--------|
| `src/classifier_protor.zig` | CREATE |
| `src/classifier_oons.zig` | CREATE |
| `src/classifier.zig` | MODIFY |

### Success Criteria

- [ ] 3種類の分類器が利用可能
- [ ] 各分類器でFreeSASAと同じ半径

---

## Phase 9.5: Config File Parser

**目標**: FreeSASA互換の設定ファイルパーサー

### Tasks

- [ ] `types:` セクションパース
- [ ] `atoms:` セクションパース
- [ ] `name:` セクションパース（オプション）
- [ ] コメント（`#`）対応
- [ ] エラーメッセージ

### API

```zig
pub fn parseConfig(allocator: Allocator, content: []const u8) !Classifier;
pub fn parseConfigFile(allocator: Allocator, path: []const u8) !Classifier;
```

### Files

| File | Action |
|------|--------|
| `src/classifier_parser.zig` | CREATE |

### Success Criteria

- [ ] FreeSASAの`naccess.config`をパースできる
- [ ] 不正な設定ファイルで適切なエラー

---

## Phase 9.6: CLI Integration

**目標**: CLIオプションとの統合

### Tasks

- [ ] `--classifier=TYPE` オプション追加（naccess/protor/oons）
- [ ] `--config=FILE` オプション追加（カスタム設定）
- [ ] 新しい入力形式の検討（座標+残基名+原子名）
- [ ] ヘルプメッセージ更新

### 新しい入力形式案

```json
{
  "atoms": [
    {"x": 1.0, "y": 2.0, "z": 3.0, "residue": "ALA", "name": "CA"},
    {"x": 2.0, "y": 3.0, "z": 4.0, "residue": "ALA", "name": "CB"},
    ...
  ]
}
```

または簡易形式:

```json
{
  "x": [1.0, 2.0, ...],
  "y": [2.0, 3.0, ...],
  "z": [3.0, 4.0, ...],
  "residue": ["ALA", "ALA", ...],
  "atom_name": ["CA", "CB", ...]
}
```

### Files

| File | Action |
|------|--------|
| `src/main.zig` | MODIFY |
| `src/json_parser.zig` | MODIFY |

### Success Criteria

- [ ] `--classifier=naccess` で分類器選択
- [ ] `--config=custom.config` でカスタム設定読み込み
- [ ] 残基/原子名入力時に自動半径割り当て

---

## Phase 9.7: Documentation & Testing

**目標**: ドキュメントと検証

### Tasks

- [ ] README更新（classifier使用法）
- [ ] docs/classifier.md 作成
- [ ] FreeSASAとの比較検証スクリプト
- [ ] 全アミノ酸の半径検証テスト

### Success Criteria

- [ ] ドキュメントが完備
- [ ] FreeSASAと同等の半径割り当て

---

## Summary

| Phase | 内容 | 依存 | 状態 |
|-------|------|------|------|
| 9.0 | Input Format Extension | - | ✅ 完了 |
| 9.1 | Core Data Structures | 9.0 | ✅ 完了 |
| 9.2 | Element-Based Guessing | 9.1 | |
| 9.3 | NACCESS Classifier | 9.1, 9.2 | |
| 9.4 | ProtOr & OONS | 9.3 | |
| 9.5 | Config Parser | 9.1 | |
| 9.6 | CLI Integration | 9.3, 9.5 | |
| 9.7 | Docs & Testing | 9.6 | |

## Files Summary

| File | Action | Phase |
|------|--------|-------|
| `scripts/cif_to_input_json.py` | MODIFY | 9.0 |
| `src/types.zig` | MODIFY | 9.0 |
| `src/json_parser.zig` | MODIFY | 9.0, 9.6 |
| `src/classifier.zig` | CREATE | 9.1-9.2 |
| `src/classifier_naccess.zig` | CREATE | 9.3 |
| `src/classifier_protor.zig` | CREATE | 9.4 |
| `src/classifier_oons.zig` | CREATE | 9.4 |
| `src/classifier_parser.zig` | CREATE | 9.5 |
| `src/main.zig` | MODIFY | 9.6 |
| `docs/classifier.md` | CREATE | 9.7 |
| `README.md` | MODIFY | 9.0, 9.7 |

---
- [ ] **DONE** - Phase 9 complete

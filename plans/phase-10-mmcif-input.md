# Phase 10: mmCIF Direct Input

## Goal

mmCIF形式の構造ファイルを直接読み込み、JSON変換なしでSASA計算できるようにする。

**スコープ**: `atom_site`データ抽出に特化した簡易CIFパーサー実装。
完全なmmCIFパーサーは別プロジェクトとして後日実装する。

## 優先度: 高

## 参考実装: gemmi

gemmiの実装を参考にする（`~/ghq/github.com/project-gemmi/gemmi/`）。

### gemmi CIF構造 (`include/gemmi/cifdoc.hpp`)

```cpp
struct Loop {
  std::vector<std::string> tags;    // タグ名リスト
  std::vector<std::string> values;  // 値（row-major順）
  // ...
};
```

### gemmi atom_site抽出 (`src/mmcif.cpp:793-823`)

```cpp
enum { kId=0, kGroupPdb, kSymbol, kLabelAtomId, kAltId, kLabelCompId,
       kLabelAsymId, kLabelEntityId, kLabelSeqId, kInsCode,
       kX, kY, kZ, kOcc, kBiso, kCharge, /* ... */ };

cif::Table atom_table = block.find("_atom_site.",
                                   {"id",
                                    "?group_PDB",
                                    "type_symbol",
                                    "?label_atom_id",
                                    "label_alt_id",
                                    "?label_comp_id",
                                    "Cartn_x",
                                    "Cartn_y",
                                    "Cartn_z",
                                    /* ... */});
```

### gemmi Element Enum (`include/gemmi/elem.hpp:13-24`)

```cpp
enum class El : unsigned char {
  X=0,  // unknown element
  H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar,  // 1-18
  K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,  // 19-36
  // ... (enum value = atomic number)
};
```

**Key insight**: Element enumの値がそのまま原子番号として使える。

## 依存関係

- Phase 9 (Radius Classifier) ✅ - 完了済み

## CIF形式の基礎

### データ構造

```
data_XXXX
#
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
C CA ALA 10.000 20.000 30.000
N N  ALA 11.000 21.000 31.000
...
```

### パース戦略（簡易版）

1. `loop_` キーワードを探す
2. `_atom_site.`で始まるタグを収集
3. タグ数を数える
4. 値を読み取り、タグ位置に基づいて割り当て

## Tasks

### Phase 10.1: Element Module

- [ ] Element enum定義（gemmi El enumと同様）
- [ ] Element symbol → atomic number変換
- [ ] atomic number → Element変換

**ファイル**: `src/element.zig`

```zig
pub const Element = enum(u8) {
    X = 0,   // Unknown
    H = 1,
    He = 2,
    Li = 3,
    // ... 118まで

    /// Element symbol (1-2文字) からElementを取得
    pub fn fromSymbol(symbol: []const u8) Element {
        // "C" -> .C, "CA" -> .Ca, "FE" -> .Fe
    }

    /// Atomic numberを取得（enum値そのもの）
    pub fn atomicNumber(self: Element) u8 {
        return @intFromEnum(self);
    }
};
```

### Phase 10.2: CIF Tokenizer

- [ ] 基本トークナイザー実装
- [ ] コメント(`#`)のスキップ
- [ ] 引用符文字列（`'...'`, `"..."`）の処理
- [ ] セミコロンテキスト（`;...;`）の処理
- [ ] 特殊値（`.`, `?`）の処理

**ファイル**: `src/cif_tokenizer.zig`

```zig
pub const Token = union(enum) {
    data_block: []const u8,  // data_XXXX
    loop,                     // loop_
    tag: []const u8,          // _category.field
    value: []const u8,        // 値
    eof,
};

pub const Tokenizer = struct {
    source: []const u8,
    pos: usize,

    pub fn next(self: *Tokenizer) Token { ... }
};
```

### Phase 10.3: AtomSite Parser

- [ ] `loop_`ブロック内の`_atom_site.*`タグ検出
- [ ] 必要なフィールドのインデックス記録
- [ ] 値の読み取りとAtomInput構築

**必要なフィールド（gemmi参照）**:

| フィールド | 用途 | 必須 |
|-----------|------|------|
| `Cartn_x` | X座標 | ✅ |
| `Cartn_y` | Y座標 | ✅ |
| `Cartn_z` | Z座標 | ✅ |
| `type_symbol` | 元素記号 → 半径推測 | ✅ |
| `label_atom_id` or `auth_atom_id` | 原子名（CA等） | ✅ |
| `label_comp_id` or `auth_comp_id` | 残基名（ALA等） | ✅ |
| `group_PDB` | ATOM/HETATM区別 | ⚪ |
| `label_alt_id` | 代替位置 | ⚪ |
| `pdbx_PDB_model_num` | モデル番号 | ⚪ |

**ファイル**: `src/mmcif_parser.zig`

```zig
pub const MmcifParser = struct {
    allocator: Allocator,

    pub fn parse(self: *MmcifParser, source: []const u8) !AtomInput {
        // 1. Tokenize
        // 2. Find atom_site loop
        // 3. Map tag positions
        // 4. Extract values
        // 5. Build AtomInput
    }

    /// Parse from file path
    pub fn parseFile(self: *MmcifParser, path: []const u8) !AtomInput {
        const file = try std.fs.cwd().openFile(path, .{});
        defer file.close();
        const source = try file.readToEndAlloc(self.allocator, 100 * 1024 * 1024);
        defer self.allocator.free(source);
        return self.parse(source);
    }
};
```

### Phase 10.4: AtomInput統合

- [ ] MmcifParser出力をAtomInput形式に変換
- [ ] element配列（atomic number）の設定
- [ ] Classifierとの連携確認

**変換フロー**:
```
mmCIF file
    ↓ parse()
atom_site loop data
    ↓ extract coordinates
    ↓ extract element (type_symbol → Element → u8)
    ↓ extract residue/atom names
AtomInput {
    x, y, z: coordinates
    r: classifier.getRadius() or element-based guess
    residue: label_comp_id or auth_comp_id
    atom_name: label_atom_id or auth_atom_id
    element: atomic numbers from type_symbol
}
```

### Phase 10.5: CLI Integration

- [ ] 入力形式の自動検出（`.cif` / `.json`）
- [ ] MmcifParserの呼び出し追加

**ファイル**: `src/main.zig`

```zig
fn detectInputFormat(path: []const u8) InputFormat {
    if (std.mem.endsWith(u8, path, ".cif")) return .mmcif;
    if (std.mem.endsWith(u8, path, ".mmcif")) return .mmcif;
    return .json;
}

// main()内で
const input = switch (detectInputFormat(input_path)) {
    .mmcif => try mmcif_parser.parseFile(input_path),
    .json => try json_parser.parseFile(input_path),
};
```

### Phase 10.6: Testing

- [ ] Element module unit tests
- [ ] Tokenizer unit tests
- [ ] Parser integration tests
- [ ] Real mmCIF file tests

**テストファイル**:
- `test/mmcif/simple.cif` - 最小限のatom_site
- `test/mmcif/1crn.cif` - 実際のPDBファイル（小さいタンパク質）
- ベンチマーク用PDBファイル（既存）

## Files

| File | Action |
|------|--------|
| `src/element.zig` | CREATE |
| `src/cif_tokenizer.zig` | CREATE |
| `src/mmcif_parser.zig` | CREATE |
| `src/main.zig` | MODIFY |
| `src/root.zig` | MODIFY |
| `test/mmcif/` | CREATE (directory) |

## Future Work（別フェーズ/別プロジェクト）

以下は今回のスコープ外:

- [ ] gzip圧縮ファイルのサポート (.cif.gz) → Phase 10.7
- [ ] 完全なCIFパーサー（save frames, nested blocks等）→ 別プロジェクト
- [ ] PDB format サポート → Phase 11
- [ ] モデル選択オプション (`--model=<n>`) → Phase 10.7
- [ ] チェーンフィルタ (`--chain=<id>`) → Phase 10.7

## Success Criteria

- [ ] `.cif` ファイルから直接SASA計算できる
- [ ] 既存のベンチマークPDBファイルをmmCIFに変換してテスト
- [ ] FreeSASAと同じ原子が選択される（atom_site座標一致）
- [ ] Classifierによる半径割り当てが正しく動作

## 実装順序

1. **Element module** - 最も独立、他に依存なし
2. **CIF Tokenizer** - Element moduleに依存なし
3. **AtomSite Parser** - Tokenizer + Element使用
4. **CLI Integration** - Parser完成後
5. **Testing** - 各フェーズで並行

---
- [ ] **DONE** - Phase 10 complete

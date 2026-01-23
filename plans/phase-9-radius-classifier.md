# Phase 9: Radius Classifier

## Goal

FreeSASAの原子半径分類アルゴリズムを移植し、PDB/mmCIFから直接SASA計算できるようにする。

## 優先度: 高

## 参考: FreeSASA実装

FreeSASAの半径分類器（`src/classifier*.c`）:
- `classifier.c` - 基本分類ロジック
- `classifier_naccess.c` - NACCESS互換半径
- `classifier_oons.c` - OONS半径セット
- `classifier_protor.c` - ProtOr半径セット

設定ファイル（`share/`）:
- `naccess.config`
- `oons.config`
- `protor.config`
- `dssp.config`

## Tasks

### Phase 9.1: Classifier Core

- [ ] 半径分類器のコア実装
- [ ] 残基名 + 原子名 → 半径 のマッピング
- [ ] 不明原子のフォールバック処理

**API設計:**
```zig
pub const Classifier = struct {
    pub fn getRadius(
        residue_name: []const u8,
        atom_name: []const u8,
    ) ?f64;

    pub fn fromConfig(allocator: Allocator, config: []const u8) !Classifier;
};

pub const ClassifierType = enum {
    naccess,  // NACCESS compatible
    oons,     // OONS
    protor,   // ProtOr
};
```

### Phase 9.2: Built-in Classifiers

- [ ] NACCESS半径テーブル実装
- [ ] OONS半径テーブル実装
- [ ] ProtOr半径テーブル実装

**NACCESSの例（一部）:**
| Residue | Atom | Radius (Å) |
|---------|------|------------|
| ALA | N | 1.65 |
| ALA | CA | 1.87 |
| ALA | C | 1.76 |
| ALA | O | 1.40 |
| ALA | CB | 1.87 |
| ... | ... | ... |

### Phase 9.3: Custom Config Parser

- [ ] FreeSASA互換の設定ファイルパーサー
- [ ] `--config=<file>` オプション追加
- [ ] カスタム半径定義のサポート

**設定ファイル形式（FreeSASA互換）:**
```
# NACCESS compatible radii
types:
N 1.65
CA 1.87
C 1.76
O 1.40
...

residues:
ALA: N CA C O CB
GLY: N CA C O
...
```

### Phase 9.4: Integration with Structure Input

- [ ] mmCIF入力との統合（Phase 10と連携）
- [ ] 残基/原子情報の保持
- [ ] `--classifier=<type>` オプション追加

## Files

| File | Action |
|------|--------|
| `src/classifier.zig` | CREATE |
| `src/classifier_naccess.zig` | CREATE |
| `src/classifier_oons.zig` | CREATE |
| `src/classifier_protor.zig` | CREATE |
| `share/naccess.config` | CREATE (optional) |
| `src/main.zig` | MODIFY |

## Success Criteria

- [ ] 3種類の組み込み分類器が動作
- [ ] FreeSASAと同じ半径が割り当てられる
- [ ] カスタム設定ファイルが読み込める

---
- [ ] **DONE** - Phase 9 complete

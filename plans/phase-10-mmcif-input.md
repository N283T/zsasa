# Phase 10: mmCIF Direct Input

## Goal

mmCIF形式の構造ファイルを直接読み込み、JSON変換なしでSASA計算できるようにする。

## 優先度: 中

## 参考: FreeSASA実装

FreeSASAの`src/cif.cc`:
- gemmi-cifライブラリを使用
- `_atom_site`カテゴリから座標を取得
- Model/Chain/Residue/Atom階層を構築

## 依存関係

- Phase 9 (Radius Classifier) - 半径割り当てに必要

## Tasks

### Phase 10.1: mmCIF Parser

- [ ] mmCIF形式の基本パーサー実装
- [ ] `_atom_site`カテゴリの読み取り
- [ ] 座標（x, y, z）の抽出

**必要なmmCIFフィールド:**
```
_atom_site.group_PDB      # ATOM/HETATM
_atom_site.label_atom_id  # 原子名
_atom_site.label_comp_id  # 残基名
_atom_site.label_asym_id  # チェーンID
_atom_site.Cartn_x        # X座標
_atom_site.Cartn_y        # Y座標
_atom_site.Cartn_z        # Z座標
_atom_site.label_seq_id   # 残基番号
_atom_site.type_symbol    # 元素記号
```

### Phase 10.2: Structure Builder

- [ ] Model/Chain/Residue/Atom階層の構築
- [ ] 代替位置（alt_id）の処理
- [ ] モデル選択オプション

**データ構造:**
```zig
pub const Structure = struct {
    models: []Model,

    pub const Model = struct {
        chains: []Chain,
    };

    pub const Chain = struct {
        id: []const u8,
        residues: []Residue,
    };

    pub const Residue = struct {
        name: []const u8,
        seq_id: i32,
        atoms: []Atom,
    };

    pub const Atom = struct {
        name: []const u8,
        element: []const u8,
        coord: Vec3,
        radius: f64,  // Classifierで割り当て
    };
};
```

### Phase 10.3: CLI Integration

- [ ] 入力形式の自動検出（.cif / .json）
- [ ] `--model=<n>` オプション追加
- [ ] `--chain=<id>` フィルタオプション

**CLI例:**
```bash
# mmCIF直接入力
freesasa_zig structure.cif output.json

# 特定モデル/チェーン
freesasa_zig --model=1 --chain=A structure.cif output.json

# 半径分類器指定
freesasa_zig --classifier=naccess structure.cif output.json
```

### Phase 10.4: Compressed File Support

- [ ] gzip圧縮ファイルのサポート (.cif.gz)
- [ ] 自動解凍処理

## Files

| File | Action |
|------|--------|
| `src/mmcif_parser.zig` | CREATE |
| `src/structure.zig` | CREATE |
| `src/main.zig` | MODIFY |

## Alternative: External Library

Zigでのmmcifパーサー実装が複雑な場合:
- gemmiのCバインディングを使用
- またはPythonスクリプトでのプリプロセスを維持

## Success Criteria

- [ ] .cif ファイルから直接SASA計算できる
- [ ] FreeSASAと同じ原子が選択される
- [ ] gzip圧縮ファイルが読み込める

---
- [ ] **DONE** - Phase 10 complete

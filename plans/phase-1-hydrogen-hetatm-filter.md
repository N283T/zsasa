# Phase 1: Hydrogen / HETATM Filtering

## Goal

zsasa の PDB/mmCIF パーサーに水素原子・HETATM のフィルタリングオプションを追加し、FreeSASA・RustSASA とデフォルト挙動を揃える。

## Background

| Parameter | zsasa (current) | FreeSASA | RustSASA |
|-----------|----------------|----------|----------|
| Hydrogen | included (no filter) | excluded (default) | excluded (default) |
| HETATM | included (atom_only=false) | excluded (default) | excluded (default) |

AlphaFold PDB には水素・HETATM が含まれないためバッチベンチマークには影響なし。
実験構造 PDB を扱う場合に差異が出る。

## Tasks

- [x] PDB パーサーに水素原子スキップオプション追加 (`skip_hydrogens: bool = true`)
- [x] mmCIF パーサーにも同様のオプション追加
- [x] CLI に `--include-hydrogens` / `--include-hetatm` フラグ追加
- [x] デフォルト: 水素除外、HETATM 除外 (FreeSASA/RustSASA と一致)
- [x] バッチモード (`batch.zig`) にもフラグ伝播
- [x] テスト追加 (水素あり/なし PDB で atom count が変わることを確認)
- [x] docs/cli.md 更新

## Notes

- `atom_only` フラグは既に `pdb_parser.zig:56` に存在するが、デフォルト false
- 水素判定: element == H or atom_name が H で始まる (PDB column 77-78 or atom name heuristic)
- FreeSASA: `freesasa.h:183-184` で `FREESASA_INCLUDE_HETATM` / `FREESASA_INCLUDE_HYDROGEN` フラグ
- RustSASA: `options.rs:86-98` で `include_hydrogens` / `include_hetatms` (both default false)

---
- [x] **DONE** - Phase complete

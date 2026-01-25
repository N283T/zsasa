# Phase 16: Analysis Features

## Goal

構造生物学で頻用される解析機能を追加する。

## 優先度: 中

## Tasks

### Phase 16.1: Relative SASA (RSA)

- [x] Gly-X-Gly参照値テーブル実装
- [x] 残基ごとのRSA計算（0-100%）
- [x] 埋没/露出判定（閾値: 25%）

**参照値 (Miller et al. 1987, Sander & Rost 1994):**
```zig
const max_asa = struct {
    // Gly-X-Gly tripeptide reference values (Å²)
    pub const ALA = 113.0;
    pub const ARG = 241.0;
    pub const ASN = 158.0;
    pub const ASP = 151.0;
    // ...
};
```

**出力例:**
```
Residue  RSA    State
ALA-1    85.2%  Exposed
VAL-2    12.3%  Buried
LYS-3    67.8%  Exposed
```

### Phase 16.2: Per-Residue Aggregation

- [x] 原子SASA → 残基SASA集約
- [x] 残基番号 + 挿入コード対応
- [ ] Backbone/Sidechain分離（オプション） - 未実装（将来拡張）

**集計ロジック:**
```
Residue SASA = Σ(atom SASA for atoms in residue)
```

### Phase 16.3: Polar/Nonpolar Classification

- [x] 原子の極性分類
- [x] Polar/Nonpolar面積内訳
- [x] 残基タイプ別集計

**分類基準:**
| Type | Atoms |
|------|-------|
| Nonpolar | C (except CO) |
| Polar | N, O (except charged) |
| Charged | Nζ(Lys), Oδ(Asp), etc. |

### Phase 16.4: Interface SASA - REMOVED

~~削除済み~~ - 計算式の定義が曖昧（ChimeraXは/2、FreeSASAは未実装）のため、
ユーザーが`--chain`オプションで個別計算する方式を推奨。

## Files

| File | Action |
|------|--------|
| `src/analysis.zig` | CREATE - RSA, aggregation, classification |
| `src/reference_asa.zig` | CREATE - Gly-X-Gly reference values |
| `src/main.zig` | MODIFY - analysis options |

## Success Criteria

- [x] RSA計算が参照実装と一致
- [x] 残基別SASAが原子SASAの合計と一致

---
- [x] **DONE** - Phase 16 complete

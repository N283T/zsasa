# Phase 8: Benchmark Dataset

## Goal

ベンチマーク用の構造ファイルとテストデータを整備する。

## 優先度: 高

## Tasks

### Phase 8.1: Structure Collection

- [ ] 異なるサイズの構造ファイル収集
- [ ] mmCIF形式で保存
- [ ] 入力JSON生成

**構造リスト:**
| PDB | Atoms | Type | Description |
|-----|-------|------|-------------|
| 1CRN | 327 | Small | Crambin |
| 1UBQ | 660 | Small | Ubiquitin |
| 1A0Q | 3,183 | Medium | Lipid transfer protein |
| 3HHB | 4,779 | Medium | Hemoglobin |
| 1AON | 58,674 | Large | GroEL-GroES complex |
| 4V6X | ~150,000 | XLarge | Ribosome |

## Benchmark Structure Selection Guide

### Size Categories

性能スケーリングを測定するため、異なるサイズの構造が必要:

| Category | Atoms | Purpose | Examples |
|----------|-------|---------|----------|
| Tiny | <500 | 正確性検証、高速テスト | 1CRN (327) |
| Small | 500-2,000 | 基本性能測定 | 1UBQ (660), 2GB1 (855) |
| Medium | 2,000-10,000 | 標準ベンチマーク | 1A0Q (3,183), 3HHB (4,779) |
| Large | 10,000-50,000 | スケーラビリティ | 1HRC (8,619), 4HHB (17,536) |
| XLarge | 50,000+ | 極限性能 | 1AON (58,674), 4V6X (~150,000) |

### Structure Types

異なる種類の構造での動作確認:

| Type | Why | Examples |
|------|-----|----------|
| Globular protein | 標準ケース | 1UBQ, 1CRN |
| Multi-chain complex | 複数chain処理 | 3HHB (4 chains), 1AON (21 chains) |
| Membrane protein | 特殊な形状 | 1OCC, 4HYT |
| Nucleic acid | DNA/RNA対応確認 | 1BNA (DNA), 4TNA (tRNA) |
| Protein-nucleic complex | 混合系 | 1EHZ (tRNA synthetase) |
| Disordered regions | 柔軟領域 | 構造によるが多いものを選択 |

### Special Cases

エッジケースのテスト:

| Case | Why | Examples |
|------|-----|----------|
| Non-standard residues | MSE, HYP等の処理 | 多数のPDBに含まれる |
| Alternate conformations | altloc処理 | 高分解能構造 |
| Missing atoms | 不完全構造 | NMR構造等 |
| Large B-factors | 柔軟領域 | B > 80 Å² |
| Hydrogen atoms | H原子の扱い | NMR構造、HADDOCK等 |

### Radius Classifier Testing

半径分類器のテストに必要な多様性:

| Requirement | Examples |
|-------------|----------|
| All 20 standard amino acids | ほとんどの構造でカバー |
| Modified residues (MSE, HYP, etc.) | 1A0Qなど |
| N/C-terminal variants | すべての構造 |
| Nucleic acids (if supported) | 1BNA, 4TNA |

### Recommended Minimum Set

最小限のベンチマークセット（6構造）:

```
benchmarks/
├── tiny/
│   └── 1crn.cif.gz      # 327 atoms - 正確性検証
├── small/
│   └── 1ubq.cif.gz      # 660 atoms - 基本性能
├── medium/
│   ├── 1a0q.cif.gz      # 3,183 atoms - 標準ベンチ
│   └── 3hhb.cif.gz      # 4,779 atoms - 複数chain
├── large/
│   └── 1aon.cif.gz      # 58,674 atoms - スケール
└── xlarge/
    └── 4v6x.cif.gz      # ~150,000 atoms - 極限
```

### Extended Set (Optional)

より詳細なベンチマーク用:

| PDB | Atoms | Purpose |
|-----|-------|---------|
| 2GB1 | 855 | Small protein, NMR structure |
| 1HRC | 8,619 | Medium-large transition |
| 1OCC | ~7,000 | Membrane protein |
| 1BNA | 486 | DNA (B-form) |
| 1EHZ | ~3,000 | Protein-RNA complex |
| 6LU7 | ~2,400 | SARS-CoV-2 Mpro (recent) |

### Phase 8.2: Reference Data Generation

- [ ] FreeSASAでの参照SASA計算
- [ ] パラメータごとの参照値（n_points: 50, 100, 200, 1000）
- [ ] プローブ半径ごとの参照値（1.2, 1.4, 1.5）

**参照データ形式:**
```json
{
  "structure": "1A0Q",
  "parameters": {
    "n_points": 100,
    "probe_radius": 1.4
  },
  "freesasa_version": "2.1.2",
  "results": {
    "total_area": 18923.28,
    "atom_areas": [...]
  }
}
```

### Phase 8.3: Validation Script

- [ ] Zig実装とFreeSASAの差分計算
- [ ] 許容誤差の設定（2%以内など）
- [ ] 全構造での自動検証

## Files

| File | Action |
|------|--------|
| `benchmarks/structures/` | CREATE (directory) |
| `benchmarks/references/` | CREATE (directory) |
| `benchmarks/inputs/` | CREATE (directory) |
| `scripts/generate_benchmark_data.py` | CREATE |
| `scripts/validate_accuracy.py` | CREATE |

## Directory Structure

```
benchmarks/
├── structures/
│   ├── 1crn.cif.gz
│   ├── 1ubq.cif.gz
│   ├── 1a0q.cif.gz
│   └── ...
├── inputs/
│   ├── 1crn.json
│   ├── 1ubq.json
│   └── ...
└── references/
    ├── 1crn_n100_p1.4.json
    ├── 1a0q_n100_p1.4.json
    └── ...
```

## Success Criteria

- [ ] 6種類以上の構造が用意されている
- [ ] 各構造に対するFreeSASA参照値がある
- [ ] 自動検証スクリプトが動作する

---
- [ ] **DONE** - Phase 8 complete

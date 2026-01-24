# Phase 18: Test Pipeline

## Goal

全件テストのための一貫したパイプラインを構築する。計算・統計・グラフ作成を自動化。

## 優先度: 高

## Tasks

### Phase 18.1: Test Dataset

- [ ] テスト構造セット定義
- [ ] PDB/mmCIF自動ダウンロード
- [ ] サイズ別カテゴリ分け

**テスト構造:**
| Category | Examples | Atoms |
|----------|----------|-------|
| Small | 1CRN, 1L2Y | < 500 |
| Medium | 1UBQ, 1AKE | 500-2000 |
| Large | 1AON, 3J3Q | > 5000 |
| Membrane | 7KQE | - |
| Nucleic | 1BNA | - |

### Phase 18.2: Reference Calculation

- [ ] FreeSASA (Python) での基準値計算
- [ ] 結果保存（JSON形式）
- [ ] バージョン管理

**ディレクトリ構造:**
```
benchmark/
├── structures/
│   ├── 1crn.cif
│   └── ...
├── reference/
│   ├── 1crn_freesasa.json
│   └── ...
└── results/
    ├── 1crn_zig.json
    └── ...
```

### Phase 18.3: Comparison Script

- [ ] Zig実装 vs FreeSASA比較
- [ ] 差分統計（平均、標準偏差、最大）
- [ ] 原子別・残基別・全体

**統計出力:**
```
Structure: 1CRN
Atoms: 327
Total SASA:
  FreeSASA: 4521.34 Å²
  Zig:      4532.18 Å²
  Diff:     +0.24%

Per-atom difference:
  Mean:   0.15%
  Std:    0.42%
  Max:    2.31% (atom 145)
```

### Phase 18.4: Visualization

- [ ] 差分ヒストグラム
- [ ] 構造サイズ vs 計算時間
- [ ] Correlation plot

**ツール:** matplotlib / plotly

**グラフ例:**
- `diff_histogram.png` - 原子別差分分布
- `time_vs_atoms.png` - スケーラビリティ
- `correlation.png` - FreeSASA vs Zig scatter

### Phase 18.5: CI Integration

- [ ] GitHub Actions での自動実行
- [ ] 結果サマリーをPR commentに投稿
- [ ] 閾値超過で警告

## Files

| File | Action |
|------|--------|
| `benchmark/run_all.py` | CREATE - メインスクリプト |
| `benchmark/download.py` | CREATE - 構造ダウンロード |
| `benchmark/compare.py` | CREATE - 比較・統計 |
| `benchmark/visualize.py` | CREATE - グラフ作成 |
| `benchmark/structures.json` | CREATE - テスト構造リスト |
| `.github/workflows/benchmark.yml` | CREATE - CI設定 |

## Success Criteria

- [ ] ワンコマンドで全件テスト実行
- [ ] 統計レポート自動生成
- [ ] グラフが見やすい
- [ ] CIで自動実行される

---
- [ ] **DONE** - Phase 18 complete

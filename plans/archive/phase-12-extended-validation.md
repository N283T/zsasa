# Phase 12: Extended Validation

## Goal

より多くの構造での検証と、他ツールとの比較を行い、実装の正確性を確認する。

## 優先度: 低

## Tasks

### Phase 12.1: Multi-Structure Validation

- [ ] 異なるサイズの構造での検証
- [ ] 膜タンパク質での検証
- [ ] 核酸（DNA/RNA）での検証
- [ ] 糖鎖での検証

**検証構造リスト:**
| PDB | Type | Description |
|-----|------|-------------|
| 1CRN | Soluble | Small protein (crambin) |
| 1A0Q | Soluble | Medium protein |
| 7KQE | Membrane | GPCR |
| 1BNA | Nucleic | DNA double helix |
| 3HYD | Glycan | Glycoprotein |

### Phase 12.2: Comparison with Other Tools

- [ ] NACCESS との比較
- [ ] DSSP との比較（relative ASA）
- [ ] PyMOL の get_area との比較

**比較対象:**
| Tool | Type | Notes |
|------|------|-------|
| FreeSASA | Reference | 主要比較対象 |
| NACCESS | Classic | 古典的ツール |
| DSSP | Secondary | 二次構造と共に計算 |
| PyMOL | Visualization | get_area コマンド |

### Phase 12.3: Edge Case Testing

- [ ] 単一原子
- [ ] 2原子（完全重複）
- [ ] 非常に大きな構造（> 100,000 atoms）
- [ ] 非標準残基

### Phase 12.4: Statistical Analysis

- [ ] 差分の統計（平均、標準偏差、最大）
- [ ] 原子タイプ別の精度分析
- [ ] 構造サイズと精度の関係

**レポート形式:**
```markdown
## Validation Report

### Summary
| Metric | Value |
|--------|-------|
| Structures tested | 10 |
| Average difference | 1.2% |
| Max difference | 3.5% |

### Per-Structure Results
| PDB | Atoms | FreeSASA | Zig | Diff |
|-----|-------|----------|-----|------|
| 1CRN | 327 | 4521.3 | 4532.1 | 0.24% |
| ... | ... | ... | ... | ... |

### Per-Atom-Type Analysis
| Atom Type | Avg Diff | Max Diff |
|-----------|----------|----------|
| C | 0.8% | 2.1% |
| N | 1.1% | 2.8% |
| O | 0.9% | 2.5% |
```

## Files

| File | Action |
|------|--------|
| `scripts/validate_multi.py` | CREATE |
| `scripts/compare_tools.py` | CREATE |
| `docs/validation_report.md` | CREATE |

## Success Criteria

- [ ] 10構造以上で検証完了
- [ ] FreeSASAとの差分が全て5%以内
- [ ] 他ツールとの比較データがある

---
- [ ] **DONE** - Phase 12 complete

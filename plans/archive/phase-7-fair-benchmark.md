# Phase 7: Fair Benchmark Comparison

## Goal

FreeSASAとの公平な比較のため、パーサー部分を統一し、純粋なアルゴリズム性能を測定する。

## 優先度: 高

## 現状の問題

現在のベンチマークは以下を含む:
1. JSONパース時間
2. SASA計算時間
3. 結果出力時間

FreeSASA (Python)との比較では条件が異なる:
- FreeSASA: mmCIF/PDBパース + SASA計算
- Zig実装: JSONパース + SASA計算

## Tasks

### Phase 7.1: Benchmark Isolation

- [ ] パーサー時間を分離計測
- [ ] SASA計算時間のみを計測するオプション追加
- [ ] メモリ使用量の計測

**計測項目:**
```
1. Input parsing time
2. Neighbor list construction time
3. SASA calculation time
4. Output writing time
Total time
Peak memory usage
```

### Phase 7.2: Unified Benchmark Script

- [ ] 同一入力データでの比較スクリプト作成
- [ ] FreeSASA (C) との比較追加
- [ ] 複数構造での自動ベンチマーク

**比較対象:**
| Tool | Language | Parser |
|------|----------|--------|
| FreeSASA CLI | C | mmCIF/PDB |
| FreeSASA Python | Python | mmCIF/PDB |
| freesasa-zig | Zig | JSON (pre-parsed) |
| freesasa-zig | Zig | mmCIF (Phase 9実装後) |

### Phase 7.3: Benchmark Dataset

- [ ] 異なるサイズの構造でのベンチマーク
- [ ] 結果の統計処理（平均、標準偏差）
- [ ] Markdown形式でのレポート生成

**データセット案:**
| PDB | Atoms | Description |
|-----|-------|-------------|
| 1CRN | 327 | Small protein (crambin) |
| 1A0Q | 3,183 | Medium protein (current) |
| 1AON | 58,674 | Large complex (GroEL) |
| 3J3Q | 200,000+ | Ribosome |

## Files

| File | Action |
|------|--------|
| `scripts/benchmark_isolated.py` | CREATE |
| `scripts/benchmark_compare.py` | CREATE |
| `src/main.zig` | MODIFY (timing output) |

## Success Criteria

- [ ] アルゴリズム部分のみの性能が明確に計測できる
- [ ] FreeSASA (C)との公平な比較データが得られる
- [ ] 複数サイズでの性能特性が把握できる

---
- [x] **DONE** - Phase 7 complete

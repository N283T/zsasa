# Phase 15: Chain Support

## Goal

mmCIF入力でのチェーン指定とチェーン別結果表示を実装する。

## 優先度: 高

## Tasks

### Phase 15.1: Chain Filtering

- [ ] `--chain=<id>` オプション追加
- [ ] label_asym_id / auth_asym_id 両対応
- [ ] 複数チェーン指定 (`--chain=A,B`)

**mmCIF カラム:**
| カラム | 説明 |
|--------|------|
| `label_asym_id` | mmCIF標準のチェーンID |
| `auth_asym_id` | PDB著者指定のチェーンID |

**CLI例:**
```bash
# 単一チェーン
freesasa-zig input.cif --chain=A

# 複数チェーン
freesasa-zig input.cif --chain=A,B

# auth_asym_id使用（デフォルトはlabel）
freesasa-zig input.cif --chain=A --auth-chain
```

### Phase 15.2: Per-Chain Results

- [ ] チェーン別SASA集計
- [ ] JSON出力にチェーン情報追加
- [ ] サマリー表示

**出力例:**
```
Total SASA: 12345.67 Å²

Per-chain:
  Chain A:  6543.21 Å² (234 atoms)
  Chain B:  5802.46 Å² (198 atoms)
```

### Phase 15.3: Model Selection

- [ ] `--model=<n>` オプション追加
- [ ] NMR構造対応（複数モデル）
- [ ] デフォルト: 最初のモデル

**mmCIF カラム:** `pdbx_PDB_model_num`

## Files

| File | Action |
|------|--------|
| `src/mmcif_parser.zig` | MODIFY - chain/model columns追加 |
| `src/main.zig` | MODIFY - CLI options追加 |
| `src/types.zig` | MODIFY - chain_id追加 |

## Success Criteria

- [ ] チェーン指定で原子フィルタリングできる
- [ ] チェーン別の面積が表示される
- [ ] 複数モデル構造で特定モデルを選択できる

---
- [ ] **DONE** - Phase 15 complete

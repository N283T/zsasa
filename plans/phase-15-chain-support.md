# Phase 15: Chain Support

## Goal

mmCIF入力でのチェーン指定とチェーン別結果表示を実装する。

## 優先度: 高

## Tasks

### Phase 15.1: Chain Filtering

- [x] `--chain=<id>` オプション追加
- [x] label_asym_id / auth_asym_id 両対応
- [x] 複数チェーン指定 (`--chain=A,B`)

**mmCIF カラム:**
| カラム | 説明 |
|--------|------|
| `label_asym_id` | mmCIF標準のチェーンID（デフォルト） |
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

- [x] チェーン別SASA集計
- [x] サマリー表示
- [ ] JSON出力にチェーン情報追加（Phase 17で対応予定）

**出力例:**
```
Calculated SASA for 41 atoms
Total area: 843.53 Å²

Per-chain SASA:
  A: 744.02 Å² (39 atoms)
  B: 99.51 Å² (2 atoms)
```

### Phase 15.3: Model Selection

- [x] `--model=<n>` オプション追加
- [x] NMR構造対応（複数モデル）
- [x] モデル番号バリデーション（>= 1）

**mmCIF カラム:** `pdbx_PDB_model_num`

## Files

| File | Action |
|------|--------|
| `src/mmcif_parser.zig` | MODIFY - chain/model columns追加 |
| `src/main.zig` | MODIFY - CLI options追加 |
| `src/types.zig` | MODIFY - chain_id追加 |

## Success Criteria

- [x] チェーン指定で原子フィルタリングできる
- [x] チェーン別の面積が表示される
- [x] 複数モデル構造で特定モデルを選択できる

---
- [x] **DONE** - Phase 15 complete

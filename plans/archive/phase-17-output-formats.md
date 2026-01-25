# Phase 17: Output Formats

## Goal

CSV形式での出力をサポートし、スプレッドシートやデータ分析ツールとの連携を容易にする。

## 優先度: 中

## Tasks

### Phase 17.1: CSV Output

- [x] `--format=csv` オプション追加
- [x] ヘッダー行付きCSV出力
- [x] 原子レベル/残基レベル選択（入力形式により自動判定）

**基本CSV出力（JSON入力時）:**
```csv
atom_index,area
0,20.778180
1,13.283910
total,2975.651128
```

**リッチCSV出力（mmCIF入力時）:**
```csv
chain,residue,resnum,atom_name,x,y,z,radius,area
A,THR,1,N,17.047,14.099,3.625,1.550,20.778180
A,THR,1,CA,16.967,12.784,4.338,1.700,13.283910
```

### Phase 17.2: Output File Option

- [x] `-o/--output=FILE` オプション
- [x] 位置引数との併用対応

**CLI例:**
```bash
# CSV to file (positional)
freesasa-zig input.cif result.csv --format=csv

# CSV to file (-o option)
freesasa-zig input.cif --format=csv -o result.csv
```

## Files

| File | Action |
|------|--------|
| `src/json_writer.zig` | MODIFY - Rich CSV formatter added |
| `src/main.zig` | MODIFY - -o option added |

## Success Criteria

- [x] CSV出力がExcel/Googleスプレッドシートで開ける
- [x] ヘッダーが適切
- [x] 数値精度が十分（小数点以下6桁）

---
- [x] **DONE** - Phase 17 complete

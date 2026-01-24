# Phase 17: Output Formats

## Goal

CSV/TSV形式での出力をサポートし、スプレッドシートやデータ分析ツールとの連携を容易にする。

## 優先度: 中

## Tasks

### Phase 17.1: CSV Output

- [ ] `--format=csv` オプション追加
- [ ] ヘッダー行付きCSV出力
- [ ] 原子レベル/残基レベル選択

**原子レベル出力:**
```csv
atom_id,chain,residue,resnum,atom_name,element,x,y,z,radius,sasa
1,A,ALA,1,N,N,10.000,20.000,30.000,1.55,15.23
2,A,ALA,1,CA,C,11.000,21.000,31.000,1.70,8.45
```

**残基レベル出力 (`--per-residue`):**
```csv
chain,residue,resnum,sasa,rsa,state
A,ALA,1,85.23,75.4,Exposed
A,VAL,2,23.45,18.2,Buried
```

### Phase 17.2: TSV Output

- [ ] `--format=tsv` オプション
- [ ] タブ区切り（スプレッドシート貼り付け用）

### Phase 17.3: Output File

- [ ] `-o/--output=<file>` オプション
- [ ] stdout vs file出力切り替え
- [ ] 拡張子による形式自動判定

**CLI例:**
```bash
# CSV to stdout
freesasa-zig input.cif --format=csv

# CSV to file
freesasa-zig input.cif -o result.csv

# TSV with per-residue
freesasa-zig input.cif --format=tsv --per-residue -o result.tsv
```

## Files

| File | Action |
|------|--------|
| `src/output.zig` | CREATE - CSV/TSV formatter |
| `src/main.zig` | MODIFY - output options |

## Success Criteria

- [ ] CSV出力がExcel/Googleスプレッドシートで開ける
- [ ] ヘッダーが適切
- [ ] 数値精度が十分（小数点以下3桁）

---
- [ ] **DONE** - Phase 17 complete

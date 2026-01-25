# Future: Batch Parallel Processing

## Status: Idea / Low Priority

## Background

RustSASAは「FreeSASAより5x高速」と謳っているが、これはバッチ処理（E. coli proteome 4,391構造）での話。
単一タンパク質では RustSASA = FreeSASA（両方4.0ms）。

freesasa-zigの強みは**シングルでも速い**こと：
- vs RustSASA: 1.7x-3.2x 高速
- vs FreeSASA C: 2.0x-4.7x 高速

## Options

### 1. External Tools (Already Works)

```bash
# GNU Parallel
find *.cif.gz | parallel ./freesasa_zig {} {.}_sasa.json

# xargs
ls *.cif | xargs -P 8 -I {} ./freesasa_zig {} {}.json
```

Pros: 実装不要、十分高速
Cons: 外部依存

### 2. Built-in Batch Mode

```bash
freesasa_zig --batch input_dir/ output_dir/ --format json
```

Pros: 統合されたUX、進捗表示可能
Cons: 実装コスト

### 3. Python multiprocessing

```python
from multiprocessing import Pool
from freesasa_zig import calculate_sasa

with Pool(8) as p:
    results = p.map(process_structure, pdb_files)
```

Pros: 柔軟、既存バインディング活用
Cons: Python依存

## Recommendation

現状は外部ツール（GNU Parallel）で十分。
需要があれば`--batch`モードを検討。

## Key Message

> freesasa-zig is fast at the algorithm level, not just batch processing.
> Single-protein SASA: 2-3x faster than RustSASA/FreeSASA.

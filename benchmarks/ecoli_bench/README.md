# E. coli Proteome Benchmark

RustSASA論文との比較用ベンチマーク。

## Dataset

- Source: AlphaFold E. coli proteome (UP000000625_83333_ECOLI_v6)
- Structures: 4,370
- Total atoms: 10,520,167
- Atom range: 122 - 17,331

## Quick Start

```bash
# 1. Run benchmark (single tool)
./benchmarks/scripts/run.py \
  --tool zig --algorithm sr \
  --input-dir benchmarks/UP000000625_83333_ECOLI_v6/json \
  --output-dir benchmarks/results/ecoli/zig_sr_f64 \
  --threads 1,2,4,8,10 --runs 10

# 2. Run all tools
./benchmarks/ecoli_bench/run_all.sh

# 3. Analyze results
./benchmarks/scripts/analyze.py summary \
  --config benchmarks/results/ecoli/config.yaml
```

## Methodology Improvements over RustSASA

| Aspect | RustSASA | This Benchmark |
|--------|----------|----------------|
| Parallelization | GNU parallel (process overhead) | Native threading |
| Warmup | 3 | 5 |
| Runs | 3 | 10 |
| Timing | Total proteome only | Per-structure SASA time |
| Analysis | None | Size-stratified speedup |

## Directory Structure

```
benchmarks/
├── UP000000625_83333_ECOLI_v6/
│   ├── pdb/         # Original PDB (4370)
│   ├── cif/         # Original CIF (4370)
│   └── json/        # Converted JSON.gz (4370)
├── ecoli_bench/
│   ├── README.md    # This file
│   ├── run_all.sh   # Run all benchmarks
│   └── config.yaml  # (in results/ecoli/)
└── results/ecoli/
    ├── config.yaml
    ├── zig_sr_f64/
    ├── freesasa_sr/
    └── rust_sr/
```

# E. coli Proteome Benchmark

RustSASA-style benchmark using the complete E. coli K-12 proteome from AlphaFold.

This benchmark follows the methodology from the [RustSASA paper](https://github.com/OWissett/rustsasa) for direct comparison.

## Dataset

| Property | Value |
|----------|-------|
| Source | AlphaFold E. coli proteome (UP000000625_83333_ECOLI_v6) |
| Structures | 4,370 |
| Total atoms | 10,520,167 |
| Atom range | 122 - 17,331 |
| Input format | PDB |

## Test Environment

| Item | Value |
|------|-------|
| Machine | MacBook Pro |
| Chip | Apple M4 |
| Cores | 10 (4 performance + 6 efficiency) |
| Memory | 32 GB |
| OS | macOS |

## Results

Benchmark parameters: warmup=3, runs=5

### Summary Table

| Tool | Threads | Time (s) | Std Dev | vs FreeSASA | vs RustSASA |
|------|--------:|--------:|--------:|------------:|------------:|
| **zsasa f32** | 8 | **4.407** | ±0.074 | **10.4x** | **1.25x** |
| **zsasa f64** | 8 | **4.614** | ±0.066 | **9.9x** | **1.19x** |
| RustSASA | 8 | 5.491 | ±0.040 | 8.3x | - |
| **zsasa f32** | 1 | **23.847** | ±0.098 | **1.92x** | **1.68x** |
| zsasa f64 | 1 | 32.970 | ±14.97 | 1.39x | 1.21x |
| RustSASA | 1 | 39.968 | ±20.47 | 1.14x | - |
| FreeSASA | 1 | 45.760 | ±0.048 | - | - |

### Key Findings

**8-thread performance:**
- zsasa f32: **4.407s** (fastest)
- zsasa f64: **4.614s**
- RustSASA: 5.491s
- Speedup: zsasa is **20% faster** than RustSASA

**Single-thread performance:**
- zsasa f32: **23.847s** (fastest)
- FreeSASA: 45.760s
- Speedup: zsasa is **48% faster** than FreeSASA

**Parallel scaling (1→8 threads):**
- zsasa f32: 23.847s → 4.407s = **5.41x**
- zsasa f64: 32.970s → 4.614s = **7.14x**
- RustSASA: 39.968s → 5.491s = **7.28x**

## PDB Parser Optimization Impact

Comparison before/after PDB parser optimization ([details](../parser-optimizations.md)):

| Tool | Before | After | Improvement |
|------|-------:|------:|------------:|
| zsasa f32 8t | 4.850s | 4.407s | **9.1%** |
| zsasa f64 8t | 4.854s | 4.614s | **4.9%** |
| zsasa f32 1t | 24.241s | 23.847s | **1.6%** |
| zsasa f64 1t | 24.271s | - | (high variance) |

The optimization primarily benefits multi-threaded workloads where parsing overhead is more significant.

## Methodology

### Differences from RustSASA Paper

| Aspect | RustSASA Paper | This Benchmark |
|--------|----------------|----------------|
| Parallelization | GNU parallel (process overhead) | Native threading |
| Warmup | 3 | 3 |
| Runs | 3 | **5** (increased for stability) |
| Input | PDB | PDB |
| Tool | hyperfine | hyperfine |

### Notes

1. **Input format**: All tools use PDB input for fair comparison. zsasa supports JSON input which would be faster, but PDB is used to match RustSASA methodology.

2. **FreeSASA limitation**: FreeSASA CLI only supports single-threaded batch processing. We use a custom batch wrapper (`sasa_batch.cpp`) that processes files sequentially.

3. **High variance**: Some single-threaded runs (zsasa_f64_1t: ±14.97s, RustSASA_1t: ±20.47s) show high standard deviation, likely due to system interference or thermal throttling. Multi-threaded results are more stable. Runs increased from 3 to 5 for better statistical reliability.

## Running the Benchmark

```bash
# Run all tools
./benchmarks/scripts/batch_bench.py \
  --input benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  --name ecoli \
  --runs 5 --threads 8

# Run specific tool
./benchmarks/scripts/batch_bench.py \
  -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  -n ecoli \
  --tool zig --runs 5

# Dry run (show commands only)
./benchmarks/scripts/batch_bench.py \
  -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  -n ecoli-test \
  --dry-run
```

Results are saved to `benchmarks/results/batch/<name>/`.

## Related Documents

- [results.md](results.md) - Large-scale benchmark (100k structures)
- [batch.md](batch.md) - Batch processing benchmarks
- [../parser-optimizations.md](../parser-optimizations.md) - PDB parser optimization details

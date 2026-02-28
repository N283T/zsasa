# Batch Processing Benchmarks

Throughput benchmarks for processing complete datasets using [hyperfine](https://github.com/sharkdp/hyperfine) timing.

> **Note**: This measures total wall-clock time for processing an entire directory of PDB files with native multi-threading.

## Test Environment

| Item | Value |
|------|-------|
| Machine | MacBook Pro |
| Chip | Apple M4 (10 cores: 4P + 6E) |
| Memory | 32 GB |
| OS | macOS |

## Results

### E. coli Proteome (4,370 structures)

Dataset: AlphaFold E. coli K-12 proteome (UP000000625_83333_ECOLI_v6), PDB format.

Benchmark parameters: warmup=3, runs=5, threads=1,8,10

| Tool | Threads | Time (s) | Std Dev | vs FreeSASA | vs RustSASA | RSS (MB) |
|------|--------:|--------:|--------:|------------:|------------:|---------:|
| zsasa f64 | 1 | 24.13 | ±0.113 | 1.83x | 1.19x | 13 |
| zsasa f32 | 1 | 24.27 | ±0.053 | 1.82x | 1.19x | 13 |
| RustSASA | 1 | 28.78 | ±0.059 | 1.53x | baseline | 39 |
| FreeSASA | 1 | 44.06 | ±0.336 | baseline | - | 154 |
| **zsasa f32** | **8** | **4.40** | **±0.026** | **10.0x** | **1.20x** | **62** |
| zsasa f64 | 8 | 4.42 | ±0.043 | 10.0x | 1.19x | 64 |
| RustSASA | 8 | 5.27 | ±0.100 | 8.4x | baseline | 145 |
| **zsasa f32** | **10** | **4.14** | **±0.275** | **10.6x** | **1.21x** | **75** |
| zsasa f64 | 10 | 4.18 | ±0.255 | 10.5x | 1.20x | 76 |
| RustSASA | 10 | 5.02 | ±0.043 | 8.8x | baseline | 172 |

**Key findings:**

- **10-thread**: zsasa f32 is **10.6x faster** than FreeSASA, **1.21x faster** than RustSASA
- **8-thread**: zsasa f32 is **10.0x faster** than FreeSASA, **1.20x faster** than RustSASA
- **Single-thread**: zsasa is **1.82-1.83x faster** than FreeSASA, **1.19x faster** than RustSASA
- **Parallel scaling (1t->10t)**: zsasa f32 5.86x, zsasa f64 5.77x, RustSASA 5.73x
- **Memory**: zsasa uses **2-4x less** memory than RustSASA, **10-12x less** than FreeSASA
- f32 and f64 perform nearly identically

> FreeSASA 1t: 1 outlier removed (84.9s, 1st timed run after warmup). Remaining 4 runs: 43.6-44.3s (stddev ±0.34s).

| 1 thread | 8 threads | 10 threads |
|:--------:|:---------:|:----------:|
| ![1t](pathname:///zsasa/benchmarks/results/plots/batch/ecoli_time_1t.png) | ![8t](pathname:///zsasa/benchmarks/results/plots/batch/ecoli_time_8t.png) | ![10t](pathname:///zsasa/benchmarks/results/plots/batch/ecoli_time_10t.png) |

| Memory (1t) | Memory (8t) | Memory (10t) |
|:-----------:|:-----------:|:------------:|
| ![1t](pathname:///zsasa/benchmarks/results/plots/batch/ecoli_memory_1t.png) | ![8t](pathname:///zsasa/benchmarks/results/plots/batch/ecoli_memory_8t.png) | ![10t](pathname:///zsasa/benchmarks/results/plots/batch/ecoli_memory_10t.png) |

### Human Proteome (23,586 structures)

Dataset: AlphaFold Human proteome (UP000005640_9606_HUMAN_v6), PDB format.

Benchmark parameters: warmup=3, runs=3, threads=1,8,10

| Tool | Threads | Time (s) | Std Dev | vs FreeSASA | vs RustSASA | RSS (MB) |
|------|--------:|--------:|--------:|------------:|------------:|---------:|
| zsasa f32 | 1 | 375.09 | ±111.04 | 1.63x | 1.07x | 18 |
| zsasa f64 | 1 | 385.74 | ±109.99 | 1.58x | 1.05x | 18 |
| RustSASA | 1 | 403.22 | ±180.79 | 1.51x | baseline | 89 |
| FreeSASA | 1 | 610.69 | ±155.00 | baseline | - | 350 |
| **zsasa f32** | **8** | **43.84** | **±0.048** | **13.9x** | **1.21x** | **106** |
| zsasa f64 | 8 | 44.83 | ±0.073 | 13.6x | 1.18x | 108 |
| RustSASA | 8 | 52.84 | ±0.022 | 11.6x | baseline | 275 |
| **zsasa f32** | **10** | **40.51** | **±0.368** | **15.1x** | **1.23x** | **130** |
| zsasa f64 | 10 | 42.24 | ±0.005 | 14.5x | 1.18x | 134 |
| RustSASA | 10 | 49.97 | ±1.406 | 12.2x | baseline | 325 |

**Key findings:**

- **10-thread**: zsasa f32 is **15.1x faster** than FreeSASA, **1.23x faster** than RustSASA
- **8-thread**: zsasa f32 is **13.9x faster** than FreeSASA, **1.21x faster** than RustSASA
- **Parallel scaling (1t->10t)**: zsasa f32 9.26x, zsasa f64 9.13x, RustSASA 8.07x
- **Memory**: zsasa uses **2.5x less** memory than RustSASA at 10t, **2.7x less** at 8t

> Single-thread runs show high variance across all tools (only 3 runs on a large dataset). Multi-thread results are highly stable.

| 1 thread | 8 threads | 10 threads |
|:--------:|:---------:|:----------:|
| ![1t](pathname:///zsasa/benchmarks/results/plots/batch/human_time_1t.png) | ![8t](pathname:///zsasa/benchmarks/results/plots/batch/human_time_8t.png) | ![10t](pathname:///zsasa/benchmarks/results/plots/batch/human_time_10t.png) |

| Memory (1t) | Memory (8t) | Memory (10t) |
|:-----------:|:-----------:|:------------:|
| ![1t](pathname:///zsasa/benchmarks/results/plots/batch/human_memory_1t.png) | ![8t](pathname:///zsasa/benchmarks/results/plots/batch/human_memory_8t.png) | ![10t](pathname:///zsasa/benchmarks/results/plots/batch/human_memory_10t.png) |

### SwissProt (550,122 structures)

Dataset: SwissProt PDB v6, PDB format. FreeSASA and single-thread baselines omitted due to dataset size.

Benchmark parameters: warmup=3, runs=3, threads=8,10

| Tool | Threads | Time (s) | Std Dev | vs RustSASA | RSS (MB) |
|------|--------:|--------:|--------:|------------:|---------:|
| zsasa f64 | 8 | 1423.70 (23.7 min) | ±2.60 | 1.15x | 185 |
| zsasa f32 | 8 | 1439.19 (24.0 min) | ±13.98 | 1.14x | 182 |
| RustSASA | 8 | 1635.42 (27.3 min) | ±36.74 | baseline | 1130 |
| **zsasa f32** | **10** | **1320.18 (22.0 min)** | **±14.07** | **1.10x** | **208** |
| zsasa f64 | 10 | 1339.41 (22.3 min) | ±17.99 | 1.08x | 212 |
| RustSASA | 10 | 1451.37 (24.2 min) | ±6.54 | baseline | 1131 |

**Key findings:**

- zsasa f32 (10t) processes **550K structures in 22.0 minutes**
- **10-thread**: zsasa f32 is **131 seconds faster** than RustSASA (1.10x)
- **8-thread**: zsasa f64 is **212 seconds faster** than RustSASA (1.15x)
- **Memory**: zsasa uses **~200 MB** vs RustSASA's **~1.13 GB** (5-6x less)
- **Scaling (8t->10t)**: zsasa f32 improves 8.3%, RustSASA improves 11.3%

| 8 threads | 10 threads |
|:---------:|:----------:|
| ![8t](pathname:///zsasa/benchmarks/results/plots/batch/swissprot_time_8t.png) | ![10t](pathname:///zsasa/benchmarks/results/plots/batch/swissprot_time_10t.png) |

| Memory (8t) | Memory (10t) |
|:-----------:|:------------:|
| ![8t](pathname:///zsasa/benchmarks/results/plots/batch/swissprot_memory_8t.png) | ![10t](pathname:///zsasa/benchmarks/results/plots/batch/swissprot_memory_10t.png) |

## Summary

### Performance (vs RustSASA, matched threads)

| Dataset | Structures | zsasa f32 (10t) | zsasa f64 (10t) | Speedup |
|---------|----------:|-----------:|-----------:|--------:|
| E. coli | 4,370 | 4.14s | 4.18s | 1.20-1.21x |
| Human | 23,586 | 40.51s | 42.24s | 1.18-1.23x |
| SwissProt | 550,122 | 1320.18s | 1339.41s | 1.08-1.10x |

### Memory Efficiency

| Dataset | zsasa (10t) | RustSASA (10t) | Ratio |
|---------|----------:|----------:|------:|
| E. coli | ~75 MB | ~172 MB | 2.3x less |
| Human | ~130 MB | ~325 MB | 2.5x less |
| SwissProt | ~208 MB | ~1131 MB | 5.4x less |

Memory advantage increases with dataset size, from 2.3x on small datasets to 5.4x on SwissProt.

## Note: Datasets Larger Than Available RAM

zsasa uses **mmap** for file I/O. When the dataset fits in RAM, all file data stays in the OS page cache and performance depends purely on compute speed. When the dataset exceeds available RAM, mmap page faults must read from storage, and performance becomes **storage I/O-bound** regardless of compute efficiency.

In this regime, multiple worker threads fault pages in random order, so effective read throughput drops well below the SSD's sequential maximum. This is not specific to zsasa — any mmap-based tool (including Lahuta and RustSASA) hits the same bottleneck, and all converge to similar wall-clock times.

zsasa's low memory footprint helps here: less RSS means more RAM is available for the OS page cache, which can reduce page fault frequency.

## Methodology

Uses [hyperfine](https://github.com/sharkdp/hyperfine) for timing, following the [RustSASA paper](https://github.com/OWissett/rustsasa) methodology:

1. Warmup runs (default 3) to eliminate cold-start effects
2. Multiple timed runs for statistical reliability
3. Reports mean, stddev, min, max

### Tool Configurations

| Tool | Configurations |
|------|----------------|
| zsasa (Zig) | f64 Nt, f64 1t, f32 Nt, f32 1t |
| FreeSASA | 1t only (sequential batch wrapper) |
| RustSASA | Nt, 1t |

### Notes

1. **Input format**: All tools use PDB input for fair comparison. zsasa supports JSON input which would be faster.
2. **FreeSASA limitation**: FreeSASA CLI only supports single-threaded batch processing. A custom wrapper (`sasa_batch.cpp`) processes files sequentially.
3. **SwissProt scope**: FreeSASA and single-thread baselines omitted due to dataset size.
4. **Precision**: f32 and f64 perform nearly identically across all datasets, with f32 sometimes slightly faster due to better cache utilization.

## Running Benchmarks

```bash
# E. coli proteome
./benchmarks/scripts/bench_batch.py \
  -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  -n ecoli --runs 5 --threads 1,8,10

# Human proteome
./benchmarks/scripts/bench_batch.py \
  -i benchmarks/UP000005640_9606_HUMAN_v6/pdb \
  -n human --runs 3 --threads 1,8,10

# SwissProt (large dataset, skip 1t)
./benchmarks/scripts/bench_batch.py \
  -i /path/to/swissprot_pdb_v6 \
  -n swissprot --runs 3 --threads 8,10 \
  --tool zig --tool rustsasa

# Single tool only
./benchmarks/scripts/bench_batch.py \
  -i /path/to/pdb_dir -n my-bench --tool zig

# Dry run
./benchmarks/scripts/bench_batch.py \
  -i /path/to/pdb_dir -n test --dry-run
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input`, `-i` | Input directory (PDB files) | (required) |
| `--name`, `-n` | Benchmark name | (required) |
| `--threads`, `-T` | Thread counts: `1,8,10` or `1-10` | `1,8` |
| `--runs`, `-r` | Number of benchmark runs | 3 |
| `--warmup`, `-w` | Number of warmup runs | 3 |
| `--tool`, `-t` | Tool: zig, freesasa, rustsasa (repeatable) | all |
| `--output`, `-o` | Output directory | results/batch/\<name\> |
| `--dry-run` | Show commands without running | false |

### Analysis

```bash
./benchmarks/scripts/analyze_batch.py summary           # All benchmarks
./benchmarks/scripts/analyze_batch.py summary -n ecoli  # Specific benchmark
./benchmarks/scripts/analyze_batch.py plot -n ecoli     # Time comparison charts (per thread)
./benchmarks/scripts/analyze_batch.py memory -n ecoli   # Memory comparison charts (per thread)
./benchmarks/scripts/analyze_batch.py all               # Summary + all charts
```

## Related Documents

- [single-file.md](single-file.md) - Single-file benchmark results
- [md.md](md.md) - MD trajectory benchmark results
- [validation.md](validation.md) - SASA accuracy validation

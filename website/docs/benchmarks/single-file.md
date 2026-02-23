# Single-File SASA Benchmarks

Large-scale benchmark results for zsasa using Shrake-Rupley algorithm.

- **Dataset**: 2,006 structures (stratified sampling from PDB + AlphaFold)
- **Precision**: Zig/FreeSASA use f64, RustSASA uses f32

> **Note**:
> - Absolute execution times are environment-dependent. **Relative speedup ratios** are the meaningful metric for comparison.
> - All implementations use identical parameters: `n_points=100`, `probe_radius=1.4Å`
> - Benchmarks measure **SASA calculation time only** (file I/O excluded). See [Methodology](#methodology) for details.
> - SASA accuracy validated: mean error <0.001% vs FreeSASA reference. See [validation.md](validation.md) for details.

## Highlights

Zig's key advantage: **Large structures + Multi-threading**

| Speedup at threads=10 | Thread Scaling (100k+ atoms) |
|:---------------:|:----------------------------:|
| ![Speedup](pathname:///zsasa/benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](pathname:///zsasa/benchmarks/results/plots/large/speedup_by_threads.png) |

**Key Results (100k+ atoms, n=506):**
- **Up to 2.93x faster** than FreeSASA (8to0: 673,884 atoms, threads=10)
- **2.3x** median speedup vs FreeSASA and RustSASA (threads=10)
- Speedup increases with thread count (parallel efficiency advantage)

---

## Methodology

### SASA-Only Timing

For fair comparison, we measure **SASA calculation time only**. File I/O is excluded.

```
Total time = File I/O + SASA calculation + Output
                        ^^^^^^^^^^^^^^^^
                        Only this is measured
```

Measurement method for each implementation:

| Implementation | Method |
|----------------|--------|
| zsasa | Internal measurement via `--timing` option (stderr output) |
| FreeSASA C | Patched binary outputs SASA calculation time to stderr |
| RustSASA | Patched binary outputs `SASA_TIME_US` to stderr |

### Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Algorithm | Shrake-Rupley | Supported by all implementations (LR: Zig/FreeSASA only) |
| n_points | 100 | Number of test points |
| probe_radius | 1.4 A | Water molecule radius |
| Runs | 3 | Average value used |

### Stratified Sampling

Stratified sampling from PDB and AlphaFold structures, 150 structures per bin (13 bins):

| Bin | Atom Range | Count |
|-----|-----------|------:|
| 0-500 | 0 - 500 | 150 |
| 500-1k | 500 - 1,000 | 150 |
| 1k-2k | 1,000 - 2,000 | 150 |
| 2k-3k | 2,000 - 3,000 | 150 |
| 3k-5k | 3,000 - 5,000 | 150 |
| 5k-10k | 5,000 - 10,000 | 150 |
| 10k-20k | 10,000 - 20,000 | 150 |
| 20k-50k | 20,000 - 50,000 | 150 |
| 50k-75k | 50,000 - 75,000 | 150 |
| 75k-100k | 75,000 - 100,000 | 150 |
| 100k-150k | 100,000 - 150,000 | 150 |
| 150k-200k | 150,000 - 200,000 | 150 |
| 200k-500k | 200,000 - 500,000 | 150 |
| 500k+ | 500,000+ | 56 |
| **Total** | | **2,006** |

Same seed produces identical samples (reproducible):

```bash
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o a.json
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o b.json
diff a.json b.json  # No differences
```

---

## Test Environment

| Item | Value |
|------|-------|
| Machine | MacBook Pro |
| Chip | Apple M4 |
| Cores | 10 (4 performance + 6 efficiency) |
| Memory | 32 GB |
| OS | macOS |

## Executive Summary

| Metric | Zig vs FreeSASA | Zig vs RustSASA |
| --- | --- | --- |
| **Best case (50k+ atoms)** | **2.93x** (8to0) | **2.60x** (8bsj) |
| **Overall (threads=10)** | **2.08x** median | **2.22x** median |
| **Large structures (100k+)** | **2.28x** | **2.32x** |
| **Largest structure (4.5M atoms)** | **2.5x** | **2.3x** |
| **Parallel efficiency (threads=10)** | **+43%** | **+112%** |

---

## Dataset

| Size Bin | Atom Range | Count |
| --- | ---: | ---: |
| 0-500 | 0-500 | 150 |
| 500-1k | 500-1,000 | 150 |
| 1k-2k | 1,000-2,000 | 150 |
| 2k-3k | 2,000-3,000 | 150 |
| 3k-5k | 3,000-5,000 | 150 |
| 5k-10k | 5,000-10,000 | 150 |
| 10k-20k | 10,000-20,000 | 150 |
| 20k-50k | 20,000-50,000 | 150 |
| 50k-75k | 50,000-75,000 | 150 |
| 75k-100k | 75,000-100,000 | 150 |
| 100k-150k | 100,000-150,000 | 150 |
| 150k-200k | 150,000-200,000 | 150 |
| 200k-500k | 200,000-500,000 | 150 |
| 500k+ | 500,000+ | 56 |
| **Total** | | **2,006** |

- **Source**: PDB and AlphaFold structures
- **Sampling**: Stratified sampling, 150 per bin (seed 42)
- **Large structures**: 500k+ bin has 56 structures (fewer available)

---

## Single-Thread Performance (threads=1)

Single-threaded comparison (excluding parallelization effects):

| Size Bin | Count | vs FreeSASA | vs RustSASA |
| --- | ---: | ---: | ---: |
| 0-500 | 150 | 1.23x | 0.84x |
| 500-1k | 150 | 1.24x | 0.89x |
| 1k-2k | 150 | 1.28x | 0.93x |
| 2k-3k | 150 | 1.35x | 1.00x |
| 3k-5k | 150 | 1.42x | 1.04x |
| 5k-10k | 150 | 1.49x | 1.05x |
| 10k-20k | 150 | 1.48x | 1.06x |
| 20k-50k | 150 | 1.46x | 1.06x |
| 50k-75k | 150 | 1.51x | 1.06x |
| 75k-100k | 150 | 1.55x | 1.05x |
| 100k-150k | 150 | **1.56x** | 1.07x |
| 150k-200k | 150 | **1.58x** | 1.07x |
| 200k-500k | 150 | **1.58x** | 1.09x |
| 500k+ | 56 | **1.59x** | 1.08x |

**Observations:**
- Zig vs FreeSASA: **1.6x** on large structures (SIMD optimization)
- Zig vs Rust: Nearly equal at 1t (both use SIMD), Zig slightly ahead on large structures

---

## Multi-Thread Performance (threads=10)

### Overall Statistics

| Tool | Structures | Median (ms) | Mean (ms) | P95 (ms) |
| --- | ---: | ---: | ---: | ---: |
| **Zig** | 2,006 | **7.51** | 62.11 | 174.93 |
| FreeSASA | 2,006 | 14.24 | 142.68 | 389.82 |
| RustSASA | 2,006 | 15.77 | 143.39 | 417.80 |

### Speedup by Structure Size

![Speedup by Size and Threads](pathname:///zsasa/benchmarks/results/plots/speedup_by_bin/grid.png)

| Size Bin | Count | vs FreeSASA | vs RustSASA |
| --- | ---: | ---: | ---: |
| 0-500 | 150 | 1.14x | 1.10x |
| 500-1k | 150 | 1.48x | 1.52x |
| 1k-2k | 150 | 1.57x | 1.66x |
| 2k-3k | 150 | 1.68x | 1.86x |
| 3k-5k | 150 | 1.82x | 1.99x |
| 5k-10k | 150 | 2.02x | 2.06x |
| 10k-20k | 150 | 2.06x | 2.17x |
| 20k-50k | 150 | 2.07x | **2.32x** |
| 50k-75k | 150 | 2.18x | **2.30x** |
| 75k-100k | 150 | **2.31x** | **2.26x** |
| 100k-150k | 150 | **2.30x** | **2.30x** |
| 150k-200k | 150 | **2.29x** | **2.31x** |
| 200k-500k | 150 | **2.23x** | **2.34x** |
| 500k+ | 56 | **2.31x** | **2.33x** |

**Observations:**
- **Small (0-500)**: Overhead dominant, modest speedup
- **Medium (1k-20k)**: Stable **1.6x-2.1x** speedup
- **Large (50k+)**: Up to **2.3x** speedup with consistent results (narrow IQR)

**Key Insight:**
- At threads=1, Zig vs Rust is nearly equal
- At threads=10, Zig takes a significant lead -> **parallel efficiency difference**

---

## Thread Scaling

### Median Execution Time by Thread Count

![Thread Scaling](pathname:///zsasa/benchmarks/results/plots/large/speedup_by_threads.png)

| Threads | Zig (ms) | FreeSASA (ms) | Rust (ms) |
| ---: | ---: | ---: | ---: |
| 1 | 20.71 | 30.34 | 22.21 |
| 2 | 12.76 | 20.87 | 18.55 |
| 4 | 9.19 | 15.96 | 16.55 |
| 8 | 7.81 | 14.48 | 15.94 |
| 10 | **7.31** | 14.19 | 15.77 |

**Speedup from threads=1 to threads=10:**
- Zig: 20.71 -> 7.31 = **2.83x**
- FreeSASA: 30.34 -> 14.19 = **2.14x**
- Rust: 22.21 -> 15.77 = **1.41x**

---

## Parallel Efficiency

### Definition

```
Parallel Efficiency = T1 / (TN x N)
```

- T1 = Single-thread execution time
- TN = N-thread execution time
- 1.0 = Ideal linear scaling

### Efficiency by Thread Count

![Parallel Efficiency](pathname:///zsasa/benchmarks/results/plots/efficiency/summary.png)

| Threads | Zig | FreeSASA | Rust | Zig vs FS | Zig vs Rust |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.000 | 1.000 | 1.000 | - | - |
| 2 | 0.816 | 0.724 | 0.596 | **+13%** | **+37%** |
| 4 | 0.597 | 0.466 | 0.333 | **+28%** | **+79%** |
| 8 | 0.350 | 0.256 | 0.173 | **+37%** | **+102%** |
| 10 | 0.295 | 0.207 | 0.139 | **+43%** | **+112%** |

### Efficiency by Size Bin (threads=10)

| Size Bin | Zig | FreeSASA | Rust |
| --- | ---: | ---: | ---: |
| 0-500 | 0.16 | 0.17 | 0.12 |
| 500-1k | 0.24 | 0.20 | 0.14 |
| 1k-2k | 0.26 | 0.21 | 0.14 |
| 2k-3k | 0.27 | 0.21 | 0.14 |
| 3k-5k | 0.28 | 0.21 | 0.14 |
| 5k-10k | 0.29 | 0.21 | 0.14 |
| 10k-20k | 0.30 | 0.22 | 0.14 |
| 20k-50k | 0.30 | 0.21 | 0.14 |
| 50k-75k | 0.30 | 0.21 | 0.14 |
| 75k-100k | 0.30 | 0.20 | 0.14 |
| 100k-150k | **0.30** | 0.20 | 0.14 |
| 150k-200k | **0.30** | 0.21 | 0.14 |
| 200k-500k | 0.29 | 0.21 | 0.14 |
| 500k+ | **0.30** | 0.20 | 0.14 |

**Observations:**
- Zig has the highest parallel efficiency across all sizes
- Efficiency peaks at ~0.30 for large structures
- Rust has low efficiency and benefits less from thread increases

---

## Large Structure Analysis

### Summary (100k+ atoms, n=506)

| Speedup at threads=10 | Thread Scaling |
|:---------------:|:--------------:|
| ![Speedup](pathname:///zsasa/benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](pathname:///zsasa/benchmarks/results/plots/large/speedup_by_threads.png) |

| Comparison | Median Speedup |
| --- | ---: |
| vs FreeSASA | **2.28x** |
| vs RustSASA | **2.32x** |

**Observations:**
- Speedup improves with thread count (1.6x->2.3x vs FreeSASA)
- vs Rust: dramatic improvement (1.1x->2.3x) due to parallel efficiency difference

### Maximum Structure: 9fqr (4,506,416 atoms)

![Max Structure Scaling](pathname:///zsasa/benchmarks/results/plots/samples/max_structure.png)

Thread scaling on the largest PDB structure (9fqr, mean of 3 runs):

| Threads | Zig (s) | FreeSASA (s) | Rust (s) | vs FS | vs Rust |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 8.89 | 14.50 | 9.52 | 1.6x | 1.1x |
| 2 | 5.41 | 10.35 | 8.00 | 1.9x | 1.5x |
| 4 | 3.68 | 8.32 | 7.22 | 2.3x | 2.0x |
| 8 | 3.08 | 7.73 | 6.94 | 2.5x | 2.3x |
| 10 | 2.95 | 7.51 | 6.87 | **2.5x** | **2.3x** |

**Key Insight**:
- Speedup ratio improves with increasing thread count
- Achieves **2.5x** vs FreeSASA, **2.3x** vs Rust at threads=10
- Rust plateaus with increasing threads (parallel efficiency issue)

### Best Speedup Structures (50k+ atoms)

![Speedup Comparison](pathname:///zsasa/benchmarks/results/plots/speedup/comparison.png)

Top 5 structures with highest speedup at any thread count:

| Rank | vs FreeSASA | Atoms | Threads | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 1 | 8to0 | 673,884 | 10 | **2.93x** |
| 2 | 8to0 | 673,884 | 8 | 2.88x |
| 3 | 8rbs | 164,605 | 10 | 2.84x |
| 4 | 8vkq | 150,620 | 10 | 2.70x |
| 5 | 8vkr | 152,898 | 10 | 2.68x |

| Rank | vs Rust | Atoms | Threads | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 1 | 8bsj | 81,467 | 10 | **2.60x** |
| 2 | 8f2n | 162,111 | 10 | 2.50x |
| 3 | 7l5q | 248,040 | 10 | 2.44x |
| 4 | 9n5l | 239,460 | 10 | 2.42x |
| 5 | 9j7k | 243,000 | 10 | 2.42x |

**Note**: threads=8 sometimes outperforms threads=10 (e.g., 8to0: 2.93x at t=10 vs 2.88x at t=8).

---

## Execution Time Distribution

![SR Scatter Plot](pathname:///zsasa/benchmarks/results/plots/scatter/sr/grid.png)

**Observations:**
- Nearly linear on log scale -> O(N) neighbor list is effective (all tools use cell list)
- Zig (green) is consistently lower (faster) across all sizes
- Gap between 3 tools widens with increasing thread count
- Few outliers -> stable performance

---

## Per-Bin Sample Results

Thread scaling details on representative structures selected from each size bin.

| Bin | Atoms Range | Sample Plot |
|-----|-------------|-------------|
| 0-500 | 0-500 | [View](pathname:///zsasa/benchmarks/results/plots/samples/0-500.png) |
| 500-1k | 500-1,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/500-1k.png) |
| 1k-2k | 1,000-2,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/1k-2k.png) |
| 2k-3k | 2,000-3,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/2k-3k.png) |
| 3k-5k | 3,000-5,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/3k-5k.png) |
| 5k-10k | 5,000-10,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/5k-10k.png) |
| 10k-20k | 10,000-20,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/10k-20k.png) |
| 20k-50k | 20,000-50,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/20k-50k.png) |
| 50k-75k | 50,000-75,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/50k-75k.png) |
| 75k-100k | 75,000-100,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/75k-100k.png) |
| 100k-150k | 100,000-150,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/100k-150k.png) |
| 150k-200k | 150,000-200,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/150k-200k.png) |
| 200k-500k | 200,000-500,000 | [View](pathname:///zsasa/benchmarks/results/plots/samples/200k-500k.png) |
| 500k+ | 500,000+ | [View](pathname:///zsasa/benchmarks/results/plots/samples/500kplus.png) |

---

## Key Takeaways

> **Why is Zig faster?** SIMD optimization (8-wide distance calculation), multi-threading with work stealing, and spatial hashing for O(1) neighbor lookup.

1. **Maximum effect on large structures**
   - **2.3x** speedup for 100k+ atoms
   - **2.5x** speedup for the largest structure (4.5M atoms)

2. **Consistent advantage**
   - Outperforms FreeSASA across all sizes (500+ atoms)
   - Fastest at all thread counts

3. **Excellent parallel efficiency**
   - **43%** higher than FreeSASA
   - **112%** higher than RustSASA (at threads=10)

4. **Accurate results**
   - See [validation.md](validation.md) for detailed accuracy analysis

---

## Running Benchmarks

### Setup

```bash
# Build Zig binary
zig build -Doptimize=ReleaseFast

# External tools setup (for comparison)
cd benchmarks/external
git clone https://github.com/N283T/freesasa-bench.git
cd freesasa-bench && ./configure --enable-threads && make && cd ..
git clone --recursive https://github.com/N283T/rustsasa-bench.git
cd rustsasa-bench && cargo build --release --features cli && cd ..
```

### Index & Sample Generation

```bash
# Create index (first time only)
./benchmarks/scripts/build_index.py benchmarks/inputs

# Check distribution
./benchmarks/scripts/sample.py benchmarks/inputs/index.json --analyze

# Generate sample
./benchmarks/scripts/sample.py benchmarks/inputs/index.json \
    --seed 42 \
    -o benchmarks/dataset/sample.json
```

### Running

```bash
# Basic usage
./benchmarks/scripts/bench.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/bench.py --tool freesasa --algorithm lr --threads 1-10

# With sample file
./benchmarks/scripts/bench.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/dataset/sample.json \
    --threads 1,2,4,8,10

# Single run for quick testing
./benchmarks/scripts/bench.py --tool zig --algorithm sr --threads 1 --runs 1

# With f32 precision (Zig only)
./benchmarks/scripts/bench.py --tool zig --algorithm sr --precision f32
```

### Analysis & Visualization

```bash
# Summary tables
./benchmarks/scripts/analyze.py summary

# Generate all plots
./benchmarks/scripts/analyze.py all

# Individual plot types
./benchmarks/scripts/analyze.py scatter      # Atoms vs time scatter
./benchmarks/scripts/analyze.py threads      # Thread scaling
./benchmarks/scripts/analyze.py grid         # Speedup grid by size/threads
./benchmarks/scripts/analyze.py validation   # SASA validation
./benchmarks/scripts/analyze.py samples      # Per-bin sample plots
./benchmarks/scripts/analyze.py large        # Large structure analysis
./benchmarks/scripts/analyze.py efficiency   # Parallel efficiency

# Export to CSV
./benchmarks/scripts/analyze.py export-csv
```

### Scripts

| Script | Purpose |
|--------|---------|
| `build_index.py` | Create atom count index from all input files |
| `sample.py` | Stratified sampling from index |
| `bench.py` | Run benchmarks (single-file mode) |
| `analyze.py` | Analyze results and generate plots |
| `generate_json.py` | Convert CIF/PDB to JSON format |

### Notes

1. **Initial runs are slow**: Due to file cache and warmup effects
2. **Thread count depends on CPU**: Optimal when matching physical core count
3. **External tools require patches**: SASA-only timing requires modified binaries

---

## Appendix: Lee-Richards (LR) Algorithm

Lee-Richards results using the same 2,006 structures. RustSASA does not support LR.

### Overall Statistics (threads=10)

| Tool | Structures | Median (ms) | Mean (ms) | P95 (ms) |
| --- | ---: | ---: | ---: | ---: |
| **Zig** | 2,006 | **20.82** | 178.39 | 510.65 |
| FreeSASA | 2,006 | 39.58 | 327.41 | 916.19 |

**Key Insight**: Zig is **1.79x** faster than FreeSASA (median per-structure speedup)

### Speedup by Structure Size (threads=10)

| Size Bin | Count | vs FreeSASA |
| --- | ---: | ---: |
| 0-500 | 150 | 1.04x |
| 500-1k | 150 | 1.41x |
| 1k-2k | 150 | 1.57x |
| 2k-3k | 150 | 1.71x |
| 3k-5k | 150 | 1.76x |
| 5k-10k | 150 | 1.86x |
| 10k-20k | 150 | 1.93x |
| 20k-50k | 150 | 1.77x |
| 50k-75k | 150 | 1.78x |
| 75k-100k | 150 | **1.82x** |
| 100k-150k | 150 | **1.82x** |
| 150k-200k | 150 | **1.82x** |
| 200k-500k | 150 | **1.80x** |
| 500k+ | 56 | **1.85x** |

### Thread Scaling

| Threads | Zig (ms) | FreeSASA (ms) |
| ---: | ---: | ---: |
| 1 | 107.08 | 172.88 |
| 2 | 58.24 | 94.63 |
| 4 | 33.45 | 55.05 |
| 8 | 22.90 | 40.13 |
| 10 | **20.82** | 39.58 |

**Speedup from threads=1 to threads=10:**
- Zig: 107.08 -> 20.82 = **5.14x**
- FreeSASA: 172.88 -> 39.58 = **4.37x**

### Execution Time Distribution

![LR Scatter Plot](pathname:///zsasa/benchmarks/results/plots/scatter/lr/grid.png)

**Observations:**
- Overall 3-4x slower than SR (slice integration cost)
- Zig's advantage is maintained in LR as well

---

## Related Documents

- [batch.md](batch.md) - Batch processing benchmarks (proteome datasets)
- [md.md](md.md) - MD trajectory benchmarks
- [validation.md](validation.md) - SASA accuracy validation
- [Algorithms](../guide/algorithms.mdx) - Algorithm comparison and usage

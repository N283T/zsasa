# Single-File SASA Benchmarks

Large-scale benchmark results for zsasa using Shrake-Rupley algorithm.

- **Dataset**: ~100k structures (stratified sampling from PDB)
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
| ![Speedup](../../benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](../../benchmarks/results/plots/large/speedup_by_threads.png) |

**Key Results (100k+ atoms, n=1,171):**
- **Up to 3.05x faster** than FreeSASA (8to0: 673,884 atoms, threads=8)
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

Stratified sampling from all PDB structures (~238K):

| Bin | Range | Strategy |
|-----|-------|----------|
| 0-500 | 0 - 500 | Proportional allocation |
| 500-2k | 500 - 2,000 | Proportional allocation |
| 2k-10k | 2,000 - 10,000 | Proportional allocation |
| 10k-50k | 10,000 - 50,000 | Proportional allocation |
| 50k-200k | 50,000 - 200,000 | **All included** |
| 200k+ | 200,000+ | **All included** |

Rare large structures (50k+ atoms) are all included; the rest use proportional allocation.

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
| **Best case (50k+ atoms)** | **3.05x** (8to0) | **2.79x** (5wk6) |
| **Overall (threads=10)** | **1.45x** median | **2.07x** median |
| **Large structures (100k+)** | **2.3x** | **2.3x** |
| **Largest structure (4.5M atoms)** | **2.9x** | **2.2x** |
| **Parallel efficiency (threads=10)** | **+26%** | **+97%** |
| **Instruction count** | **2.4x fewer** | Comparable |

---

## Dataset

![Dataset Distribution](../../benchmarks/results/plots/dataset/stratified_100k.png)

| Sampling Bin | Size Bin | Atoms | Count | Percentage |
| :---: | --- | ---: | ---: | ---: |
| 0-500 | 0-500 | 0-500 | 2,506 | 2.5% |
| 500-2k | 500-1k | 500-1,000 | 5,744 | 5.7% |
| | 1k-2k | 1,000-2,000 | 15,922 | 15.9% |
| 2k-10k | 2k-5k | 2,000-5,000 | 36,123 | 36.1% |
| | 5k-10k | 5,000-10,000 | 19,835 | 19.8% |
| 10k-50k | 10k-20k | 10,000-20,000 | 10,187 | 10.2% |
| | 20k-50k | 20,000-50,000 | 5,377 | 5.4% |
| 50k-200k | 50k-100k | 50,000-100,000 | 3,133 | 3.1% |
| | 100k-200k | 100,000-200,000 | 900 | 0.9% |
| 200k+ | 200k+ | 200,000+ | 271 | 0.3% |
|  | **Total** |  | **99,998** |  |

- **Source**: All structures from the Protein Data Bank (PDB)
- **Sampling**: Stratified sampling (seed 42)
- **Large structures**: 50k+ atoms are all included (due to rarity)

---

## Single-Thread Performance (threads=1)

Single-threaded comparison (excluding parallelization effects):

| Size Bin | vs FreeSASA | vs RustSASA |
| --- | ---: | ---: |
| 0-500 | 0.95x | 0.74x |
| 500-1k | 0.99x | 0.79x |
| 1k-2k | 1.04x | 0.82x |
| 2k-5k | 1.13x | 0.86x |
| 5k-10k | 1.24x | 0.93x |
| 10k-20k | 1.35x | 1.01x |
| 20k-50k | 1.48x | 1.09x |
| 50k-100k | **1.56x** | 1.10x |
| 100k-200k | **1.60x** | 1.11x |
| 200k+ | **1.60x** | 1.13x |

**Observations:**
- Zig vs FreeSASA: **1.6x** on large structures (SIMD optimization)
- Zig vs Rust: Nearly equal (both use SIMD)

---

## Multi-Thread Performance (threads=10)

### Overall Statistics

| Tool | Structures | Median (ms) | Mean (ms) | P95 (ms) |
| --- | ---: | ---: | ---: | ---: |
| **Zig** | 99,998 | **3.24** | 7.78 | 28.91 |
| FreeSASA | 99,998 | 4.69 | 14.53 | 57.06 |
| RustSASA | 99,998 | 5.60 | 15.68 | 60.39 |

### Speedup by Structure Size

![Speedup by Size and Threads](../../benchmarks/results/plots/speedup_by_bin/grid.png)

| Size Bin | Count | vs FreeSASA | vs RustSASA |
| --- | ---: | ---: | ---: |
| 0-500 | 2,506 | 0.92x | 0.97x |
| 500-1k | 5,744 | 1.18x | 1.36x |
| 1k-2k | 15,922 | 1.26x | 1.54x |
| 2k-5k | 36,123 | 1.42x | 1.70x |
| 5k-10k | 19,835 | 1.56x | 1.84x |
| 10k-20k | 10,187 | 1.68x | 1.95x |
| 20k-50k | 5,377 | 1.93x | 2.11x |
| 50k-100k | 3,133 | **2.22x** | **2.25x** |
| 100k-200k | 900 | **2.31x** | **2.30x** |
| 200k+ | 271 | **2.28x** | **2.34x** |

**Observations:**
- **Small (0-500)**: Overhead dominant, no speedup
- **Medium (1k-20k)**: Stable **1.3x-1.9x** speedup
- **Large (50k+)**: Up to **2.3x** speedup with consistent results (narrow IQR)

**Key Insight:**
- At threads=1, Zig vs Rust is nearly equal
- At threads=10, Zig takes a significant lead -> **parallel efficiency difference**

---

## Thread Scaling

### Median Execution Time by Thread Count

![Thread Scaling](../../benchmarks/results/plots/large/speedup_by_threads.png)

| Threads | Zig (ms) | FreeSASA (ms) | Rust (ms) |
| ---: | ---: | ---: | ---: |
| 1 | 8.94 | 10.11 | 7.71 |
| 2 | 5.64 | 6.89 | 6.46 |
| 4 | 4.01 | 5.20 | 5.84 |
| 8 | 3.37 | 4.79 | 5.63 |
| 10 | **3.24** | 4.69 | 5.60 |

**Speedup from threads=1 to threads=10:**
- Zig: 8.94 -> 3.24 = **2.76x**
- FreeSASA: 10.11 -> 4.69 = **2.16x**
- Rust: 7.71 -> 5.60 = **1.38x**

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

![Parallel Efficiency](../../benchmarks/results/plots/efficiency/summary.png)

| Threads | Zig | FreeSASA | Rust | Zig vs FS | Zig vs Rust |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.000 | 1.000 | 1.000 | - | - |
| 2 | 0.791 | 0.734 | 0.595 | **+8%** | **+33%** |
| 4 | 0.556 | 0.486 | 0.331 | **+14%** | **+68%** |
| 8 | 0.326 | 0.263 | 0.171 | **+24%** | **+90%** |
| 10 | 0.271 | 0.214 | 0.138 | **+26%** | **+97%** |

### Efficiency by Size Bin (threads=10)

| Size Bin | Zig | FreeSASA | Rust |
| --- | ---: | ---: | ---: |
| 0-500 | 0.160 | 0.165 | 0.123 |
| 500-1k | 0.241 | 0.201 | 0.139 |
| 1k-2k | 0.261 | 0.214 | 0.139 |
| 2k-5k | 0.274 | 0.216 | 0.138 |
| 5k-10k | 0.274 | 0.216 | 0.137 |
| 10k-20k | 0.269 | 0.214 | 0.137 |
| 20k-50k | 0.273 | 0.209 | 0.139 |
| 50k-100k | 0.288 | 0.202 | 0.141 |
| 100k-200k | **0.291** | 0.200 | 0.141 |
| 200k+ | **0.285** | 0.200 | 0.138 |

**Observations:**
- Zig has the highest parallel efficiency across all sizes
- Efficiency improves for large structures (achieves 0.93)
- Rust has low efficiency and benefits less from thread increases

---

## Large Structure Analysis

### Summary (100k+ atoms, n=1,171)

| Speedup at threads=10 | Thread Scaling |
|:---------------:|:--------------:|
| ![Speedup](../../benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](../../benchmarks/results/plots/large/speedup_by_threads.png) |

| Comparison | Median Speedup |
| --- | ---: |
| vs FreeSASA | **2.3x** |
| vs RustSASA | **2.3x** |

**Observations:**
- Speedup improves with thread count (1.6x->2.3x vs FreeSASA)
- vs Rust: dramatic improvement (1.1x->2.3x) due to parallel efficiency difference

### Maximum Structure: 9fqr (4,506,416 atoms)

![Max Structure Scaling](../../benchmarks/results/plots/samples/max_structure.png)

Thread scaling on the largest PDB structure (9fqr, mean of 3 runs):

| Threads | Zig (s) | FreeSASA (s) | Rust (s) | vs FS | vs Rust |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 8.74 | 16.61 | 9.13 | 1.9x | 1.0x |
| 2 | 5.46 | 12.23 | 7.90 | 2.2x | 1.4x |
| 4 | 3.71 | 9.68 | 7.18 | 2.6x | 1.9x |
| 8 | 3.16 | 8.85 | 6.94 | 2.8x | 2.2x |
| 10 | 3.10 | 9.01 | 6.88 | **2.9x** | **2.2x** |

**Key Insight**:
- Speedup ratio improves with increasing thread count
- Achieves **2.9x** vs FreeSASA, **2.2x** vs Rust at threads=10
- Rust plateaus with increasing threads (parallel efficiency issue)

### Best Speedup Structures (50k+ atoms)

![Speedup Comparison](../../benchmarks/results/plots/speedup/comparison.png)

Top 5 structures with highest speedup at any thread count:

| Rank | vs FreeSASA | Atoms | Threads | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 1 | 8to0 | 673,884 | 8 | **3.05x** |
| 2 | 8to0 | 673,884 | 10 | 3.00x |
| 3 | 7zts | 280,200 | 10 | 2.99x |
| 4 | 6pwb | 157,148 | 10 | 2.97x |
| 5 | 7zts | 280,200 | 8 | 2.87x |

| Rank | vs Rust | Atoms | Threads | Speedup |
| ---: | --- | ---: | ---: | ---: |
| 1 | 5wk6 | 81,139 | 10 | **2.79x** |
| 2 | 5wlc | 161,064 | 10 | 2.69x |
| 3 | 8tqe | 100,035 | 10 | 2.58x |
| 4 | 8sep | 142,736 | 10 | 2.57x |
| 5 | 8uxk | 72,372 | 10 | 2.53x |

**Note**: threads=8 sometimes outperforms threads=10 (e.g., 8to0: 3.05x vs 3.00x).

**Rust outliers**: Small structures (<50k atoms) show anomalous speedup values (e.g., 17x) likely due to measurement noise. Large structure comparisons are more reliable.

---

## Execution Time Distribution

![SR Scatter Plot](../../benchmarks/results/plots/scatter/sr/grid.png)

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
| 0-500 | 0-500 | [View](../../benchmarks/results/plots/samples/0-500.png) |
| 500-1k | 500-1,000 | [View](../../benchmarks/results/plots/samples/500-1k.png) |
| 1k-2k | 1,000-2,000 | [View](../../benchmarks/results/plots/samples/1k-2k.png) |
| 2k-5k | 2,000-5,000 | [View](../../benchmarks/results/plots/samples/2k-5k.png) |
| 5k-10k | 5,000-10,000 | [View](../../benchmarks/results/plots/samples/5k-10k.png) |
| 10k-20k | 10,000-20,000 | [View](../../benchmarks/results/plots/samples/10k-20k.png) |
| 20k-50k | 20,000-50,000 | [View](../../benchmarks/results/plots/samples/20k-50k.png) |
| 50k-100k | 50,000-100,000 | [View](../../benchmarks/results/plots/samples/50k-100k.png) |
| 100k-200k | 100,000-200,000 | [View](../../benchmarks/results/plots/samples/100k-200k.png) |
| 200k+ | 200,000+ | [View](../../benchmarks/results/plots/samples/200kplus.png) |

---

---

## Key Takeaways

> **Why is Zig faster?** -> See [optimizations.md](../optimizations.md)

1. **Maximum effect on large structures**
   - **2.3x** speedup for 100k+ atoms
   - **2.9x** speedup for the largest structure (4.5M atoms)

2. **Consistent advantage**
   - Outperforms FreeSASA across all sizes (1k+ atoms)
   - Fastest at all thread counts

3. **Excellent parallel efficiency**
   - **26%** higher than FreeSASA
   - **97%** higher than RustSASA (at threads=10)

4. **Efficient code**
   - Same computation with **2.4x fewer instructions** than FreeSASA
   - Result of SIMD optimization

5. **Accurate results**
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
    --target 100000 --seed 42 \
    -o benchmarks/samples/stratified_100k.json
```

### Running

```bash
# Basic usage
./benchmarks/scripts/bench.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/bench.py --tool freesasa --algorithm lr --threads 1-10

# With sample file
./benchmarks/scripts/bench.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_100k.json \
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

Lee-Richards results using ~30k structures. RustSASA does not support LR.

### Dataset

![LR Dataset](../../benchmarks/results/plots/dataset/stratified_30k.png)

### Overall Statistics (threads=10)

| Tool | Structures | Median (ms) | Mean (ms) | P95 (ms) |
| --- | ---: | ---: | ---: | ---: |
| **Zig** | 29,998 | **9.47** | 40.11 | 166.88 |
| FreeSASA | 29,998 | 15.51 | 71.25 | 302.61 |

**Key Insight**: Zig is **1.64x** faster than FreeSASA (median)

### Speedup by Structure Size (threads=10)

| Size Bin | Count | vs FreeSASA |
| --- | ---: | ---: |
| 0-500 | 673 | 0.92x |
| 500-1k | 1,543 | 1.42x |
| 1k-2k | 4,274 | 1.53x |
| 2k-5k | 9,629 | 1.66x |
| 5k-10k | 5,396 | 1.71x |
| 10k-20k | 2,729 | 1.72x |
| 20k-50k | 1,450 | 1.74x |
| 50k-100k | 3,133 | **1.79x** |
| 100k-200k | 900 | **1.82x** |
| 200k+ | 271 | **1.82x** |

### Thread Scaling

| Threads | Zig (ms) | FreeSASA (ms) |
| ---: | ---: | ---: |
| 1 | 46.61 | 72.20 |
| 2 | 25.94 | 39.11 |
| 4 | 15.02 | 22.31 |
| 8 | 10.59 | 16.98 |
| 10 | **9.47** | 15.51 |

**Speedup from threads=1 to threads=10:**
- Zig: 46.61 -> 9.47 = **4.92x**
- FreeSASA: 72.20 -> 15.51 = **4.66x**

### Execution Time Distribution

![LR Scatter Plot](../../benchmarks/results/plots/scatter/lr/grid.png)

**Observations:**
- Overall 3-4x slower than SR (slice integration cost)
- Zig's advantage is maintained in LR as well

---

## Related Documents

- [batch.md](batch.md) - Batch processing benchmarks (proteome datasets)
- [md.md](md.md) - MD trajectory benchmarks
- [validation.md](validation.md) - SASA accuracy validation
- [optimizations.md](../optimizations.md) - Detailed optimization techniques
- [algorithm.md](../algorithm.md) - Algorithm details

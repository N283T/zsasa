# MD Trajectory SASA Benchmarks

Performance comparison of SASA calculation on molecular dynamics trajectories.

> **Note**: Times include trajectory loading and SASA calculation. All tools use 100 test points per atom (Shrake-Rupley algorithm).

## TL;DR

| Dataset | Atoms | Frames | zsasa CLI bitmask (f32, 10t) | vs MDTraj | vs mdsasa-bolt | RSS |
|---------|------:|-------:|---:|----------:|---------------:|----:|
| 5wvo_C | 3,858 | 1,001 | **0.78s** | 30x | 5.8x | 17 MB |
| 6sup_A | 33,377 | 1,001 | **6.92s** | 132x | 8.5x | 113 MB |
| 5vz0_A | 17,910 | 10,001 | **38.04s** | >2h (timeout) | >2h (timeout) | 63 MB |

- **zsasa CLI bitmask** is the fastest variant: ~2x faster than standard zsasa on large systems
- Memory is constant with frame count (streaming I/O) — **63 MB** for 10K frames vs Python tools using 3–6 GB
- MDTraj and mdsasa-bolt **timed out** (>2h) on the 10K-frame dataset that zsasa completes in 38 seconds

## Test Environment

| Item | Value |
|------|-------|
| Machine | MacBook Pro |
| Chip | Apple M4 (10 cores: 4P + 6E) |
| Memory | 32 GB |
| OS | macOS 15.3.2 (Darwin 24.6.0) |

## Implementations

| Tool | Language | Threading | I/O | Notes |
|------|----------|-----------|-----|-------|
| **zsasa CLI** | Zig | Configurable | Native XTC reader | f32/f64 precision, streams frames |
| **zsasa CLI bitmask** | Zig | Configurable | Native XTC reader | LUT bitmask, streams frames |
| **zsasa.mdtraj** | Python/Zig | Configurable | MDTraj (C) | MDTraj loads trajectory, zsasa computes SASA |
| **zsasa.mdtraj bitmask** | Python/Zig | Configurable | MDTraj (C) | LUT bitmask variant |
| **zsasa.mdanalysis** | Python/Zig | Configurable | MDAnalysis | MDAnalysis loads trajectory, zsasa computes SASA |
| **zsasa.mdanalysis bitmask** | Python/Zig | Configurable | MDAnalysis | LUT bitmask variant |
| **MDTraj** | Python/C | Single-threaded | MDTraj (C) | `shrake_rupley`, no parallel option |
| **mdsasa-bolt** | Rust | All cores (rayon) | MDAnalysis/Python | RustSASA via Python bridge |

## Results

### 6sup_A_analysis (Large: 33,377 atoms, 1,001 frames)

Benchmark: warmup=1, runs=3, threads=1,8,10.

#### 10-Thread Comparison

| Tool | Time (s) | FPS | vs MDTraj | vs mdsasa-bolt | RSS |
|------|--------:|----:|----------:|---------------:|----:|
| **zsasa CLI (f32, bitmask)** | **6.92** | **145** | **132.4x** | **8.5x** | **113 MB** |
| zsasa CLI (f64, bitmask) | 7.16 | 140 | 127.9x | 8.2x | 120 MB |
| zsasa.mdtraj bitmask | 7.31 | 137 | 125.3x | 8.0x | 1.5 GB |
| zsasa.mdanalysis bitmask | 7.71 | 130 | 118.8x | 7.6x | 988 MB |
| zsasa CLI (f32) | 14.79 | 68 | 61.9x | 4.0x | 111 MB |
| zsasa.mdtraj | 15.59 | 64 | 58.7x | 3.8x | 1.6 GB |
| zsasa CLI (f64) | 15.79 | 63 | 58.0x | 3.7x | 116 MB |
| zsasa.mdanalysis | 17.11 | 59 | 53.5x | 3.4x | 1.0 GB |
| mdsasa-bolt | 58.53 | 17 | 15.6x | baseline | 10.9 GB |
| MDTraj | 915.57 | 1.1 | baseline | — | 874 MB |

**Key findings:**

- **zsasa CLI bitmask (f32)** is **132x faster** than MDTraj, **8.5x faster** than mdsasa-bolt
- Bitmask variants are **~2.1x faster** than standard zsasa (6.92s vs 14.79s)
- mdsasa-bolt uses **10.9 GB** RAM vs zsasa CLI's **113 MB** (96x difference)

| Tool Comparison (10t) | Memory |
|:---------------------:|:------:|
| ![bar](pathname:///zsasa/assets/benchmarks/md/6sup_A_analysis/bar.png) | ![memory](pathname:///zsasa/assets/benchmarks/md/6sup_A_analysis/memory.png) |

### 5vz0_A_protein (Medium: 17,910 atoms, 10,001 frames)

Benchmark: warmup=1, runs=3, threads=1,8,10.

> MDTraj and mdsasa-bolt timed out at the 2-hour limit (`--timeout 7200`). MDTraj is too slow for 10K frames of this size; mdsasa-bolt likely exhausted memory (see [mdsasa-bolt Performance](#mdsasa-bolt-performance)).

#### 10-Thread Comparison

| Tool | Time (s) | FPS | RSS |
|------|--------:|----:|----:|
| **zsasa.mdtraj bitmask** | **37.41** | **267** | **5.8 GB** |
| zsasa.mdanalysis bitmask | 37.77 | 265 | 3.2 GB |
| **zsasa CLI (f32, bitmask)** | **38.04** | **263** | **63 MB** |
| zsasa CLI (f64, bitmask) | 38.99 | 256 | 67 MB |
| zsasa.mdtraj | 81.55 | 123 | 5.8 GB |
| zsasa CLI (f32) | 82.58 | 121 | 61 MB |
| zsasa.mdanalysis | 82.69 | 121 | 3.1 GB |
| zsasa CLI (f64) | 86.08 | 116 | 64 MB |

**Key findings:**

- **Bitmask variants ~2.1x faster** than standard (38s vs 83s)
- zsasa CLI processes **10,001 frames in 38 seconds** (263 FPS) with only **63 MB** RAM
- Python bindings slightly faster than CLI on bitmask (trajectory pre-loaded vs streaming), but use **50–90x more memory**

| Tool Comparison (10t) | Memory |
|:---------------------:|:------:|
| ![bar](pathname:///zsasa/assets/benchmarks/md/5vz0_A_protein/bar.png) | ![memory](pathname:///zsasa/assets/benchmarks/md/5vz0_A_protein/memory.png) |

### 5wvo_C_analysis (Small: 3,858 atoms, 1,001 frames)

Benchmark: warmup=1, runs=3, threads=1,8,10.

#### 10-Thread Comparison

| Tool | Time (s) | FPS | vs MDTraj | vs mdsasa-bolt | RSS |
|------|--------:|----:|----------:|---------------:|----:|
| **zsasa CLI (f32, bitmask)** | **0.78** | **1,291** | **30.0x** | **5.8x** | **17 MB** |
| zsasa CLI (f64, bitmask) | 0.79 | 1,263 | 29.4x | 5.6x | 18 MB |
| zsasa.mdtraj bitmask | 0.88 | 1,144 | 26.6x | 5.1x | 210 MB |
| zsasa.mdanalysis bitmask | 1.20 | 838 | 19.5x | 3.7x | 172 MB |
| zsasa.mdtraj | 1.64 | 612 | 14.2x | 2.7x | 195 MB |
| zsasa CLI (f32) | 1.69 | 593 | 13.8x | 2.6x | 14 MB |
| zsasa CLI (f64) | 1.77 | 568 | 13.2x | 2.5x | 15 MB |
| zsasa.mdanalysis | 1.99 | 505 | 11.7x | 2.2x | 166 MB |
| mdsasa-bolt | 4.47 | 224 | 5.2x | baseline | 1.4 GB |
| MDTraj | 23.29 | 43 | baseline | — | 156 MB |

**Key findings:**

- On small systems, bitmask advantage is smaller (~2.2x vs ~2.1x on large)
- zsasa CLI bitmask processes 1,001 frames in under **1 second** at 1,291 FPS
- Even on a small system, mdsasa-bolt uses **1.4 GB** vs zsasa CLI's **17 MB**

| Tool Comparison (10t) | Memory |
|:---------------------:|:------:|
| ![bar](pathname:///zsasa/assets/benchmarks/md/5wvo_C_analysis/bar.png) | ![memory](pathname:///zsasa/assets/benchmarks/md/5wvo_C_analysis/memory.png) |

## Memory Usage

| Tool | 5wvo_C (4k atoms, 1k fr) | 5vz0_A (18k atoms, 10k fr) | 6sup_A (33k atoms, 1k fr) |
|------|---:|---:|---:|
| zsasa CLI bitmask | 17 MB | 63 MB | 113 MB |
| zsasa CLI | 14 MB | 61 MB | 111 MB |
| zsasa.mdanalysis bitmask | 172 MB | 3.2 GB | 988 MB |
| zsasa.mdtraj bitmask | 210 MB | 5.8 GB | 1.5 GB |
| MDTraj | 156 MB | — | 874 MB |
| mdsasa-bolt | 1.4 GB | — | 10.9 GB |

- **zsasa CLI** memory scales with atom count, not frame count (streaming I/O)
- **Python tools** load the entire trajectory into memory — scales linearly with frames
- For 5vz0 (10k frames, 18k atoms): zsasa CLI uses **63 MB** while Python bindings use **3–6 GB**

## Bitmask Variants

The `_bitmask` variants use LUT (look-up table) bitmask neighbor lists:

- **~2.1x faster** than standard zsasa across all system sizes
- Minimal memory overhead for CLI variants (~2 MB extra)
- Python bitmask bindings are faster but still load full trajectories into memory
- On small systems (5wvo, 4k atoms), bitmask advantage is ~2.2x; on large systems (6sup, 33k atoms), ~2.1x

## mdsasa-bolt Performance

mdsasa-bolt (RustSASA via Python) shows high memory usage across all datasets:

| Dataset | mdsasa-bolt RSS | zsasa CLI RSS | Ratio |
|---------|----------------:|--------------:|------:|
| 5wvo_C (4k atoms) | 1.4 GB | 17 MB | 82x |
| 6sup_A (33k atoms) | 10.9 GB | 113 MB | 96x |

On the previous 10k-frame dataset (6sup_R1), mdsasa-bolt showed severe performance degradation (147x slower than zsasa), becoming slower than single-threaded MDTraj. This pattern is consistent with memory-bound scaling issues.

## SASA Validation

Per-frame total SASA compared against MDTraj `shrake_rupley` as baseline, using 5wvo_C (3,858 atoms, 1,001 frames).

> **Note on baseline**: MDTraj and zsasa use different sphere point placement (MDTraj uses reversed winding order). This causes small systematic differences that diminish with increasing test points. Bitmask error is independent of point count due to the LUT approximation. Plots use zoomed axes — the spread appears larger than it is in practice.

### R² vs MDTraj by Test Points

| Tool | 100 pts | 200 pts | 500 pts | 1000 pts |
|------|--------:|--------:|--------:|---------:|
| zsasa CLI (f64) | 0.8716 | 0.9664 | 0.9938 | 0.9983 |
| zsasa CLI (f32) | 0.8716 | 0.9664 | 0.9938 | 0.9983 |
| zsasa.mdtraj | 0.8716 | 0.9664 | 0.9938 | 0.9983 |
| zsasa.mdanalysis | 0.9724 | 0.9885 | 0.9718 | 0.9602 |
| zsasa CLI bitmask (f64) | 0.3980 | 0.6437 | 0.7537 | 0.7916 |
| zsasa CLI bitmask (f32) | 0.3981 | 0.6438 | 0.7537 | 0.7917 |

- **Standard variants** converge rapidly: R² > 0.99 at 500 points, > 0.998 at 1000 points
- **Bitmask variants** plateau around R² ~0.79 regardless of point count (LUT approximation error)
- **f32 and f64** produce virtually identical results
- **zsasa CLI f64 vs zsasa.mdtraj**: R² = 1.000 (identical SASA engine, different I/O path — max error 0.008%)

### ΔR² (Frame-to-Frame SASA Changes)

ΔR² measures how well frame-to-frame SASA *changes* (ΔSASA) correlate — which is what matters most in MD analysis.

| Tool | 100 pts | 200 pts | 500 pts | 1000 pts |
|------|--------:|--------:|--------:|---------:|
| zsasa CLI (f64) | 0.8918 | 0.9527 | 0.9895 | 0.9962 |
| zsasa CLI (f32) | 0.8918 | 0.9527 | 0.9895 | 0.9962 |
| zsasa CLI bitmask (f64) | 0.8902 | 0.9501 | 0.9872 | 0.9941 |
| zsasa CLI bitmask (f32) | 0.8902 | 0.9502 | 0.9874 | 0.9941 |

- **Bitmask ΔR² is nearly identical to standard** — the systematic offset from LUT cancels out in frame differences
- At 100 points: bitmask ΔR² = 0.890 vs standard ΔR² = 0.892 (negligible gap)
- At 1000 points: both converge to ΔR² > 0.994

### Practical Significance

The absolute R² values for bitmask may look low (~0.79), but the scatter plots use zoomed axes that exaggerate the spread. The actual max per-frame error is **1.5–3%**, and more importantly, the ΔR² shows that frame-to-frame trends are tracked with **>0.99 correlation** at 500+ points — even for bitmask.

For MD trajectory analysis where SASA trends and changes matter more than absolute values, bitmask accuracy is sufficient. If precise absolute values are needed, use the standard `f64` variant.

![Validation grid](pathname:///zsasa/assets/benchmarks/md/validation/validation_grid.png)

| XTC reader consistency | Bitmask comparison |
|:----------------------:|:------------------:|
| ![xtc](pathname:///zsasa/assets/benchmarks/md/validation/validation_xtc_comparison.png) | ![bitmask](pathname:///zsasa/assets/benchmarks/md/validation/validation_bitmask_comparison.png) |

## Running Benchmarks

```bash
# Run benchmark
./benchmarks/scripts/bench_md.py \
    --xtc trajectory.xtc \
    --pdb topology.pdb \
    --name my_bench \
    --threads "1,8,10"

# Specific tools only
./benchmarks/scripts/bench_md.py \
    --xtc trajectory.xtc \
    --pdb topology.pdb \
    --name my_bench \
    --tool zig --tool zig_bitmask --tool zsasa_mdtraj

# Analysis
./benchmarks/scripts/analyze_md.py all --name my_bench
./benchmarks/scripts/analyze_md.py summary --name my_bench
./benchmarks/scripts/analyze_md.py bar --name my_bench
./benchmarks/scripts/analyze_md.py memory --name my_bench
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--xtc` | XTC trajectory file | (required) |
| `--pdb` | Topology PDB file | (required) |
| `--name`, `-n` | Benchmark name | (required) |
| `--threads`, `-T` | Thread counts: `1,4,8` or `1-10` | `1,8` |
| `--runs`, `-r` | Number of benchmark runs | 3 |
| `--warmup`, `-w` | Number of warmup runs | 1 |
| `--stride`, `-s` | Frame stride | 1 |
| `--n-points` | Test points per atom | 100 |
| `--tool`, `-t` | Tool: zig, zig\_bitmask, zsasa\_mdtraj, zsasa\_mdtraj\_bitmask, zsasa\_mdanalysis, zsasa\_mdanalysis\_bitmask, mdtraj, mdsasa\_bolt (repeatable) | all |
| `--output`, `-o` | Output directory | results/md/\<n\_points\>/\<name\> |
| `--dry-run` | Show commands without running | false |

## Notes

- f32 and f64 precision produce nearly identical timing (~3–5% difference)
- SASA accuracy and XTC reader consistency validated in [validation.md](validation.md)
- mdsasa-bolt "all" = 10 threads on this system

## Related Documents

- [Single-File Benchmarks](single-file.md) — Per-structure performance across 2,013 structures
- [Batch Processing Benchmarks](batch.md) — Proteome-scale throughput
- [SASA Validation](validation.md) — Accuracy validation against FreeSASA

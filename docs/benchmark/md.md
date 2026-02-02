# MD Trajectory SASA Benchmark

Comparison of SASA calculation performance for molecular dynamics trajectories.

## Implementations

| Tool | Language | Threading | Notes |
|------|----------|-----------|-------|
| **zsasa** | Zig | Configurable | Native XTC reader, parallel SASA |
| **MDTraj** | Python/C | Single-threaded | `shrake_rupley`, no parallel option |
| **mdsasa-bolt** | Rust | All cores (rayon) | RustSASA via MDAnalysis |

## Benchmark Results

### 6QFK (Medium: 20k atoms)

**System**: 6QFK protein (20,391 atoms, 1,001 frames, 90.5 MB XTC)

```
┏━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━━━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Method      ┃ Threads ┃ Time (s) ┃ fps ┃ vs MDTraj ┃ vs mdsasa-bolt ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━━━━━━━╇━━━━━━━━━━━━━━━━┩
│ zsasa       │       1 │    55.29 │  18 │      6.7x │           0.5x │
│ zsasa       │       2 │    28.66 │  35 │     12.9x │           0.9x │
│ zsasa       │       4 │    15.04 │  67 │     24.5x │           1.7x │
│ zsasa       │       8 │    10.43 │  96 │     35.3x │           2.4x │
│ zsasa       │      10 │     9.25 │ 108 │     39.8x │           2.7x │
│ MDTraj      │       1 │   368.37 │   3 │         - │           0.1x │
│ mdsasa-bolt │     all │    24.97 │  40 │     14.8x │              - │
└─────────────┴─────────┴──────────┴─────┴───────────┴────────────────┘
```

### 5LTJ (Large: 11k atoms, 10k frames)

**System**: 5LTJ protein (11,487 atoms, 10,001 frames, 535 MB XTC)

```
┏━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━┳━━━━━┳━━━━━━━━━━━━━━━━┓
┃ Method      ┃ Threads ┃ Time (s) ┃ fps ┃ vs mdsasa-bolt ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━╇━━━━━╇━━━━━━━━━━━━━━━━┩
│ zsasa       │       8 │    56.91 │ 176 │          11.3x │
│ zsasa       │      10 │    49.91 │ 200 │          12.8x │
│ mdsasa-bolt │     all │   640.89 │  16 │              - │
└─────────────┴─────────┴──────────┴─────┴────────────────┘
```

MDTraj skipped (too slow for 10k frames).

**Environment**: macOS, Apple M2 Pro (10 cores)

## Key Findings

### zsasa vs MDTraj

- **Single-threaded**: zsasa is **6.7x faster** than MDTraj
- **10 threads**: zsasa is **40x faster** than MDTraj

MDTraj's `shrake_rupley` is single-threaded only (no `n_jobs` parameter).

### zsasa vs mdsasa-bolt

**6QFK (20k atoms, 1k frames)**:
- **Single-threaded zsasa vs all-core mdsasa-bolt**: mdsasa-bolt is 2x faster
- **10-thread zsasa vs all-core mdsasa-bolt**: zsasa is **2.7x faster**

**5LTJ (11k atoms, 10k frames)**:
- zsasa (10 threads) is **12.8x faster** than mdsasa-bolt (all cores)

mdsasa-bolt performance varies significantly between datasets. See [Parallelization Strategy](#parallelization-strategy) for details.

### Thread Scaling

| Threads | Time (s) | Speedup vs 1 thread |
|---------|----------|---------------------|
| 1 | 55.29 | 1.0x |
| 2 | 28.66 | 1.9x |
| 4 | 15.04 | 3.7x |
| 8 | 10.43 | 5.3x |
| 10 | 9.25 | 6.0x |

Near-linear scaling up to 4 threads, with diminishing returns beyond 8 threads.

## Running the Benchmark

```bash
./benchmarks/scripts/benchmark_md.py \
    trajectory.xtc \
    topology.pdb \
    --threads "1,2,4,8,10"
```

Options:
- `--threads`, `-t`: Comma-separated thread counts for zsasa (default: "1,4,8")
- `--runs`, `-n`: Number of runs per configuration (default: 3)

## Parallelization Strategy

Both zsasa and mdsasa-bolt use the same parallelization approach:

| Level | zsasa | mdsasa-bolt |
|-------|-------|-------------|
| Frame-level | Parallel (configurable threads) | Parallel (rayon, all cores) |
| Per-frame SASA | Single-threaded | Single-threaded (`n_threads=1` hardcoded) |

The key difference is **per-frame SASA implementation**:
- **zsasa**: Zig with SIMD optimizations (fast)
- **mdsasa-bolt**: rust_sasa library (slower)

### mdsasa-bolt Observations

- CPU utilization during mdsasa-bolt runs was lower than expected for "all cores"
- Performance varies significantly between datasets:
  - 6QFK (1k frames): 2.7x slower than zsasa
  - 5LTJ (10k frames): 12.8x slower than zsasa

**Hypothesis**: The performance degradation with more frames may be due to Python/Rust boundary overhead:

```python
# mdsasa-bolt collects all frames into Python lists before Rust call
input_atoms_per_frame = []
for _ in trajectory:
    input_atoms = [(tuple(pos), radius, resnum) for ...]  # Nested tuples
    input_atoms_per_frame.append(input_atoms)
```

- 10k frames × 11k atoms = 110M Python tuples created
- Inefficient memory layout (nested lists vs contiguous arrays)
- zsasa passes numpy arrays directly via CFFI (minimal overhead)

## Running the Benchmark

```bash
./benchmarks/scripts/benchmark_md.py \
    trajectory.xtc \
    topology.pdb \
    --threads "1,2,4,8,10"
```

Options:
- `--threads`, `-t`: Comma-separated thread counts for zsasa (default: "1,4,8")
- `--runs`, `-n`: Number of runs per configuration (default: 3)
- `--skip-mdtraj`: Skip MDTraj benchmark (for large structures)

## Notes

- All implementations use 100 test points per atom (Shrake-Rupley algorithm)
- Times include trajectory loading and SASA calculation
- mdsasa-bolt "all" = 10 threads on this system

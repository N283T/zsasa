# MD Trajectory SASA Benchmark

Comparison of SASA calculation performance for molecular dynamics trajectories.

## Implementations

| Tool | Language | Threading | Notes |
|------|----------|-----------|-------|
| **zsasa** | Zig | Configurable | Native XTC reader, parallel SASA |
| **MDTraj** | Python/C | Single-threaded | `shrake_rupley`, no parallel option |
| **mdsasa-bolt** | Rust | All cores (rayon) | RustSASA via MDAnalysis |

## Benchmark Results

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

**Environment**: macOS, Apple M2 Pro (10 cores)

## Key Findings

### zsasa vs MDTraj

- **Single-threaded**: zsasa is **6.7x faster** than MDTraj
- **10 threads**: zsasa is **40x faster** than MDTraj

MDTraj's `shrake_rupley` is single-threaded only (no `n_jobs` parameter).

### zsasa vs mdsasa-bolt

- **Single-threaded zsasa vs all-core mdsasa-bolt**: mdsasa-bolt is 2x faster
- **10-thread zsasa vs all-core mdsasa-bolt**: zsasa is **2.7x faster**

mdsasa-bolt uses Rust's rayon for parallelization but has no thread control - it always uses all available cores.

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

## Notes

- All implementations use 100 test points per atom (Shrake-Rupley algorithm)
- Times include trajectory loading and SASA calculation
- mdsasa-bolt "all" = 10 threads on this system

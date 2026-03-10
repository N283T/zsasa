---
sidebar_position: 0
sidebar_label: "Overview"
---

# Benchmarks

Performance and accuracy benchmarks for zsasa compared with FreeSASA, RustSASA, and Lahuta.

## Benchmark Suite

| Benchmark | What it measures | Key result |
|-----------|-----------------|------------|
| [Single-File](single-file.md) | Per-structure SASA speed across 2,013 structures | **1.88x** faster than FreeSASA, **1.84x** faster than RustSASA (f64, t=10) |
| [Batch](batch.md) | Proteome-scale directory processing throughput | **3.7x** faster than RustSASA, **7x** less memory at SwissProt scale |
| [MD Trajectory](md.md) | Molecular dynamics trajectory SASA performance | **82–96x** less memory than mdsasa-bolt |
| [Validation](validation.md) | Accuracy against FreeSASA reference | Mean error **<0.001%** (f64), **<0.05%** (f32) |

## Environment

All benchmarks were run on the same machine. Absolute times are environment-dependent — **relative speedup ratios** are the meaningful metric. All tools use identical parameters (`n_points=100`, `probe_radius=1.4Å`) unless otherwise noted.

## See Also

- [Comparison with Other Tools](/docs/comparison) — feature-level comparison with source code references

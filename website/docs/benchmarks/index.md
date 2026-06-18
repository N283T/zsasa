---
sidebar_position: 0
sidebar_label: "Overview"
---

# Benchmarks

This section summarizes the benchmark evidence used by the `zsasa` manuscript. The current headline results come from the pinned `zsasa` v0.6.0 benchmark harness in [`zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks), with figures copied from the paper and benchmark repositories.

The older exploratory benchmark pages have been replaced with the current pinned suite layout. The only pre-pinned result retained here is the SwissProt-scale batch benchmark, labeled as legacy context on the [Batch](batch.md#legacy-swissprot-benchmark-pre-pinned) page.

## Benchmark suites

| Suite | Purpose | Primary dataset | Main readouts |
| --- | --- | --- | --- |
| [Validation](validation.md) | Agreement with established SASA implementations | *E. coli* AFDB and 5wvo_C trajectory | R², relative error |
| [Batch throughput](batch.md) | Proteome-scale directory processing | *E. coli* and Human AFDB | Runtime, structures/s, RSS, scaling |
| [Single-file stress tests](single-file.md) | Large structures and parser-heavy cases | 8 curated structures up to 4.5M atoms | Runtime, RSS, parse/SASA timing |
| [MD trajectories](md.md) | Streaming frame-wise trajectory SASA | 5wvo_C, 6sup_A, 5vz0_A | Frames/s, RSS, comparator speedup |

## TL;DR

The clearest batch-performance view is throughput versus peak RSS: points in the upper-left are faster and more memory-efficient.

![Batch throughput vs peak RSS](pathname:///zsasa/assets/benchmarks/paper/batch/batch_tldr_throughput_vs_rss_2grid.png)

**Figure 1. Batch throughput vs peak memory.** The `zsasa` bitmask modes sit in the high-throughput, low-RSS region for both E. coli and Human AFDB.

## Test environment

The pinned current pinned benchmarks were run on one consumer laptop:

| Item | Value |
| --- | --- |
| Machine | MacBook Pro (`Mac16,1`) |
| Chip | Apple M4 |
| Cores | 10 total: 4 performance + 6 efficiency |
| Memory | 32 GB |
| OS | macOS 26.2 |
| Timing | hyperfine 1.20.0, with warmups and measured runs per suite |
| Tool pinning | Nix for native tools, uv for Python dependencies |
| `zsasa` version | v0.6.0 |

## Interpreting the numbers

- Speedups are comparator runtime divided by `zsasa` runtime; higher is better.
- RSS is peak resident set size; lower memory means a higher RSS-reduction ratio.
- Exact f64/f32 modes target numerical continuity with matched Shrake--Rupley outputs.
- Bitmask mode is a throughput-oriented approximation with an explicit validation envelope.
- Batch FreeSASA timings use a `freesasa_batch` wrapper because upstream FreeSASA has no native directory mode.

## Evidence sources

- Manuscript repository: [`N283T/zsasa-paper`](https://github.com/N283T/zsasa-paper)
- Benchmark harness and result tables: [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks)
- Feature comparison: [Comparison with Other Tools](/docs/comparison)

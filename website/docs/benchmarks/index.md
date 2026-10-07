---
sidebar_position: 0
sidebar_label: "Overview"
---

# Benchmarks

These pages report the `zsasa` 0.9.0 benchmark suite: agreement with established SASA implementations, proteome-scale batch throughput, large single structures, and MD trajectories. Every chart and table is generated from the result tables of [`zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks), archived at [doi:10.5281/zenodo.23184113](https://doi.org/10.5281/zenodo.23184113).

:::info[Numbers differ from the preprint]
The current [bioRxiv preprint](https://doi.org/10.64898/2026.06.29.733683) reports the `zsasa` 0.6.0 suite. These pages use the 0.9.0 rerun, so headline values differ slightly from that version; a revision of the preprint based on the 0.9.0 suite is in preparation. The 0.6.0 archive remains available at [doi:10.5281/zenodo.20577561](https://doi.org/10.5281/zenodo.20577561).
:::

<div data-chart="batch/map"></div>

Hover or focus a mark to read its values. Each chart also has a data table underneath it.

## Benchmark suites

| Suite | Question | Data | Readouts |
| --- | --- | --- | --- |
| [Validation](validation.md) | Does `zsasa` agree with established implementations? | 4,370 *E. coli* AFDB structures; 1,001 frames of 5wvo_C | Signed relative difference, R² |
| [Batch throughput](batch.md) | How fast is a whole directory? | *E. coli*, Human and SwissProt AFDB | Structures/s, peak RSS, thread scaling |
| [Single-file stress tests](single-file.md) | How does one large structure behave? | 8 structures, 10,919 to 4,506,416 atoms | Runtime, peak RSS, thread scaling |
| [MD trajectories](md.md) | How fast is frame-by-frame SASA? | 5wvo_C, 6sup_A, 5vz0_A | Frames/s, peak RSS |

## Reading the charts

- **Colour is the tool.** Yellow is `zsasa`. Blue is the reference implementation of each suite (FreeSASA, MDTraj), magenta the Rust implementation (RustSASA, mdsasa-bolt), and green the remaining comparator (Lahuta, PDBTools.jl).
- **Marker shape is the mode.** Circles are exact Shrake–Rupley runs and diamonds are bitmask runs; dashed lines also mark bitmask. Hollow markers are f32. Squares and triangles mark the Python integrations and the mmCIF parser paths, as named in each legend.
- **Whiskers** show one standard deviation over three measured runs.
- **Speedup** is comparator runtime divided by `zsasa` runtime, so higher is better.
- **Peak RSS** is the peak resident set size of the whole process.

Exact f64 and f32 modes aim to reproduce matched Shrake–Rupley outputs. Bitmask mode is a throughput-oriented approximation; its measured difference from FreeSASA is on the [validation page](validation.md).

## Test environment

All runs were collected on one consumer laptop:

| Item | Value |
| --- | --- |
| Machine | MacBook Pro (`Mac16,1`) |
| Chip | Apple M4 |
| Cores | 10 total: 4 performance + 6 efficiency |
| Memory | 32 GB |
| OS | macOS 26.2 |
| Timing | hyperfine 1.20.0, three measured runs per condition, `--prepare sync` |
| Tool pinning | Nix for native tools, uv for Python dependencies |
| `zsasa` version | 0.9.0 |

Batch FreeSASA timings use a `freesasa_batch` wrapper, because upstream FreeSASA has no native directory mode. The pinned comparator builds did not change between `zsasa` releases, so their batch timings are retained from the 0.6.0 session on the same machine; validation, single-file and trajectory runs were repeated for all tools.

## Evidence sources

- Benchmark harness and result tables: [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks)
- Archived results: [doi:10.5281/zenodo.23184113](https://doi.org/10.5281/zenodo.23184113)
- Preprint: [Nagae and Tomii, *bioRxiv* (2026)](https://doi.org/10.64898/2026.06.29.733683)
- Feature comparison: [Comparison with Other Tools](/docs/comparison)

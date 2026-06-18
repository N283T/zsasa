---
sidebar_position: 3
---

# Comparison with Other Tools

How `zsasa` compares with FreeSASA, RustSASA, and Lahuta for SASA-focused workflows. This page emphasizes factual feature differences and points to the benchmark pages for paper-era performance numbers.

## Overview

| Feature | zsasa | FreeSASA | RustSASA | Lahuta |
| --- | --- | --- | --- | --- |
| Language | Zig | C/C++ | Rust | C++ |
| Algorithms | Shrake-Rupley + Lee-Richards | Shrake-Rupley + Lee-Richards | Shrake-Rupley | Shrake-Rupley |
| Precision | f64/f32 selectable | f64 | f32 | f64 |
| Bitmask/LUT mode | ✅ 1-1024 points | — | — | ✅ 64/128/256 only |
| Structure input for SASA | PDB, mmCIF, BinaryCIF, SDF/MOL, JSON | PDB, mmCIF | PDB, mmCIF | AF2 model CIF/PDB |
| Directory / batch mode | ✅ native workflow/batch | ❌ single file only | ✅ native | ✅ native |
| Multi-chain SASA | ✅ | ✅ | ✅ | ❌ chain `A` hardcoded in SASA kernel |
| MD trajectory SASA | ✅ CLI + Python | — | △ via mdsasa-bolt | — |
| Trajectory formats | XTC, TRR, DCD, AMBER NetCDF | — | via MDAnalysis bridge | — |
| Python SASA API | ✅ multi-threaded | ✅ single-threaded | △ separate package | ❌ no SASA API |
| External dependencies | Zig + first-party `ztraj` module | none | several Rust crates | 13+ C++ libraries |

## Algorithms and precision

`zsasa` and FreeSASA support both Shrake-Rupley (SR) and Lee-Richards (LR). RustSASA and Lahuta support SR only for the SASA paths evaluated here.

`zsasa` exposes f64 and f32 modes. f64 is the default continuity path for matched FreeSASA-style SR output; f32 is useful when users want lower memory bandwidth and compatibility with f32-oriented tools. RustSASA uses f32 throughout its computation pipeline.

## Bitmask LUT mode

`zsasa` includes a bitmask lookup-table mode for high-throughput SR calculations. It is an approximation, not a numerically identical replacement for exact SR mode. In the paper-era static validation set, bitmask f64 at 128 points reached R² = 0.999811, mean relative difference = 0.662%, and max relative difference = 2.02% versus FreeSASA.

| Tool | Bitmask support | Point-count support |
| --- | --- | --- |
| zsasa | ✅ | 1-1024 |
| Lahuta | ✅ | 64, 128, 256 |
| RustSASA | — | — |
| FreeSASA | — | — |

See [SASA Validation](benchmarks/validation.md) and [Algorithms](guide/algorithms.mdx#bitmask-lut-optimization) for details.

## Input and parser behavior

`zsasa` uses native parsers that extract the fields needed for SASA calculation. This avoids rejecting structures because of unrelated fields and helps with very large assemblies.

RustSASA relies on `pdbtbx`; the benchmark input preparation had to normalize large PDB records and other parser-sensitive cases for comparator compatibility. Lahuta's evaluated SASA command is intended for AlphaFold2-model inputs and treats the SASA structure as chain `A`, which excludes it from the multi-chain single-file stress suite.

## Batch processing

FreeSASA has no native directory mode. The paper-era batch benchmarks therefore use `freesasa_batch`, a thin wrapper around the pinned FreeSASA C API, so that comparisons are multi-file workloads rather than shell loops.

Representative 10-thread batch results at 128 points:

| Workload | zsasa mode | Runtime | Throughput | RSS | Speedup vs FreeSASA batch |
| --- | --- | ---: | ---: | ---: | ---: |
| *E. coli* AFDB, 4,370 structures | f64 | 4.411 s | 991 str/s | 45.5 MiB | 2.94× |
| *E. coli* AFDB, 4,370 structures | bitmask f32 | 1.481 s | 2,951 str/s | 45.1 MiB | 8.77× |
| Human AFDB, 23,586 structures | f64 | 45.508 s | 518 str/s | 82.3 MiB | 2.94× |
| Human AFDB, 23,586 structures | bitmask f32 | 13.814 s | 1,707 str/s | 79.5 MiB | 9.70× |

See [Batch Processing Benchmarks](benchmarks/batch.md) for full tables and the legacy SwissProt note.

## Single-file stress behavior

The paper-era single-file suite uses eight curated structures up to 4,506,416 atoms. On the largest assembly, `zsasa` completed in 4.696 s in f64 mode and 3.788 s in bitmask f64 mode at 100 points and 10 threads. The same normalized input took 191.876 s with FreeSASA and 8.731 s with RustSASA.

| Workload | Mode | Runtime | RSS | Speedup vs FreeSASA |
| --- | --- | ---: | ---: | ---: |
| 9fqr, 4.5M atoms | zsasa f64 | 4.696 s | 1,615 MiB | 40.9× |
| 9fqr, 4.5M atoms | zsasa bitmask f64 | 3.788 s | 1,616 MiB | 50.6× |

See [Single-File Stress Benchmarks](benchmarks/single-file.md).

## MD trajectory support

`zsasa` provides native CLI trajectory processing and Python integrations. The benchmarked CLI path streams frames, keeping peak memory close to the current-frame working set. RustSASA trajectory support is represented by mdsasa-bolt, which uses an MDAnalysis front-end and can materialize much more trajectory data in memory.

Representative 10-thread trajectory results at 100 points:

| Workload | zsasa mode | Runtime | Frames/s | RSS | Speedup |
| --- | --- | ---: | ---: | ---: | --- |
| 5wvo_C, 1,001 frames | CLI bitmask f32 | 0.839 s | 1,194 | 22.6 MiB | 27.8× vs MDTraj |
| 6sup_A, 1,001 frames | CLI bitmask f32 | 6.949 s | 144 | 115.9 MiB | 132× vs MDTraj |
| 5vz0_A, 10,001 frames | CLI bitmask f32 | 38.056 s | 263 | 64.6 MiB | 86.5× vs mdsasa-bolt |

See [MD Trajectory Benchmarks](benchmarks/md.md).

## Build complexity

| Tool | Build command | Dependencies |
| --- | --- | --- |
| zsasa | `zig build -Doptimize=ReleaseFast` | Zig and first-party `ztraj` module |
| FreeSASA | `./configure && make` or CMake | none |
| RustSASA | `cargo build --release` | pdbtbx, pulp, mimalloc, rayon, and other crates |
| Lahuta | CMake/vcpkg | Boost, Eigen, gemmi, Highway, RDKit, lmdb, and more |

## Known limitations of zsasa

### Zig language maturity

Zig has not yet reached version 1.0. The language may introduce breaking changes between releases. The Python package communicates with the Zig core through a C ABI, which helps insulate Python users from Zig-internal changes as long as the C ABI is maintained.

### Benchmark scope

The current paper-era benchmark claims come from one consumer laptop and pinned comparator versions. Absolute runtimes will vary across hardware. Use the relative comparisons and the documented benchmark settings when interpreting the results.

## Links

| Tool | Repository |
| --- | --- |
| zsasa | [github.com/N283T/zsasa](https://github.com/N283T/zsasa) |
| FreeSASA | [github.com/mittinatten/freesasa](https://github.com/mittinatten/freesasa) |
| RustSASA | [github.com/maxall41/RustSASA](https://github.com/maxall41/RustSASA) |
| Lahuta | [github.com/bisejdiu/lahuta](https://github.com/bisejdiu/lahuta) |

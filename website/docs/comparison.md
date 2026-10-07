---
sidebar_position: 3
---

# Comparison with Other Tools

How `zsasa` compares with FreeSASA, RustSASA, and Lahuta for SASA-focused workflows. This page emphasizes factual feature differences and points to the benchmark pages for current pinned performance numbers.

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

`zsasa` includes a bitmask lookup-table mode for high-throughput SR calculations. It is an approximation, not a numerically identical replacement for exact SR mode. In the current pinned static validation set, bitmask f64 at 128 points reached R² = 0.999811, mean relative difference = 0.662%, and max relative difference = 2.02% versus FreeSASA.

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

FreeSASA has no native directory mode. The current pinned batch benchmarks therefore use `freesasa_batch`, a thin wrapper around the pinned FreeSASA C API, so that comparisons are multi-file workloads rather than shell loops.

Results at 10 threads and 128 sphere points, from the `zsasa` 0.9.0 suite:

<div data-chart="batch/bars"></div>

See [Batch Processing Benchmarks](benchmarks/batch.md) for full tables, thread scaling, and the mmCIF and SwissProt runs.

## Single-file stress behavior

The single-file suite uses eight curated structures up to 4,506,416 atoms, at 100 sphere points and 10 threads.

<div data-table="single/t10-pdb"></div>

Two of the large ratios are parser outliers rather than kernel scaling: FreeSASA on 8rbs and RustSASA on 5vyc. See [Single-File Stress Benchmarks](benchmarks/single-file.md) for the caveats.

## MD trajectory support

`zsasa` provides native CLI trajectory processing and Python integrations. The benchmarked CLI path streams frames, keeping peak memory close to the current-frame working set. RustSASA trajectory support is represented by mdsasa-bolt, which uses an MDAnalysis front-end and can materialize much more trajectory data in memory.

Results at 10 threads and 128 sphere points:

<div data-chart="md/bars"></div>

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

The benchmark claims come from one consumer laptop and pinned comparator versions. Absolute runtimes will vary across hardware. Use the relative comparisons and the documented benchmark settings when interpreting the results.

## Links

| Tool | Repository |
| --- | --- |
| zsasa | [github.com/N283T/zsasa](https://github.com/N283T/zsasa) |
| FreeSASA | [github.com/mittinatten/freesasa](https://github.com/mittinatten/freesasa) |
| RustSASA | [github.com/maxall41/RustSASA](https://github.com/maxall41/RustSASA) |
| Lahuta | [github.com/bisejdiu/lahuta](https://github.com/bisejdiu/lahuta) |

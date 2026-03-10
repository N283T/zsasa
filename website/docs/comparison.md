---
sidebar_position: 3
---

# Comparison with Other Tools

How zsasa compares to other modern SASA calculation tools.

This page focuses on **factual, verifiable differences** with source code references.
For performance numbers, see [Benchmarks](benchmarks/single-file.md).

## Overview

| Feature | zsasa | FreeSASA | RustSASA | Lahuta |
|---------|-------|----------|----------|--------|
| Language | Zig | C | Rust | C++ |
| Algorithms | SR + LR | SR + LR | SR only | SR only |
| Input formats | mmCIF, PDB, JSON | PDB, mmCIF | PDB (via pdbtbx) | AF2 model (CIF/PDB) |
| MD trajectory | Native (XTC, DCD) | — | △ (mdsasa-bolt) | — |
| Python API | ✅ (multi-threaded) | ✅ (single-threaded) | △ (separate package) | ❌ |
| Bitmask range | 1–1024 | — | — | 64/128/256 |
| External deps | None | None | pdbtbx, pulp, mimalloc | 13+ (Boost, Eigen, gemmi, Highway, RDKit, lmdb, ...) |
| Batch processing | ✅ Native | ❌ | ✅ | ✅ |
| Multi-chain | ✅ | ✅ | ✅ | ❌ (chain "A" hardcoded) |
| Precision | f64/f32 (selectable) | f64 | f32 only | f64 |
| Build system | `zig build` | autotools/CMake | Cargo | CMake + vcpkg |

## Algorithms

zsasa and FreeSASA support both Shrake-Rupley (SR) and Lee-Richards (LR) algorithms. RustSASA and Lahuta support SR only.

- RustSASA: [lib.rs](https://github.com/maxall41/RustSASA/blob/main/src/lib.rs) — "Shrake-Rupley algorithm"
- Lahuta: [sasa.hpp](https://github.com/bisejdiu/lahuta/blob/main/core/src/analysis/sasa/sasa.hpp) — SR computation only

## Input & Parsing

### zsasa: Native parser, zero external dependencies

zsasa includes a native mmCIF/PDB parser written in Zig. It extracts only the fields needed for SASA calculation (coordinates, element, residue name, atom name), avoiding errors from irrelevant data fields. JSON input is also supported for maximum flexibility.

### RustSASA: pdbtbx dependency causes parse errors

RustSASA uses [pdbtbx](https://github.com/maxall41/RustSASA/blob/main/Cargo.toml) for PDB parsing. This external parser attempts to validate the entire PDB structure, causing failures on:

- **Multi-character chain IDs** (common in large complexes from mmCIF → PDB conversion)
- **Hybrid-36 residue/serial numbers** (structures with >9999 residues or >99999 atoms)
- **Missing CRYST1 Z values** (common in EM/NMR structures)

Our benchmarking required a [dedicated preprocessing script](https://github.com/N283T/zsasa/blob/main/benchmarks/scripts/generate_pdb.py) to work around these pdbtbx limitations — reassigning chain IDs, wrapping serial/residue numbers, and fixing CRYST1 records.

### Lahuta: AlphaFold2 model input required

Lahuta's `sasa-sr` command requires AlphaFold2 model inputs. File-based sources must pass the `--is_af2_model` flag, otherwise an error is thrown:

> *"sasa-sr expects AlphaFold2 model inputs. For file-based sources, pass --is_af2_model"*

Source: [sasa_sr_spec.cpp L211–214](https://github.com/bisejdiu/lahuta/blob/main/cli/specs/sasa_sr_spec.cpp#L211-L214)

## Chain Support

Lahuta's SASA kernel hardcodes `chain_id = "A"`:

```cpp
constexpr std::string_view chain_id = "A";
```

Source: [kernel.hpp L164, L192](https://github.com/bisejdiu/lahuta/blob/main/core/src/analysis/sasa/kernel.hpp#L164)

This means multi-chain structures are treated as single-chain. zsasa, FreeSASA, and RustSASA all support multi-chain structures natively.

## MD Trajectory Support

| Tool | XTC | DCD | Integration |
|------|-----|-----|-------------|
| zsasa | ✅ Native CLI + Python | ✅ Native CLI + Python | MDTraj, MDAnalysis |
| FreeSASA | — | — | — |
| RustSASA | △ via mdsasa-bolt | — | MDAnalysis (Python bridge) |
| Lahuta | — | — | — |

Lahuta's `sasa-sr` has no trajectory options — grep for "trajectory", "xtc", or "dcd" in [sasa_sr_spec.cpp](https://github.com/bisejdiu/lahuta/blob/main/cli/specs/sasa_sr_spec.cpp) returns no matches.

### mdsasa-bolt (RustSASA) memory issues

RustSASA's trajectory support is via [mdsasa-bolt](https://github.com/maxall41/RustSASA#mdsasa-bolt), a separate Python package that bridges MDAnalysis to RustSASA. It has significant memory overhead:

| Dataset | mdsasa-bolt RSS | zsasa CLI RSS | Ratio |
|---------|----------------|---------------|-------|
| 5vz0 (1K frames, 28K atoms) | 1.4 GB | 17 MB | **82x** |
| 5wvo (1K frames, 97K atoms) | 10.9 GB | 113 MB | **96x** |

On 10K-frame datasets, mdsasa-bolt **timed out** (>2 hours) where zsasa completed in 38 seconds.

See [MD Trajectory Benchmarks](benchmarks/md.md) for full results.

## Bitmask LUT Optimization

| Tool | Bitmask support | Supported n_points | Range |
|------|----------------|--------------------|-------|
| zsasa | ✅ | 1–1024 (any value) | Continuous |
| Lahuta | ✅ | 64, 128, 256 | 3 discrete values |
| RustSASA | — | — | — |
| FreeSASA | — | — | — |

zsasa's bitmask supports any `n_points` from 1 to 1024 ([bitmask_lut.zig](https://github.com/N283T/zsasa/blob/main/src/bitmask_lut.zig#L8)), while Lahuta is limited to 64/128/256 ([sasa_sr_spec.cpp L172–173](https://github.com/bisejdiu/lahuta/blob/main/cli/specs/sasa_sr_spec.cpp#L172-L173)).

See [Bitmask LUT Algorithm](guide/algorithms.mdx#bitmask-lut-optimization) for how it works.

## Python Bindings

| Tool | Python SASA API | How |
|------|----------------|-----|
| zsasa | ✅ `from zsasa import calculate_sasa` | Native C extension via Zig, multi-threaded |
| FreeSASA | ✅ `import freesasa` | Native C extension, **single-threaded only** |
| RustSASA | △ `rust-sasa-python` (separate package) | PyO3 bridge |
| Lahuta | ❌ | pybind11 bindings exist but expose pipeline framework only — no SASA API |

Lahuta has pybind11 bindings, but they expose the pipeline framework — not SASA calculation. The Python code contains no reference to "sasa" at all.

Source: [interop/python/](https://github.com/bisejdiu/lahuta/tree/main/interop/python)

## Precision

| Tool | Precision | Selectable |
|------|-----------|------------|
| zsasa | f64 or f32 | ✅ `--precision=f32/f64` |
| FreeSASA | f64 | — |
| RustSASA | f32 | — |
| Lahuta | f64 | — |

RustSASA uses f32 throughout its entire computation pipeline — coordinates, sphere points, and SASA output are all single-precision.

Source: [lib.rs](https://github.com/maxall41/RustSASA/blob/main/src/lib.rs) — `Vec<f32>` for coordinates, `f32` for all SIMD operations

zsasa defaults to f64 for maximum accuracy, with f32 available as an option for trajectory processing or when comparing with f32-based tools. See [Validation Benchmarks](benchmarks/validation.md) for accuracy comparison.

## Batch Processing

| Tool | Directory batch | Parallelism |
|------|----------------|-------------|
| zsasa | ✅ Native `zsasa batch` | File-level multi-threading |
| FreeSASA | ❌ Single file only | Single file at a time |
| RustSASA | ✅ `rust-sasa dir/ out/` | File-level multi-threading |
| Lahuta | ✅ `--directory` / `--files` | Pipeline-based |

FreeSASA processes one file at a time and has no built-in batch/directory mode. Processing a proteome-scale dataset requires external scripting with sequential invocations.

See [Batch Benchmarks](benchmarks/batch.md) for performance comparison across proteome-scale datasets.

## Build Complexity

| Tool | Build command | Dependencies |
|------|--------------|--------------|
| zsasa | `zig build -Doptimize=ReleaseFast` | None (Zig only) |
| FreeSASA | `./configure && make` | None (C only) |
| RustSASA | `cargo build --release` | pdbtbx, pulp, mimalloc, rayon |
| Lahuta | CMake | Boost, Eigen, gemmi, Google Highway (via distopia), RDKit, lmdb, mimalloc, spdlog, libdivide, simde, ZLIB, ... |

zsasa and FreeSASA have zero external dependencies. RustSASA pulls in several crates. Lahuta bundles [13+ external libraries](https://github.com/bisejdiu/lahuta/tree/main/core/external) — building from source requires resolving all of them.

### SIMD Architecture

| Tool | SIMD approach | Architecture independence |
|------|--------------|--------------------------|
| zsasa | Zig `@Vector` | ✅ Compiler auto-targets native ISA (x86 SSE/AVX, ARM NEON, etc.) |
| RustSASA | `pulp` crate | ✅ Runtime dispatch via `ARCH.dispatch` ([lib.rs](https://github.com/maxall41/RustSASA/blob/main/src/lib.rs)) |
| Lahuta | SIMDe (AVX2 API) | △ Written against AVX2 intrinsics — SIMDe emulates on non-x86 |
| FreeSASA | None | — |

zsasa uses Zig's built-in `@Vector` primitives — the compiler maps them directly to the target architecture's SIMD instructions without any abstraction layer.

RustSASA uses [pulp](https://crates.io/crates/pulp) for runtime SIMD dispatch — it detects available instruction sets at runtime and selects the best implementation.

Lahuta uses [SIMDe](https://github.com/simd-everywhere/simde) but with explicit AVX2 API calls (`#include <simde/x86/avx2.h>`). All SIMD functions are named with `_avx2` suffix (e.g., `fast_rsqrt_avx2`, `octa_encode_avx2`). SIMDe can **emulate** AVX2 on non-x86 platforms, but the code is fundamentally written against the x86 AVX2 instruction set.

Source: [bitmask_simd.hpp](https://github.com/bisejdiu/lahuta/blob/main/core/src/analysis/sasa/methods/bitmask_simd.hpp)

## Performance Summary

### Single-file (100 test points, median across 128 structures)

| Metric | zsasa f64 (t=10) | vs FreeSASA | vs RustSASA |
|--------|-----------------|-------------|-------------|
| Median speedup | — | **1.88x** | **1.84x** |
| Thread scaling (t=1→10) | **2.71x** | 1.99x | 1.39x |

RustSASA shows [pathological slowdowns](benchmarks/single-file.md#rustsasa-pathological-cases) on certain structures (up to **280x slower** than zsasa).

See [Single-file Benchmarks](benchmarks/single-file.md) for full results.

### Batch (128 test points, E. coli K-12 proteome, 4,370 structures)

| Tool | Time | RSS | vs zsasa |
|------|------|-----|----------|
| zsasa bitmask (f32) | **1.42s** | 43 MB | baseline |
| Lahuta bitmask | 2.01s | 291 MB | 1.4x slower, 6.8x more memory |
| RustSASA | 5.24s | 169 MB | 3.7x slower, 3.9x more memory |

At SwissProt scale (550K structures), zsasa uses **157 MB** vs RustSASA **1.1 GB** (7.2x less) and Lahuta bitmask **2.2 GB** (14x less).

See [Batch Benchmarks](benchmarks/batch.md) for full results including memory scaling charts.

## Known Limitations of zsasa

### Zig language maturity

Zig has not yet reached version 1.0. The language may introduce breaking changes between releases, which could require updates to the zsasa codebase. Unreported compiler bugs are also possible in pre-1.0 software.

**Mitigation:** The Python package (`zsasa`) communicates with the Zig core through the C ABI, which is stable and standardized. This means Python users are insulated from Zig-specific changes — even if the Zig internals are updated, the Python API remains unaffected as long as the C ABI interface is maintained.

## Links

| Tool | Repository |
|------|-----------|
| zsasa | [github.com/N283T/zsasa](https://github.com/N283T/zsasa) |
| FreeSASA | [github.com/mittinatten/freesasa](https://github.com/mittinatten/freesasa) |
| RustSASA | [github.com/maxall41/RustSASA](https://github.com/maxall41/RustSASA) |
| Lahuta | [github.com/bisejdiu/lahuta](https://github.com/bisejdiu/lahuta) |

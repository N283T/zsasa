---
sidebar_position: 10
---

# Changelog

All notable changes to zsasa. See [GitHub Releases](https://github.com/N283T/zsasa/releases) for full details.

## [Unreleased]

Changes since v0.2.1 (2026-02-26 – 2026-03-10).

### Added

- **Bitmask variants for Python MD wrappers**: `use_bitmask` option for MDTraj, MDAnalysis, XTC, and DCD integrations (#275)
- **JSONL streaming batch output**: `--format=jsonl` for memory-efficient batch results (#227, #235)
- **Comparison page**: factual feature comparison with FreeSASA, RustSASA, and Lahuta, backed by source code references (#292)
- **Landing page**: hero section with highlight cards and navigation buttons (#291)
- **Python autodoc**: pdoc-generated API documentation deployed to `/python-autodoc/` (#290)
- **Benchmarks overview page**: summary table linking all benchmark categories (#294)
- **Changelog page**: release history on the documentation site (#294)

### Changed

- **Documentation site overhaul**:
  - Split CLI reference into Commands, Input, and Output pages (#289)
  - Rewrote single-file, batch, MD trajectory, and validation benchmark pages with new results (#280, #282, #283, #284)
  - Enhanced bitmask LUT documentation with implementation details (#287, #288)
  - Streamlined docs index to hub layout, removed feature duplication with landing page (#294)

### Performance

- **mmap file reading**: replaced `readToEndAlloc` with memory-mapped I/O for structure files (#229)
- **Flat buffer NeighborList**: replaced dynamic `ArrayList` with pre-allocated flat buffers (#230)
- **Trajectory parallel workers**: aligned with batch allocator pattern for lower overhead (#231)
- **64KB write buffer** for JSONL output (#228)

### Fixed

- Corrected Lahuta placeholder URL and language (Zig → C++) in batch benchmarks (#293)
- Corrected RustSASA precision (f64 → f32) in single-file benchmarks (#293)
- Fixed `pathname://` image paths for Docusaurus compatibility (#284)
- Various PDB generation fixes: chain name shortening, serial number wrapping, CRYST1 Z value handling (#251, #252, #253)

## [v0.2.1](https://github.com/N283T/zsasa/releases/tag/v0.2.1) — 2026-02-25

### Changed

- Relaxed bitmask `n_points` from fixed 64/128/256 to any value 1–1024, with internal storage expanded from `[4]u64` to `[16]u64`

## [v0.2.0](https://github.com/N283T/zsasa/releases/tag/v0.2.0) — 2026-02-25

### Added

- **Bitmask-optimized Shrake-Rupley** (`--use-bitmask`): precomputed occlusion bitmask LUT with SIMD 4-neighbor batching
- Batch and trajectory mode bitmask support (shared LUT, built once)
- Python `use_bitmask=True` parameter for `calculate_sasa()` and `calculate_sasa_batch()`
- C API bitmask exports (`zsasa_calc_sr_bitmask`, etc.)

### Changed

- **BREAKING**: CLI now requires subcommands: `zsasa calc`, `zsasa batch`, `zsasa traj`
- **BREAKING**: Removed `--parallelism` option
- Trajectory mode: `--include-hydrogens` is now the default

### Removed

- Pipeline parallelism mode

### Performance

- E. coli batch (f32, 10 threads, 128 points): **4.92s → 2.57s** with bitmask SR (1.9x speedup)

## [v0.1.3](https://github.com/N283T/zsasa/releases/tag/v0.1.3) — 2026-02-25

### Added

- Directory batch processing C API and Python `process_directory()`
- Docusaurus documentation site

### Fixed

- json_writer performance regression: reverted to in-memory string building (~8x speedup)

### Removed

- Streaming output (`--stream`, `--stream-format`, `--stream-output`)

## [v0.1.2](https://github.com/N283T/zsasa/releases/tag/v0.1.2) — 2026-02-22

### Added

- **DCD trajectory reader** (native Zig): NAMD/CHARMM DCD support with endianness auto-detection
- JSON streaming output for batch processing
- Zig package manager distribution (`zig fetch`)
- TOML format support for custom classifier configs
- Fuzz tests for parsers

## [v0.1.1](https://github.com/N283T/zsasa/releases/tag/v0.1.1) — 2026-02-22

### Added

- Pre-built wheels for Linux, macOS, and Windows via cibuildwheel
- Windows support (Zig build, tests, CLI, Python bindings)
- PyPI publish workflow with OIDC trusted publishing
- `CODE_OF_CONDUCT.md`, `CITATION.cff`, issue/PR templates

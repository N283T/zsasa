---
sidebar_position: 10
---

# Changelog

All notable changes to zsasa. See [GitHub Releases](https://github.com/N283T/zsasa/releases) for full details.

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

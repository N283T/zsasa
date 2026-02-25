# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-02-25

### Added

- **Bitmask-optimized Shrake-Rupley algorithm** (`--use-bitmask`): precomputed occlusion bitmask LUT with O(1) octahedral encoding for direction lookup, replacing per-point neighbor testing (#197)
- **SIMD 4-neighbor batching** for bitmask SR: processes 4 neighbors simultaneously with `@Vector(4, T)`, branchless octahedral encoding via `@select`, and combined mask accumulation (#198)
- **Batch-mode `--use-bitmask` support**: shared LUT across all files in batch processing, avoiding redundant ~20ms LUT construction per file (#198)
- **Trajectory-mode `--use-bitmask` support**: build LUT once, reuse across all frames in both sequential and batch-parallel paths
- **Python bindings `use_bitmask` parameter**: added `use_bitmask=True` option to `calculate_sasa()` and `calculate_sasa_batch()`, with pass-through to MDTraj, MDAnalysis, XTC, and DCD integrations
- **C API bitmask exports**: `zsasa_calc_sr_bitmask`, `zsasa_calc_sr_batch_bitmask`, `zsasa_calc_sr_batch_bitmask_f32` for FFI access to bitmask LUT optimization

### Changed

- **BREAKING**: CLI now requires subcommands: `zsasa calc`, `zsasa batch`, `zsasa traj`
- **BREAKING**: Removed `--parallelism` option (`calc` uses atom-level, `batch` uses file-level parallelism)
- **Trajectory mode**: `--include-hydrogens` is now the default (hydrogen atoms included). Use `--no-hydrogens` to exclude. MD trajectories typically include all atoms.
- CI: removed Windows from PR checks (linux + macOS only); Windows builds remain in release workflow (#199)

### Removed

- Pipeline parallelism mode (`--parallelism=pipeline`)
- Atom-level batch parallelism (`--parallelism=atom`)

### Performance

- E.coli proteome batch (f32, 10 threads, 128 points): **4.92s → 2.57s** with bitmask SR (1.9x speedup, within 12% of lahuta reference)

## [0.1.3] - 2026-02-25

### Added

- **Directory batch processing C API**: `zsasa_batch_dir_*` functions for processing all structure files in a directory (#191)
- **Python bindings for directory batch processing**: `process_directory()` function and `BatchDirResult` dataclass wrapping the C API (#193)
  - Process all supported structure files in a directory from Python
  - Support for SR/LR algorithms, all classifiers, threading, output directory
  - Per-file results: filename, atom count, total SASA, status
  - Error mapping: `ValueError`, `FileNotFoundError`, `MemoryError`, `RuntimeError`
- Docusaurus documentation site with GitHub Pages deployment (#186)

### Changed

- Slimmed down README, added uv install instructions (#190)
- CI: skip workflow for website-only changes (#189)

### Fixed

- **json_writer performance regression**: Reverted from unbuffered streaming writes to in-memory string building + single writeAll, fixing ~8x slowdown in batch processing caused by millions of write syscalls (#157 regression, #194)
- Updated benchmark docs for new dataset (#188)
- Fixed absolute paths for benchmark images (#187)

### Removed

- **Streaming output** (`--stream`, `--stream-format`, `--stream-output`): Removed StreamWriter module and CLI options to reduce code complexity (#194)

## [0.1.2] - 2026-02-22

### Added

- JSON streaming output for batch processing (`--stream`): stream results as NDJSON or JSON array as each file completes (#157)
- Stream format selection (`--stream-format`): choose between `ndjson` (default) and `json` array format (#157)
- Stream output destination (`--stream-output`): write stream to file instead of stdout (#157)
- Writer-based streaming JSON output for per-file results, reducing memory usage (#157)
- Zig package manager (zon) distribution: zsasa can now be used as a library dependency via `zig fetch` (#160)
- Public library API in `root.zig`: `shrake_rupley`, `lee_richards`, `types`, `pdb_parser`, `mmcif_parser`, `json_parser`, `classifier`, `analysis`
- Fuzz tests for CIF tokenizer, PDB parser, and mmCIF parser using Zig's built-in `std.testing.fuzz()` (#161)
- TOML format support for custom classifier configs (`--config=file.toml`): human-friendly alternative to FreeSASA format with auto-detection by file extension (#158)
- **DCD trajectory reader** (native Zig): read NAMD/CHARMM DCD binary trajectories without external dependencies (#154)
  - Zig DCD reader (`src/dcd.zig`) with endianness auto-detection and CHARMM extension support
  - C API: `zsasa_dcd_open`, `zsasa_dcd_close`, `zsasa_dcd_read_frame`, `zsasa_dcd_get_natoms`
  - Python DCD reader (`zsasa.dcd`): `DcdReader` class and `compute_sasa_trajectory()` function
  - CLI: `zsasa traj` now supports `.dcd` files with auto-detection by file extension
- Zig library API reference documentation (`docs/zig-api/`): types, algorithms, parsers, classifier, analysis (#156)
- Python pdoc auto-generated API documentation (`scripts/generate-python-docs.sh`) (#156)
- Documentation site with MkDocs Material and GitHub Pages deployment (#184)
- `zig build docs` step for interactive Zig autodoc generation (#184)

### Changed

- `build.zig`: Removed boilerplate template comments (176 → 63 lines)
- `build.zig.zon`: Synced version with `build.zig`
- Homepage and Documentation URLs now point to GitHub Pages (#184)

### Fixed

- `simd.zig`: Fixed `std.math.atan2` comptime_float errors with explicit `@as(f64, ...)` casts
- `simd.zig`: Corrected `fastAtan2` test tolerance for negative quadrant inputs (0.005 → 0.07)
- `root.zig`: Re-enabled `shrake_rupley` and `lee_richards` in test block (previously excluded due to transitive simd test failure)

## [0.1.1] - 2026-02-22

### Added

- `CODE_OF_CONDUCT.md` (Contributor Covenant v2.1)
- `CITATION.cff` (CFF 1.2.0 format for academic citation)
- GitHub issue templates (bug report, feature request) in YAML form
- Pull request template
- `speedup_by_threads.png` plot generation in `analyze.py large` (thread scaling for 50k+ atoms)
- CI status, license, Zig, and Python badges to READMEs
- Pre-built wheel distribution via cibuildwheel (Linux x86_64/aarch64, macOS x86_64/arm64, Windows x86_64)
- Windows support: Zig build, tests, CLI, Python bindings
- PyPI publish workflow (`.github/workflows/publish.yml`) with OIDC trusted publishing
- `python -m ziglang` fallback in `hatch_build.py` for Zig discovery

### Changed

- `python/pyproject.toml`: Updated author name, added Documentation/Issues/Changelog URLs

### Fixed

- README: Corrected MD trajectory benchmark data (6sup_A_analysis: 4.3x speedup, verified against actual data)
- README: Fixed broken documentation link (`docs/python.md` → `docs/python-api/`)
- README: Fixed broken image reference for thread scaling plot

## [0.1.0] - 2026-01-31

### Added

- **PDB file format support**
  - Fixed-width PDB parser (`src/pdb_parser.zig`)
  - Auto-detection of `.pdb` and `.ent` files
  - MODEL/ENDMDL, alternate location, chain filtering
  - Element inference from atom names

- **mmCIF parser** (internal)
  - Native mmCIF parser (`src/mmcif_parser.zig`)
  - No external dependencies (replaces gemmi for CLI)

- **Gzip support**
  - Transparent decompression for `.json.gz` and `.cif.gz`
  - Streaming decompression with zlib

- **Batch processing** (directory input)
  - Process entire directories: `zsasa ./input_dir/ ./output_dir/`
  - File-level parallelism with work stealing
  - Per-thread arena allocators for memory efficiency
  - Progress bar with file count
  - `--parallelism` option for concurrent file processing
  - Duplicate coordinate detection and warning

- **f32 precision option** (`--precision=f32`)
  - Single-precision mode for reduced memory usage
  - Comptime generics for zero-cost abstraction
  - ~Same speed, slightly lower accuracy

- **AVX-512 auto-optimization**
  - 16-wide SIMD on supported CPUs
  - Automatic detection and fallback (16 → 8 → 4 → scalar)

- **Benchmark infrastructure**
  - `benchmarks/scripts/run.py` - Unified benchmark runner
  - `benchmarks/scripts/analyze.py` - Results analysis and plotting
  - `benchmarks/scripts/sample.py` - Stratified sampling for large datasets
  - `benchmarks/scripts/build_index.py` - Dataset indexing
  - Full PDB dataset benchmark (238,124 structures)
  - Batch benchmark: Zig +7% faster than RustSASA

- **Example files** (`examples/`)
  - Sample PDB and mmCIF files for quick testing
  - README with usage examples

- **8-wide SIMD optimization** for Shrake-Rupley algorithm
  - `@Vector(8, f64)` for processing 8 atoms in parallel
  - Tiered processing: 8-wide → 4-wide → scalar for remaining atoms
  - ~16% speedup on large structures (4V6X: 237k atoms)

- **Fast trigonometry** for Lee-Richards algorithm
  - Polynomial approximations for `acos` and `atan2`
  - ~37% speedup on large structures (4V6X: 1021ms → 743ms)
  - Accuracy within 0.3% of reference (well within 2% tolerance)

- **Area difference column** in benchmark output
  - Shows percentage difference between Zig and FreeSASA C results

- **Analysis options** (CLI)
  - `--per-residue` - Per-residue SASA aggregation
  - `--rsa` - Relative Solvent Accessibility calculation
  - `--polar` - Polar/nonpolar SASA classification

- **Python bindings** (`python/zsasa`)
  - C ABI shared library (`libzsasa.dylib/.so/.dll`)
  - NumPy-based Python API with ctypes bindings
  - Both SR and LR algorithms supported
  - `calculate_sasa(coords, radii, algorithm="sr"|"lr", ...)` function
  - RSA functions: `calculate_rsa()`, `calculate_rsa_batch()`, `get_max_sasa()`
  - Per-residue aggregation: `aggregate_by_residue()`, `ResidueResult` class
  - 161 unit tests with pytest

- **Python integrations** (optional dependencies)
  - Gemmi integration for mmCIF/PDB file loading
  - BioPython integration for structure file support
  - Biotite integration for structure analysis workflows
  - MDTraj integration (`zsasa.mdtraj`) - drop-in replacement for `mdtraj.shrake_rupley()`
  - MDAnalysis integration (`zsasa.mdanalysis`) - `SASAAnalysis` class compatible with `AnalysisBase`

- **XTC trajectory reader** (native Zig)
  - Zig port of GROMACS libxdrfile (BSD-2-Clause)
  - Python XTC reader (`zsasa.xtc`) - no MDTraj/MDAnalysis dependency required

- **Trajectory subcommand** (`zsasa traj`)
  - `zsasa traj trajectory.xtc topology.pdb` - CLI trajectory mode
  - Frame-level batch parallelism with work-stealing
  - `--stride=N`, `--start=N`, `--end=N` frame filtering
  - `--batch-size=N` for controlling parallel batch size
  - Default f32 precision for speed in trajectory mode

- **Hydrogen and HETATM filtering**
  - `--include-hydrogens` - Include H/D atoms (default: excluded)
  - `--include-hetatm` - Include HETATM records (default: excluded)
  - Applied to both PDB and mmCIF parsers

- **SASA validation infrastructure**
  - `benchmarks/scripts/validation.py` - Accuracy validation vs FreeSASA C
  - `benchmarks/scripts/validation_md.py` - MD trajectory validation across implementations
  - Lee-Richards validation with E. coli proteome (R²=1.0)

- **MD trajectory benchmarks**
  - `benchmarks/scripts/bench_md.py` - Hyperfine-based MD benchmark
  - `benchmarks/scripts/analyze_md.py` - MD analysis and plots
  - E. coli proteome batch benchmark

- **Timing breakdown** (`--timing` flag)
  - Reports detailed timing for each phase: parsing, classification, SASA calculation, output
  - Enables fair performance comparison by measuring SASA-only time

- **Benchmark dataset** (6 structures from tiny to xlarge)
  - 1CRN (327 atoms), 1UBQ (602), 1A0Q (3,183), 3HHB (4,384), 1AON (58,674), 4V6X (237,685)
  - `benchmarks/inputs_protor/` - Pre-generated inputs with ProtOr radii
  - `scripts/data/generate_protor.py` - Generate inputs with ProtOr radii
  - `scripts/benchmark.py` - Unified benchmark comparing Zig vs FreeSASA

- **Lee-Richards algorithm** (`--algorithm=lr`)
  - Slice-based method with exact arc integration
  - `--n-slices=N` option (default: 20)
  - Multi-threading and SIMD support
  - 1.1x-1.7x faster than FreeSASA C

- **Atom classifier module** with CLI integration
  - `classifier.zig` - Core data structures, element-based radius guessing, ClassifierType enum
  - `classifier_naccess.zig` - NACCESS-compatible built-in classifier
  - `classifier_protor.zig` - ProtOr classifier (hybridization-based, Tsai et al. 1999)
  - `classifier_oons.zig` - OONS classifier (older FreeSASA default)
  - O(1) compile-time hash lookup using `StaticStringMap`
  - Support for 20 standard amino acids + SEC/MSE/PYL/ASX/GLX
  - Support for RNA/DNA nucleotides (A, C, G, I, T, U, DA, DC, DG, DI, DT, DU)
  - ANY fallback for backbone atoms (NACCESS/OONS)
  - Element-based radius guessing from atom names
  - `--classifier=naccess|protor|oons` - Use built-in classifier
  - `--config=FILE` - Use custom config file (FreeSASA format)

- Extended input format with optional `residue` and `atom_name` fields

### Changed

- **Updated benchmark**: Now compares against FreeSASA C (native binary) instead of Python
  - SR: 1.2x-2.3x faster than FreeSASA C
  - LR: 1.1x-1.7x faster than FreeSASA C
- **Scripts reorganization**
  - Created `scripts/data/` subdirectory for data preparation scripts
  - Renamed: `benchmark_all.py` → `benchmark.py`, `validate_accuracy.py` → `validate.py`
  - Moved data scripts to `scripts/data/` with shorter names
- **Project renamed**: `freesasa-zig` → `zsasa` (repository, binary, Python package)
- **Default ProtOr classifier** for PDB/mmCIF input (no `--classifier` flag needed)
- **Internal**: Refactored `AtomInput.r` from `[]const f64` to `[]f64` to properly support classifier mutations
- **Internal**: `FixedString5` for mmCIF 5-character `comp_id` (modified residue support)
- Added LICENSE (MIT) and CONTRIBUTING.md
- Enabled GitHub Actions CI/CD (format, build, test, Python)

## [0.0.5] - 2025-01-23

### Added

- **Extended CLI options**
  - `--probe-radius=R` - Configure probe radius (default: 1.4 Å)
  - `--n-points=N` - Configure test points per atom (default: 100)
  - `--quiet` / `-q` - Suppress progress output
  - `--help` / `-h` - Show help message
  - `--version` / `-V` - Show version

- **Output format options**
  - `--format=json` - Pretty-printed JSON (default)
  - `--format=compact` - Single-line JSON
  - `--format=csv` - CSV format with header

- **Input validation**
  - `--validate` - Validate input without calculation
  - Array length consistency check
  - Coordinate finiteness check (NaN/Inf detection)
  - Radius range validation (positive, ≤ 100 Å)
  - Detailed error messages with atom index and value

## [0.0.4] - 2025-01-22

### Added

- Multi-threading support with configurable thread count
- `--threads=N` CLI option (auto-detect by default)
- Generic thread pool implementation with work-stealing

### Changed

- Default execution mode is now multi-threaded
- 6.4x faster than FreeSASA (Python) on 3,183 atoms

## [0.0.3] - 2025-01-21

### Added

- SIMD optimization using `@Vector(4, f64)` for batch distance calculations
- 4.5x faster than FreeSASA (Python) single-threaded

## [0.0.2] - 2025-01-20

### Added

- Neighbor list optimization for O(N) neighbor lookup
- Spatial hashing with configurable cell size
- 3.9x faster than FreeSASA (Python) single-threaded

### Changed

- Reduced algorithmic complexity from O(N²) to O(N)

## [0.0.1] - 2025-01-19

### Added

- Initial Shrake-Rupley algorithm implementation
- Golden Section Spiral test point generation
- JSON input/output format
- Basic CLI with input/output file arguments
- Python scripts for structure conversion and benchmarking
  - `cif_to_input_json.py` - Convert mmCIF/PDB to input JSON
  - `calc_reference_sasa.py` - Generate reference SASA
  - `benchmark.py` - Performance benchmarking

[Unreleased]: https://github.com/N283T/zsasa/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/N283T/zsasa/compare/v0.1.3...v0.2.0
[0.1.3]: https://github.com/N283T/zsasa/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/N283T/zsasa/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/N283T/zsasa/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/N283T/zsasa/compare/v0.0.5...v0.1.0
[0.0.5]: https://github.com/N283T/zsasa/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/N283T/zsasa/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/N283T/zsasa/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/N283T/zsasa/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/N283T/zsasa/releases/tag/v0.0.1

# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

- Removed Windows from CI matrix (WSL recommended for Windows users)
- **Updated benchmark**: Now compares against FreeSASA C (native binary) instead of Python
  - SR: 1.2x-2.3x faster than FreeSASA C
  - LR: 1.1x-1.7x faster than FreeSASA C
- **Scripts reorganization**
  - Created `scripts/data/` subdirectory for data preparation scripts
  - Renamed: `benchmark_all.py` → `benchmark.py`, `validate_accuracy.py` → `validate.py`
  - Moved data scripts to `scripts/data/` with shorter names
- **Internal**: Refactored `AtomInput.r` from `[]const f64` to `[]f64` to properly support classifier mutations
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

[Unreleased]: https://github.com/N283T/zsasa/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/N283T/zsasa/compare/v0.0.5...v0.1.0
[0.0.5]: https://github.com/N283T/zsasa/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/N283T/zsasa/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/N283T/zsasa/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/N283T/zsasa/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/N283T/zsasa/releases/tag/v0.0.1

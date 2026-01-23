# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Lee-Richards algorithm** (`--algorithm=lr`)
  - Slice-based method with exact arc integration
  - `--n-slices=N` option (default: 20)
  - Multi-threading and SIMD support
  - 2.0x faster than FreeSASA (Python)

- **Atom classifier module** (library only, CLI integration planned)
  - `classifier.zig` - Core data structures, element-based radius guessing
  - `classifier_naccess.zig` - NACCESS-compatible built-in classifier
  - O(1) compile-time hash lookup using `StaticStringMap`
  - Support for 20 standard amino acids + SEC/MSE
  - Support for RNA/DNA nucleotides (A, C, G, I, T, U, DA, DC, DG, DI, DT, DU)
  - ANY fallback for backbone atoms
  - Element-based radius guessing from atom names

- Extended input format with optional `residue` and `atom_name` fields

### Changed

- Removed Windows from CI matrix (WSL recommended for Windows users)

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

[Unreleased]: https://github.com/N283T/freesasa-zig/compare/v0.0.5...HEAD
[0.0.5]: https://github.com/N283T/freesasa-zig/compare/v0.0.4...v0.0.5
[0.0.4]: https://github.com/N283T/freesasa-zig/compare/v0.0.3...v0.0.4
[0.0.3]: https://github.com/N283T/freesasa-zig/compare/v0.0.2...v0.0.3
[0.0.2]: https://github.com/N283T/freesasa-zig/compare/v0.0.1...v0.0.2
[0.0.1]: https://github.com/N283T/freesasa-zig/releases/tag/v0.0.1

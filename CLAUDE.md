# Project Rules

zsasa: High-performance SASA (Solvent Accessible Surface Area) calculation in Zig.

## Quick Reference

```bash
# Build
zig build -Doptimize=ReleaseFast

# Test
zig build test

# Format
zig fmt src/

# Python
cd python && uv pip install -e ".[dev]" && pytest tests/ -v
```

## Build

**ALWAYS use ReleaseFast for Zig builds:**

```bash
zig build -Doptimize=ReleaseFast
```

Debug builds are 3-10x slower. Use ReleaseFast for:
- Benchmarking
- Performance testing
- Production builds

Output:
- `zig-out/bin/zsasa` - CLI binary
- `zig-out/lib/libzsasa.{so,dylib}` - Shared library for Python

## Testing

### Zig Tests

```bash
zig build test                    # Run all tests
zig build test 2>&1 | head -100   # Truncate verbose output
```

### Python Tests

```bash
cd python
uv pip install -e ".[dev]"        # Install with dev dependencies
pytest tests/ -v                  # Run all tests
pytest tests/test_sasa.py -v      # Run specific test file
```

## Code Style

### Zig

```bash
zig fmt src/                      # Format all source files
zig fmt --check src/              # Check formatting (CI)
```

### Python

```bash
cd python
ruff check zsasa/          # Lint
ruff format zsasa/         # Format
ruff check zsasa/ --fix    # Auto-fix lint issues
```

## Project Structure

```
src/                    # Zig source code
├── main.zig           # CLI entry point
├── shrake_rupley.zig  # SR algorithm
├── lee_richards.zig   # LR algorithm
├── simd.zig           # SIMD optimizations
└── ...

python/                 # Python bindings
├── zsasa/      # Python package
├── tests/             # Python tests
└── pyproject.toml     # Package config

docs/                   # Documentation
├── benchmark/         # Benchmark methodology & results
├── ja/                # Japanese translations
└── *.md               # English docs

benchmarks/             # Benchmark infrastructure
├── scripts/           # Benchmark scripts (PEP 723)
├── results/           # Output data & plots
└── external/          # Comparison tools (FreeSASA, RustSASA)
```

## Benchmark Scripts

**DO NOT run benchmarks.** The user runs benchmarks manually. Only modify benchmark scripts when requested.

Located in `benchmarks/scripts/`. All scripts use PEP 723 metadata (typer + rich).

### Running Benchmarks

```bash
# Single-file mode
./benchmarks/scripts/bench.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/bench.py --tool freesasa --algorithm lr --threads 1-10

# With stratified sampling
./benchmarks/scripts/build_index.py benchmarks/inputs
./benchmarks/scripts/sample.py benchmarks/inputs/index.json \
    --target 100000 --seed 42 -o benchmarks/samples/stratified_100k.json
./benchmarks/scripts/bench.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_100k.json

# Batch mode (hyperfine-based, RustSASA paper style)
./benchmarks/scripts/bench_batch.py -i benchmarks/UP000000625_83333_ECOLI_v6/pdb -n ecoli --threads 8
```

### Analysis

```bash
./benchmarks/scripts/analyze.py summary      # Summary tables
./benchmarks/scripts/analyze.py all          # Generate all plots
./benchmarks/scripts/analyze.py scatter      # Atoms vs time
./benchmarks/scripts/analyze.py threads      # Thread scaling
./benchmarks/scripts/analyze.py grid         # Speedup by size/threads
./benchmarks/scripts/analyze.py validation   # SASA validation
./benchmarks/scripts/analyze.py samples      # Per-bin sample plots
./benchmarks/scripts/analyze.py large        # Large structure (100k+) analysis
./benchmarks/scripts/analyze.py efficiency   # Parallel efficiency
./benchmarks/scripts/analyze.py export-csv   # Export to CSV
```

### MD Trajectory Analysis

```bash
./benchmarks/scripts/analyze_md.py summary --name 6sup_R1   # Summary table
./benchmarks/scripts/analyze_md.py all --name 6sup_R1       # All plots + summary
./benchmarks/scripts/analyze_md.py bar --name 6sup_R1       # Tool comparison bar chart
./benchmarks/scripts/analyze_md.py memory --name 6sup_R1    # Memory usage
```

### Other Scripts

```bash
./benchmarks/scripts/generate_json.py /path/to/cif /path/to/output  # CIF→JSON
```

## CI/CD

GitHub Actions runs on push/PR to main:

| Job | Purpose |
|-----|---------|
| fmt | `zig fmt --check src/` |
| build | Build & test on Linux/macOS |
| python | Python tests & lint |
| validate | Integration tests with real data |

CI skips for: `*.md`, `docs/**`, `plans/**`, `benchmarks/**`, `LICENSE`

## Documentation

| Doc | Content |
|-----|---------|
| `docs/cli.md` | CLI reference |
| `docs/python-api/` | Python API |
| `docs/algorithm.md` | SR/LR algorithms |
| `docs/architecture.md` | Code structure |
| `docs/optimizations.md` | Performance techniques |
| `docs/benchmark/` | Benchmark methodology & results |
| `docs/ja/` | Japanese translations |

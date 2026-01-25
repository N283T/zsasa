# Project Rules

## Build

**ALWAYS use ReleaseFast for Zig builds:**

```bash
zig build -Doptimize=ReleaseFast
```

Never use `zig build` without `-Doptimize=ReleaseFast` - debug builds are 3-10x slower and will give misleading benchmark results.

## Benchmark Scripts

Benchmark scripts use typer + rich with PEP 723 metadata. Located in `benchmarks/scripts/`.

```bash
# Run benchmark (single tool, single algorithm)
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/run.py --tool freesasa --algorithm lr --threads 1-10

# Analyze results
./benchmarks/scripts/analyze.py summary    # Show summary tables
./benchmarks/scripts/analyze.py validate   # Validate SASA values
./benchmarks/scripts/analyze.py plot       # Generate graphs
./benchmarks/scripts/analyze.py all        # All of the above

# Generate JSON input from CIF files
./benchmarks/scripts/generate_json.py /path/to/cif /path/to/output

# Lint/Type check
ruff check benchmarks/scripts/*.py
```

## Benchmark

Before benchmarking, ensure:
1. Zig binary is built with ReleaseFast
2. FreeSASA C is built with thread support

```bash
zig build -Doptimize=ReleaseFast
./benchmarks/scripts/run.py --tool zig --algorithm sr
```

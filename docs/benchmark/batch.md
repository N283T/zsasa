# Batch Processing Benchmarks

Throughput benchmarks for processing multiple files in parallel. This measures total wall-clock time for processing a large set of structures.

> **Note**: This is supplementary to the main [single-file benchmarks](results.md) which measure per-structure performance.

## Overview

Unlike single-file benchmarks that measure SASA calculation time per structure, batch benchmarks measure **total throughput** when processing many files concurrently.

| Mode | Measures | Use Case |
|------|----------|----------|
| Single-file | Per-structure SASA time | Algorithm comparison |
| **Batch** | Total wall-clock time | Production throughput |

## Batch Implementation

Each tool uses its native batch processing capability:

| Tool | Implementation |
|------|----------------|
| Zig | Native directory input (auto-detected) |
| Rust | `--json-dir` flag with rayon parallelism |
| FreeSASA | Shell script wrapper with background jobs |

## Running Batch Benchmarks

### Basic Usage

```bash
# Run batch benchmark
./benchmarks/scripts/run_batch.py --tool zig --algorithm sr --threads 1,4,8

# With sample file
./benchmarks/scripts/run_batch.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_1k.json \
    --threads 1,4,8 --runs 3
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--tool` | Tool: zig, rust, freesasa | (required) |
| `--algorithm` | Algorithm: sr, lr | sr |
| `--threads` | Thread counts | 1,4,8 |
| `--runs` | Runs per configuration | 3 |
| `--precision` | f32 or f64 (Zig only) | f64 |
| `--parallelism` | file, atom, pipeline (Zig only) | file |

### Parallelism Modes (Zig only)

| Mode | Description |
|------|-------------|
| `file` | Process multiple files in parallel |
| `atom` | Parallelize within single file (default for single-file) |
| `pipeline` | Combination of file and atom parallelism |

## Analysis

```bash
# Summary table
./benchmarks/scripts/analyze_batch.py summary

# Generate comparison plots
./benchmarks/scripts/analyze_batch.py plot

# Both
./benchmarks/scripts/analyze_batch.py all
```

## Output

Results are saved to:

```
benchmarks/results/batch_{tool}_{algorithm}/
├── config.json   # System info and parameters
└── results.csv   # Benchmark results (per-run timing)
```

## Results

Processing all **238,124 PDB structures** with 10 threads:

| Tool | Precision | Total Time | Throughput | vs Rust |
|------|-----------|------------|------------|---------|
| **Zig** | f32 | **724.6s** | **328.6 files/s** | **+7%** |
| Zig | f64 | 747.9s | 318.4 files/s | +4% |
| Rust | f32 | 774.5s | 307.5 files/s | baseline |

> FreeSASA C excluded as single-file benchmarks show it is slower than Rust.

## Notes

1. **Different from single-file benchmarks**: Batch benchmarks include file I/O overhead
2. **Parallelism strategy matters**: File-level parallelism vs atom-level parallelism
3. **Memory usage**: Batch mode may use more memory for concurrent file processing

## Related Documents

- [methodology.md](methodology.md) - Benchmark methodology
- [results.md](results.md) - Single-file benchmark results

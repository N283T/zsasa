# Batch Processing Benchmarks

Throughput benchmarks for processing a complete dataset using hyperfine timing.

> **Note**: This measures total wall-clock time for processing an entire directory of PDB files.

## Overview

Batch benchmarks measure **total throughput** when processing many files with native multi-threading.

| Benchmark | Measures | Tool |
|-----------|----------|------|
| [Single-file](results.md) | Per-structure SASA time | run.py |
| **Batch** | Total wall-clock time | batch_bench.py |
| [E. coli](ecoli-proteome.md) | Proteome throughput | batch_bench.py |

## Running Batch Benchmarks

### Basic Usage

```bash
# E. coli proteome benchmark
./benchmarks/scripts/batch_bench.py \
  --input benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  --name ecoli \
  --runs 5 --threads 8

# Custom dataset
./benchmarks/scripts/batch_bench.py \
  -i /path/to/pdb_dir \
  -n my-benchmark \
  --runs 3

# Single tool only
./benchmarks/scripts/batch_bench.py \
  -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  -n ecoli \
  --tool zig --runs 3

# Dry run (show commands)
./benchmarks/scripts/batch_bench.py \
  -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
  -n test \
  --dry-run
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input`, `-i` | Input directory (PDB files) | (required) |
| `--name`, `-n` | Benchmark name | (required) |
| `--threads`, `-T` | Thread count for multi-threaded runs | 8 |
| `--runs`, `-r` | Number of benchmark runs | 3 |
| `--warmup`, `-w` | Number of warmup runs | 3 |
| `--tool`, `-t` | Tool: zig, freesasa, rustsasa, all | all |
| `--output`, `-o` | Output directory | results/batch/<name> |
| `--dry-run` | Show commands without running | false |

## Benchmark Configurations

For each tool, the following configurations are tested:

| Tool | Configurations |
|------|----------------|
| zsasa (Zig) | f64 Nt, f64 1t, f32 Nt, f32 1t |
| FreeSASA | 1t only (sequential) |
| RustSASA | Nt, 1t |

Where N = `--threads` value (default 8).

## Output

Results are saved to:

```
benchmarks/results/batch/<name>/
├── bench_zsasa_f64_8t.json
├── bench_zsasa_f64_1t.json
├── bench_zsasa_f32_8t.json
├── bench_zsasa_f32_1t.json
├── bench_freesasa_1t.json
├── bench_rustsasa_8t.json
└── bench_rustsasa_1t.json
```

Each JSON file contains hyperfine output:

```json
{
  "results": [{
    "command": "...",
    "mean": 4.614,
    "stddev": 0.066,
    "min": 4.511,
    "max": 4.688,
    "times": [4.51, 4.65, ...],
    "memory_usage_byte": [...],
    "exit_codes": [...]
  }]
}
```

## Methodology

Uses [hyperfine](https://github.com/sharkdp/hyperfine) for timing:

1. Warmup runs (default 3) to eliminate cold-start effects
2. Multiple timed runs (default 3-5) for statistical reliability
3. Reports mean, stddev, min, max

This methodology matches the [RustSASA paper](https://github.com/OWissett/rustsasa) for direct comparison.

## Results

See [E. coli Proteome Benchmark](ecoli-proteome.md) for detailed results.

**Key findings (E. coli proteome, 4,370 structures):**

| Tool | Threads | Time | vs FreeSASA |
|------|--------:|-----:|------------:|
| zsasa f32 | 8 | 4.4s | **10.4x** |
| zsasa f64 | 8 | 4.6s | **9.9x** |
| RustSASA | 8 | 5.5s | 8.3x |
| FreeSASA | 1 | 45.8s | baseline |

## Related Documents

- [methodology.md](methodology.md) - Benchmark methodology
- [results.md](results.md) - Single-file benchmark results
- [ecoli-proteome.md](ecoli-proteome.md) - E. coli proteome benchmark details

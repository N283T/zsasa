# Benchmark Methodology

Measurement methods and execution procedures for benchmarks. See [results.md](results.md) for results.

## SASA-Only Timing

For fair comparison, we measure **SASA calculation time only**. File I/O is excluded.

```
Total time = File I/O + SASA calculation + Output
                        ^^^^^^^^^^^^^^^^
                        Only this is measured
```

Measurement method for each implementation:

| Implementation | Method |
|----------------|--------|
| zsasa | Internal measurement via `--timing` option (stderr output) |
| FreeSASA C | Patched binary outputs SASA calculation time to stderr |
| RustSASA | Patched binary outputs `SASA_TIME_US` to stderr |

## Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Algorithm | Shrake-Rupley | Supported by all implementations (LR: Zig/FreeSASA only) |
| n_points | 100 | Number of test points |
| probe_radius | 1.4 Å | Water molecule radius |
| Runs | 3 | Average value used |

## Stratified Sampling

Stratified sampling from all PDB structures (~238K):

| Bin | Range | Strategy |
|-----|-------|----------|
| 0-500 | 0 - 500 | Proportional allocation |
| 500-2k | 500 - 2,000 | Proportional allocation |
| 2k-10k | 2,000 - 10,000 | Proportional allocation |
| 10k-50k | 10,000 - 50,000 | Proportional allocation |
| 50k-200k | 50,000 - 200,000 | **All included** |
| 200k+ | 200,000+ | **All included** |

Rare large structures (50k+ atoms) are all included; the rest use proportional allocation.

## Running Benchmarks

### Setup

```bash
# Build Zig binary
zig build -Doptimize=ReleaseFast

# External tools setup (for comparison)
cd benchmarks/external
git clone https://github.com/N283T/freesasa-bench.git
cd freesasa-bench && ./configure --enable-threads && make && cd ..
git clone --recursive https://github.com/N283T/rustsasa-bench.git
cd rustsasa-bench && cargo build --release --features cli && cd ..
```

### Index & Sample Generation

```bash
# Create index (first time only)
./benchmarks/scripts/build_index.py benchmarks/inputs

# Check distribution
./benchmarks/scripts/sample.py benchmarks/inputs/index.json --analyze

# Generate sample
./benchmarks/scripts/sample.py benchmarks/inputs/index.json \
    --target 100000 --seed 42 \
    -o benchmarks/samples/stratified_100k.json
```

### Running Benchmarks

```bash
# Basic usage
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/run.py --tool freesasa --algorithm lr --threads 1-10

# With sample file
./benchmarks/scripts/run.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_100k.json \
    --threads 1,2,4,8,10

# Single run for quick testing
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1 --runs 1

# With f32 precision (Zig only)
./benchmarks/scripts/run.py --tool zig --algorithm sr --precision f32
```

### Analysis & Visualization

```bash
# Summary tables
./benchmarks/scripts/analyze.py summary

# Generate all plots
./benchmarks/scripts/analyze.py all

# Individual plot types
./benchmarks/scripts/analyze.py scatter      # Atoms vs time scatter
./benchmarks/scripts/analyze.py threads      # Thread scaling
./benchmarks/scripts/analyze.py grid         # Speedup grid by size/threads
./benchmarks/scripts/analyze.py validation   # SASA validation
./benchmarks/scripts/analyze.py samples      # Per-bin sample plots
./benchmarks/scripts/analyze.py large        # Large structure analysis
./benchmarks/scripts/analyze.py efficiency   # Parallel efficiency

# Export to CSV
./benchmarks/scripts/analyze.py export-csv
```

## Scripts

| Script | Purpose |
|--------|---------|
| `build_index.py` | Create atom count index from all input files |
| `sample.py` | Stratified sampling from index |
| `bench.py` | Run benchmarks (single-file mode) |
| `bench_batch.py` | Batch benchmarks (hyperfine-based) |
| `analyze.py` | Analyze results and generate plots |
| `generate_json.py` | Convert CIF/PDB to JSON format |

## Reproducibility

Same seed produces identical samples:

```bash
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o a.json
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o b.json
diff a.json b.json  # No differences
```

System information displayed during benchmark execution:

```
=== System Info ===
Platform: Darwin (arm64)
CPU: Apple M4
Cores: 10
Memory: 24.0 GB

=== Benchmark Config ===
Threads: 4
N-points: 100
Probe radius: 1.4 Å
Runs: 5
```

## Notes

1. **Initial runs are slow**: Due to file cache and warmup effects
2. **Thread count depends on CPU**: Optimal when matching physical core count
3. **External tools require patches**: SASA-only timing requires modified binaries

## Related Documents

- [results.md](results.md) - Benchmark results
- [batch.md](batch.md) - Batch processing benchmarks

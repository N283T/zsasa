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
| freesasa-zig | Internal measurement via `--timing` option (stderr output) |
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

### Single-File Mode

Process files sequentially to test intra-algorithm parallelization:

```bash
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/run.py --tool freesasa --algorithm lr --threads 1-10

# With sample file
./benchmarks/scripts/run.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_100k.json \
    --threads 1,2,4,8,10
```

### Batch Mode

Process multiple files in parallel to test throughput:

```bash
./benchmarks/scripts/run_batch.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_1k.json \
    --threads 1,4,8,10 --runs 3
```

| Tool | Batch Implementation |
|------|---------------------|
| Zig | Native directory input |
| Rust | `--json-dir` + rayon |
| FreeSASA | Shell script + background jobs |

### Analysis

```bash
./benchmarks/scripts/analyze.py summary    # Summary tables
./benchmarks/scripts/analyze.py validate   # SASA validation
./benchmarks/scripts/analyze.py plot       # Generate graphs
./benchmarks/scripts/analyze.py all        # All of the above

# CSV export
./benchmarks/scripts/analyze.py export_csv
```

## Scripts

| Script | Purpose |
|--------|---------|
| `build_index.py` | Create atom count index from all files |
| `sample.py` | Stratified sampling from index |
| `run.py` | Run benchmarks in single-file mode |
| `run_batch.py` | Run benchmarks in batch mode |
| `analyze.py` | Analyze results and generate graphs |
| `generate_json.py` | Convert CIF to JSON |

## Reproducibility

Same seed produces identical samples:

```bash
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o a.json
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o b.json
diff a.json b.json  # 差分なし
```

ベンチマーク実行時に表示されるシステム情報：

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

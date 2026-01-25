# External Benchmark Dependencies

Optional dependencies for running comparative benchmarks against FreeSASA C and RustSASA.

## Setup

### FreeSASA (C implementation)

```bash
git clone https://github.com/N283T/freesasa-bench.git
cd freesasa-bench
./configure --enable-threads
make
```

Binary: `freesasa-bench/src/freesasa`

### RustSASA (Rust implementation)

```bash
git clone https://github.com/N283T/rustsasa-bench.git
cd rustsasa-bench
cargo build --release --features cli
```

Binary: `rustsasa-bench/target/release/rust-sasa`

## Timing Patches

Both forks include timing patches that output SASA calculation time to stderr:

- **FreeSASA**: `SASA calculation time: X.XX ms`
- **RustSASA**: `SASA_TIME_US:XXXXX`

This allows fair SASA-only timing comparison (excluding file I/O).

## Running Benchmarks

After building both dependencies:

```bash
cd /path/to/freesasa-zig
./scripts/benchmark.py
```

The benchmark script automatically detects binaries in this directory.

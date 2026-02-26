# External Benchmark Dependencies

Optional dependencies for running comparative benchmarks.

## Directory Layout

```
external/
├── freesasa/        ← vanilla upstream (single-file bench via hyperfine)
├── freesasa-bench/  ← fork with freesasa_batch (batch bench)
├── lahuta/          ← as-is
└── rustsasa/        ← vanilla upstream (replaces rustsasa-bench)
```

## Setup

### FreeSASA — vanilla (single-file benchmark)

```bash
cd benchmarks/external
git clone https://github.com/mittinatten/freesasa.git
cd freesasa
./configure --enable-threads
make
```

Binary: `freesasa/src/freesasa`

### FreeSASA — fork with freesasa_batch (batch benchmark)

```bash
cd benchmarks/external
git clone https://github.com/N283T/freesasa-bench.git
cd freesasa-bench
./configure --enable-threads
make
cd src && make freesasa_batch
```

Binary: `freesasa-bench/src/freesasa_batch`

### RustSASA — vanilla (single-file & batch benchmark)

```bash
cd benchmarks/external
git clone --recursive https://github.com/mcisb/rustsasa.git
cd rustsasa
cargo build --release --features cli
```

Binary: `rustsasa/target/release/rust-sasa`

### Lahuta

```bash
cd benchmarks/external
git clone https://github.com/DominikSko/lahuta.git
cd lahuta
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Binary: `lahuta/build/cli/lahuta`

## Usage

Single-file benchmarks use **hyperfine** for wall-clock timing:

```bash
# SR benchmark (Shrake-Rupley)
./benchmarks/scripts/bench.py --tool zig_f64 --threads 1,4,8

# LR benchmark (Lee-Richards)
./benchmarks/scripts/bench_lr.py --tool zig_f64 --threads 1,4,8

# Batch benchmark
./benchmarks/scripts/bench_batch.py -i benchmarks/dataset/pdb -n test
```

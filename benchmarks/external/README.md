# External Benchmark Dependencies

Optional dependencies for running comparative benchmarks.

## Directory Layout

```
external/
├── freesasa/        ← vanilla upstream (single-file & batch build dependency)
├── freesasa_batch/  ← batch runner (source tracked, builds against vanilla freesasa)
├── lahuta/          ← as-is
└── rustsasa/        ← vanilla upstream
```

## Setup

### 1. FreeSASA (vanilla)

Required for both single-file benchmarks and building `freesasa_batch`.

```bash
cd benchmarks/external
git clone https://github.com/mittinatten/freesasa.git
cd freesasa
./configure --enable-threads
make
```

Binary: `freesasa/src/freesasa`

### 2. freesasa_batch (batch benchmark)

Builds against the vanilla FreeSASA library from step 1.

```bash
cd benchmarks/external/freesasa_batch
make
```

Binary: `freesasa_batch/freesasa_batch`

### 3. RustSASA (vanilla)

```bash
cd benchmarks/external
git clone --recursive https://github.com/mcisb/rustsasa.git
cd rustsasa
cargo build --release --features cli
```

Binary: `rustsasa/target/release/rust-sasa`

### 4. Lahuta

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

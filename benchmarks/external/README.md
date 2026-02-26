# External Benchmark Dependencies

Optional dependencies for running comparative benchmarks.

## Directory Layout

```
external/
├── bin/              ← symlinks to all tool binaries (created by setup.sh)
├── freesasa/         ← vanilla upstream (cloned by setup.sh)
├── freesasa_batch/   ← batch runner source (tracked, builds against vanilla freesasa)
├── lahuta/           ← cloned by setup.sh
├── rustsasa/         ← vanilla upstream (cloned by setup.sh)
└── setup.sh          ← one-command setup: clone, build, symlink
```

## Quick Start

```bash
cd benchmarks/external
./setup.sh
```

This will clone, build, and symlink all tool binaries into `bin/`.

To build a single tool:

```bash
./setup.sh freesasa
./setup.sh rustsasa
./setup.sh lahuta
./setup.sh freesasa_batch
./setup.sh zsasa          # symlink only (run 'zig build' first)
```

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

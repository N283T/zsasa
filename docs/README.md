# Documentation

Technical documentation for zsasa.

## Contents

| File | Description |
|------|-------------|
| [zig-api/](zig-api/) | Zig Library API Reference |
| [architecture.md](architecture.md) | Architecture overview |
| [algorithm.md](algorithm.md) | SASA algorithms (SR/LR) |
| [optimizations.md](optimizations.md) | Optimization techniques |
| [benchmark/](benchmark/) | Benchmark documentation |
| ├─ [single-file.md](benchmark/single-file.md) | Single-file benchmark results with plots |
| ├─ [batch.md](benchmark/batch.md) | Batch processing benchmarks |
| └─ [md.md](benchmark/md.md) | MD trajectory benchmarks |
| [cli.md](cli.md) | CLI Reference |
| [python-api/](python-api/) | Python API Reference |
| [xtc.md](xtc.md) | XTC trajectory file reader |
| [classifier.md](classifier.md) | Atom radius classifiers |
| [ci.md](ci.md) | CI/CD configuration |

## Benchmark Environment

Performance numbers in this documentation were measured on:

- **Hardware**: Apple M4 (10 cores: 4 performance + 6 efficiency), 32 GB RAM
- **Test data**: ~100k PDB structures (stratified sampling)
- **Zig version**: 0.15.2

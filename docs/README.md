# Documentation

Technical documentation for freesasa-zig.

## Contents

| File | Description |
|------|-------------|
| [architecture.md](architecture.md) | Architecture overview |
| [algorithm.md](algorithm.md) | SASA algorithms (SR/LR) |
| [optimizations.md](optimizations.md) | Optimization techniques |
| [benchmark/](benchmark/) | Benchmark documentation |
| ├─ [results.md](benchmark/results.md) | Large-scale benchmark results with plots |
| ├─ [methodology.md](benchmark/methodology.md) | Benchmark methodology |
| └─ [batch.md](benchmark/batch.md) | Batch processing benchmarks |
| [cpu-efficiency.md](cpu-efficiency.md) | CPU efficiency analysis (IPC, instructions) |
| [cli.md](cli.md) | CLI Reference |
| [python.md](python.md) | Python API Reference |
| [classifier.md](classifier.md) | Atom radius classifiers |
| [ci.md](ci.md) | CI/CD configuration |
| [ja/](ja/) | Japanese documentation |

## Benchmark Environment

Performance numbers in this documentation were measured on:

- **Hardware**: Apple M4 (10 cores: 4 performance + 6 efficiency), 32 GB RAM
- **Test data**: ~100k PDB structures (stratified sampling)
- **Zig version**: 0.15.2

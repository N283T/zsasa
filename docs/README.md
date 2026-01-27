# Documentation / ドキュメント

Technical documentation for freesasa-zig (Japanese).

freesasa-zig の技術ドキュメント（日本語）。

## Contents / 目次

| File | Description | 内容 |
|------|-------------|------|
| [architecture.md](architecture.md) | Architecture overview | アーキテクチャ概要 |
| [algorithm.md](algorithm.md) | SASA algorithms (SR/LR) | アルゴリズム詳解 |
| [optimizations.md](optimizations.md) | Optimization techniques | 最適化技術 |
| [benchmark.md](benchmark.md) | Benchmark methodology | ベンチマーク手法 |
| [benchmark-results.md](benchmark-results.md) | Large-scale benchmark results with plots | 大規模ベンチマーク結果 |
| [cpu-efficiency.md](cpu-efficiency.md) | CPU efficiency analysis (IPC, instructions) | CPU効率解析 |
| [cli-io.md](cli-io.md) | CLI, I/O, and analysis | CLI・入出力・解析 |
| [classifier.md](classifier.md) | Atom radius classifiers | 原子分類器詳解 |
| [ci.md](ci.md) | CI/CD configuration | CI/CD構成 |

## Benchmark Environment / ベンチマーク環境

Performance numbers in this documentation were measured on:

- **Hardware**: Apple M2 (8 cores)
- **Test data**: PDB 1A0Q (3,183 atoms)
- **Zig version**: 0.15.2

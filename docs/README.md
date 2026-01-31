# Documentation / ドキュメント

Technical documentation for freesasa-zig (Japanese).

freesasa-zig の技術ドキュメント（日本語）。

## Contents / 目次

| File | Description | 内容 |
|------|-------------|------|
| [architecture.md](architecture.md) | Architecture overview | アーキテクチャ概要 |
| [algorithm.md](algorithm.md) | SASA algorithms (SR/LR) | アルゴリズム詳解 |
| [optimizations.md](optimizations.md) | Optimization techniques | 最適化技術 |
| [benchmark/](benchmark/) | Benchmark documentation | ベンチマーク |
| ├─ [methodology.md](benchmark/methodology.md) | Benchmark methodology | 測定手法 |
| └─ [results.md](benchmark/results.md) | Large-scale benchmark results with plots | 大規模結果 |
| [cpu-efficiency.md](cpu-efficiency.md) | CPU efficiency analysis (IPC, instructions) | CPU効率解析 |
| [cli.md](cli.md) | CLI Reference (English) | CLIリファレンス（英語） |
| [cli-io.md](cli-io.md) | CLI, I/O (Japanese) | CLI・入出力（日本語） |
| [classifier.md](classifier.md) | Atom radius classifiers | 原子分類器詳解 |
| [ci.md](ci.md) | CI/CD configuration | CI/CD構成 |

## Benchmark Environment / ベンチマーク環境

Performance numbers in this documentation were measured on:

- **Hardware**: Apple M4 (10 cores: 4 performance + 6 efficiency), 32 GB RAM
- **Test data**: ~100k PDB structures (stratified sampling)
- **Zig version**: 0.15.2

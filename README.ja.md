# zsasa

[![CI](https://github.com/N283T/zsasa/actions/workflows/ci.yml/badge.svg)](https://github.com/N283T/zsasa/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/zsasa?color=blue)](https://pypi.org/project/zsasa/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Zig](https://img.shields.io/badge/Zig-0.15.2+-f7a41d?logo=zig&logoColor=white)](https://ziglang.org/)
[![Python](https://img.shields.io/badge/Python-3.11+-3776ab?logo=python&logoColor=white)](https://www.python.org/)

[English](README.md) | 日本語

Zig で実装された高性能 Solvent Accessible Surface Area (SASA) 計算ツール。
FreeSASA C より**最大3倍高速**、f64 精度を維持。

**[ドキュメント](https://n283t.github.io/zsasa/)**

## 特徴

- **2つのアルゴリズム**: Shrake-Rupley（高速）と Lee-Richards（高精度）
- **複数の入力形式**: mmCIF, PDB, JSON
- **解析機能**: 残基単位集計、RSA、極性/非極性分類
- **高性能**: SIMD最適化、マルチスレッド、近傍リスト O(N)
- **クロスプラットフォーム**: Linux, macOS, Windows（`pip install zsasa` でビルド済みホイール利用可）
- **Python バインディング**: NumPy 連携、BioPython/Biotite/Gemmi 対応
- **MD トラジェクトリ解析**: ネイティブ XTC リーダー、MDTraj / MDAnalysis 連携

## クイックスタート

### Python

```bash
pip install zsasa
# or
uv add zsasa
```

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
```

### CLI

[Zig 0.15.2+](https://ziglang.org/download/) が必要です。

```bash
git clone https://github.com/N283T/zsasa.git
cd zsasa && zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa calc structure.cif output.json
```

## ベンチマーク

### 単一ファイル性能

| スピードアップ (threads=10) | スレッドスケーリング (50k+ atoms) |
|:--------------------:|:----------------------------:|
| ![Speedup](benchmarks/results/plots/large/speedup_bar.png) | ![Thread Scaling](benchmarks/results/plots/large/speedup_by_threads.png) |

**主要結果 (50k+ atoms, threads=10):**
- FreeSASA/RustSASA 比で**2.3倍**の中央値スピードアップ (SR)
- スレッド数増加に伴いスピードアップ向上（優れた並列効率）

> **注**: zsasa/FreeSASA は f64、RustSASA は f32 を使用。

### MD トラジェクトリ性能

実際の MD トラジェクトリデータで mdsasa-bolt (RustSASA) より**4.3倍高速**。

![MD Trajectory Benchmark](benchmarks/results/md/6sup_A_analysis/plots/bar.png)

*33,377 atoms, 1,001 frames, n_points=100*

詳細は[ベンチマーク結果](https://n283t.github.io/zsasa/benchmarks/single-file)を参照。

## ドキュメント

完全なドキュメントは **[n283t.github.io/zsasa](https://n283t.github.io/zsasa/)** で公開しています。

| セクション | 説明 |
|---------|-------------|
| [Getting Started](https://n283t.github.io/zsasa/getting-started) | インストールと最初の計算 |
| [CLI Reference](https://n283t.github.io/zsasa/cli) | CLI オプションと出力形式 |
| [Python API](https://n283t.github.io/zsasa/python-api) | コア API、インテグレーション、MD トラジェクトリ |
| [Benchmarks](https://n283t.github.io/zsasa/benchmarks/single-file) | 方法論と結果 |

## コントリビュート

開発セットアップとガイドラインは [CONTRIBUTING.md](CONTRIBUTING.md) を参照。

## ライセンス

[MIT](LICENSE)

## 参考文献

- Shrake, A.; Rupley, J. A. Environment and Exposure to Solvent of Protein Atoms. *J. Mol. Biol.* 1973, 79(2), 351-371. [doi:10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
- Lee, B.; Richards, F. M. The Interpretation of Protein Structures: Estimation of Static Accessibility. *J. Mol. Biol.* 1971, 55(3), 379-400. [doi:10.1016/0022-2836(71)90324-x](https://doi.org/10.1016/0022-2836(71)90324-x)
- Mitternacht, S. FreeSASA: An Open Source C Library for Solvent Accessible Surface Area Calculations. *F1000Res.* 2016, 5, 189. [doi:10.12688/f1000research.7931.1](https://doi.org/10.12688/f1000research.7931.1)
- Campbell, M. J. RustSASA: A Rust Crate for Accelerated Solvent Accessible Surface Area Calculations. *J. Open Source Softw.* 2026, 11(117), 9537. [doi:10.21105/joss.09537](https://doi.org/10.21105/joss.09537)

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

- **2つのアルゴリズム**: Shrake-Rupley（高速）と Lee-Richards（高精度）、ビットマスク LUT 最適化対応
- **複数の入力形式**: mmCIF, PDB, JSON
- **MD トラジェクトリ解析**: ネイティブ XTC/DCD リーダー、MDTraj / MDAnalysis 連携
- **バッチ処理**: プロテオームスケールのディレクトリ一括処理
- **解析機能**: 残基単位集計、RSA、極性/非極性分類
- **高性能**: SIMD最適化、マルチスレッド、f64/f32 選択可能
- **クロスプラットフォーム**: Linux, macOS, Windows（`pip install zsasa` でビルド済みホイール利用可）
- **Python バインディング**: NumPy 連携、Gemmi/BioPython/Biotite 対応

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

### 単一ファイル（2,013構造、n_points=100、threads=10）

| 指標 | vs FreeSASA | vs RustSASA |
|------|-------------|-------------|
| 中央値スピードアップ | **1.88x** | **1.84x** |
| スレッドスケーリング (t=1→10) | 2.71x | 1.39x |

### バッチ（E. coli K-12 プロテオーム、4,370構造）

| ツール | 時間 | メモリ |
|--------|------|--------|
| zsasa bitmask (f32) | **1.42s** | 43 MB |
| Lahuta bitmask | 2.01s | 291 MB |
| RustSASA | 5.24s | 169 MB |

### MD トラジェクトリ

mdsasa-bolt (RustSASA) 比で**82〜96倍少ないメモリ使用量**。10Kフレームのデータセットでは mdsasa-bolt が2時間超でタイムアウトする中、zsasa は38秒で完了。

詳細は[ベンチマーク結果](https://n283t.github.io/zsasa/docs/benchmarks)を参照。

## ドキュメント

完全なドキュメントは **[n283t.github.io/zsasa](https://n283t.github.io/zsasa/)** で公開しています。

| セクション | 説明 |
|---------|-------------|
| [Getting Started](https://n283t.github.io/zsasa/docs/getting-started) | インストールと最初の計算 |
| [Comparison](https://n283t.github.io/zsasa/docs/comparison) | FreeSASA, RustSASA, Lahuta との比較 |
| [CLI Reference](https://n283t.github.io/zsasa/docs/cli/commands) | CLI オプションと出力形式 |
| [Python API](https://n283t.github.io/zsasa/docs/python-api) | コア API、インテグレーション、MD トラジェクトリ |
| [Benchmarks](https://n283t.github.io/zsasa/docs/benchmarks) | 性能・精度ベンチマーク |

## コントリビュート

開発セットアップとガイドラインは [CONTRIBUTING.md](CONTRIBUTING.md) を参照。

## ライセンス

[MIT](LICENSE)

## 参考文献

- Shrake, A.; Rupley, J. A. Environment and Exposure to Solvent of Protein Atoms. *J. Mol. Biol.* 1973, 79(2), 351-371. [doi:10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
- Lee, B.; Richards, F. M. The Interpretation of Protein Structures: Estimation of Static Accessibility. *J. Mol. Biol.* 1971, 55(3), 379-400. [doi:10.1016/0022-2836(71)90324-x](https://doi.org/10.1016/0022-2836(71)90324-x)
- Mitternacht, S. FreeSASA: An Open Source C Library for Solvent Accessible Surface Area Calculations. *F1000Res.* 2016, 5, 189. [doi:10.12688/f1000research.7931.1](https://doi.org/10.12688/f1000research.7931.1)
- Campbell, M. J. RustSASA: A Rust Crate for Accelerated Solvent Accessible Surface Area Calculations. *J. Open Source Softw.* 2026, 11(117), 9537. [doi:10.21105/joss.09537](https://doi.org/10.21105/joss.09537)

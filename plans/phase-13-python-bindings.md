# Phase 13: Python Bindings

## Goal

ZigライブラリをPythonから利用できるようにし、科学計算ワークフローへの統合を容易にする。

## 優先度: 低

## アプローチ

### Option A: C ABI + ctypes

ZigのC ABI互換性を利用してPythonから呼び出す。

```zig
// Zig側 (lib.zig)
export fn freesasa_calculate(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    r: [*]const f64,
    n: usize,
    result: [*]f64,
) callconv(.C) c_int {
    // ...
}
```

```python
# Python側
import ctypes

lib = ctypes.CDLL("libfreesasa_zig.so")
lib.freesasa_calculate.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # x
    ctypes.POINTER(ctypes.c_double),  # y
    ctypes.POINTER(ctypes.c_double),  # z
    ctypes.POINTER(ctypes.c_double),  # r
    ctypes.c_size_t,                   # n
    ctypes.POINTER(ctypes.c_double),  # result
]
```

### Option B: PyO3/Maturin (via C)

より洗練されたPythonインターフェース。

## Tasks

### Phase 13.1: C ABI Export

- [x] 共有ライブラリビルド設定
- [x] C互換APIの設計
- [x] エラーハンドリング（戻り値）

**API設計:**
```zig
// エラーコード
pub const FREESASA_OK: c_int = 0;
pub const FREESASA_ERROR_INVALID_INPUT: c_int = -1;
pub const FREESASA_ERROR_OUT_OF_MEMORY: c_int = -2;

// メイン関数
export fn freesasa_calculate_sasa(
    x: [*]const f64,
    y: [*]const f64,
    z: [*]const f64,
    r: [*]const f64,
    n_atoms: usize,
    n_points: u32,
    probe_radius: f64,
    n_threads: usize,
    atom_areas: [*]f64,
    total_area: *f64,
) callconv(.C) c_int;

// バージョン取得
export fn freesasa_version() callconv(.C) [*:0]const u8;
```

### Phase 13.2: Python Wrapper

- [x] ctypesラッパー実装
- [x] NumPy配列サポート
- [x] Pythonic API設計

**Python API:**
```python
import numpy as np
from freesasa_zig import calculate_sasa

# NumPy配列で入力
coords = np.array([[x1, y1, z1], [x2, y2, z2], ...])
radii = np.array([r1, r2, ...])

# SASA計算
result = calculate_sasa(
    coords,
    radii,
    n_points=100,
    probe_radius=1.4,
    n_threads=0  # auto
)

print(f"Total SASA: {result.total_area}")
print(f"Atom areas: {result.atom_areas}")
```

### Phase 13.3: Package Distribution

- [x] PyPI用パッケージング (pyproject.toml)
- [ ] wheel ビルド（manylinux, macosx）- 将来の課題
- [x] ドキュメント

## Files

| File | Action |
|------|--------|
| `src/c_api.zig` | CREATE |
| `python/freesasa_zig/__init__.py` | CREATE |
| `python/freesasa_zig/core.py` | CREATE |
| `python/pyproject.toml` | CREATE |
| `python/tests/test_sasa.py` | CREATE |
| `scripts/benchmark_python.py` | CREATE |
| `build.zig` | MODIFY (shared lib) |

## Success Criteria

- [x] Pythonから計算が実行できる
- [x] NumPy配列が直接渡せる
- [x] pip install可能なパッケージ

## Benchmark Results

Python bindings vs FreeSASA Python bindings:

| PDB | Atoms | Zig SR | FS SR | LR Speedup |
|-----|-------|--------|-------|------------|
| 1crn | 327 | 1.9ms | 0.7ms | 0.8x |
| 1ubq | 602 | 2.3ms | 1.2ms | 1.2x |
| 1a0q | 3183 | 10.6ms | 7.5ms | 1.4x |
| 3hhb | 4384 | 14.4ms | 11.0ms | 1.4x |
| 1aon | 58674 | 201ms | 163ms | 1.4x |

- SR: FreeSASA is slightly faster (C extension optimized)
- LR: Zig is 1.2-1.4x faster
- Accuracy: 0.00 Å² difference (identical results)

---
- [x] **DONE** - Phase 13 complete

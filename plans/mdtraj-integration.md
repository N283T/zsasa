# MDTraj Integration Plan

## Goal

freesasa-zigをMDTrajと統合し、MD trajectoryのSASA計算を高速化する。

## Background

### MDTraj現状
- `md.shrake_rupley(traj)` でSASA計算
- 内部はC++実装（OpenMP対応だがドキュメントなし）
- Biotiteより約2倍遅い

### freesasa-zigの強み
- SIMD最適化（AVX2/NEON）
- マルチスレッド対応
- Lee-Richards/Shrake-Rupley両対応

### MDの特性（活用できるポイント）
- トポロジー固定（原子数、半径は変わらない）
- フレーム間で独立（embarrassingly parallel）
- 座標のみ変化

## API Design

### Option A: MDTraj互換API

```python
from freesasa_zig import compute_sasa_trajectory

# MDTraj Trajectoryを直接受け取る
sasa = compute_sasa_trajectory(
    traj,                    # mdtraj.Trajectory
    algorithm='sr',          # 'sr' or 'lr'
    n_threads=None,          # None = auto
    probe_radius=1.4,        # Angstrom
    n_points=100,            # sphere points (SR only)
    mode='atom',             # 'atom' or 'residue'
)
# Returns: np.ndarray shape (n_frames, n_atoms) or (n_frames, n_residues)
```

### Option B: NumPy配列API（汎用）

```python
from freesasa_zig import compute_sasa_batch

# 生のNumPy配列を受け取る（MDAnalysisにも対応可能）
sasa = compute_sasa_batch(
    coordinates,             # np.ndarray (n_frames, n_atoms, 3)
    radii,                   # np.ndarray (n_atoms,)
    algorithm='sr',
    n_threads=None,
    probe_radius=1.4,
    n_points=100,
)
# Returns: np.ndarray shape (n_frames, n_atoms)
```

### Recommendation

**両方実装**:
- `compute_sasa_batch()` - 低レベルAPI（NumPy配列）
- `compute_sasa_trajectory()` - 高レベルAPI（MDTraj wrapper）

## Implementation Phases

### Phase 1: Zig側バッチ処理API

現状の`calculate_sasa`は単一構造用。複数フレーム対応を追加。

```zig
// 新規追加
pub fn calculateSasaBatch(
    coordinates: [][*]const f32,  // n_frames x (n_atoms * 3)
    n_frames: usize,
    n_atoms: usize,
    radii: []const f32,
    config: Config,
    results: [][]f32,             // n_frames x n_atoms (output)
) !void
```

最適化ポイント:
- [ ] フレーム間並列（既存スレッドプール活用）
- [ ] バッファ再利用（隣接リスト、点群など）
- [ ] 原子半径は1回だけ設定

### Phase 2: C API追加

```c
// freesasa_zig.h に追加
int freesasa_calculate_batch(
    const float* coordinates,    // n_frames * n_atoms * 3
    int n_frames,
    int n_atoms,
    const float* radii,
    FreeSasaConfig config,
    float* results               // n_frames * n_atoms (output)
);
```

### Phase 3: Python binding更新

```python
# python/freesasa_zig/_core.py

def compute_sasa_batch(
    coordinates: np.ndarray,
    radii: np.ndarray,
    algorithm: str = 'sr',
    n_threads: int | None = None,
    probe_radius: float = 1.4,
    n_points: int = 100,
) -> np.ndarray:
    """
    Compute SASA for multiple frames (batch processing).

    Parameters
    ----------
    coordinates : np.ndarray, shape (n_frames, n_atoms, 3)
        Atomic coordinates in Angstrom
    radii : np.ndarray, shape (n_atoms,)
        Van der Waals radii in Angstrom
    ...

    Returns
    -------
    sasa : np.ndarray, shape (n_frames, n_atoms)
        SASA values in Angstrom^2
    """
    ...
```

### Phase 4: MDTraj wrapper

```python
# python/freesasa_zig/mdtraj.py

def compute_sasa_trajectory(
    traj,  # mdtraj.Trajectory
    algorithm: str = 'sr',
    n_threads: int | None = None,
    probe_radius: float = 1.4,
    n_points: int = 100,
    mode: str = 'atom',
) -> np.ndarray:
    """
    Compute SASA for MDTraj trajectory.

    Drop-in replacement for mdtraj.shrake_rupley() with better performance.
    """
    # Extract coordinates (nm -> Angstrom)
    coords = traj.xyz * 10.0  # MDTraj uses nm

    # Get radii from topology
    radii = np.array([
        _ATOMIC_RADII[atom.element.symbol]
        for atom in traj.topology.atoms
    ])

    # Call batch API
    sasa = compute_sasa_batch(
        coords, radii, algorithm, n_threads, probe_radius, n_points
    )

    # Aggregate by residue if needed
    if mode == 'residue':
        sasa = _aggregate_by_residue(sasa, traj.topology)

    return sasa
```

## Optimization Strategy

### Memory Reuse (Key for MD)

```
First frame:
  - Allocate buffers (neighbor list, sphere points, etc.)
  - Build radii array from topology

Subsequent frames:
  - Reuse all buffers
  - Only update coordinates
```

### Parallelization

```
Option 1: Frame-level parallelism
  Thread 0: Frame 0, 4, 8, ...
  Thread 1: Frame 1, 5, 9, ...
  Thread 2: Frame 2, 6, 10, ...
  Thread 3: Frame 3, 7, 11, ...

Option 2: Hybrid (for large systems)
  Frame-level + atom-level parallelism within frame
```

## Benchmark Plan

### Test Cases

1. **Small protein**: Lysozyme (~1,000 atoms) x 10,000 frames
2. **Medium protein**: ~10,000 atoms x 1,000 frames
3. **Large system**: ~100,000 atoms x 100 frames

### Comparison

- MDTraj `shrake_rupley()`
- FreeSASA (if trajectory support added)
- freesasa-zig batch

### Metrics

- Total time for trajectory
- Per-frame time
- Memory usage
- Thread scaling

## File Structure

```
python/
├── freesasa_zig/
│   ├── __init__.py
│   ├── _core.py          # Low-level binding (existing)
│   ├── batch.py          # NEW: Batch processing API
│   ├── mdtraj.py         # NEW: MDTraj integration
│   └── mdanalysis.py     # FUTURE: MDAnalysis integration
└── tests/
    ├── test_sasa.py      # Existing tests
    └── test_trajectory.py # NEW: Trajectory tests
```

## Dependencies

### Required
- numpy (already)

### Optional (for high-level API)
- mdtraj (for compute_sasa_trajectory)
- mdanalysis (future)

```toml
# pyproject.toml
[project.optional-dependencies]
mdtraj = ["mdtraj>=1.9"]
mdanalysis = ["MDAnalysis>=2.0"]
all = ["mdtraj>=1.9", "MDAnalysis>=2.0"]
```

## Questions to Resolve

1. **Unit conversion**: MDTraj uses nm, freesasa-zig uses Angstrom. Convert in Python or Zig?
   - Recommendation: Python側で変換（明示的で分かりやすい）

2. **Streaming mode**: 巨大トラジェクトリ用にチャンク処理が必要か？
   - MDTrajの`iterload()`と組み合わせて使えばOK

3. **Residue aggregation**: Zig側でやるかPython側でやるか？
   - Recommendation: Python側（柔軟性重視）

## Success Criteria

- [ ] MDTrajより2倍以上高速
- [ ] 同じ結果（numerical tolerance内）
- [ ] メモリ使用量が線形（フレーム数に対して）
- [ ] 使いやすいAPI

---
- [ ] **DONE** - Phase complete

# Shrake-Rupleyアルゴリズム詳解

## 概要

Shrake-Rupley法は、タンパク質などの生体分子の溶媒接触可能表面積（SASA: Solvent Accessible Surface Area）を計算するための数値的手法。1973年にShrakeとRupleyによって提案された。

## 理論的背景

### 溶媒接触可能表面積とは

- 分子表面を「プローブ球」（通常は水分子を表す半径1.4Åの球）でなぞったときの軌跡
- 各原子の「拡張半径」= van der Waals半径 + プローブ半径
- 生化学的に重要：タンパク質の折り畳み自由エネルギー、結合エネルギーの推定に使用

### アルゴリズムの基本原理

1. 各原子の表面に均一にテストポイントを配置
2. 各テストポイントが他の原子の内部にあるかチェック
3. 露出しているポイントの割合からSASAを計算

```
SASA_i = 4π × r_i² × (露出ポイント数 / 総ポイント数)
```

## 実装詳細

### ファイル: `src/shrake_rupley.zig`

#### メインAPI

```zig
pub fn calculateSasa(
    allocator: Allocator,
    atoms: types.AtomInput,
    options: SasaOptions,
) !types.SasaResult
```

**パラメータ:**
- `allocator`: メモリアロケータ
- `atoms`: 原子座標と半径（SoA形式）
- `options`: 計算オプション

**オプション構造体:**
```zig
pub const SasaOptions = struct {
    n_points: u32 = 100,        // テストポイント数（デフォルト: 100）
    probe_radius: f64 = 1.4,    // プローブ半径（デフォルト: 1.4 Å）
    n_threads: usize = 0,       // スレッド数（0 = 自動検出）
};
```

#### 計算フロー

```zig
pub fn calculateSasa(...) !types.SasaResult {
    // 1. テストポイント生成（単位球上）
    const test_points = try test_points_mod.generateTestPoints(
        allocator,
        options.n_points
    );
    defer allocator.free(test_points);

    // 2. 近傍リスト構築
    const max_radius = findMaxRadius(atoms.r);
    const cutoff = 2.0 * (max_radius + options.probe_radius);
    var neighbor_list = try NeighborList.init(allocator, atoms, cutoff);
    defer neighbor_list.deinit();

    // 3. SASA計算（並列または逐次）
    if (options.n_threads > 1) {
        return calculateSasaParallel(...);
    } else {
        return calculateSasaSequential(...);
    }
}
```

### 単一原子のSASA計算

```zig
fn calculateAtomSasa(
    atom_idx: usize,
    atoms: types.AtomInput,
    test_points: []const types.Vec3,
    neighbor_list: *const NeighborList,
    probe_radius: f64,
) f64 {
    const center = types.Vec3{
        .x = atoms.x[atom_idx],
        .y = atoms.y[atom_idx],
        .z = atoms.z[atom_idx],
    };
    const radius = atoms.r[atom_idx] + probe_radius;

    // 近傍原子を取得（O(1)）
    const neighbors = neighbor_list.getNeighbors(atom_idx);

    var exposed_count: u32 = 0;

    // 各テストポイントをチェック
    for (test_points) |tp| {
        // テストポイントを原子表面に配置
        const point = types.Vec3{
            .x = center.x + tp.x * radius,
            .y = center.y + tp.y * radius,
            .z = center.z + tp.z * radius,
        };

        // 近傍原子との衝突チェック（SIMD最適化）
        if (!isPointBuried(point, neighbors, atoms, probe_radius)) {
            exposed_count += 1;
        }
    }

    // SASA = 4πr² × (露出率)
    const total_points: f64 = @floatFromInt(test_points.len);
    const exposed_ratio = @as(f64, @floatFromInt(exposed_count)) / total_points;
    return 4.0 * std.math.pi * radius * radius * exposed_ratio;
}
```

### 埋没判定（SIMD最適化版）

```zig
fn isPointBuried(
    point: types.Vec3,
    neighbors: []const usize,
    atoms: types.AtomInput,
    probe_radius: f64,
) bool {
    // 4原子ずつSIMD処理
    var i: usize = 0;
    while (i + 4 <= neighbors.len) : (i += 4) {
        const indices = neighbors[i..][0..4];

        // SIMD距離計算
        if (simd.isAnyWithinRadius(
            point,
            atoms,
            indices,
            probe_radius
        )) {
            return true;  // 埋没している
        }
    }

    // 残りを逐次処理
    while (i < neighbors.len) : (i += 1) {
        const j = neighbors[i];
        const neighbor_center = types.Vec3{
            .x = atoms.x[j],
            .y = atoms.y[j],
            .z = atoms.z[j],
        };
        const neighbor_radius = atoms.r[j] + probe_radius;

        const diff = point.sub(neighbor_center);
        if (diff.lengthSquared() < neighbor_radius * neighbor_radius) {
            return true;  // 埋没している
        }
    }

    return false;  // 露出している
}
```

## テストポイント生成

### ファイル: `src/test_points.zig`

### Golden Section Spiralアルゴリズム

球面上に均一にポイントを分布させる効率的な方法。

```zig
pub fn generateTestPoints(allocator: Allocator, n: u32) ![]types.Vec3 {
    const points = try allocator.alloc(types.Vec3, n);
    errdefer allocator.free(points);

    const golden_angle = std.math.pi * (3.0 - @sqrt(5.0));  // ≈ 2.399963...

    for (0..n) |i| {
        const i_f: f64 = @floatFromInt(i);
        const n_f: f64 = @floatFromInt(n);

        // y座標: -1から1へ均等に分布
        const y = 1.0 - (2.0 * i_f + 1.0) / n_f;

        // xz平面での半径
        const r = @sqrt(1.0 - y * y);

        // 黄金角で回転
        const theta = golden_angle * i_f;

        points[i] = types.Vec3{
            .x = r * @cos(theta),
            .y = y,
            .z = r * @sin(theta),
        };
    }

    return points;
}
```

**黄金角の数学的背景:**
- 黄金角 = π(3 - √5) ≈ 137.5°
- 連続するポイントがこの角度で回転することで、螺旋パターンが形成
- ひまわりの種の配置と同じ原理
- どの視点から見ても均一に分布

**ポイント数と精度のトレードオフ:**

| ポイント数 | 精度 | 計算時間（相対） |
|-----------|------|-----------------|
| 50 | 低 | 0.5x |
| 100 | 中（デフォルト） | 1.0x |
| 200 | 高 | 2.0x |
| 1000 | 非常に高 | 10.0x |

## 精度検証

### FreeSASAとの比較

FreeSASA（C実装の参照実装）との比較結果:

| 構造 | 原子数 | FreeSASA (Å²) | Zig実装 (Å²) | 差分 |
|------|--------|---------------|--------------|------|
| 1A0Q | 3,183 | 18,923.28 | 19,211.19 | 1.52% |

**差分の原因:**
1. テストポイント生成アルゴリズムの違い（Golden Spiral vs Fibonacci Lattice）
2. 浮動小数点演算の順序の違い
3. 近傍リストのカットオフ距離の違い

1.52%の差分は、SASA計算の用途（相対比較、トレンド分析）では許容範囲内。

## 計算量解析

| 処理 | 計算量 | 備考 |
|------|--------|------|
| テストポイント生成 | O(P) | P = ポイント数、事前計算 |
| 近傍リスト構築 | O(N) | N = 原子数 |
| 近傍探索（1原子） | O(1) | 空間ハッシュにより定数時間 |
| 埋没判定（1ポイント） | O(K) | K = 平均近傍数（通常 < 20） |
| 全体 | O(N × P × K) | 実質 O(N × P) |

**従来のナイーブ実装との比較:**
- ナイーブ: O(N² × P) - 全原子ペアをチェック
- 最適化版: O(N × P × K) - 近傍のみチェック
- 大規模分子（N > 1000）では10倍以上の高速化

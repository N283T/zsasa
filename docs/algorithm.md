# SASAアルゴリズム詳解

## 概要

本実装では2つのSASA計算アルゴリズムを提供する：

1. **Shrake-Rupley (SR)**: テストポイント法（1973年）
2. **Lee-Richards (LR)**: スライス法（1971年）

## アルゴリズム比較

| 特性 | Shrake-Rupley | Lee-Richards |
|------|---------------|--------------|
| 発表年 | 1973 | 1971 |
| 手法 | テストポイント | スライス/円弧積分 |
| 精度制御 | ポイント数 (`--n-points`) | スライス数 (`--n-slices`) |
| 計算量 | O(N × P × K) | O(N × S × K) |
| 長所 | 実装が単純、高速 | 数学的に厳密 |
| 短所 | 統計的誤差 | 計算コスト高 |
| 推奨用途 | 大規模構造、迅速な解析 | 高精度が必要な場合 |

**凡例:** N=原子数, P=ポイント数, S=スライス数, K=平均近傍数

### パフォーマンス比較（1A0Q, 3183原子, 4スレッド, ReleaseFast）

| アルゴリズム | Zig | FreeSASA C | 高速化 |
|-------------|-----|------------|--------|
| Shrake-Rupley (100点) | 2.6ms | 4.4ms | 1.7x |
| Lee-Richards (20スライス) | 9.9ms | 16.1ms | 1.6x |

詳細な最適化手法については [docs/optimizations.md](optimizations.md) を参照。

### どちらを選ぶべきか

- **Shrake-Rupley**: デフォルト。ほとんどの用途に適切
- **Lee-Richards**: 以下の場合に推奨
  - 絶対値の精度が重要な場合
  - 他のLee-Richards実装との比較が必要な場合
  - 円弧の幾何学的な厳密性が必要な場合

---

# Shrake-Rupleyアルゴリズム

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
    var i: usize = 0;

    // 8原子ずつSIMD処理（8-wide SIMD）
    while (i + 8 <= neighbors.len) : (i += 8) {
        const indices = neighbors[i..][0..8];
        if (simd.isAnyWithinRadiusBatch8(point, atoms, indices, probe_radius)) {
            return true;
        }
    }

    // 4原子ずつSIMD処理（残り用）
    while (i + 4 <= neighbors.len) : (i += 4) {
        const indices = neighbors[i..][0..4];
        if (simd.isAnyWithinRadius(point, atoms, indices, probe_radius)) {
            return true;
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
            return true;
        }
    }

    return false;
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

---

# Lee-Richardsアルゴリズム

## 概要

Lee-Richards法は、1971年にLeeとRichardsによって提案されたSASA計算手法。原子球をZ軸方向にスライスし、各スライスでの露出円弧を積分することで表面積を計算する。

## 理論的背景

### スライスベースのアプローチ

```
        ╭─────╮
       ╱   ●   ╲      ← 原子球
      │    │    │
   ───┼────┼────┼───  ← スライス（Z軸に垂直）
      │    │    │
       ╲       ╱
        ╰─────╯

スライス断面:
    ╭───╮
   ╱     ╲
  │   ●   │  ← 原子断面（円）
   ╲     ╱
    ╰───╯
     ↑
    露出円弧を計算
```

各スライスで：
1. 原子の断面円を計算
2. 近傍原子による遮蔽円弧を特定
3. 露出円弧の長さを計算
4. スライス厚と掛け合わせて表面積に寄与

### 数学的定式化

スライスzでの原子iの断面半径：
```
R'_i(z) = √(R_i² - (z - z_i)²)
```

ここで R_i = 原子半径 + プローブ半径

総表面積：
```
SASA_i = ∫ R_i × L_exposed(z) dz
```

ここで L_exposed(z) はスライスzでの露出円弧長。

## 実装詳細

### ファイル: `src/lee_richards.zig`

#### メインAPI

```zig
pub fn calculateSasa(
    allocator: Allocator,
    input: AtomInput,
    config: LeeRichardsConfig,
) !SasaResult
```

**設定構造体:**
```zig
pub const LeeRichardsConfig = struct {
    n_slices: u32 = 20,      // スライス数（デフォルト: 20）
    probe_radius: f64 = 1.4, // プローブ半径（デフォルト: 1.4 Å）
};
```

### 円弧計算

2つの円の交差から、遮蔽される円弧を計算：

```
2つの円の交点から、原子1の遮蔽円弧を計算

    ╭───╮
   ╱  1  ╲╲
  │   ●───●│ 2  ← 円2が円1を遮蔽
   ╲     ╱╱
    ╰───╯

遮蔽角度 α = acos((R'_i² + d_ij² - R'_j²) / (2 × R'_i × d_ij))
遮蔽円弧の中心角 β = atan2(dy, dx) + π
```

### SIMD最適化

Lee-Richardsでも8近傍ずつのSIMD処理を行う：

```zig
// 8-wide SIMD: スライス半径とxy距離を計算
const rj_primes = simd.sliceRadiiBatch8(slice_z, batch_z, batch_r);
const dij_batch = simd.xyDistanceBatch8(xi, yi, batch_x, batch_y);

// 8-wide SIMD: 円の重なりをチェック
const overlap_mask = simd.circlesOverlapBatch8(dij_batch, Ri_prime, rj_primes);
```

### 高速三角関数

円弧角度計算では、標準ライブラリの`acos`/`atan2`の代わりに多項式近似を使用：

```zig
// 高速近似（最大誤差: acos ~0.02°, atan2 ~0.09°）
const alpha = simd.fastAcos(cos_alpha);
const beta = simd.fastAtan2(dy, dx) + std.math.pi;
```

これにより約27%の高速化を達成。精度は0.3%以内で実用上問題なし。

### 円弧のマージ

複数の近傍原子からの遮蔽円弧が重なる場合、マージ処理が必要：

```zig
fn exposedArcLength(arcs: []Arc) f64 {
    // 円弧を開始角度でソート
    sortArcs(arcs);

    // 重なりをマージして露出長を計算
    var sum: f64 = arcs[0].start;
    var sup: f64 = arcs[0].end;

    for (arcs[1..]) |arc| {
        if (sup < arc.start) {
            sum += arc.start - sup;  // ギャップ（露出部分）
        }
        sup = @max(sup, arc.end);
    }

    sum += TWOPI - sup;  // 最後の円弧以降の露出
    return sum;
}
```

## スライス数と精度のトレードオフ

| スライス数 | 精度 | 計算時間（相対） |
|-----------|------|-----------------|
| 10 | 低 | 0.5x |
| 20 | 中（デフォルト） | 1.0x |
| 50 | 高 | 2.5x |
| 100 | 非常に高 | 5.0x |

## Shrake-Rupleyとの精度比較

同じ構造（1A0Q）での比較：

| アルゴリズム | パラメータ | 結果 (Å²) |
|-------------|-----------|-----------|
| Shrake-Rupley | 100ポイント | 19,211.19 |
| Lee-Richards | 20スライス | 19,201.26 |
| 差分 | | 0.05% |

両アルゴリズムは非常に近い結果を出力する。

## 計算量解析

| 処理 | 計算量 | 備考 |
|------|--------|------|
| 近傍リスト構築 | O(N) | N = 原子数 |
| スライス処理（1原子） | O(S × K) | S = スライス数, K = 近傍数 |
| 円弧ソート | O(K log K) | 挿入ソート（小配列に効率的） |
| 全体 | O(N × S × K) | 実質 O(N × S) |

## 参考文献

- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.

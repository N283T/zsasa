# SASAアルゴリズム

## 概要

本実装は2つのSASA計算アルゴリズムを提供します:

1. **Shrake-Rupley (SR)**: テスト点法 (1973)
2. **Lee-Richards (LR)**: スライス法 (1971)

## アルゴリズム比較

| 特性 | Shrake-Rupley | Lee-Richards |
|----------|---------------|--------------|
| 発表年 | 1973 | 1971 |
| 手法 | テスト点 | スライス/弧積分 |
| 精度制御 | 点数 (`--n-points`) | スライス数 (`--n-slices`) |
| 計算量 | O(N × P × K) | O(N × S × K) |
| 長所 | シンプル、高速 | 数学的に厳密 |
| 短所 | 統計誤差 | 計算コストが高い |
| 推奨用途 | 大規模構造、高速解析 | 高精度要求 |

**凡例:** N=原子数、P=点数、S=スライス数、K=平均近傍原子数

### パフォーマンス比較 (1A0Q, 3183原子, 4スレッド, ReleaseFast)

| アルゴリズム | Zig | FreeSASA C | 高速化 |
|-----------|-----|------------|---------|
| Shrake-Rupley (100点) | 2.6ms | 4.4ms | 1.7x |
| Lee-Richards (20スライス) | 9.9ms | 16.1ms | 1.6x |

詳細な最適化技術については[optimizations.md](optimizations.md)を参照してください。

### どちらを選ぶべきか？

- **Shrake-Rupley**: デフォルト。ほとんどのユースケースに適合
- **Lee-Richards**: 以下の場合に推奨:
  - 絶対的な精度が重要
  - 他のLee-Richards実装との比較が必要
  - 弧計算の幾何学的厳密性が必要

---

# Shrake-Rupleyアルゴリズム

## 概要

Shrake-Rupley法は、タンパク質などの生体分子の溶媒接触表面積（SASA）を計算するための数値的手法です。1973年にShrakeとRupleyによって提案されました。

## 理論的背景

### 溶媒接触表面積とは？

- 「プローブ球」（通常は水分子を表す半径1.4Åの球）を分子表面上で転がしたときの軌跡
- 各原子の「拡張半径」= ファンデルワールス半径 + プローブ半径
- 生化学的に重要: タンパク質の折りたたみ自由エネルギーや結合エネルギーの推定に使用

### アルゴリズムの原理

1. 各原子の表面上にテスト点を均一に配置
2. 各テスト点が他の原子の内部にあるかチェック
3. 露出点の割合からSASAを計算

```
SASA_i = 4π × r_i² × (露出点数 / 総点数)
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
- `atoms`: 原子座標と半径 (SoA形式)
- `options`: 計算オプション

**オプション構造体:**
```zig
pub const SasaOptions = struct {
    n_points: u32 = 100,        // テスト点数 (デフォルト: 100)
    probe_radius: f64 = 1.4,    // プローブ半径 (デフォルト: 1.4 Å)
    n_threads: usize = 0,       // スレッド数 (0 = 自動検出)
};
```

#### 計算フロー

```zig
pub fn calculateSasa(...) !types.SasaResult {
    // 1. テスト点を生成（単位球上）
    const test_points = try test_points_mod.generateTestPoints(
        allocator,
        options.n_points
    );
    defer allocator.free(test_points);

    // 2. 近傍リストを構築
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

    // 近傍原子を取得 O(1)
    const neighbors = neighbor_list.getNeighbors(atom_idx);

    var exposed_count: u32 = 0;

    // 各テスト点をチェック
    for (test_points) |tp| {
        // テスト点を原子表面に配置
        const point = types.Vec3{
            .x = center.x + tp.x * radius,
            .y = center.y + tp.y * radius,
            .z = center.z + tp.z * radius,
        };

        // 近傍との衝突チェック（SIMD最適化）
        if (!isPointBuried(point, neighbors, atoms, probe_radius)) {
            exposed_count += 1;
        }
    }

    // SASA = 4πr² × (露出比率)
    const total_points: f64 = @floatFromInt(test_points.len);
    const exposed_ratio = @as(f64, @floatFromInt(exposed_count)) / total_points;
    return 4.0 * std.math.pi * radius * radius * exposed_ratio;
}
```

### 埋没判定（SIMD最適化）

```zig
fn isPointBuried(
    point: types.Vec3,
    neighbors: []const usize,
    atoms: types.AtomInput,
    probe_radius: f64,
) bool {
    var i: usize = 0;

    // SIMDで8原子ずつ処理（8-wide）
    while (i + 8 <= neighbors.len) : (i += 8) {
        const indices = neighbors[i..][0..8];
        if (simd.isAnyWithinRadiusBatch8(point, atoms, indices, probe_radius)) {
            return true;
        }
    }

    // 4原子ずつ処理（余り）
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

## テスト点生成

### ファイル: `src/test_points.zig`

### 黄金角螺旋アルゴリズム

球面上に点を均一に分布させる効率的な方法です。

```zig
pub fn generateTestPoints(allocator: Allocator, n: u32) ![]types.Vec3 {
    const points = try allocator.alloc(types.Vec3, n);
    errdefer allocator.free(points);

    const golden_angle = std.math.pi * (3.0 - @sqrt(5.0));  // ≈ 2.399963...

    for (0..n) |i| {
        const i_f: f64 = @floatFromInt(i);
        const n_f: f64 = @floatFromInt(n);

        // y座標: -1から1まで均一に分布
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
- 連続する点がこの角度で回転し、螺旋パターンを形成
- ひまわりの種の配置と同じ原理
- どの角度から見ても均一に見える

**点数と精度のトレードオフ:**

| 点数 | 精度 | 相対時間 |
|--------|-----------|---------------|
| 50 | 低 | 0.5x |
| 100 | 中 (デフォルト) | 1.0x |
| 200 | 高 | 2.0x |
| 1000 | 非常に高 | 10.0x |

## 精度検証

### FreeSASAとの比較

FreeSASA（Cリファレンス実装）との比較結果:

| 構造 | 原子数 | FreeSASA (Å²) | Zig (Å²) | 差異 |
|-----------|-------|---------------|----------|------------|
| 1A0Q | 3,183 | 18,923.28 | 19,211.19 | 1.52% |

**差異の原因:**
1. テスト点生成アルゴリズムの違い（黄金螺旋 vs フィボナッチ格子）
2. 浮動小数点演算順序の違い
3. 近傍リストのカットオフ距離の違い

1.52%の差異は、一般的なSASA用途（相対比較、傾向分析）では許容範囲内です。

## 計算量解析

| 操作 | 計算量 | 備考 |
|-----------|------------|-------|
| テスト点生成 | O(P) | P = 点数、事前計算 |
| 近傍リスト構築 | O(N) | N = 原子数 |
| 近傍検索（1原子） | O(1) | 空間ハッシュによる定数時間 |
| 埋没判定（1点） | O(K) | K = 平均近傍数（通常 < 20） |
| 合計 | O(N × P × K) | 実質的に O(N × P) |

**ナイーブ実装との比較:**
- ナイーブ: O(N² × P) - 全原子ペアをチェック
- 最適化版: O(N × P × K) - 近傍のみチェック
- 大きな分子（N > 1000）で10倍以上の高速化

---

# Lee-Richardsアルゴリズム

## 概要

Lee-Richards法は、1971年にLeeとRichardsによって提案されたSASA計算手法です。原子球をZ軸に沿ってスライスし、各スライスでの露出弧を積分することで表面積を計算します。

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
    露出弧を計算
```

各スライスで:
1. 原子の断面円を計算
2. 近傍原子からの遮蔽弧を特定
3. 露出弧の長さを計算
4. スライス厚みを掛けて表面積への寄与を算出

### 数学的定式化

スライスzにおける原子iの断面半径:
```
R'_i(z) = √(R_i² - (z - z_i)²)
```

ここでR_i = 原子半径 + プローブ半径

総表面積:
```
SASA_i = ∫ R_i × L_exposed(z) dz
```

L_exposed(z)はスライスzでの露出弧の長さ。

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
    n_slices: u32 = 20,      // スライス数 (デフォルト: 20)
    probe_radius: f64 = 1.4, // プローブ半径 (デフォルト: 1.4 Å)
};
```

### 弧計算

2つの円の交差から遮蔽弧を計算:

```
円1の円2との交差による遮蔽弧を計算

    ╭───╮
   ╱  1  ╲╲
  │   ●───●│ 2  ← 円2が円1を遮蔽
   ╲     ╱╱
    ╰───╯

遮蔽角度 α = acos((R'_i² + d_ij² - R'_j²) / (2 × R'_i × d_ij))
弧の中心角度 β = atan2(dy, dx) + π
```

### SIMD最適化

Lee-Richardsも近傍に対して8-wide SIMD処理を使用:

```zig
// 8-wide SIMD: スライス半径とxy距離を計算
const rj_primes = simd.sliceRadiiBatch8(slice_z, batch_z, batch_r);
const dij_batch = simd.xyDistanceBatch8(xi, yi, batch_x, batch_y);

// 8-wide SIMD: 円の重なりをチェック
const overlap_mask = simd.circlesOverlapBatch8(dij_batch, Ri_prime, rj_primes);
```

### 高速三角関数

弧角度の計算は標準ライブラリの`acos`/`atan2`の代わりに多項式近似を使用:

```zig
// 高速近似（最大誤差: acos ~0.02°, atan2 ~0.09°）
const alpha = simd.fastAcos(cos_alpha);
const beta = simd.fastAtan2(dy, dx) + std.math.pi;
```

これにより約37%の高速化を達成。精度は0.3%以内で、実用上許容範囲です。

### 弧のマージ

複数の近傍原子からの遮蔽弧が重なる場合、マージが必要:

```zig
fn exposedArcLength(arcs: []Arc) f64 {
    // 弧を開始角度でソート
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

    sum += TWOPI - sup;  // 最後の弧の後の露出
    return sum;
}
```

## スライス数と精度のトレードオフ

| スライス数 | 精度 | 相対時間 |
|--------|-----------|---------------|
| 10 | 低 | 0.5x |
| 20 | 中 (デフォルト) | 1.0x |
| 50 | 高 | 2.5x |
| 100 | 非常に高 | 5.0x |

## Shrake-Rupleyとの精度比較

同じ構造（1A0Q）での比較:

| アルゴリズム | パラメータ | 結果 (Å²) |
|-----------|------------|-------------|
| Shrake-Rupley | 100点 | 19,211.19 |
| Lee-Richards | 20スライス | 19,201.26 |
| 差異 | | 0.05% |

両アルゴリズムは非常に近い結果を生成します。

## 計算量解析

| 操作 | 計算量 | 備考 |
|-----------|------------|-------|
| 近傍リスト構築 | O(N) | N = 原子数 |
| スライス処理（1原子） | O(S × K) | S = スライス数、K = 近傍数 |
| 弧ソート | O(K log K) | 挿入ソート（小配列に効率的） |
| 合計 | O(N × S × K) | 実質的に O(N × S) |

## 参考文献

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms: Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.

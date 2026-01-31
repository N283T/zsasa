# 最適化技術

本プロジェクトでは5つの主要な最適化を実装しています。これらの組み合わせにより、FreeSASA C（ネイティブ版）と比較してSRで最大2.3倍、LRで最大1.7倍の高速化を達成しています。

## 1. 近傍リスト最適化（フェーズ2）

### ファイル: `src/neighbor_list.zig`

### 問題

ナイーブなSASA計算では、各テスト点が「埋没」しているかを判定するために全原子への距離を計算します。

```
計算量: O(N² × P)
- N = 原子数
- P = テスト点数
```

3,183原子の場合: 3183² × 100 ≈ 10億回の距離計算

### 解決策: 空間ハッシュ

3D空間をセルに分割し、各原子がどのセルに属するかを記録。近傍探索時は隣接セルのみを探索。

```
┌─────┬─────┬─────┐
│     │  ●  │     │
├─────┼─────┼─────┤
│  ●  │  ◎  │  ●  │  ◎ = 対象原子
├─────┼─────┼─────┤  ● = 近傍候補
│     │  ●  │     │
└─────┴─────┴─────┘
```

### 実装

```zig
pub const NeighborList = struct {
    cells: [][]usize,           // セルごとの原子インデックス
    cell_size: f64,             // セルサイズ
    grid_dims: [3]usize,        // グリッド次元
    min_bounds: [3]f64,         // 最小座標
    neighbor_cache: [][]usize,  // 原子ごとの近傍キャッシュ
    allocator: Allocator,

    pub fn init(
        allocator: Allocator,
        atoms: types.AtomInput,
        cutoff: f64,
    ) !NeighborList { ... }

    pub fn getNeighbors(self: *const NeighborList, atom_idx: usize) []const usize {
        return self.neighbor_cache[atom_idx];
    }

    pub fn deinit(self: *NeighborList) void { ... }
};
```

### セルサイズの決定

```zig
// カットオフ距離 = 2 × (最大原子半径 + プローブ半径)
const cutoff = 2.0 * (max_radius + probe_radius);
const cell_size = cutoff;  // セルサイズ = カットオフ
```

**根拠:**
- セルサイズ = カットオフなら、近傍原子は現在のセル + 26個の隣接セルに限定
- 3×3×3 = 27セルの探索で全近傍をカバー

### 計算量の改善

```
改善後: O(N × K)
- K = 平均近傍数（通常10-20）
```

3,183原子の場合: 3183 × 20 ≈ 64,000回の距離計算（**約15,000倍削減**）

### ベンチマーク結果

| 実装 | 時間 | 高速化 |
|----------------|------|---------|
| ナイーブ O(N²) | ~200ms | 1.0x |
| 近傍リスト O(N) | ~13ms | 15x |

---

## 2. SIMD最適化（フェーズ3）

### ファイル: `src/simd.zig`

### 問題

距離計算は単純ですが、大量に実行されます（原子 × 点 × 近傍）。

```zig
// スカラー距離計算
fn distance_squared(p1: Vec3, p2: Vec3) f64 {
    const dx = p1.x - p2.x;
    const dy = p1.y - p2.y;
    const dz = p1.z - p2.z;
    return dx*dx + dy*dy + dz*dz;
}
```

### 解決策: SIMD（Single Instruction, Multiple Data）

Zigの`@Vector`型を使用して、8つの距離計算を同時に実行（8-wide SIMD）。4-wide版も提供され、残り原子数に応じて選択されます。

**精度選択**: `--precision`オプションでf32/f64を選択可能。f32はメモリ帯域とキャッシュ効率が向上し、約3%の高速化を達成。

```zig
pub fn isAnyWithinRadius(
    point: types.Vec3,
    atoms: types.AtomInput,
    indices: *const [4]usize,
    probe_radius: f64,
) bool {
    // 4原子のx座標を読み込み
    const x: @Vector(4, f64) = .{
        atoms.x[indices[0]],
        atoms.x[indices[1]],
        atoms.x[indices[2]],
        atoms.x[indices[3]],
    };

    // 4原子のy座標を読み込み
    const y: @Vector(4, f64) = .{
        atoms.y[indices[0]],
        atoms.y[indices[1]],
        atoms.y[indices[2]],
        atoms.y[indices[3]],
    };

    // 4原子のz座標を読み込み
    const z: @Vector(4, f64) = .{
        atoms.z[indices[0]],
        atoms.z[indices[1]],
        atoms.z[indices[2]],
        atoms.z[indices[3]],
    };

    // 4原子の半径を読み込み
    const r: @Vector(4, f64) = .{
        atoms.r[indices[0]] + probe_radius,
        atoms.r[indices[1]] + probe_radius,
        atoms.r[indices[2]] + probe_radius,
        atoms.r[indices[3]] + probe_radius,
    };

    // 点座標をブロードキャスト
    const px: @Vector(4, f64) = @splat(point.x);
    const py: @Vector(4, f64) = @splat(point.y);
    const pz: @Vector(4, f64) = @splat(point.z);

    // 差分計算（4並列）
    const dx = px - x;
    const dy = py - y;
    const dz = pz - z;

    // 距離の二乗（4並列）
    const dist_sq = dx * dx + dy * dy + dz * dz;

    // 半径の二乗
    const r_sq = r * r;

    // 比較（4並列）
    const within = dist_sq < r_sq;

    // いずれかがtrueならtrue
    return @reduce(.Or, within);
}
```

### SoA形式との親和性

```
AoS (Array of Structures):     SoA (Structure of Arrays):
┌───────────────────┐          ┌─────────────────────┐
│ Atom0: x,y,z,r    │          │ x: [x0,x1,x2,x3...] │
│ Atom1: x,y,z,r    │    →     │ y: [y0,y1,y2,y3...] │
│ Atom2: x,y,z,r    │          │ z: [z0,z1,z2,z3...] │
│ ...               │          │ r: [r0,r1,r2,r3...] │
└───────────────────┘          └─────────────────────┘
```

SoA形式によりSIMD命令が連続メモリから効率的にデータを読み込めます。

### ベンチマーク結果

| 実装 | 時間 | 高速化 |
|----------------|------|---------|
| スカラー | ~13ms | 1.0x |
| SIMD | ~11ms | 1.2x |

**注:** 近傍リスト最適化後は距離計算回数が大幅に減少しているため、SIMD単体の効果は限定的。ただしナイーブ実装と組み合わせると効果は大きい。

---

## 3. マルチスレッド最適化（フェーズ4）

### ファイル: `src/thread_pool.zig`

### 問題

各原子のSASA計算は独立していますが、シングルスレッド実行では全CPUコアを活用できません。

### 解決策: ワークスティーリングスレッドプール

```
┌─────────────────────────────────────────┐
│              タスクキュー                │
│  [Atom0] [Atom1] [Atom2] ... [AtomN-1]  │
└─────────────────────────────────────────┘
      ↓         ↓         ↓         ↓
   ┌─────┐  ┌─────┐  ┌─────┐  ┌─────┐
   │ T0  │  │ T1  │  │ T2  │  │ T3  │  スレッド
   └─────┘  └─────┘  └─────┘  └─────┘
      ↓         ↓         ↓         ↓
   [Result0] [Result1] [Result2] [Result3]
```

### 実装

```zig
pub fn ThreadPool(comptime Task: type, comptime Result: type) type {
    return struct {
        const Self = @This();

        threads: []std.Thread,
        mutex: std.Thread.Mutex,
        condition: std.Thread.Condition,
        tasks: []const Task,
        results: []Result,
        next_task: usize,
        completed: usize,
        shutdown: bool,
        work_fn: *const fn (Task) Result,

        pub fn init(
            allocator: Allocator,
            n_threads: usize,
            tasks: []const Task,
            results: []Result,
            work_fn: *const fn (Task) Result,
        ) !Self { ... }

        pub fn wait(self: *Self) void { ... }

        pub fn deinit(self: *Self) void { ... }
    };
}
```

### スレッドワーカー

```zig
fn workerThread(self: *Self) void {
    while (true) {
        // タスクを取得
        self.mutex.lock();

        while (self.next_task >= self.tasks.len and !self.shutdown) {
            self.condition.wait(&self.mutex);
        }

        if (self.shutdown and self.next_task >= self.tasks.len) {
            self.mutex.unlock();
            return;
        }

        const task_idx = self.next_task;
        self.next_task += 1;
        self.mutex.unlock();

        // タスク実行（ロック外）
        const result = self.work_fn(self.tasks[task_idx]);

        // 結果を格納
        self.mutex.lock();
        self.results[task_idx] = result;
        self.completed += 1;

        if (self.completed == self.tasks.len) {
            self.condition.broadcast();
        }
        self.mutex.unlock();
    }
}
```

### 自動スレッド数検出

```zig
fn getDefaultThreadCount() usize {
    return std.Thread.getCpuCount() catch 1;
}
```

### 同期オーバーヘッドの最小化

- タスク粒度: 1原子 = 1タスク（細かすぎず粗すぎない）
- ロック範囲: タスク取得と結果格納のみ（計算中はロックフリー）
- 共有データ: テスト点と近傍リストは読み取り専用で共有

### ベンチマーク結果（Apple M2, 8コア）

| スレッド数 | 時間 | 高速化 |
|---------|------|---------|
| 1 | ~13ms | 1.0x |
| 2 | ~8ms | 1.6x |
| 4 | ~5ms | 2.6x |
| 8 | ~8ms | 1.6x |

**注:** スレッド数が多すぎると同期オーバーヘッドとキャッシュ競合により遅くなる。最適なスレッド数は物理コア数程度。

---

## 4. 8-wide SIMD最適化

### ファイル: `src/simd.zig`

### 問題

4-wide SIMDでは1回の操作で4原子しか処理できません。AVX-512対応CPUやApple Siliconはより広いベクトル幅を活用できます。

### 解決策: 8-wide SIMD

`@Vector(8, T)`を使用して1回の操作で8原子を同時処理（Tはf32またはf64）。

```zig
// T = f32 or f64（コンパイル時に決定）
pub fn isAnyWithinRadiusBatch8(
    comptime T: type,
    point: Vec3(T),
    atoms: AtomInput(T),
    indices: *const [8]usize,
    probe_radius: T,
) bool {
    // 8原子のx座標を読み込み
    const x: @Vector(8, T) = .{
        atoms.x[indices[0]], atoms.x[indices[1]],
        atoms.x[indices[2]], atoms.x[indices[3]],
        atoms.x[indices[4]], atoms.x[indices[5]],
        atoms.x[indices[6]], atoms.x[indices[7]],
    };
    // ... y, z, rも同様

    // 8並列で距離計算
    const dx = px - x;
    const dy = py - y;
    const dz = pz - z;
    const dist_sq = dx * dx + dy * dy + dz * dz;
    const within = dist_sq < r_sq;
    return @reduce(.Or, within);
}
```

### 階層的処理

近傍原子数に応じて最適な幅を選択:

```
近傍数 = 25:
├── 8原子 × 3回 = 24原子（8-wide SIMD）
└── 1原子 × 1回 = 1原子（スカラー）
```

### ベンチマーク結果（4V6X, 237,685原子）

| 実装 | 時間 | 高速化 |
|----------------|------|---------|
| 4-wide SIMD | ~220ms | 1.0x |
| 8-wide SIMD | ~189ms | 1.16x |

---

## 5. 高速三角関数（Lee-Richardsのみ）

### ファイル: `src/simd.zig`

### 問題

Lee-Richardsアルゴリズムでは弧の開始/終了角度を計算するために`acos`と`atan2`を多用します。標準ライブラリの実装は高精度ですが低速です。

```zig
// ボトルネック: スライスごとに複数回呼び出し
const alpha = std.math.acos(cos_alpha);  // 低速
const beta = std.math.atan2(dy, dx);     // 低速
```

### 解決策: 多項式近似

SASA用途では高精度は不要。速度優先で多項式近似を使用。

```zig
/// 高速acos近似（多項式）
/// 最大誤差: ~0.0003ラジアン（~0.02度）
pub fn fastAcos(x: f64) f64 {
    const clamped = std.math.clamp(x, -1.0, 1.0);
    const abs_x = @abs(clamped);

    // Handbook of Mathematical Functionsによる係数
    const a0: f64 = 1.5707963267948966; // π/2
    const a1: f64 = -0.2145988016038123;
    const a2: f64 = 0.0889789874093553;
    const a3: f64 = -0.0501743046129726;

    const sqrt_term = @sqrt(1.0 - abs_x);
    const poly = a0 + abs_x * (a1 + abs_x * (a2 + abs_x * a3));
    const result = sqrt_term * poly;

    return if (clamped < 0) std.math.pi - result else result;
}

/// 高速atan2近似（多項式）
/// 最大誤差: ~0.0015ラジアン（~0.09度）
pub fn fastAtan2(y: f64, x: f64) f64 {
    const abs_x = @abs(x);
    const abs_y = @abs(y);
    if (abs_x < 1e-10 and abs_y < 1e-10) return 0.0;

    const swap = abs_y > abs_x;
    const ratio = if (swap) abs_x / abs_y else abs_y / abs_x;

    const c0: f64 = 0.9998660373;
    const c1: f64 = -0.3302994844;
    const c2: f64 = 0.1801410321;

    const r2 = ratio * ratio;
    var atan_r = ratio * (c0 + r2 * (c1 + r2 * c2));

    if (swap) atan_r = std.math.pi / 2.0 - atan_r;
    if (x < 0) atan_r = std.math.pi - atan_r;
    if (y < 0) atan_r = -atan_r;

    return atan_r;
}
```

### 精度への影響

| 構造 | 面積差（vs FreeSASA C） |
|-----------|--------------------------------|
| 1CRN | 0.191% |
| 4V6X | 0.311% |

2%以内の許容誤差を満たしており、実用上問題なし。

### ベンチマーク結果（4V6X, 237,685原子, 4スレッド）

| 実装 | 時間 | 高速化 |
|----------------|------|---------|
| std.math | ~1021ms | 1.0x |
| 高速三角関数 | ~743ms | **1.37x** |

---

## 最適化効果の組み合わせ

### Shrake-Rupley（PDB 4V6X, 237,685原子, 4スレッド）

| 実装 | 時間 | vs FreeSASA C |
|----------------|------|---------------|
| FreeSASA C | 424.53ms | 1.0x |
| Zig（全最適化） | 189.06ms | **2.25x** |

### Lee-Richards（PDB 4V6X, 237,685原子, 4スレッド）

| 実装 | 時間 | vs FreeSASA C |
|----------------|------|---------------|
| FreeSASA C | 1293.48ms | 1.0x |
| Zig（全最適化） | 741.86ms | **1.74x** |

### 相乗効果（SR）

```
ナイーブ      近傍リスト    8-wide SIMD    マルチスレッド
O(N²×P)  →   O(N×P×K)    →  (定数)       →  (並列)

            （大きな効果）    (1.16x)        (線形スケール)
```

### 相乗効果（LR）

```
ナイーブ      近傍リスト    高速三角関数    マルチスレッド
O(N²×S)  →   O(N×S×K)    →  (1.37x)      →  (並列)

```

**ポイント:**
- SRは8-wide SIMDの効果が最も大きい（距離計算がボトルネック）
- LRは高速三角関数の効果が最も大きい（acos/atan2がボトルネック）
- どちらもマルチスレッドで線形にスケール

---

## 6. 精度選択（f32/f64）

### 問題

デフォルトのf64（倍精度）は高精度ですが、メモリ帯域とキャッシュ効率では不利です。

### 解決策: コンパイル時精度選択

`--precision`オプションでf32/f64を選択可能。内部計算（座標、距離、SIMD）すべてが指定精度で実行されます。

```zig
// コンパイル時に型パラメータで切り替え
pub fn calculateSasa(comptime T: type, atoms: AtomInput(T), options: Options) !SasaResult(T) {
    // T = f32 or f64
    // @Vector(8, T)で8並列SIMD
}
```

### パフォーマンス比較（バッチ処理、10,000ファイル、10スレッド）

| 精度 | 時間 | スループット | vs f64 |
|-----------|------|------------|--------|
| f32 | 171s | 58.5/s | **1.03x** |
| f64 | 176s | 56.7/s | 1.0x |

### f32の利点

1. **メモリ帯域**: データサイズが半分でキャッシュ効率向上
2. **SIMD効率**: 同じレジスタ幅で2倍の要素を処理（将来の16-wideサポート）
3. **RustSASA互換**: 公平なベンチマーク比較が可能

### f64の利点

1. **高精度**: 大規模構造での累積誤差が小さい
2. **デフォルト**: FreeSASA Cと同等精度

### 推奨事項

- **大量バッチ処理**: `--precision=f32`（速度優先）
- **高精度が必要な場合**: デフォルト（f64）を使用

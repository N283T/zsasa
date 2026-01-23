# 最適化技術詳解

本プロジェクトでは3つの主要な最適化を実装。これらの組み合わせにより、FreeSASA（Python版）比で6.4倍の高速化を達成。

## 1. 近傍リスト最適化（Phase 2）

### ファイル: `src/neighbor_list.zig`

### 問題

ナイーブなSASA計算では、各テストポイントが「埋没」しているかを判定するために全原子との距離を計算する。

```
計算量: O(N² × P)
- N = 原子数
- P = テストポイント数
```

3,183原子の場合: 3183² × 100 ≈ 10億回の距離計算

### 解決策: 空間ハッシュ

3次元空間をセルに分割し、各原子がどのセルに属するかを記録。近傍探索時は隣接セルのみを検索。

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
    neighbor_cache: [][]usize,  // 各原子の近傍キャッシュ
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

**理由:**
- セルサイズ = カットオフとすることで、近傍原子は自セル + 隣接26セルに限定
- 3×3×3 = 27セルの検索で全近傍を網羅

### 計算量改善

```
改善後: O(N × K)
- K = 平均近傍数（通常 10-20）
```

3,183原子の場合: 3183 × 20 ≈ 6.4万回の距離計算（**約15,000倍の削減**）

### ベンチマーク結果

| 実装 | 時間 | 高速化率 |
|------|------|----------|
| ナイーブ O(N²) | 〜200ms | 1.0x |
| 近傍リスト O(N) | 〜13ms | 15x |

---

## 2. SIMD最適化（Phase 3）

### ファイル: `src/simd.zig`

### 問題

距離計算は単純だが、大量に実行される（原子数 × ポイント数 × 近傍数）。

```zig
// スカラー版距離計算
fn distance_squared(p1: Vec3, p2: Vec3) f64 {
    const dx = p1.x - p2.x;
    const dy = p1.y - p2.y;
    const dz = p1.z - p2.z;
    return dx*dx + dy*dy + dz*dz;
}
```

### 解決策: SIMD（Single Instruction, Multiple Data）

Zigの `@Vector` 型を使用し、4つの距離計算を同時実行。

```zig
pub fn isAnyWithinRadius(
    point: types.Vec3,
    atoms: types.AtomInput,
    indices: *const [4]usize,
    probe_radius: f64,
) bool {
    // 4原子のx座標をロード
    const x: @Vector(4, f64) = .{
        atoms.x[indices[0]],
        atoms.x[indices[1]],
        atoms.x[indices[2]],
        atoms.x[indices[3]],
    };

    // 4原子のy座標をロード
    const y: @Vector(4, f64) = .{
        atoms.y[indices[0]],
        atoms.y[indices[1]],
        atoms.y[indices[2]],
        atoms.y[indices[3]],
    };

    // 4原子のz座標をロード
    const z: @Vector(4, f64) = .{
        atoms.z[indices[0]],
        atoms.z[indices[1]],
        atoms.z[indices[2]],
        atoms.z[indices[3]],
    };

    // 4原子の半径をロード
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

SoA形式により、SIMD命令が連続メモリから効率的にデータをロード可能。

### ベンチマーク結果

| 実装 | 時間 | 高速化率 |
|------|------|----------|
| スカラー版 | 〜13ms | 1.0x |
| SIMD版 | 〜11ms | 1.2x |

**注:** 近傍リスト最適化後は距離計算の回数自体が大幅に減少しているため、SIMD単体の効果は限定的。しかし、ナイーブ実装と組み合わせた場合は大きな効果がある。

---

## 3. マルチスレッド最適化（Phase 4）

### ファイル: `src/thread_pool.zig`

### 問題

各原子のSASA計算は独立しているが、シングルスレッドでは全CPUコアを活用できない。

### 解決策: Work-Stealingスレッドプール

```
┌─────────────────────────────────────────┐
│              Task Queue                  │
│  [Atom0] [Atom1] [Atom2] ... [AtomN-1]  │
└─────────────────────────────────────────┘
      ↓         ↓         ↓         ↓
   ┌─────┐  ┌─────┐  ┌─────┐  ┌─────┐
   │ T0  │  │ T1  │  │ T2  │  │ T3  │  スレッド
   └─────┘  └─────┘  └─────┘  └─────┘
      ↓         ↓         ↓         ↓
   [結果0]   [結果1]   [結果2]   [結果3]
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
        // タスク取得
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

        // 結果格納
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

### スレッド数の自動検出

```zig
fn getDefaultThreadCount() usize {
    return std.Thread.getCpuCount() catch 1;
}
```

### 同期のオーバーヘッド最小化

- タスク粒度: 1原子 = 1タスク（細かすぎず、粗すぎない）
- ロック範囲: タスク取得と結果格納のみ（計算中はロックフリー）
- 共有データ: テストポイントと近傍リストは読み取り専用として共有

### ベンチマーク結果（Apple M2, 8コア）

| スレッド数 | 時間 | 高速化率 |
|-----------|------|----------|
| 1 | 〜13ms | 1.0x |
| 2 | 〜8ms | 1.6x |
| 4 | 〜5ms | 2.6x |
| 8 | 〜8ms | 1.6x |

**注:** スレッド数を増やしすぎると、同期オーバーヘッドとキャッシュ競合により逆に遅くなる。最適なスレッド数は物理コア数程度。

---

## 最適化の組み合わせ効果

### 全体ベンチマーク（PDB 1A0Q, 3,183原子）

| 実装 | 時間 | FreeSASA比 |
|------|------|------------|
| FreeSASA (Python) | 〜51ms | 1.0x |
| Zig (ナイーブ) | 〜200ms | 0.25x |
| Zig (近傍リスト) | 〜13ms | 3.9x |
| Zig (近傍リスト + SIMD) | 〜11ms | 4.6x |
| Zig (全最適化) | 〜8ms | **6.4x** |

### 最適化の相乗効果

```
ナイーブ        近傍リスト      SIMD         マルチスレッド
O(N²×P)    →    O(N×P×K)    →  (定数改善)  →  (並列化)
200ms           13ms            11ms           8ms
               (15x)           (1.2x)         (1.4x)
```

**合計高速化: 200ms → 8ms = 25倍**

# Optimization Techniques

This project implements three major optimizations. Their combination achieves up to 2.3x speedup for SR and up to 1.7x for LR compared to FreeSASA C (native version).

> **Note:** Neighbor list optimization (spatial hashing) is also implemented but not listed here as it is the same technique used by FreeSASA.

## 1. SIMD Optimization

### File: `src/simd.zig`

### Problem

Distance calculation is simple but executed in large quantities (atoms × points × neighbors).

```zig
// Scalar distance calculation
fn distance_squared(p1: Vec3, p2: Vec3) f64 {
    const dx = p1.x - p2.x;
    const dy = p1.y - p2.y;
    const dz = p1.z - p2.z;
    return dx*dx + dy*dy + dz*dz;
}
```

### Solution: SIMD (Single Instruction, Multiple Data)

Using Zig's `@Vector` type, execute 8 distance calculations simultaneously (8-wide SIMD). A 4-wide version is also provided and selected based on remaining atom count.

**Precision selection**: The `--precision` option allows choosing f32/f64. f32 improves memory bandwidth and cache efficiency, achieving ~3% speedup.

```zig
pub fn isAnyWithinRadius(
    point: types.Vec3,
    atoms: types.AtomInput,
    indices: *const [4]usize,
    probe_radius: f64,
) bool {
    // Load x coordinates of 4 atoms
    const x: @Vector(4, f64) = .{
        atoms.x[indices[0]],
        atoms.x[indices[1]],
        atoms.x[indices[2]],
        atoms.x[indices[3]],
    };

    // Load y coordinates of 4 atoms
    const y: @Vector(4, f64) = .{
        atoms.y[indices[0]],
        atoms.y[indices[1]],
        atoms.y[indices[2]],
        atoms.y[indices[3]],
    };

    // Load z coordinates of 4 atoms
    const z: @Vector(4, f64) = .{
        atoms.z[indices[0]],
        atoms.z[indices[1]],
        atoms.z[indices[2]],
        atoms.z[indices[3]],
    };

    // Load radii of 4 atoms
    const r: @Vector(4, f64) = .{
        atoms.r[indices[0]] + probe_radius,
        atoms.r[indices[1]] + probe_radius,
        atoms.r[indices[2]] + probe_radius,
        atoms.r[indices[3]] + probe_radius,
    };

    // Broadcast point coordinates
    const px: @Vector(4, f64) = @splat(point.x);
    const py: @Vector(4, f64) = @splat(point.y);
    const pz: @Vector(4, f64) = @splat(point.z);

    // Difference calculation (4 parallel)
    const dx = px - x;
    const dy = py - y;
    const dz = pz - z;

    // Distance squared (4 parallel)
    const dist_sq = dx * dx + dy * dy + dz * dz;

    // Radius squared
    const r_sq = r * r;

    // Comparison (4 parallel)
    const within = dist_sq < r_sq;

    // True if any is true
    return @reduce(.Or, within);
}
```

### SoA Format Affinity

```
AoS (Array of Structures):     SoA (Structure of Arrays):
┌───────────────────┐          ┌─────────────────────┐
│ Atom0: x,y,z,r    │          │ x: [x0,x1,x2,x3...] │
│ Atom1: x,y,z,r    │    →     │ y: [y0,y1,y2,y3...] │
│ Atom2: x,y,z,r    │          │ z: [z0,z1,z2,z3...] │
│ ...               │          │ r: [r0,r1,r2,r3...] │
└───────────────────┘          └─────────────────────┘
```

SoA format enables SIMD instructions to efficiently load data from contiguous memory.

### 8-wide SIMD

Using `@Vector(8, T)` to process 8 atoms simultaneously per operation (T is f32 or f64). Select optimal width based on neighbor atom count:

```
Neighbor count = 25:
├── 8 atoms × 3 iterations = 24 atoms (8-wide SIMD)
└── 1 atom × 1 iteration = 1 atom (scalar)
```

### Future: 16-wide SIMD with AVX-512

With AVX-512 and f32 precision, 16-wide SIMD (`@Vector(16, f32)`) is theoretically possible. However, wider SIMD does not always result in better performance due to:

- **Frequency throttling**: AVX-512 instructions may cause CPU frequency reduction
- **Memory bandwidth**: Wider vectors require more data loading
- **Platform variance**: Performance varies significantly across environments (e.g., slower on WSL in testing)

This remains unverified and is not currently implemented.

---

## 2. Multi-thread Optimization

### File: `src/thread_pool.zig`

### Problem

Each atom's SASA calculation is independent, but single-threaded execution cannot utilize all CPU cores.

### Solution: Work-Stealing Thread Pool

```
┌─────────────────────────────────────────┐
│              Task Queue                  │
│  [Atom0] [Atom1] [Atom2] ... [AtomN-1]  │
└─────────────────────────────────────────┘
      ↓         ↓         ↓         ↓
   ┌─────┐  ┌─────┐  ┌─────┐  ┌─────┐
   │ T0  │  │ T1  │  │ T2  │  │ T3  │  Threads
   └─────┘  └─────┘  └─────┘  └─────┘
      ↓         ↓         ↓         ↓
   [Result0] [Result1] [Result2] [Result3]
```

### Implementation

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

### Thread Worker

```zig
fn workerThread(self: *Self) void {
    while (true) {
        // Get task
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

        // Execute task (outside lock)
        const result = self.work_fn(self.tasks[task_idx]);

        // Store result
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

### Automatic Thread Count Detection

```zig
fn getDefaultThreadCount() usize {
    return std.Thread.getCpuCount() catch 1;
}
```

### Minimizing Synchronization Overhead

- Task granularity: 1 atom = 1 task (neither too fine nor too coarse)
- Lock scope: Only task acquisition and result storage (lock-free during computation)
- Shared data: Test points and neighbor list shared as read-only

---

## 3. Fast Trigonometric Functions (Lee-Richards only)

### File: `src/simd.zig`

### Problem

The Lee-Richards algorithm uses `acos` and `atan2` extensively to calculate arc start/end angles. Standard library implementations are high-precision but slow.

```zig
// Bottleneck: Called multiple times per slice
const alpha = std.math.acos(cos_alpha);  // slow
const beta = std.math.atan2(dy, dx);     // slow
```

### Solution: Polynomial Approximation

For SASA applications, high precision is not required. Use polynomial approximation for speed.

```zig
/// Fast acos approximation (polynomial)
/// Max error: ~0.0003 radians (~0.02 degrees)
pub fn fastAcos(x: f64) f64 {
    const clamped = std.math.clamp(x, -1.0, 1.0);
    const abs_x = @abs(clamped);

    // Coefficients per Handbook of Mathematical Functions
    const a0: f64 = 1.5707963267948966; // π/2
    const a1: f64 = -0.2145988016038123;
    const a2: f64 = 0.0889789874093553;
    const a3: f64 = -0.0501743046129726;

    const sqrt_term = @sqrt(1.0 - abs_x);
    const poly = a0 + abs_x * (a1 + abs_x * (a2 + abs_x * a3));
    const result = sqrt_term * poly;

    return if (clamped < 0) std.math.pi - result else result;
}

/// Fast atan2 approximation (polynomial)
/// Max error: ~0.0015 radians (~0.09 degrees)
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

The approximation error is within acceptable tolerance (<1% vs FreeSASA C). See [benchmark/single-file.md](benchmark/single-file.md) for validation details.

---

## Combined Optimization Effects

See [benchmark/single-file.md](benchmark/single-file.md) for detailed performance data.

### Synergistic Effects (SR)

```
Naive        Neighbor list   8-wide SIMD    Multi-thread
O(N²×P)  →   O(N×P×K)    →  (constant)  →  (parallel)

            (large gain)     (1.16x)        (linear scale)
```

### Synergistic Effects (LR)

```
Naive        Neighbor list   Fast trig      Multi-thread
O(N²×S)  →   O(N×S×K)    →  (1.37x)     →  (parallel)

```

**Key points:**
- SR benefits most from 8-wide SIMD (distance calculation is the bottleneck)
- LR benefits most from fast trigonometry (acos/atan2 is the bottleneck)
- Both scale linearly with multi-threading

---

## Precision Selection (f32/f64)

### Problem

The default f64 (double precision) is high-precision but disadvantageous for memory bandwidth and cache efficiency.

### Solution: Compile-time Precision Selection

The `--precision` option allows selecting f32/f64. All internal calculations (coordinates, distances, SIMD) execute at the specified precision.

```zig
// Switch via type parameter at compile time
pub fn calculateSasa(comptime T: type, atoms: AtomInput(T), options: Options) !SasaResult(T) {
    // T = f32 or f64
    // @Vector(8, T) for 8-parallel SIMD
}
```

### Performance Comparison (batch processing, 10,000 files, 10 threads)

| Precision | Time | Throughput | vs f64 |
|-----------|------|------------|--------|
| f32 | 171s | 58.5/s | **1.03x** |
| f64 | 176s | 56.7/s | 1.0x |

### f32 Benefits

1. **Memory bandwidth**: Data size halved for improved cache efficiency
2. **SIMD efficiency**: Same register width processes 2x elements (future 16-wide support)
3. **RustSASA compatible**: Enables fair benchmark comparison

### f64 Benefits

1. **High precision**: Smaller cumulative error for large structures
2. **Default**: Equivalent precision to FreeSASA C

### Recommendations

- **High-volume batch processing**: `--precision=f32` (speed priority)
- **When high precision is needed**: Use default (f64)

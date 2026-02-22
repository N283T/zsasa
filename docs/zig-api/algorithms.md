# Algorithms

Two SASA algorithms are available: **Shrake-Rupley** (test point method) and **Lee-Richards** (slice method).

Both support:
- Single-threaded and parallel execution
- f64 (default) and f32 precision
- Automatic SIMD optimization (AVX-512/AVX2/NEON/scalar)

---

## Shrake-Rupley

Import via `zsasa.shrake_rupley`.

The Shrake-Rupley algorithm distributes test points on a sphere around each atom and counts how many are accessible (not buried inside neighboring atoms).

### calculateSasa

```zig
pub fn calculateSasa(
    allocator: std.mem.Allocator,
    input: AtomInput,
    config: Config,
) !SasaResult
```

Single-threaded SASA calculation.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `allocator` | `Allocator` | Memory allocator |
| `input` | `AtomInput` | Atom coordinates and radii |
| `config` | `Config` | Algorithm parameters (default: `n_points=100`, `probe_radius=1.4`) |

**Returns:** `SasaResult` with `total_area` and per-atom `atom_areas`.

**Example:**

```zig
const zsasa = @import("zsasa");

var atoms = try parser.parseFile("protein.pdb");
defer atoms.deinit();

var result = try zsasa.shrake_rupley.calculateSasa(allocator, atoms, .{});
defer result.deinit();

std.debug.print("Total: {d:.2} A^2\n", .{result.total_area});
for (result.atom_areas, 0..) |area, i| {
    std.debug.print("Atom {d}: {d:.2} A^2\n", .{ i, area });
}
```

### calculateSasaParallel

```zig
pub fn calculateSasaParallel(
    allocator: std.mem.Allocator,
    input: AtomInput,
    config: Config,
    n_threads: usize,
) !SasaResult
```

Multi-threaded SASA calculation with work-stealing thread pool.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `allocator` | `Allocator` | Memory allocator |
| `input` | `AtomInput` | Atom coordinates and radii |
| `config` | `Config` | Algorithm parameters |
| `n_threads` | `usize` | Thread count (`0` = auto-detect CPU cores) |

**Returns:** `SasaResult`

**Example:**

```zig
// Auto-detect thread count
var result = try zsasa.shrake_rupley.calculateSasaParallel(allocator, atoms, .{}, 0);
defer result.deinit();

// Explicit 4 threads
var result = try zsasa.shrake_rupley.calculateSasaParallel(allocator, atoms, .{}, 4);
defer result.deinit();
```

### f32 Variants

```zig
pub fn calculateSasaf32(
    allocator: Allocator,
    input: AtomInput,
    config: Configf32,
) !SasaResultf32

pub fn calculateSasaParallelf32(
    allocator: Allocator,
    input: AtomInput,
    config: Configf32,
    n_threads: usize,
) !SasaResultf32
```

Single-precision variants. Use less memory with slightly lower accuracy. Input coordinates are cast from f64 internally.

### Generic API

For advanced use, `ShrakeRupleyGen(T)` provides the full generic implementation:

```zig
const SR = zsasa.shrake_rupley.ShrakeRupleyGen(f32);
// Access internal methods for custom pipelines
```

---

## Lee-Richards

Import via `zsasa.lee_richards`.

The Lee-Richards algorithm slices each atom into parallel planes and computes exposed arc lengths on each slice circle.

### LeeRichardsConfig

```zig
pub const LeeRichardsConfig = struct {
    n_slices: u32 = 20,
    probe_radius: f64 = 1.4,
};
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `n_slices` | `u32` | `20` | Slices per atom diameter (higher = more accurate, slower) |
| `probe_radius` | `f64` | `1.4` | Water probe radius in Angstroms |

### calculateSasa

```zig
pub fn calculateSasa(
    allocator: std.mem.Allocator,
    input: AtomInput,
    config: LeeRichardsConfig,
) !SasaResult
```

Single-threaded Lee-Richards SASA calculation.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `allocator` | `Allocator` | Memory allocator |
| `input` | `AtomInput` | Atom coordinates and radii |
| `config` | `LeeRichardsConfig` | Algorithm parameters |

**Returns:** `SasaResult`

### calculateSasaParallel

```zig
pub fn calculateSasaParallel(
    allocator: std.mem.Allocator,
    input: AtomInput,
    config: LeeRichardsConfig,
    n_threads: usize,
) !SasaResult
```

Multi-threaded Lee-Richards SASA calculation.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `allocator` | `Allocator` | Memory allocator |
| `input` | `AtomInput` | Atom coordinates and radii |
| `config` | `LeeRichardsConfig` | Algorithm parameters |
| `n_threads` | `usize` | Thread count (`0` = auto-detect) |

**Returns:** `SasaResult`

**Example:**

```zig
var result = try zsasa.lee_richards.calculateSasaParallel(
    allocator, atoms, .{ .n_slices = 40 }, 0,
);
defer result.deinit();
```

### f32 Variants

```zig
pub fn calculateSasaf32(
    allocator: Allocator,
    input: AtomInput,
    config: LeeRichardsConfigf32,
) !SasaResultf32

pub fn calculateSasaParallelf32(
    allocator: Allocator,
    input: AtomInput,
    config: LeeRichardsConfigf32,
    n_threads: usize,
) !SasaResultf32
```

### Generic API

`LeeRichardsGen(T)` and `LeeRichardsConfigGen(T)` provide the full generic implementation.

---

## Choosing an Algorithm

| Aspect | Shrake-Rupley | Lee-Richards |
|--------|---------------|--------------|
| Method | Test points on sphere | Slicing into planes |
| Accuracy control | `n_points` (default 100) | `n_slices` (default 20) |
| Speed | Faster for small structures | Competitive for large structures |
| SIMD | 8-wide vectorized | Arc-based with fast trig |
| Recommended for | General use, trajectory analysis | High-accuracy requirements |

Both algorithms produce results within 0.01% of FreeSASA C reference implementation.

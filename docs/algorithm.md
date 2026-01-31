# SASA Algorithms

## Overview

This implementation provides two SASA calculation algorithms:

1. **Shrake-Rupley (SR)**: Test point method (1973)
2. **Lee-Richards (LR)**: Slice method (1971)

## Algorithm Comparison

| Property | Shrake-Rupley | Lee-Richards |
|----------|---------------|--------------|
| Published | 1973 | 1971 |
| Method | Test points | Slice/arc integration |
| Precision control | Point count (`--n-points`) | Slice count (`--n-slices`) |
| Complexity | O(N × P × K) | O(N × S × K) |
| Strengths | Simple, fast | Mathematically rigorous |
| Weaknesses | Statistical error | Higher computational cost |
| Recommended for | Large structures, rapid analysis | High-precision requirements |

**Legend:** N=atom count, P=point count, S=slice count, K=average neighbor count

### Performance Comparison (1A0Q, 3183 atoms, 4 threads, ReleaseFast)

| Algorithm | Zig | FreeSASA C | Speedup |
|-----------|-----|------------|---------|
| Shrake-Rupley (100 points) | 2.6ms | 4.4ms | 1.7x |
| Lee-Richards (20 slices) | 9.9ms | 16.1ms | 1.6x |

For detailed optimization techniques, see [optimizations.md](optimizations.md).

### Which to Choose?

- **Shrake-Rupley**: Default. Suitable for most use cases
- **Lee-Richards**: Recommended when:
  - Absolute precision is critical
  - Comparison with other Lee-Richards implementations is needed
  - Geometric rigor of arc calculations is required

---

# Shrake-Rupley Algorithm

## Overview

The Shrake-Rupley method is a numerical approach for calculating the Solvent Accessible Surface Area (SASA) of biomolecules such as proteins. It was proposed by Shrake and Rupley in 1973.

## Theoretical Background

### What is Solvent Accessible Surface Area?

- The trajectory traced when rolling a "probe sphere" (typically a 1.4Å radius sphere representing a water molecule) over the molecular surface
- Each atom's "extended radius" = van der Waals radius + probe radius
- Biochemically important: Used for estimating protein folding free energy and binding energy

### Algorithm Principles

1. Place test points uniformly on each atom's surface
2. Check whether each test point is inside any other atom
3. Calculate SASA from the proportion of exposed points

```
SASA_i = 4π × r_i² × (exposed points / total points)
```

## Implementation Details

### File: `src/shrake_rupley.zig`

#### Main API

```zig
pub fn calculateSasa(
    allocator: Allocator,
    atoms: types.AtomInput,
    options: SasaOptions,
) !types.SasaResult
```

**Parameters:**
- `allocator`: Memory allocator
- `atoms`: Atom coordinates and radii (SoA format)
- `options`: Calculation options

**Options struct:**
```zig
pub const SasaOptions = struct {
    n_points: u32 = 100,        // Test points (default: 100)
    probe_radius: f64 = 1.4,    // Probe radius (default: 1.4 Å)
    n_threads: usize = 0,       // Threads (0 = auto-detect)
};
```

#### Calculation Flow

```zig
pub fn calculateSasa(...) !types.SasaResult {
    // 1. Generate test points (on unit sphere)
    const test_points = try test_points_mod.generateTestPoints(
        allocator,
        options.n_points
    );
    defer allocator.free(test_points);

    // 2. Build neighbor list
    const max_radius = findMaxRadius(atoms.r);
    const cutoff = 2.0 * (max_radius + options.probe_radius);
    var neighbor_list = try NeighborList.init(allocator, atoms, cutoff);
    defer neighbor_list.deinit();

    // 3. Calculate SASA (parallel or sequential)
    if (options.n_threads > 1) {
        return calculateSasaParallel(...);
    } else {
        return calculateSasaSequential(...);
    }
}
```

### Single Atom SASA Calculation

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

    // Get neighbor atoms O(1)
    const neighbors = neighbor_list.getNeighbors(atom_idx);

    var exposed_count: u32 = 0;

    // Check each test point
    for (test_points) |tp| {
        // Place test point on atom surface
        const point = types.Vec3{
            .x = center.x + tp.x * radius,
            .y = center.y + tp.y * radius,
            .z = center.z + tp.z * radius,
        };

        // Collision check with neighbors (SIMD optimized)
        if (!isPointBuried(point, neighbors, atoms, probe_radius)) {
            exposed_count += 1;
        }
    }

    // SASA = 4πr² × (exposure ratio)
    const total_points: f64 = @floatFromInt(test_points.len);
    const exposed_ratio = @as(f64, @floatFromInt(exposed_count)) / total_points;
    return 4.0 * std.math.pi * radius * radius * exposed_ratio;
}
```

### Burial Detection (SIMD Optimized)

```zig
fn isPointBuried(
    point: types.Vec3,
    neighbors: []const usize,
    atoms: types.AtomInput,
    probe_radius: f64,
) bool {
    var i: usize = 0;

    // Process 8 atoms at a time with SIMD (8-wide)
    while (i + 8 <= neighbors.len) : (i += 8) {
        const indices = neighbors[i..][0..8];
        if (simd.isAnyWithinRadiusBatch8(point, atoms, indices, probe_radius)) {
            return true;
        }
    }

    // Process 4 atoms at a time (remainder)
    while (i + 4 <= neighbors.len) : (i += 4) {
        const indices = neighbors[i..][0..4];
        if (simd.isAnyWithinRadius(point, atoms, indices, probe_radius)) {
            return true;
        }
    }

    // Sequential processing for remaining atoms
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

## Test Point Generation

### File: `src/test_points.zig`

### Golden Section Spiral Algorithm

An efficient method for uniformly distributing points on a sphere.

```zig
pub fn generateTestPoints(allocator: Allocator, n: u32) ![]types.Vec3 {
    const points = try allocator.alloc(types.Vec3, n);
    errdefer allocator.free(points);

    const golden_angle = std.math.pi * (3.0 - @sqrt(5.0));  // ≈ 2.399963...

    for (0..n) |i| {
        const i_f: f64 = @floatFromInt(i);
        const n_f: f64 = @floatFromInt(n);

        // y coordinate: uniformly distributed from -1 to 1
        const y = 1.0 - (2.0 * i_f + 1.0) / n_f;

        // Radius in xz plane
        const r = @sqrt(1.0 - y * y);

        // Rotate by golden angle
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

**Mathematical Background of Golden Angle:**
- Golden angle = π(3 - √5) ≈ 137.5°
- Consecutive points rotate by this angle, forming a spiral pattern
- Same principle as sunflower seed arrangement
- Appears uniform from any viewing angle

**Point Count vs Precision Tradeoff:**

| Points | Precision | Relative Time |
|--------|-----------|---------------|
| 50 | Low | 0.5x |
| 100 | Medium (default) | 1.0x |
| 200 | High | 2.0x |
| 1000 | Very high | 10.0x |

## Precision Validation

### Comparison with FreeSASA

Comparison results with FreeSASA (C reference implementation):

| Structure | Atoms | FreeSASA (Å²) | Zig (Å²) | Difference |
|-----------|-------|---------------|----------|------------|
| 1A0Q | 3,183 | 18,923.28 | 19,211.19 | 1.52% |

**Causes of difference:**
1. Different test point generation algorithms (Golden Spiral vs Fibonacci Lattice)
2. Floating-point operation order differences
3. Neighbor list cutoff distance differences

The 1.52% difference is acceptable for typical SASA applications (relative comparison, trend analysis).

## Complexity Analysis

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Test point generation | O(P) | P = point count, precomputed |
| Neighbor list construction | O(N) | N = atom count |
| Neighbor lookup (1 atom) | O(1) | Constant time via spatial hash |
| Burial check (1 point) | O(K) | K = average neighbors (typically < 20) |
| Total | O(N × P × K) | Effectively O(N × P) |

**Comparison with naive implementation:**
- Naive: O(N² × P) - checks all atom pairs
- Optimized: O(N × P × K) - checks only neighbors
- 10x+ speedup for large molecules (N > 1000)

---

# Lee-Richards Algorithm

## Overview

The Lee-Richards method is a SASA calculation approach proposed by Lee and Richards in 1971. It slices atomic spheres along the Z-axis and calculates surface area by integrating exposed arcs at each slice.

## Theoretical Background

### Slice-Based Approach

```
        ╭─────╮
       ╱   ●   ╲      ← Atomic sphere
      │    │    │
   ───┼────┼────┼───  ← Slice (perpendicular to Z-axis)
      │    │    │
       ╲       ╱
        ╰─────╯

Slice cross-section:
    ╭───╮
   ╱     ╲
  │   ●   │  ← Atom cross-section (circle)
   ╲     ╱
    ╰───╯
     ↑
    Calculate exposed arc
```

At each slice:
1. Calculate the atom's cross-sectional circle
2. Identify obscuring arcs from neighbor atoms
3. Calculate exposed arc length
4. Multiply by slice thickness for surface area contribution

### Mathematical Formulation

Cross-sectional radius of atom i at slice z:
```
R'_i(z) = √(R_i² - (z - z_i)²)
```

Where R_i = atom radius + probe radius

Total surface area:
```
SASA_i = ∫ R_i × L_exposed(z) dz
```

Where L_exposed(z) is the exposed arc length at slice z.

## Implementation Details

### File: `src/lee_richards.zig`

#### Main API

```zig
pub fn calculateSasa(
    allocator: Allocator,
    input: AtomInput,
    config: LeeRichardsConfig,
) !SasaResult
```

**Config struct:**
```zig
pub const LeeRichardsConfig = struct {
    n_slices: u32 = 20,      // Slice count (default: 20)
    probe_radius: f64 = 1.4, // Probe radius (default: 1.4 Å)
};
```

### Arc Calculation

Calculate obscured arc from intersection of two circles:

```
Calculate obscured arc of circle 1 from intersection with circle 2

    ╭───╮
   ╱  1  ╲╲
  │   ●───●│ 2  ← Circle 2 obscures circle 1
   ╲     ╱╱
    ╰───╯

Obscuring angle α = acos((R'_i² + d_ij² - R'_j²) / (2 × R'_i × d_ij))
Arc center angle β = atan2(dy, dx) + π
```

### SIMD Optimization

Lee-Richards also uses 8-wide SIMD processing for neighbors:

```zig
// 8-wide SIMD: Calculate slice radii and xy distances
const rj_primes = simd.sliceRadiiBatch8(slice_z, batch_z, batch_r);
const dij_batch = simd.xyDistanceBatch8(xi, yi, batch_x, batch_y);

// 8-wide SIMD: Check circle overlap
const overlap_mask = simd.circlesOverlapBatch8(dij_batch, Ri_prime, rj_primes);
```

### Fast Trigonometric Functions

Arc angle calculations use polynomial approximations instead of standard library `acos`/`atan2`:

```zig
// Fast approximation (max error: acos ~0.02°, atan2 ~0.09°)
const alpha = simd.fastAcos(cos_alpha);
const beta = simd.fastAtan2(dy, dx) + std.math.pi;
```

This achieves approximately 37% speedup. Precision is within 0.3%, acceptable for practical use.

### Arc Merging

When obscuring arcs from multiple neighbor atoms overlap, merging is required:

```zig
fn exposedArcLength(arcs: []Arc) f64 {
    // Sort arcs by start angle
    sortArcs(arcs);

    // Merge overlaps and calculate exposed length
    var sum: f64 = arcs[0].start;
    var sup: f64 = arcs[0].end;

    for (arcs[1..]) |arc| {
        if (sup < arc.start) {
            sum += arc.start - sup;  // Gap (exposed portion)
        }
        sup = @max(sup, arc.end);
    }

    sum += TWOPI - sup;  // Exposure after last arc
    return sum;
}
```

## Slice Count vs Precision Tradeoff

| Slices | Precision | Relative Time |
|--------|-----------|---------------|
| 10 | Low | 0.5x |
| 20 | Medium (default) | 1.0x |
| 50 | High | 2.5x |
| 100 | Very high | 5.0x |

## Precision Comparison with Shrake-Rupley

Comparison on same structure (1A0Q):

| Algorithm | Parameters | Result (Å²) |
|-----------|------------|-------------|
| Shrake-Rupley | 100 points | 19,211.19 |
| Lee-Richards | 20 slices | 19,201.26 |
| Difference | | 0.05% |

Both algorithms produce very similar results.

## Complexity Analysis

| Operation | Complexity | Notes |
|-----------|------------|-------|
| Neighbor list construction | O(N) | N = atom count |
| Slice processing (1 atom) | O(S × K) | S = slice count, K = neighbor count |
| Arc sorting | O(K log K) | Insertion sort (efficient for small arrays) |
| Total | O(N × S × K) | Effectively O(N × S) |

## References

- Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms: Lysozyme and insulin. *Journal of Molecular Biology*, 79(2), 351-371.
- Lee, B., & Richards, F. M. (1971). The interpretation of protein structures: estimation of static accessibility. *Journal of Molecular Biology*, 55(3), 379-400.

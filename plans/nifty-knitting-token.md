# Phase 3: SIMD Optimization

## Goal

Use Zig's `@Vector` to parallelize distance calculations and achieve further performance improvement.

**Current**: ~53ms (optimized with neighbor list)
**Target**: ~30ms (1.5-2x improvement)

## Analysis Summary

### Hot Paths Identified

| Location | Function | Lines | Priority |
|----------|----------|-------|----------|
| `shrake_rupley.zig` | `atomSasaWithNeighbors` | 130-133 | **HIGHEST** |
| `shrake_rupley.zig` | `atomSasa` | 64-67 | Keep as reference |
| `neighbor_list.zig` | `addNeighborPairs` | 276-279 | LOW (runs once) |

### Current Distance Calculation
```zig
const dx = point.x - other_pos.x;
const dy = point.y - other_pos.y;
const dz = point.z - other_pos.z;
const dist_sq = dx * dx + dy * dy + dz * dz;
```

This runs millions of times per SASA calculation.

## Implementation Plan

### Phase 3.1: Add SIMD Vec3 Helper (`src/simd.zig`)

**New File: `src/simd.zig`**

```zig
const std = @import("std");
const types = @import("types.zig");
const Vec3 = types.Vec3;

/// SIMD-optimized batch distance squared calculation
/// Process 4 positions simultaneously using @Vector(4, f64)
pub fn distanceSquaredBatch4(
    point: Vec3,
    positions: [4]Vec3,
) [4]f64 {
    // Load into vectors
    const px: @Vector(4, f64) = @splat(point.x);
    const py: @Vector(4, f64) = @splat(point.y);
    const pz: @Vector(4, f64) = @splat(point.z);

    const ox = @Vector(4, f64){ positions[0].x, positions[1].x, positions[2].x, positions[3].x };
    const oy = @Vector(4, f64){ positions[0].y, positions[1].y, positions[2].y, positions[3].y };
    const oz = @Vector(4, f64){ positions[0].z, positions[1].z, positions[2].z, positions[3].z };

    // Calculate differences
    const dx = px - ox;
    const dy = py - oy;
    const dz = pz - oz;

    // Calculate squared distances
    const dist_sq = dx * dx + dy * dy + dz * dz;

    return dist_sq;
}

/// Check if point is buried by any of 4 atoms (returns true if buried)
pub fn isPointBuriedBatch4(
    point: Vec3,
    positions: [4]Vec3,
    radii_sq: [4]f64,
) bool {
    const dist_sq = distanceSquaredBatch4(point, positions);
    const radii_v: @Vector(4, f64) = radii_sq;

    // Check if any distance < radius
    const inside = dist_sq < radii_v;
    return @reduce(.Or, inside);
}
```

### Phase 3.2: Integrate into `atomSasaWithNeighbors`

**Modify: `src/shrake_rupley.zig`**

Change the inner neighbor loop to process 4 neighbors at a time:

```zig
fn atomSasaWithNeighbors(...) f64 {
    // ...
    for (test_points_array) |test_point| {
        const point = atom_pos.add(test_point.scale(atom_radius_probe));

        var is_buried = false;
        var i: usize = 0;

        // Process 4 neighbors at a time with SIMD
        while (i + 4 <= neighbors.len) : (i += 4) {
            const batch_positions = [4]Vec3{
                positions[neighbors[i]],
                positions[neighbors[i + 1]],
                positions[neighbors[i + 2]],
                positions[neighbors[i + 3]],
            };
            const batch_radii = [4]f64{
                radii_with_probe_sq[neighbors[i]],
                radii_with_probe_sq[neighbors[i + 1]],
                radii_with_probe_sq[neighbors[i + 2]],
                radii_with_probe_sq[neighbors[i + 3]],
            };

            if (simd.isPointBuriedBatch4(point, batch_positions, batch_radii)) {
                is_buried = true;
                break;
            }
        }

        // Handle remaining neighbors (0-3) with scalar fallback
        if (!is_buried) {
            while (i < neighbors.len) : (i += 1) {
                // ... existing scalar code ...
            }
        }

        if (!is_buried) n_exposed += 1;
    }
    // ...
}
```

### Phase 3.3: Unit Tests for SIMD Functions

```zig
test "distanceSquaredBatch4 - correctness" {
    const point = Vec3{ .x = 0, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 1, .y = 0, .z = 0 },  // dist² = 1
        Vec3{ .x = 0, .y = 2, .z = 0 },  // dist² = 4
        Vec3{ .x = 0, .y = 0, .z = 3 },  // dist² = 9
        Vec3{ .x = 1, .y = 1, .z = 1 },  // dist² = 3
    };

    const result = distanceSquaredBatch4(point, positions);
    try std.testing.expectApproxEqAbs(1.0, result[0], 1e-10);
    try std.testing.expectApproxEqAbs(4.0, result[1], 1e-10);
    try std.testing.expectApproxEqAbs(9.0, result[2], 1e-10);
    try std.testing.expectApproxEqAbs(3.0, result[3], 1e-10);
}

test "isPointBuriedBatch4 - one inside" {
    // Point is inside atom[1] only
    const point = Vec3{ .x = 0.5, .y = 0, .z = 0 };
    const positions = [4]Vec3{
        Vec3{ .x = 10, .y = 0, .z = 0 },  // far
        Vec3{ .x = 0, .y = 0, .z = 0 },   // close (radius 1.0)
        Vec3{ .x = 10, .y = 0, .z = 0 },  // far
        Vec3{ .x = 10, .y = 0, .z = 0 },  // far
    };
    const radii_sq = [4]f64{ 1.0, 1.0, 1.0, 1.0 };

    try std.testing.expect(isPointBuriedBatch4(point, positions, radii_sq));
}
```

### Phase 3.4: Verification

1. All existing tests must pass
2. SASA values must be identical (< 0.01% difference)
3. Benchmark should show improvement

## Files to Modify

| File | Action |
|------|--------|
| `src/simd.zig` | **CREATE** - SIMD helper functions |
| `src/shrake_rupley.zig` | **MODIFY** - Use SIMD in atomSasaWithNeighbors |

## Success Criteria

- [x] All existing tests pass
- [x] SASA values unchanged (0% difference - identical output: 19211.19 Ų)
- [x] Performance improvement (achieved ~11ms, 4.6x faster than FreeSASA baseline)
- [x] No memory leaks

## Verification Commands

```bash
# Build and run tests
zig build test

# Benchmark
time zig build run -- examples/input_1a0q.json /tmp/output.json

# Compare with FreeSASA
uv run scripts/benchmark.py examples/1A0Q.cif.gz --runs 3
```

## Notes

- Start with 4-wide SIMD (`@Vector(4, f64)`) as it's well-supported across platforms
- AVX-512 (8-wide) can be explored later if beneficial
- Keep scalar fallback for remaining elements (neighbors.len % 4)

---
- [x] **DONE** - Phase complete

## Results

**Performance achieved**: ~11ms (4.6x faster than FreeSASA's ~50ms baseline)
**Target was**: ~30ms (1.5-2x improvement)
**Actually achieved**: ~11ms (4.6x improvement)

### Implementation Summary

1. Created `src/simd.zig` with:
   - `distanceSquaredBatch4()` - SIMD batch distance calculation using `@Vector(4, f64)`
   - `isPointBuriedBatch4()` - SIMD batch burial check with early exit

2. Modified `src/shrake_rupley.zig`:
   - Updated `atomSasaWithNeighbors()` to process 4 neighbors at a time
   - Added scalar fallback for remaining neighbors (0-3)

3. All tests pass, SASA values unchanged.

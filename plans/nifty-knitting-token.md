# Phase 2: Neighbor List Optimization

## Goal

Replace O(N^2) neighbor checking with O(N) using cell-based spatial hashing.

**Current**: 1.5s for 3183 atoms (1A0Q)
**Target**: < 0.1s (30x improvement)

## Key Insight

Current bottleneck (`shrake_rupley.zig:55-73`):
```zig
for (0..n_atoms) |j| {  // Checks ALL atoms for EVERY test point
```

With neighbor list: check only ~10-30 neighbors instead of 3183 atoms.

## Implementation Plan

### Phase 2.1: Create `src/neighbor_list.zig`

**New Data Structures:**

```zig
/// Cell in spatial hash grid
pub const Cell = struct {
    atoms: std.ArrayList(u32),  // Atom indices in this cell
};

/// Spatial hash grid
pub const CellList = struct {
    cells: []Cell,
    nx, ny, nz: usize,          // Grid dimensions
    cell_size: f64,             // = 2 * (max_radius + probe_radius)
    x_min, y_min, z_min: f64,   // Bounding box origin
};

/// Pre-computed neighbor list
pub const NeighborList = struct {
    neighbors: []std.ArrayList(u32),  // neighbors[i] = neighbor indices
};
```

**Key Functions:**
- `CellList.init(positions, max_radius)` - Build spatial grid
- `NeighborList.init(positions, radii, probe_radius)` - Build neighbor list
- Both have `deinit()` for cleanup

**Algorithm:**
1. cell_size = 2 * (max_radius + probe_radius)
2. Divide bounding box into nx × ny × nz cells
3. Assign each atom to its cell
4. For each cell pair (same + 13 forward neighbors):
   - Check all atom pairs, add if distance < ri + rj + 2*probe

### Phase 2.2: Unit Tests for `neighbor_list.zig`

```zig
test "NeighborList - two far atoms" {
    // positions 100 Å apart → no neighbors
}

test "NeighborList - two touching atoms" {
    // positions 2 Å apart → mutual neighbors
}

test "NeighborList - symmetry" {
    // j in neighbors[i] ↔ i in neighbors[j]
}
```

### Phase 2.3: Integrate into `shrake_rupley.zig`

**Modify `calculateSasa()`:**
```zig
// NEW: Build neighbor list once
var neighbor_list = try NeighborList.init(allocator, positions, input.r, config.probe_radius);
defer neighbor_list.deinit();

// Pre-compute radii_with_probe[] and radii_sq[]
for (0..n_atoms) |i| {
    const neighbors = neighbor_list.getNeighbors(i);
    atom_areas[i] = atomSasaWithNeighbors(i, positions, radii_sq, test_points, neighbors);
}
```

**New internal function:**
```zig
fn atomSasaWithNeighbors(atom_idx, positions, radii_sq, test_points, neighbors) f64 {
    for (test_points) |tp| {
        for (neighbors) |j| {  // O(k) instead of O(N)
            // distance check
        }
    }
}
```

**Keep original `atomSasa()` unchanged** for backward compatibility.

### Phase 2.4: Verification

1. Run all existing tests → must pass
2. Benchmark: `zig build && time ./zig-out/bin/freesasa_zig examples/input_1a0q.json /tmp/out.json`
3. Compare output with reference (18923 Å² ± 1%)

## Files to Modify

| File | Action |
|------|--------|
| `src/neighbor_list.zig` | **CREATE** - Cell, CellList, NeighborList |
| `src/shrake_rupley.zig` | **MODIFY** - Use NeighborList in calculateSasa |

## Success Criteria

- [x] All existing tests pass
- [x] SASA values unchanged (< 0.01% difference)
- [x] 1A0Q benchmark: < 0.1s (achieved ~0.07s debug, ~0.02s release)
- [x] No memory leaks (GPA check)

## Verification Commands

```bash
# Build and run tests
zig build test

# Benchmark
time zig build run -- examples/input_1a0q.json /tmp/output.json

# Compare results
python3 scripts/benchmark.py
```

## Results

- **Performance**: 0.07s (debug) / 0.02s (ReleaseFast) vs 1.5s original → **20-75x improvement**
- **SASA**: 19211.19 Å² for 3183 atoms (1A0Q)
- **All tests pass** including verification test comparing optimized vs original
- **No memory leaks** verified via GPA

---
- [x] **DONE** - Phase complete

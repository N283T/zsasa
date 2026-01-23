# Phase 1: Shrake-Rupley SASA Algorithm Implementation

## Overview

Port the Shrake-Rupley Solvent Accessible Surface Area (SASA) algorithm from FreeSASA (C) to Zig.
Focus on correctness first, optimization later.

## Reference

- FreeSASA source: `~/ghq/github.com/mittinatten/freesasa/src/sasa_sr.c`
- Algorithm: Golden Section Spiral test points + per-atom exposure check

## File Structure

| File | Purpose |
|------|---------|
| `src/types.zig` | Shared data structures (Vec3, AtomInput, SasaResult, Config) |
| `src/json_parser.zig` | Parse input JSON |
| `src/json_writer.zig` | Write output JSON |
| `src/test_points.zig` | Generate Golden Section Spiral points |
| `src/shrake_rupley.zig` | Core SASA calculation |
| `src/main.zig` | CLI entry point |

## Data Structures

```zig
// src/types.zig

pub const Vec3 = struct {
    x: f64,
    y: f64,
    z: f64,
};

pub const AtomInput = struct {
    x: []const f64,
    y: []const f64,
    z: []const f64,
    r: []const f64,
};

pub const SasaResult = struct {
    total_area: f64,
    atom_areas: []f64,
};

pub const Config = struct {
    n_points: u32 = 100,
    probe_radius: f64 = 1.4,
};
```

## Implementation Steps

### Step 1: Types Module
- [ ] Create `src/types.zig` with Vec3, AtomInput, SasaResult, Config

### Step 2: JSON Parser
- [ ] Create `src/json_parser.zig`
- [ ] Parse x, y, z, r arrays from JSON using `std.json`
- [ ] Test: Parse `examples/input_1a0q.json` -> 3183 atoms

### Step 3: JSON Writer
- [ ] Create `src/json_writer.zig`
- [ ] Write SasaResult to JSON format

### Step 4: Test Points Generation
- [ ] Create `src/test_points.zig`
- [ ] Implement Golden Section Spiral algorithm
- [ ] Test: All points on unit sphere (|p| = 1.0)

Algorithm:
```zig
const dlong = PI * (3.0 - sqrt(5.0));
const dz = 2.0 / n;
var z = 1.0 - dz / 2.0;
var longitude = 0.0;

for each point:
    r = sqrt(1 - z*z)
    p = (cos(longitude)*r, sin(longitude)*r, z)
    z -= dz
    longitude += dlong
```

### Step 5: Single-Atom SASA
- [ ] Create `src/shrake_rupley.zig`
- [ ] Implement `atomSasa()` - SASA for one atom
- [ ] Test: Isolated atom -> SASA = 4π(r + probe)²

Algorithm:
```
For atom i:
1. Scale test points by (r[i] + probe_radius)
2. Translate to atom position
3. For each test point:
   - Check if inside any other atom j
   - If not, count as exposed
4. SASA = 4π × r² × (exposed / total)
```

### Step 6: Full SASA Calculation
- [ ] Implement `calculateSasa()` - all atoms
- [ ] Test: Compare with reference (total ~18923 Å²)

### Step 7: CLI Entry Point
- [ ] Update `src/main.zig`
- [ ] Usage: `freesasa_zig <input.json> [output.json]`
- [ ] End-to-end test with input_1a0q.json

## Acceptance Criteria

- [ ] Parse input_1a0q.json (3183 atoms)
- [ ] Total SASA within 1% of reference (18923.28 Å²)
- [ ] All tests pass
- [ ] No memory leaks (GPA check)

## Future Optimization (Phase 2+)

- Neighbor list (cell-based spatial hashing) - O(N²) → O(N)
- SIMD with `@Vector`
- Comptime test points
- Multi-threading with `std.Thread`

---
- [ ] **DONE** - Phase complete

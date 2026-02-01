# XTC Trajectory Reader

## Overview

The XTC reader provides native support for reading GROMACS XTC trajectory files without external dependencies. This is a pure Zig implementation ported from the GROMACS xdrfile library.

## Features

- Pure Zig implementation (no C dependencies)
- Sequential frame reading
- Bit-for-bit identical output to C xdrfile library
- Memory-efficient streaming design

## XTC Format

XTC (eXternal compressed Trajectory) is GROMACS's compressed trajectory format:

| Component | Description |
|-----------|-------------|
| XDR encoding | Big-endian binary (platform-independent) |
| Compression | Custom 3D lossy compression with delta encoding |
| Precision | User-specified (typically 1000 = 0.001 nm resolution) |

### Frame Structure

Each frame contains:
- Magic number (1995)
- Number of atoms
- Step number
- Simulation time
- Box matrix (3×3)
- Compressed coordinates

### Compression Algorithm

The XTC compression exploits spatial coherence in molecular dynamics:

1. **Bounding box**: Coordinates are bounded by min/max integers
2. **Delta encoding**: Most atoms store offset from previous position
3. **Run-length encoding**: Repeated patterns (e.g., water molecules) are compressed
4. **Variable-bit packing**: Smaller deltas use fewer bits

The `run` variable persists across iterations when `flag=0`, indicating "reuse previous run length" - this is a key optimization for water molecules that often have similar coordinate patterns.

## Precision and Accuracy

### Bit-for-Bit Identical Results

The Zig implementation produces **exactly the same floating-point values** as the C xdrfile library. This has been verified by comparing hex representations:

```
C:   atom[0]: 0xBF63DD98 (-0.8901000023)
Zig: atom[0]: 0xBF63DD98 (-0.8901000023)
```

The values are not just "close" - they are identical at the binary level.

### Why Tolerance in Tests?

Despite bit-for-bit identical results, the test code uses tolerance-based comparison:

```zig
const tolerance: f32 = 0.0001;
try std.testing.expectApproxEqAbs(@as(f32, -0.8901), frame.coords[0], tolerance);
```

This is **defensive programming**, not because values differ. Reasons:

1. **Best practice**: Floating-point comparisons should generally use tolerance
2. **Future-proofing**: Implementation changes could introduce minor differences
3. **Readability**: `expectApproxEqAbs` clearly shows intent

For strict verification, bit-exact comparison is possible:

```zig
try std.testing.expectEqual(
    @as(u32, @bitCast(frame.coords[0])),
    @as(u32, 0xBF63DD98)
);
```

### Impact on SASA Calculations

Since coordinates are identical to C xdrfile:
- **SASA values will be identical** when using the same input trajectory
- No precision loss from the XTC reader
- Any differences in SASA output come from algorithm differences, not coordinate reading

## Usage

### Zig API

```zig
const xtc = @import("xtc.zig");

// Open XTC file
var reader = try xtc.XtcReader.open(allocator, "trajectory.xtc");
defer reader.close();

// Get atom count
const natoms = reader.natoms;

// Read frames
while (true) {
    var frame = reader.readFrame() catch |err| {
        if (err == xtc.XtcError.EndOfFile) break;
        return err;
    };
    defer frame.deinit(allocator);

    // Access frame data
    const step = frame.step;
    const time = frame.time;
    const box = frame.box;        // 3x3 matrix
    const coords = frame.coords;  // flat [x0,y0,z0, x1,y1,z1, ...]
    const precision = frame.precision;

    // Process coordinates...
    for (0..natoms) |i| {
        const x = coords[i * 3 + 0];
        const y = coords[i * 3 + 1];
        const z = coords[i * 3 + 2];
        // ...
    }
}
```

### XtcFrame Fields

| Field | Type | Description |
|-------|------|-------------|
| `step` | `i32` | Simulation step number |
| `time` | `f32` | Simulation time (ps) |
| `box` | `[3][3]f32` | Unit cell box matrix (nm) |
| `coords` | `[]f32` | Coordinates as flat array (nm) |
| `precision` | `f32` | Compression precision (typically 1000) |

## Performance

Tested with real MD trajectories:

| File | Atoms | Frames | Size | Result |
|------|-------|--------|------|--------|
| 1l2y.xtc | 304 | 38 | 62 KB | ✓ |
| 6qfk_A_R1.xtc | 20,391 | 1,001 | 90 MB | ✓ |
| 5ltj_A_prod_R1_fit.xtc | 11,487 | 10,001 | 511 MB | ✓ |

## Limitations

- **Read-only**: Write/compression not yet implemented
- **Sequential access**: No random frame seeking (XTC frames are variable-length)
- **Single-file**: No trajectory concatenation

## References

- GROMACS xdrfile library (BSD-2-Clause): Erik Lindahl & David van der Spoel
- XDR specification: RFC 1014 (External Data Representation)
- XTC format: GROMACS manual

# Architecture Overview

## Project Structure

```
src/
├── main.zig           # CLI entry point
├── root.zig           # Library root module
├── types.zig          # Data structure definitions
├── json_parser.zig    # JSON input parsing and validation
├── json_writer.zig    # Output formats (JSON/CSV)
├── test_points.zig    # Test point generation via Golden Section Spiral
├── neighbor_list.zig  # Neighbor search with spatial hashing
├── simd.zig           # SIMD optimization (8-wide distance calculation, fast trigonometry)
├── thread_pool.zig    # Generic thread pool
├── shrake_rupley.zig  # Shrake-Rupley algorithm
└── lee_richards.zig   # Lee-Richards algorithm
```

## Module Dependencies

```
main.zig
    │
    ├── json_parser.zig ──► types.zig
    │
    ├── shrake_rupley.zig
    │       │
    │       ├── test_points.zig ──► types.zig
    │       │
    │       ├── neighbor_list.zig ──► types.zig
    │       │
    │       ├── simd.zig ──► types.zig
    │       │
    │       └── thread_pool.zig
    │
    └── json_writer.zig ──► types.zig
```

## Data Flow

```
┌─────────────────┐
│   Input JSON    │
│ (x, y, z, r)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  json_parser    │  Validation
│  - Array length │  - Coordinates are finite
│  - Radius check │  - Radius > 0, ≤ 100
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   AtomInput     │  Internal data structure
│   - x: []f64    │
│   - y: []f64    │
│   - z: []f64    │
│   - r: []f64    │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────────────┐
│           shrake_rupley                  │
│  ┌─────────────────┐                    │
│  │  test_points    │ Test point generation │
│  │  (precomputed)  │ Golden Section     │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │  neighbor_list  │ Build neighbor list │
│  │  (O(N) build)   │ Spatial hashing    │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │  thread_pool    │ Parallel processing │
│  │  (per atom)     │ Work-stealing      │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │     simd        │ Distance calculation │
│  │  (@Vector 8)    │ Process 8 atoms    │
│  └─────────────────┘                    │
└────────────────────┬────────────────────┘
                     │
                     ▼
┌─────────────────┐
│   SasaResult    │
│   - total_area  │
│   - atom_areas  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  json_writer    │  Output formats
│  - JSON (pretty)│
│  - JSON (compact)│
│  - CSV          │
└─────────────────┘
```

## Key Type Definitions (types.zig)

### Vec3

3D vector. Used for atom coordinates and test points. Supports precision parameters (f32/f64).

```zig
pub fn Vec3(comptime T: type) type {
    return struct {
        x: T,
        y: T,
        z: T,

        pub fn sub(self: @This(), other: @This()) @This() { ... }
        pub fn scale(self: @This(), s: T) @This() { ... }
        pub fn add(self: @This(), other: @This()) @This() { ... }
        pub fn lengthSquared(self: @This()) T { ... }
    };
}

// Usage examples
const Vec3f32 = Vec3(f32);  // f32 precision
const Vec3f64 = Vec3(f64);  // f64 precision (default)
```

### AtomInput

Input data structure. Uses SoA (Structure of Arrays) format for optimized memory and cache efficiency. Supports precision parameters (f32/f64).

```zig
pub fn AtomInput(comptime T: type) type {
    return struct {
        x: []const T,
        y: []const T,
        z: []const T,
        r: []const T,

        pub fn len(self: @This()) usize {
            return self.x.len;
        }
    };
}

// Usage examples
const AtomInputf32 = AtomInput(f32);  // f32 precision (fast)
const AtomInputf64 = AtomInput(f64);  // f64 precision (default)
```

**Benefits of SoA format:**
- High affinity with SIMD processing
- Contiguous memory access patterns
- Good cache line efficiency
- With f32, data size is halved for further efficiency gains

### SasaResult

Holds computation results.

```zig
pub const SasaResult = struct {
    total_area: f64,
    atom_areas: []f64,
    allocator: Allocator,

    pub fn deinit(self: *SasaResult) void {
        self.allocator.free(self.atom_areas);
    }
};
```

## Memory Management

Uses Zig's explicit allocator pattern. All dynamic memory allocations are done through an `Allocator`.

```zig
// Typical pattern
pub fn calculate(allocator: Allocator, input: AtomInput) !SasaResult {
    const atom_areas = try allocator.alloc(f64, input.len());
    errdefer allocator.free(atom_areas);  // Free on error

    // ... calculation ...

    return SasaResult{
        .total_area = total,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}
```

**Benefits:**
- Memory leak prevention (trackable at compile time)
- Leak detection with `std.testing.allocator` during tests
- Easy substitution with custom allocators (arena, etc.)

## Error Handling

Uses Zig's error union types. All errors can propagate to the caller.

```zig
pub const ParseError = error{
    InvalidJson,
    MissingField,
    ArrayLengthMismatch,
    EmptyInput,
};

pub const ValidationError = error{
    InvalidRadius,
    InvalidCoordinate,
};

pub fn parseAtomInput(allocator: Allocator, json: []const u8) !AtomInput {
    // Errors automatically propagate to caller
    const parsed = try std.json.parseFromSlice(...);
    // ...
}
```

## Thread Safety

- `thread_pool.zig`: Inter-thread synchronization implemented with `std.Thread.Mutex` and `std.Thread.Condition`
- `shrake_rupley.zig`: Each atom's calculation is independent, enabling parallel execution
- Shared state: Test point array and neighbor list are shared as read-only

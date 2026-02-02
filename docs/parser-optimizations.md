# PDB Parser Optimizations

This document describes the optimizations applied to the PDB parser that achieved a **32% reduction in parse time**.

## Overview

The PDB parser was optimized through three complementary techniques:

| Technique | Description | Impact |
|-----------|-------------|--------|
| FixedString4 | Fixed-size string type for atom names | Eliminates per-atom string allocations |
| Fast coordinate parser | Custom parser for PDB coordinate fields | Avoids `std.fmt.parseFloat` overhead |
| Pre-allocation | ArrayList capacity reservation | Reduces reallocation during parsing |

### Benchmark Results

Tested with a 17,331 atom structure:

| Version | Parse Time | Improvement |
|---------|------------|-------------|
| Before | 2.78 ms | - |
| After | 1.90 ms | **32% faster** |

---

## 1. FixedString4 Type

### File: `src/types.zig`

### Problem

PDB atom names are at most 4 characters (e.g., "CA", "N", "OG1"). Using `[]const u8` slices requires:
- String interning or allocation per atom
- Pointer indirection for access
- Memory fragmentation

### Solution

A fixed-size 5-byte struct that stores the string inline:

```zig
pub const FixedString4 = struct {
    data: [4]u8 = .{ 0, 0, 0, 0 },
    len: u8 = 0,

    pub fn fromSlice(s: []const u8) FixedString4 {
        var result = FixedString4{};
        const copy_len: u8 = @intCast(@min(s.len, 4));
        @memcpy(result.data[0..copy_len], s[0..copy_len]);
        result.len = copy_len;
        return result;
    }

    pub fn slice(self: *const FixedString4) []const u8 {
        return self.data[0..self.len];
    }
};
```

### Benefits

1. **Zero allocation**: String data stored inline in the struct
2. **Cache-friendly**: Contiguous memory layout
3. **Copy semantics**: No pointer management needed
4. **5 bytes vs 16 bytes**: Smaller than a slice (pointer + length)

### Usage

```zig
// AtomInput fields use FixedString4
pub const AtomInput = struct {
    atom_name: []FixedString4,  // was: [][]const u8
    res_name: []FixedString4,   // was: [][]const u8
    chain_id: []FixedString4,   // was: [][]const u8
    // ... coordinates, radii, etc.
};

// Accessing the string value
const name = atom.atom_name[i].slice();  // returns []const u8
```

---

## 2. Fast Coordinate Parser

### File: `src/pdb_parser.zig`

### Problem

PDB coordinates are in fixed-width format (8.3):
```
ATOM      1  N   ALA A   1      12.345  67.890 -12.345  1.00 20.00
                                ^^^^^^^^^^^^^^^^^^^^^^^^
                                Columns 31-38, 39-46, 47-54
```

Using `std.fmt.parseFloat` for every coordinate is slow because:
- Generic parser handles all float formats (scientific notation, NaN, etc.)
- Unnecessary validation for well-formed PDB files
- Function call overhead per coordinate

### Solution

A custom parser optimized for PDB's fixed-width decimal format:

```zig
fn parseCoordinate(field: []const u8) ?f64 {
    // Skip leading whitespace
    var start: usize = 0;
    const len = field.len;
    while (start < len and field[start] == ' ') : (start += 1) {}
    if (start >= len) return null;

    // Handle sign
    var negative = false;
    if (field[start] == '-') {
        negative = true;
        start += 1;
    } else if (field[start] == '+') {
        start += 1;
    }
    if (start >= len) return null;

    // Parse integer part with overflow detection
    var int_part: i64 = 0;
    var has_int_digits = false;
    while (start < len and field[start] >= '0' and field[start] <= '9') : (start += 1) {
        has_int_digits = true;
        const mul_result = @mulWithOverflow(int_part, 10);
        if (mul_result[1] != 0) return null;  // overflow
        const add_result = @addWithOverflow(mul_result[0], @as(i64, field[start] - '0'));
        if (add_result[1] != 0) return null;  // overflow
        int_part = add_result[0];
    }

    // Parse fractional part
    var frac_part: i64 = 0;
    var frac_divisor: f64 = 1.0;
    var has_frac_digits = false;

    if (start < len and field[start] == '.') {
        start += 1;
        while (start < len and field[start] >= '0' and field[start] <= '9') : (start += 1) {
            has_frac_digits = true;
            const mul_result = @mulWithOverflow(frac_part, 10);
            if (mul_result[1] != 0) break;
            const add_result = @addWithOverflow(mul_result[0], @as(i64, field[start] - '0'));
            if (add_result[1] != 0) break;
            frac_part = add_result[0];
            frac_divisor *= 10.0;
        }
    }

    // Reject sign-only input (e.g., "-" or "+")
    if (!has_int_digits and !has_frac_digits) return null;

    // Combine parts
    var result = @as(f64, @floatFromInt(int_part)) + @as(f64, @floatFromInt(frac_part)) / frac_divisor;
    if (negative) result = -result;

    return result;
}
```

### Key Features

1. **Overflow detection**: Uses `@mulWithOverflow` and `@addWithOverflow` for safe integer arithmetic
2. **Early termination**: Stops at first non-digit character
3. **Sign-only rejection**: Returns `null` for invalid input like "-" or "+"
4. **No heap allocation**: All computation on stack

---

## 3. Pre-allocation

### File: `src/pdb_parser.zig`

### Problem

`ArrayList.append()` may trigger reallocation and copying when capacity is exceeded. For large PDB files, this happens multiple times during parsing.

### Solution

Estimate atom count from file size and pre-allocate:

```zig
pub fn parse(allocator: Allocator, source: []const u8) !AtomInput {
    // Estimate: ~80 bytes per ATOM line in PDB format
    const estimated_atoms = source.len / 80;

    var x_list = std.ArrayList(f64).init(allocator);
    try x_list.ensureTotalCapacity(estimated_atoms);

    var y_list = std.ArrayList(f64).init(allocator);
    try y_list.ensureTotalCapacity(estimated_atoms);

    // ... repeat for all arrays
}
```

### Benefits

1. **Single allocation**: Arrays sized correctly from the start
2. **No copying**: No reallocation during parsing
3. **Predictable performance**: Consistent parse times regardless of file size

---

## Combined Effect

The three optimizations work synergistically:

```
File read → Pre-allocate arrays → Parse lines → Store FixedString4 + coordinates
                    ↓                  ↓                    ↓
              No reallocation    Fast coordinate      Zero string
                                    parser            allocation
```

For a typical PDB file:
- **FixedString4**: Eliminates ~3 string allocations per atom
- **Fast parser**: ~3 coordinates parsed per atom
- **Pre-allocation**: Eliminates O(log N) reallocations

The 32% improvement primarily comes from reduced allocation overhead and the faster coordinate parser.

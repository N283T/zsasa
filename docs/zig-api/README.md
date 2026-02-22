# Zig API Reference

API reference for using zsasa as a Zig library dependency.

## Installation

Add zsasa to your `build.zig.zon`:

```bash
zig fetch --save git+https://github.com/N283T/zsasa.git
```

Then in `build.zig`:

```zig
const zsasa_dep = b.dependency("zsasa", .{
    .target = target,
    .optimize = optimize,
});
exe.root_module.addImport("zsasa", zsasa_dep.module("zsasa"));
```

## Quick Start

```zig
const std = @import("std");
const zsasa = @import("zsasa");

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // Parse a PDB file
    var parser = zsasa.pdb_parser.PdbParser.init(allocator);
    var atoms = try parser.parseFile("protein.pdb");
    defer atoms.deinit();

    // Calculate SASA (Shrake-Rupley, parallel)
    var result = try zsasa.shrake_rupley.calculateSasaParallel(
        allocator, atoms, .{}, 0,  // 0 = auto-detect threads
    );
    defer result.deinit();

    std.debug.print("Total SASA: {d:.2} A^2\n", .{result.total_area});
}
```

## Modules

| Module | Description | Reference |
|--------|-------------|-----------|
| `types` | Core types: `AtomInput`, `Config`, `SasaResult`, `Vec3` | [types.md](types.md) |
| `shrake_rupley` | Shrake-Rupley algorithm | [algorithms.md](algorithms.md) |
| `lee_richards` | Lee-Richards algorithm | [algorithms.md](algorithms.md) |
| `pdb_parser` | PDB format parser | [parsers.md](parsers.md) |
| `mmcif_parser` | mmCIF format parser | [parsers.md](parsers.md) |
| `json_parser` | JSON format parser | [parsers.md](parsers.md) |
| `classifier` | Atom radius classifier | [classifier.md](classifier.md) |
| `analysis` | Per-residue aggregation, RSA, polar summary | [analysis.md](analysis.md) |

## Typical Workflow

```
Parse structure file  →  (Optional) Reclassify  →  Calculate SASA  →  Analyze results
  (pdb_parser)            (classifier)              (shrake_rupley)     (analysis)
  (mmcif_parser)                                    (lee_richards)
  (json_parser)
```

1. **Parse** a structure file into `AtomInput` (PDB, mmCIF, or JSON). PDB/mmCIF parsers assign element-based vdW radii automatically.
2. **(Optional) Reclassify** atoms with a specialized classifier (ProtOr, NACCESS, OONS) for more accurate radii
3. **Calculate** SASA using Shrake-Rupley or Lee-Richards
4. **Analyze** results: per-residue aggregation, RSA, polar/nonpolar summary

## Precision

zsasa supports both f64 (default) and f32 precision via comptime generics. All core types have f32 variants:

| f64 (default) | f32 variant |
|----------------|-------------|
| `Config` | `Configf32` |
| `LeeRichardsConfig` | `LeeRichardsConfigf32` |
| `SasaResult` | `SasaResultf32` |
| `Vec3` | `Vec3f32` |
| `calculateSasa()` | `calculateSasaf32()` |
| `calculateSasaParallel()` | `calculateSasaParallelf32()` |

Both Shrake-Rupley and Lee-Richards modules provide these f32 function variants.

f32 uses less memory with slightly lower accuracy. Use f64 unless memory is a constraint.

## Memory Management

zsasa follows Zig conventions: the caller provides an allocator and is responsible for freeing returned resources via `deinit()`.

```zig
var result = try zsasa.shrake_rupley.calculateSasa(allocator, atoms, .{});
defer result.deinit();  // Always defer deinit

var atoms = try parser.parseFile("protein.pdb");
defer atoms.deinit();   // Free parsed atom data
```

## Error Handling

All fallible functions return Zig error unions. Common patterns:

```zig
// Parser errors
const atoms = parser.parseFile("missing.pdb") catch |err| switch (err) {
    error.FileReadError => std.debug.print("File not found\n", .{}),
    error.NoAtomsFound => std.debug.print("No atoms in file\n", .{}),
    else => return err,
};

// Analysis errors (missing metadata)
const residues = zsasa.analysis.aggregateByResidue(allocator, atoms, areas) catch |err| switch (err) {
    error.MissingChainInfo => std.debug.print("No chain IDs\n", .{}),
    error.MissingResidueInfo => std.debug.print("No residue names\n", .{}),
    else => return err,
};
```

const std = @import("std");

/// Input file format
pub const InputFormat = enum {
    json,
    mmcif,
    pdb,
    sdf,
};

/// Supported file extensions (lowercase, with .gz compression)
pub const supported_extensions = [_][]const u8{
    ".json",
    ".json.gz",
    ".pdb",
    ".pdb.gz",
    ".cif",
    ".cif.gz",
    ".mmcif",
    ".mmcif.gz",
    ".ent",
    ".ent.gz",
    ".sdf",
    ".sdf.gz",
    ".mol",
    ".mol.gz",
};

/// Check if filename has a supported structure file extension
pub fn isSupportedFile(name: []const u8) bool {
    // Check lowercase extensions
    for (supported_extensions) |ext| {
        if (std.mem.endsWith(u8, name, ext)) return true;
    }
    // Check uppercase extensions (.PDB, .CIF, .ENT)
    if (std.mem.endsWith(u8, name, ".PDB")) return true;
    if (std.mem.endsWith(u8, name, ".CIF")) return true;
    if (std.mem.endsWith(u8, name, ".ENT")) return true;
    if (std.mem.endsWith(u8, name, ".JSON")) return true;
    if (std.mem.endsWith(u8, name, ".MMCIF")) return true;
    if (std.mem.endsWith(u8, name, ".mmCIF")) return true;
    if (std.mem.endsWith(u8, name, ".SDF")) return true;
    if (std.mem.endsWith(u8, name, ".MOL")) return true;
    return false;
}

/// Strip .gz extension from path
fn stripGzExt(path: []const u8) []const u8 {
    if (std.mem.endsWith(u8, path, ".gz")) {
        return path[0 .. path.len - 3];
    }
    return path;
}

/// Detect input file format from extension (handles .gz compression)
pub fn detectInputFormat(path: []const u8) InputFormat {
    const base = stripGzExt(path);

    // mmCIF formats
    if (std.mem.endsWith(u8, base, ".cif")) return .mmcif;
    if (std.mem.endsWith(u8, base, ".mmcif")) return .mmcif;
    if (std.mem.endsWith(u8, base, ".CIF")) return .mmcif;
    if (std.mem.endsWith(u8, base, ".MMCIF")) return .mmcif;
    if (std.mem.endsWith(u8, base, ".mmCIF")) return .mmcif;

    // PDB formats
    if (std.mem.endsWith(u8, base, ".pdb")) return .pdb;
    if (std.mem.endsWith(u8, base, ".PDB")) return .pdb;
    if (std.mem.endsWith(u8, base, ".ent")) return .pdb;
    if (std.mem.endsWith(u8, base, ".ENT")) return .pdb;

    // SDF/MOL formats
    if (std.mem.endsWith(u8, base, ".sdf")) return .sdf;
    if (std.mem.endsWith(u8, base, ".SDF")) return .sdf;
    if (std.mem.endsWith(u8, base, ".mol")) return .sdf;
    if (std.mem.endsWith(u8, base, ".MOL")) return .sdf;

    // Default to JSON
    return .json;
}

// =============================================================================
// Tests
// =============================================================================

test "isSupportedFile accepts structure files" {
    // PDB
    try std.testing.expect(isSupportedFile("structure.pdb"));
    try std.testing.expect(isSupportedFile("structure.pdb.gz"));
    try std.testing.expect(isSupportedFile("structure.PDB"));

    // mmCIF
    try std.testing.expect(isSupportedFile("structure.cif"));
    try std.testing.expect(isSupportedFile("structure.cif.gz"));
    try std.testing.expect(isSupportedFile("structure.mmcif"));
    try std.testing.expect(isSupportedFile("structure.CIF"));
    try std.testing.expect(isSupportedFile("structure.mmCIF"));

    // ENT
    try std.testing.expect(isSupportedFile("structure.ent"));
    try std.testing.expect(isSupportedFile("structure.ENT"));

    // JSON
    try std.testing.expect(isSupportedFile("structure.json"));
    try std.testing.expect(isSupportedFile("structure.json.gz"));

    // Unsupported
    try std.testing.expect(!isSupportedFile("structure.txt"));
    try std.testing.expect(!isSupportedFile("structure.xyz"));
    try std.testing.expect(!isSupportedFile("structure.mol2"));
}

test "detectInputFormat handles plain files" {
    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("file.pdb"));
    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("file.PDB"));
    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("file.ent"));
    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("file.ENT"));

    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.cif"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.CIF"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.mmcif"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.mmCIF"));

    try std.testing.expectEqual(InputFormat.json, detectInputFormat("file.json"));
    try std.testing.expectEqual(InputFormat.json, detectInputFormat("unknown.xyz"));
}

test "detectInputFormat handles .gz compressed files" {
    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("file.pdb.gz"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.cif.gz"));
    try std.testing.expectEqual(InputFormat.json, detectInputFormat("file.json.gz"));
}

test "isSupportedFile accepts SDF/MOL files" {
    try std.testing.expect(isSupportedFile("molecule.sdf"));
    try std.testing.expect(isSupportedFile("molecule.sdf.gz"));
    try std.testing.expect(isSupportedFile("molecule.SDF"));
    try std.testing.expect(isSupportedFile("molecule.mol"));
    try std.testing.expect(isSupportedFile("molecule.mol.gz"));
    try std.testing.expect(isSupportedFile("molecule.MOL"));
}

test "detectInputFormat handles SDF/MOL files" {
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.SDF"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.MOL"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf.gz"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol.gz"));
}

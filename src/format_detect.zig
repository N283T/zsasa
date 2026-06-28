const std = @import("std");

/// Input file format
pub const InputFormat = enum {
    bcif,
    json,
    mmcif,
    pdb,
    sdf,
};

/// Supported file extensions (lowercase, with .gz/.zst compression)
pub const supported_extensions = [_][]const u8{
    ".json",
    ".json.gz",
    ".json.zst",
    ".pdb",
    ".pdb.gz",
    ".pdb.zst",
    ".cif",
    ".cif.gz",
    ".cif.zst",
    ".bcif",
    ".bcif.gz",
    ".bcif.zst",
    ".mmcif",
    ".mmcif.gz",
    ".mmcif.zst",
    ".ent",
    ".ent.gz",
    ".ent.zst",
    ".sdf",
    ".sdf.gz",
    ".sdf.zst",
    ".mol",
    ".mol.gz",
    ".mol.zst",
};

/// Check if filename has a supported structure file extension
pub fn isSupportedFile(name: []const u8) bool {
    for (supported_extensions) |ext| {
        if (std.ascii.endsWithIgnoreCase(name, ext)) return true;
    }
    return false;
}

/// Strip supported compression extension from path.
fn stripCompressionExt(path: []const u8) []const u8 {
    if (std.ascii.endsWithIgnoreCase(path, ".gz")) {
        return path[0 .. path.len - 3];
    }
    if (std.ascii.endsWithIgnoreCase(path, ".zst")) {
        return path[0 .. path.len - 4];
    }
    return path;
}

/// Detect input file format from extension (handles .gz/.zst compression)
pub fn detectInputFormat(path: []const u8) InputFormat {
    const base = stripCompressionExt(path);

    // BinaryCIF formats (detect before mmCIF)
    if (std.ascii.endsWithIgnoreCase(base, ".bcif")) return .bcif;

    // mmCIF formats
    if (std.ascii.endsWithIgnoreCase(base, ".cif")) return .mmcif;
    if (std.ascii.endsWithIgnoreCase(base, ".mmcif")) return .mmcif;

    // PDB formats
    if (std.ascii.endsWithIgnoreCase(base, ".pdb")) return .pdb;
    if (std.ascii.endsWithIgnoreCase(base, ".ent")) return .pdb;

    // SDF/MOL formats
    if (std.ascii.endsWithIgnoreCase(base, ".sdf")) return .sdf;
    if (std.ascii.endsWithIgnoreCase(base, ".mol")) return .sdf;

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
    try std.testing.expect(isSupportedFile("structure.pdb.zst"));
    try std.testing.expect(isSupportedFile("structure.PDB"));

    // mmCIF
    try std.testing.expect(isSupportedFile("structure.cif"));
    try std.testing.expect(isSupportedFile("structure.cif.gz"));
    try std.testing.expect(isSupportedFile("structure.cif.zst"));
    try std.testing.expect(isSupportedFile("structure.mmcif"));
    try std.testing.expect(isSupportedFile("structure.CIF"));
    try std.testing.expect(isSupportedFile("structure.mmCIF"));

    // ENT
    try std.testing.expect(isSupportedFile("structure.ent"));
    try std.testing.expect(isSupportedFile("structure.ENT"));

    // JSON
    try std.testing.expect(isSupportedFile("structure.json"));
    try std.testing.expect(isSupportedFile("structure.json.gz"));
    try std.testing.expect(isSupportedFile("structure.json.zst"));

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
    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("file.pdb.zst"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.cif.gz"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("file.cif.zst"));
    try std.testing.expectEqual(InputFormat.json, detectInputFormat("file.json.gz"));
    try std.testing.expectEqual(InputFormat.json, detectInputFormat("file.json.zst"));
}

test "format detection handles uppercase compressed suffixes" {
    try std.testing.expect(isSupportedFile("MODEL.PDB.GZ"));
    try std.testing.expect(isSupportedFile("MODEL.CIF.ZST"));
    try std.testing.expect(isSupportedFile("MODEL.BCIF.GZ"));
    try std.testing.expect(isSupportedFile("MODEL.SDF.ZST"));

    try std.testing.expectEqual(InputFormat.pdb, detectInputFormat("MODEL.PDB.GZ"));
    try std.testing.expectEqual(InputFormat.mmcif, detectInputFormat("MODEL.CIF.ZST"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("MODEL.BCIF.GZ"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("MODEL.SDF.ZST"));
}

test "isSupportedFile accepts SDF/MOL files" {
    try std.testing.expect(isSupportedFile("molecule.sdf"));
    try std.testing.expect(isSupportedFile("molecule.sdf.gz"));
    try std.testing.expect(isSupportedFile("molecule.sdf.zst"));
    try std.testing.expect(isSupportedFile("molecule.SDF"));
    try std.testing.expect(isSupportedFile("molecule.mol"));
    try std.testing.expect(isSupportedFile("molecule.mol.gz"));
    try std.testing.expect(isSupportedFile("molecule.mol.zst"));
    try std.testing.expect(isSupportedFile("molecule.MOL"));
}

test "isSupportedFile accepts BinaryCIF files" {
    try std.testing.expect(isSupportedFile("structure.bcif"));
    try std.testing.expect(isSupportedFile("structure.bcif.gz"));
    try std.testing.expect(isSupportedFile("structure.bcif.zst"));
    try std.testing.expect(isSupportedFile("structure.BCIF"));
}

test "detectInputFormat handles SDF/MOL files" {
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.SDF"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.MOL"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf.gz"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.sdf.zst"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol.gz"));
    try std.testing.expectEqual(InputFormat.sdf, detectInputFormat("file.mol.zst"));
}

test "detectInputFormat handles BinaryCIF files" {
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.bcif"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.bcif.gz"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.bcif.zst"));
    try std.testing.expectEqual(InputFormat.bcif, detectInputFormat("file.BCIF"));
}

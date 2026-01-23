//! OONS classifier - Ooi, Oobatake, Nemethy, Scheraga radii
//!
//! This was the default classifier in older versions of FreeSASA.
//! OONS uses larger radii for aliphatic carbons (2.00 Å) compared to NACCESS (1.87 Å).

const std = @import("std");
const classifier = @import("classifier.zig");
const AtomClass = classifier.AtomClass;
const AtomProperties = classifier.AtomProperties;

/// Atom type with radius and polarity class
const AtomType = struct {
    radius: f64,
    class: AtomClass,
};

/// OONS atom types
const atom_types = struct {
    const C_ALI = AtomType{ .radius = 2.00, .class = .apolar }; // Aliphatic carbon
    const C_ARO = AtomType{ .radius = 1.75, .class = .apolar }; // Aromatic carbon
    const C_CAR = AtomType{ .radius = 1.55, .class = .polar }; // Carbonyl carbon (polar in OONS)
    const N = AtomType{ .radius = 1.55, .class = .polar }; // Nitrogen
    const O = AtomType{ .radius = 1.40, .class = .polar }; // Oxygen
    const S = AtomType{ .radius = 2.00, .class = .polar }; // Sulfur (polar in OONS)
    const P = AtomType{ .radius = 1.80, .class = .polar }; // Phosphorus
    const SE = AtomType{ .radius = 1.90, .class = .polar }; // Selenium
    const U_POL = AtomType{ .radius = 1.50, .class = .polar }; // Unknown polar (ASX, GLX)
    const WATER = AtomType{ .radius = 1.40, .class = .polar }; // Water
};

/// ANY atoms (fallback for all residues)
/// Key format: "ANY :ATOM" (9 chars)
const any_atoms = std.StaticStringMap(AtomType).initComptime(.{
    // Backbone
    .{ "ANY :C   ", atom_types.C_CAR },
    .{ "ANY :O   ", atom_types.O },
    .{ "ANY :CA  ", atom_types.C_ALI },
    .{ "ANY :N   ", atom_types.N },
    .{ "ANY :CB  ", atom_types.C_ALI },
    .{ "ANY :OXT ", atom_types.O },

    // DNA/RNA phosphate
    .{ "ANY :P   ", atom_types.P },
    .{ "ANY :OP1 ", atom_types.O },
    .{ "ANY :OP2 ", atom_types.O },
    .{ "ANY :OP3 ", atom_types.O },
    .{ "ANY :O5' ", atom_types.O },

    // DNA/RNA sugar
    .{ "ANY :C5' ", atom_types.C_ALI },
    .{ "ANY :C4' ", atom_types.C_ARO },
    .{ "ANY :O4' ", atom_types.O },
    .{ "ANY :C3' ", atom_types.C_ARO },
    .{ "ANY :O3' ", atom_types.O },
    .{ "ANY :C2' ", atom_types.C_ARO },
    .{ "ANY :O2' ", atom_types.O },
    .{ "ANY :C1' ", atom_types.C_ARO },

    // DNA/RNA base
    .{ "ANY :N1  ", atom_types.N },
    .{ "ANY :N2  ", atom_types.N },
    .{ "ANY :N3  ", atom_types.N },
    .{ "ANY :N4  ", atom_types.N },
    .{ "ANY :N6  ", atom_types.N },
    .{ "ANY :N7  ", atom_types.N },
    .{ "ANY :N9  ", atom_types.N },
    .{ "ANY :C2  ", atom_types.C_ARO },
    .{ "ANY :C4  ", atom_types.C_ARO },
    .{ "ANY :C5  ", atom_types.C_ARO },
    .{ "ANY :C6  ", atom_types.C_ARO },
    .{ "ANY :C7  ", atom_types.C_ARO },
    .{ "ANY :C8  ", atom_types.C_ARO },
    .{ "ANY :O2  ", atom_types.O },
    .{ "ANY :O4  ", atom_types.O },
    .{ "ANY :O6  ", atom_types.O },

    // Methylation
    .{ "ANY :CM2 ", atom_types.C_ALI },
});

/// Residue-specific atom lookup
/// Key format: "RES :ATOM" (9 chars, space-padded)
const residue_atoms = std.StaticStringMap(AtomType).initComptime(.{
    // === Amino Acid Side Chains ===

    // ARG
    .{ "ARG :CG  ", atom_types.C_ALI },
    .{ "ARG :CD  ", atom_types.C_ALI },
    .{ "ARG :NE  ", atom_types.N },
    .{ "ARG :CZ  ", atom_types.C_ALI },
    .{ "ARG :NH1 ", atom_types.N },
    .{ "ARG :NH2 ", atom_types.N },

    // ASN
    .{ "ASN :CG  ", atom_types.C_CAR },
    .{ "ASN :OD1 ", atom_types.O },
    .{ "ASN :ND2 ", atom_types.N },

    // ASP
    .{ "ASP :CG  ", atom_types.C_CAR },
    .{ "ASP :OD1 ", atom_types.O },
    .{ "ASP :OD2 ", atom_types.O },

    // CYS
    .{ "CYS :SG  ", atom_types.S },

    // GLN
    .{ "GLN :CG  ", atom_types.C_ALI },
    .{ "GLN :CD  ", atom_types.C_CAR },
    .{ "GLN :OE1 ", atom_types.O },
    .{ "GLN :NE2 ", atom_types.N },

    // GLU
    .{ "GLU :CG  ", atom_types.C_ALI },
    .{ "GLU :CD  ", atom_types.C_CAR },
    .{ "GLU :OE1 ", atom_types.O },
    .{ "GLU :OE2 ", atom_types.O },

    // HIS
    .{ "HIS :CG  ", atom_types.C_ARO },
    .{ "HIS :ND1 ", atom_types.N },
    .{ "HIS :CD2 ", atom_types.C_ARO },
    .{ "HIS :NE2 ", atom_types.N },
    .{ "HIS :CE1 ", atom_types.C_ARO },

    // ILE
    .{ "ILE :CG1 ", atom_types.C_ALI },
    .{ "ILE :CG2 ", atom_types.C_ALI },
    .{ "ILE :CD1 ", atom_types.C_ALI },

    // LEU
    .{ "LEU :CG  ", atom_types.C_ALI },
    .{ "LEU :CD1 ", atom_types.C_ALI },
    .{ "LEU :CD2 ", atom_types.C_ALI },

    // LYS
    .{ "LYS :CG  ", atom_types.C_ALI },
    .{ "LYS :CD  ", atom_types.C_ALI },
    .{ "LYS :CE  ", atom_types.C_ALI },
    .{ "LYS :NZ  ", atom_types.N },

    // MET
    .{ "MET :CG  ", atom_types.C_ALI },
    .{ "MET :SD  ", atom_types.S },
    .{ "MET :CE  ", atom_types.C_ALI },

    // PHE
    .{ "PHE :CG  ", atom_types.C_ARO },
    .{ "PHE :CD1 ", atom_types.C_ARO },
    .{ "PHE :CD2 ", atom_types.C_ARO },
    .{ "PHE :CE1 ", atom_types.C_ARO },
    .{ "PHE :CE2 ", atom_types.C_ARO },
    .{ "PHE :CZ  ", atom_types.C_ARO },

    // PRO (ring carbons are aromatic-like)
    .{ "PRO :CB  ", atom_types.C_ARO },
    .{ "PRO :CG  ", atom_types.C_ARO },
    .{ "PRO :CD  ", atom_types.C_ARO },

    // SER
    .{ "SER :OG  ", atom_types.O },

    // THR
    .{ "THR :OG1 ", atom_types.O },
    .{ "THR :CG2 ", atom_types.C_ALI },

    // TRP
    .{ "TRP :CG  ", atom_types.C_ARO },
    .{ "TRP :CD1 ", atom_types.C_ARO },
    .{ "TRP :CD2 ", atom_types.C_ARO },
    .{ "TRP :NE1 ", atom_types.N },
    .{ "TRP :CE2 ", atom_types.C_ARO },
    .{ "TRP :CE3 ", atom_types.C_ARO },
    .{ "TRP :CZ2 ", atom_types.C_ARO },
    .{ "TRP :CZ3 ", atom_types.C_ARO },
    .{ "TRP :CH2 ", atom_types.C_ARO },

    // TYR
    .{ "TYR :CG  ", atom_types.C_ARO },
    .{ "TYR :CD1 ", atom_types.C_ARO },
    .{ "TYR :CD2 ", atom_types.C_ARO },
    .{ "TYR :CE1 ", atom_types.C_ARO },
    .{ "TYR :CE2 ", atom_types.C_ARO },
    .{ "TYR :CZ  ", atom_types.C_ARO },
    .{ "TYR :OH  ", atom_types.O },

    // VAL
    .{ "VAL :CG1 ", atom_types.C_ALI },
    .{ "VAL :CG2 ", atom_types.C_ALI },

    // === Non-standard Amino Acids ===

    // ASX (ambiguous ASN/ASP)
    .{ "ASX :CG  ", atom_types.C_CAR },
    .{ "ASX :XD1 ", atom_types.U_POL },
    .{ "ASX :XD2 ", atom_types.U_POL },
    .{ "ASX :AD1 ", atom_types.U_POL },
    .{ "ASX :AD2 ", atom_types.U_POL },

    // GLX (ambiguous GLN/GLU)
    .{ "GLX :CG  ", atom_types.C_ALI },
    .{ "GLX :CD  ", atom_types.C_CAR },
    .{ "GLX :XE1 ", atom_types.U_POL },
    .{ "GLX :XE2 ", atom_types.U_POL },
    .{ "GLX :AE1 ", atom_types.U_POL },
    .{ "GLX :AE2 ", atom_types.U_POL },

    // SEC (selenocysteine)
    .{ "SEC :SE  ", atom_types.SE },
    .{ "CSE :SE  ", atom_types.SE }, // alternate name

    // MSE (selenomethionine)
    .{ "MSE :SE  ", atom_types.SE },
    .{ "MSE :CG  ", atom_types.C_ALI },
    .{ "MSE :CE  ", atom_types.C_ALI },

    // PYL (pyrrolysine)
    .{ "PYL :CG  ", atom_types.C_ALI },
    .{ "PYL :CD  ", atom_types.C_ALI },
    .{ "PYL :CE  ", atom_types.C_ALI },
    .{ "PYL :NZ  ", atom_types.N },
    .{ "PYL :O2  ", atom_types.O },
    .{ "PYL :C2  ", atom_types.C_CAR },
    .{ "PYL :CA2 ", atom_types.C_ARO },
    .{ "PYL :CB2 ", atom_types.C_ALI },
    .{ "PYL :CG2 ", atom_types.C_ARO },
    .{ "PYL :CD2 ", atom_types.C_ARO },
    .{ "PYL :CE2 ", atom_types.C_ARO },
    .{ "PYL :N2  ", atom_types.N },

    // === Capping Groups ===

    // ACE (acetyl)
    .{ "ACE :CH3 ", atom_types.C_ALI },

    // NH2 (amide)
    .{ "NH2 :NH2 ", atom_types.N },

    // === Water ===
    .{ "HOH :O   ", atom_types.WATER },
});

/// Create lookup key at runtime
fn makeKeyRuntime(residue: []const u8, atom: []const u8) [9]u8 {
    var key: [9]u8 = .{ ' ', ' ', ' ', ' ', ':', ' ', ' ', ' ', ' ' };

    // Copy residue (up to 4 chars, left-justified)
    const res_len = @min(residue.len, 4);
    for (residue[0..res_len], 0..) |c, i| {
        key[i] = c;
    }

    // Copy atom name (up to 4 chars)
    const atom_len = @min(atom.len, 4);
    for (atom[0..atom_len], 0..) |c, i| {
        key[5 + i] = c;
    }

    return key;
}

/// Get radius for a (residue, atom) pair
/// Lookup order: residue-specific → ANY fallback → element guess
pub fn getRadius(residue: []const u8, atom: []const u8) ?f64 {
    // Try residue-specific first
    const key = makeKeyRuntime(residue, atom);
    if (residue_atoms.get(&key)) |t| {
        return t.radius;
    }

    // Try ANY fallback
    const any_key = makeKeyRuntime("ANY", atom);
    if (any_atoms.get(&any_key)) |t| {
        return t.radius;
    }

    // Fall back to element guess
    return classifier.guessRadiusFromAtomName(atom);
}

/// Get polarity class for a (residue, atom) pair
pub fn getClass(residue: []const u8, atom: []const u8) AtomClass {
    // Try residue-specific first
    const key = makeKeyRuntime(residue, atom);
    if (residue_atoms.get(&key)) |t| {
        return t.class;
    }

    // Try ANY fallback
    const any_key = makeKeyRuntime("ANY", atom);
    if (any_atoms.get(&any_key)) |t| {
        return t.class;
    }

    return .unknown;
}

/// Get both radius and class
pub fn getProperties(residue: []const u8, atom: []const u8) ?AtomProperties {
    // Try residue-specific first
    const key = makeKeyRuntime(residue, atom);
    if (residue_atoms.get(&key)) |t| {
        return AtomProperties{ .radius = t.radius, .class = t.class };
    }

    // Try ANY fallback
    const any_key = makeKeyRuntime("ANY", atom);
    if (any_atoms.get(&any_key)) |t| {
        return AtomProperties{ .radius = t.radius, .class = t.class };
    }

    // Try element guess for radius only
    if (classifier.guessRadiusFromAtomName(atom)) |r| {
        return AtomProperties{ .radius = r, .class = .unknown };
    }

    return null;
}

// ============================================================================
// Tests
// ============================================================================

test "OONS backbone atoms via ANY" {
    // All residues should get backbone atoms from ANY
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("ALA", "CA").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.55), getRadius("ALA", "C").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.55), getRadius("ALA", "N").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.40), getRadius("ALA", "O").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("ALA", "CB").?, 0.001);

    // ANY works for any residue
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("GLY", "CA").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("VAL", "CA").?, 0.001);
}

test "OONS side chain atoms" {
    // PHE aromatic ring
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), getRadius("PHE", "CG").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), getRadius("PHE", "CD1").?, 0.001);

    // LEU aliphatic
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("LEU", "CG").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("LEU", "CD1").?, 0.001);

    // CYS sulfur
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("CYS", "SG").?, 0.001);
}

test "OONS polarity classes" {
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "N"));
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "O"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("ALA", "CA"));
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "C")); // C_CAR is polar in OONS
    try std.testing.expectEqual(AtomClass.polar, getClass("CYS", "SG")); // S is polar in OONS
}

test "OONS vs NACCESS radii differences" {
    // OONS has larger aliphatic carbon radius
    // NACCESS C_ALI = 1.87, OONS C_ALI = 2.00
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("ALA", "CA").?, 0.001);

    // OONS aromatic carbon is smaller
    // NACCESS C_CAR = 1.76, OONS C_ARO = 1.75
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), getRadius("PHE", "CG").?, 0.001);
}

test "OONS selenocysteine and selenomethionine" {
    // SEC
    try std.testing.expectApproxEqAbs(@as(f64, 1.90), getRadius("SEC", "SE").?, 0.001);
    // CB comes from ANY
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("SEC", "CB").?, 0.001);

    // MSE
    try std.testing.expectApproxEqAbs(@as(f64, 1.90), getRadius("MSE", "SE").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("MSE", "CG").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 2.00), getRadius("MSE", "CE").?, 0.001);
}

test "OONS nucleotides via ANY" {
    // Phosphate
    try std.testing.expectApproxEqAbs(@as(f64, 1.80), getRadius("A", "P").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.40), getRadius("A", "OP1").?, 0.001);

    // Sugar
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), getRadius("A", "C1'").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.40), getRadius("A", "O4'").?, 0.001);

    // Base
    try std.testing.expectApproxEqAbs(@as(f64, 1.55), getRadius("A", "N9").?, 0.001);
}

test "OONS water" {
    try std.testing.expectApproxEqAbs(@as(f64, 1.40), getRadius("HOH", "O").?, 0.001);
    try std.testing.expectEqual(AtomClass.polar, getClass("HOH", "O"));
}

test "OONS element-based fallback" {
    // Unknown residue - should use element guess
    const r = getRadius("UNK", "FE");
    try std.testing.expect(r != null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.26), r.?, 0.001); // Fe element
}

test "OONS getProperties" {
    const props = getProperties("PHE", "CG");
    try std.testing.expect(props != null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.75), props.?.radius, 0.001);
    try std.testing.expectEqual(AtomClass.apolar, props.?.class);
}

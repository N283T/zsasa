//! NACCESS-compatible built-in classifier.
//!
//! This module provides van der Waals radii and polarity classes compatible
//! with the NACCESS program. Data is embedded at compile time for zero-allocation
//! O(1) lookup.
//!
//! Based on FreeSASA's naccess.config (contributed by João Rodrigues).

const std = @import("std");
const classifier = @import("classifier.zig");

pub const AtomClass = classifier.AtomClass;
pub const AtomProperties = classifier.AtomProperties;

// =============================================================================
// Atom Type Definitions
// =============================================================================

/// NACCESS atom type radii and classes
const AtomType = struct {
    radius: f64,
    class: AtomClass,
};

/// All NACCESS atom types
const atom_types = struct {
    const C_ALI = AtomType{ .radius = 1.87, .class = .apolar };
    const C_CAR = AtomType{ .radius = 1.76, .class = .apolar };
    const C_NUC = AtomType{ .radius = 1.80, .class = .apolar };
    const N_AMN = AtomType{ .radius = 1.50, .class = .polar };
    const N_AMD = AtomType{ .radius = 1.65, .class = .polar };
    const N_NUC = AtomType{ .radius = 1.60, .class = .polar };
    const O = AtomType{ .radius = 1.40, .class = .polar };
    const S = AtomType{ .radius = 1.85, .class = .apolar };
    const SE = AtomType{ .radius = 1.80, .class = .apolar };
    const P = AtomType{ .radius = 1.90, .class = .apolar };
};

// =============================================================================
// Lookup Key Utilities
// =============================================================================

/// Runtime key generation for lookup (9 chars: RES:ATOM with padding)
fn makeKeyRuntime(residue: []const u8, atom: []const u8) [9]u8 {
    var key: [9]u8 = .{ ' ', ' ', ' ', ' ', ':', ' ', ' ', ' ', ' ' };
    const res_len = @min(residue.len, 4);
    const atm_len = @min(atom.len, 4);
    for (residue[0..res_len], 0..) |c, i| {
        key[i] = c;
    }
    for (atom[0..atm_len], 0..) |c, i| {
        key[5 + i] = c;
    }
    return key;
}

// =============================================================================
// Atom Lookup Maps
// Keys are formatted as "RES :ATOM" (4 char residue + colon + 4 char atom)
// =============================================================================

/// ANY residue atoms (fallback for all residues)
const any_atoms = std.StaticStringMap(AtomType).initComptime(.{
    // Backbone atoms
    .{ "ANY :C   ", atom_types.C_CAR },
    .{ "ANY :O   ", atom_types.O },
    .{ "ANY :CA  ", atom_types.C_ALI },
    .{ "ANY :N   ", atom_types.N_AMD },
    .{ "ANY :CB  ", atom_types.C_ALI },
    .{ "ANY :OXT ", atom_types.O },
    // Nucleic acid backbone
    .{ "ANY :P   ", atom_types.P },
    .{ "ANY :OP1 ", atom_types.O },
    .{ "ANY :OP2 ", atom_types.O },
    .{ "ANY :OP3 ", atom_types.O },
    .{ "ANY :O5' ", atom_types.O },
    .{ "ANY :O4' ", atom_types.O },
    .{ "ANY :O3' ", atom_types.O },
    .{ "ANY :O2' ", atom_types.O },
    .{ "ANY :C5' ", atom_types.C_NUC },
    .{ "ANY :C4' ", atom_types.C_NUC },
    .{ "ANY :C3' ", atom_types.C_NUC },
    .{ "ANY :C2' ", atom_types.C_NUC },
    .{ "ANY :C1' ", atom_types.C_NUC },
});

/// Residue-specific atoms (overrides ANY)
const residue_atoms = std.StaticStringMap(AtomType).initComptime(.{
    // ALA
    .{ "ALA :CB  ", atom_types.C_ALI },

    // ARG
    .{ "ARG :CG  ", atom_types.C_ALI },
    .{ "ARG :CD  ", atom_types.C_ALI },
    .{ "ARG :NE  ", atom_types.N_AMD },
    .{ "ARG :CZ  ", atom_types.C_CAR },
    .{ "ARG :NH1 ", atom_types.N_AMD },
    .{ "ARG :NH2 ", atom_types.N_AMD },

    // ASN
    .{ "ASN :CG  ", atom_types.C_CAR },
    .{ "ASN :OD1 ", atom_types.O },
    .{ "ASN :ND2 ", atom_types.N_AMD },

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
    .{ "GLN :NE2 ", atom_types.N_AMD },

    // GLU
    .{ "GLU :CG  ", atom_types.C_ALI },
    .{ "GLU :CD  ", atom_types.C_CAR },
    .{ "GLU :OE1 ", atom_types.O },
    .{ "GLU :OE2 ", atom_types.O },

    // GLY
    .{ "GLY :CA  ", atom_types.C_ALI },

    // HIS
    .{ "HIS :CG  ", atom_types.C_CAR },
    .{ "HIS :ND1 ", atom_types.N_AMD },
    .{ "HIS :CD2 ", atom_types.C_CAR },
    .{ "HIS :NE2 ", atom_types.N_AMD },
    .{ "HIS :CE1 ", atom_types.C_CAR },

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
    .{ "LYS :NZ  ", atom_types.N_AMN },

    // MET
    .{ "MET :CG  ", atom_types.C_ALI },
    .{ "MET :SD  ", atom_types.S },
    .{ "MET :CE  ", atom_types.C_ALI },

    // PHE
    .{ "PHE :CG  ", atom_types.C_CAR },
    .{ "PHE :CD1 ", atom_types.C_CAR },
    .{ "PHE :CD2 ", atom_types.C_CAR },
    .{ "PHE :CE1 ", atom_types.C_CAR },
    .{ "PHE :CE2 ", atom_types.C_CAR },
    .{ "PHE :CZ  ", atom_types.C_CAR },

    // PRO
    .{ "PRO :CG  ", atom_types.C_ALI },
    .{ "PRO :CD  ", atom_types.C_ALI },

    // SEC (selenocysteine) - like CYS but with SE instead of S
    .{ "SEC :CB  ", atom_types.C_ALI },
    .{ "SEC :SE  ", atom_types.SE },

    // MSE (selenomethionine) - like MET but with SE instead of S
    .{ "MSE :CG  ", atom_types.C_ALI },
    .{ "MSE :SE  ", atom_types.SE },
    .{ "MSE :CE  ", atom_types.C_ALI },

    // SER
    .{ "SER :OG  ", atom_types.O },

    // THR
    .{ "THR :OG1 ", atom_types.O },
    .{ "THR :CG2 ", atom_types.C_ALI },

    // TRP
    .{ "TRP :CG  ", atom_types.C_CAR },
    .{ "TRP :CD1 ", atom_types.C_CAR },
    .{ "TRP :CD2 ", atom_types.C_CAR },
    .{ "TRP :NE1 ", atom_types.N_AMD },
    .{ "TRP :CE2 ", atom_types.C_CAR },
    .{ "TRP :CE3 ", atom_types.C_CAR },
    .{ "TRP :CZ2 ", atom_types.C_CAR },
    .{ "TRP :CZ3 ", atom_types.C_CAR },
    .{ "TRP :CH2 ", atom_types.C_CAR },

    // TYR
    .{ "TYR :CG  ", atom_types.C_CAR },
    .{ "TYR :CD1 ", atom_types.C_CAR },
    .{ "TYR :CD2 ", atom_types.C_CAR },
    .{ "TYR :CE1 ", atom_types.C_CAR },
    .{ "TYR :CE2 ", atom_types.C_CAR },
    .{ "TYR :CZ  ", atom_types.C_CAR },
    .{ "TYR :OH  ", atom_types.O },

    // VAL
    .{ "VAL :CG1 ", atom_types.C_ALI },
    .{ "VAL :CG2 ", atom_types.C_ALI },

    // RNA Nucleotides
    // A (Adenine)
    .{ "A   :N9  ", atom_types.N_NUC },
    .{ "A   :C8  ", atom_types.C_NUC },
    .{ "A   :N7  ", atom_types.N_NUC },
    .{ "A   :C5  ", atom_types.C_NUC },
    .{ "A   :C6  ", atom_types.C_NUC },
    .{ "A   :N6  ", atom_types.N_NUC },
    .{ "A   :N1  ", atom_types.N_NUC },
    .{ "A   :C2  ", atom_types.C_NUC },
    .{ "A   :N3  ", atom_types.N_NUC },
    .{ "A   :C4  ", atom_types.C_NUC },

    // C (Cytosine)
    .{ "C   :N1  ", atom_types.N_NUC },
    .{ "C   :C2  ", atom_types.C_NUC },
    .{ "C   :O2  ", atom_types.O },
    .{ "C   :N3  ", atom_types.N_NUC },
    .{ "C   :C4  ", atom_types.C_NUC },
    .{ "C   :N4  ", atom_types.N_NUC },
    .{ "C   :C5  ", atom_types.C_NUC },
    .{ "C   :C6  ", atom_types.C_NUC },

    // G (Guanine)
    .{ "G   :N9  ", atom_types.N_NUC },
    .{ "G   :C8  ", atom_types.C_NUC },
    .{ "G   :N7  ", atom_types.N_NUC },
    .{ "G   :C5  ", atom_types.C_NUC },
    .{ "G   :C6  ", atom_types.C_NUC },
    .{ "G   :O6  ", atom_types.O },
    .{ "G   :N1  ", atom_types.N_NUC },
    .{ "G   :C2  ", atom_types.C_NUC },
    .{ "G   :N2  ", atom_types.N_NUC },
    .{ "G   :N3  ", atom_types.N_NUC },
    .{ "G   :C4  ", atom_types.C_NUC },

    // I (Inosine)
    .{ "I   :N9  ", atom_types.N_NUC },
    .{ "I   :C8  ", atom_types.C_NUC },
    .{ "I   :N7  ", atom_types.N_NUC },
    .{ "I   :C5  ", atom_types.C_NUC },
    .{ "I   :C6  ", atom_types.C_NUC },
    .{ "I   :O6  ", atom_types.O },
    .{ "I   :N1  ", atom_types.N_NUC },
    .{ "I   :C2  ", atom_types.C_NUC },
    .{ "I   :N3  ", atom_types.N_NUC },
    .{ "I   :C4  ", atom_types.C_NUC },

    // T (Thymine - RNA)
    .{ "T   :N1  ", atom_types.N_NUC },
    .{ "T   :C2  ", atom_types.C_NUC },
    .{ "T   :O2  ", atom_types.O },
    .{ "T   :N3  ", atom_types.N_NUC },
    .{ "T   :C4  ", atom_types.C_NUC },
    .{ "T   :O4  ", atom_types.O },
    .{ "T   :C5  ", atom_types.C_NUC },
    .{ "T   :C7  ", atom_types.C_NUC },
    .{ "T   :C6  ", atom_types.C_NUC },

    // U (Uracil)
    .{ "U   :N1  ", atom_types.N_NUC },
    .{ "U   :C2  ", atom_types.C_NUC },
    .{ "U   :O2  ", atom_types.O },
    .{ "U   :N3  ", atom_types.N_NUC },
    .{ "U   :C4  ", atom_types.C_NUC },
    .{ "U   :O4  ", atom_types.O },
    .{ "U   :C5  ", atom_types.C_NUC },
    .{ "U   :C6  ", atom_types.C_NUC },

    // DNA Nucleotides
    // DA (Deoxyadenosine)
    .{ "DA  :N9  ", atom_types.N_NUC },
    .{ "DA  :C8  ", atom_types.C_NUC },
    .{ "DA  :N7  ", atom_types.N_NUC },
    .{ "DA  :C5  ", atom_types.C_NUC },
    .{ "DA  :C6  ", atom_types.C_NUC },
    .{ "DA  :N6  ", atom_types.N_NUC },
    .{ "DA  :N1  ", atom_types.N_NUC },
    .{ "DA  :C2  ", atom_types.C_NUC },
    .{ "DA  :N3  ", atom_types.N_NUC },
    .{ "DA  :C4  ", atom_types.C_NUC },

    // DC (Deoxycytidine)
    .{ "DC  :N1  ", atom_types.N_NUC },
    .{ "DC  :C2  ", atom_types.C_NUC },
    .{ "DC  :O2  ", atom_types.O },
    .{ "DC  :N3  ", atom_types.N_NUC },
    .{ "DC  :C4  ", atom_types.C_NUC },
    .{ "DC  :N4  ", atom_types.N_NUC },
    .{ "DC  :C5  ", atom_types.C_NUC },
    .{ "DC  :C6  ", atom_types.C_NUC },

    // DG (Deoxyguanosine)
    .{ "DG  :N9  ", atom_types.N_NUC },
    .{ "DG  :C8  ", atom_types.C_NUC },
    .{ "DG  :N7  ", atom_types.N_NUC },
    .{ "DG  :C5  ", atom_types.C_NUC },
    .{ "DG  :C6  ", atom_types.C_NUC },
    .{ "DG  :O6  ", atom_types.O },
    .{ "DG  :N1  ", atom_types.N_NUC },
    .{ "DG  :C2  ", atom_types.C_NUC },
    .{ "DG  :N2  ", atom_types.N_NUC },
    .{ "DG  :N3  ", atom_types.N_NUC },
    .{ "DG  :C4  ", atom_types.C_NUC },

    // DI (Deoxyinosine)
    .{ "DI  :N9  ", atom_types.N_NUC },
    .{ "DI  :C8  ", atom_types.C_NUC },
    .{ "DI  :N7  ", atom_types.N_NUC },
    .{ "DI  :C5  ", atom_types.C_NUC },
    .{ "DI  :C6  ", atom_types.C_NUC },
    .{ "DI  :O6  ", atom_types.O },
    .{ "DI  :N1  ", atom_types.N_NUC },
    .{ "DI  :C2  ", atom_types.C_NUC },
    .{ "DI  :N3  ", atom_types.N_NUC },
    .{ "DI  :C4  ", atom_types.C_NUC },

    // DT (Deoxythymidine)
    .{ "DT  :N1  ", atom_types.N_NUC },
    .{ "DT  :C2  ", atom_types.C_NUC },
    .{ "DT  :O2  ", atom_types.O },
    .{ "DT  :N3  ", atom_types.N_NUC },
    .{ "DT  :C4  ", atom_types.C_NUC },
    .{ "DT  :O4  ", atom_types.O },
    .{ "DT  :C5  ", atom_types.C_NUC },
    .{ "DT  :C7  ", atom_types.C_NUC },
    .{ "DT  :C6  ", atom_types.C_NUC },

    // DU (Deoxyuridine)
    .{ "DU  :N1  ", atom_types.N_NUC },
    .{ "DU  :C2  ", atom_types.C_NUC },
    .{ "DU  :O2  ", atom_types.O },
    .{ "DU  :N3  ", atom_types.N_NUC },
    .{ "DU  :C4  ", atom_types.C_NUC },
    .{ "DU  :O4  ", atom_types.O },
    .{ "DU  :C5  ", atom_types.C_NUC },
    .{ "DU  :C6  ", atom_types.C_NUC },
});

// =============================================================================
// Public API
// =============================================================================

/// Get radius for an atom using NACCESS classification.
/// Lookup order: residue-specific -> ANY fallback -> element guess -> null
pub fn getRadius(residue: []const u8, atom: []const u8) ?f64 {
    const key = makeKeyRuntime(residue, atom);

    // Try residue-specific first
    if (residue_atoms.get(&key)) |t| {
        return t.radius;
    }

    // Fall back to ANY
    const any_key = makeKeyRuntime("ANY", atom);
    if (any_atoms.get(&any_key)) |t| {
        return t.radius;
    }

    // Fall back to element-based guessing
    return classifier.guessRadiusFromAtomName(atom);
}

/// Get class for an atom using NACCESS classification.
/// Lookup order: residue-specific -> ANY fallback -> unknown
pub fn getClass(residue: []const u8, atom: []const u8) AtomClass {
    const key = makeKeyRuntime(residue, atom);

    // Try residue-specific first
    if (residue_atoms.get(&key)) |t| {
        return t.class;
    }

    // Fall back to ANY
    const any_key = makeKeyRuntime("ANY", atom);
    if (any_atoms.get(&any_key)) |t| {
        return t.class;
    }

    return .unknown;
}

/// Get both radius and class for an atom.
/// Lookup order: residue-specific -> ANY fallback -> element guess (radius only)
pub fn getProperties(residue: []const u8, atom: []const u8) ?AtomProperties {
    const key = makeKeyRuntime(residue, atom);

    // Try residue-specific first
    if (residue_atoms.get(&key)) |t| {
        return AtomProperties{ .radius = t.radius, .class = t.class };
    }

    // Fall back to ANY
    const any_key = makeKeyRuntime("ANY", atom);
    if (any_atoms.get(&any_key)) |t| {
        return AtomProperties{ .radius = t.radius, .class = t.class };
    }

    // Fall back to element-based guessing (class will be unknown)
    if (classifier.guessRadiusFromAtomName(atom)) |radius| {
        return AtomProperties{ .radius = radius, .class = .unknown };
    }

    return null;
}

/// Classifier name
pub const name = "NACCESS";

// =============================================================================
// Tests
// =============================================================================

test "NACCESS backbone atoms (ANY)" {
    // Standard backbone atoms should work for any residue
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("ALA", "C"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("ALA", "O"));
    try std.testing.expectEqual(@as(?f64, 1.87), getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("ALA", "N"));

    // Same for different residues
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("GLY", "C"));
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("TRP", "C"));
}

test "NACCESS backbone classes" {
    try std.testing.expectEqual(AtomClass.apolar, getClass("ALA", "C"));
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "O"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("ALA", "CA"));
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "N"));
}

test "NACCESS residue-specific atoms" {
    // ARG side chain
    try std.testing.expectEqual(@as(?f64, 1.87), getRadius("ARG", "CG"));
    try std.testing.expectEqual(@as(?f64, 1.87), getRadius("ARG", "CD"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("ARG", "NE"));
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("ARG", "CZ"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("ARG", "NH1"));

    // CYS sulfur
    try std.testing.expectEqual(@as(?f64, 1.85), getRadius("CYS", "SG"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("CYS", "SG"));

    // LYS amino group
    try std.testing.expectEqual(@as(?f64, 1.50), getRadius("LYS", "NZ"));
    try std.testing.expectEqual(AtomClass.polar, getClass("LYS", "NZ"));

    // TYR hydroxyl
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("TYR", "OH"));
    try std.testing.expectEqual(AtomClass.polar, getClass("TYR", "OH"));
}

test "NACCESS aromatic residues" {
    // PHE ring carbons
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("PHE", "CG"));
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("PHE", "CD1"));
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("PHE", "CZ"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("PHE", "CG"));

    // TRP ring
    try std.testing.expectEqual(@as(?f64, 1.76), getRadius("TRP", "CG"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("TRP", "NE1"));
}

test "NACCESS selenocysteine and selenomethionine" {
    // SEC (like CYS but with SE)
    try std.testing.expectEqual(@as(?f64, 1.87), getRadius("SEC", "CB"));
    try std.testing.expectEqual(@as(?f64, 1.80), getRadius("SEC", "SE"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("SEC", "SE"));

    // MSE (like MET but with SE)
    try std.testing.expectEqual(@as(?f64, 1.87), getRadius("MSE", "CG"));
    try std.testing.expectEqual(@as(?f64, 1.80), getRadius("MSE", "SE"));
    try std.testing.expectEqual(@as(?f64, 1.87), getRadius("MSE", "CE"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("MSE", "SE"));
}

test "NACCESS nucleic acid backbone" {
    // Phosphate
    try std.testing.expectEqual(@as(?f64, 1.90), getRadius("DA", "P"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("DA", "OP1"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("DA", "OP2"));

    // Sugar
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("DA", "O5'"));
    try std.testing.expectEqual(@as(?f64, 1.80), getRadius("DA", "C5'"));
}

test "NACCESS nucleotide bases" {
    // Adenine (purine)
    try std.testing.expectEqual(@as(?f64, 1.60), getRadius("A", "N9"));
    try std.testing.expectEqual(@as(?f64, 1.80), getRadius("A", "C8"));
    try std.testing.expectEqual(@as(?f64, 1.60), getRadius("A", "N6"));

    // Guanine
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("G", "O6"));
    try std.testing.expectEqual(@as(?f64, 1.60), getRadius("G", "N2"));

    // Thymine
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("DT", "O2"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("DT", "O4"));
    try std.testing.expectEqual(@as(?f64, 1.80), getRadius("DT", "C7"));
}

test "NACCESS element fallback" {
    // Unknown atom should fall back to element-based guessing
    const radius = getRadius("ALA", "XX");
    // XX doesn't match any known atom, and X isn't a known element
    try std.testing.expectEqual(@as(?f64, null), radius);

    // Hydrogen should fall back to element guess
    const h_radius = getRadius("ALA", "H");
    try std.testing.expect(h_radius != null);
    try std.testing.expectEqual(@as(f64, 1.10), h_radius.?);
}

test "NACCESS getProperties" {
    const props = getProperties("ALA", "CA");
    try std.testing.expect(props != null);
    try std.testing.expectEqual(@as(f64, 1.87), props.?.radius);
    try std.testing.expectEqual(AtomClass.apolar, props.?.class);

    // With element fallback
    const h_props = getProperties("ALA", "H");
    try std.testing.expect(h_props != null);
    try std.testing.expectEqual(@as(f64, 1.10), h_props.?.radius);
    try std.testing.expectEqual(AtomClass.unknown, h_props.?.class);

    // Unknown
    const unknown = getProperties("XXX", "ZZZ");
    try std.testing.expectEqual(@as(?AtomProperties, null), unknown);
}

// =============================================================================
// Comprehensive Amino Acid Coverage Tests
// =============================================================================

test "NACCESS all 20 standard amino acids have backbone radii" {
    const residues = [_][]const u8{
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    };
    const backbone_atoms = [_][]const u8{ "N", "CA", "C", "O" };

    for (residues) |res| {
        for (backbone_atoms) |atom| {
            const radius = getRadius(res, atom);
            try std.testing.expect(radius != null);
            // Backbone radii should be consistent
            if (std.mem.eql(u8, atom, "N")) {
                try std.testing.expectEqual(@as(f64, 1.65), radius.?);
            } else if (std.mem.eql(u8, atom, "CA")) {
                try std.testing.expectEqual(@as(f64, 1.87), radius.?);
            } else if (std.mem.eql(u8, atom, "C")) {
                try std.testing.expectEqual(@as(f64, 1.76), radius.?);
            } else if (std.mem.eql(u8, atom, "O")) {
                try std.testing.expectEqual(@as(f64, 1.40), radius.?);
            }
        }
    }
}

test "NACCESS all amino acids with CB have CB radius" {
    // All except GLY have CB
    const residues_with_cb = [_][]const u8{
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    };

    for (residues_with_cb) |res| {
        const radius = getRadius(res, "CB");
        try std.testing.expect(radius != null);
        try std.testing.expectEqual(@as(f64, 1.87), radius.?);
    }

    // GLY has no CB - should fall back to ANY CB
    const gly_cb = getRadius("GLY", "CB");
    try std.testing.expectEqual(@as(?f64, 1.87), gly_cb); // ANY fallback
}

test "NACCESS sulfur-containing amino acids" {
    // CYS: SG
    try std.testing.expectEqual(@as(?f64, 1.85), getRadius("CYS", "SG"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("CYS", "SG"));

    // MET: SD
    try std.testing.expectEqual(@as(?f64, 1.85), getRadius("MET", "SD"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("MET", "SD"));
}

test "NACCESS all charged amino acids have correct terminal atoms" {
    // Positively charged
    try std.testing.expectEqual(@as(?f64, 1.50), getRadius("LYS", "NZ")); // amine
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("ARG", "NH1")); // guanidinium
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("ARG", "NH2"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("HIS", "ND1")); // imidazole

    // Negatively charged
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("ASP", "OD1"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("ASP", "OD2"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("GLU", "OE1"));
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("GLU", "OE2"));
}

test "NACCESS hydroxyl-containing amino acids" {
    // SER
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("SER", "OG"));
    try std.testing.expectEqual(AtomClass.polar, getClass("SER", "OG"));

    // THR
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("THR", "OG1"));
    try std.testing.expectEqual(AtomClass.polar, getClass("THR", "OG1"));

    // TYR
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("TYR", "OH"));
    try std.testing.expectEqual(AtomClass.polar, getClass("TYR", "OH"));
}

test "NACCESS amide-containing amino acids" {
    // ASN
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("ASN", "OD1"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("ASN", "ND2"));
    try std.testing.expectEqual(AtomClass.polar, getClass("ASN", "ND2"));

    // GLN
    try std.testing.expectEqual(@as(?f64, 1.40), getRadius("GLN", "OE1"));
    try std.testing.expectEqual(@as(?f64, 1.65), getRadius("GLN", "NE2"));
    try std.testing.expectEqual(AtomClass.polar, getClass("GLN", "NE2"));
}

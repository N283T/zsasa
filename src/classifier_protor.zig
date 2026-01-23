//! ProtOr classifier - Hybridization-based radii from Tsai et al. 1999
//!
//! Reference:
//! Tsai, J., Taylor, R., Chothia, C., & Gerstein, M. (1999).
//! The packing density in proteins: standard radii and volumes.
//! Journal of molecular biology, 290(1), 253-266.
//!
//! Unlike NACCESS, ProtOr does not have ANY fallback - all atoms are
//! specified per residue with hybridization-based types (C3H0, C4H1, etc.).

const std = @import("std");
const classifier = @import("classifier.zig");
const AtomClass = classifier.AtomClass;
const AtomProperties = classifier.AtomProperties;

/// Atom type with radius and polarity class
const AtomType = struct {
    radius: f64,
    class: AtomClass,
};

/// ProtOr atom types based on hybridization state
const atom_types = struct {
    // Carbon types
    const C3H0 = AtomType{ .radius = 1.61, .class = .apolar }; // sp2, no H
    const C3H1 = AtomType{ .radius = 1.76, .class = .apolar }; // sp2, 1 H
    const C4H1 = AtomType{ .radius = 1.88, .class = .apolar }; // sp3, 1 H
    const C4H2 = AtomType{ .radius = 1.88, .class = .apolar }; // sp3, 2 H
    const C4H3 = AtomType{ .radius = 1.88, .class = .apolar }; // sp3, 3 H (methyl)

    // Nitrogen types
    const N2H0 = AtomType{ .radius = 1.64, .class = .polar }; // sp, no H
    const N2H2 = AtomType{ .radius = 1.64, .class = .polar }; // (really N3H2)
    const N3H0 = AtomType{ .radius = 1.64, .class = .polar }; // sp2, no H
    const N3H1 = AtomType{ .radius = 1.64, .class = .polar }; // sp2, 1 H
    const N3H2 = AtomType{ .radius = 1.64, .class = .polar }; // sp2, 2 H (amide)
    const N4H3 = AtomType{ .radius = 1.64, .class = .polar }; // sp3, 3 H (amino)

    // Oxygen types
    const O1H0 = AtomType{ .radius = 1.42, .class = .polar }; // carbonyl
    const O2H0 = AtomType{ .radius = 1.46, .class = .polar }; // ether/ester
    const O2H1 = AtomType{ .radius = 1.46, .class = .polar }; // hydroxyl
    const O2H2 = AtomType{ .radius = 1.46, .class = .polar }; // water

    // Sulfur types
    const S2H0 = AtomType{ .radius = 1.77, .class = .polar }; // thioether
    const S2H1 = AtomType{ .radius = 1.77, .class = .polar }; // thiol

    // Selenium types
    const SE2H0 = AtomType{ .radius = 1.90, .class = .polar }; // selenoether (MSE)
    const SE2H1 = AtomType{ .radius = 1.90, .class = .polar }; // selenol (SEC)

    // Phosphorus
    const P4H0 = AtomType{ .radius = 1.80, .class = .polar };

    // Unknown polar (ASX, GLX)
    const X1H0 = AtomType{ .radius = 1.50, .class = .polar };
};

/// Residue-specific atom lookup (no ANY fallback in ProtOr)
/// Key format: "RES :ATOM" (9 chars, space-padded)
const residue_atoms = std.StaticStringMap(AtomType).initComptime(.{
    // === Standard Amino Acids ===

    // ALA
    .{ "ALA :N   ", atom_types.N3H2 },
    .{ "ALA :CA  ", atom_types.C4H1 },
    .{ "ALA :C   ", atom_types.C3H0 },
    .{ "ALA :O   ", atom_types.O1H0 },
    .{ "ALA :CB  ", atom_types.C4H3 },
    .{ "ALA :OXT ", atom_types.O2H1 },

    // ARG
    .{ "ARG :N   ", atom_types.N3H2 },
    .{ "ARG :CA  ", atom_types.C4H1 },
    .{ "ARG :C   ", atom_types.C3H0 },
    .{ "ARG :O   ", atom_types.O1H0 },
    .{ "ARG :CB  ", atom_types.C4H2 },
    .{ "ARG :CG  ", atom_types.C4H2 },
    .{ "ARG :CD  ", atom_types.C4H2 },
    .{ "ARG :NE  ", atom_types.N3H1 },
    .{ "ARG :CZ  ", atom_types.C3H0 },
    .{ "ARG :NH1 ", atom_types.N3H2 },
    .{ "ARG :NH2 ", atom_types.N3H2 },
    .{ "ARG :OXT ", atom_types.O2H1 },

    // ASN
    .{ "ASN :N   ", atom_types.N3H2 },
    .{ "ASN :CA  ", atom_types.C4H1 },
    .{ "ASN :C   ", atom_types.C3H0 },
    .{ "ASN :O   ", atom_types.O1H0 },
    .{ "ASN :CB  ", atom_types.C4H2 },
    .{ "ASN :CG  ", atom_types.C3H0 },
    .{ "ASN :OD1 ", atom_types.O1H0 },
    .{ "ASN :ND2 ", atom_types.N3H2 },
    .{ "ASN :OXT ", atom_types.O2H1 },

    // ASP
    .{ "ASP :N   ", atom_types.N3H2 },
    .{ "ASP :CA  ", atom_types.C4H1 },
    .{ "ASP :C   ", atom_types.C3H0 },
    .{ "ASP :O   ", atom_types.O1H0 },
    .{ "ASP :CB  ", atom_types.C4H2 },
    .{ "ASP :CG  ", atom_types.C3H0 },
    .{ "ASP :OD1 ", atom_types.O1H0 },
    .{ "ASP :OD2 ", atom_types.O2H1 },
    .{ "ASP :OXT ", atom_types.O2H1 },

    // CYS
    .{ "CYS :N   ", atom_types.N3H2 },
    .{ "CYS :CA  ", atom_types.C4H1 },
    .{ "CYS :C   ", atom_types.C3H0 },
    .{ "CYS :O   ", atom_types.O1H0 },
    .{ "CYS :CB  ", atom_types.C4H2 },
    .{ "CYS :SG  ", atom_types.S2H1 },
    .{ "CYS :OXT ", atom_types.O2H1 },

    // GLN
    .{ "GLN :N   ", atom_types.N3H2 },
    .{ "GLN :CA  ", atom_types.C4H1 },
    .{ "GLN :C   ", atom_types.C3H0 },
    .{ "GLN :O   ", atom_types.O1H0 },
    .{ "GLN :CB  ", atom_types.C4H2 },
    .{ "GLN :CG  ", atom_types.C4H2 },
    .{ "GLN :CD  ", atom_types.C3H0 },
    .{ "GLN :OE1 ", atom_types.O1H0 },
    .{ "GLN :NE2 ", atom_types.N3H2 },
    .{ "GLN :OXT ", atom_types.O2H1 },

    // GLU
    .{ "GLU :N   ", atom_types.N3H2 },
    .{ "GLU :CA  ", atom_types.C4H1 },
    .{ "GLU :C   ", atom_types.C3H0 },
    .{ "GLU :O   ", atom_types.O1H0 },
    .{ "GLU :CB  ", atom_types.C4H2 },
    .{ "GLU :CG  ", atom_types.C4H2 },
    .{ "GLU :CD  ", atom_types.C3H0 },
    .{ "GLU :OE1 ", atom_types.O1H0 },
    .{ "GLU :OE2 ", atom_types.O2H1 },
    .{ "GLU :OXT ", atom_types.O2H1 },

    // GLY
    .{ "GLY :N   ", atom_types.N3H2 },
    .{ "GLY :CA  ", atom_types.C4H2 },
    .{ "GLY :C   ", atom_types.C3H0 },
    .{ "GLY :O   ", atom_types.O1H0 },
    .{ "GLY :OXT ", atom_types.O2H1 },

    // HIS
    .{ "HIS :N   ", atom_types.N3H2 },
    .{ "HIS :CA  ", atom_types.C4H1 },
    .{ "HIS :C   ", atom_types.C3H0 },
    .{ "HIS :O   ", atom_types.O1H0 },
    .{ "HIS :CB  ", atom_types.C4H2 },
    .{ "HIS :CG  ", atom_types.C3H0 },
    .{ "HIS :ND1 ", atom_types.N3H1 },
    .{ "HIS :CD2 ", atom_types.C3H1 },
    .{ "HIS :CE1 ", atom_types.C3H1 },
    .{ "HIS :NE2 ", atom_types.N3H1 },
    .{ "HIS :OXT ", atom_types.O2H1 },

    // ILE
    .{ "ILE :N   ", atom_types.N3H2 },
    .{ "ILE :CA  ", atom_types.C4H1 },
    .{ "ILE :C   ", atom_types.C3H0 },
    .{ "ILE :O   ", atom_types.O1H0 },
    .{ "ILE :CB  ", atom_types.C4H1 },
    .{ "ILE :CG1 ", atom_types.C4H2 },
    .{ "ILE :CG2 ", atom_types.C4H3 },
    .{ "ILE :CD1 ", atom_types.C4H3 },
    .{ "ILE :OXT ", atom_types.O2H1 },

    // LEU
    .{ "LEU :N   ", atom_types.N3H2 },
    .{ "LEU :CA  ", atom_types.C4H1 },
    .{ "LEU :C   ", atom_types.C3H0 },
    .{ "LEU :O   ", atom_types.O1H0 },
    .{ "LEU :CB  ", atom_types.C4H2 },
    .{ "LEU :CG  ", atom_types.C4H1 },
    .{ "LEU :CD1 ", atom_types.C4H3 },
    .{ "LEU :CD2 ", atom_types.C4H3 },
    .{ "LEU :OXT ", atom_types.O2H1 },

    // LYS
    .{ "LYS :N   ", atom_types.N3H2 },
    .{ "LYS :CA  ", atom_types.C4H1 },
    .{ "LYS :C   ", atom_types.C3H0 },
    .{ "LYS :O   ", atom_types.O1H0 },
    .{ "LYS :CB  ", atom_types.C4H2 },
    .{ "LYS :CG  ", atom_types.C4H2 },
    .{ "LYS :CD  ", atom_types.C4H2 },
    .{ "LYS :CE  ", atom_types.C4H2 },
    .{ "LYS :NZ  ", atom_types.N4H3 },
    .{ "LYS :OXT ", atom_types.O2H1 },

    // MET
    .{ "MET :N   ", atom_types.N3H2 },
    .{ "MET :CA  ", atom_types.C4H1 },
    .{ "MET :C   ", atom_types.C3H0 },
    .{ "MET :O   ", atom_types.O1H0 },
    .{ "MET :CB  ", atom_types.C4H2 },
    .{ "MET :CG  ", atom_types.C4H2 },
    .{ "MET :SD  ", atom_types.S2H0 },
    .{ "MET :CE  ", atom_types.C4H3 },
    .{ "MET :OXT ", atom_types.O2H1 },

    // PHE
    .{ "PHE :N   ", atom_types.N3H2 },
    .{ "PHE :CA  ", atom_types.C4H1 },
    .{ "PHE :C   ", atom_types.C3H0 },
    .{ "PHE :O   ", atom_types.O1H0 },
    .{ "PHE :CB  ", atom_types.C4H2 },
    .{ "PHE :CG  ", atom_types.C3H0 },
    .{ "PHE :CD1 ", atom_types.C3H1 },
    .{ "PHE :CD2 ", atom_types.C3H1 },
    .{ "PHE :CE1 ", atom_types.C3H1 },
    .{ "PHE :CE2 ", atom_types.C3H1 },
    .{ "PHE :CZ  ", atom_types.C3H1 },
    .{ "PHE :OXT ", atom_types.O2H1 },

    // PRO
    .{ "PRO :N   ", atom_types.N3H1 },
    .{ "PRO :CA  ", atom_types.C4H1 },
    .{ "PRO :C   ", atom_types.C3H0 },
    .{ "PRO :O   ", atom_types.O1H0 },
    .{ "PRO :CB  ", atom_types.C4H2 },
    .{ "PRO :CG  ", atom_types.C4H2 },
    .{ "PRO :CD  ", atom_types.C4H2 },
    .{ "PRO :OXT ", atom_types.O2H1 },

    // SER
    .{ "SER :N   ", atom_types.N3H2 },
    .{ "SER :CA  ", atom_types.C4H1 },
    .{ "SER :C   ", atom_types.C3H0 },
    .{ "SER :O   ", atom_types.O1H0 },
    .{ "SER :CB  ", atom_types.C4H2 },
    .{ "SER :OG  ", atom_types.O2H1 },
    .{ "SER :OXT ", atom_types.O2H1 },

    // THR
    .{ "THR :N   ", atom_types.N3H2 },
    .{ "THR :CA  ", atom_types.C4H1 },
    .{ "THR :C   ", atom_types.C3H0 },
    .{ "THR :O   ", atom_types.O1H0 },
    .{ "THR :CB  ", atom_types.C4H1 },
    .{ "THR :OG1 ", atom_types.O2H1 },
    .{ "THR :CG2 ", atom_types.C4H3 },
    .{ "THR :OXT ", atom_types.O2H1 },

    // TRP
    .{ "TRP :N   ", atom_types.N3H2 },
    .{ "TRP :CA  ", atom_types.C4H1 },
    .{ "TRP :C   ", atom_types.C3H0 },
    .{ "TRP :O   ", atom_types.O1H0 },
    .{ "TRP :CB  ", atom_types.C4H2 },
    .{ "TRP :CG  ", atom_types.C3H0 },
    .{ "TRP :CD1 ", atom_types.C3H1 },
    .{ "TRP :CD2 ", atom_types.C3H0 },
    .{ "TRP :NE1 ", atom_types.N3H1 },
    .{ "TRP :CE2 ", atom_types.C3H0 },
    .{ "TRP :CE3 ", atom_types.C3H1 },
    .{ "TRP :CZ2 ", atom_types.C3H1 },
    .{ "TRP :CZ3 ", atom_types.C3H1 },
    .{ "TRP :CH2 ", atom_types.C3H1 },
    .{ "TRP :OXT ", atom_types.O2H1 },

    // TYR
    .{ "TYR :N   ", atom_types.N3H2 },
    .{ "TYR :CA  ", atom_types.C4H1 },
    .{ "TYR :C   ", atom_types.C3H0 },
    .{ "TYR :O   ", atom_types.O1H0 },
    .{ "TYR :CB  ", atom_types.C4H2 },
    .{ "TYR :CG  ", atom_types.C3H0 },
    .{ "TYR :CD1 ", atom_types.C3H1 },
    .{ "TYR :CD2 ", atom_types.C3H1 },
    .{ "TYR :CE1 ", atom_types.C3H1 },
    .{ "TYR :CE2 ", atom_types.C3H1 },
    .{ "TYR :CZ  ", atom_types.C3H0 },
    .{ "TYR :OH  ", atom_types.O2H1 },
    .{ "TYR :OXT ", atom_types.O2H1 },

    // VAL
    .{ "VAL :N   ", atom_types.N3H2 },
    .{ "VAL :CA  ", atom_types.C4H1 },
    .{ "VAL :C   ", atom_types.C3H0 },
    .{ "VAL :O   ", atom_types.O1H0 },
    .{ "VAL :CB  ", atom_types.C4H1 },
    .{ "VAL :CG1 ", atom_types.C4H3 },
    .{ "VAL :CG2 ", atom_types.C4H3 },
    .{ "VAL :OXT ", atom_types.O2H1 },

    // === Non-standard Amino Acids ===

    // ASX (ambiguous ASN/ASP)
    .{ "ASX :N   ", atom_types.N3H2 },
    .{ "ASX :CA  ", atom_types.C4H1 },
    .{ "ASX :C   ", atom_types.C3H0 },
    .{ "ASX :O   ", atom_types.O1H0 },
    .{ "ASX :CB  ", atom_types.C4H2 },
    .{ "ASX :CG  ", atom_types.C3H0 },
    .{ "ASX :XD1 ", atom_types.X1H0 },
    .{ "ASX :XD2 ", atom_types.X1H0 },
    .{ "ASX :OXT ", atom_types.O2H1 },

    // GLX (ambiguous GLN/GLU)
    .{ "GLX :N   ", atom_types.N3H2 },
    .{ "GLX :CA  ", atom_types.C4H1 },
    .{ "GLX :C   ", atom_types.C3H0 },
    .{ "GLX :O   ", atom_types.O1H0 },
    .{ "GLX :CB  ", atom_types.C4H2 },
    .{ "GLX :CG  ", atom_types.C4H2 },
    .{ "GLX :CD  ", atom_types.C3H0 },
    .{ "GLX :XE1 ", atom_types.X1H0 },
    .{ "GLX :XE2 ", atom_types.X1H0 },
    .{ "GLX :OXT ", atom_types.O2H1 },

    // SEC (selenocysteine)
    .{ "SEC :N   ", atom_types.N3H2 },
    .{ "SEC :CA  ", atom_types.C4H1 },
    .{ "SEC :CB  ", atom_types.C4H2 },
    .{ "SEC :SE  ", atom_types.SE2H1 },
    .{ "SEC :C   ", atom_types.C3H0 },
    .{ "SEC :O   ", atom_types.O1H0 },
    .{ "SEC :OXT ", atom_types.O2H1 },

    // MSE (selenomethionine)
    .{ "MSE :N   ", atom_types.N3H2 },
    .{ "MSE :CA  ", atom_types.C4H1 },
    .{ "MSE :C   ", atom_types.C3H0 },
    .{ "MSE :O   ", atom_types.O1H0 },
    .{ "MSE :OXT ", atom_types.O2H1 },
    .{ "MSE :CB  ", atom_types.C4H2 },
    .{ "MSE :CG  ", atom_types.C4H2 },
    .{ "MSE :SE  ", atom_types.SE2H0 },
    .{ "MSE :CE  ", atom_types.C4H3 },

    // PYL (pyrrolysine)
    .{ "PYL :N   ", atom_types.N3H2 },
    .{ "PYL :CA  ", atom_types.C4H1 },
    .{ "PYL :C   ", atom_types.C3H0 },
    .{ "PYL :O   ", atom_types.O1H0 },
    .{ "PYL :OXT ", atom_types.O2H1 },
    .{ "PYL :CB  ", atom_types.C4H2 },
    .{ "PYL :CG  ", atom_types.C4H2 },
    .{ "PYL :CD  ", atom_types.C4H2 },
    .{ "PYL :CE  ", atom_types.C4H2 },
    .{ "PYL :NZ  ", atom_types.N3H1 },
    .{ "PYL :C2  ", atom_types.C3H0 },
    .{ "PYL :O2  ", atom_types.O1H0 },
    .{ "PYL :CA2 ", atom_types.C4H1 },
    .{ "PYL :N2  ", atom_types.N2H0 },
    .{ "PYL :CB2 ", atom_types.C4H3 },
    .{ "PYL :CG2 ", atom_types.C4H1 },
    .{ "PYL :CD2 ", atom_types.C4H2 },
    .{ "PYL :CE2 ", atom_types.C3H1 },

    // === Capping Groups ===

    // ACE (acetyl)
    .{ "ACE :C   ", atom_types.C3H1 },
    .{ "ACE :O   ", atom_types.O1H0 },
    .{ "ACE :CH3 ", atom_types.C4H3 },

    // NH2 (amide)
    .{ "NH2 :N   ", atom_types.N2H2 },

    // HOH (water)
    .{ "HOH :O   ", atom_types.O2H2 },

    // === RNA Nucleotides ===

    // A (adenine)
    .{ "A   :OP3 ", atom_types.O2H1 },
    .{ "A   :P   ", atom_types.P4H0 },
    .{ "A   :OP1 ", atom_types.O1H0 },
    .{ "A   :OP2 ", atom_types.O2H1 },
    .{ "A   :O5' ", atom_types.O2H0 },
    .{ "A   :C5' ", atom_types.C4H2 },
    .{ "A   :C4' ", atom_types.C4H1 },
    .{ "A   :O4' ", atom_types.O2H0 },
    .{ "A   :C3' ", atom_types.C4H1 },
    .{ "A   :O3' ", atom_types.O2H1 },
    .{ "A   :C2' ", atom_types.C4H1 },
    .{ "A   :O2' ", atom_types.O2H1 },
    .{ "A   :C1' ", atom_types.C4H1 },
    .{ "A   :N9  ", atom_types.N3H0 },
    .{ "A   :C8  ", atom_types.C3H1 },
    .{ "A   :N7  ", atom_types.N2H0 },
    .{ "A   :C5  ", atom_types.C3H0 },
    .{ "A   :C6  ", atom_types.C3H0 },
    .{ "A   :N6  ", atom_types.N3H2 },
    .{ "A   :N1  ", atom_types.N2H0 },
    .{ "A   :C2  ", atom_types.C3H1 },
    .{ "A   :N3  ", atom_types.N2H0 },
    .{ "A   :C4  ", atom_types.C3H0 },

    // C (cytosine)
    .{ "C   :OP3 ", atom_types.O2H1 },
    .{ "C   :P   ", atom_types.P4H0 },
    .{ "C   :OP1 ", atom_types.O1H0 },
    .{ "C   :OP2 ", atom_types.O2H1 },
    .{ "C   :O5' ", atom_types.O2H0 },
    .{ "C   :C5' ", atom_types.C4H2 },
    .{ "C   :C4' ", atom_types.C4H1 },
    .{ "C   :O4' ", atom_types.O2H0 },
    .{ "C   :C3' ", atom_types.C4H1 },
    .{ "C   :O3' ", atom_types.O2H1 },
    .{ "C   :C2' ", atom_types.C4H1 },
    .{ "C   :O2' ", atom_types.O2H1 },
    .{ "C   :C1' ", atom_types.C4H1 },
    .{ "C   :N1  ", atom_types.N3H0 },
    .{ "C   :C2  ", atom_types.C3H0 },
    .{ "C   :O2  ", atom_types.O1H0 },
    .{ "C   :N3  ", atom_types.N2H0 },
    .{ "C   :C4  ", atom_types.C3H0 },
    .{ "C   :N4  ", atom_types.N3H2 },
    .{ "C   :C5  ", atom_types.C3H1 },
    .{ "C   :C6  ", atom_types.C3H1 },

    // G (guanine)
    .{ "G   :OP3 ", atom_types.O2H1 },
    .{ "G   :P   ", atom_types.P4H0 },
    .{ "G   :OP1 ", atom_types.O1H0 },
    .{ "G   :OP2 ", atom_types.O2H1 },
    .{ "G   :O5' ", atom_types.O2H0 },
    .{ "G   :C5' ", atom_types.C4H2 },
    .{ "G   :C4' ", atom_types.C4H1 },
    .{ "G   :O4' ", atom_types.O2H0 },
    .{ "G   :C3' ", atom_types.C4H1 },
    .{ "G   :O3' ", atom_types.O2H1 },
    .{ "G   :C2' ", atom_types.C4H1 },
    .{ "G   :O2' ", atom_types.O2H1 },
    .{ "G   :C1' ", atom_types.C4H1 },
    .{ "G   :N9  ", atom_types.N3H0 },
    .{ "G   :C8  ", atom_types.C3H1 },
    .{ "G   :N7  ", atom_types.N2H0 },
    .{ "G   :C5  ", atom_types.C3H0 },
    .{ "G   :C6  ", atom_types.C3H0 },
    .{ "G   :O6  ", atom_types.O1H0 },
    .{ "G   :N1  ", atom_types.N3H1 },
    .{ "G   :C2  ", atom_types.C3H0 },
    .{ "G   :N2  ", atom_types.N3H2 },
    .{ "G   :N3  ", atom_types.N2H0 },
    .{ "G   :C4  ", atom_types.C3H0 },

    // I (inosine)
    .{ "I   :OP3 ", atom_types.O2H1 },
    .{ "I   :P   ", atom_types.P4H0 },
    .{ "I   :OP1 ", atom_types.O1H0 },
    .{ "I   :OP2 ", atom_types.O2H1 },
    .{ "I   :O5' ", atom_types.O2H0 },
    .{ "I   :C5' ", atom_types.C4H2 },
    .{ "I   :C4' ", atom_types.C4H1 },
    .{ "I   :O4' ", atom_types.O2H0 },
    .{ "I   :C3' ", atom_types.C4H1 },
    .{ "I   :O3' ", atom_types.O2H1 },
    .{ "I   :C2' ", atom_types.C4H1 },
    .{ "I   :O2' ", atom_types.O2H1 },
    .{ "I   :C1' ", atom_types.C4H1 },
    .{ "I   :N9  ", atom_types.N3H0 },
    .{ "I   :C8  ", atom_types.C3H1 },
    .{ "I   :N7  ", atom_types.N2H0 },
    .{ "I   :C5  ", atom_types.C3H0 },
    .{ "I   :C6  ", atom_types.C3H0 },
    .{ "I   :O6  ", atom_types.O1H0 },
    .{ "I   :N1  ", atom_types.N3H1 },
    .{ "I   :C2  ", atom_types.C3H1 },
    .{ "I   :N3  ", atom_types.N2H0 },
    .{ "I   :C4  ", atom_types.C3H0 },

    // T (thymine - for RNA)
    .{ "T   :OP3 ", atom_types.O2H1 },
    .{ "T   :P   ", atom_types.P4H0 },
    .{ "T   :OP1 ", atom_types.O1H0 },
    .{ "T   :OP2 ", atom_types.O2H1 },
    .{ "T   :O5' ", atom_types.O2H0 },
    .{ "T   :C5' ", atom_types.C4H2 },
    .{ "T   :C4' ", atom_types.C4H1 },
    .{ "T   :O4' ", atom_types.O2H0 },
    .{ "T   :C3' ", atom_types.C4H1 },
    .{ "T   :O3' ", atom_types.O2H1 },
    .{ "T   :C2' ", atom_types.C4H2 },
    .{ "T   :C1' ", atom_types.C4H1 },
    .{ "T   :N1  ", atom_types.N3H0 },
    .{ "T   :C2  ", atom_types.C3H0 },
    .{ "T   :O2  ", atom_types.O1H0 },
    .{ "T   :N3  ", atom_types.N3H1 },
    .{ "T   :C4  ", atom_types.C3H0 },
    .{ "T   :O4  ", atom_types.O1H0 },
    .{ "T   :C5  ", atom_types.C3H0 },
    .{ "T   :C7  ", atom_types.C4H3 },
    .{ "T   :C6  ", atom_types.C3H1 },

    // U (uracil)
    .{ "U   :OP3 ", atom_types.O2H1 },
    .{ "U   :P   ", atom_types.P4H0 },
    .{ "U   :OP1 ", atom_types.O1H0 },
    .{ "U   :OP2 ", atom_types.O2H1 },
    .{ "U   :O5' ", atom_types.O2H0 },
    .{ "U   :C5' ", atom_types.C4H2 },
    .{ "U   :C4' ", atom_types.C4H1 },
    .{ "U   :O4' ", atom_types.O2H0 },
    .{ "U   :C3' ", atom_types.C4H1 },
    .{ "U   :O3' ", atom_types.O2H1 },
    .{ "U   :C2' ", atom_types.C4H1 },
    .{ "U   :O2' ", atom_types.O2H1 },
    .{ "U   :C1' ", atom_types.C4H1 },
    .{ "U   :N1  ", atom_types.N3H0 },
    .{ "U   :C2  ", atom_types.C3H0 },
    .{ "U   :O2  ", atom_types.O1H0 },
    .{ "U   :N3  ", atom_types.N3H1 },
    .{ "U   :C4  ", atom_types.C3H0 },
    .{ "U   :O4  ", atom_types.O1H0 },
    .{ "U   :C5  ", atom_types.C3H1 },
    .{ "U   :C6  ", atom_types.C3H1 },

    // === DNA Nucleotides ===

    // DA (deoxyadenosine)
    .{ "DA  :OP3 ", atom_types.O2H1 },
    .{ "DA  :P   ", atom_types.P4H0 },
    .{ "DA  :OP1 ", atom_types.O1H0 },
    .{ "DA  :OP2 ", atom_types.O2H1 },
    .{ "DA  :O5' ", atom_types.O2H0 },
    .{ "DA  :C5' ", atom_types.C4H2 },
    .{ "DA  :C4' ", atom_types.C4H1 },
    .{ "DA  :O4' ", atom_types.O2H0 },
    .{ "DA  :C3' ", atom_types.C4H1 },
    .{ "DA  :O3' ", atom_types.O2H1 },
    .{ "DA  :C2' ", atom_types.C4H2 },
    .{ "DA  :C1' ", atom_types.C4H1 },
    .{ "DA  :N9  ", atom_types.N3H0 },
    .{ "DA  :C8  ", atom_types.C3H1 },
    .{ "DA  :N7  ", atom_types.N2H0 },
    .{ "DA  :C5  ", atom_types.C3H0 },
    .{ "DA  :C6  ", atom_types.C3H0 },
    .{ "DA  :N6  ", atom_types.N3H2 },
    .{ "DA  :N1  ", atom_types.N2H0 },
    .{ "DA  :C2  ", atom_types.C3H1 },
    .{ "DA  :N3  ", atom_types.N2H0 },
    .{ "DA  :C4  ", atom_types.C3H0 },

    // DC (deoxycytidine)
    .{ "DC  :OP3 ", atom_types.O2H1 },
    .{ "DC  :P   ", atom_types.P4H0 },
    .{ "DC  :OP1 ", atom_types.O1H0 },
    .{ "DC  :OP2 ", atom_types.O2H1 },
    .{ "DC  :O5' ", atom_types.O2H0 },
    .{ "DC  :C5' ", atom_types.C4H2 },
    .{ "DC  :C4' ", atom_types.C4H1 },
    .{ "DC  :O4' ", atom_types.O2H0 },
    .{ "DC  :C3' ", atom_types.C4H1 },
    .{ "DC  :O3' ", atom_types.O2H1 },
    .{ "DC  :C2' ", atom_types.C4H2 },
    .{ "DC  :C1' ", atom_types.C4H1 },
    .{ "DC  :N1  ", atom_types.N3H0 },
    .{ "DC  :C2  ", atom_types.C3H0 },
    .{ "DC  :O2  ", atom_types.O1H0 },
    .{ "DC  :N3  ", atom_types.N2H0 },
    .{ "DC  :C4  ", atom_types.C3H0 },
    .{ "DC  :N4  ", atom_types.N3H2 },
    .{ "DC  :C5  ", atom_types.C3H1 },
    .{ "DC  :C6  ", atom_types.C3H1 },

    // DG (deoxyguanosine)
    .{ "DG  :OP3 ", atom_types.O2H1 },
    .{ "DG  :P   ", atom_types.P4H0 },
    .{ "DG  :OP1 ", atom_types.O1H0 },
    .{ "DG  :OP2 ", atom_types.O2H1 },
    .{ "DG  :O5' ", atom_types.O2H0 },
    .{ "DG  :C5' ", atom_types.C4H2 },
    .{ "DG  :C4' ", atom_types.C4H1 },
    .{ "DG  :O4' ", atom_types.O2H0 },
    .{ "DG  :C3' ", atom_types.C4H1 },
    .{ "DG  :O3' ", atom_types.O2H1 },
    .{ "DG  :C2' ", atom_types.C4H2 },
    .{ "DG  :C1' ", atom_types.C4H1 },
    .{ "DG  :N9  ", atom_types.N3H0 },
    .{ "DG  :C8  ", atom_types.C3H1 },
    .{ "DG  :N7  ", atom_types.N2H0 },
    .{ "DG  :C5  ", atom_types.C3H0 },
    .{ "DG  :C6  ", atom_types.C3H0 },
    .{ "DG  :O6  ", atom_types.O1H0 },
    .{ "DG  :N1  ", atom_types.N3H1 },
    .{ "DG  :C2  ", atom_types.C3H0 },
    .{ "DG  :N2  ", atom_types.N3H2 },
    .{ "DG  :N3  ", atom_types.N2H0 },
    .{ "DG  :C4  ", atom_types.C3H0 },

    // DI (deoxyinosine)
    .{ "DI  :OP3 ", atom_types.O2H1 },
    .{ "DI  :P   ", atom_types.P4H0 },
    .{ "DI  :OP1 ", atom_types.O1H0 },
    .{ "DI  :OP2 ", atom_types.O2H1 },
    .{ "DI  :O5' ", atom_types.O2H0 },
    .{ "DI  :C5' ", atom_types.C4H2 },
    .{ "DI  :C4' ", atom_types.C4H1 },
    .{ "DI  :O4' ", atom_types.O2H0 },
    .{ "DI  :C3' ", atom_types.C4H1 },
    .{ "DI  :O3' ", atom_types.O2H1 },
    .{ "DI  :C2' ", atom_types.C4H2 },
    .{ "DI  :C1' ", atom_types.C4H1 },
    .{ "DI  :N9  ", atom_types.N3H0 },
    .{ "DI  :C8  ", atom_types.C3H1 },
    .{ "DI  :N7  ", atom_types.N2H0 },
    .{ "DI  :C5  ", atom_types.C3H0 },
    .{ "DI  :C6  ", atom_types.C3H0 },
    .{ "DI  :O6  ", atom_types.O1H0 },
    .{ "DI  :N1  ", atom_types.N3H1 },
    .{ "DI  :C2  ", atom_types.C3H1 },
    .{ "DI  :N3  ", atom_types.N2H0 },
    .{ "DI  :C4  ", atom_types.C3H0 },

    // DT (deoxythymidine)
    .{ "DT  :OP3 ", atom_types.O2H1 },
    .{ "DT  :P   ", atom_types.P4H0 },
    .{ "DT  :OP1 ", atom_types.O1H0 },
    .{ "DT  :OP2 ", atom_types.O2H1 },
    .{ "DT  :O5' ", atom_types.O2H0 },
    .{ "DT  :C5' ", atom_types.C4H2 },
    .{ "DT  :C4' ", atom_types.C4H1 },
    .{ "DT  :O4' ", atom_types.O2H0 },
    .{ "DT  :C3' ", atom_types.C4H1 },
    .{ "DT  :O3' ", atom_types.O2H1 },
    .{ "DT  :C2' ", atom_types.C4H2 },
    .{ "DT  :C1' ", atom_types.C4H1 },
    .{ "DT  :N1  ", atom_types.N3H0 },
    .{ "DT  :C2  ", atom_types.C3H0 },
    .{ "DT  :O2  ", atom_types.O1H0 },
    .{ "DT  :N3  ", atom_types.N3H1 },
    .{ "DT  :C4  ", atom_types.C3H0 },
    .{ "DT  :O4  ", atom_types.O1H0 },
    .{ "DT  :C5  ", atom_types.C3H0 },
    .{ "DT  :C7  ", atom_types.C4H3 },
    .{ "DT  :C6  ", atom_types.C3H1 },

    // DU (deoxyuridine)
    .{ "DU  :OP3 ", atom_types.O2H1 },
    .{ "DU  :P   ", atom_types.P4H0 },
    .{ "DU  :OP1 ", atom_types.O1H0 },
    .{ "DU  :OP2 ", atom_types.O2H1 },
    .{ "DU  :O5' ", atom_types.O2H0 },
    .{ "DU  :C5' ", atom_types.C4H2 },
    .{ "DU  :C4' ", atom_types.C4H1 },
    .{ "DU  :O4' ", atom_types.O2H0 },
    .{ "DU  :C3' ", atom_types.C4H1 },
    .{ "DU  :O3' ", atom_types.O2H1 },
    .{ "DU  :C2' ", atom_types.C4H2 },
    .{ "DU  :C1' ", atom_types.C4H1 },
    .{ "DU  :N1  ", atom_types.N3H0 },
    .{ "DU  :C2  ", atom_types.C3H0 },
    .{ "DU  :O2  ", atom_types.O1H0 },
    .{ "DU  :N3  ", atom_types.N3H1 },
    .{ "DU  :C4  ", atom_types.C3H0 },
    .{ "DU  :O4  ", atom_types.O1H0 },
    .{ "DU  :C5  ", atom_types.C3H1 },
    .{ "DU  :C6  ", atom_types.C3H1 },
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
/// ProtOr does not have ANY fallback, only element-based guessing
pub fn getRadius(residue: []const u8, atom: []const u8) ?f64 {
    const key = makeKeyRuntime(residue, atom);
    if (residue_atoms.get(&key)) |t| {
        return t.radius;
    }
    // No ANY fallback in ProtOr, go directly to element guess
    return classifier.guessRadiusFromAtomName(atom);
}

/// Get polarity class for a (residue, atom) pair
pub fn getClass(residue: []const u8, atom: []const u8) AtomClass {
    const key = makeKeyRuntime(residue, atom);
    if (residue_atoms.get(&key)) |t| {
        return t.class;
    }
    return .unknown;
}

/// Get both radius and class
pub fn getProperties(residue: []const u8, atom: []const u8) ?AtomProperties {
    const key = makeKeyRuntime(residue, atom);
    if (residue_atoms.get(&key)) |t| {
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

test "ProtOr backbone atoms" {
    // ALA backbone
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), getRadius("ALA", "N").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("ALA", "CA").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.61), getRadius("ALA", "C").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.42), getRadius("ALA", "O").?, 0.001);

    // GLY has different CA (C4H2 instead of C4H1)
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("GLY", "CA").?, 0.001);
}

test "ProtOr side chain atoms" {
    // ALA CB (methyl)
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("ALA", "CB").?, 0.001);

    // PHE aromatic ring
    try std.testing.expectApproxEqAbs(@as(f64, 1.61), getRadius("PHE", "CG").?, 0.001); // C3H0
    try std.testing.expectApproxEqAbs(@as(f64, 1.76), getRadius("PHE", "CD1").?, 0.001); // C3H1

    // CYS thiol
    try std.testing.expectApproxEqAbs(@as(f64, 1.77), getRadius("CYS", "SG").?, 0.001);

    // MET thioether
    try std.testing.expectApproxEqAbs(@as(f64, 1.77), getRadius("MET", "SD").?, 0.001);

    // LYS amino
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), getRadius("LYS", "NZ").?, 0.001);
}

test "ProtOr polarity classes" {
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "N"));
    try std.testing.expectEqual(AtomClass.polar, getClass("ALA", "O"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("ALA", "CA"));
    try std.testing.expectEqual(AtomClass.apolar, getClass("ALA", "CB"));
    try std.testing.expectEqual(AtomClass.polar, getClass("CYS", "SG")); // S is polar in ProtOr
}

test "ProtOr selenocysteine and selenomethionine" {
    // SEC
    try std.testing.expectApproxEqAbs(@as(f64, 1.90), getRadius("SEC", "SE").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("SEC", "CB").?, 0.001);

    // MSE
    try std.testing.expectApproxEqAbs(@as(f64, 1.90), getRadius("MSE", "SE").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("MSE", "CG").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("MSE", "CE").?, 0.001);
}

test "ProtOr nucleotides" {
    // RNA adenine
    try std.testing.expectApproxEqAbs(@as(f64, 1.80), getRadius("A", "P").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("A", "C1'").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), getRadius("A", "N9").?, 0.001);

    // DNA thymine
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), getRadius("DT", "C7").?, 0.001); // methyl
}

test "ProtOr element-based fallback" {
    // Unknown residue - should use element guess
    // Note: " CA " (with leading space) -> Carbon (1.70)
    // "CA" (without leading space) -> Calcium (2.31)
    const r = getRadius("UNK", " CA ");
    try std.testing.expect(r != null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.70), r.?, 0.001); // C element

    // Without leading space, "CA" is Calcium
    const r2 = getRadius("UNK", "CA");
    try std.testing.expect(r2 != null);
    try std.testing.expectApproxEqAbs(@as(f64, 2.31), r2.?, 0.001); // CA element (Calcium)
}

test "ProtOr getProperties" {
    const props = getProperties("ALA", "CB");
    try std.testing.expect(props != null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), props.?.radius, 0.001);
    try std.testing.expectEqual(AtomClass.apolar, props.?.class);
}

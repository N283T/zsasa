//! CCD-based classifier with hardcoded StaticStringMap for standard residues
//! and runtime component support for non-standard residues.
//!
//! This classifier is the primary lookup for CCD-based SASA calculation.
//! Standard residues use a compile-time StaticStringMap for O(1) lookup.
//! Non-standard residues are added at runtime via `addComponent`.
//!
//! The hardcoded table produces the same radii as ProtOr for standard residues.
//! Reference: Tsai et al. 1999, J. Mol. Biol. 290:253-266

const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const AtomClass = classifier.AtomClass;
const AtomProperties = classifier.AtomProperties;
const hybridization = @import("hybridization.zig");
const Component = hybridization.Component;

// =============================================================================
// Hardcoded table (same entries and values as classifier_protor.zig)
// =============================================================================

/// Key format: "RES :ATOM" — 4 chars residue (space-padded) + ':' + 4 chars atom (space-padded) = 9 chars total
const hardcoded_table = std.StaticStringMap(AtomProperties).initComptime(.{
    // === Standard Amino Acids ===

    // ALA
    .{ "ALA :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ALA :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ALA :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ALA :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ALA :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ALA :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // ARG
    .{ "ARG :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ARG :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ARG :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ARG :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ARG :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ARG :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ARG :CD  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ARG :NE  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ARG :CZ  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ARG :NH1 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ARG :NH2 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ARG :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // ASN
    .{ "ASN :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ASN :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ASN :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ASN :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ASN :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ASN :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ASN :OD1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ASN :ND2 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ASN :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // ASP
    .{ "ASP :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ASP :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ASP :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ASP :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ASP :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ASP :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ASP :OD1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ASP :OD2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "ASP :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // CYS
    .{ "CYS :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "CYS :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "CYS :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "CYS :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "CYS :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "CYS :SG  ", AtomProperties{ .radius = 1.77, .class = .polar } },
    .{ "CYS :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // GLN
    .{ "GLN :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "GLN :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLN :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLN :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "GLN :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLN :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLN :CD  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLN :OE1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "GLN :NE2 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "GLN :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // GLU
    .{ "GLU :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "GLU :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLU :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLU :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "GLU :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLU :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLU :CD  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLU :OE1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "GLU :OE2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "GLU :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // GLY
    .{ "GLY :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "GLY :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLY :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLY :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "GLY :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // HIS
    .{ "HIS :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "HIS :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "HIS :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "HIS :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "HIS :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "HIS :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "HIS :ND1 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "HIS :CD2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "HIS :CE1 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "HIS :NE2 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "HIS :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // ILE
    .{ "ILE :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ILE :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ILE :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ILE :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ILE :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ILE :CG1 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ILE :CG2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ILE :CD1 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ILE :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // LEU
    .{ "LEU :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "LEU :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LEU :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "LEU :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "LEU :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LEU :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LEU :CD1 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LEU :CD2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LEU :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // LYS
    .{ "LYS :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "LYS :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LYS :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "LYS :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "LYS :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LYS :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LYS :CD  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LYS :CE  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "LYS :NZ  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "LYS :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // MET
    .{ "MET :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "MET :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MET :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "MET :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "MET :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MET :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MET :SD  ", AtomProperties{ .radius = 1.77, .class = .polar } },
    .{ "MET :CE  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MET :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // PHE
    .{ "PHE :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PHE :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PHE :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PHE :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PHE :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PHE :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PHE :CD1 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "PHE :CD2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "PHE :CE1 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "PHE :CE2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "PHE :CZ  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "PHE :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // PRO
    .{ "PRO :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PRO :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PRO :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PRO :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PRO :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PRO :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PRO :CD  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PRO :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // SER
    .{ "SER :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "SER :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "SER :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "SER :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "SER :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "SER :OG  ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "SER :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // THR
    .{ "THR :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "THR :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "THR :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "THR :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "THR :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "THR :OG1 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "THR :CG2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "THR :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // TRP
    .{ "TRP :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "TRP :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TRP :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TRP :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "TRP :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TRP :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TRP :CD1 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TRP :CD2 ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TRP :NE1 ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "TRP :CE2 ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TRP :CE3 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TRP :CZ2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TRP :CZ3 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TRP :CH2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TRP :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // TYR
    .{ "TYR :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "TYR :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TYR :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TYR :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "TYR :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TYR :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TYR :CD1 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TYR :CD2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TYR :CE1 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TYR :CE2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "TYR :CZ  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TYR :OH  ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "TYR :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // VAL
    .{ "VAL :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "VAL :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "VAL :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "VAL :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "VAL :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "VAL :CG1 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "VAL :CG2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "VAL :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // === Non-standard Amino Acids ===

    // ASX (ambiguous ASN/ASP)
    .{ "ASX :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "ASX :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ASX :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ASX :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "ASX :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "ASX :CG  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "ASX :XD1 ", AtomProperties{ .radius = 1.50, .class = .polar } },
    .{ "ASX :XD2 ", AtomProperties{ .radius = 1.50, .class = .polar } },
    .{ "ASX :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // GLX (ambiguous GLN/GLU)
    .{ "GLX :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "GLX :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLX :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLX :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "GLX :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLX :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "GLX :CD  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "GLX :XE1 ", AtomProperties{ .radius = 1.50, .class = .polar } },
    .{ "GLX :XE2 ", AtomProperties{ .radius = 1.50, .class = .polar } },
    .{ "GLX :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // SEC (selenocysteine)
    .{ "SEC :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "SEC :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "SEC :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "SEC :SE  ", AtomProperties{ .radius = 1.90, .class = .polar } },
    .{ "SEC :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "SEC :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "SEC :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // MSE (selenomethionine)
    .{ "MSE :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "MSE :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MSE :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "MSE :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "MSE :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "MSE :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MSE :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MSE :SE  ", AtomProperties{ .radius = 1.90, .class = .polar } },
    .{ "MSE :CE  ", AtomProperties{ .radius = 1.88, .class = .apolar } },

    // PYL (pyrrolysine)
    .{ "PYL :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PYL :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PYL :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PYL :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PYL :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :CD  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :CE  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :NZ  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PYL :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PYL :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PYL :CA2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :N2  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PYL :CB2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :CG2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :CD2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PYL :CE2 ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // === New residues (not in ProtOr) ===

    // HYP (hydroxyproline) — like PRO but with hydroxyl on CG
    .{ "HYP :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "HYP :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "HYP :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "HYP :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "HYP :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } }, // sp3, 2H
    .{ "HYP :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } }, // sp3, 1H
    .{ "HYP :CD  ", AtomProperties{ .radius = 1.88, .class = .apolar } }, // sp3, 2H
    .{ "HYP :OD1 ", AtomProperties{ .radius = 1.46, .class = .polar } }, // sp3, 1H (hydroxyl)
    .{ "HYP :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // MLY (dimethyllysine) — like LYS but NZ has 0H (methylated)
    .{ "MLY :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "MLY :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MLY :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "MLY :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "MLY :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MLY :CG  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MLY :CD  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MLY :CE  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "MLY :NZ  ", AtomProperties{ .radius = 1.64, .class = .polar } }, // sp3, 0H
    .{ "MLY :CH1 ", AtomProperties{ .radius = 1.88, .class = .apolar } }, // methyl
    .{ "MLY :CH2 ", AtomProperties{ .radius = 1.88, .class = .apolar } }, // methyl
    .{ "MLY :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // SEP (phosphoserine) — like SER but OG bonded to phosphate
    .{ "SEP :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "SEP :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "SEP :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "SEP :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "SEP :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "SEP :OG  ", AtomProperties{ .radius = 1.46, .class = .polar } }, // sp3, 0H (bonded to P)
    .{ "SEP :P   ", AtomProperties{ .radius = 1.80, .class = .apolar } },
    .{ "SEP :O1P ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "SEP :O2P ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "SEP :O3P ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "SEP :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // TPO (phosphothreonine) — like THR but OG1 bonded to phosphate
    .{ "TPO :N   ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "TPO :CA  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TPO :C   ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "TPO :O   ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "TPO :CB  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TPO :CG2 ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "TPO :OG1 ", AtomProperties{ .radius = 1.46, .class = .polar } }, // sp3, 0H (bonded to P)
    .{ "TPO :P   ", AtomProperties{ .radius = 1.80, .class = .apolar } },
    .{ "TPO :O1P ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "TPO :O2P ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "TPO :O3P ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "TPO :OXT ", AtomProperties{ .radius = 1.46, .class = .polar } },

    // PSU (pseudouridine) — like U but C-glycosidic bond (C1'-C5 instead of N1)
    .{ "PSU :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PSU :P   ", AtomProperties{ .radius = 1.80, .class = .apolar } },
    .{ "PSU :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PSU :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PSU :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PSU :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PSU :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PSU :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PSU :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PSU :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PSU :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PSU :O2' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "PSU :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "PSU :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PSU :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PSU :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PSU :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "PSU :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "PSU :O4  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "PSU :C5  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "PSU :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // === Capping Groups ===

    // ACE (acetyl cap) — CH3-C(=O)
    .{ "ACE :C   ", AtomProperties{ .radius = 1.76, .class = .apolar } }, // sp2, 1H (C3H1 in ProtOr)
    .{ "ACE :O   ", AtomProperties{ .radius = 1.42, .class = .polar } }, // sp2, 0H
    .{ "ACE :CH3 ", AtomProperties{ .radius = 1.88, .class = .apolar } }, // sp3, 3H

    // NH2 (amino cap)
    .{ "NH2 :N   ", AtomProperties{ .radius = 1.64, .class = .polar } }, // sp3, 2H (N2H2 in ProtOr)

    // HOH (water)
    .{ "HOH :O   ", AtomProperties{ .radius = 1.46, .class = .polar } }, // sp3, 2H

    // === RNA Nucleotides ===

    // A (adenine)
    .{ "A   :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "A   :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "A   :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "A   :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "A   :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "A   :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "A   :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "A   :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "A   :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "A   :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "A   :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "A   :O2' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "A   :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "A   :N9  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "A   :C8  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "A   :N7  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "A   :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "A   :C6  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "A   :N6  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "A   :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "A   :C2  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "A   :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "A   :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },

    // C (cytosine)
    .{ "C   :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "C   :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "C   :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "C   :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "C   :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "C   :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "C   :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "C   :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "C   :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "C   :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "C   :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "C   :O2' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "C   :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "C   :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "C   :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "C   :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "C   :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "C   :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "C   :N4  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "C   :C5  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "C   :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // G (guanine)
    .{ "G   :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "G   :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "G   :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "G   :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "G   :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "G   :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "G   :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "G   :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "G   :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "G   :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "G   :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "G   :O2' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "G   :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "G   :N9  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "G   :C8  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "G   :N7  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "G   :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "G   :C6  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "G   :O6  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "G   :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "G   :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "G   :N2  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "G   :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "G   :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },

    // I (inosine)
    .{ "I   :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "I   :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "I   :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "I   :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "I   :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "I   :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "I   :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "I   :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "I   :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "I   :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "I   :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "I   :O2' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "I   :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "I   :N9  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "I   :C8  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "I   :N7  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "I   :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "I   :C6  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "I   :O6  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "I   :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "I   :C2  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "I   :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "I   :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },

    // T (thymine - for RNA)
    .{ "T   :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "T   :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "T   :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "T   :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "T   :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "T   :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "T   :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "T   :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "T   :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "T   :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "T   :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "T   :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "T   :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "T   :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "T   :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "T   :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "T   :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "T   :O4  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "T   :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "T   :C7  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "T   :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // U (uracil)
    .{ "U   :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "U   :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "U   :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "U   :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "U   :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "U   :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "U   :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "U   :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "U   :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "U   :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "U   :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "U   :O2' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "U   :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "U   :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "U   :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "U   :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "U   :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "U   :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "U   :O4  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "U   :C5  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "U   :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // === DNA Nucleotides ===

    // DA (deoxyadenosine)
    .{ "DA  :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DA  :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "DA  :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DA  :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DA  :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DA  :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DA  :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DA  :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DA  :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DA  :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DA  :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DA  :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DA  :N9  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DA  :C8  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DA  :N7  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DA  :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DA  :C6  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DA  :N6  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DA  :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DA  :C2  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DA  :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DA  :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },

    // DC (deoxycytidine)
    .{ "DC  :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DC  :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "DC  :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DC  :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DC  :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DC  :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DC  :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DC  :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DC  :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DC  :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DC  :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DC  :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DC  :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DC  :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DC  :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DC  :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DC  :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DC  :N4  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DC  :C5  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DC  :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // DG (deoxyguanosine)
    .{ "DG  :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DG  :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "DG  :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DG  :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DG  :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DG  :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DG  :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DG  :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DG  :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DG  :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DG  :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DG  :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DG  :N9  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DG  :C8  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DG  :N7  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DG  :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DG  :C6  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DG  :O6  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DG  :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DG  :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DG  :N2  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DG  :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DG  :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },

    // DI (deoxyinosine)
    .{ "DI  :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DI  :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "DI  :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DI  :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DI  :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DI  :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DI  :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DI  :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DI  :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DI  :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DI  :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DI  :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DI  :N9  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DI  :C8  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DI  :N7  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DI  :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DI  :C6  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DI  :O6  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DI  :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DI  :C2  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DI  :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DI  :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },

    // DT (deoxythymidine)
    .{ "DT  :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DT  :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "DT  :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DT  :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DT  :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DT  :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DT  :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DT  :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DT  :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DT  :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DT  :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DT  :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DT  :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DT  :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DT  :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DT  :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DT  :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DT  :O4  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DT  :C5  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DT  :C7  ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DT  :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },

    // DU (deoxyuridine)
    .{ "DU  :OP3 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DU  :P   ", AtomProperties{ .radius = 1.80, .class = .polar } },
    .{ "DU  :OP1 ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DU  :OP2 ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DU  :O5' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DU  :C5' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DU  :C4' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DU  :O4' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DU  :C3' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DU  :O3' ", AtomProperties{ .radius = 1.46, .class = .polar } },
    .{ "DU  :C2' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DU  :C1' ", AtomProperties{ .radius = 1.88, .class = .apolar } },
    .{ "DU  :N1  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DU  :C2  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DU  :O2  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DU  :N3  ", AtomProperties{ .radius = 1.64, .class = .polar } },
    .{ "DU  :C4  ", AtomProperties{ .radius = 1.61, .class = .apolar } },
    .{ "DU  :O4  ", AtomProperties{ .radius = 1.42, .class = .polar } },
    .{ "DU  :C5  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
    .{ "DU  :C6  ", AtomProperties{ .radius = 1.76, .class = .apolar } },
});

/// List of hardcoded residue names for isHardcoded check
const hardcoded_residues = std.StaticStringMap(void).initComptime(.{
    // Standard amino acids
    .{ "ALA", {} },
    .{ "ARG", {} },
    .{ "ASN", {} },
    .{ "ASP", {} },
    .{ "CYS", {} },
    .{ "GLN", {} },
    .{ "GLU", {} },
    .{ "GLY", {} },
    .{ "HIS", {} },
    .{ "ILE", {} },
    .{ "LEU", {} },
    .{ "LYS", {} },
    .{ "MET", {} },
    .{ "PHE", {} },
    .{ "PRO", {} },
    .{ "SER", {} },
    .{ "THR", {} },
    .{ "TRP", {} },
    .{ "TYR", {} },
    .{ "VAL", {} },
    // Non-standard amino acids
    .{ "ASX", {} },
    .{ "GLX", {} },
    .{ "SEC", {} },
    .{ "MSE", {} },
    .{ "PYL", {} },
    // New residues
    .{ "HYP", {} },
    .{ "MLY", {} },
    .{ "SEP", {} },
    .{ "TPO", {} },
    .{ "PSU", {} },
    // Capping groups
    .{ "ACE", {} },
    .{ "NH2", {} },
    // Solvent
    .{ "HOH", {} },
    // Nucleic acids
    .{ "A", {} },
    .{ "C", {} },
    .{ "G", {} },
    .{ "I", {} },
    .{ "T", {} },
    .{ "U", {} },
    .{ "DA", {} },
    .{ "DC", {} },
    .{ "DG", {} },
    .{ "DI", {} },
    .{ "DT", {} },
    .{ "DU", {} },
});

// =============================================================================
// Key formatting
// =============================================================================

/// Convert (residue, atom) to 9-char key for StaticStringMap lookup.
/// Format: "RES :ATOM" — 4 chars residue (left-justified, space-padded) + ':' + 4 chars atom (left-justified, space-padded)
pub fn formatKey(residue: []const u8, atom: []const u8) [9]u8 {
    var key: [9]u8 = .{ ' ', ' ', ' ', ' ', ':', ' ', ' ', ' ', ' ' };

    const res_len = @min(residue.len, 4);
    for (residue[0..res_len], 0..) |c, i| {
        key[i] = c;
    }

    const atom_len = @min(atom.len, 4);
    for (atom[0..atom_len], 0..) |c, i| {
        key[5 + i] = c;
    }

    return key;
}

// =============================================================================
// Runtime entry for non-standard residues
// =============================================================================

/// Runtime entry stored per atom for a component added via addComponent.
const RuntimeEntry = struct {
    props: AtomProperties,
};

/// Runtime component storage: maps formatted 9-char keys to properties.
const RuntimeMap = std.StringHashMap(RuntimeEntry);

// =============================================================================
// CcdClassifier
// =============================================================================

pub const CcdClassifier = struct {
    runtime_components: RuntimeMap,
    allocator: Allocator,

    const Self = @This();

    pub fn init(allocator: Allocator) CcdClassifier {
        return .{
            .runtime_components = RuntimeMap.init(allocator),
            .allocator = allocator,
        };
    }

    pub fn deinit(self: *CcdClassifier) void {
        // Free all allocated keys
        var it = self.runtime_components.keyIterator();
        while (it.next()) |key_ptr| {
            self.allocator.free(key_ptr.*);
        }
        self.runtime_components.deinit();
    }

    /// Get radius for a (residue, atom) pair.
    /// Lookup priority: hardcoded table -> runtime map -> null
    pub fn getRadius(self: *const CcdClassifier, residue: []const u8, atom: []const u8) ?f64 {
        const props = self.getProperties(residue, atom) orelse return null;
        return props.radius;
    }

    /// Get polarity class for a (residue, atom) pair.
    pub fn getClass(self: *const CcdClassifier, residue: []const u8, atom: []const u8) AtomClass {
        const props = self.getProperties(residue, atom) orelse return .unknown;
        return props.class;
    }

    /// Get both radius and class for a (residue, atom) pair.
    /// Lookup priority:
    /// 1. Hardcoded StaticStringMap (O(1), compile-time)
    /// 2. Runtime HashMap (populated via addComponent)
    /// 3. Return null (caller uses element fallback)
    pub fn getProperties(self: *const CcdClassifier, residue: []const u8, atom: []const u8) ?AtomProperties {
        // 1. Hardcoded lookup
        const key = formatKey(residue, atom);
        if (hardcoded_table.get(&key)) |props| {
            return props;
        }

        // 2. Runtime lookup
        if (self.runtime_components.get(&key)) |entry| {
            return entry.props;
        }

        // 3. Not found
        return null;
    }

    /// Check if a residue name has hardcoded entries.
    pub fn isHardcoded(residue: []const u8) bool {
        return hardcoded_residues.get(residue) != null;
    }

    /// Add a component from hybridization analysis to the runtime map.
    /// Calls deriveComponentProperties and stores results keyed by (comp_id, atom_id).
    pub fn addComponent(self: *CcdClassifier, component: *const Component) !void {
        const entries = try hybridization.deriveComponentProperties(self.allocator, component);
        defer self.allocator.free(entries);

        const comp_id = component.compIdSlice();

        for (entries) |entry| {
            const atom_id = entry.atomIdSlice();
            const key = formatKey(comp_id, atom_id);

            // Allocate a persistent copy of the key on the heap
            const key_copy = try self.allocator.alloc(u8, 9);
            @memcpy(key_copy, &key);

            try self.runtime_components.put(key_copy, RuntimeEntry{ .props = entry.props });
        }
    }
};

// =============================================================================
// Tests
// =============================================================================

const protor = @import("classifier_protor.zig");
const CompAtom = hybridization.CompAtom;
const CompBond = hybridization.CompBond;

test "CCD hardcoded ALA lookup — getRadius and getClass" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    // N
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), ccd.getRadius("ALA", "N").?, 0.001);
    try std.testing.expectEqual(AtomClass.polar, ccd.getClass("ALA", "N"));

    // CA
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), ccd.getRadius("ALA", "CA").?, 0.001);
    try std.testing.expectEqual(AtomClass.apolar, ccd.getClass("ALA", "CA"));

    // C
    try std.testing.expectApproxEqAbs(@as(f64, 1.61), ccd.getRadius("ALA", "C").?, 0.001);
    try std.testing.expectEqual(AtomClass.apolar, ccd.getClass("ALA", "C"));

    // O
    try std.testing.expectApproxEqAbs(@as(f64, 1.42), ccd.getRadius("ALA", "O").?, 0.001);
    try std.testing.expectEqual(AtomClass.polar, ccd.getClass("ALA", "O"));

    // CB
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), ccd.getRadius("ALA", "CB").?, 0.001);
    try std.testing.expectEqual(AtomClass.apolar, ccd.getClass("ALA", "CB"));
}

test "CCD unknown residue returns null" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    try std.testing.expectEqual(@as(?f64, null), ccd.getRadius("ZZZ", "N"));
    try std.testing.expectEqual(@as(?f64, null), ccd.getRadius("ZZZ", "CA"));
    try std.testing.expectEqual(AtomClass.unknown, ccd.getClass("ZZZ", "N"));
    try std.testing.expectEqual(@as(?AtomProperties, null), ccd.getProperties("ZZZ", "N"));
}

test "CCD isHardcoded" {
    try std.testing.expect(CcdClassifier.isHardcoded("ALA"));
    try std.testing.expect(CcdClassifier.isHardcoded("GLY"));
    try std.testing.expect(CcdClassifier.isHardcoded("PHE"));
    try std.testing.expect(CcdClassifier.isHardcoded("HOH"));
    try std.testing.expect(CcdClassifier.isHardcoded("MSE"));
    try std.testing.expect(CcdClassifier.isHardcoded("A"));
    try std.testing.expect(CcdClassifier.isHardcoded("DA"));
    try std.testing.expect(CcdClassifier.isHardcoded("HYP"));
    try std.testing.expect(CcdClassifier.isHardcoded("SEP"));

    try std.testing.expect(!CcdClassifier.isHardcoded("ZZZ"));
    try std.testing.expect(!CcdClassifier.isHardcoded("LIG"));
    try std.testing.expect(!CcdClassifier.isHardcoded(""));
}

test "CCD runtime component via addComponent" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    // Create a fake ligand "LIG" with 2 atoms: C1(C) and N1(N)
    const atoms = [_]CompAtom{
        CompAtom.init("C1", "C"), // 0
        CompAtom.init("N1", "N"), // 1
        CompAtom.init("H1", "H"), // 2 (hydrogen)
    };
    const bonds = [_]CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single, .aromatic = false },
        .{ .atom_idx_1 = 0, .atom_idx_2 = 2, .order = .single, .aromatic = false },
    };
    const comp = Component{
        .comp_id = .{ 'L', 'I', 'G', 0, 0 },
        .comp_id_len = 3,
        .atoms = &atoms,
        .bonds = &bonds,
    };

    // Before adding, should return null
    try std.testing.expectEqual(@as(?f64, null), ccd.getRadius("LIG", "C1"));

    try ccd.addComponent(&comp);

    // After adding, should find it
    const c1_radius = ccd.getRadius("LIG", "C1");
    try std.testing.expect(c1_radius != null);

    const n1_radius = ccd.getRadius("LIG", "N1");
    try std.testing.expect(n1_radius != null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), n1_radius.?, 0.001); // N always 1.64

    // H1 should not be added (hydrogen is skipped)
    try std.testing.expectEqual(@as(?f64, null), ccd.getRadius("LIG", "H1"));
}

test "CCD ProtOr compatibility — backbone atoms for all 20 standard amino acids" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    const residues = [_][]const u8{
        "ALA", "ARG", "ASN", "ASP", "CYS",
        "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO",
        "SER", "THR", "TRP", "TYR", "VAL",
    };

    const backbone_atoms = [_][]const u8{ "N", "CA", "C", "O" };

    for (residues) |res| {
        for (backbone_atoms) |atom| {
            const ccd_radius = ccd.getRadius(res, atom);
            const protor_radius = protor.getRadius(res, atom);

            // Both must return a value
            try std.testing.expect(ccd_radius != null);
            try std.testing.expect(protor_radius != null);

            // Values must match
            try std.testing.expectApproxEqAbs(protor_radius.?, ccd_radius.?, 0.001);
        }
    }
}

test "CCD HOH lookup" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    try std.testing.expectApproxEqAbs(@as(f64, 1.46), ccd.getRadius("HOH", "O").?, 0.001);
    try std.testing.expectEqual(AtomClass.polar, ccd.getClass("HOH", "O"));
}

test "CCD formatKey" {
    // Standard 3-letter residue + short atom name
    const key1 = formatKey("ALA", "N");
    try std.testing.expectEqualStrings("ALA :N   ", &key1);

    // 3-letter residue + 2-char atom
    const key2 = formatKey("ALA", "CA");
    try std.testing.expectEqualStrings("ALA :CA  ", &key2);

    // 1-letter residue (nucleic acid)
    const key3 = formatKey("A", "P");
    try std.testing.expectEqualStrings("A   :P   ", &key3);

    // 2-letter residue (DNA)
    const key4 = formatKey("DA", "C1'");
    try std.testing.expectEqualStrings("DA  :C1' ", &key4);

    // 4-char atom name
    const key5 = formatKey("ALA", "OXT");
    try std.testing.expectEqualStrings("ALA :OXT ", &key5);
}

test "CCD new residues — HYP, SEP, TPO, HOH" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    // HYP backbone matches standard AA
    try std.testing.expectApproxEqAbs(@as(f64, 1.64), ccd.getRadius("HYP", "N").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), ccd.getRadius("HYP", "CA").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.61), ccd.getRadius("HYP", "C").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.42), ccd.getRadius("HYP", "O").?, 0.001);
    // HYP specific: OD1 hydroxyl
    try std.testing.expectApproxEqAbs(@as(f64, 1.46), ccd.getRadius("HYP", "OD1").?, 0.001);

    // SEP phosphate
    try std.testing.expectApproxEqAbs(@as(f64, 1.80), ccd.getRadius("SEP", "P").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.42), ccd.getRadius("SEP", "O1P").?, 0.001);

    // TPO phosphate
    try std.testing.expectApproxEqAbs(@as(f64, 1.80), ccd.getRadius("TPO", "P").?, 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 1.42), ccd.getRadius("TPO", "O1P").?, 0.001);
}

test "CCD getProperties returns correct struct" {
    var ccd = CcdClassifier.init(std.testing.allocator);
    defer ccd.deinit();

    const props = ccd.getProperties("ALA", "CB");
    try std.testing.expect(props != null);
    try std.testing.expectApproxEqAbs(@as(f64, 1.88), props.?.radius, 0.001);
    try std.testing.expectEqual(AtomClass.apolar, props.?.class);
}

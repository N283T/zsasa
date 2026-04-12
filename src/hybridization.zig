//! Hybridization derivation module for CCD-based atom classification.
//!
//! Provides pure functions for deriving ProtOr-compatible VdW radii from
//! CCD (Chemical Component Dictionary) bond topology. Given a component's
//! atom and bond tables, this module analyzes bond orders to determine
//! hybridization states, then maps those to ProtOr radii and polarity classes.
//!
//! Reference: Tsai et al. 1999, J. Mol. Biol. 290:253-266

const std = @import("std");
const Allocator = std.mem.Allocator;
const classifier = @import("classifier.zig");
const AtomClass = classifier.AtomClass;
const AtomProperties = classifier.AtomProperties;

// =============================================================================
// Types
// =============================================================================

/// Bond order parsed from CCD `_chem_comp_bond.value_order` strings.
pub const BondOrder = enum {
    single,
    double,
    triple,
    aromatic,
    delocalized,
    unknown,

    /// Parse a CCD value_order string (case-insensitive).
    /// "SING" -> single, "DOUB" -> double, "TRIP" -> triple,
    /// "AROM" -> aromatic, "DELO" / "POLY" -> delocalized,
    /// anything else -> unknown.
    pub fn fromString(s: []const u8) BondOrder {
        if (s.len == 4) {
            var buf: [4]u8 = undefined;
            for (s[0..4], 0..) |c, i| {
                buf[i] = std.ascii.toUpper(c);
            }
            if (std.mem.eql(u8, &buf, "SING")) return .single;
            if (std.mem.eql(u8, &buf, "DOUB")) return .double;
            if (std.mem.eql(u8, &buf, "TRIP")) return .triple;
            if (std.mem.eql(u8, &buf, "AROM")) return .aromatic;
            if (std.mem.eql(u8, &buf, "DELO")) return .delocalized;
            if (std.mem.eql(u8, &buf, "POLY")) return .delocalized;
        }
        return .unknown;
    }
};

/// Hybridization state of an atom.
pub const Hybridization = enum {
    sp3,
    sp2,
    sp,
    unknown,
};

/// Fixed-size atom representation from a CCD component.
pub const CompAtom = struct {
    atom_id: [4]u8,
    atom_id_len: u3,
    type_symbol: [4]u8,
    type_symbol_len: u3,
    aromatic: bool,
    leaving: bool,

    /// Construct a CompAtom from slices (truncates to max field sizes).
    pub fn init(atom_id: []const u8, type_symbol: []const u8) CompAtom {
        var self: CompAtom = undefined;
        self.aromatic = false;
        self.leaving = false;

        const aid_len: usize = @min(atom_id.len, 4);
        self.atom_id = .{ 0, 0, 0, 0 };
        self.atom_id_len = @intCast(aid_len);
        for (atom_id[0..aid_len], 0..) |c, i| {
            self.atom_id[i] = c;
        }

        const ts_len: usize = @min(type_symbol.len, 4);
        self.type_symbol = .{ 0, 0, 0, 0 };
        self.type_symbol_len = @intCast(ts_len);
        for (type_symbol[0..ts_len], 0..) |c, i| {
            self.type_symbol[i] = c;
        }

        return self;
    }

    /// Get the atom_id as a slice.
    pub fn atomIdSlice(self: *const CompAtom) []const u8 {
        return self.atom_id[0..self.atom_id_len];
    }

    /// Get the type_symbol as a slice.
    pub fn typeSymbolSlice(self: *const CompAtom) []const u8 {
        return self.type_symbol[0..self.type_symbol_len];
    }
};

/// Bond between two atoms in a CCD component.
pub const CompBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: BondOrder,
    aromatic: bool,
};

/// A CCD component (residue/ligand) with atoms and bonds.
pub const Component = struct {
    comp_id: [5]u8,
    comp_id_len: u3,
    atoms: []const CompAtom,
    bonds: []const CompBond,

    /// Get the comp_id as a slice.
    pub fn compIdSlice(self: *const Component) []const u8 {
        return self.comp_id[0..self.comp_id_len];
    }
};

/// Result of analyzing bonds around a single atom.
pub const BondAnalysis = struct {
    hybridization: Hybridization,
    heavy_bond_count: u8,
    h_bond_count: u8,
    has_double: bool,
    has_triple: bool,
    has_aromatic: bool,
};

/// Result entry from deriveComponentProperties.
pub const DerivedAtomEntry = struct {
    atom_id: [4]u8,
    atom_id_len: u3,
    props: AtomProperties,

    pub fn atomIdSlice(self: *const DerivedAtomEntry) []const u8 {
        return self.atom_id[0..self.atom_id_len];
    }
};

// =============================================================================
// Bond Analysis
// =============================================================================

/// Check if a type_symbol represents hydrogen (case-insensitive).
fn isHydrogen(type_symbol: []const u8) bool {
    if (type_symbol.len == 1) {
        return type_symbol[0] == 'H' or type_symbol[0] == 'h';
    }
    return false;
}

/// Normalize an element symbol to uppercase into a fixed buffer.
/// Returns the length of the normalized symbol.
fn normalizeElement(type_symbol: []const u8, buf: *[4]u8) u3 {
    const trimmed = std.mem.trim(u8, type_symbol, " ");
    const len: usize = @min(trimmed.len, 4);
    for (trimmed[0..len], 0..) |c, i| {
        buf[i] = std.ascii.toUpper(c);
    }
    return @intCast(len);
}

/// Analyze bonds around a given atom in a component.
///
/// Iterates all bonds in the component to find those involving `atom_idx`,
/// counts heavy-atom vs hydrogen neighbors, and detects double/triple/aromatic
/// bonds to determine hybridization.
pub fn analyzeBonds(component: *const Component, atom_idx: u16) BondAnalysis {
    var result = BondAnalysis{
        .hybridization = .unknown,
        .heavy_bond_count = 0,
        .h_bond_count = 0,
        .has_double = false,
        .has_triple = false,
        .has_aromatic = false,
    };

    for (component.bonds) |bond| {
        const neighbor_idx: ?u16 = if (bond.atom_idx_1 == atom_idx)
            bond.atom_idx_2
        else if (bond.atom_idx_2 == atom_idx)
            bond.atom_idx_1
        else
            null;

        if (neighbor_idx) |nidx| {
            if (nidx >= component.atoms.len) continue;

            const neighbor = &component.atoms[nidx];
            if (isHydrogen(neighbor.typeSymbolSlice())) {
                result.h_bond_count += 1;
            } else {
                result.heavy_bond_count += 1;
            }

            switch (bond.order) {
                .double => result.has_double = true,
                .triple => result.has_triple = true,
                .aromatic => result.has_aromatic = true,
                .delocalized => result.has_double = true,
                else => {},
            }

            if (bond.aromatic) {
                result.has_aromatic = true;
            }
        }
    }

    // Determine hybridization from bond analysis
    if (result.has_triple) {
        result.hybridization = .sp;
    } else if (result.has_double or result.has_aromatic) {
        result.hybridization = .sp2;
    } else if (result.heavy_bond_count > 0 or result.h_bond_count > 0) {
        result.hybridization = .sp3;
    }

    return result;
}

// =============================================================================
// Radius Derivation (ProtOr table)
// =============================================================================

/// Derive ProtOr-compatible VdW radius from element, hybridization, and
/// implicit hydrogen count.
///
/// Returns null for elements not in the ProtOr table.
pub fn deriveRadius(type_symbol: []const u8, hybridization: Hybridization, implicit_h: u8) ?f64 {
    var buf: [4]u8 = undefined;
    const len = normalizeElement(type_symbol, &buf);
    const elem = buf[0..len];

    if (std.mem.eql(u8, elem, "C")) {
        return switch (hybridization) {
            .sp2 => if (implicit_h == 0) 1.61 else 1.76,
            .sp3 => if (implicit_h == 0) 1.61 else 1.88,
            .sp => 1.61,
            .unknown => 1.88,
        };
    }
    if (std.mem.eql(u8, elem, "N")) return 1.64;
    if (std.mem.eql(u8, elem, "O")) {
        return switch (hybridization) {
            .sp2 => 1.42,
            .sp3 => 1.46,
            .sp => 1.42,
            .unknown => 1.42,
        };
    }
    if (std.mem.eql(u8, elem, "S")) return 1.77;
    if (std.mem.eql(u8, elem, "P")) return 1.80;
    if (std.mem.eql(u8, elem, "SE")) return 1.90;

    return null;
}

// =============================================================================
// Atom Class Derivation
// =============================================================================

/// Derive polarity class from element symbol.
///
/// N, O -> polar; C, S, P, SE -> apolar; others -> unknown.
pub fn deriveClass(type_symbol: []const u8) AtomClass {
    var buf: [4]u8 = undefined;
    const len = normalizeElement(type_symbol, &buf);
    const elem = buf[0..len];

    if (std.mem.eql(u8, elem, "N") or std.mem.eql(u8, elem, "O")) return .polar;
    if (std.mem.eql(u8, elem, "C") or
        std.mem.eql(u8, elem, "S") or
        std.mem.eql(u8, elem, "P") or
        std.mem.eql(u8, elem, "SE")) return .apolar;

    return .unknown;
}

// =============================================================================
// Valence / Implicit Hydrogen
// =============================================================================

/// Typical valence for common biological elements.
fn typicalValence(type_symbol: []const u8) ?u8 {
    var buf: [4]u8 = undefined;
    const len = normalizeElement(type_symbol, &buf);
    const elem = buf[0..len];

    if (std.mem.eql(u8, elem, "C")) return 4;
    if (std.mem.eql(u8, elem, "N")) return 3;
    if (std.mem.eql(u8, elem, "O")) return 2;
    if (std.mem.eql(u8, elem, "S")) return 2;
    if (std.mem.eql(u8, elem, "P")) return 5;
    if (std.mem.eql(u8, elem, "SE")) return 2;

    return null;
}

/// Compute the number of implicit hydrogens for an atom.
///
/// implicit_h = typical_valence - (heavy_bond_count + explicit_h_count),
/// clamped to 0. Returns 0 if the element has no known typical valence.
pub fn implicitHCount(type_symbol: []const u8, heavy_bond_count: u8, explicit_h_count: u8) u8 {
    const valence = typicalValence(type_symbol) orelse return 0;
    const total_bonds = @as(u16, heavy_bond_count) + @as(u16, explicit_h_count);
    if (total_bonds >= valence) return 0;
    return @intCast(valence - total_bonds);
}

// =============================================================================
// Component-Level Derivation
// =============================================================================

/// Derive ProtOr-compatible atom properties for all non-hydrogen atoms
/// in a component.
///
/// For each non-hydrogen atom:
/// 1. Analyze bonds to determine hybridization
/// 2. Compute implicit hydrogen count
/// 3. Derive ProtOr radius and polarity class
/// 4. Skip atoms where deriveRadius returns null (unsupported elements)
///
/// Returns an owned slice that the caller must free with the same allocator.
pub fn deriveComponentProperties(
    allocator: Allocator,
    component: *const Component,
) ![]const DerivedAtomEntry {
    var results = std.ArrayListUnmanaged(DerivedAtomEntry){};
    errdefer results.deinit(allocator);

    for (component.atoms, 0..) |atom, idx| {
        // Skip hydrogen atoms
        if (isHydrogen(atom.typeSymbolSlice())) continue;

        const atom_idx: u16 = @intCast(idx);
        const analysis = analyzeBonds(component, atom_idx);

        const implicit_h = implicitHCount(
            atom.typeSymbolSlice(),
            analysis.heavy_bond_count,
            analysis.h_bond_count,
        );

        const radius = deriveRadius(
            atom.typeSymbolSlice(),
            analysis.hybridization,
            implicit_h,
        ) orelse continue; // skip unsupported elements

        const class = deriveClass(atom.typeSymbolSlice());

        try results.append(allocator, .{
            .atom_id = atom.atom_id,
            .atom_id_len = atom.atom_id_len,
            .props = .{
                .radius = radius,
                .class = class,
            },
        });
    }

    return results.toOwnedSlice(allocator);
}

// =============================================================================
// Tests
// =============================================================================

test "BondOrder.fromString — all variants and case insensitivity" {
    // Standard uppercase
    try std.testing.expectEqual(BondOrder.single, BondOrder.fromString("SING"));
    try std.testing.expectEqual(BondOrder.double, BondOrder.fromString("DOUB"));
    try std.testing.expectEqual(BondOrder.triple, BondOrder.fromString("TRIP"));
    try std.testing.expectEqual(BondOrder.aromatic, BondOrder.fromString("AROM"));
    try std.testing.expectEqual(BondOrder.delocalized, BondOrder.fromString("DELO"));
    try std.testing.expectEqual(BondOrder.delocalized, BondOrder.fromString("POLY"));

    // Case insensitivity
    try std.testing.expectEqual(BondOrder.single, BondOrder.fromString("sing"));
    try std.testing.expectEqual(BondOrder.double, BondOrder.fromString("Doub"));
    try std.testing.expectEqual(BondOrder.triple, BondOrder.fromString("tRiP"));
    try std.testing.expectEqual(BondOrder.aromatic, BondOrder.fromString("arom"));

    // Unknown
    try std.testing.expectEqual(BondOrder.unknown, BondOrder.fromString("XXXX"));
    try std.testing.expectEqual(BondOrder.unknown, BondOrder.fromString(""));
    try std.testing.expectEqual(BondOrder.unknown, BondOrder.fromString("SI"));
}

test "analyzeBonds — sp3 atom (ALA CA)" {
    // ALA CA: bonded to N(single), C(single), CB(single), HA(single)
    // => sp3, 3 heavy, 1 H, no double/triple/aromatic
    const atoms = [_]CompAtom{
        CompAtom.init("N", "N"), // 0
        CompAtom.init("CA", "C"), // 1
        CompAtom.init("C", "C"), // 2
        CompAtom.init("CB", "C"), // 3
        CompAtom.init("HA", "H"), // 4
    };
    const bonds = [_]CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single, .aromatic = false }, // N-CA
        .{ .atom_idx_1 = 1, .atom_idx_2 = 2, .order = .single, .aromatic = false }, // CA-C
        .{ .atom_idx_1 = 1, .atom_idx_2 = 3, .order = .single, .aromatic = false }, // CA-CB
        .{ .atom_idx_1 = 1, .atom_idx_2 = 4, .order = .single, .aromatic = false }, // CA-HA
    };
    const comp = Component{
        .comp_id = .{ 'A', 'L', 'A', 0, 0 },
        .comp_id_len = 3,
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const result = analyzeBonds(&comp, 1); // CA
    try std.testing.expectEqual(Hybridization.sp3, result.hybridization);
    try std.testing.expectEqual(@as(u8, 3), result.heavy_bond_count);
    try std.testing.expectEqual(@as(u8, 1), result.h_bond_count);
    try std.testing.expect(!result.has_double);
    try std.testing.expect(!result.has_triple);
    try std.testing.expect(!result.has_aromatic);
}

test "analyzeBonds — sp2 atom (carbonyl C)" {
    // Carbonyl C: bonded to CA(single), O(double), N_next(single)
    // => sp2, 3 heavy, 0 H, has_double=true
    const atoms = [_]CompAtom{
        CompAtom.init("CA", "C"), // 0
        CompAtom.init("C", "C"), // 1
        CompAtom.init("O", "O"), // 2
        CompAtom.init("N", "N"), // 3 (next residue's N, modeled in same component)
    };
    const bonds = [_]CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single, .aromatic = false }, // CA-C
        .{ .atom_idx_1 = 1, .atom_idx_2 = 2, .order = .double, .aromatic = false }, // C=O
        .{ .atom_idx_1 = 1, .atom_idx_2 = 3, .order = .single, .aromatic = false }, // C-N
    };
    const comp = Component{
        .comp_id = .{ 'A', 'L', 'A', 0, 0 },
        .comp_id_len = 3,
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const result = analyzeBonds(&comp, 1); // C
    try std.testing.expectEqual(Hybridization.sp2, result.hybridization);
    try std.testing.expectEqual(@as(u8, 3), result.heavy_bond_count);
    try std.testing.expectEqual(@as(u8, 0), result.h_bond_count);
    try std.testing.expect(result.has_double);
    try std.testing.expect(!result.has_triple);
    try std.testing.expect(!result.has_aromatic);
}

test "analyzeBonds — aromatic atom" {
    // PHE CG: aromatic ring atom bonded to CB(single), CD1(arom), CD2(arom)
    const atoms = [_]CompAtom{
        CompAtom.init("CB", "C"), // 0
        CompAtom.init("CG", "C"), // 1
        CompAtom.init("CD1", "C"), // 2
        CompAtom.init("CD2", "C"), // 3
    };
    const bonds = [_]CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single, .aromatic = false }, // CB-CG
        .{ .atom_idx_1 = 1, .atom_idx_2 = 2, .order = .aromatic, .aromatic = true }, // CG-CD1
        .{ .atom_idx_1 = 1, .atom_idx_2 = 3, .order = .aromatic, .aromatic = true }, // CG-CD2
    };
    const comp = Component{
        .comp_id = .{ 'P', 'H', 'E', 0, 0 },
        .comp_id_len = 3,
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const result = analyzeBonds(&comp, 1); // CG
    try std.testing.expectEqual(Hybridization.sp2, result.hybridization);
    try std.testing.expectEqual(@as(u8, 3), result.heavy_bond_count);
    try std.testing.expectEqual(@as(u8, 0), result.h_bond_count);
    try std.testing.expect(!result.has_double);
    try std.testing.expect(!result.has_triple);
    try std.testing.expect(result.has_aromatic);
}

test "deriveRadius — all element/hybridization combos" {
    // Carbon
    try std.testing.expectEqual(@as(?f64, 1.61), deriveRadius("C", .sp2, 0));
    try std.testing.expectEqual(@as(?f64, 1.76), deriveRadius("C", .sp2, 1));
    try std.testing.expectEqual(@as(?f64, 1.76), deriveRadius("C", .sp2, 2));
    try std.testing.expectEqual(@as(?f64, 1.61), deriveRadius("C", .sp3, 0));
    try std.testing.expectEqual(@as(?f64, 1.88), deriveRadius("C", .sp3, 1));
    try std.testing.expectEqual(@as(?f64, 1.88), deriveRadius("C", .sp3, 2));
    try std.testing.expectEqual(@as(?f64, 1.88), deriveRadius("C", .sp3, 3));
    try std.testing.expectEqual(@as(?f64, 1.61), deriveRadius("C", .sp, 0));
    try std.testing.expectEqual(@as(?f64, 1.88), deriveRadius("C", .unknown, 0));

    // Nitrogen (always 1.64)
    try std.testing.expectEqual(@as(?f64, 1.64), deriveRadius("N", .sp2, 0));
    try std.testing.expectEqual(@as(?f64, 1.64), deriveRadius("N", .sp3, 1));
    try std.testing.expectEqual(@as(?f64, 1.64), deriveRadius("N", .sp, 0));
    try std.testing.expectEqual(@as(?f64, 1.64), deriveRadius("N", .unknown, 0));

    // Oxygen
    try std.testing.expectEqual(@as(?f64, 1.42), deriveRadius("O", .sp2, 0));
    try std.testing.expectEqual(@as(?f64, 1.46), deriveRadius("O", .sp3, 0));
    try std.testing.expectEqual(@as(?f64, 1.42), deriveRadius("O", .sp, 0));
    try std.testing.expectEqual(@as(?f64, 1.42), deriveRadius("O", .unknown, 0));

    // Sulfur (always 1.77)
    try std.testing.expectEqual(@as(?f64, 1.77), deriveRadius("S", .sp3, 0));
    try std.testing.expectEqual(@as(?f64, 1.77), deriveRadius("S", .sp2, 0));

    // Phosphorus (always 1.80)
    try std.testing.expectEqual(@as(?f64, 1.80), deriveRadius("P", .sp3, 0));

    // Selenium (always 1.90)
    try std.testing.expectEqual(@as(?f64, 1.90), deriveRadius("SE", .sp3, 0));

    // Case insensitivity
    try std.testing.expectEqual(@as(?f64, 1.61), deriveRadius("c", .sp2, 0));
    try std.testing.expectEqual(@as(?f64, 1.64), deriveRadius("n", .sp3, 0));
    try std.testing.expectEqual(@as(?f64, 1.90), deriveRadius("Se", .sp3, 0));
}

test "deriveRadius — unknown element returns null" {
    try std.testing.expectEqual(@as(?f64, null), deriveRadius("FE", .sp3, 0));
    try std.testing.expectEqual(@as(?f64, null), deriveRadius("ZN", .sp3, 0));
    try std.testing.expectEqual(@as(?f64, null), deriveRadius("XX", .sp3, 0));
}

test "deriveClass — all elements" {
    // Polar
    try std.testing.expectEqual(AtomClass.polar, deriveClass("N"));
    try std.testing.expectEqual(AtomClass.polar, deriveClass("O"));
    try std.testing.expectEqual(AtomClass.polar, deriveClass("n"));
    try std.testing.expectEqual(AtomClass.polar, deriveClass("o"));

    // Apolar
    try std.testing.expectEqual(AtomClass.apolar, deriveClass("C"));
    try std.testing.expectEqual(AtomClass.apolar, deriveClass("S"));
    try std.testing.expectEqual(AtomClass.apolar, deriveClass("P"));
    try std.testing.expectEqual(AtomClass.apolar, deriveClass("SE"));
    try std.testing.expectEqual(AtomClass.apolar, deriveClass("c"));
    try std.testing.expectEqual(AtomClass.apolar, deriveClass("Se"));

    // Unknown
    try std.testing.expectEqual(AtomClass.unknown, deriveClass("FE"));
    try std.testing.expectEqual(AtomClass.unknown, deriveClass("ZN"));
    try std.testing.expectEqual(AtomClass.unknown, deriveClass("H"));
}

test "implicitHCount — C with 3 heavy bonds" {
    // C valence=4, 3 heavy bonds, 0 explicit H => 1 implicit H
    try std.testing.expectEqual(@as(u8, 1), implicitHCount("C", 3, 0));
}

test "implicitHCount — N with 1 heavy bond" {
    // N valence=3, 1 heavy bond, 0 explicit H => 2 implicit H
    try std.testing.expectEqual(@as(u8, 2), implicitHCount("N", 1, 0));
}

test "implicitHCount — with explicit H (=0)" {
    // N valence=3, 1 heavy bond, 2 explicit H => 0 implicit H
    try std.testing.expectEqual(@as(u8, 0), implicitHCount("N", 1, 2));
}

test "implicitHCount — clamped to zero" {
    // If bonds exceed valence, result is 0
    try std.testing.expectEqual(@as(u8, 0), implicitHCount("C", 5, 0));
    try std.testing.expectEqual(@as(u8, 0), implicitHCount("O", 3, 0));
}

test "implicitHCount — unknown element returns 0" {
    try std.testing.expectEqual(@as(u8, 0), implicitHCount("FE", 2, 0));
}

test "deriveComponentProperties — simplified ALA" {
    // Simplified ALA without H atoms: N, CA, C, O, CB
    // Bonds: N-CA(single), CA-C(single), C=O(double), CA-CB(single)
    // Plus backbone N-C bond for completeness: not included since
    // this is intra-residue only.
    const atoms = [_]CompAtom{
        CompAtom.init("N", "N"), // 0: N
        CompAtom.init("CA", "C"), // 1: CA
        CompAtom.init("C", "C"), // 2: C (carbonyl)
        CompAtom.init("O", "O"), // 3: O
        CompAtom.init("CB", "C"), // 4: CB
    };
    const bonds = [_]CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single, .aromatic = false }, // N-CA
        .{ .atom_idx_1 = 1, .atom_idx_2 = 2, .order = .single, .aromatic = false }, // CA-C
        .{ .atom_idx_1 = 2, .atom_idx_2 = 3, .order = .double, .aromatic = false }, // C=O
        .{ .atom_idx_1 = 1, .atom_idx_2 = 4, .order = .single, .aromatic = false }, // CA-CB
    };
    const comp = Component{
        .comp_id = .{ 'A', 'L', 'A', 0, 0 },
        .comp_id_len = 3,
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const results = try deriveComponentProperties(std.testing.allocator, &comp);
    defer std.testing.allocator.free(results);

    // Should have 5 non-hydrogen atoms, all with known radii
    try std.testing.expectEqual(@as(usize, 5), results.len);

    // Build a lookup by atom_id for easier assertions
    var found_n = false;
    var found_ca = false;
    var found_c = false;
    var found_o = false;
    var found_cb = false;

    for (results) |entry| {
        const name = entry.atomIdSlice();
        if (std.mem.eql(u8, name, "N")) {
            // N: 1 heavy bond (CA), 0 H bonds
            // implicit_h = 3 - (1+0) = 2
            // radius = 1.64 (N always 1.64)
            try std.testing.expectEqual(@as(f64, 1.64), entry.props.radius);
            try std.testing.expectEqual(AtomClass.polar, entry.props.class);
            found_n = true;
        } else if (std.mem.eql(u8, name, "CA")) {
            // CA: 3 heavy bonds (N, C, CB), 0 H bonds
            // implicit_h = 4 - (3+0) = 1
            // hybridization = sp3 (all single bonds)
            // radius = sp3/H1+ = 1.88
            try std.testing.expectEqual(@as(f64, 1.88), entry.props.radius);
            try std.testing.expectEqual(AtomClass.apolar, entry.props.class);
            found_ca = true;
        } else if (std.mem.eql(u8, name, "C")) {
            // C (carbonyl): 2 heavy bonds (CA, O), 0 H bonds
            // implicit_h = 4 - (2+0) = 2
            // hybridization = sp2 (has double bond C=O)
            // radius = sp2/H1+ = 1.76
            try std.testing.expectEqual(@as(f64, 1.76), entry.props.radius);
            try std.testing.expectEqual(AtomClass.apolar, entry.props.class);
            found_c = true;
        } else if (std.mem.eql(u8, name, "O")) {
            // O: 1 heavy bond (C), 0 H bonds
            // implicit_h = 2 - (1+0) = 1
            // hybridization = sp2 (double bond C=O)
            // radius = sp2 = 1.42
            try std.testing.expectEqual(@as(f64, 1.42), entry.props.radius);
            try std.testing.expectEqual(AtomClass.polar, entry.props.class);
            found_o = true;
        } else if (std.mem.eql(u8, name, "CB")) {
            // CB: 1 heavy bond (CA), 0 H bonds
            // implicit_h = 4 - (1+0) = 3
            // hybridization = sp3 (single bonds only)
            // radius = sp3/H1+ = 1.88
            try std.testing.expectEqual(@as(f64, 1.88), entry.props.radius);
            try std.testing.expectEqual(AtomClass.apolar, entry.props.class);
            found_cb = true;
        }
    }

    try std.testing.expect(found_n);
    try std.testing.expect(found_ca);
    try std.testing.expect(found_c);
    try std.testing.expect(found_o);
    try std.testing.expect(found_cb);
}

test "deriveComponentProperties — skips hydrogen atoms and unknown elements" {
    const atoms = [_]CompAtom{
        CompAtom.init("C1", "C"), // 0
        CompAtom.init("H1", "H"), // 1 (hydrogen, should be skipped)
        CompAtom.init("FE", "FE"), // 2 (unknown element, should be skipped)
    };
    const bonds = [_]CompBond{
        .{ .atom_idx_1 = 0, .atom_idx_2 = 1, .order = .single, .aromatic = false },
    };
    const comp = Component{
        .comp_id = .{ 'T', 'S', 'T', 0, 0 },
        .comp_id_len = 3,
        .atoms = &atoms,
        .bonds = &bonds,
    };

    const results = try deriveComponentProperties(std.testing.allocator, &comp);
    defer std.testing.allocator.free(results);

    // Only C1 should appear (H1 is hydrogen, FE has no ProtOr radius)
    try std.testing.expectEqual(@as(usize, 1), results.len);
    try std.testing.expectEqualStrings("C1", results[0].atomIdSlice());
}

test "CompAtom init and accessors" {
    const atom = CompAtom.init("CA", "C");
    try std.testing.expectEqualStrings("CA", atom.atomIdSlice());
    try std.testing.expectEqualStrings("C", atom.typeSymbolSlice());
    try std.testing.expect(!atom.aromatic);
    try std.testing.expect(!atom.leaving);
}

test "CompAtom truncation" {
    // atom_id max is 4, type_symbol max is 4
    const atom = CompAtom.init("TOOLONG", "VERYLONGELEMENT");
    try std.testing.expectEqualStrings("TOOL", atom.atomIdSlice());
    try std.testing.expectEqualStrings("VERY", atom.typeSymbolSlice());
}

test "Component compIdSlice" {
    const comp = Component{
        .comp_id = .{ 'A', 'L', 'A', 0, 0 },
        .comp_id_len = 3,
        .atoms = &.{},
        .bonds = &.{},
    };
    try std.testing.expectEqualStrings("ALA", comp.compIdSlice());
}

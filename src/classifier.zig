//! Atom radius classifier for SASA calculation.
//!
//! This module provides classification of atoms based on residue name and atom name,
//! returning appropriate van der Waals radii and polarity classes.
//!
//! The classifier follows FreeSASA's approach:
//! 1. Look up by (residue_name, atom_name)
//! 2. Fall back to ("ANY", atom_name)
//! 3. Fall back to element-based radius guessing
//!
//! ## Available Built-in Classifiers
//!
//! - **CCD**: Default. Derives radii from CCD bond topology (ProtOr-compatible).
//!   ProtOr is an alias for CCD — they use the same hardcoded radii, but CCD
//!   additionally supports runtime analysis of non-standard residues.
//! - **NACCESS**: NACCESS-compatible radii (João Rodrigues)
//! - **OONS**: Ooi, Oobatake, Nemethy, Scheraga radii (older FreeSASA default)
//!
//! ## Usage
//!
//! ```zig
//! const classifier = @import("classifier.zig");
//! const naccess = @import("classifier_naccess.zig");
//! const ccd = @import("classifier_ccd.zig");
//! const oons = @import("classifier_oons.zig");
//!
//! // Get classifier by type
//! const classifier_type = classifier.ClassifierType.ccd;
//!
//! // Or use built-in classifiers directly
//! const r = naccess.getRadius("ALA", "CA");
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;

/// Built-in classifier types
pub const ClassifierType = enum {
    /// NACCESS-compatible radii
    /// Based on naccess.config from FreeSASA (João Rodrigues)
    naccess,

    /// Alias for CCD — kept for backward compatibility
    /// Uses the same CCD classifier internally
    protor,

    /// OONS radii (older FreeSASA default)
    /// Reference: Ooi, Oobatake, Nemethy, Scheraga
    oons,

    /// CCD-based radii derived from bond topology (default)
    /// Hardcoded ProtOr radii + runtime CCD analysis for non-standard residues
    ccd,

    pub fn name(self: ClassifierType) []const u8 {
        return switch (self) {
            .naccess => "NACCESS",
            .protor => "ProtOr",
            .oons => "OONS",
            .ccd => "CCD",
        };
    }

    pub fn fromString(s: []const u8) ?ClassifierType {
        if (std.ascii.eqlIgnoreCase(s, "naccess")) return .naccess;
        if (std.ascii.eqlIgnoreCase(s, "protor")) return .protor;
        if (std.ascii.eqlIgnoreCase(s, "oons")) return .oons;
        if (std.ascii.eqlIgnoreCase(s, "ccd")) return .ccd;
        return null;
    }
};

/// Atom polarity classification
pub const AtomClass = enum {
    polar,
    apolar,
    unknown,

    pub fn fromString(s: []const u8) AtomClass {
        if (std.mem.eql(u8, s, "polar")) return .polar;
        if (std.mem.eql(u8, s, "apolar")) return .apolar;
        return .unknown;
    }

    pub fn toString(self: AtomClass) []const u8 {
        return switch (self) {
            .polar => "polar",
            .apolar => "apolar",
            .unknown => "unknown",
        };
    }
};

/// Properties of an atom type (radius and class)
pub const AtomProperties = struct {
    radius: f64,
    class: AtomClass,
};

/// Key for atom lookup: (residue_name, atom_name)
/// Uses fixed-size arrays to avoid allocation for keys
pub const AtomKey = struct {
    residue: [5]u8, // mmCIF comp_id can be up to 5 chars
    residue_len: u8,
    atom: [4]u8, // PDB atom names are max 4 chars
    atom_len: u8,

    pub fn init(residue: []const u8, atom: []const u8) AtomKey {
        var self: AtomKey = undefined;
        self.initInPlace(residue, atom);
        return self;
    }

    pub fn initInPlace(self: *AtomKey, residue: []const u8, atom: []const u8) void {
        const res_len: usize = @min(residue.len, 5);
        const atm_len: usize = @min(atom.len, 4);

        self.residue = .{ 0, 0, 0, 0, 0 };
        self.atom = .{ 0, 0, 0, 0 };
        self.residue_len = @intCast(res_len);
        self.atom_len = @intCast(atm_len);

        for (residue[0..res_len], 0..) |c, i| {
            self.residue[i] = c;
        }
        for (atom[0..atm_len], 0..) |c, i| {
            self.atom[i] = c;
        }
    }

    pub fn residueName(self: *const AtomKey) []const u8 {
        return self.residue[0..self.residue_len];
    }

    pub fn atomName(self: *const AtomKey) []const u8 {
        return self.atom[0..self.atom_len];
    }

    pub fn eql(a: AtomKey, b: AtomKey) bool {
        return a.residue_len == b.residue_len and
            a.atom_len == b.atom_len and
            std.mem.eql(u8, a.residue[0..a.residue_len], b.residue[0..b.residue_len]) and
            std.mem.eql(u8, a.atom[0..a.atom_len], b.atom[0..b.atom_len]);
    }

    pub fn hash(self: AtomKey) u64 {
        var h = std.hash.Wyhash.init(0);
        h.update(self.residue[0..self.residue_len]);
        h.update(&[_]u8{0}); // separator
        h.update(self.atom[0..self.atom_len]);
        return h.final();
    }
};

/// Hash map context for AtomKey
pub const AtomKeyContext = struct {
    pub fn hash(_: AtomKeyContext, key: AtomKey) u64 {
        return key.hash();
    }

    pub fn eql(_: AtomKeyContext, a: AtomKey, b: AtomKey) bool {
        return a.eql(b);
    }
};

/// Atom radius and class classifier
pub const Classifier = struct {
    name: []const u8,
    /// Map from (residue, atom) -> AtomProperties
    atoms: std.HashMap(AtomKey, AtomProperties, AtomKeyContext, 80),
    allocator: Allocator,

    const Self = @This();

    /// Initialize an empty classifier
    pub fn init(allocator: Allocator, name: []const u8) !Self {
        const name_copy = try allocator.dupe(u8, name);
        return Self{
            .name = name_copy,
            .atoms = std.HashMap(AtomKey, AtomProperties, AtomKeyContext, 80).init(allocator),
            .allocator = allocator,
        };
    }

    /// Free all resources
    pub fn deinit(self: *Self) void {
        self.allocator.free(self.name);
        self.atoms.deinit();
    }

    /// Add an atom type to the classifier
    pub fn addAtom(
        self: *Self,
        residue: []const u8,
        atom: []const u8,
        radius: f64,
        class: AtomClass,
    ) !void {
        const key = AtomKey.init(residue, atom);
        try self.atoms.put(key, AtomProperties{ .radius = radius, .class = class });
    }

    /// Get radius for an atom, with ANY fallback
    /// Returns null if not found
    pub fn getRadius(self: *const Self, residue: []const u8, atom: []const u8) ?f64 {
        // Try exact match first
        const key = AtomKey.init(residue, atom);
        if (self.atoms.get(key)) |props| {
            return props.radius;
        }

        // Fall back to ANY
        const any_key = AtomKey.init("ANY", atom);
        if (self.atoms.get(any_key)) |props| {
            return props.radius;
        }

        return null;
    }

    /// Get class for an atom, with ANY fallback
    pub fn getClass(self: *const Self, residue: []const u8, atom: []const u8) AtomClass {
        // Try exact match first
        const key = AtomKey.init(residue, atom);
        if (self.atoms.get(key)) |props| {
            return props.class;
        }

        // Fall back to ANY
        const any_key = AtomKey.init("ANY", atom);
        if (self.atoms.get(any_key)) |props| {
            return props.class;
        }

        return .unknown;
    }

    /// Get both radius and class for an atom
    pub fn getProperties(self: *const Self, residue: []const u8, atom: []const u8) ?AtomProperties {
        // Try exact match first
        const key = AtomKey.init(residue, atom);
        if (self.atoms.get(key)) |props| {
            return props;
        }

        // Fall back to ANY
        const any_key = AtomKey.init("ANY", atom);
        if (self.atoms.get(any_key)) |props| {
            return props;
        }

        return null;
    }

    /// Get the number of atom entries
    pub fn count(self: *const Self) usize {
        return self.atoms.count();
    }
};

// =============================================================================
// Element-Based Radius Guessing
// =============================================================================

/// van der Waals radii by atomic number (O(1) array lookup)
/// Sources: Mantina et al. 2009 for common elements, gemmi for others
/// Index = atomic number, Value = radius in Angstroms (0 = unknown)
const atomic_number_radii = blk: {
    var radii: [119]f64 = .{0} ** 119;
    // Common biological elements
    radii[1] = 1.10; // H
    radii[6] = 1.70; // C
    radii[7] = 1.55; // N
    radii[8] = 1.52; // O
    radii[15] = 1.80; // P
    radii[16] = 1.80; // S
    radii[34] = 1.90; // Se
    // Halogens
    radii[9] = 1.47; // F
    radii[17] = 1.75; // Cl
    radii[35] = 1.83; // Br
    radii[53] = 1.98; // I
    // Alkali and Alkali Earth metals
    radii[3] = 1.81; // Li
    radii[4] = 1.53; // Be
    radii[11] = 2.27; // Na
    radii[12] = 1.73; // Mg
    radii[19] = 2.75; // K
    radii[20] = 2.31; // Ca
    // Transition metals (common in proteins)
    radii[25] = 1.19; // Mn
    radii[26] = 1.26; // Fe
    radii[27] = 1.13; // Co
    radii[28] = 1.63; // Ni
    radii[29] = 1.40; // Cu
    radii[30] = 1.39; // Zn
    // Other metals
    radii[33] = 1.85; // As
    radii[48] = 1.58; // Cd
    radii[80] = 1.55; // Hg
    radii[82] = 2.02; // Pb
    break :blk radii;
};

/// Get van der Waals radius from atomic number
/// Returns null if atomic number is not recognized or out of range
pub fn guessRadiusFromAtomicNumber(atomic_number: u8) ?f64 {
    if (atomic_number == 0 or atomic_number >= atomic_number_radii.len) return null;
    const radius = atomic_number_radii[atomic_number];
    if (radius == 0) return null;
    return radius;
}

/// van der Waals radii map (O(1) lookup)
/// Sources: Mantina et al. 2009 for common elements, gemmi for others
/// Keys are uppercase element symbols
const element_radii_map = std.StaticStringMap(f64).initComptime(.{
    // Common biological elements
    .{ "H", 1.10 },
    .{ "C", 1.70 },
    .{ "N", 1.55 },
    .{ "O", 1.52 },
    .{ "P", 1.80 },
    .{ "S", 1.80 },
    .{ "SE", 1.90 },
    // Halogens
    .{ "F", 1.47 },
    .{ "CL", 1.75 },
    .{ "BR", 1.83 },
    .{ "I", 1.98 },
    // Alkali and Alkali Earth metals
    .{ "LI", 1.81 },
    .{ "BE", 1.53 },
    .{ "NA", 2.27 },
    .{ "MG", 1.73 },
    .{ "K", 2.75 },
    .{ "CA", 2.31 },
    // Transition metals (common in proteins)
    .{ "FE", 1.26 },
    .{ "CO", 1.13 },
    .{ "NI", 1.63 },
    .{ "CU", 1.40 },
    .{ "ZN", 1.39 },
    .{ "MN", 1.19 },
    // Other metals
    .{ "CD", 1.58 },
    .{ "HG", 1.55 },
    .{ "PB", 2.02 },
    .{ "AS", 1.85 },
});

/// Guess van der Waals radius from element symbol
/// Returns null if element is not recognized
/// Lookup is case-insensitive and whitespace is trimmed
pub fn guessRadius(element: []const u8) ?f64 {
    if (element.len == 0) return null;

    // Normalize: trim whitespace
    const trimmed = std.mem.trim(u8, element, " ");
    if (trimmed.len == 0 or trimmed.len > 2) return null;

    // Convert to uppercase for map lookup
    var upper: [2]u8 = undefined;
    for (trimmed, 0..) |c, i| {
        upper[i] = std.ascii.toUpper(c);
    }

    return element_radii_map.get(upper[0..trimmed.len]);
}

/// Extract element symbol from PDB atom name
/// PDB atom names are 4 characters with specific conventions:
/// - Columns 13-16 in PDB format
/// - For most atoms: first char is space, element is inferred from name
/// - For 2-char elements (Fe, Ca, etc.): both chars are the element at start
///
/// Key rule: If original input starts with a space, element is single-character.
/// Only check for 2-char elements if input starts with non-space.
///
/// **Lifetime note**: Returns a slice into the input `atom_name`.
/// The returned slice is only valid as long as `atom_name` is valid.
///
/// Examples:
/// " CA " -> "C" (carbon alpha, not calcium - leading space!)
/// " NE2" -> "N" (nitrogen)
/// "FE  " -> "FE" (iron - no leading space)
/// "CA  " -> "CA" (calcium - no leading space, 2-char element)
/// " OG1" -> "O" (oxygen)
/// "SE  " -> "SE" (selenium)
pub fn extractElement(atom_name: []const u8) []const u8 {
    if (atom_name.len == 0) return "";

    const trimmed = std.mem.trim(u8, atom_name, " ");
    if (trimmed.len == 0) return "";

    // PDB convention: if original input starts with space, element is single-char
    // Only check for 2-char elements if input starts with a letter (no leading space)
    const has_leading_space = atom_name[0] == ' ';

    if (!has_leading_space and trimmed.len >= 2) {
        const first_two = trimmed[0..2];
        // Check if it's a known 2-character element using O(1) map lookup
        var upper: [2]u8 = .{ std.ascii.toUpper(first_two[0]), std.ascii.toUpper(first_two[1]) };
        if (element_radii_map.get(&upper) != null) {
            return first_two;
        }
    }

    // Default: first character is the element (C, N, O, H, S, P)
    return trimmed[0..1];
}

/// Get radius by guessing from atom name (extracts element first)
pub fn guessRadiusFromAtomName(atom_name: []const u8) ?f64 {
    const element = extractElement(atom_name);
    return guessRadius(element);
}

// =============================================================================
// Tests
// =============================================================================

test "ClassifierType fromString" {
    try std.testing.expectEqual(ClassifierType.naccess, ClassifierType.fromString("naccess").?);
    try std.testing.expectEqual(ClassifierType.protor, ClassifierType.fromString("protor").?);
    try std.testing.expectEqual(ClassifierType.oons, ClassifierType.fromString("oons").?);
    try std.testing.expectEqual(ClassifierType.ccd, ClassifierType.fromString("ccd").?);
    try std.testing.expectEqual(ClassifierType.ccd, ClassifierType.fromString("CCD").?);
    try std.testing.expectEqual(@as(?ClassifierType, null), ClassifierType.fromString("invalid"));
}

test "ClassifierType name" {
    try std.testing.expectEqualStrings("NACCESS", ClassifierType.naccess.name());
    try std.testing.expectEqualStrings("ProtOr", ClassifierType.protor.name());
    try std.testing.expectEqualStrings("OONS", ClassifierType.oons.name());
    try std.testing.expectEqualStrings("CCD", ClassifierType.ccd.name());
}

test "AtomClass fromString" {
    try std.testing.expectEqual(AtomClass.polar, AtomClass.fromString("polar"));
    try std.testing.expectEqual(AtomClass.apolar, AtomClass.fromString("apolar"));
    try std.testing.expectEqual(AtomClass.unknown, AtomClass.fromString("invalid"));
    try std.testing.expectEqual(AtomClass.unknown, AtomClass.fromString(""));
}

test "AtomClass toString" {
    try std.testing.expectEqualStrings("polar", AtomClass.polar.toString());
    try std.testing.expectEqualStrings("apolar", AtomClass.apolar.toString());
    try std.testing.expectEqualStrings("unknown", AtomClass.unknown.toString());
}

test "AtomKey init and accessors" {
    const key = AtomKey.init("ALA", "CA");
    try std.testing.expectEqualStrings("ALA", key.residueName());
    try std.testing.expectEqualStrings("CA", key.atomName());
}

test "AtomKey equality" {
    const key1 = AtomKey.init("ALA", "CA");
    const key2 = AtomKey.init("ALA", "CA");
    const key3 = AtomKey.init("GLY", "CA");
    const key4 = AtomKey.init("ALA", "CB");

    try std.testing.expect(key1.eql(key2));
    try std.testing.expect(!key1.eql(key3));
    try std.testing.expect(!key1.eql(key4));
}

test "AtomKey hash consistency" {
    const key1 = AtomKey.init("ALA", "CA");
    const key2 = AtomKey.init("ALA", "CA");
    const key3 = AtomKey.init("GLY", "CA");

    try std.testing.expectEqual(key1.hash(), key2.hash());
    try std.testing.expect(key1.hash() != key3.hash());
}

test "Classifier init and deinit" {
    const allocator = std.testing.allocator;
    var classifier = try Classifier.init(allocator, "test");
    defer classifier.deinit();

    try std.testing.expectEqualStrings("test", classifier.name);
    try std.testing.expectEqual(@as(usize, 0), classifier.count());
}

test "Classifier addAtom and getRadius" {
    const allocator = std.testing.allocator;
    var classifier = try Classifier.init(allocator, "test");
    defer classifier.deinit();

    try classifier.addAtom("ALA", "CA", 1.87, .apolar);
    try classifier.addAtom("ALA", "CB", 1.87, .apolar);
    try classifier.addAtom("ANY", "C", 1.76, .apolar);
    try classifier.addAtom("ANY", "O", 1.40, .polar);

    // Exact match
    try std.testing.expectEqual(@as(?f64, 1.87), classifier.getRadius("ALA", "CA"));
    try std.testing.expectEqual(@as(?f64, 1.87), classifier.getRadius("ALA", "CB"));

    // ANY fallback
    try std.testing.expectEqual(@as(?f64, 1.76), classifier.getRadius("ALA", "C"));
    try std.testing.expectEqual(@as(?f64, 1.76), classifier.getRadius("GLY", "C"));
    try std.testing.expectEqual(@as(?f64, 1.40), classifier.getRadius("ALA", "O"));

    // Not found
    try std.testing.expectEqual(@as(?f64, null), classifier.getRadius("ALA", "XX"));
}

test "Classifier getClass with fallback" {
    const allocator = std.testing.allocator;
    var classifier = try Classifier.init(allocator, "test");
    defer classifier.deinit();

    try classifier.addAtom("ALA", "CA", 1.87, .apolar);
    try classifier.addAtom("ANY", "N", 1.55, .polar);

    // Exact match
    try std.testing.expectEqual(AtomClass.apolar, classifier.getClass("ALA", "CA"));

    // ANY fallback
    try std.testing.expectEqual(AtomClass.polar, classifier.getClass("ALA", "N"));
    try std.testing.expectEqual(AtomClass.polar, classifier.getClass("GLY", "N"));

    // Not found
    try std.testing.expectEqual(AtomClass.unknown, classifier.getClass("ALA", "XX"));
}

test "Classifier getProperties" {
    const allocator = std.testing.allocator;
    var classifier = try Classifier.init(allocator, "test");
    defer classifier.deinit();

    try classifier.addAtom("ALA", "CA", 1.87, .apolar);
    try classifier.addAtom("ANY", "O", 1.40, .polar);

    // Exact match
    const props1 = classifier.getProperties("ALA", "CA");
    try std.testing.expect(props1 != null);
    try std.testing.expectEqual(@as(f64, 1.87), props1.?.radius);
    try std.testing.expectEqual(AtomClass.apolar, props1.?.class);

    // ANY fallback
    const props2 = classifier.getProperties("GLY", "O");
    try std.testing.expect(props2 != null);
    try std.testing.expectEqual(@as(f64, 1.40), props2.?.radius);
    try std.testing.expectEqual(AtomClass.polar, props2.?.class);

    // Not found
    try std.testing.expectEqual(@as(?AtomProperties, null), classifier.getProperties("ALA", "XX"));
}

test "Classifier residue-specific override" {
    const allocator = std.testing.allocator;
    var classifier = try Classifier.init(allocator, "test");
    defer classifier.deinit();

    // Generic ANY entry
    try classifier.addAtom("ANY", "CB", 1.70, .apolar);
    // Residue-specific override
    try classifier.addAtom("CYS", "CB", 1.80, .apolar);

    // CYS should use specific value
    try std.testing.expectEqual(@as(?f64, 1.80), classifier.getRadius("CYS", "CB"));

    // Other residues should use ANY
    try std.testing.expectEqual(@as(?f64, 1.70), classifier.getRadius("ALA", "CB"));
    try std.testing.expectEqual(@as(?f64, 1.70), classifier.getRadius("GLY", "CB"));
}

test "AtomKey truncation behavior" {
    // Residue names longer than 5 chars are truncated (mmCIF comp_id support)
    // Atom names longer than 4 chars are truncated
    const key = AtomKey.init("TOOLONG", "VERYLONGNAME");
    try std.testing.expectEqualStrings("TOOLO", key.residueName());
    try std.testing.expectEqualStrings("VERY", key.atomName());
}

test "AtomKey hash distinguishes similar keys" {
    // Ensure separator prevents collision between ("A", "BCD") and ("AB", "CD")
    const key1 = AtomKey.init("A", "BCD");
    const key2 = AtomKey.init("AB", "CD");
    try std.testing.expect(key1.hash() != key2.hash());

    // Also test ("ABC", "D") vs ("A", "BCD")
    const key3 = AtomKey.init("ABC", "D");
    try std.testing.expect(key1.hash() != key3.hash());
    try std.testing.expect(key2.hash() != key3.hash());
}

// =============================================================================
// Element-Based Radius Tests
// =============================================================================

test "guessRadius common elements" {
    // Common biological elements
    try std.testing.expectEqual(@as(?f64, 1.10), guessRadius("H"));
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadius("C"));
    try std.testing.expectEqual(@as(?f64, 1.55), guessRadius("N"));
    try std.testing.expectEqual(@as(?f64, 1.52), guessRadius("O"));
    try std.testing.expectEqual(@as(?f64, 1.80), guessRadius("P"));
    try std.testing.expectEqual(@as(?f64, 1.80), guessRadius("S"));
    try std.testing.expectEqual(@as(?f64, 1.90), guessRadius("SE"));
}

test "guessRadius case insensitive" {
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadius("c"));
    try std.testing.expectEqual(@as(?f64, 1.55), guessRadius("n"));
    try std.testing.expectEqual(@as(?f64, 1.90), guessRadius("se"));
    try std.testing.expectEqual(@as(?f64, 1.90), guessRadius("Se"));
}

test "guessRadius with whitespace" {
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadius(" C"));
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadius("C "));
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadius(" C "));
    try std.testing.expectEqual(@as(?f64, 1.39), guessRadius(" ZN "));
}

test "guessRadius metals" {
    try std.testing.expectEqual(@as(?f64, 1.26), guessRadius("FE"));
    try std.testing.expectEqual(@as(?f64, 1.39), guessRadius("ZN"));
    try std.testing.expectEqual(@as(?f64, 1.40), guessRadius("CU"));
    try std.testing.expectEqual(@as(?f64, 1.73), guessRadius("MG"));
    try std.testing.expectEqual(@as(?f64, 2.31), guessRadius("CA"));
}

test "guessRadius unknown element" {
    try std.testing.expectEqual(@as(?f64, null), guessRadius("XX"));
    try std.testing.expectEqual(@as(?f64, null), guessRadius(""));
    try std.testing.expectEqual(@as(?f64, null), guessRadius("   "));
}

test "extractElement single char elements" {
    // Standard amino acid atoms - first char after space is element
    try std.testing.expectEqualStrings("C", extractElement(" CA "));
    try std.testing.expectEqualStrings("C", extractElement(" CB "));
    try std.testing.expectEqualStrings("N", extractElement(" N  "));
    try std.testing.expectEqualStrings("N", extractElement(" NE2"));
    try std.testing.expectEqualStrings("O", extractElement(" O  "));
    try std.testing.expectEqualStrings("O", extractElement(" OG1"));
    try std.testing.expectEqualStrings("H", extractElement(" H  "));
    try std.testing.expectEqualStrings("S", extractElement(" SG "));
}

test "extractElement two char elements" {
    // Metal ions and selenium
    try std.testing.expectEqualStrings("FE", extractElement("FE  "));
    try std.testing.expectEqualStrings("ZN", extractElement("ZN  "));
    try std.testing.expectEqualStrings("CA", extractElement("CA  "));
    try std.testing.expectEqualStrings("MG", extractElement("MG  "));
    try std.testing.expectEqualStrings("SE", extractElement("SE  "));
    try std.testing.expectEqualStrings("CU", extractElement("CU  "));
}

test "extractElement trimmed input" {
    // Without leading space, 2-char elements are recognized
    try std.testing.expectEqualStrings("CA", extractElement("CA")); // Calcium (no leading space)
    try std.testing.expectEqualStrings("N", extractElement("N")); // Nitrogen (1-char)
    try std.testing.expectEqualStrings("FE", extractElement("FE")); // Iron (2-char)
    try std.testing.expectEqualStrings("C", extractElement("C")); // Carbon (1-char)
}

test "extractElement empty input" {
    try std.testing.expectEqualStrings("", extractElement(""));
    try std.testing.expectEqualStrings("", extractElement("    "));
}

test "guessRadiusFromAtomName" {
    // Standard atoms
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadiusFromAtomName(" CA "));
    try std.testing.expectEqual(@as(?f64, 1.55), guessRadiusFromAtomName(" N  "));
    try std.testing.expectEqual(@as(?f64, 1.52), guessRadiusFromAtomName(" O  "));
    try std.testing.expectEqual(@as(?f64, 1.80), guessRadiusFromAtomName(" SG "));

    // Metals
    try std.testing.expectEqual(@as(?f64, 1.26), guessRadiusFromAtomName("FE  "));
    try std.testing.expectEqual(@as(?f64, 1.39), guessRadiusFromAtomName("ZN  "));
    try std.testing.expectEqual(@as(?f64, 1.90), guessRadiusFromAtomName("SE  "));
}

test "guessRadiusFromAtomicNumber common elements" {
    // Common biological elements
    try std.testing.expectEqual(@as(?f64, 1.10), guessRadiusFromAtomicNumber(1)); // H
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadiusFromAtomicNumber(6)); // C
    try std.testing.expectEqual(@as(?f64, 1.55), guessRadiusFromAtomicNumber(7)); // N
    try std.testing.expectEqual(@as(?f64, 1.52), guessRadiusFromAtomicNumber(8)); // O
    try std.testing.expectEqual(@as(?f64, 1.80), guessRadiusFromAtomicNumber(15)); // P
    try std.testing.expectEqual(@as(?f64, 1.80), guessRadiusFromAtomicNumber(16)); // S
    try std.testing.expectEqual(@as(?f64, 1.90), guessRadiusFromAtomicNumber(34)); // Se
}

test "guessRadiusFromAtomicNumber metals" {
    try std.testing.expectEqual(@as(?f64, 1.26), guessRadiusFromAtomicNumber(26)); // Fe
    try std.testing.expectEqual(@as(?f64, 1.39), guessRadiusFromAtomicNumber(30)); // Zn
    try std.testing.expectEqual(@as(?f64, 1.40), guessRadiusFromAtomicNumber(29)); // Cu
    try std.testing.expectEqual(@as(?f64, 1.73), guessRadiusFromAtomicNumber(12)); // Mg
    try std.testing.expectEqual(@as(?f64, 2.31), guessRadiusFromAtomicNumber(20)); // Ca
}

test "guessRadiusFromAtomicNumber distinguishes CA (Carbon) from Ca (Calcium)" {
    // This is the key test: atomic numbers unambiguously identify elements
    // Carbon (atomic number 6) vs Calcium (atomic number 20)
    try std.testing.expectEqual(@as(?f64, 1.70), guessRadiusFromAtomicNumber(6)); // Carbon
    try std.testing.expectEqual(@as(?f64, 2.31), guessRadiusFromAtomicNumber(20)); // Calcium
}

test "guessRadiusFromAtomicNumber unknown" {
    try std.testing.expectEqual(@as(?f64, null), guessRadiusFromAtomicNumber(0)); // Invalid
    try std.testing.expectEqual(@as(?f64, null), guessRadiusFromAtomicNumber(2)); // He (not in our map)
    try std.testing.expectEqual(@as(?f64, null), guessRadiusFromAtomicNumber(118)); // Oganesson
}

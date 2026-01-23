//! Atom radius classifier for SASA calculation.
//!
//! This module provides classification of atoms based on residue name and atom name,
//! returning appropriate van der Waals radii and polarity classes.
//!
//! The classifier follows FreeSASA's approach:
//! 1. Look up by (residue_name, atom_name)
//! 2. Fall back to ("ANY", atom_name)
//! 3. Fall back to element-based radius guessing

const std = @import("std");
const Allocator = std.mem.Allocator;

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
    residue: [4]u8, // PDB residue names are max 3-4 chars
    residue_len: u8,
    atom: [4]u8, // PDB atom names are max 4 chars
    atom_len: u8,

    pub fn init(residue: []const u8, atom: []const u8) AtomKey {
        var key = AtomKey{
            .residue = [_]u8{0} ** 4,
            .residue_len = @intCast(@min(residue.len, 4)),
            .atom = [_]u8{0} ** 4,
            .atom_len = @intCast(@min(atom.len, 4)),
        };
        @memcpy(key.residue[0..key.residue_len], residue[0..key.residue_len]);
        @memcpy(key.atom[0..key.atom_len], atom[0..key.atom_len]);
        return key;
    }

    pub fn residueName(self: AtomKey) []const u8 {
        return self.residue[0..self.residue_len];
    }

    pub fn atomName(self: AtomKey) []const u8 {
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
// Tests
// =============================================================================

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
    // Names longer than 4 chars are truncated
    const key = AtomKey.init("TOOLONG", "VERYLONGNAME");
    try std.testing.expectEqualStrings("TOOL", key.residueName());
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

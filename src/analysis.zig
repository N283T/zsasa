//! Analysis module for SASA aggregation and classification.
//!
//! Provides functions for:
//! - Per-residue SASA aggregation
//! - RSA (Relative Solvent Accessibility) calculation
//! - Polar/nonpolar classification

const std = @import("std");
const classifier = @import("classifier.zig");
const types = @import("types.zig");

/// Polarity class that a classifier gives an atom
pub const AtomClass = classifier.AtomClass;

/// Maximum SASA values for standard amino acids (in Å²).
/// Values from Tien et al. (2013) "Maximum allowed solvent accessibilities
/// of residues in proteins" - empirical values.
/// These represent the theoretical maximum exposure for each residue type.
pub const MaxSASA = struct {
    // Standard 20 amino acids
    pub const ALA: f64 = 129.0;
    pub const ARG: f64 = 274.0;
    pub const ASN: f64 = 195.0;
    pub const ASP: f64 = 193.0;
    pub const CYS: f64 = 167.0;
    pub const GLN: f64 = 225.0;
    pub const GLU: f64 = 223.0;
    pub const GLY: f64 = 104.0;
    pub const HIS: f64 = 224.0;
    pub const ILE: f64 = 197.0;
    pub const LEU: f64 = 201.0;
    pub const LYS: f64 = 236.0;
    pub const MET: f64 = 224.0;
    pub const PHE: f64 = 240.0;
    pub const PRO: f64 = 159.0;
    pub const SER: f64 = 155.0;
    pub const THR: f64 = 172.0;
    pub const TRP: f64 = 285.0;
    pub const TYR: f64 = 263.0;
    pub const VAL: f64 = 174.0;

    /// Get MaxSASA for a given 3-letter residue code.
    /// Returns null for unknown residues.
    pub fn get(residue_name: []const u8) ?f64 {
        if (std.mem.eql(u8, residue_name, "ALA")) return ALA;
        if (std.mem.eql(u8, residue_name, "ARG")) return ARG;
        if (std.mem.eql(u8, residue_name, "ASN")) return ASN;
        if (std.mem.eql(u8, residue_name, "ASP")) return ASP;
        if (std.mem.eql(u8, residue_name, "CYS")) return CYS;
        if (std.mem.eql(u8, residue_name, "GLN")) return GLN;
        if (std.mem.eql(u8, residue_name, "GLU")) return GLU;
        if (std.mem.eql(u8, residue_name, "GLY")) return GLY;
        if (std.mem.eql(u8, residue_name, "HIS")) return HIS;
        if (std.mem.eql(u8, residue_name, "ILE")) return ILE;
        if (std.mem.eql(u8, residue_name, "LEU")) return LEU;
        if (std.mem.eql(u8, residue_name, "LYS")) return LYS;
        if (std.mem.eql(u8, residue_name, "MET")) return MET;
        if (std.mem.eql(u8, residue_name, "PHE")) return PHE;
        if (std.mem.eql(u8, residue_name, "PRO")) return PRO;
        if (std.mem.eql(u8, residue_name, "SER")) return SER;
        if (std.mem.eql(u8, residue_name, "THR")) return THR;
        if (std.mem.eql(u8, residue_name, "TRP")) return TRP;
        if (std.mem.eql(u8, residue_name, "TYR")) return TYR;
        if (std.mem.eql(u8, residue_name, "VAL")) return VAL;
        return null;
    }
};

/// Residue classification for polar/nonpolar analysis
pub const ResidueClass = enum {
    polar,
    nonpolar,
    unknown,

    /// Classify a residue by its 3-letter code.
    /// Polar: charged (D, E, H, K, R) or H-bonding (N, Q, S, T, Y)
    /// Nonpolar: hydrophobic (A, C, F, G, I, L, M, P, V, W)
    pub fn fromResidueName(residue_name: []const u8) ResidueClass {
        // Polar residues (charged or H-bonding capable)
        if (std.mem.eql(u8, residue_name, "ARG")) return .polar; // R - charged
        if (std.mem.eql(u8, residue_name, "ASN")) return .polar; // N - H-bonding
        if (std.mem.eql(u8, residue_name, "ASP")) return .polar; // D - charged
        if (std.mem.eql(u8, residue_name, "GLN")) return .polar; // Q - H-bonding
        if (std.mem.eql(u8, residue_name, "GLU")) return .polar; // E - charged
        if (std.mem.eql(u8, residue_name, "HIS")) return .polar; // H - charged
        if (std.mem.eql(u8, residue_name, "LYS")) return .polar; // K - charged
        if (std.mem.eql(u8, residue_name, "SER")) return .polar; // S - H-bonding
        if (std.mem.eql(u8, residue_name, "THR")) return .polar; // T - H-bonding
        if (std.mem.eql(u8, residue_name, "TYR")) return .polar; // Y - H-bonding

        // Nonpolar residues (hydrophobic)
        if (std.mem.eql(u8, residue_name, "ALA")) return .nonpolar; // A
        if (std.mem.eql(u8, residue_name, "CYS")) return .nonpolar; // C
        if (std.mem.eql(u8, residue_name, "PHE")) return .nonpolar; // F
        if (std.mem.eql(u8, residue_name, "GLY")) return .nonpolar; // G
        if (std.mem.eql(u8, residue_name, "ILE")) return .nonpolar; // I
        if (std.mem.eql(u8, residue_name, "LEU")) return .nonpolar; // L
        if (std.mem.eql(u8, residue_name, "MET")) return .nonpolar; // M
        if (std.mem.eql(u8, residue_name, "PRO")) return .nonpolar; // P
        if (std.mem.eql(u8, residue_name, "TRP")) return .nonpolar; // W
        if (std.mem.eql(u8, residue_name, "VAL")) return .nonpolar; // V

        return .unknown;
    }
};

/// Summary of polar/nonpolar SASA
pub const PolarSummary = struct {
    polar_sasa: f64,
    nonpolar_sasa: f64,
    unknown_sasa: f64,
    polar_residue_count: usize,
    nonpolar_residue_count: usize,
    unknown_residue_count: usize,

    pub fn polarFraction(self: PolarSummary) f64 {
        const total = self.polar_sasa + self.nonpolar_sasa;
        if (total > 0) {
            return self.polar_sasa / total;
        }
        return 0;
    }

    pub fn nonpolarFraction(self: PolarSummary) f64 {
        const total = self.polar_sasa + self.nonpolar_sasa;
        if (total > 0) {
            return self.nonpolar_sasa / total;
        }
        return 0;
    }
};

/// Calculate polar/nonpolar SASA summary from per-residue data
pub fn calculatePolarSummary(residues: []const ResidueSasa) PolarSummary {
    var summary = PolarSummary{
        .polar_sasa = 0,
        .nonpolar_sasa = 0,
        .unknown_sasa = 0,
        .polar_residue_count = 0,
        .nonpolar_residue_count = 0,
        .unknown_residue_count = 0,
    };

    for (residues) |res| {
        switch (ResidueClass.fromResidueName(res.residue_name.slice())) {
            .polar => {
                summary.polar_sasa += res.sasa;
                summary.polar_residue_count += 1;
            },
            .nonpolar => {
                summary.nonpolar_sasa += res.sasa;
                summary.nonpolar_residue_count += 1;
            },
            .unknown => {
                summary.unknown_sasa += res.sasa;
                summary.unknown_residue_count += 1;
            },
        }
    }

    return summary;
}

/// Print polar/nonpolar SASA summary.
/// Note: Percentages are calculated excluding unknown residues.
pub fn printPolarSummary(summary: PolarSummary) void {
    std.debug.print("\nPolar/Nonpolar SASA:\n", .{});
    std.debug.print("  Polar:    {d:>10.2} Å² ({d:>5.1}%) - {d} residues\n", .{
        summary.polar_sasa,
        summary.polarFraction() * 100,
        summary.polar_residue_count,
    });
    std.debug.print("  Nonpolar: {d:>10.2} Å² ({d:>5.1}%) - {d} residues\n", .{
        summary.nonpolar_sasa,
        summary.nonpolarFraction() * 100,
        summary.nonpolar_residue_count,
    });
    if (summary.unknown_sasa > 0) {
        std.debug.print("  Unknown:  {d:>10.2} Å² - {d} residues (excluded from %)\n", .{
            summary.unknown_sasa,
            summary.unknown_residue_count,
        });
    }
}

/// Whether an atom counts as polar in the polar/non-polar partition by atom:
/// the `Non-polar` and `All polar` columns of the RSA file and the atom
/// summary of `--polar`.
///
/// `class` is the class that the active classifier gives the atom, the same
/// classifier that set its radius. Classifiers disagree on some atoms
/// (NACCESS classes sulfur as apolar, OONS classes carbonyl carbon as polar),
/// so the partition follows the classifier. An atom that the classifier does
/// not class (`.unknown`: hydrogens, ligands outside its tables), and every
/// atom when no classifier ran, falls back on the element: N, O, P and S are
/// polar, everything else is apolar. The element is the atom's entry in the
/// input's element column or, without such a column, the first letter of the
/// atom name; an atom without a name is apolar.
pub fn isPolarAtom(class: AtomClass, atom_name: ?[]const u8, element: ?u8) bool {
    switch (class) {
        .polar => return true,
        .apolar => return false,
        .unknown => {},
    }
    const name = atom_name orelse return false;
    if (element) |atomic_number| {
        return atomic_number == 7 or atomic_number == 8 or atomic_number == 15 or atomic_number == 16;
    }
    const trimmed = std.mem.trim(u8, name, " ");
    if (trimmed.len == 0) return false;
    const c = std.ascii.toUpper(trimmed[0]);
    return c == 'N' or c == 'O' or c == 'P' or c == 'S';
}

/// Class of atom `i` in `atom_classes`, or `.unknown` when no classifier ran.
pub fn atomClassAt(atom_classes: ?[]const AtomClass, i: usize) AtomClass {
    return if (atom_classes) |classes| classes[i] else .unknown;
}

/// Polar/non-polar SASA by atom, as partitioned by `isPolarAtom`.
pub const AtomPolarSummary = struct {
    polar_sasa: f64 = 0,
    apolar_sasa: f64 = 0,
    polar_atom_count: usize = 0,
    apolar_atom_count: usize = 0,
    /// Atoms without a class from the classifier, classed by their element
    fallback_atom_count: usize = 0,

    pub fn polarFraction(self: AtomPolarSummary) f64 {
        const total = self.polar_sasa + self.apolar_sasa;
        return if (total > 0) self.polar_sasa / total else 0;
    }

    pub fn apolarFraction(self: AtomPolarSummary) f64 {
        const total = self.polar_sasa + self.apolar_sasa;
        return if (total > 0) self.apolar_sasa / total else 0;
    }
};

/// Sum atom areas by polarity. `atom_classes` holds the class of every atom
/// from the active classifier, or is null when no classifier ran.
pub fn calculateAtomPolarSummary(
    input: types.AtomInput,
    atom_areas: []const f64,
    atom_classes: ?[]const AtomClass,
) !AtomPolarSummary {
    const n = input.atomCount();
    if (atom_areas.len != n) return error.LengthMismatch;
    if (atom_classes) |classes| {
        if (classes.len != n) return error.LengthMismatch;
    }

    var summary = AtomPolarSummary{};
    for (atom_areas, 0..) |area, i| {
        const class = atomClassAt(atom_classes, i);
        if (class == .unknown) summary.fallback_atom_count += 1;
        const polar = isPolarAtom(
            class,
            if (input.atom_name) |names| names[i].slice() else null,
            if (input.element) |elements| elements[i] else null,
        );
        if (polar) {
            summary.polar_sasa += area;
            summary.polar_atom_count += 1;
        } else {
            summary.apolar_sasa += area;
            summary.apolar_atom_count += 1;
        }
    }
    return summary;
}

/// Print the polar/non-polar SASA by atom class below the summary by residue
/// type. The areas are those of the `TOTAL` row of the RSA file.
pub fn printAtomPolarSummary(summary: AtomPolarSummary, classifier_name: []const u8) void {
    std.debug.print("\nPolar/Nonpolar SASA by atom class (classifier: {s}):\n", .{classifier_name});
    std.debug.print("  Polar:    {d:>10.2} Å² ({d:>5.1}%) - {d} atoms\n", .{
        summary.polar_sasa,
        summary.polarFraction() * 100,
        summary.polar_atom_count,
    });
    std.debug.print("  Nonpolar: {d:>10.2} Å² ({d:>5.1}%) - {d} atoms\n", .{
        summary.apolar_sasa,
        summary.apolarFraction() * 100,
        summary.apolar_atom_count,
    });
    if (summary.fallback_atom_count > 0) {
        std.debug.print("  ({d} atoms without a class from the classifier are classed by element)\n", .{
            summary.fallback_atom_count,
        });
    }
}

/// Per-residue SASA data
pub const ResidueSasa = struct {
    chain_id: types.FixedString4,
    chain_id_full: ?[]const u8 = null,
    residue_name: types.FixedString5,
    residue_num: i32,
    insertion_code: types.FixedString4,
    sasa: f64,
    atom_count: usize,
    /// Relative Solvent Accessibility (0.0-1.0+), null if MaxSASA unknown
    rsa: ?f64 = null,

    /// Calculate RSA from SASA and residue name.
    /// RSA values > 1.0 are possible for exposed terminal residues.
    pub fn calculateRsa(self: *ResidueSasa) void {
        if (MaxSASA.get(self.residue_name.slice())) |max_sasa| {
            if (max_sasa > 0) {
                self.rsa = self.sasa / max_sasa;
            } else {
                self.rsa = null;
            }
        } else {
            self.rsa = null;
        }
    }

    pub fn chainLabel(self: ResidueSasa) []const u8 {
        return self.chain_id_full orelse self.chain_id.slice();
    }
};

/// Result of per-residue aggregation
pub const ResidueResult = struct {
    residues: []ResidueSasa,
    allocator: std.mem.Allocator,

    pub fn deinit(self: *ResidueResult) void {
        for (self.residues) |res| {
            if (res.chain_id_full) |chain_id_full| {
                self.allocator.free(chain_id_full);
            }
        }
        self.allocator.free(self.residues);
    }
};

fn freeResidueFullChainLabels(allocator: std.mem.Allocator, residues: []ResidueSasa) void {
    for (residues) |res| {
        if (res.chain_id_full) |chain_id_full| {
            allocator.free(chain_id_full);
        }
    }
}

/// Residue identity, shared by every output that reports residues: the
/// `--per-residue` / `--rsa` table (`aggregateByResidue`), the RSA file
/// (`json_writer.sasaResultToRsa`), and the JSONL residue map
/// (`json_writer.buildResidueMap`) with the BSA residue arrays built from it.
///
/// Two atoms belong to the same residue when they are adjacent in the input
/// and agree in all of
///
/// - the chain ID (the full ID from `chain_id_full` where the parser kept
///   one, otherwise `chain_id`),
/// - the residue number,
/// - the insertion code, and
/// - the residue name.
///
/// A residue is thus a maximal run of consecutive atoms with one identity,
/// reported in input order. Two consequences are deliberate:
///
/// - Atoms of one residue that are not contiguous in the input give one
///   entry per run, with the same labels. The JSONL residue map describes a
///   residue as an atom range (`residue_atom_start`, `residue_atom_count`),
///   which cannot hold a residue that is scattered over the file, and the
///   other outputs follow it so that all of them report the same entries.
///   FreeSASA starts a new residue the same way.
/// - A multi-model file read with all models superimposed (the default)
///   repeats every residue once per model. Each repetition is its own
///   entry, with the area that this copy has inside the superimposed
///   structure; the copies are not summed. Only when a model ends with the
///   identity that the next one begins with (a model of a single residue)
///   are the two runs adjacent and reported as one entry.
pub const ResidueIdentity = struct {
    chain_ids: []const types.FixedString4,
    chain_ids_full: ?[]const []const u8,
    residue_names: []const types.FixedString5,
    residue_nums: []const i32,
    insertion_codes: []const types.FixedString4,

    /// Atom index range `[start, end)` of one residue.
    pub const Range = struct {
        start: usize,
        end: usize,

        pub fn atomCount(self: Range) usize {
            return self.end - self.start;
        }
    };

    /// Iterates over the residues of an input in input order.
    pub const Iterator = struct {
        identity: ResidueIdentity,
        next_atom: usize = 0,

        pub fn next(self: *Iterator) ?Range {
            const start = self.next_atom;
            if (start >= self.identity.atomCount()) return null;
            var end = start + 1;
            while (end < self.identity.atomCount() and self.identity.sameResidue(start, end)) : (end += 1) {}
            self.next_atom = end;
            return .{ .start = start, .end = end };
        }
    };

    /// Borrows the identity columns of `input`, which must outlive the result.
    pub fn init(input: types.AtomInput) !ResidueIdentity {
        return .{
            .chain_ids = input.chain_id orelse return error.MissingChainInfo,
            .chain_ids_full = input.chain_id_full,
            .residue_names = input.residue orelse return error.MissingResidueInfo,
            .residue_nums = input.residue_num orelse return error.MissingResidueNumInfo,
            .insertion_codes = input.insertion_code orelse return error.MissingInsertionCodeInfo,
        };
    }

    pub fn atomCount(self: ResidueIdentity) usize {
        return self.chain_ids.len;
    }

    /// Chain ID of an atom as it is written to output: the full ID where the
    /// parser kept one, otherwise the (at most four-character) `chain_id`.
    pub fn chainLabel(self: ResidueIdentity, atom: usize) []const u8 {
        return if (self.chain_ids_full) |full| full[atom] else self.chain_ids[atom].slice();
    }

    pub fn sameChain(self: ResidueIdentity, a: usize, b: usize) bool {
        return std.mem.eql(u8, self.chainLabel(a), self.chainLabel(b));
    }

    /// Whether atoms `a` and `b` have the same residue identity. Atoms of one
    /// residue are also adjacent; see the type's documentation.
    pub fn sameResidue(self: ResidueIdentity, a: usize, b: usize) bool {
        return self.residue_nums[a] == self.residue_nums[b] and
            self.sameChain(a, b) and
            std.mem.eql(u8, self.insertion_codes[a].slice(), self.insertion_codes[b].slice()) and
            std.mem.eql(u8, self.residue_names[a].slice(), self.residue_names[b].slice());
    }

    pub fn residues(self: ResidueIdentity) Iterator {
        return .{ .identity = self };
    }

    pub fn residueCount(self: ResidueIdentity) usize {
        var count: usize = 0;
        var it = self.residues();
        while (it.next()) |_| count += 1;
        return count;
    }
};

/// Aggregate atom SASA values to per-residue SASA.
/// Atoms are grouped into residues by `ResidueIdentity`.
pub fn aggregateByResidue(
    allocator: std.mem.Allocator,
    input: types.AtomInput,
    atom_areas: []const f64,
) !ResidueResult {
    // Check if we have the required residue info
    const identity = try ResidueIdentity.init(input);

    const n = input.atomCount();
    if (atom_areas.len != n) {
        return error.LengthMismatch;
    }

    var residue_list = std.ArrayListUnmanaged(ResidueSasa).empty;
    errdefer freeResidueFullChainLabels(allocator, residue_list.items);
    defer residue_list.deinit(allocator);
    try residue_list.ensureTotalCapacity(allocator, identity.residueCount());

    var it = identity.residues();
    while (it.next()) |range| {
        var sasa = atom_areas[range.start];
        for (atom_areas[range.start + 1 .. range.end]) |area| sasa += area;

        var residue = ResidueSasa{
            .chain_id = identity.chain_ids[range.start],
            .chain_id_full = if (identity.chain_ids_full) |full|
                try allocator.dupe(u8, full[range.start])
            else
                null,
            .residue_name = identity.residue_names[range.start],
            .residue_num = identity.residue_nums[range.start],
            .insertion_code = identity.insertion_codes[range.start],
            .sasa = sasa,
            .atom_count = range.atomCount(),
        };
        residue.calculateRsa();
        residue_list.appendAssumeCapacity(residue);
    }

    return ResidueResult{
        .residues = try residue_list.toOwnedSlice(allocator),
        .allocator = allocator,
    };
}

/// Print per-residue SASA results
pub fn printResidueResults(residues: []const ResidueSasa) void {
    std.debug.print("\nPer-residue SASA:\n", .{});
    std.debug.print("{s:>5} {s:>4} {s:>6} {s:>10} {s:>6}\n", .{
        "Chain", "Res", "Num", "SASA", "Atoms",
    });
    std.debug.print("{s:->5} {s:->4} {s:->6} {s:->10} {s:->6}\n", .{
        "", "", "", "", "",
    });

    for (residues) |res| {
        // Format residue number as string to avoid Zig's "+" prefix for positive integers
        var num_buf: [16]u8 = undefined;
        const num_str = std.fmt.bufPrint(&num_buf, "{d}", .{res.residue_num}) catch "?";

        if (res.insertion_code.len > 0) {
            std.debug.print("{s:>5} {s:>4} {s:>5}{s:<1} {d:>10.2} {d:>6}\n", .{
                res.chainLabel(),
                res.residue_name.slice(),
                num_str,
                res.insertion_code.slice(),
                res.sasa,
                res.atom_count,
            });
        } else {
            std.debug.print("{s:>5} {s:>4} {s:>6} {d:>10.2} {d:>6}\n", .{
                res.chainLabel(),
                res.residue_name.slice(),
                num_str,
                res.sasa,
                res.atom_count,
            });
        }
    }
}

/// Print per-residue SASA results with RSA (Relative Solvent Accessibility)
pub fn printResidueResultsWithRsa(residues: []const ResidueSasa) void {
    std.debug.print("\nPer-residue SASA with RSA:\n", .{});
    std.debug.print("{s:>5} {s:>4} {s:>6} {s:>10} {s:>6} {s:>6}\n", .{
        "Chain", "Res", "Num", "SASA", "RSA", "Atoms",
    });
    std.debug.print("{s:->5} {s:->4} {s:->6} {s:->10} {s:->6} {s:->6}\n", .{
        "", "", "", "", "", "",
    });

    for (residues) |res| {
        // Format residue number as string to avoid Zig's "+" prefix for positive integers
        var num_buf: [16]u8 = undefined;
        const num_str = std.fmt.bufPrint(&num_buf, "{d}", .{res.residue_num}) catch "?";

        // Format RSA as percentage or N/A
        var rsa_buf: [8]u8 = undefined;
        const rsa_str = if (res.rsa) |rsa|
            std.fmt.bufPrint(&rsa_buf, "{d:.2}", .{rsa}) catch "?"
        else
            "N/A";

        if (res.insertion_code.len > 0) {
            std.debug.print("{s:>5} {s:>4} {s:>5}{s:<1} {d:>10.2} {s:>6} {d:>6}\n", .{
                res.chainLabel(),
                res.residue_name.slice(),
                num_str,
                res.insertion_code.slice(),
                res.sasa,
                rsa_str,
                res.atom_count,
            });
        } else {
            std.debug.print("{s:>5} {s:>4} {s:>6} {d:>10.2} {s:>6} {d:>6}\n", .{
                res.chainLabel(),
                res.residue_name.slice(),
                num_str,
                res.sasa,
                rsa_str,
                res.atom_count,
            });
        }
    }
}

// Tests
test "aggregateByResidue basic" {
    const allocator = std.testing.allocator;

    // Create mock input with 4 atoms in 2 residues
    const x = try allocator.alloc(f64, 4);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 4);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 4);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 4);
    defer allocator.free(r);

    // Chain IDs
    var chain_ids = try allocator.alloc(types.FixedString4, 4);
    defer allocator.free(chain_ids);
    chain_ids[0] = types.FixedString4.fromSlice("A");
    chain_ids[1] = types.FixedString4.fromSlice("A");
    chain_ids[2] = types.FixedString4.fromSlice("A");
    chain_ids[3] = types.FixedString4.fromSlice("A");

    // Residue names
    var residue_names = try allocator.alloc(types.FixedString5, 4);
    defer allocator.free(residue_names);
    residue_names[0] = types.FixedString5.fromSlice("ALA");
    residue_names[1] = types.FixedString5.fromSlice("ALA");
    residue_names[2] = types.FixedString5.fromSlice("GLY");
    residue_names[3] = types.FixedString5.fromSlice("GLY");

    // Residue numbers
    var residue_nums = try allocator.alloc(i32, 4);
    defer allocator.free(residue_nums);
    residue_nums[0] = 1;
    residue_nums[1] = 1;
    residue_nums[2] = 2;
    residue_nums[3] = 2;

    // Insertion codes (all empty)
    var insertion_codes = try allocator.alloc(types.FixedString4, 4);
    defer allocator.free(insertion_codes);
    insertion_codes[0] = types.FixedString4.fromSlice("");
    insertion_codes[1] = types.FixedString4.fromSlice("");
    insertion_codes[2] = types.FixedString4.fromSlice("");
    insertion_codes[3] = types.FixedString4.fromSlice("");

    const input = types.AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .chain_id = chain_ids,
        .residue = residue_names,
        .residue_num = residue_nums,
        .insertion_code = insertion_codes,
        .allocator = allocator,
    };

    // Mock atom areas
    const atom_areas = [_]f64{ 10.0, 15.0, 20.0, 25.0 };

    var result = try aggregateByResidue(allocator, input, &atom_areas);
    defer result.deinit();

    try std.testing.expectEqual(@as(usize, 2), result.residues.len);

    // First residue: ALA-1 with SASA = 10 + 15 = 25
    try std.testing.expectEqualStrings("ALA", result.residues[0].residue_name.slice());
    try std.testing.expectEqual(@as(i32, 1), result.residues[0].residue_num);
    try std.testing.expectEqual(@as(f64, 25.0), result.residues[0].sasa);
    try std.testing.expectEqual(@as(usize, 2), result.residues[0].atom_count);

    // Second residue: GLY-2 with SASA = 20 + 25 = 45
    try std.testing.expectEqualStrings("GLY", result.residues[1].residue_name.slice());
    try std.testing.expectEqual(@as(i32, 2), result.residues[1].residue_num);
    try std.testing.expectEqual(@as(f64, 45.0), result.residues[1].sasa);
    try std.testing.expectEqual(@as(usize, 2), result.residues[1].atom_count);

    // Check RSA values are calculated
    // ALA: RSA = 25.0 / 129.0 ≈ 0.194
    try std.testing.expect(result.residues[0].rsa != null);
    try std.testing.expectApproxEqRel(25.0 / 129.0, result.residues[0].rsa.?, 0.001);

    // GLY: RSA = 45.0 / 104.0 ≈ 0.433
    try std.testing.expect(result.residues[1].rsa != null);
    try std.testing.expectApproxEqRel(45.0 / 104.0, result.residues[1].rsa.?, 0.001);
}

test "aggregateByResidue groups by full chain IDs when present" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 2);
    defer allocator.free(x);
    const y = try allocator.alloc(f64, 2);
    defer allocator.free(y);
    const z = try allocator.alloc(f64, 2);
    defer allocator.free(z);
    const r = try allocator.alloc(f64, 2);
    defer allocator.free(r);
    @memset(x, 0);
    @memset(y, 0);
    @memset(z, 0);
    @memset(r, 1);

    var chain_ids = try allocator.alloc(types.FixedString4, 2);
    defer allocator.free(chain_ids);
    chain_ids[0] = types.FixedString4.fromSlice("ABCD");
    chain_ids[1] = types.FixedString4.fromSlice("ABCD");
    const chain_ids_full = [_][]const u8{ "ABCD1", "ABCD2" };

    var residue_names = try allocator.alloc(types.FixedString5, 2);
    defer allocator.free(residue_names);
    residue_names[0] = types.FixedString5.fromSlice("ALA");
    residue_names[1] = types.FixedString5.fromSlice("ALA");

    var residue_nums = try allocator.alloc(i32, 2);
    defer allocator.free(residue_nums);
    residue_nums[0] = 1;
    residue_nums[1] = 1;

    var insertion_codes = try allocator.alloc(types.FixedString4, 2);
    defer allocator.free(insertion_codes);
    insertion_codes[0] = types.FixedString4.fromSlice("");
    insertion_codes[1] = types.FixedString4.fromSlice("");

    const input = types.AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .chain_id = chain_ids,
        .chain_id_full = chain_ids_full[0..],
        .residue = residue_names,
        .residue_num = residue_nums,
        .insertion_code = insertion_codes,
        .allocator = allocator,
    };
    const atom_areas = [_]f64{ 10.0, 20.0 };

    var result = try aggregateByResidue(allocator, input, atom_areas[0..]);
    defer result.deinit();

    try std.testing.expectEqual(@as(usize, 2), result.residues.len);
    try std.testing.expectEqualStrings("ABCD1", result.residues[0].chainLabel());
    try std.testing.expectEqualStrings("ABCD2", result.residues[1].chainLabel());
}

test "aggregateByResidue owns full chain labels independently of input" {
    const allocator = std.testing.allocator;

    const x = try allocator.alloc(f64, 1);
    const y = try allocator.alloc(f64, 1);
    const z = try allocator.alloc(f64, 1);
    const r = try allocator.alloc(f64, 1);
    const chain_id = try allocator.alloc(types.FixedString4, 1);
    const chain_id_full = try allocator.alloc([]const u8, 1);
    const residue = try allocator.alloc(types.FixedString5, 1);
    const residue_num = try allocator.alloc(i32, 1);
    const insertion_code = try allocator.alloc(types.FixedString4, 1);

    x[0] = 0.0;
    y[0] = 0.0;
    z[0] = 0.0;
    r[0] = 1.0;
    chain_id[0] = types.FixedString4.fromSlice("ABCD");
    chain_id_full[0] = try allocator.dupe(u8, "ABCD1");
    residue[0] = types.FixedString5.fromSlice("ALA");
    residue_num[0] = 1;
    insertion_code[0] = types.FixedString4.fromSlice("");

    var input = types.AtomInput{
        .x = x,
        .y = y,
        .z = z,
        .r = r,
        .chain_id = chain_id,
        .chain_id_full = chain_id_full,
        .residue = residue,
        .residue_num = residue_num,
        .insertion_code = insertion_code,
        .allocator = allocator,
    };

    const atom_areas = [_]f64{10.0};
    var result = try aggregateByResidue(allocator, input, atom_areas[0..]);
    defer result.deinit();
    input.deinit();

    try std.testing.expectEqualStrings("ABCD1", result.residues[0].chainLabel());
}

test "MaxSASA lookup" {
    // Test known amino acids
    try std.testing.expectEqual(@as(f64, 129.0), MaxSASA.get("ALA").?);
    try std.testing.expectEqual(@as(f64, 104.0), MaxSASA.get("GLY").?);
    try std.testing.expectEqual(@as(f64, 285.0), MaxSASA.get("TRP").?);
    try std.testing.expectEqual(@as(f64, 274.0), MaxSASA.get("ARG").?);

    // Unknown residue returns null
    try std.testing.expect(MaxSASA.get("XXX") == null);
    try std.testing.expect(MaxSASA.get("HOH") == null);
}

test "ResidueSasa calculateRsa" {
    // Known residue
    var res_ala = ResidueSasa{
        .chain_id = types.FixedString4.fromSlice("A"),
        .residue_name = types.FixedString5.fromSlice("ALA"),
        .residue_num = 1,
        .insertion_code = types.FixedString4.fromSlice(""),
        .sasa = 64.5, // 50% of MaxSASA
        .atom_count = 5,
    };
    res_ala.calculateRsa();
    try std.testing.expect(res_ala.rsa != null);
    try std.testing.expectApproxEqRel(64.5 / 129.0, res_ala.rsa.?, 0.001);

    // Unknown residue
    var res_unk = ResidueSasa{
        .chain_id = types.FixedString4.fromSlice("A"),
        .residue_name = types.FixedString5.fromSlice("UNK"),
        .residue_num = 1,
        .insertion_code = types.FixedString4.fromSlice(""),
        .sasa = 100.0,
        .atom_count = 10,
    };
    res_unk.calculateRsa();
    try std.testing.expect(res_unk.rsa == null);

    // RSA > 1.0 is possible for exposed terminal residues
    var res_gly = ResidueSasa{
        .chain_id = types.FixedString4.fromSlice("A"),
        .residue_name = types.FixedString5.fromSlice("GLY"),
        .residue_num = 1,
        .insertion_code = types.FixedString4.fromSlice(""),
        .sasa = 150.0, // Exceeds MaxSASA of 104.0
        .atom_count = 4,
    };
    res_gly.calculateRsa();
    try std.testing.expect(res_gly.rsa != null);
    try std.testing.expect(res_gly.rsa.? > 1.0);
    try std.testing.expectApproxEqRel(150.0 / 104.0, res_gly.rsa.?, 0.001);
}

test "ResidueClass classification" {
    // Polar residues
    try std.testing.expectEqual(ResidueClass.polar, ResidueClass.fromResidueName("ARG"));
    try std.testing.expectEqual(ResidueClass.polar, ResidueClass.fromResidueName("ASP"));
    try std.testing.expectEqual(ResidueClass.polar, ResidueClass.fromResidueName("GLU"));
    try std.testing.expectEqual(ResidueClass.polar, ResidueClass.fromResidueName("LYS"));
    try std.testing.expectEqual(ResidueClass.polar, ResidueClass.fromResidueName("SER"));
    try std.testing.expectEqual(ResidueClass.polar, ResidueClass.fromResidueName("THR"));

    // Nonpolar residues
    try std.testing.expectEqual(ResidueClass.nonpolar, ResidueClass.fromResidueName("ALA"));
    try std.testing.expectEqual(ResidueClass.nonpolar, ResidueClass.fromResidueName("ILE"));
    try std.testing.expectEqual(ResidueClass.nonpolar, ResidueClass.fromResidueName("LEU"));
    try std.testing.expectEqual(ResidueClass.nonpolar, ResidueClass.fromResidueName("PHE"));
    try std.testing.expectEqual(ResidueClass.nonpolar, ResidueClass.fromResidueName("VAL"));

    // Unknown
    try std.testing.expectEqual(ResidueClass.unknown, ResidueClass.fromResidueName("UNK"));
    try std.testing.expectEqual(ResidueClass.unknown, ResidueClass.fromResidueName("HOH"));
}

test "calculatePolarSummary" {
    const residues = [_]ResidueSasa{
        .{ .chain_id = types.FixedString4.fromSlice("A"), .residue_name = types.FixedString5.fromSlice("ALA"), .residue_num = 1, .insertion_code = types.FixedString4.fromSlice(""), .sasa = 50.0, .atom_count = 5 },
        .{ .chain_id = types.FixedString4.fromSlice("A"), .residue_name = types.FixedString5.fromSlice("SER"), .residue_num = 2, .insertion_code = types.FixedString4.fromSlice(""), .sasa = 30.0, .atom_count = 6 },
        .{ .chain_id = types.FixedString4.fromSlice("A"), .residue_name = types.FixedString5.fromSlice("LEU"), .residue_num = 3, .insertion_code = types.FixedString4.fromSlice(""), .sasa = 40.0, .atom_count = 8 },
        .{ .chain_id = types.FixedString4.fromSlice("A"), .residue_name = types.FixedString5.fromSlice("ASP"), .residue_num = 4, .insertion_code = types.FixedString4.fromSlice(""), .sasa = 20.0, .atom_count = 8 },
    };

    const summary = calculatePolarSummary(&residues);

    // Polar: SER (30) + ASP (20) = 50
    try std.testing.expectEqual(@as(f64, 50.0), summary.polar_sasa);
    try std.testing.expectEqual(@as(usize, 2), summary.polar_residue_count);

    // Nonpolar: ALA (50) + LEU (40) = 90
    try std.testing.expectEqual(@as(f64, 90.0), summary.nonpolar_sasa);
    try std.testing.expectEqual(@as(usize, 2), summary.nonpolar_residue_count);

    // No unknown
    try std.testing.expectEqual(@as(f64, 0.0), summary.unknown_sasa);
    try std.testing.expectEqual(@as(usize, 0), summary.unknown_residue_count);

    // Fractions: polar = 50/140, nonpolar = 90/140
    try std.testing.expectApproxEqRel(50.0 / 140.0, summary.polarFraction(), 0.001);
    try std.testing.expectApproxEqRel(90.0 / 140.0, summary.nonpolarFraction(), 0.001);
}

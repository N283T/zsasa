//! Analysis module for SASA aggregation and classification.
//!
//! Provides functions for:
//! - Per-residue SASA aggregation
//! - RSA (Relative Solvent Accessibility) calculation
//! - Polar/nonpolar classification

const std = @import("std");
const types = @import("types.zig");

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
        switch (ResidueClass.fromResidueName(res.residue_name)) {
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

/// Per-chain SASA data for interface calculation
pub const ChainSasa = struct {
    chain_id: []const u8,
    complex_sasa: f64, // SASA when in complex
    isolated_sasa: f64, // SASA when isolated
    buried_sasa: f64, // Difference (buried at interface)
    atom_count: usize,
};

/// Interface SASA summary
pub const InterfaceSummary = struct {
    chains: []ChainSasa,
    total_complex_sasa: f64,
    total_isolated_sasa: f64,
    total_buried_sasa: f64,
    allocator: std.mem.Allocator,

    pub fn deinit(self: *InterfaceSummary) void {
        self.allocator.free(self.chains);
    }
};

/// Calculate per-chain SASA from atom areas and chain IDs.
/// Returns a map of chain_id -> total SASA for that chain.
pub fn calculatePerChainSasa(
    allocator: std.mem.Allocator,
    chain_ids: []const []const u8,
    atom_areas: []const f64,
) ![]ChainSasa {
    var chain_map = std.ArrayListUnmanaged(ChainSasa){};
    defer chain_map.deinit(allocator);

    for (chain_ids, 0..) |chain, i| {
        // Find existing chain or create new
        var found = false;
        for (chain_map.items) |*c| {
            if (std.mem.eql(u8, c.chain_id, chain)) {
                c.complex_sasa += atom_areas[i];
                c.atom_count += 1;
                found = true;
                break;
            }
        }
        if (!found) {
            try chain_map.append(allocator, ChainSasa{
                .chain_id = chain,
                .complex_sasa = atom_areas[i],
                .isolated_sasa = 0, // Will be set later
                .buried_sasa = 0, // Will be calculated later
                .atom_count = 1,
            });
        }
    }

    return try chain_map.toOwnedSlice(allocator);
}

/// Print interface SASA summary
pub fn printInterfaceSummary(summary: InterfaceSummary) void {
    std.debug.print("\nInterface SASA Analysis:\n", .{});
    std.debug.print("{s:>5} {s:>12} {s:>12} {s:>12} {s:>6}\n", .{
        "Chain", "Complex", "Isolated", "Buried", "Atoms",
    });
    std.debug.print("{s:->5} {s:->12} {s:->12} {s:->12} {s:->6}\n", .{
        "", "", "", "", "",
    });

    for (summary.chains) |chain| {
        std.debug.print("{s:>5} {d:>12.2} {d:>12.2} {d:>12.2} {d:>6}\n", .{
            chain.chain_id,
            chain.complex_sasa,
            chain.isolated_sasa,
            chain.buried_sasa,
            chain.atom_count,
        });
    }

    std.debug.print("{s:->5} {s:->12} {s:->12} {s:->12} {s:->6}\n", .{
        "", "", "", "", "",
    });
    std.debug.print("{s:>5} {d:>12.2} {d:>12.2} {d:>12.2}\n", .{
        "Total",
        summary.total_complex_sasa,
        summary.total_isolated_sasa,
        summary.total_buried_sasa,
    });

    if (summary.total_isolated_sasa > 0) {
        const buried_pct = summary.total_buried_sasa / summary.total_isolated_sasa * 100;
        std.debug.print("\nBuried surface: {d:.1}% of total isolated surface\n", .{buried_pct});
    }
}

/// Per-residue SASA data
pub const ResidueSasa = struct {
    chain_id: []const u8,
    residue_name: []const u8,
    residue_num: i32,
    insertion_code: []const u8,
    sasa: f64,
    atom_count: usize,
    /// Relative Solvent Accessibility (0.0-1.0+), null if MaxSASA unknown
    rsa: ?f64 = null,

    /// Calculate RSA from SASA and residue name.
    /// RSA values > 1.0 are possible for exposed terminal residues.
    pub fn calculateRsa(self: *ResidueSasa) void {
        if (MaxSASA.get(self.residue_name)) |max_sasa| {
            if (max_sasa > 0) {
                self.rsa = self.sasa / max_sasa;
            } else {
                self.rsa = null;
            }
        } else {
            self.rsa = null;
        }
    }
};

/// Result of per-residue aggregation
pub const ResidueResult = struct {
    residues: []ResidueSasa,
    allocator: std.mem.Allocator,

    pub fn deinit(self: *ResidueResult) void {
        self.allocator.free(self.residues);
    }
};

/// Aggregate atom SASA values to per-residue SASA.
/// Atoms are grouped by (chain_id, residue_num, insertion_code).
pub fn aggregateByResidue(
    allocator: std.mem.Allocator,
    input: types.AtomInput,
    atom_areas: []const f64,
) !ResidueResult {
    // Check if we have the required residue info
    const chain_ids = input.chain_id orelse return error.MissingChainInfo;
    const residue_names = input.residue orelse return error.MissingResidueInfo;
    const residue_nums = input.residue_num orelse return error.MissingResidueNumInfo;
    const insertion_codes = input.insertion_code orelse return error.MissingInsertionCodeInfo;

    const n = input.atomCount();
    if (atom_areas.len != n) {
        return error.LengthMismatch;
    }
    if (n == 0) {
        return ResidueResult{
            .residues = try allocator.alloc(ResidueSasa, 0),
            .allocator = allocator,
        };
    }

    // Use a simple approach: collect unique residues and sum areas
    // For efficiency, we use a fixed-size buffer then convert to owned slice
    var residue_list = std.ArrayListUnmanaged(ResidueSasa){};
    defer residue_list.deinit(allocator);

    for (0..n) |i| {
        const chain = chain_ids[i];
        const res_num = residue_nums[i];
        const ins_code = insertion_codes[i];

        // Find if this residue already exists
        var found_idx: ?usize = null;
        for (residue_list.items, 0..) |*res, j| {
            if (res.residue_num == res_num and
                std.mem.eql(u8, res.chain_id, chain) and
                std.mem.eql(u8, res.insertion_code, ins_code))
            {
                found_idx = j;
                break;
            }
        }

        if (found_idx) |idx| {
            // Add to existing residue
            residue_list.items[idx].sasa += atom_areas[i];
            residue_list.items[idx].atom_count += 1;
        } else {
            // Add new residue
            try residue_list.append(allocator, ResidueSasa{
                .chain_id = chain,
                .residue_name = residue_names[i],
                .residue_num = res_num,
                .insertion_code = ins_code,
                .sasa = atom_areas[i],
                .atom_count = 1,
            });
        }
    }

    // Calculate RSA for each residue
    for (residue_list.items) |*res| {
        res.calculateRsa();
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
                res.chain_id,
                res.residue_name,
                num_str,
                res.insertion_code,
                res.sasa,
                res.atom_count,
            });
        } else {
            std.debug.print("{s:>5} {s:>4} {s:>6} {d:>10.2} {d:>6}\n", .{
                res.chain_id,
                res.residue_name,
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
                res.chain_id,
                res.residue_name,
                num_str,
                res.insertion_code,
                res.sasa,
                rsa_str,
                res.atom_count,
            });
        } else {
            std.debug.print("{s:>5} {s:>4} {s:>6} {d:>10.2} {s:>6} {d:>6}\n", .{
                res.chain_id,
                res.residue_name,
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
    var chain_ids = try allocator.alloc([]const u8, 4);
    defer allocator.free(chain_ids);
    chain_ids[0] = "A";
    chain_ids[1] = "A";
    chain_ids[2] = "A";
    chain_ids[3] = "A";

    // Residue names
    var residue_names = try allocator.alloc([]const u8, 4);
    defer allocator.free(residue_names);
    residue_names[0] = "ALA";
    residue_names[1] = "ALA";
    residue_names[2] = "GLY";
    residue_names[3] = "GLY";

    // Residue numbers
    var residue_nums = try allocator.alloc(i32, 4);
    defer allocator.free(residue_nums);
    residue_nums[0] = 1;
    residue_nums[1] = 1;
    residue_nums[2] = 2;
    residue_nums[3] = 2;

    // Insertion codes (all empty)
    var insertion_codes = try allocator.alloc([]const u8, 4);
    defer allocator.free(insertion_codes);
    insertion_codes[0] = "";
    insertion_codes[1] = "";
    insertion_codes[2] = "";
    insertion_codes[3] = "";

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
    try std.testing.expectEqualStrings("ALA", result.residues[0].residue_name);
    try std.testing.expectEqual(@as(i32, 1), result.residues[0].residue_num);
    try std.testing.expectEqual(@as(f64, 25.0), result.residues[0].sasa);
    try std.testing.expectEqual(@as(usize, 2), result.residues[0].atom_count);

    // Second residue: GLY-2 with SASA = 20 + 25 = 45
    try std.testing.expectEqualStrings("GLY", result.residues[1].residue_name);
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
        .chain_id = "A",
        .residue_name = "ALA",
        .residue_num = 1,
        .insertion_code = "",
        .sasa = 64.5, // 50% of MaxSASA
        .atom_count = 5,
    };
    res_ala.calculateRsa();
    try std.testing.expect(res_ala.rsa != null);
    try std.testing.expectApproxEqRel(64.5 / 129.0, res_ala.rsa.?, 0.001);

    // Unknown residue
    var res_unk = ResidueSasa{
        .chain_id = "A",
        .residue_name = "UNK",
        .residue_num = 1,
        .insertion_code = "",
        .sasa = 100.0,
        .atom_count = 10,
    };
    res_unk.calculateRsa();
    try std.testing.expect(res_unk.rsa == null);

    // RSA > 1.0 is possible for exposed terminal residues
    var res_gly = ResidueSasa{
        .chain_id = "A",
        .residue_name = "GLY",
        .residue_num = 1,
        .insertion_code = "",
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
        .{ .chain_id = "A", .residue_name = "ALA", .residue_num = 1, .insertion_code = "", .sasa = 50.0, .atom_count = 5 },
        .{ .chain_id = "A", .residue_name = "SER", .residue_num = 2, .insertion_code = "", .sasa = 30.0, .atom_count = 6 },
        .{ .chain_id = "A", .residue_name = "LEU", .residue_num = 3, .insertion_code = "", .sasa = 40.0, .atom_count = 8 },
        .{ .chain_id = "A", .residue_name = "ASP", .residue_num = 4, .insertion_code = "", .sasa = 20.0, .atom_count = 8 },
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

//! Analysis module for SASA aggregation and classification.
//!
//! Provides functions for:
//! - Per-residue SASA aggregation
//! - RSA (Relative Solvent Accessibility) calculation
//! - Polar/nonpolar classification (future)

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

    /// Calculate RSA from SASA and residue name
    pub fn calculateRsa(self: *ResidueSasa) void {
        if (MaxSASA.get(self.residue_name)) |max_sasa| {
            self.rsa = self.sasa / max_sasa;
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
}

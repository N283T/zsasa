//! Analysis module for SASA aggregation and classification.
//!
//! Provides functions for:
//! - Per-residue SASA aggregation
//! - RSA (Relative Solvent Accessibility) calculation (future)
//! - Polar/nonpolar classification (future)

const std = @import("std");
const types = @import("types.zig");

/// Per-residue SASA data
pub const ResidueSasa = struct {
    chain_id: []const u8,
    residue_name: []const u8,
    residue_num: i32,
    insertion_code: []const u8,
    sasa: f64,
    atom_count: usize,
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
}

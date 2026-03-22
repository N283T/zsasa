//! mmCIF Parser for extracting atom_site data.
//!
//! This module provides a simplified mmCIF parser focused on extracting
//! atom coordinates from the _atom_site category for SASA calculation.
//!
//! ## Extracted Fields
//!
//! Required:
//! - Cartn_x, Cartn_y, Cartn_z: Cartesian coordinates
//! - type_symbol: Element symbol (for atomic number)
//!
//! Optional (for classification):
//! - label_atom_id / auth_atom_id: Atom name (e.g., CA, CB)
//! - label_comp_id / auth_comp_id: Residue name (e.g., ALA, GLY)
//! - group_PDB: Record type (ATOM/HETATM)
//! - label_alt_id: Alternate location indicator
//!
//! ## Usage
//!
//! ```zig
//! const parser = @import("mmcif_parser.zig");
//!
//! var mmcif = parser.MmcifParser.init(allocator);
//! // Note: Parser itself has no deinit - only the result needs cleanup
//!
//! const input = try mmcif.parseFile("structure.cif");
//! defer input.deinit();
//! ```

const std = @import("std");
const Allocator = std.mem.Allocator;
const cif = @import("cif_tokenizer.zig");
const elem = @import("element.zig");
const mmap_reader = @import("mmap_reader.zig");
const gzip = @import("gzip.zig");
const types = @import("types.zig");
const AtomInput = types.AtomInput;

/// Error types for mmCIF parsing
pub const ParseError = error{
    /// No atom_site loop found in the file
    NoAtomSiteLoop,
    /// Missing required coordinate field (Cartn_x, Cartn_y, or Cartn_z)
    MissingCoordinateField,
    /// Invalid coordinate value (not a valid number)
    InvalidCoordinate,
};

/// Column indices for atom_site fields
const AtomSiteColumns = struct {
    cartn_x: ?usize = null,
    cartn_y: ?usize = null,
    cartn_z: ?usize = null,
    type_symbol: ?usize = null,
    label_atom_id: ?usize = null,
    auth_atom_id: ?usize = null,
    label_comp_id: ?usize = null,
    auth_comp_id: ?usize = null,
    label_asym_id: ?usize = null,
    auth_asym_id: ?usize = null,
    label_seq_id: ?usize = null,
    auth_seq_id: ?usize = null,
    pdbx_pdb_ins_code: ?usize = null,
    group_pdb: ?usize = null,
    label_alt_id: ?usize = null,
    pdbx_pdb_model_num: ?usize = null,

    /// Check if required coordinate fields are present
    fn hasRequiredFields(self: AtomSiteColumns) bool {
        return self.cartn_x != null and self.cartn_y != null and self.cartn_z != null;
    }

    /// Get atom name column (prefer label over auth)
    fn getAtomNameCol(self: AtomSiteColumns) ?usize {
        return self.label_atom_id orelse self.auth_atom_id;
    }

    /// Get residue name column (prefer label over auth)
    fn getResNameCol(self: AtomSiteColumns) ?usize {
        return self.label_comp_id orelse self.auth_comp_id;
    }

    /// Get chain ID column (prefer label over auth by default)
    fn getChainCol(self: AtomSiteColumns, use_auth: bool) ?usize {
        if (use_auth) {
            return self.auth_asym_id orelse self.label_asym_id;
        }
        return self.label_asym_id orelse self.auth_asym_id;
    }

    /// Get residue sequence number column (prefer label over auth)
    fn getResSeqCol(self: AtomSiteColumns) ?usize {
        return self.label_seq_id orelse self.auth_seq_id;
    }

    /// Get insertion code column
    fn getInsCodeCol(self: AtomSiteColumns) ?usize {
        return self.pdbx_pdb_ins_code;
    }
};

/// mmCIF Parser
pub const MmcifParser = struct {
    allocator: Allocator,
    /// Filter to include only ATOM records (exclude HETATM)
    /// Default: true (matches FreeSASA/RustSASA behavior)
    atom_only: bool = true,
    /// Skip hydrogen atoms
    /// Default: true (matches FreeSASA/RustSASA behavior)
    skip_hydrogens: bool = true,
    /// Filter to include only first alternate location
    first_alt_loc_only: bool = true,
    /// Model number to extract (null = first model or all)
    model_num: ?u32 = null,
    /// Chain IDs to include (null = all chains)
    chain_filter: ?[]const []const u8 = null,
    /// Use auth_asym_id instead of label_asym_id for chain
    use_auth_chain: bool = false,

    pub fn init(allocator: Allocator) MmcifParser {
        return .{ .allocator = allocator };
    }

    /// Parse mmCIF from a string
    pub fn parse(self: *MmcifParser, source: []const u8) !AtomInput {
        var tokenizer = cif.Tokenizer.init(source);

        // Skip to atom_site loop
        const loop_info = try self.findAtomSiteLoop(&tokenizer);

        if (!loop_info.columns.hasRequiredFields()) {
            return ParseError.MissingCoordinateField;
        }

        // Parse atom data
        return try self.parseAtomData(&tokenizer, loop_info.columns, loop_info.num_cols);
    }

    /// Parse mmCIF from a file (handles both plain and .gz compressed)
    pub fn parseFile(self: *MmcifParser, path: []const u8) !AtomInput {
        if (std.mem.endsWith(u8, path, ".gz")) {
            const data = try gzip.readGzip(self.allocator, path);
            defer self.allocator.free(data);
            return self.parse(data);
        }
        const mapped = try mmap_reader.mmapFile(self.allocator, path);
        defer mapped.deinit();
        return self.parse(mapped.data);
    }

    /// Result from findAtomSiteLoop containing columns and their count
    const LoopInfo = struct {
        columns: AtomSiteColumns,
        num_cols: usize,
    };

    /// Find the atom_site loop and return column indices
    fn findAtomSiteLoop(self: *MmcifParser, tokenizer: *cif.Tokenizer) !LoopInfo {
        _ = self;
        var columns = AtomSiteColumns{};
        var in_atom_site_loop = false;
        var col_index: usize = 0;

        while (true) {
            // Save position before reading token
            const saved_pos = tokenizer.pos;
            const saved_line = tokenizer.line;
            const saved_col = tokenizer.col;

            const token = tokenizer.next();

            switch (token) {
                .eof => {
                    if (in_atom_site_loop and columns.hasRequiredFields()) {
                        return LoopInfo{ .columns = columns, .num_cols = col_index };
                    }
                    return ParseError.NoAtomSiteLoop;
                },
                .loop => {
                    // Start of a new loop - reset state
                    in_atom_site_loop = false;
                    col_index = 0;
                    columns = AtomSiteColumns{};
                },
                .tag => |tag| {
                    if (startsWithIgnoreCase(tag, "_atom_site.")) {
                        in_atom_site_loop = true;
                        const field = tag["_atom_site.".len..];

                        // Map field names to column indices (case-insensitive)
                        if (eqlIgnoreCase(field, "Cartn_x")) {
                            columns.cartn_x = col_index;
                        } else if (eqlIgnoreCase(field, "Cartn_y")) {
                            columns.cartn_y = col_index;
                        } else if (eqlIgnoreCase(field, "Cartn_z")) {
                            columns.cartn_z = col_index;
                        } else if (eqlIgnoreCase(field, "type_symbol")) {
                            columns.type_symbol = col_index;
                        } else if (eqlIgnoreCase(field, "label_atom_id")) {
                            columns.label_atom_id = col_index;
                        } else if (eqlIgnoreCase(field, "auth_atom_id")) {
                            columns.auth_atom_id = col_index;
                        } else if (eqlIgnoreCase(field, "label_comp_id")) {
                            columns.label_comp_id = col_index;
                        } else if (eqlIgnoreCase(field, "auth_comp_id")) {
                            columns.auth_comp_id = col_index;
                        } else if (eqlIgnoreCase(field, "label_asym_id")) {
                            columns.label_asym_id = col_index;
                        } else if (eqlIgnoreCase(field, "auth_asym_id")) {
                            columns.auth_asym_id = col_index;
                        } else if (eqlIgnoreCase(field, "label_seq_id")) {
                            columns.label_seq_id = col_index;
                        } else if (eqlIgnoreCase(field, "auth_seq_id")) {
                            columns.auth_seq_id = col_index;
                        } else if (eqlIgnoreCase(field, "pdbx_PDB_ins_code")) {
                            columns.pdbx_pdb_ins_code = col_index;
                        } else if (eqlIgnoreCase(field, "group_PDB")) {
                            columns.group_pdb = col_index;
                        } else if (eqlIgnoreCase(field, "label_alt_id")) {
                            columns.label_alt_id = col_index;
                        } else if (eqlIgnoreCase(field, "pdbx_PDB_model_num")) {
                            columns.pdbx_pdb_model_num = col_index;
                        }

                        col_index += 1;
                    } else if (in_atom_site_loop) {
                        // New category tag - we're done with atom_site
                        // Restore position so parseAtomData doesn't miss this
                        tokenizer.pos = saved_pos;
                        tokenizer.line = saved_line;
                        tokenizer.col = saved_col;
                        return LoopInfo{ .columns = columns, .num_cols = col_index };
                    }
                },
                .value => {
                    if (in_atom_site_loop) {
                        // We've hit values - columns are complete
                        // Restore position so parseAtomData reads from the first value
                        tokenizer.pos = saved_pos;
                        tokenizer.line = saved_line;
                        tokenizer.col = saved_col;
                        return LoopInfo{ .columns = columns, .num_cols = col_index };
                    }
                },
                .data_block => {
                    // New data block - if we were in atom_site, we're done
                    if (in_atom_site_loop and columns.hasRequiredFields()) {
                        tokenizer.pos = saved_pos;
                        tokenizer.line = saved_line;
                        tokenizer.col = saved_col;
                        return LoopInfo{ .columns = columns, .num_cols = col_index };
                    }
                },
                else => {},
            }
        }
    }

    /// Parse atom data from the loop values
    fn parseAtomData(
        self: *MmcifParser,
        tokenizer: *cif.Tokenizer,
        columns: AtomSiteColumns,
        num_cols: usize,
    ) !AtomInput {
        // Dynamic arrays for collecting data
        var x_list = std.ArrayListUnmanaged(f64){};
        defer x_list.deinit(self.allocator);
        var y_list = std.ArrayListUnmanaged(f64){};
        defer y_list.deinit(self.allocator);
        var z_list = std.ArrayListUnmanaged(f64){};
        defer z_list.deinit(self.allocator);
        var r_list = std.ArrayListUnmanaged(f64){};
        defer r_list.deinit(self.allocator);
        var residue_list = std.ArrayListUnmanaged(types.FixedString5){};
        defer residue_list.deinit(self.allocator);
        var atom_name_list = std.ArrayListUnmanaged(types.FixedString4){};
        defer atom_name_list.deinit(self.allocator);
        var element_list = std.ArrayListUnmanaged(u8){};
        defer element_list.deinit(self.allocator);
        var chain_id_list = std.ArrayListUnmanaged(types.FixedString4){};
        defer chain_id_list.deinit(self.allocator);
        var residue_num_list = std.ArrayListUnmanaged(i32){};
        defer residue_num_list.deinit(self.allocator);
        var insertion_code_list = std.ArrayListUnmanaged(types.FixedString4){};
        defer insertion_code_list.deinit(self.allocator);

        // Buffer for current row values
        var row_values = try self.allocator.alloc([]const u8, num_cols);
        defer self.allocator.free(row_values);

        var first_alt_loc: ?u8 = null;
        var col: usize = 0;

        // We need to re-read from the beginning of values since findAtomSiteLoop
        // consumed the first value. Let's backtrack by using a separate tokenizer
        // Actually, we should handle this differently - the tokenizer already consumed one value
        // For now, let's read values one by one

        while (true) {
            const token = tokenizer.next();

            switch (token) {
                .value => |value| {
                    row_values[col] = value;
                    col += 1;

                    if (col >= num_cols) {
                        // Complete row - process it
                        const should_include = try self.shouldIncludeAtom(
                            row_values,
                            columns,
                            &first_alt_loc,
                        );

                        if (should_include) {
                            // Parse coordinates
                            const x = try parseFloat(row_values[columns.cartn_x.?]);
                            const y = try parseFloat(row_values[columns.cartn_y.?]);
                            const z = try parseFloat(row_values[columns.cartn_z.?]);

                            try x_list.append(self.allocator, x);
                            try y_list.append(self.allocator, y);
                            try z_list.append(self.allocator, z);

                            // Get element and radius
                            var element_enum = elem.Element.X;
                            if (columns.type_symbol) |type_col| {
                                const symbol = row_values[type_col];
                                if (!cif.isNull(symbol)) {
                                    element_enum = elem.fromSymbol(symbol);
                                }
                            }

                            // Use VdW radius as default (classifier will override)
                            try r_list.append(self.allocator, element_enum.vdwRadius());
                            try element_list.append(self.allocator, element_enum.atomicNumber());

                            // Get residue name
                            if (columns.getResNameCol()) |res_col| {
                                const res = row_values[res_col];
                                if (cif.isNull(res)) {
                                    try residue_list.append(self.allocator, types.FixedString5.fromSlice("UNK"));
                                } else {
                                    try residue_list.append(self.allocator, types.FixedString5.fromSlice(res));
                                }
                            } else {
                                try residue_list.append(self.allocator, types.FixedString5.fromSlice("UNK"));
                            }

                            // Get atom name
                            if (columns.getAtomNameCol()) |atom_col| {
                                const name = row_values[atom_col];
                                if (cif.isNull(name)) {
                                    try atom_name_list.append(self.allocator, types.FixedString4.fromSlice("X"));
                                } else {
                                    try atom_name_list.append(self.allocator, types.FixedString4.fromSlice(name));
                                }
                            } else {
                                try atom_name_list.append(self.allocator, types.FixedString4.fromSlice("X"));
                            }

                            // Get chain ID
                            if (columns.getChainCol(self.use_auth_chain)) |chain_col| {
                                const chain = row_values[chain_col];
                                if (cif.isNull(chain)) {
                                    try chain_id_list.append(self.allocator, types.FixedString4.fromSlice(""));
                                } else {
                                    try chain_id_list.append(self.allocator, types.FixedString4.fromSlice(chain));
                                }
                            } else {
                                try chain_id_list.append(self.allocator, types.FixedString4.fromSlice(""));
                            }

                            // Get residue sequence number
                            if (columns.getResSeqCol()) |seq_col| {
                                const seq_str = row_values[seq_col];
                                if (cif.isNull(seq_str)) {
                                    try residue_num_list.append(self.allocator, 0);
                                } else {
                                    const seq_num = std.fmt.parseInt(i32, seq_str, 10) catch 0;
                                    try residue_num_list.append(self.allocator, seq_num);
                                }
                            } else {
                                try residue_num_list.append(self.allocator, 0);
                            }

                            // Get insertion code
                            if (columns.getInsCodeCol()) |ins_col| {
                                const ins_code = row_values[ins_col];
                                if (cif.isNull(ins_code)) {
                                    try insertion_code_list.append(self.allocator, types.FixedString4.fromSlice(""));
                                } else {
                                    try insertion_code_list.append(self.allocator, types.FixedString4.fromSlice(ins_code));
                                }
                            } else {
                                try insertion_code_list.append(self.allocator, types.FixedString4.fromSlice(""));
                            }
                        }

                        col = 0;
                    }
                },
                .eof, .loop, .data_block => {
                    // End of atom_site data
                    break;
                },
                .tag => {
                    // New category - end of atom_site
                    break;
                },
                else => {},
            }
        }

        // Convert to AtomInput
        const n = x_list.items.len;
        if (n == 0) {
            return ParseError.NoAtomSiteLoop;
        }

        return AtomInput{
            .x = try x_list.toOwnedSlice(self.allocator),
            .y = try y_list.toOwnedSlice(self.allocator),
            .z = try z_list.toOwnedSlice(self.allocator),
            .r = try r_list.toOwnedSlice(self.allocator),
            .residue = try residue_list.toOwnedSlice(self.allocator),
            .atom_name = try atom_name_list.toOwnedSlice(self.allocator),
            .element = try element_list.toOwnedSlice(self.allocator),
            .chain_id = try chain_id_list.toOwnedSlice(self.allocator),
            .residue_num = try residue_num_list.toOwnedSlice(self.allocator),
            .insertion_code = try insertion_code_list.toOwnedSlice(self.allocator),
            .allocator = self.allocator,
        };
    }

    /// Check if an atom should be included based on filters
    fn shouldIncludeAtom(
        self: *MmcifParser,
        row_values: []const []const u8,
        columns: AtomSiteColumns,
        first_alt_loc: *?u8,
    ) !bool {
        // Check group_PDB filter (ATOM vs HETATM)
        if (self.atom_only) {
            if (columns.group_pdb) |col| {
                const group = row_values[col];
                if (!std.mem.eql(u8, group, "ATOM")) {
                    return false;
                }
            }
        }

        // Check hydrogen filter
        if (self.skip_hydrogens) {
            if (columns.type_symbol) |col| {
                const symbol = row_values[col];
                if (!cif.isNull(symbol) and (std.mem.eql(u8, symbol, "H") or std.mem.eql(u8, symbol, "D"))) {
                    return false;
                }
            }
        }

        // Check alternate location filter
        if (self.first_alt_loc_only) {
            if (columns.label_alt_id) |col| {
                const alt_id = row_values[col];
                if (!cif.isNull(alt_id) and alt_id.len > 0) {
                    const alt_char = alt_id[0];
                    if (first_alt_loc.*) |first| {
                        if (alt_char != first) {
                            return false;
                        }
                    } else {
                        first_alt_loc.* = alt_char;
                    }
                }
            }
        }

        // Check model number filter
        if (self.model_num) |target_model| {
            if (columns.pdbx_pdb_model_num) |col| {
                const model_str = row_values[col];
                if (!cif.isNull(model_str)) {
                    const model = std.fmt.parseInt(u32, model_str, 10) catch 1;
                    if (model != target_model) {
                        return false;
                    }
                }
            }
        }

        // Check chain filter
        if (self.chain_filter) |chains| {
            if (columns.getChainCol(self.use_auth_chain)) |col| {
                const chain = row_values[col];
                if (cif.isNull(chain)) {
                    return false;
                }
                // Check if chain is in the filter list
                var found = false;
                for (chains) |target_chain| {
                    if (std.mem.eql(u8, chain, target_chain)) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    return false;
                }
            }
        }

        return true;
    }
};

/// Parse a float from a string, handling CIF null values
fn parseFloat(s: []const u8) !f64 {
    if (cif.isNull(s)) {
        return ParseError.InvalidCoordinate;
    }

    // Handle parenthetical uncertainty notation: "1.234(5)" -> "1.234"
    var end = s.len;
    for (s, 0..) |c, i| {
        if (c == '(') {
            end = i;
            break;
        }
    }

    return std.fmt.parseFloat(f64, s[0..end]) catch {
        return ParseError.InvalidCoordinate;
    };
}

/// Case-insensitive string comparison
fn eqlIgnoreCase(a: []const u8, b: []const u8) bool {
    if (a.len != b.len) return false;
    for (a, b) |ca, cb| {
        if (std.ascii.toLower(ca) != std.ascii.toLower(cb)) {
            return false;
        }
    }
    return true;
}

/// Case-insensitive prefix check
fn startsWithIgnoreCase(haystack: []const u8, prefix: []const u8) bool {
    if (haystack.len < prefix.len) return false;
    for (haystack[0..prefix.len], prefix) |h, p| {
        if (std.ascii.toLower(h) != std.ascii.toLower(p)) {
            return false;
        }
    }
    return true;
}

// ============================================================================
// Tests
// ============================================================================

test "parse simple mmCIF" {
    const source =
        \\data_TEST
        \\#
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA 10.000 20.000 30.000
        \\2 N N  ALA 11.000 21.000 31.000
        \\3 O O  ALA 12.000 22.000 32.000
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());

    // Check coordinates
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 20.0), input.y[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 30.0), input.z[0], 0.001);

    try std.testing.expectApproxEqAbs(@as(f64, 11.0), input.x[1], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 21.0), input.y[1], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 31.0), input.z[1], 0.001);

    // Check residue names
    try std.testing.expectEqualStrings("ALA", input.residue.?[0].slice());
    try std.testing.expectEqualStrings("ALA", input.residue.?[1].slice());

    // Check atom names
    try std.testing.expectEqualStrings("CA", input.atom_name.?[0].slice());
    try std.testing.expectEqualStrings("N", input.atom_name.?[1].slice());

    // Check elements
    try std.testing.expectEqual(@as(u8, 6), input.element.?[0]); // C
    try std.testing.expectEqual(@as(u8, 7), input.element.?[1]); // N
    try std.testing.expectEqual(@as(u8, 8), input.element.?[2]); // O
}

test "parse with alternate locations" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_alt_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA . 10.0 20.0 30.0
        \\2 C CB ALA A 11.0 21.0 31.0
        \\3 C CB ALA B 11.5 21.5 31.5
        \\4 N N  ALA . 12.0 22.0 32.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    parser.first_alt_loc_only = true;
    var input = try parser.parse(source);
    defer input.deinit();

    // Should have 3 atoms (CB with alt B should be excluded)
    try std.testing.expectEqual(@as(usize, 3), input.atomCount());

    try std.testing.expectEqualStrings("CA", input.atom_name.?[0].slice());
    try std.testing.expectEqualStrings("CB", input.atom_name.?[1].slice()); // alt A
    try std.testing.expectEqualStrings("N", input.atom_name.?[2].slice());
}

test "parse with parenthetical uncertainty" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C 10.123(45) 20.456(67) 30.789(89)
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.123), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 20.456), input.y[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 30.789), input.z[0], 0.001);
}

test "case insensitive field names" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.ID
        \\_atom_site.TYPE_SYMBOL
        \\_atom_site.CARTN_X
        \\_atom_site.CARTN_Y
        \\_atom_site.CARTN_Z
        \\1 C 10.0 20.0 30.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
}

test "missing required fields" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\1 C 10.0 20.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    const result = parser.parse(source);
    try std.testing.expectError(ParseError.MissingCoordinateField, result);
}

test "no atom_site loop" {
    const source =
        \\data_TEST
        \\loop_
        \\_cell.length_a
        \\_cell.length_b
        \\10.0 20.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    const result = parser.parse(source);
    try std.testing.expectError(ParseError.NoAtomSiteLoop, result);
}

test "parse with HETATM filter (default atom_only=true)" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.group_PDB
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 N N   ALA ATOM   10.0 20.0 30.0
        \\2 C CA  ALA ATOM   11.0 21.0 31.0
        \\3 O O   HOH HETATM 12.0 22.0 32.0
        \\#
    ;

    // Default: atom_only=true -> HETATM excluded
    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());

    // Include HETATM
    var parser2 = MmcifParser.init(std.testing.allocator);
    parser2.atom_only = false;
    var input2 = try parser2.parse(source);
    defer input2.deinit();

    try std.testing.expectEqual(@as(usize, 3), input2.atomCount());
}

test "parse with hydrogen filter (default skip_hydrogens=true)" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.group_PDB
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 N N   ALA ATOM 10.0 20.0 30.0
        \\2 C CA  ALA ATOM 11.0 21.0 31.0
        \\3 H H   ALA ATOM 12.0 22.0 32.0
        \\4 H HB  ALA ATOM 13.0 23.0 33.0
        \\5 O O   ALA ATOM 14.0 24.0 34.0
        \\#
    ;

    // Default: skip_hydrogens=true -> H atoms excluded
    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount()); // N, CA, O

    // Include hydrogens
    var parser2 = MmcifParser.init(std.testing.allocator);
    parser2.skip_hydrogens = false;
    var input2 = try parser2.parse(source);
    defer input2.deinit();

    try std.testing.expectEqual(@as(usize, 5), input2.atomCount()); // All atoms
}

test "parse with deuterium filter" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 N N  ALA 10.0 20.0 30.0
        \\2 D D  ALA 11.0 21.0 31.0
        \\3 C CA ALA 12.0 22.0 32.0
        \\#
    ;

    // Default: skip_hydrogens=true also skips deuterium (D)
    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount()); // N, CA
}

test "eqlIgnoreCase" {
    try std.testing.expect(eqlIgnoreCase("Cartn_x", "CARTN_X"));
    try std.testing.expect(eqlIgnoreCase("cartn_x", "Cartn_x"));
    try std.testing.expect(!eqlIgnoreCase("Cartn_x", "Cartn_y"));
    try std.testing.expect(!eqlIgnoreCase("short", "longer"));
}

test "startsWithIgnoreCase" {
    try std.testing.expect(startsWithIgnoreCase("_atom_site.Cartn_x", "_atom_site."));
    try std.testing.expect(startsWithIgnoreCase("_ATOM_SITE.Cartn_x", "_atom_site."));
    try std.testing.expect(startsWithIgnoreCase("_Atom_Site.Cartn_x", "_atom_site."));
    try std.testing.expect(!startsWithIgnoreCase("_cell.length_a", "_atom_site."));
    try std.testing.expect(!startsWithIgnoreCase("_atom", "_atom_site."));
}

test "fuzz mmcif parser" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, input: []const u8) !void {
            var parser = MmcifParser.init(std.testing.allocator);
            var result = parser.parse(input) catch return;
            result.deinit();
        }
    }.testOne, .{
        .corpus = &.{
            "data_TEST\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C 10.0 20.0 30.0\n#\n",
            "data_TEST\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C CA ALA 10.0 20.0 30.0\n2 N N ALA 11.0 21.0 31.0\n#\n",
            "data_TEST\nloop_\n_cell.length_a\n_cell.length_b\n10.0 20.0\n#\n",
        },
    });
}

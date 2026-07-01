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
const compressed = @import("compressed.zig");
const types = @import("types.zig");
const AtomInput = types.AtomInput;
const ccd_parser = @import("ccd_parser.zig");

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
    occupancy: ?usize = null,
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
    /// Parse inline CCD data from `_chem_comp_atom`/`_chem_comp_bond` loops.
    /// Disable this when the caller will not use inline CCD resources.
    parse_inline_ccd: bool = true,
    /// Inline CCD data parsed from `_chem_comp_atom`/`_chem_comp_bond` loops
    inline_ccd: ?ccd_parser.ComponentDict = null,

    pub fn init(allocator: Allocator) MmcifParser {
        return .{ .allocator = allocator };
    }

    /// Return a pointer to the inline CCD data, or null if none was parsed.
    pub fn getInlineCcd(self: *MmcifParser) ?*ccd_parser.ComponentDict {
        if (self.inline_ccd != null) return &self.inline_ccd.?;
        return null;
    }

    /// Move inline CCD data out of the parser. The caller owns the returned
    /// dictionary and must call `.deinit()` on it.
    pub fn takeInlineCcd(self: *MmcifParser) ?ccd_parser.ComponentDict {
        const dict = self.inline_ccd;
        self.inline_ccd = null;
        return dict;
    }

    /// Clean up inline CCD data. Must be called before the parser goes out of
    /// scope when `parse()` has been called and `getInlineCcd()` is non-null.
    pub fn deinitCcd(self: *MmcifParser) void {
        if (self.inline_ccd) |*dict| {
            dict.deinit();
            self.inline_ccd = null;
        }
    }

    /// Parse mmCIF from a string
    pub fn parse(self: *MmcifParser, source: []const u8) !AtomInput {
        self.deinitCcd();

        // First pass: extract any inline CCD data (_chem_comp_atom / _chem_comp_bond)
        // only when the caller can use it. Most protein-only batch workflows do
        // not need this full-file scan.
        if (self.parse_inline_ccd and hasInlineCcdAtomTag(source)) {
            var ccd_dict = try ccd_parser.parseCcdData(self.allocator, source, null);
            if (ccd_dict.count() == 0) {
                ccd_dict.deinit();
                self.inline_ccd = null;
            } else {
                self.inline_ccd = ccd_dict;
            }
        } else {
            self.inline_ccd = null;
        }

        // Second pass: parse _atom_site loop for coordinates
        var tokenizer = cif.Tokenizer.init(source);

        // Skip to atom_site loop
        const loop_info = self.findAtomSiteLoop(&tokenizer) catch |err| {
            // Clean up CCD data on failure
            self.deinitCcd();
            return err;
        };

        if (!loop_info.columns.hasRequiredFields()) {
            self.deinitCcd();
            return ParseError.MissingCoordinateField;
        }

        // Parse atom data
        return self.parseAtomData(&tokenizer, loop_info.columns, loop_info.num_cols) catch |err| {
            self.deinitCcd();
            return err;
        };
    }

    /// Parse mmCIF from a file (handles plain, .gz, and .zst compressed)
    pub fn parseFile(self: *MmcifParser, io: std.Io, path: []const u8) !AtomInput {
        if (compressed.isCompressed(path)) {
            const data = try compressed.read(self.allocator, path);
            defer self.allocator.free(data);
            return self.parse(data);
        }
        const mapped = try mmap_reader.mmapFile(self.allocator, io, path);
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
                        } else if (eqlIgnoreCase(field, "occupancy")) {
                            columns.occupancy = col_index;
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
        var x_list = std.ArrayListUnmanaged(f64).empty;
        defer x_list.deinit(self.allocator);
        var y_list = std.ArrayListUnmanaged(f64).empty;
        defer y_list.deinit(self.allocator);
        var z_list = std.ArrayListUnmanaged(f64).empty;
        defer z_list.deinit(self.allocator);
        var r_list = std.ArrayListUnmanaged(f64).empty;
        defer r_list.deinit(self.allocator);
        var residue_list = std.ArrayListUnmanaged(types.FixedString5).empty;
        defer residue_list.deinit(self.allocator);
        var atom_name_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer atom_name_list.deinit(self.allocator);
        var element_list = std.ArrayListUnmanaged(u8).empty;
        defer element_list.deinit(self.allocator);
        var chain_id_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer chain_id_list.deinit(self.allocator);
        var chain_id_full_list = std.ArrayListUnmanaged([]const u8).empty;
        defer chain_id_full_list.deinit(self.allocator);
        errdefer freeStringItems(self.allocator, chain_id_full_list.items);
        var has_extended_chain = false;
        var residue_num_list = std.ArrayListUnmanaged(i32).empty;
        defer residue_num_list.deinit(self.allocator);
        var insertion_code_list = std.ArrayListUnmanaged(types.FixedString4).empty;
        defer insertion_code_list.deinit(self.allocator);
        var atom_records = std.ArrayListUnmanaged(AtomRecord).empty;
        defer atom_records.deinit(self.allocator);
        var has_non_blank_alt_loc = false;

        // Buffer for current row values
        var row_values = try self.allocator.alloc([]const u8, num_cols);
        defer self.allocator.free(row_values);
        const estimated_rows = estimateAtomSiteRows(tokenizer.source.len - tokenizer.pos, num_cols);
        try atom_records.ensureTotalCapacity(self.allocator, estimated_rows);

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
                        const should_include = try self.shouldIncludeAtom(row_values, columns);

                        if (should_include) {
                            const atom = try self.atomRecordFromRow(row_values, columns);
                            has_non_blank_alt_loc = has_non_blank_alt_loc or atom.alt_loc != ' ';
                            try atom_records.append(self.allocator, atom);
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

        const max_output_atoms = atom_records.items.len;
        try x_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try y_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try z_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try r_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try residue_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try atom_name_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try element_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try chain_id_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try residue_num_list.ensureTotalCapacity(self.allocator, max_output_atoms);
        try insertion_code_list.ensureTotalCapacity(self.allocator, max_output_atoms);

        for (atom_records.items, 0..) |atom, i| {
            if (has_non_blank_alt_loc and !self.shouldKeepAltLoc(atom_records.items, i)) continue;

            try x_list.append(self.allocator, atom.x);
            try y_list.append(self.allocator, atom.y);
            try z_list.append(self.allocator, atom.z);
            try r_list.append(self.allocator, atom.radius);
            try element_list.append(self.allocator, atom.element.atomicNumber());
            try residue_list.append(self.allocator, types.FixedString5.fromSlice(atom.residue));
            try atom_name_list.append(self.allocator, types.FixedString4.fromSlice(atom.atom_name));
            try chain_id_list.append(self.allocator, types.FixedString4.fromSlice(atom.chain_id));
            if (!has_extended_chain and atom.chain_id.len > 4) {
                try chain_id_full_list.ensureTotalCapacity(self.allocator, max_output_atoms);
                for (chain_id_list.items[0 .. chain_id_list.items.len - 1]) |chain_id| {
                    chain_id_full_list.appendAssumeCapacity(try self.allocator.dupe(u8, chain_id.slice()));
                }
                has_extended_chain = true;
            }
            if (has_extended_chain) {
                chain_id_full_list.appendAssumeCapacity(try self.allocator.dupe(u8, atom.chain_id));
            }
            try residue_num_list.append(self.allocator, atom.residue_num);
            try insertion_code_list.append(self.allocator, types.FixedString4.fromSlice(atom.insertion_code));
        }

        // Convert to AtomInput
        const n = x_list.items.len;
        if (n == 0) {
            return ParseError.NoAtomSiteLoop;
        }

        const x = try x_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(x);
        const y = try y_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(y);
        const z = try z_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(z);
        const r = try r_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(r);
        const residue = try residue_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(residue);
        const atom_name = try atom_name_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(atom_name);
        const element = try element_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(element);
        const chain_id = try chain_id_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(chain_id);
        const chain_id_full: ?[]const []const u8 = if (has_extended_chain)
            try chain_id_full_list.toOwnedSlice(self.allocator)
        else
            null;
        errdefer if (chain_id_full) |chains| {
            freeStringItems(self.allocator, chains);
            self.allocator.free(chains);
        };
        const residue_num = try residue_num_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(residue_num);
        const insertion_code = try insertion_code_list.toOwnedSlice(self.allocator);
        errdefer self.allocator.free(insertion_code);

        return AtomInput{
            .x = x,
            .y = y,
            .z = z,
            .r = r,
            .residue = residue,
            .atom_name = atom_name,
            .element = element,
            .chain_id = chain_id,
            .chain_id_full = chain_id_full,
            .residue_num = residue_num,
            .insertion_code = insertion_code,
            .allocator = self.allocator,
        };
    }

    const AtomRecord = struct {
        x: f64,
        y: f64,
        z: f64,
        radius: f64,
        element: elem.Element,
        atom_name: []const u8,
        residue: []const u8,
        chain_id: []const u8,
        residue_num: i32,
        insertion_code: []const u8,
        alt_loc: u8,
        occupancy: f64,
        model_num: ?u32,
    };

    fn atomRecordFromRow(self: *MmcifParser, row_values: []const []const u8, columns: AtomSiteColumns) !AtomRecord {
        const x = try parseFloat(row_values[columns.cartn_x.?]);
        const y = try parseFloat(row_values[columns.cartn_y.?]);
        const z = try parseFloat(row_values[columns.cartn_z.?]);

        var element_enum = elem.Element.X;
        if (columns.type_symbol) |type_col| {
            const symbol = row_values[type_col];
            if (!cif.isNull(symbol)) {
                element_enum = elem.fromSymbol(symbol);
            }
        }

        const residue = if (columns.getResNameCol()) |res_col| blk: {
            const res = row_values[res_col];
            break :blk if (cif.isNull(res)) "UNK" else res;
        } else "UNK";

        const atom_name = if (columns.getAtomNameCol()) |atom_col| blk: {
            const name = row_values[atom_col];
            break :blk if (cif.isNull(name)) "X" else name;
        } else "X";

        const chain_id = if (columns.getChainCol(self.use_auth_chain)) |chain_col| blk: {
            const chain = row_values[chain_col];
            break :blk if (cif.isNull(chain)) "" else chain;
        } else "";

        const residue_num = if (columns.getResSeqCol()) |seq_col| blk: {
            const seq_str = row_values[seq_col];
            break :blk if (cif.isNull(seq_str)) 0 else std.fmt.parseInt(i32, seq_str, 10) catch 0;
        } else 0;

        const insertion_code = if (columns.getInsCodeCol()) |ins_col| blk: {
            const ins_code = row_values[ins_col];
            break :blk if (cif.isNull(ins_code)) "" else ins_code;
        } else "";

        const alt_loc: u8 = if (columns.label_alt_id) |alt_col| blk: {
            const alt_id = row_values[alt_col];
            break :blk if (cif.isNull(alt_id) or alt_id.len == 0) ' ' else alt_id[0];
        } else ' ';

        const occupancy = if (columns.occupancy) |occ_col| blk: {
            const occ = row_values[occ_col];
            break :blk if (cif.isNull(occ)) 0.0 else std.fmt.parseFloat(f64, occ) catch 0.0;
        } else 0.0;
        const model_num = if (columns.pdbx_pdb_model_num) |model_col| blk: {
            const model = row_values[model_col];
            break :blk if (cif.isNull(model)) null else std.fmt.parseInt(u32, model, 10) catch null;
        } else null;

        return .{
            .x = x,
            .y = y,
            .z = z,
            .radius = element_enum.vdwRadius(),
            .element = element_enum,
            .atom_name = atom_name,
            .residue = residue,
            .chain_id = chain_id,
            .residue_num = residue_num,
            .insertion_code = insertion_code,
            .alt_loc = alt_loc,
            .occupancy = occupancy,
            .model_num = model_num,
        };
    }

    fn sameAltLocSite(a: AtomRecord, b: AtomRecord) bool {
        return a.model_num == b.model_num and
            a.residue_num == b.residue_num and
            std.mem.eql(u8, a.chain_id, b.chain_id) and
            std.mem.eql(u8, a.residue, b.residue) and
            std.mem.eql(u8, a.insertion_code, b.insertion_code) and
            std.mem.eql(u8, a.atom_name, b.atom_name);
    }

    fn shouldKeepAltLoc(self: *MmcifParser, atoms: []const AtomRecord, index: usize) bool {
        if (!self.first_alt_loc_only) return true;

        const atom = atoms[index];
        if (atom.alt_loc == ' ') return true;

        var best_non_preferred: ?usize = null;
        for (atoms, 0..) |other, other_index| {
            if (!sameAltLocSite(atom, other)) continue;
            if (other.alt_loc == ' ') return false;
            if (other.alt_loc == 'A') return atom.alt_loc == 'A';
            if (best_non_preferred) |best_index| {
                if (other.occupancy > atoms[best_index].occupancy) {
                    best_non_preferred = other_index;
                }
            } else {
                best_non_preferred = other_index;
            }
        }
        return best_non_preferred == index;
    }

    /// Check if an atom should be included based on filters
    fn shouldIncludeAtom(
        self: *MmcifParser,
        row_values: []const []const u8,
        columns: AtomSiteColumns,
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

fn hasInlineCcdAtomTag(source: []const u8) bool {
    return std.mem.indexOf(u8, source, "_chem_comp_atom.") != null;
}

fn freeStringItems(allocator: Allocator, items: []const []const u8) void {
    for (items) |item| allocator.free(item);
}

fn estimateAtomSiteRows(remaining_bytes: usize, num_cols: usize) usize {
    if (remaining_bytes == 0) return 16;
    const bytes_per_col = 4;
    const min_row_bytes = @max(num_cols * bytes_per_col, 96);
    return @max(@as(usize, 16), @min(remaining_bytes / min_row_bytes, 128 * 1024));
}

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

test "parse mmCIF keeps extended chain IDs when label_asym_id exceeds four characters" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg 10.000 20.000 30.000
        \\2 N N  ALA ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef  11.000 21.000 31.000
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expect(input.chain_id_full != null);
    try std.testing.expectEqualStrings("ABCD", input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefg", input.chain_id_full.?[0]);
    try std.testing.expectEqualStrings("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdef", input.chain_id_full.?[1]);
}

test "parse mmCIF backfills extended chain IDs after short chains" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA A         10.000 20.000 30.000
        \\2 N N  GLY LONGCHAIN 11.000 21.000 31.000
        \\3 O O  SER B         12.000 22.000 32.000
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
    try std.testing.expect(input.chain_id_full != null);
    try std.testing.expectEqual(input.atomCount(), input.chain_id_full.?.len);
    try std.testing.expectEqualStrings("A", input.chain_id_full.?[0]);
    try std.testing.expectEqualStrings("LONGCHAIN", input.chain_id_full.?[1]);
    try std.testing.expectEqualStrings("B", input.chain_id_full.?[2]);
    try std.testing.expectEqualStrings("LONG", input.chain_id.?[1].slice());
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

test "parse mmCIF default model selection includes all models" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_model_num
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA A 1 1 10.0 20.0 30.0
        \\2 C CA GLY B 2 2 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectEqualStrings("A", input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("B", input.chain_id.?[1].slice());
}

test "parse mmCIF explicit model selection filters requested model" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.pdbx_PDB_model_num
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA A 1 1 10.0 20.0 30.0
        \\2 C CA GLY B 2 2 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    parser.model_num = 2;
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectEqualStrings("B", input.chain_id.?[0].slice());
}

test "parse mmCIF handles quoted atom_site values" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C "C A" "M ET" "A#1" 10.0 20.0 30.0
        \\2 N 'N''X' GLY B 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectEqualStrings("C A", input.atom_name.?[0].slice());
    try std.testing.expectEqualStrings("M ET", input.residue.?[0].slice());
    try std.testing.expectEqualStrings("A#1", input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("N''X", input.atom_name.?[1].slice());
}

test "parse mmCIF skips semicolon text in ignored atom_site column" {
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
        \\_atom_site.pdbx_description
        \\1 C CA ALA 10.0 20.0 30.0
        \\;
        \\ignored atom-site text
        \\still ignored
        \\;
        \\2 N N GLY 11.0 21.0 31.0 plain
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 11.0), input.x[1], 0.001);
}

test "parse mmCIF atom_site stops at new data block" {
    const source =
        \\data_ONE
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C 10.0 20.0 30.0
        \\data_TWO
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\2 N 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
}

test "parse mmCIF null coordinate is invalid" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C . 20.0 30.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    const result = parser.parse(source);
    try std.testing.expectError(ParseError.InvalidCoordinate, result);
}

test "parse mmCIF model filter keeps null model values" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_asym_id
        \\_atom_site.pdbx_PDB_model_num
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C A ? 10.0 20.0 30.0
        \\2 C B 2 11.0 21.0 31.0
        \\3 C C 3 12.0 22.0 32.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    parser.model_num = 2;
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectEqualStrings("A", input.chain_id.?[0].slice());
    try std.testing.expectEqualStrings("B", input.chain_id.?[1].slice());
}

test "parse mmCIF auth chain filter uses auth_asym_id" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_asym_id
        \\_atom_site.auth_asym_id
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C L1 A 10.0 20.0 30.0
        \\2 C L2 B 11.0 21.0 31.0
        \\#
    ;
    const chains = [_][]const u8{"A"};

    var parser = MmcifParser.init(std.testing.allocator);
    parser.use_auth_chain = true;
    parser.chain_filter = chains[0..];
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectEqualStrings("A", input.chain_id.?[0].slice());
}

test "parse mmCIF altLoc selection keeps highest occupancy without A" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.label_alt_id
        \\_atom_site.occupancy
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CB ALA A 1 B 0.40 10.0 20.0 30.0
        \\2 C CB ALA A 1 C 0.70 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 11.0), input.x[0], 0.001);
}

test "parse mmCIF null altLoc is treated as blank" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.label_alt_id
        \\_atom_site.occupancy
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CB ALA A 1 ? 1.00 10.0 20.0 30.0
        \\2 C CB ALA A 1 B 0.50 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 1), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
}

test "parse mmCIF altLoc selection is per atom site and keeps later B-only sites" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.label_alt_id
        \\_atom_site.occupancy
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA A 1 A 0.60 10.0 20.0 30.0
        \\2 C CA ALA A 1 B 0.40 12.0 22.0 32.0
        \\3 C CA GLY A 2 B 0.50 14.0 24.0 34.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 14.0), input.x[1], 0.001);
    try std.testing.expectEqual(@as(i32, 2), input.residue_num.?[1]);
}

test "parse mmCIF altLoc selection is scoped by model" {
    const source =
        \\data_TEST
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.label_asym_id
        \\_atom_site.label_seq_id
        \\_atom_site.label_alt_id
        \\_atom_site.occupancy
        \\_atom_site.pdbx_PDB_model_num
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C CA ALA A 1 A 0.60 1 10.0 20.0 30.0
        \\2 C CA ALA A 1 B 0.40 1 12.0 22.0 32.0
        \\3 C CA ALA A 1 B 0.50 2 14.0 24.0 34.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);
    try std.testing.expectApproxEqAbs(@as(f64, 14.0), input.x[1], 0.001);
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
        fn testOne(_: void, smith: *std.testing.Smith) !void {
            const input = smith.in orelse return;
            var parser = MmcifParser.init(std.testing.allocator);
            defer parser.deinitCcd();
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

test "parse mmCIF with inline CCD data" {
    // A mmCIF file containing both _chem_comp_atom/_chem_comp_bond AND _atom_site loops
    const source =
        \\data_TEST
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.pdbx_aromatic_flag
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\ALA N   N N N
        \\ALA CA  C N N
        \\ALA C   C N N
        \\ALA O   O N N
        \\ALA CB  C N N
        \\ALA OXT O N Y
        \\#
        \\loop_
        \\_chem_comp_bond.comp_id
        \\_chem_comp_bond.atom_id_1
        \\_chem_comp_bond.atom_id_2
        \\_chem_comp_bond.value_order
        \\ALA N   CA  SING
        \\ALA CA  C   SING
        \\ALA CA  CB  SING
        \\ALA C   O   DOUB
        \\ALA C   OXT SING
        \\#
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.group_PDB
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 N N  ALA ATOM 10.0 20.0 30.0
        \\2 C CA ALA ATOM 11.0 21.0 31.0
        \\3 C C  ALA ATOM 12.0 22.0 32.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    defer parser.deinitCcd();
    var input = try parser.parse(source);
    defer input.deinit();

    // Verify atom_site parsing works as normal
    try std.testing.expectEqual(@as(usize, 3), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);

    // Verify inline CCD was parsed
    const ccd = parser.getInlineCcd();
    try std.testing.expect(ccd != null);
    try std.testing.expectEqual(@as(usize, 1), ccd.?.count());

    const comp = ccd.?.get("ALA");
    try std.testing.expect(comp != null);
    try std.testing.expectEqual(@as(usize, 6), comp.?.atoms.len);
    try std.testing.expectEqual(@as(usize, 5), comp.?.bonds.len);

    // Verify OXT is a leaving atom
    var found_oxt = false;
    for (comp.?.atoms) |atom| {
        if (std.mem.eql(u8, atom.atomIdSlice(), "OXT")) {
            try std.testing.expect(atom.leaving);
            found_oxt = true;
        }
    }
    try std.testing.expect(found_oxt);
}

test "parse mmCIF without inline CCD — getInlineCcd returns null" {
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
        \\1 C CA ALA 10.0 20.0 30.0
        \\2 N N  ALA 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    defer parser.deinitCcd();
    var input = try parser.parse(source);
    defer input.deinit();

    // No CCD loops present, so inline_ccd should be null
    try std.testing.expect(parser.getInlineCcd() == null);

    // Normal parsing still works
    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
}

test "parse mmCIF with inline CCD after atom_site" {
    // CCD loops come AFTER the _atom_site loop
    const source =
        \\data_TEST
        \\#
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.group_PDB
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 N N  GLY ATOM 10.0 20.0 30.0
        \\2 C CA GLY ATOM 11.0 21.0 31.0
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\GLY N  N
        \\GLY CA C
        \\GLY C  C
        \\GLY O  O
        \\#
        \\loop_
        \\_chem_comp_bond.comp_id
        \\_chem_comp_bond.atom_id_1
        \\_chem_comp_bond.atom_id_2
        \\_chem_comp_bond.value_order
        \\GLY N  CA SING
        \\GLY CA C  SING
        \\GLY C  O  DOUB
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    defer parser.deinitCcd();
    var input = try parser.parse(source);
    defer input.deinit();

    // Atom coordinates parsed correctly
    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expectApproxEqAbs(@as(f64, 10.0), input.x[0], 0.001);

    // CCD data parsed even though it was after atom_site
    const ccd = parser.getInlineCcd();
    try std.testing.expect(ccd != null);
    try std.testing.expectEqual(@as(usize, 1), ccd.?.count());

    const comp = ccd.?.get("GLY");
    try std.testing.expect(comp != null);
    try std.testing.expectEqual(@as(usize, 4), comp.?.atoms.len);
    try std.testing.expectEqual(@as(usize, 3), comp.?.bonds.len);
}

test "parse mmCIF can skip inline CCD extraction" {
    const source =
        \\data_TEST
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\LIG C1 C
        \\LIG O1 O
        \\#
        \\loop_
        \\_atom_site.id
        \\_atom_site.type_symbol
        \\_atom_site.label_atom_id
        \\_atom_site.label_comp_id
        \\_atom_site.group_PDB
        \\_atom_site.Cartn_x
        \\_atom_site.Cartn_y
        \\_atom_site.Cartn_z
        \\1 C C1 LIG ATOM 10.0 20.0 30.0
        \\2 O O1 LIG ATOM 11.0 21.0 31.0
        \\#
    ;

    var parser = MmcifParser.init(std.testing.allocator);
    parser.parse_inline_ccd = false;
    defer parser.deinitCcd();

    var input = try parser.parse(source);
    defer input.deinit();

    try std.testing.expectEqual(@as(usize, 2), input.atomCount());
    try std.testing.expect(parser.getInlineCcd() == null);
}

//! CCD (Chemical Component Dictionary) streaming parser.
//!
//! Parses `_chem_comp_atom` and `_chem_comp_bond` loops from CIF data
//! using the existing `cif_tokenizer.zig`. Returns a `ComponentDict` keyed
//! by comp_id, where each entry owns its atom and bond arrays.

const std = @import("std");
const Allocator = std.mem.Allocator;
const cif = @import("cif_tokenizer.zig");
const Token = cif.Token;
const Tokenizer = cif.Tokenizer;
const hyb = @import("hybridization.zig");
const CompAtom = hyb.CompAtom;
const CompBond = hyb.CompBond;
const Component = hyb.Component;
const BondOrder = hyb.BondOrder;

// =============================================================================
// StoredComponent — owns allocated atom/bond arrays
// =============================================================================

/// A component that owns its atom and bond data (heap-allocated).
pub const StoredComponent = struct {
    comp_id: [5]u8,
    comp_id_len: u3,
    atoms: []CompAtom,
    bonds: []CompBond,
    allocator: Allocator,

    pub fn deinit(self: *StoredComponent) void {
        self.allocator.free(self.atoms);
        self.allocator.free(self.bonds);
    }

    /// Return a non-owning `Component` view.
    pub fn view(self: *const StoredComponent) Component {
        return .{
            .comp_id = self.comp_id,
            .comp_id_len = self.comp_id_len,
            .atoms = self.atoms,
            .bonds = self.bonds,
        };
    }
};

// =============================================================================
// ComponentDict — stores parsed components keyed by comp_id
// =============================================================================

/// Dictionary of parsed CCD components.
pub const ComponentDict = struct {
    components: std.StringHashMap(StoredComponent),
    allocator: Allocator,
    /// Keys that we allocated and must free.
    owned_keys: std.ArrayListUnmanaged([]const u8),

    pub fn init(allocator: Allocator) ComponentDict {
        return .{
            .components = std.StringHashMap(StoredComponent).init(allocator),
            .allocator = allocator,
            .owned_keys = .empty,
        };
    }

    pub fn deinit(self: *ComponentDict) void {
        var it = self.components.iterator();
        while (it.next()) |entry| {
            var stored = entry.value_ptr;
            stored.deinit();
        }
        // Free all owned keys
        for (self.owned_keys.items) |key| {
            self.allocator.free(key);
        }
        self.owned_keys.deinit(self.allocator);
        self.components.deinit();
    }

    /// Look up a component by comp_id string, returning a non-owning view.
    pub fn get(self: *const ComponentDict, comp_id: []const u8) ?Component {
        const stored = self.components.get(comp_id) orelse return null;
        return stored.view();
    }

    /// Number of components stored.
    pub fn count(self: *const ComponentDict) usize {
        return self.components.count();
    }
};

// =============================================================================
// Column index helpers
// =============================================================================

/// Indices of relevant columns within a `_chem_comp_atom` loop.
const AtomColumns = struct {
    comp_id: ?usize = null,
    atom_id: ?usize = null,
    type_symbol: ?usize = null,
    pdbx_aromatic_flag: ?usize = null,
    pdbx_leaving_atom_flag: ?usize = null,
    total: usize = 0,

    fn isValid(self: *const AtomColumns) bool {
        return self.comp_id != null and self.atom_id != null and self.type_symbol != null;
    }
};

/// Indices of relevant columns within a `_chem_comp_bond` loop.
const BondColumns = struct {
    comp_id: ?usize = null,
    atom_id_1: ?usize = null,
    atom_id_2: ?usize = null,
    value_order: ?usize = null,
    pdbx_aromatic_flag: ?usize = null,
    total: usize = 0,

    fn isValid(self: *const BondColumns) bool {
        return self.comp_id != null and self.atom_id_1 != null and self.atom_id_2 != null;
    }
};

// =============================================================================
// Temp builders
// =============================================================================

/// Temporary per-component data accumulated during parsing.
const ComponentBuilder = struct {
    atoms: std.ArrayListUnmanaged(CompAtom),
    bonds: std.ArrayListUnmanaged(CompBond),
    /// Map from atom_id (fixed-size key) -> index in atoms list.
    atom_name_map: std.StringHashMap(u16),

    fn init(allocator: Allocator) ComponentBuilder {
        return .{
            .atoms = .empty,
            .bonds = .empty,
            .atom_name_map = std.StringHashMap(u16).init(allocator),
        };
    }

    fn deinit(self: *ComponentBuilder, allocator: Allocator) void {
        self.atoms.deinit(allocator);
        self.bonds.deinit(allocator);
        // Free owned keys in atom_name_map
        var it = self.atom_name_map.iterator();
        while (it.next()) |entry| {
            allocator.free(entry.key_ptr.*);
        }
        self.atom_name_map.deinit();
    }
};

// =============================================================================
// Public API
// =============================================================================

/// Parse CCD data from a CIF source string.
///
/// If `filter` is non-null, only components whose comp_id appears in the
/// filter list are included in the result. All others are skipped.
///
/// The returned `ComponentDict` owns all allocated data; call `.deinit()`
/// when done.
pub fn parseCcdData(
    allocator: Allocator,
    source: []const u8,
    filter: ?[]const []const u8,
) !ComponentDict {
    var dict = ComponentDict.init(allocator);
    errdefer dict.deinit();

    // Temporary builders keyed by comp_id (duped key).
    var builders = std.StringHashMap(ComponentBuilder).init(allocator);
    defer {
        var bit = builders.iterator();
        while (bit.next()) |entry| {
            allocator.free(entry.key_ptr.*);
            entry.value_ptr.deinit(allocator);
        }
        builders.deinit();
    }

    var tokenizer = Tokenizer.init(source);

    while (true) {
        const tok = tokenizer.next();
        switch (tok) {
            .eof => break,
            .loop => {
                // Peek at the first tag to determine which loop this is.
                const first_tag_tok = tokenizer.next();
                switch (first_tag_tok) {
                    .tag => |tag_str| {
                        if (isAtomTag(tag_str)) {
                            try parseAtomLoop(allocator, &tokenizer, first_tag_tok, &builders, filter);
                        } else if (isBondTag(tag_str)) {
                            try parseBondLoop(allocator, &tokenizer, first_tag_tok, &builders, filter);
                        }
                        // else: unrelated loop, skip tokens until next structural token
                        // (the parse*Loop functions consume the loop; for unrelated loops
                        // we just let the outer loop continue — the tokenizer will
                        // naturally skip values until the next data_block/loop/eof).
                    },
                    .eof => break,
                    else => {},
                }
            },
            else => {},
        }
    }

    // Convert builders into StoredComponents and insert into dict.
    var bit = builders.iterator();
    while (bit.next()) |entry| {
        const comp_id_str = entry.key_ptr.*;
        const builder = entry.value_ptr;

        var stored = StoredComponent{
            .comp_id = .{ 0, 0, 0, 0, 0 },
            .comp_id_len = 0,
            .atoms = try builder.atoms.toOwnedSlice(allocator),
            .bonds = try builder.bonds.toOwnedSlice(allocator),
            .allocator = allocator,
        };

        const cid_len: usize = @min(comp_id_str.len, 5);
        stored.comp_id_len = @intCast(cid_len);
        for (comp_id_str[0..cid_len], 0..) |c, i| {
            stored.comp_id[i] = c;
        }

        // Dupe the key for the dict
        const dict_key = try allocator.dupe(u8, comp_id_str);
        dict.owned_keys.append(allocator, dict_key) catch |e| {
            allocator.free(dict_key);
            return e;
        };
        try dict.components.put(dict_key, stored);
    }

    return dict;
}

// =============================================================================
// Tag classification
// =============================================================================

fn isAtomTag(tag: []const u8) bool {
    return std.mem.startsWith(u8, tag, "_chem_comp_atom.");
}

fn isBondTag(tag: []const u8) bool {
    return std.mem.startsWith(u8, tag, "_chem_comp_bond.");
}

// =============================================================================
// Atom loop parsing
// =============================================================================

fn parseAtomLoop(
    allocator: Allocator,
    tokenizer: *Tokenizer,
    first_tag: Token,
    builders: *std.StringHashMap(ComponentBuilder),
    filter: ?[]const []const u8,
) !void {
    var cols = AtomColumns{};

    // Process the first tag we already consumed.
    classifyAtomTag(first_tag.tag, 0, &cols);
    cols.total = 1;

    // Read remaining tags.
    while (true) {
        // Save position in case we need to "unread" the token.
        const saved_pos = tokenizer.pos;
        const saved_line = tokenizer.line;
        const saved_col = tokenizer.col;

        const tok = tokenizer.next();
        switch (tok) {
            .tag => |tag_str| {
                classifyAtomTag(tag_str, cols.total, &cols);
                cols.total += 1;
            },
            else => {
                // Not a tag — this is the first value. Restore position
                // so the value-reading loop picks it up.
                tokenizer.pos = saved_pos;
                tokenizer.line = saved_line;
                tokenizer.col = saved_col;
                break;
            },
        }
    }

    if (!cols.isValid()) return;

    // Read row values.
    var col_idx: usize = 0;
    var row_comp_id: ?[]const u8 = null;
    var row_atom_id: ?[]const u8 = null;
    var row_type_symbol: ?[]const u8 = null;
    var row_aromatic: bool = false;
    var row_leaving: bool = false;

    while (true) {
        const saved_pos = tokenizer.pos;
        const saved_line = tokenizer.line;
        const saved_col = tokenizer.col;

        const tok = tokenizer.next();
        switch (tok) {
            .value => |val| {
                if (col_idx == cols.comp_id.?) {
                    row_comp_id = val;
                } else if (col_idx == cols.atom_id.?) {
                    row_atom_id = val;
                } else if (col_idx == cols.type_symbol.?) {
                    row_type_symbol = val;
                } else if (cols.pdbx_aromatic_flag != null and col_idx == cols.pdbx_aromatic_flag.?) {
                    row_aromatic = (val.len == 1 and (val[0] == 'Y' or val[0] == 'y'));
                } else if (cols.pdbx_leaving_atom_flag != null and col_idx == cols.pdbx_leaving_atom_flag.?) {
                    row_leaving = (val.len == 1 and (val[0] == 'Y' or val[0] == 'y'));
                }

                col_idx += 1;
                if (col_idx >= cols.total) {
                    // End of row — emit atom.
                    if (row_comp_id) |cid| {
                        if (row_atom_id) |aid| {
                            if (row_type_symbol) |ts| {
                                if (!isFiltered(cid, filter)) {
                                    const builder = try getOrCreateBuilder(allocator, builders, cid);
                                    var atom = CompAtom.init(aid, ts);
                                    atom.aromatic = row_aromatic;
                                    atom.leaving = row_leaving;

                                    const atom_idx: u16 = @intCast(builder.atoms.items.len);
                                    try builder.atoms.append(allocator, atom);

                                    // Store atom name -> index mapping for bond resolution.
                                    const name_key = try allocator.dupe(u8, aid);
                                    builder.atom_name_map.put(name_key, atom_idx) catch |e| {
                                        allocator.free(name_key);
                                        return e;
                                    };
                                }
                            }
                        }
                    }

                    // Reset for next row.
                    col_idx = 0;
                    row_comp_id = null;
                    row_atom_id = null;
                    row_type_symbol = null;
                    row_aromatic = false;
                    row_leaving = false;
                }
            },
            .loop, .data_block => {
                // New structural token — put it back and return.
                tokenizer.pos = saved_pos;
                tokenizer.line = saved_line;
                tokenizer.col = saved_col;
                return;
            },
            .tag => {
                // A tag outside a loop signals a key-value pair at block level.
                tokenizer.pos = saved_pos;
                tokenizer.line = saved_line;
                tokenizer.col = saved_col;
                return;
            },
            .eof => return,
            else => {},
        }
    }
}

fn classifyAtomTag(tag: []const u8, idx: usize, cols: *AtomColumns) void {
    const prefix = "_chem_comp_atom.";
    if (tag.len <= prefix.len) return;
    const field = tag[prefix.len..];

    if (std.mem.eql(u8, field, "comp_id")) {
        cols.comp_id = idx;
    } else if (std.mem.eql(u8, field, "atom_id")) {
        cols.atom_id = idx;
    } else if (std.mem.eql(u8, field, "type_symbol")) {
        cols.type_symbol = idx;
    } else if (std.mem.eql(u8, field, "pdbx_aromatic_flag")) {
        cols.pdbx_aromatic_flag = idx;
    } else if (std.mem.eql(u8, field, "pdbx_leaving_atom_flag")) {
        cols.pdbx_leaving_atom_flag = idx;
    }
}

// =============================================================================
// Bond loop parsing
// =============================================================================

fn parseBondLoop(
    allocator: Allocator,
    tokenizer: *Tokenizer,
    first_tag: Token,
    builders: *std.StringHashMap(ComponentBuilder),
    filter: ?[]const []const u8,
) !void {
    var cols = BondColumns{};

    classifyBondTag(first_tag.tag, 0, &cols);
    cols.total = 1;

    // Read remaining tags.
    while (true) {
        const saved_pos = tokenizer.pos;
        const saved_line = tokenizer.line;
        const saved_col = tokenizer.col;

        const tok = tokenizer.next();
        switch (tok) {
            .tag => |tag_str| {
                classifyBondTag(tag_str, cols.total, &cols);
                cols.total += 1;
            },
            else => {
                tokenizer.pos = saved_pos;
                tokenizer.line = saved_line;
                tokenizer.col = saved_col;
                break;
            },
        }
    }

    if (!cols.isValid()) return;

    // Read row values.
    var col_idx: usize = 0;
    var row_comp_id: ?[]const u8 = null;
    var row_atom_id_1: ?[]const u8 = null;
    var row_atom_id_2: ?[]const u8 = null;
    var row_value_order: ?[]const u8 = null;
    var row_aromatic: bool = false;

    while (true) {
        const saved_pos = tokenizer.pos;
        const saved_line = tokenizer.line;
        const saved_col = tokenizer.col;

        const tok = tokenizer.next();
        switch (tok) {
            .value => |val| {
                if (col_idx == cols.comp_id.?) {
                    row_comp_id = val;
                } else if (col_idx == cols.atom_id_1.?) {
                    row_atom_id_1 = val;
                } else if (col_idx == cols.atom_id_2.?) {
                    row_atom_id_2 = val;
                } else if (cols.value_order != null and col_idx == cols.value_order.?) {
                    row_value_order = val;
                } else if (cols.pdbx_aromatic_flag != null and col_idx == cols.pdbx_aromatic_flag.?) {
                    row_aromatic = (val.len == 1 and (val[0] == 'Y' or val[0] == 'y'));
                }

                col_idx += 1;
                if (col_idx >= cols.total) {
                    // End of row — emit bond.
                    if (row_comp_id) |cid| {
                        if (row_atom_id_1) |aid1| {
                            if (row_atom_id_2) |aid2| {
                                if (!isFiltered(cid, filter)) {
                                    const builder = try getOrCreateBuilder(allocator, builders, cid);

                                    // Resolve atom names to indices.
                                    const idx1 = builder.atom_name_map.get(aid1);
                                    const idx2 = builder.atom_name_map.get(aid2);

                                    if (idx1 != null and idx2 != null) {
                                        const order = if (row_value_order) |vo|
                                            BondOrder.fromString(vo)
                                        else
                                            BondOrder.unknown;

                                        try builder.bonds.append(allocator, .{
                                            .atom_idx_1 = idx1.?,
                                            .atom_idx_2 = idx2.?,
                                            .order = order,
                                            .aromatic = row_aromatic,
                                        });
                                    }
                                }
                            }
                        }
                    }

                    // Reset for next row.
                    col_idx = 0;
                    row_comp_id = null;
                    row_atom_id_1 = null;
                    row_atom_id_2 = null;
                    row_value_order = null;
                    row_aromatic = false;
                }
            },
            .loop, .data_block => {
                tokenizer.pos = saved_pos;
                tokenizer.line = saved_line;
                tokenizer.col = saved_col;
                return;
            },
            .tag => {
                tokenizer.pos = saved_pos;
                tokenizer.line = saved_line;
                tokenizer.col = saved_col;
                return;
            },
            .eof => return,
            else => {},
        }
    }
}

fn classifyBondTag(tag: []const u8, idx: usize, cols: *BondColumns) void {
    const prefix = "_chem_comp_bond.";
    if (tag.len <= prefix.len) return;
    const field = tag[prefix.len..];

    if (std.mem.eql(u8, field, "comp_id")) {
        cols.comp_id = idx;
    } else if (std.mem.eql(u8, field, "atom_id_1")) {
        cols.atom_id_1 = idx;
    } else if (std.mem.eql(u8, field, "atom_id_2")) {
        cols.atom_id_2 = idx;
    } else if (std.mem.eql(u8, field, "value_order")) {
        cols.value_order = idx;
    } else if (std.mem.eql(u8, field, "pdbx_aromatic_flag")) {
        cols.pdbx_aromatic_flag = idx;
    }
}

// =============================================================================
// Helpers
// =============================================================================

/// Return true if `comp_id` should be SKIPPED (not in filter list).
fn isFiltered(comp_id: []const u8, filter: ?[]const []const u8) bool {
    const f = filter orelse return false; // no filter => include everything
    for (f) |allowed| {
        if (std.mem.eql(u8, comp_id, allowed)) return false; // found in filter => include
    }
    return true; // not found in filter => skip
}

/// Get or create a ComponentBuilder for a given comp_id.
fn getOrCreateBuilder(
    allocator: Allocator,
    builders: *std.StringHashMap(ComponentBuilder),
    comp_id: []const u8,
) !*ComponentBuilder {
    const gop = try builders.getOrPut(comp_id);
    if (!gop.found_existing) {
        // Dupe the key so the builder map owns it.
        const key_copy = try allocator.dupe(u8, comp_id);
        gop.key_ptr.* = key_copy;
        gop.value_ptr.* = ComponentBuilder.init(allocator);
    }
    return gop.value_ptr;
}

// =============================================================================
// Tests
// =============================================================================

test "parseCcdData — single ALA component" {
    const source =
        \\data_ALA
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
        \\ALA H   H N N
        \\ALA HA  H N N
        \\#
        \\loop_
        \\_chem_comp_bond.comp_id
        \\_chem_comp_bond.atom_id_1
        \\_chem_comp_bond.atom_id_2
        \\_chem_comp_bond.value_order
        \\_chem_comp_bond.pdbx_aromatic_flag
        \\ALA N   CA  SING N
        \\ALA CA  C   SING N
        \\ALA CA  CB  SING N
        \\ALA CA  HA  SING N
        \\ALA C   O   DOUB N
        \\ALA C   OXT SING N
        \\ALA N   H   SING N
        \\#
    ;

    var dict = try parseCcdData(std.testing.allocator, source, null);
    defer dict.deinit();

    try std.testing.expectEqual(@as(usize, 1), dict.count());

    const comp = dict.get("ALA") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqualStrings("ALA", comp.compIdSlice());

    // 8 atoms: N, CA, C, O, CB, OXT, H, HA
    try std.testing.expectEqual(@as(usize, 8), comp.atoms.len);

    // 7 bonds
    try std.testing.expectEqual(@as(usize, 7), comp.bonds.len);

    // Verify C=O is double bond. Find the bond where atom names are C and O.
    // C is index 2, O is index 3 in our atom list.
    var found_co_double = false;
    for (comp.bonds) |bond| {
        const a1 = comp.atoms[bond.atom_idx_1].atomIdSlice();
        const a2 = comp.atoms[bond.atom_idx_2].atomIdSlice();
        if (std.mem.eql(u8, a1, "C") and std.mem.eql(u8, a2, "O")) {
            try std.testing.expectEqual(BondOrder.double, bond.order);
            found_co_double = true;
        }
    }
    try std.testing.expect(found_co_double);

    // Verify OXT is a leaving atom
    var found_oxt_leaving = false;
    for (comp.atoms) |atom| {
        if (std.mem.eql(u8, atom.atomIdSlice(), "OXT")) {
            try std.testing.expect(atom.leaving);
            found_oxt_leaving = true;
        }
    }
    try std.testing.expect(found_oxt_leaving);
}

test "parseCcdData — filter by comp_id" {
    const source =
        \\data_ALA
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\ALA N  N
        \\ALA CA C
        \\ALA C  C
        \\ALA O  O
        \\ALA CB C
        \\#
        \\data_GLY
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
    ;

    const filter = &[_][]const u8{"GLY"};
    var dict = try parseCcdData(std.testing.allocator, source, filter);
    defer dict.deinit();

    // Only GLY should be present.
    try std.testing.expectEqual(@as(usize, 1), dict.count());
    try std.testing.expect(dict.get("GLY") != null);
    try std.testing.expect(dict.get("ALA") == null);

    const gly = dict.get("GLY").?;
    try std.testing.expectEqual(@as(usize, 4), gly.atoms.len);
}

test "parseCcdData — minimal columns (HOH, no aromatic/leaving flags, no bonds)" {
    const source =
        \\data_HOH
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\HOH O  O
        \\HOH H1 H
        \\HOH H2 H
        \\#
    ;

    var dict = try parseCcdData(std.testing.allocator, source, null);
    defer dict.deinit();

    try std.testing.expectEqual(@as(usize, 1), dict.count());

    const comp = dict.get("HOH") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 3), comp.atoms.len);
    try std.testing.expectEqual(@as(usize, 0), comp.bonds.len);

    // Aromatic/leaving should default to false.
    for (comp.atoms) |atom| {
        try std.testing.expect(!atom.aromatic);
        try std.testing.expect(!atom.leaving);
    }
}

test "parseCcdData — aromatic flags" {
    const source =
        \\data_PHE
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.pdbx_aromatic_flag
        \\PHE CB  C N
        \\PHE CG  C Y
        \\PHE CD1 C Y
        \\PHE CD2 C Y
        \\PHE CE1 C Y
        \\PHE CE2 C Y
        \\PHE CZ  C Y
        \\#
        \\loop_
        \\_chem_comp_bond.comp_id
        \\_chem_comp_bond.atom_id_1
        \\_chem_comp_bond.atom_id_2
        \\_chem_comp_bond.value_order
        \\_chem_comp_bond.pdbx_aromatic_flag
        \\PHE CB  CG  SING N
        \\PHE CG  CD1 DOUB Y
        \\PHE CG  CD2 SING Y
        \\PHE CD1 CE1 SING Y
        \\PHE CD2 CE2 DOUB Y
        \\PHE CE1 CZ  DOUB Y
        \\PHE CE2 CZ  SING Y
        \\#
    ;

    var dict = try parseCcdData(std.testing.allocator, source, null);
    defer dict.deinit();

    const comp = dict.get("PHE") orelse return error.TestUnexpectedResult;

    // CB should NOT be aromatic; CG, CD1, CD2, CE1, CE2, CZ should be aromatic.
    for (comp.atoms) |atom| {
        const name = atom.atomIdSlice();
        if (std.mem.eql(u8, name, "CB")) {
            try std.testing.expect(!atom.aromatic);
        } else {
            try std.testing.expect(atom.aromatic);
        }
    }

    // Bond CG-CD1 should be aromatic.
    var found_cg_cd1 = false;
    for (comp.bonds) |bond| {
        const a1 = comp.atoms[bond.atom_idx_1].atomIdSlice();
        const a2 = comp.atoms[bond.atom_idx_2].atomIdSlice();
        if (std.mem.eql(u8, a1, "CG") and std.mem.eql(u8, a2, "CD1")) {
            try std.testing.expect(bond.aromatic);
            try std.testing.expectEqual(BondOrder.double, bond.order);
            found_cg_cd1 = true;
        }
    }
    try std.testing.expect(found_cg_cd1);

    // Bond CB-CG should NOT be aromatic.
    var found_cb_cg = false;
    for (comp.bonds) |bond| {
        const a1 = comp.atoms[bond.atom_idx_1].atomIdSlice();
        const a2 = comp.atoms[bond.atom_idx_2].atomIdSlice();
        if (std.mem.eql(u8, a1, "CB") and std.mem.eql(u8, a2, "CG")) {
            try std.testing.expect(!bond.aromatic);
            try std.testing.expectEqual(BondOrder.single, bond.order);
            found_cb_cg = true;
        }
    }
    try std.testing.expect(found_cb_cg);
}

test "parseCcdData — extra columns are skipped" {
    // Real CCD files have many more columns. Verify we handle extra columns.
    const source =
        \\data_ALA
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.alt_atom_id
        \\_chem_comp_atom.type_symbol
        \\_chem_comp_atom.charge
        \\_chem_comp_atom.pdbx_align
        \\_chem_comp_atom.pdbx_aromatic_flag
        \\_chem_comp_atom.pdbx_leaving_atom_flag
        \\_chem_comp_atom.pdbx_stereo_config
        \\_chem_comp_atom.model_Cartn_x
        \\_chem_comp_atom.model_Cartn_y
        \\_chem_comp_atom.model_Cartn_z
        \\ALA N   N   N 0 1 N N N 0.000 0.000 0.000
        \\ALA CA  CA  C 0 1 N N S 1.458 0.000 0.000
        \\ALA C   C   C 0 1 N N N 2.009 1.420 0.000
        \\#
    ;

    var dict = try parseCcdData(std.testing.allocator, source, null);
    defer dict.deinit();

    const comp = dict.get("ALA") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 3), comp.atoms.len);

    // Verify atom data was correctly extracted despite extra columns.
    try std.testing.expectEqualStrings("N", comp.atoms[0].atomIdSlice());
    try std.testing.expectEqualStrings("N", comp.atoms[0].typeSymbolSlice());
    try std.testing.expectEqualStrings("CA", comp.atoms[1].atomIdSlice());
    try std.testing.expectEqualStrings("C", comp.atoms[1].typeSymbolSlice());
    try std.testing.expectEqualStrings("C", comp.atoms[2].atomIdSlice());
    try std.testing.expectEqualStrings("C", comp.atoms[2].typeSymbolSlice());
}

test "parseCcdData — multiple comp_ids in single loop" {
    // CCD files may have all components in a single loop.
    const source =
        \\data_all
        \\#
        \\loop_
        \\_chem_comp_atom.comp_id
        \\_chem_comp_atom.atom_id
        \\_chem_comp_atom.type_symbol
        \\ALA N  N
        \\ALA CA C
        \\GLY N  N
        \\GLY CA C
        \\GLY C  C
        \\#
    ;

    var dict = try parseCcdData(std.testing.allocator, source, null);
    defer dict.deinit();

    try std.testing.expectEqual(@as(usize, 2), dict.count());

    const ala = dict.get("ALA") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 2), ala.atoms.len);

    const gly = dict.get("GLY") orelse return error.TestUnexpectedResult;
    try std.testing.expectEqual(@as(usize, 3), gly.atoms.len);
}

test "ComponentDict — empty source" {
    var dict = try parseCcdData(std.testing.allocator, "", null);
    defer dict.deinit();
    try std.testing.expectEqual(@as(usize, 0), dict.count());
}

//! Alternate-location (altLoc) handling shared by the PDB, mmCIF and
//! BinaryCIF parsers: the `--altloc` setting and the rules that decide which
//! alternates survive.
//!
//! ## Terms
//!
//! - An *alternate* is an atom with a non-blank altLoc ID.
//! - A *residue position* is a model, chain, residue number and insertion
//!   code.
//! - An *atom site* is an atom name within a residue name at a residue
//!   position.
//!
//! ## Rules
//!
//! Atoms without an altLoc ID are never dropped. `all` keeps every alternate
//! too, and `none` expects the caller to reject alternates while it reads
//! them. The other modes decide in two steps.
//!
//! 1. One residue per residue position. With microheterogeneity the
//!    alternates of a position belong to different residues (SER as altLoc A
//!    and PRO as altLoc B), and keeping an alternate of each atom name would
//!    superimpose them. So when the alternates of a position carry more than
//!    one residue name, one residue survives and the alternates of the others
//!    are dropped:
//!    - `auto`: the residue of the first alternate A. Without an alternate A,
//!      the residue with the highest occupancy.
//!    - `highest_occupancy`: the residue with the highest occupancy.
//!    - `selected`: only alternates with the selected ID count, so the
//!      residue that carries the ID survives. If several do, the first one.
//!
//!    The occupancy of a residue is the occupancy of its first alternate in
//!    the file.
//!
//! 2. One alternate per atom site, among the alternates of the surviving
//!    residue:
//!    - `auto`: whichever comes first in the file, an atom without an altLoc
//!      or an alternate A, decides. The first drops every alternate of the
//!      site, the second keeps the alternates A. A site with neither keeps
//!      the alternate with the highest occupancy.
//!    - `highest_occupancy`: the alternate with the highest occupancy is
//!      kept, unless the site also has an atom without an altLoc, which drops
//!      every alternate.
//!    - `selected`: the alternates with the selected ID are kept.
//!
//! Equal occupancies are a tie, and the residue or alternate that comes first
//! in the file wins it.

const std = @import("std");
const Allocator = std.mem.Allocator;

pub const AltLocMode = enum {
    /// Historical behavior: blank altLoc wins, then A, then highest occupancy.
    auto,
    /// Assume there are no non-blank altLocs; fail fast if one is encountered.
    none,
    /// Keep all alternate locations.
    all,
    /// Keep blank altLocs and the selected altLoc ID.
    selected,
    /// Keep the highest-occupancy residue for each residue position and the
    /// highest-occupancy alternate for each of its atom sites.
    highest_occupancy,
};

pub const AltLocSetting = struct {
    mode: AltLocMode = .auto,
    id: u8 = 'A',
};

pub fn parseSetting(value: []const u8) ?AltLocSetting {
    if (std.mem.eql(u8, value, "auto")) return .{ .mode = .auto };
    if (std.mem.eql(u8, value, "none")) return .{ .mode = .none };
    if (std.mem.eql(u8, value, "all")) return .{ .mode = .all };
    if (std.mem.eql(u8, value, "highest-occupancy") or
        std.mem.eql(u8, value, "highest_occupancy"))
    {
        return .{ .mode = .highest_occupancy };
    }
    if (value.len == 1 and value[0] != '.' and value[0] != '?') {
        return .{ .mode = .selected, .id = value[0] };
    }
    return null;
}

/// What altLoc resolution reads from one atom record. The slices are borrowed
/// from the record.
pub const Site = struct {
    model_num: ?u32,
    chain_id: []const u8,
    /// Residue number that identifies the residue within its chain.
    seq: i32,
    /// True when `seq` belongs to another numbering than the `seq` of rows
    /// where this is false (mmCIF and BinaryCIF number non-polymer rows by
    /// auth_seq_id and polymer rows by label_seq_id). Equal numbers of the
    /// two numberings are different residue positions.
    seq_is_auth: bool = false,
    insertion_code: []const u8,
    residue: []const u8,
    atom_name: []const u8,
    /// ' ' for an atom without an altLoc ID.
    alt_loc: u8,
    occupancy: f64,
};

/// Decide which atom records survive altLoc resolution (see the rules in the
/// module documentation).
///
/// `Record` is the atom record type of a parser and has a method
/// `pub fn altLocSite(self: Record) Site`. `records` are in file order.
///
/// Returns one flag per record, true for the records to keep, or null when
/// every record is kept. The caller owns the returned slice.
///
/// Runs in time linear in the number of records: residue positions, residues
/// and atom sites are found through hash maps, and their records do not have
/// to be contiguous.
pub fn resolve(
    comptime Record: type,
    allocator: Allocator,
    records: []const Record,
    setting: AltLocSetting,
) Allocator.Error!?[]bool {
    switch (setting.mode) {
        .all, .none => return null,
        .auto, .selected, .highest_occupancy => {},
    }

    const has_alternates = for (records) |record| {
        if (record.altLocSite().alt_loc != ' ') break true;
    } else false;
    if (!has_alternates) return null;

    const keep = try allocator.alloc(bool, records.len);
    errdefer allocator.free(keep);
    @memset(keep, true);

    var positions = PositionIndex{};
    defer positions.deinit(allocator);
    var position_states = std.ArrayListUnmanaged(PositionState).empty;
    defer position_states.deinit(allocator);
    var residues = ChildMap.empty;
    defer residues.deinit(allocator);
    var sites = ChildMap.empty;
    defer sites.deinit(allocator);
    var site_states = std.ArrayListUnmanaged(SiteState).empty;
    defer site_states.deinit(allocator);

    // Per record: the residue of an alternate after pass 1, and the atom
    // site of an alternate after pass 2. Not set for atoms without an altLoc.
    const group_of = try allocator.alloc(u32, records.len);
    defer allocator.free(group_of);

    // Pass 1: the residue positions and residues that hold alternates, and
    // the residue that survives at each position.
    for (records, group_of, keep) |record, *group, *flag| {
        const site = record.altLocSite();
        if (site.alt_loc == ' ') continue;
        if (setting.mode == .selected and site.alt_loc != setting.id) {
            flag.* = false;
            continue;
        }

        const position = try positions.getOrPut(allocator, PositionKey.of(site));
        if (position == position_states.items.len) try position_states.append(allocator, .{});
        const n_residues = residues.count();
        const residue = try getOrPutChild(allocator, &residues, .{ .parent = position, .name = site.residue });
        position_states.items[position].add(setting.mode, site, residue, residue == n_residues);
        group.* = residue;
    }

    const survives = try allocator.alloc(bool, residues.count());
    defer allocator.free(survives);
    @memset(survives, false);
    for (position_states.items) |state| survives[state.residue] = true;

    if (setting.mode == .selected) {
        for (records, group_of, keep) |record, residue, *flag| {
            if (record.altLocSite().alt_loc == ' ' or !flag.*) continue;
            flag.* = survives[residue];
        }
        return keep;
    }

    // Pass 2: what each atom site of a surviving residue holds. Atoms
    // without an altLoc take part when their residue has alternates, and
    // they can come before them.
    for (records, group_of, keep, 0..) |record, *group, *flag, i| {
        const site = record.altLocSite();
        const residue = if (site.alt_loc != ' ') blk: {
            if (!survives[group.*]) {
                flag.* = false;
                continue;
            }
            break :blk group.*;
        } else blk: {
            const position = positions.get(PositionKey.of(site)) orelse continue;
            break :blk residues.get(.{ .parent = position, .name = site.residue }) orelse continue;
        };

        const site_index = try getOrPutChild(allocator, &sites, .{ .parent = residue, .name = site.atom_name });
        if (site_index == site_states.items.len) try site_states.append(allocator, .{});
        site_states.items[site_index].add(setting.mode, site, @intCast(i));
        group.* = site_index;
    }

    // Pass 3: keep the alternate that its atom site settled on.
    for (records, group_of, keep, 0..) |record, site_index, *flag, i| {
        const alt_loc = record.altLocSite().alt_loc;
        if (alt_loc == ' ' or !flag.*) continue;
        flag.* = site_states.items[site_index].keeps(alt_loc, @intCast(i));
    }

    return keep;
}

const none_index = std.math.maxInt(u32);

/// The residue that survives at one residue position.
const PositionState = struct {
    residue: u32 = none_index,
    /// Occupancy of `residue`: that of its first alternate in the file.
    occupancy: f64 = 0.0,
    /// `auto` mode: `residue` holds the first alternate A of the position,
    /// and no other residue displaces it.
    has_alt_a: bool = false,

    /// `site` is an alternate of `residue`, the first one in the file when
    /// `is_first_of_residue`. In `selected` mode it has the selected ID.
    fn add(self: *PositionState, mode: AltLocMode, site: Site, residue: u32, is_first_of_residue: bool) void {
        if (self.has_alt_a) return;
        if (mode == .auto and site.alt_loc == 'A') {
            self.residue = residue;
            self.has_alt_a = true;
            return;
        }
        if (!is_first_of_residue) return;
        if (self.residue == none_index or (mode != .selected and site.occupancy > self.occupancy)) {
            self.residue = residue;
            self.occupancy = site.occupancy;
        }
    }
};

/// What the atoms of one atom site amount to.
const SiteState = struct {
    decider: Decider = .undecided,
    /// Alternate with the highest occupancy, the first one in the file on a
    /// tie. In `auto` mode the alternates A are not candidates.
    best: u32 = none_index,
    best_occupancy: f64 = 0.0,

    const Decider = enum {
        undecided,
        /// The site has an atom without an altLoc, which drops its
        /// alternates. In `auto` mode only when it comes before an
        /// alternate A.
        blank,
        /// `auto` mode: an alternate A came before any atom without an
        /// altLoc.
        alt_a,
    };

    fn add(self: *SiteState, mode: AltLocMode, site: Site, index: u32) void {
        if (site.alt_loc == ' ') {
            if (self.decider == .undecided) self.decider = .blank;
            return;
        }
        if (mode == .auto and site.alt_loc == 'A') {
            if (self.decider == .undecided) self.decider = .alt_a;
            return;
        }
        if (self.best == none_index or site.occupancy > self.best_occupancy) {
            self.best = index;
            self.best_occupancy = site.occupancy;
        }
    }

    fn keeps(self: SiteState, alt_loc: u8, index: u32) bool {
        return switch (self.decider) {
            .blank => false,
            .alt_a => alt_loc == 'A',
            .undecided => self.best == index,
        };
    }
};

const PositionKey = struct {
    model_num: ?u32,
    seq: i32,
    seq_is_auth: bool,
    chain_id: []const u8,
    insertion_code: []const u8,

    fn of(site: Site) PositionKey {
        return .{
            .model_num = site.model_num,
            .seq = site.seq,
            .seq_is_auth = site.seq_is_auth,
            .chain_id = site.chain_id,
            .insertion_code = site.insertion_code,
        };
    }

    fn eql(a: PositionKey, b: PositionKey) bool {
        return a.seq == b.seq and
            a.model_num == b.model_num and
            a.seq_is_auth == b.seq_is_auth and
            std.mem.eql(u8, a.chain_id, b.chain_id) and
            std.mem.eql(u8, a.insertion_code, b.insertion_code);
    }

    fn hash(self: PositionKey) u64 {
        var hasher = std.hash.Wyhash.init(0);
        std.hash.autoHash(&hasher, self.model_num);
        std.hash.autoHash(&hasher, self.seq);
        std.hash.autoHash(&hasher, self.seq_is_auth);
        hasher.update(self.chain_id);
        hasher.update(&.{0}); // separator
        hasher.update(self.insertion_code);
        return hasher.final();
    }
};

const PositionContext = struct {
    pub fn hash(_: PositionContext, key: PositionKey) u64 {
        return key.hash();
    }

    pub fn eql(_: PositionContext, a: PositionKey, b: PositionKey) bool {
        return a.eql(b);
    }
};

/// Numbers the residue positions that hold alternates.
const PositionIndex = struct {
    map: std.HashMapUnmanaged(PositionKey, u32, PositionContext, std.hash_map.default_max_load_percentage) = .empty,
    /// Key and result of the previous lookup. The atoms of a residue are
    /// contiguous in real files, so most lookups repeat the previous one.
    /// This only saves time: positions are found wherever their atoms are.
    last_key: ?PositionKey = null,
    last: ?u32 = null,

    fn deinit(self: *PositionIndex, allocator: Allocator) void {
        self.map.deinit(allocator);
    }

    fn getOrPut(self: *PositionIndex, allocator: Allocator, key: PositionKey) Allocator.Error!u32 {
        if (self.last_key) |last_key| {
            if (self.last != null and last_key.eql(key)) return self.last.?;
        }
        const entry = try self.map.getOrPut(allocator, key);
        if (!entry.found_existing) entry.value_ptr.* = self.map.count() - 1;
        self.last_key = key;
        self.last = entry.value_ptr.*;
        return entry.value_ptr.*;
    }

    fn get(self: *PositionIndex, key: PositionKey) ?u32 {
        if (self.last_key) |last_key| {
            if (last_key.eql(key)) return self.last;
        }
        self.last_key = key;
        self.last = self.map.get(key);
        return self.last;
    }
};

/// A name within a numbered parent: a residue name at a residue position, or
/// an atom name in a residue.
const ChildKey = struct {
    parent: u32,
    name: []const u8,
};

const ChildContext = struct {
    pub fn hash(_: ChildContext, key: ChildKey) u64 {
        return std.hash.Wyhash.hash(key.parent, key.name);
    }

    pub fn eql(_: ChildContext, a: ChildKey, b: ChildKey) bool {
        return a.parent == b.parent and std.mem.eql(u8, a.name, b.name);
    }
};

const ChildMap = std.HashMapUnmanaged(ChildKey, u32, ChildContext, std.hash_map.default_max_load_percentage);

/// Number of `key` in `map`, in order of first appearance.
fn getOrPutChild(allocator: Allocator, map: *ChildMap, key: ChildKey) Allocator.Error!u32 {
    const entry = try map.getOrPut(allocator, key);
    if (!entry.found_existing) entry.value_ptr.* = map.count() - 1;
    return entry.value_ptr.*;
}

// ============================================================================
// Tests
// ============================================================================

const TestRecord = struct {
    site: Site,

    pub fn altLocSite(self: TestRecord) Site {
        return self.site;
    }
};

fn testRecord(seq: i32, residue: []const u8, atom_name: []const u8, alt_loc: u8, occupancy: f64) TestRecord {
    return .{ .site = .{
        .model_num = 1,
        .chain_id = "A",
        .seq = seq,
        .insertion_code = "",
        .residue = residue,
        .atom_name = atom_name,
        .alt_loc = alt_loc,
        .occupancy = occupancy,
    } };
}

fn expectKept(setting: AltLocSetting, records: []const TestRecord, expected: []const bool) !void {
    const keep = try resolve(TestRecord, std.testing.allocator, records, setting);
    defer if (keep) |flags| std.testing.allocator.free(flags);
    try std.testing.expectEqualSlices(bool, expected, keep orelse return error.TestExpectedFlags);
}

fn sameAtomSite(a: Site, b: Site) bool {
    return PositionKey.of(a).eql(PositionKey.of(b)) and
        std.mem.eql(u8, a.residue, b.residue) and
        std.mem.eql(u8, a.atom_name, b.atom_name);
}

/// Name of the residue that survives at the residue position of `atom`,
/// found by a scan over all records.
fn referenceResidue(records: []const TestRecord, setting: AltLocSetting, atom: Site) []const u8 {
    var seen: [8][]const u8 = undefined;
    var n_seen: usize = 0;
    var best: ?Site = null;
    for (records) |record| {
        const other = record.site;
        if (other.alt_loc == ' ' or !PositionKey.of(atom).eql(PositionKey.of(other))) continue;
        if (setting.mode == .selected and other.alt_loc != setting.id) continue;
        if (setting.mode == .auto and other.alt_loc == 'A') return other.residue;

        // Only the first alternate of a residue counts
        const is_first = for (seen[0..n_seen]) |name| {
            if (std.mem.eql(u8, name, other.residue)) break false;
        } else true;
        if (!is_first) continue;
        seen[n_seen] = other.residue;
        n_seen += 1;

        if (best == null or (setting.mode != .selected and other.occupancy > best.?.occupancy)) best = other;
    }
    return best.?.residue;
}

/// The rules of the module documentation, written as scans over all records
/// for every alternate. Quadratic, and independent of the hash maps.
fn referenceKeeps(records: []const TestRecord, setting: AltLocSetting, index: usize) bool {
    const atom = records[index].site;
    if (atom.alt_loc == ' ') return true;
    switch (setting.mode) {
        .all, .none => return true,
        .selected => if (atom.alt_loc != setting.id) return false,
        .auto, .highest_occupancy => {},
    }
    if (!std.mem.eql(u8, atom.residue, referenceResidue(records, setting, atom))) return false;
    if (setting.mode == .selected) return true;

    var best: ?usize = null;
    for (records, 0..) |other_record, other_index| {
        const other = other_record.site;
        if (!sameAtomSite(atom, other)) continue;
        if (other.alt_loc == ' ') return false;
        if (setting.mode == .auto and other.alt_loc == 'A') return atom.alt_loc == 'A';
        if (best == null or other.occupancy > records[best.?].site.occupancy) best = other_index;
    }
    return best == index;
}

test "parseSetting" {
    try std.testing.expectEqual(AltLocMode.auto, parseSetting("auto").?.mode);
    try std.testing.expectEqual(AltLocMode.none, parseSetting("none").?.mode);
    try std.testing.expectEqual(AltLocMode.all, parseSetting("all").?.mode);
    try std.testing.expectEqual(AltLocMode.highest_occupancy, parseSetting("highest-occupancy").?.mode);
    try std.testing.expectEqual(AltLocMode.highest_occupancy, parseSetting("highest_occupancy").?.mode);

    const selected = parseSetting("B").?;
    try std.testing.expectEqual(AltLocMode.selected, selected.mode);
    try std.testing.expectEqual(@as(u8, 'B'), selected.id);

    try std.testing.expect(parseSetting("") == null);
    try std.testing.expect(parseSetting(".") == null);
    try std.testing.expect(parseSetting("?") == null);
    try std.testing.expect(parseSetting("first") == null);
}

test "resolve returns null when nothing is dropped" {
    const allocator = std.testing.allocator;
    const plain = [_]TestRecord{
        testRecord(1, "ALA", "N", ' ', 1.0),
        testRecord(1, "ALA", "CA", ' ', 1.0),
    };
    const alternates = [_]TestRecord{
        testRecord(1, "ALA", "CA", 'A', 0.5),
        testRecord(1, "ALA", "CA", 'B', 0.5),
    };

    for ([_]AltLocMode{ .auto, .none, .all, .selected, .highest_occupancy }) |mode| {
        try std.testing.expect(try resolve(TestRecord, allocator, &plain, .{ .mode = mode }) == null);
    }
    try std.testing.expect(try resolve(TestRecord, allocator, &alternates, .{ .mode = .all }) == null);
    try std.testing.expect(try resolve(TestRecord, allocator, &alternates, .{ .mode = .none }) == null);
    try std.testing.expect(try resolve(TestRecord, allocator, &.{}, .{}) == null);
}

test "resolve auto prefers a blank altLoc, then A, then the highest occupancy" {
    const records = [_]TestRecord{
        // A wins whatever its occupancy and position
        testRecord(1, "ALA", "CA", 'B', 0.9),
        testRecord(1, "ALA", "CA", 'A', 0.1),
        // Without A the highest occupancy wins
        testRecord(1, "ALA", "CB", 'B', 0.3),
        testRecord(1, "ALA", "CB", 'C', 0.7),
        // An atom without an altLoc drops the alternates of its site
        testRecord(1, "ALA", "N", ' ', 1.0),
        testRecord(1, "ALA", "N", 'B', 0.5),
        // A tie goes to the first alternate
        testRecord(1, "ALA", "O", 'C', 0.5),
        testRecord(1, "ALA", "O", 'B', 0.5),
    };
    try expectKept(.{ .mode = .auto }, &records, &.{ false, true, false, true, true, false, true, false });
}

test "resolve highest occupancy keeps exactly one alternate per site" {
    const records = [_]TestRecord{
        testRecord(1, "ALA", "CA", 'A', 0.4),
        testRecord(1, "ALA", "CA", 'B', 0.6),
        // 0.50/0.50 and three-way ties go to the first alternate
        testRecord(1, "ALA", "CB", 'B', 0.5),
        testRecord(1, "ALA", "CB", 'A', 0.5),
        testRecord(1, "ALA", "CG", 'A', 0.33),
        testRecord(1, "ALA", "CG", 'B', 0.33),
        testRecord(1, "ALA", "CG", 'C', 0.33),
        // No occupancy at all
        testRecord(1, "ALA", "CD", 'A', 0.0),
        testRecord(1, "ALA", "CD", 'B', 0.0),
        // An atom without an altLoc is kept, and drops the alternates
        testRecord(1, "ALA", "N", 'A', 0.9),
        testRecord(1, "ALA", "N", ' ', 0.1),
    };
    try expectKept(.{ .mode = .highest_occupancy }, &records, &.{
        false, true,
        true,  false,
        true,  false,
        false, true,
        false, false,
        true,
    });
}

test "resolve selected keeps the requested ID and the atoms without an altLoc" {
    const records = [_]TestRecord{
        testRecord(1, "ALA", "N", ' ', 1.0),
        testRecord(1, "ALA", "CA", 'A', 0.6),
        testRecord(1, "ALA", "CA", 'B', 0.4),
        // A site without the requested ID loses its atom
        testRecord(2, "GLY", "CA", 'A', 1.0),
    };
    try expectKept(.{ .mode = .selected, .id = 'B' }, &records, &.{ true, false, true, false });
}

test "resolve finds the atoms of a site that are not contiguous" {
    // The alternates B are listed after the rest of the chain, and the atom
    // without an altLoc that drops the alternates of N comes last.
    const records = [_]TestRecord{
        testRecord(1, "ALA", "CA", 'A', 0.5),
        testRecord(1, "ALA", "N", 'B', 0.5),
        testRecord(2, "GLY", "CA", 'C', 0.2),
        testRecord(1, "ALA", "CA", 'B', 0.5),
        testRecord(2, "GLY", "CA", 'B', 0.8),
        testRecord(1, "ALA", "N", ' ', 1.0),
    };
    try expectKept(.{ .mode = .auto }, &records, &.{ true, false, false, false, true, true });
    try expectKept(.{ .mode = .highest_occupancy }, &records, &.{ true, false, false, false, true, true });
}

test "resolve keeps every atom without an altLoc that shares a site" {
    // Two chains without chain IDs repeat the residue number
    const records = [_]TestRecord{
        testRecord(1, "ALA", "CA", ' ', 1.0),
        testRecord(1, "ALA", "CB", 'A', 0.5),
        testRecord(1, "ALA", "CB", 'B', 0.5),
        testRecord(1, "ALA", "CA", ' ', 0.8),
    };
    for ([_]AltLocMode{ .auto, .highest_occupancy }) |mode| {
        try expectKept(.{ .mode = mode }, &records, &.{ true, true, false, true });
    }
}

test "resolve keeps one residue where alternates are different residues" {
    // Position 2 is PRO as altLoc A and SER as altLoc B, with a shared N that
    // has no altLoc. Position 3 is ILE as B and VAL as C, with no alternate A.
    const records = [_]TestRecord{
        testRecord(2, "PRO", "N", ' ', 1.0),
        testRecord(2, "PRO", "CA", 'A', 0.4),
        testRecord(2, "PRO", "CD", 'A', 0.4),
        testRecord(2, "SER", "CA", 'B', 0.6),
        testRecord(2, "SER", "OG", 'B', 0.6),
        testRecord(3, "ILE", "CA", 'B', 0.3),
        testRecord(3, "ILE", "CD1", 'B', 0.3),
        testRecord(3, "VAL", "CA", 'C', 0.7),
        testRecord(3, "VAL", "CG1", 'C', 0.7),
    };

    // A, then the highest occupancy
    try expectKept(.{ .mode = .auto }, &records, &.{ true, true, true, false, false, false, false, true, true });
    try expectKept(.{ .mode = .highest_occupancy }, &records, &.{ true, false, false, true, true, false, false, true, true });
    // The residue that carries the ID, or none
    try expectKept(.{ .mode = .selected, .id = 'A' }, &records, &.{ true, true, true, false, false, false, false, false, false });
    try expectKept(.{ .mode = .selected, .id = 'B' }, &records, &.{ true, false, false, true, true, true, true, false, false });
    try expectKept(.{ .mode = .selected, .id = 'C' }, &records, &.{ true, false, false, false, false, false, false, true, true });
}

test "resolve compares residues by their first alternate and keeps the first on a tie" {
    const tie = [_]TestRecord{
        testRecord(1, "LEU", "CA", 'B', 0.5),
        testRecord(1, "LEU", "CB", 'B', 0.5),
        testRecord(1, "ILE", "CA", 'C', 0.5),
        testRecord(1, "ILE", "CB", 'C', 0.5),
    };
    // Only the first alternate of a residue counts: 0.4 against 0.5
    const first_alternate = [_]TestRecord{
        testRecord(1, "LEU", "CA", 'B', 0.4),
        testRecord(1, "LEU", "CB", 'B', 0.9),
        testRecord(1, "ILE", "CA", 'C', 0.5),
        testRecord(1, "ILE", "CB", 'C', 0.1),
    };
    for ([_]AltLocMode{ .auto, .highest_occupancy }) |mode| {
        try expectKept(.{ .mode = mode }, &tie, &.{ true, true, false, false });
        try expectKept(.{ .mode = mode }, &first_alternate, &.{ false, false, true, true });
    }
}

test "resolve keeps one residue whose atoms are not contiguous" {
    // The atoms of the two residues alternate, as in PDB-format files, and
    // the alternate A that decides for PRO comes last
    const records = [_]TestRecord{
        testRecord(2, "SER", "CA", 'B', 0.6),
        testRecord(2, "PRO", "CA", 'C', 0.2),
        testRecord(2, "SER", "CB", 'B', 0.6),
        testRecord(2, "PRO", "CB", 'C', 0.2),
        testRecord(2, "SER", "OG", 'B', 0.6),
        testRecord(2, "PRO", "CG", 'A', 0.2),
    };
    try expectKept(.{ .mode = .auto }, &records, &.{ false, true, false, true, false, true });
    try expectKept(.{ .mode = .highest_occupancy }, &records, &.{ true, false, true, false, true, false });
}

test "resolve chooses among the alternates of the surviving residue" {
    // SER has two conformers of its own next to PRO
    const records = [_]TestRecord{
        testRecord(2, "PRO", "CA", 'A', 0.5),
        testRecord(2, "SER", "CA", 'B', 0.2),
        testRecord(2, "SER", "CA", 'C', 0.3),
        testRecord(2, "SER", "OG", 'B', 0.2),
        testRecord(2, "SER", "OG", 'C', 0.3),
    };
    try expectKept(.{ .mode = .auto }, &records, &.{ true, false, false, false, false });
    try expectKept(.{ .mode = .highest_occupancy }, &records, &.{ true, false, false, false, false });
    try expectKept(.{ .mode = .selected, .id = 'C' }, &records, &.{ false, false, true, false, true });

    // Without PRO, each atom site of SER keeps its best alternate
    try expectKept(.{ .mode = .auto }, records[1..], &.{ false, true, false, true });
}

test "resolve leaves residues without alternates alone" {
    // Two residues without altLocs share chain and number (a file without
    // chain IDs), next to a residue with alternates at the same position
    const records = [_]TestRecord{
        testRecord(1, "ALA", "CA", ' ', 1.0),
        testRecord(1, "GLY", "CA", ' ', 1.0),
        testRecord(1, "SER", "CA", 'A', 0.5),
        testRecord(1, "PRO", "CA", 'B', 0.5),
    };
    for ([_]AltLocMode{ .auto, .highest_occupancy }) |mode| {
        try expectKept(.{ .mode = mode }, &records, &.{ true, true, true, false });
    }
}

test "resolve tells residue positions apart" {
    const base = testRecord(1, "ALA", "CA", 'B', 0.5).site;
    var other_model = base;
    other_model.model_num = 2;
    var other_chain = base;
    other_chain.chain_id = "B";
    var other_seq = base;
    other_seq.seq = 2;
    var other_numbering = base;
    other_numbering.seq_is_auth = true;
    var other_insertion = base;
    other_insertion.insertion_code = "A";
    var same = base;
    same.alt_loc = 'C';

    const records = [_]TestRecord{
        .{ .site = base },
        .{ .site = other_model },
        .{ .site = other_chain },
        .{ .site = other_seq },
        .{ .site = other_numbering },
        .{ .site = other_insertion },
        .{ .site = same },
    };
    try expectKept(.{ .mode = .auto }, &records, &.{ true, true, true, true, true, true, false });
}

test "resolve matches the reference rules on random records" {
    const allocator = std.testing.allocator;
    var prng = std.Random.DefaultPrng.init(0x425);
    const random = prng.random();

    // Few distinct values, so that sites collide in every way
    const chains = [_][]const u8{ "", "A", "B" };
    const insertion_codes = [_][]const u8{ "", "A" };
    const residues = [_][]const u8{ "ALA", "SER", "PRO" };
    const atom_names = [_][]const u8{ "N", "CA", "CB" };
    const alt_locs = [_]u8{ ' ', ' ', 'A', 'B', 'C' };
    const occupancies = [_]f64{ 0.0, 0.3, 0.5, 0.7 };
    const settings = [_]AltLocSetting{
        .{ .mode = .auto },
        .{ .mode = .highest_occupancy },
        .{ .mode = .selected, .id = 'B' },
    };

    var records: [48]TestRecord = undefined;
    for (0..400) |_| {
        const n = random.intRangeAtMost(usize, 1, records.len);
        for (records[0..n]) |*record| {
            record.* = .{ .site = .{
                .model_num = if (random.boolean()) null else 1,
                .chain_id = chains[random.uintLessThan(usize, chains.len)],
                .seq = random.intRangeAtMost(i32, 1, 3),
                .seq_is_auth = random.uintLessThan(u8, 8) == 0,
                .insertion_code = insertion_codes[random.uintLessThan(usize, insertion_codes.len)],
                .residue = residues[random.uintLessThan(usize, residues.len)],
                .atom_name = atom_names[random.uintLessThan(usize, atom_names.len)],
                .alt_loc = alt_locs[random.uintLessThan(usize, alt_locs.len)],
                .occupancy = occupancies[random.uintLessThan(usize, occupancies.len)],
            } };
        }

        for (settings) |setting| {
            const keep = try resolve(TestRecord, allocator, records[0..n], setting);
            defer if (keep) |flags| allocator.free(flags);
            for (0..n) |i| {
                const kept = if (keep) |flags| flags[i] else true;
                try std.testing.expectEqual(referenceKeeps(records[0..n], setting, i), kept);
            }
        }
    }
}

/// Counts how often `resolve` reads a record. Reading is the only way to
/// compare two records, so the count bounds the work whatever the machine.
const CountingRecord = struct {
    site: Site,

    var reads: usize = 0;

    pub fn altLocSite(self: CountingRecord) Site {
        reads += 1;
        return self.site;
    }
};

test "resolve reads every record a constant number of times" {
    const allocator = std.testing.allocator;
    const n = 60_000;

    const records = try allocator.alloc(CountingRecord, n);
    defer allocator.free(records);
    const names = try allocator.alloc([6]u8, n);
    defer allocator.free(names);
    for (names, 0..) |*name, i| _ = std.fmt.bufPrint(name, "{d:0>6}", .{i}) catch unreachable;

    const Shape = enum {
        /// Every atom is a residue of its own with altLoc A
        one_atom_per_residue,
        /// One residue with n/2 atom sites of two alternates each
        one_residue,
        /// One atom site with n alternates
        one_site,
        /// 100 residues whose atoms are interleaved, with four alternates each
        interleaved,
        /// One residue position with n/2 residue names of two atoms each
        one_position,
    };

    for (std.enums.values(Shape)) |shape| {
        for (records, 0..) |*record, i| {
            record.* = .{ .site = .{
                .model_num = 1,
                .chain_id = "A",
                .seq = 1,
                .insertion_code = "",
                .residue = "ALA",
                .atom_name = "CA",
                .alt_loc = 'A',
                .occupancy = 0.5,
            } };
            switch (shape) {
                .one_atom_per_residue => record.site.seq = @intCast(i),
                .one_residue => {
                    record.site.atom_name = &names[i / 2];
                    record.site.alt_loc = if (i % 2 == 0) 'B' else 'C';
                },
                .one_site => record.site.alt_loc = if (i % 2 == 0) 'B' else 'C',
                .interleaved => {
                    record.site.seq = @intCast(i % 100);
                    record.site.atom_name = &names[i / 400];
                    record.site.alt_loc = "BCDE"[(i / 100) % 4];
                },
                .one_position => {
                    record.site.residue = &names[i / 2];
                    record.site.atom_name = if (i % 2 == 0) "CA" else "CB";
                    record.site.alt_loc = 'B';
                },
            }
        }

        for ([_]AltLocMode{ .auto, .highest_occupancy }) |mode| {
            CountingRecord.reads = 0;
            const keep = (try resolve(CountingRecord, allocator, records, .{ .mode = mode })).?;
            defer allocator.free(keep);

            // The quadratic scan this replaces read every record once for
            // each alternate
            try std.testing.expect(CountingRecord.reads <= 4 * n);

            const n_kept = std.mem.count(bool, keep, &.{true});
            try std.testing.expectEqual(@as(usize, switch (shape) {
                .one_atom_per_residue => n,
                .one_residue => n / 2,
                .one_site => 1,
                .interleaved => n / 4,
                .one_position => 2,
            }), n_kept);
        }
    }
}

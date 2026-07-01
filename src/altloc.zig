const std = @import("std");

pub const AltLocMode = enum {
    /// Historical behavior: blank altLoc wins, then A, then highest occupancy.
    auto,
    /// Assume there are no non-blank altLocs; fail fast if one is encountered.
    none,
    /// Keep all alternate locations.
    all,
    /// Keep blank altLocs and the selected altLoc ID.
    selected,
    /// Keep the highest-occupancy alternate for each atom site.
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

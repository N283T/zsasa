//! zsasa: High-performance SASA (Solvent Accessible Surface Area) calculation library.
//!
//! ## Modules
//!
//! - `shrake_rupley` - Shrake-Rupley algorithm
//! - `lee_richards` - Lee-Richards algorithm
//! - `types` - Common types (AtomInput, Config, SasaResult)
//! - `pdb_parser` - PDB format parser
//! - `mmcif_parser` - mmCIF format parser
//! - `json_parser` - JSON format parser
//! - `classifier` - Atom classifier (FreeSASA-compatible)
//! - `analysis` - Result aggregation and RSA calculation
//!
//! ## Usage
//!
//! ```zig
//! const zsasa = @import("zsasa");
//! var parser = zsasa.pdb_parser.PdbParser.init(allocator);
//! const atoms = try parser.parse(pdb_data);
//! const result = try zsasa.shrake_rupley.calculateSasa(allocator, atoms, .{});
//! ```

// Public API
pub const shrake_rupley = @import("shrake_rupley.zig");
pub const lee_richards = @import("lee_richards.zig");
pub const types = @import("types.zig");
pub const pdb_parser = @import("pdb_parser.zig");
pub const mmcif_parser = @import("mmcif_parser.zig");
pub const json_parser = @import("json_parser.zig");
pub const classifier = @import("classifier.zig");
pub const analysis = @import("analysis.zig");

test {
    // Algorithm modules transitively pull in all internal modules' tests,
    // which conflicts with the exe_tests target. Only reference leaf modules.
    _ = types;
    _ = pdb_parser;
    _ = mmcif_parser;
    _ = json_parser;
    _ = classifier;
    _ = analysis;
}

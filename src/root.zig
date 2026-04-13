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
//! - `toml_parser` - Minimal TOML subset parser
//! - `toml_classifier_parser` - TOML-to-Classifier converter
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
pub const toml_parser = @import("toml_parser.zig");
pub const toml_classifier_parser = @import("toml_classifier_parser.zig");
pub const dcd = @import("dcd.zig");
pub const bitmask_lut = @import("bitmask_lut.zig");
pub const shrake_rupley_bitmask = @import("shrake_rupley_bitmask.zig");
pub const gzip = @import("gzip.zig");
pub const hybridization = @import("hybridization.zig");
pub const ccd_parser = @import("ccd_parser.zig");
pub const classifier_ccd = @import("classifier_ccd.zig");
pub const ccd_binary = @import("ccd_binary.zig");
pub const compile_dict = @import("compile_dict.zig");
pub const sdf_parser = @import("sdf_parser.zig");

test {
    _ = shrake_rupley;
    _ = lee_richards;
    _ = types;
    _ = pdb_parser;
    _ = mmcif_parser;
    _ = json_parser;
    _ = classifier;
    _ = analysis;
    _ = toml_parser;
    _ = toml_classifier_parser;
    _ = dcd;
    _ = bitmask_lut;
    _ = shrake_rupley_bitmask;
    _ = gzip;
    _ = hybridization;
    _ = ccd_parser;
    _ = classifier_ccd;
    _ = ccd_binary;
    _ = compile_dict;
    _ = sdf_parser;
}

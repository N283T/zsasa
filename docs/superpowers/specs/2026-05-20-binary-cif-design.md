# BinaryCIF Input Support Design

Date: 2026-05-20

## Summary

Add native BinaryCIF (`.bcif`) input support to `zsasa` for single-structure and batch SASA calculations. The first implementation focuses on decoding the `_atom_site` category into the existing `AtomInput` representation and deliberately avoids broader CIF object-model support.

## Goals

- Accept `.bcif`, `.bcif.gz`, and `.bcif.zst` inputs in `zsasa calc`.
- Include `.bcif`, `.bcif.gz`, and `.bcif.zst` in `zsasa batch` file discovery and processing.
- Preserve existing mmCIF filtering semantics for BinaryCIF input:
  - `--model`
  - `--chain`
  - `--auth-chain`
  - `--include-hydrogens`
  - `--include-hetatm`
- Decode enough BinaryCIF to read the `_atom_site` fields required for SASA calculation and atom classification.
- Keep the implementation native Zig so CLI binaries and Python wheels remain self-contained.

## Non-goals

- Do not add BinaryCIF writing.
- Do not expose a full generic BinaryCIF or CIF object API.
- Do not implement inline CCD extraction from BinaryCIF in the first cut.
- Do not explicitly add BinaryCIF topology support to `traj`; users rarely need this path, and it can be added later if needed.
- Do not change existing mmCIF, PDB, SDF, JSON, C ABI, or Python public behavior except where format discovery lists supported files.

## Background

BinaryCIF is not coordinates-only. It stores CIF data blocks, categories, columns, data encodings, and optional masks in a MessagePack container. A BinaryCIF file can therefore contain `_atom_site` and other categories such as `_chem_comp_atom` or `_chem_comp_bond`. For `zsasa`, the useful first slice is `_atom_site`, because the current mmCIF parser already maps that category into `AtomInput` and applies the relevant filters.

The existing mmCIF parser in `src/mmcif_parser.zig` extracts `_atom_site` text tokens and builds arrays for coordinates, radii, element, residue, atom name, chain ID, residue number, and insertion code. BinaryCIF support should mirror this behavior rather than introduce a parallel output schema.

## Architecture

### New module: `src/bcif_parser.zig`

Add a native BinaryCIF parser with a public API shaped like the existing mmCIF parser:

```zig
pub const BcifParser = struct {
    allocator: Allocator,
    atom_only: bool = true,
    skip_hydrogens: bool = true,
    first_alt_loc_only: bool = true,
    model_num: ?u32 = null,
    chain_filter: ?[]const []const u8 = null,
    use_auth_chain: bool = false,

    pub fn init(allocator: Allocator) BcifParser;
    pub fn parse(self: *BcifParser, source: []const u8) !AtomInput;
    pub fn parseFile(self: *BcifParser, io: std.Io, path: []const u8) !AtomInput;
};
```

`parseFile()` reads plain, gzip, and zstd-compressed files using the same compressed-file helpers used by existing parsers.

### Format detection

Extend `src/format_detect.zig`:

- Add `InputFormat.bcif`.
- Add supported extensions:
  - `.bcif`
  - `.bcif.gz`
  - `.bcif.zst`
  - `.BCIF` for plain uppercase consistency.
- Make `detectInputFormat()` return `.bcif` for these paths.

### CLI and batch integration

Update `src/calc.zig` and `src/batch.zig` to route `.bcif` to `BcifParser`. Pass through the same filter settings currently passed to `MmcifParser`.

`traj` is intentionally left unchanged in the first implementation. If format detection is shared in a way that affects topology parsing, ensure the new `.bcif` case reports an explicit unsupported topology-format error rather than silently falling through.

### Root module

Export the parser from `src/root.zig` as `bcif_parser` so Zig package consumers can use it directly, matching the existing parser modules.

## BinaryCIF decoding scope

The parser only needs enough generic decoding to reach and decode `_atom_site` columns. It should validate the top-level shape but does not need to materialize all categories.

Expected top-level structure:

- `version`
- `encoder`
- `dataBlocks`
  - `header`
  - `categories`
    - `name`
    - `rowCount`
    - `columns`
      - `name`
      - `data`
      - optional or nullable `mask`

The implementation should scan categories until it finds `_atom_site` or `atom_site`, then decode only columns that may be needed.

### Required `_atom_site` fields

Required for coordinates:

- `Cartn_x`
- `Cartn_y`
- `Cartn_z`

Optional but important for classification and filtering:

- `type_symbol`
- `label_atom_id`
- `auth_atom_id`
- `label_comp_id`
- `auth_comp_id`
- `label_asym_id`
- `auth_asym_id`
- `label_seq_id`
- `auth_seq_id`
- `pdbx_PDB_ins_code`
- `group_PDB`
- `label_alt_id`
- `pdbx_PDB_model_num`

Column selection should reuse the same preference rules as `MmcifParser`: label values are preferred over auth values except when `use_auth_chain` is enabled.

### Supported encodings

Implement the BinaryCIF encodings needed by RCSB/Mol* files and by fixtures:

- `ByteArray`
- `FixedPoint`
- `IntervalQuantization`
- `RunLength`
- `Delta`
- `IntegerPacking`
- `StringArray`

Decoding applies encodings in reverse order, as specified by BinaryCIF. Byte arrays are little-endian.

### Masks

Support column masks:

- no mask: all values are present
- mask value `0`: present
- mask value `1`: CIF `.`
- mask value `2`: CIF `?`

For nullable string values, map mask values `1` and `2` to the same semantics currently handled by text mmCIF null tokens. For required coordinates, masked values produce `InvalidCoordinate`.

## Error handling

Add BinaryCIF-specific parse errors where needed while mapping common failures to existing behavior where possible:

- missing or malformed top-level MessagePack object
- missing `dataBlocks`
- no `_atom_site` category
- missing required coordinate field
- unsupported BinaryCIF encoding kind
- invalid numeric type for a required field
- invalid coordinate value
- inconsistent decoded column length versus category `rowCount`

Errors should be explicit enough for CLI users to distinguish unsupported or corrupt BinaryCIF from a valid structure with no usable atoms.

## Data flow

1. `calc` or `batch` calls `format_detect.detectInputFormat(path)`.
2. `.bcif` routes to `BcifParser.parseFile()`.
3. `parseFile()` reads/decompresses bytes.
4. The MessagePack reader parses the BinaryCIF container.
5. The parser finds `_atom_site`.
6. Needed columns and masks are decoded.
7. Rows are filtered using the same logic as mmCIF.
8. Included rows are appended into `AtomInput` arrays.
9. The existing classifier and SASA calculation pipeline runs unchanged.

## Testing plan

### Unit tests

Add focused tests in `src/bcif_parser.zig` or helper modules for:

- MessagePack primitives used by BinaryCIF: maps, arrays, strings, binary blobs, integers, booleans, floats, nil.
- `ByteArray` decoding for signed/unsigned integers and floats.
- `IntegerPacking` signed and unsigned expansion.
- `Delta`, `RunLength`, `FixedPoint`, and `IntervalQuantization` decoding.
- `StringArray` decoding.
- mask handling for present, `.`, and `?` values.

### Parser tests

- Generate a tiny synthetic BCIF with `_atom_site` and assert parsed coordinates and metadata.
- Compare a generated BCIF against an equivalent text mmCIF for atom count, coordinates, elements, residue names, atom names, chains, residue numbers, and insertion codes.
- Include filter coverage for HETATM, hydrogens, chain filter, auth chain selection, alternate locations, and model selection.

### Integration tests

- `zsasa calc <fixture>.bcif /tmp/out.json`
- `zsasa calc <fixture>.bcif.gz /tmp/out.json`
- `zsasa calc <fixture>.bcif.zst /tmp/out.json`
- batch directory discovery includes `.bcif`, `.bcif.gz`, and `.bcif.zst`.

Fixtures may be generated using `py-mmcif` or `rcsb.utils.io`, but generated files committed to the repository should be small and deterministic.

## Documentation updates

Update user-facing docs that list supported input formats:

- `README.md`
- `website/docs/cli/input.md`
- `website/docs/python-api/core.md` if batch-supported formats are listed there
- `website/docs/changelog.md`

Documentation should state that BinaryCIF support currently reads `_atom_site` for SASA calculation and does not yet use inline CCD data from BinaryCIF.

## Compatibility and performance notes

- Existing output schemas remain unchanged.
- Existing mmCIF behavior remains unchanged.
- Native decoding avoids adding runtime Python or external binary dependencies.
- The first implementation can allocate decoded selected columns before row conversion. Later optimization can stream or selectively decode fewer values if large BCIF files expose memory pressure.
- Coordinate precision should preserve the decoded values supplied by the file. Fixed-point and interval-quantized columns decode to `f64` for `AtomInput`.

## Open follow-ups

- Add BinaryCIF topology support to `traj` if users ask for it.
- Add inline CCD extraction from BCIF once `_atom_site` support is stable.
- Consider reusing the MessagePack subset reader for other compact binary formats if future features need it.

---
sidebar_position: 2
---

# Input Formats

The input format is auto-detected from the file extension.

| Extension | Format |
|-----------|--------|
| `.json`, `.json.gz`, `.json.zst` | JSON |
| `.cif`, `.cif.gz`, `.cif.zst`, `.mmcif`, `.mmcif.gz`, `.mmcif.zst`, `.CIF`, `.mmCIF` | mmCIF |
| `.bcif`, `.bcif.gz`, `.bcif.zst`, `.BCIF` | BinaryCIF (`_atom_site` support) |
| `.pdb`, `.pdb.gz`, `.pdb.zst`, `.PDB`, `.ent`, `.ent.gz`, `.ent.zst`, `.ENT` | PDB |
| `.sdf`, `.sdf.gz`, `.sdf.zst`, `.mol`, `.mol.gz`, `.mol.zst` | SDF/MOL small molecules |

`.gz` inputs may hold several concatenated gzip members, as produced by `cat a.gz b.gz` or by `bgzip` (BGZF). All members are decompressed and joined, and the CRC32 and size of each member are verified. A `.gz` file with trailing bytes that are not a gzip member, including zero padding, is rejected as corrupt, and so is a file that ends inside a member, such as an incomplete download.

## JSON Format

Minimal JSON input with coordinates and radii:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [4.0, 5.0, 6.0],
  "z": [7.0, 8.0, 9.0],
  "r": [1.7, 1.55, 1.52]
}
```

Extended JSON with classification info (required for `--classifier`):

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [4.0, 5.0, 6.0],
  "z": [7.0, 8.0, 9.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["N", "CA", "C"],
  "element": [7, 6, 6]
}
```

### Fields

| Field | Required | Description |
|-------|----------|-------------|
| `x`, `y`, `z` | Yes | Atom coordinates in Å |
| `r` | Yes | Van der Waals radii in Å |
| `residue` | For classifier | 3-letter residue code (e.g., "ALA") |
| `atom_name` | For classifier | Atom name (e.g., "CA", "N") |
| `element` | Optional | Atomic numbers (e.g., 6=C, 7=N, 8=O) |

### Validation Rules

- All arrays must have the same length
- Arrays cannot be empty
- Coordinates must be finite (no NaN or Inf)
- Radii must be positive and ≤ 100 Å

## mmCIF Format

Standard mmCIF files are supported. The parser extracts:

- `_atom_site.Cartn_x/y/z` - Coordinates
- `_atom_site.type_symbol` - Element (for VdW radius)
- `_atom_site.label_atom_id` / `auth_atom_id` - Atom name
- `_atom_site.label_comp_id` / `auth_comp_id` - Residue name
- `_atom_site.label_asym_id` / `auth_asym_id` - Chain ID
- `_atom_site.label_seq_id` / `auth_seq_id` - Residue number
- `_atom_site.pdbx_PDB_ins_code` - Insertion code
- `_atom_site.pdbx_PDB_model_num` - Model number
- `_atom_site.label_alt_id` - Alternate location (`--altloc=auto` by default)

Residue numbers come from `label_seq_id`. Non-polymer residues (waters, ligands, glycans) have no `label_seq_id` and are numbered by `auth_seq_id`, the number they have in the PDB-format file. With `--auth-chain`, every residue is numbered by `auth_seq_id`, so chain IDs, residue numbers and insertion codes together match the PDB-format file; a row without a usable `auth_seq_id` falls back to `label_seq_id`. A residue with neither value is numbered 0. BinaryCIF input follows the same rules.

Alternate-location handling can be controlled with `--altloc=MODE` for mmCIF and BinaryCIF input. `auto` preserves the historical behavior (blank alternate location first, then `A`, then highest occupancy), `none` assumes no non-blank alternate locations and errors if one is found, `all` keeps every alternate, `highest-occupancy` keeps the highest-occupancy atom for each site, and a single ID such as `--altloc=A` keeps blank atoms plus that alternate ID.

## BinaryCIF Format

BinaryCIF input decodes `_atom_site` for SASA calculation, supports the same `--altloc` policy as mmCIF, and uses embedded `_chem_comp_atom` / `_chem_comp_bond` inline CCD data when `--classifier=ccd` (or the `ccd` default) needs bond topology for non-standard compounds. You can still provide external CCD or SDF topology when the BinaryCIF file does not include component topology.

## PDB Format

Standard PDB format files are supported with ATOM and HETATM records.

The element of each atom is read from columns 77-78. It decides which atoms are hydrogens (removed unless `--include-hydrogens` is given) and which generic radius an atom gets when the classifier has no entry for it. When those columns are blank or do not hold an element symbol (files written before the element column existed can carry an ID code and line number there), the element is inferred from the atom name in columns 13-16:

- An ion is a residue named after its atom: `CA` in residue `CA` is calcium, `HG` in residue `HG` is mercury.
- Otherwise the columns decide. A name that starts in column 14 is a one-letter element (` CA ` is an alpha carbon, ` NA ` in a heme is nitrogen, ` HG ` is a hydrogen), and a name that starts in column 13 begins with a two-letter element (`FE  `, `ZN  `, `CL1 `, `SE  `).
- Four-character names fill column 13 whatever their element, so `HG21` or `HD11` is a hydrogen.

Some programs left-justify or center atom names instead (`CA  ` for an alpha carbon). zsasa detects such files and then relies on the names alone: a name starting with H, C, N, O, P or S is that element. In those files a metal or halogen inside a larger residue (`CL1` in a ligand) cannot be told from carbon, so write the element column if you can.

## SDF/MOL Format

SDF and MOL files are supported for small-molecule SASA. V2000 and V3000 records are accepted, and batch mode expands multi-molecule SDF files so each molecule is calculated independently. Use `--mol=NAME_OR_INDEX` to select one molecule from a multi-molecule SDF.

The first line of a record is the molecule title. It may be blank, which is what RDKit writes for a molecule without a name; select such a molecule by its 1-based index. Blank lines after a `$$$$` separator and at the end of the file are skipped.

Hydrogens are excluded unless `--include-hydrogens` is given. Deuterium and tritium atoms (symbols `D` and `T`) are hydrogens: they are excluded and included together with `H` atoms.

In batch mode each molecule is named `stem_title` after the file stem and the molecule title, or `stem_N` when the title is blank; molecules of one file that share a title get their position appended (`stem_title_N`). See [SDF and MOL Output Names](../guide/batch.md#sdf-and-mol-output-names) for the full rules and the characters replaced in output file names.

## Trajectory Formats

For the `traj` subcommand, the following trajectory formats are supported:

| Extension | Format | Coordinates |
|-----------|--------|-------------|
| `.xtc` | GROMACS XTC | nm (auto-converted to Å) |
| `.trr` | GROMACS TRR | nm (auto-converted to Å) |
| `.dcd` | NAMD/CHARMM DCD | Å |
| `.nc`, `.ncdf` | AMBER NetCDF | Å |

A topology file (PDB or mmCIF) is required alongside the trajectory to provide atom names and radii classification.

## Classifiers

Built-in classifiers assign atom radii based on residue and atom names. See [Classifiers](../guide/classifiers.mdx) for detailed documentation.

### Built-in Classifiers

| Classifier | Description |
|------------|-------------|
| `naccess` | NACCESS-compatible radii (Hubbard & Thornton 1993) |
| `ccd` | CCD bond-topology radii — **default for calc/batch PDB/mmCIF** |
| `protor` | Static ProtOr-compatible radii without runtime CCD resource parsing |
| `oons` | OONS radii (Ooi et al. 1987) |

### Custom Config (`--config=FILE`)

Custom classifiers are an advanced feature for studies that require a specific non-standard radius table. The CCD classifier remains the recommended default for most structures.

Custom classifier files must use TOML and the `.toml` extension.

#### TOML Format

```toml
name = "my-classifier"

[types]
C_ALI = { radius = 1.87, class = "apolar" }
C_CAR = { radius = 1.76, class = "apolar" }
N     = { radius = 1.65, class = "polar" }
O     = { radius = 1.40, class = "polar" }
S     = { radius = 1.85, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"

[[atoms]]
residue = "ALA"
atom = "CB"
type = "C_ALI"
```

- `name` - Classifier name (optional, default: "custom")
- `[types]` - Define atom types with radius (angstrom) and class (`"polar"` or `"apolar"`)
- `[[atoms]]` - Map (residue, atom) pairs to defined types. Use `"ANY"` for fallback entries.

#### Migrating from legacy FreeSASA-style text

Older custom classifier text files are no longer loaded directly. Convert each type row and atom mapping to TOML:

```text
# Old FreeSASA-style text
C_ALI 1.87 apolar
ANY CA C_ALI
```

```toml
[types]
C_ALI = { radius = 1.87, class = "apolar" }

[[atoms]]
residue = "ANY"
atom = "CA"
type = "C_ALI"
```

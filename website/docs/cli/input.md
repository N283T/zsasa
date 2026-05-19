---
sidebar_position: 2
---

# Input Formats

The input format is auto-detected from the file extension.

| Extension | Format |
|-----------|--------|
| `.json`, `.json.gz`, `.json.zst` | JSON |
| `.cif`, `.cif.gz`, `.cif.zst`, `.mmcif`, `.mmcif.gz`, `.mmcif.zst`, `.CIF`, `.mmCIF` | mmCIF |
| `.pdb`, `.pdb.gz`, `.pdb.zst`, `.PDB`, `.ent`, `.ent.gz`, `.ent.zst`, `.ENT` | PDB |

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
- `_atom_site.label_alt_id` - Alternate location (first kept by default)

## PDB Format

Standard PDB format files are supported with ATOM and HETATM records.

## Trajectory Formats

For the `traj` subcommand, the following trajectory formats are supported:

| Extension | Format | Coordinates |
|-----------|--------|-------------|
| `.xtc` | GROMACS XTC | nm (auto-converted to Å) |
| `.dcd` | NAMD/CHARMM DCD | Å |

A topology file (PDB or mmCIF) is required alongside the trajectory to provide atom names and radii classification.

## Classifiers

Built-in classifiers assign atom radii based on residue and atom names. See [Classifiers](../guide/classifiers.mdx) for detailed documentation.

### Built-in Classifiers

| Classifier | Description |
|------------|-------------|
| `naccess` | NACCESS-compatible radii (Hubbard & Thornton 1993) |
| `ccd` | CCD bond-topology radii — **default for calc/batch PDB/mmCIF** |
| `protor` | Alias for CCD, accepted for compatibility |
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

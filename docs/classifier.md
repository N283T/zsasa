# Atom Classifier

A module that assigns van der Waals radii and polarity classes based on residue names and atom names.

## Overview

SASA calculation requires the van der Waals radius of each atom. The classifier follows FreeSASA's configuration file format and determines radii in this order:

1. **Residue-specific match**: Exact match of (residue name, atom name)
2. **ANY fallback**: Atom definitions common to all residues
3. **Element estimation**: Extract element symbol from atom name and estimate van der Waals radius

## Module Structure

| File | Description |
|------|-------------|
| `classifier.zig` | Core data structures, element-based estimation, ClassifierType |
| `classifier_naccess.zig` | NACCESS-compatible built-in classifier |
| `classifier_protor.zig` | ProtOr classifier (hybridization-based) |
| `classifier_oons.zig` | OONS classifier (former FreeSASA default) |
| `classifier_parser.zig` | FreeSASA-compatible configuration file parser |

## CLI Usage

### Built-in Classifiers

```bash
# NACCESS classifier (recommended default)
./zig-out/bin/freesasa_zig --classifier=naccess input.json output.json

# ProtOr classifier (hybridization-based)
./zig-out/bin/freesasa_zig --classifier=protor input.json output.json

# OONS classifier (former FreeSASA default)
./zig-out/bin/freesasa_zig --classifier=oons input.json output.json
```

### Custom Configuration File

```bash
# Use FreeSASA format configuration file
./zig-out/bin/freesasa_zig --config=my_radii.config input.json output.json
```

### Input Requirements

To use a classifier, the input JSON must contain `residue` and `atom_name` fields:

```json
{
  "x": [1.0, 2.0, 3.0],
  "y": [1.0, 2.0, 3.0],
  "z": [1.0, 2.0, 3.0],
  "r": [1.7, 1.55, 1.52],
  "residue": ["ALA", "ALA", "ALA"],
  "atom_name": ["CA", "CB", "C"]
}
```

The classifier overwrites the original `r` values and assigns radii based on residue/atom names.

### Priority

When both `--config` and `--classifier` are specified, `--config` takes precedence.

## API

### Classifier Selection

```zig
const classifier = @import("classifier.zig");

// Classifier types
const ClassifierType = classifier.ClassifierType;
// .naccess - NACCESS-compatible (default)
// .protor  - ProtOr (hybridization-based)
// .oons    - OONS (former FreeSASA default)

// Convert from string
const ct = ClassifierType.fromString("naccess");  // .naccess
const name = ClassifierType.naccess.name();       // "NACCESS"
```

### NACCESS Classifier

```zig
const naccess = @import("classifier_naccess.zig");

// Get radius (residue-specific → ANY → element estimation)
const radius: ?f64 = naccess.getRadius("ALA", "CA");  // 1.87

// Get polarity class (residue-specific → ANY → unknown)
const class: AtomClass = naccess.getClass("ALA", "O");  // .polar

// Get both
const props: ?AtomProperties = naccess.getProperties("CYS", "SG");
// props.radius = 1.85, props.class = .apolar
```

### ProtOr Classifier

Radii based on hybridization state (Tsai et al. 1999).

```zig
const protor = @import("classifier_protor.zig");

// Get radius (residue-specific → element estimation, no ANY fallback)
const radius: ?f64 = protor.getRadius("ALA", "CA");  // 1.88

// In ProtOr, sulfur is polar
const class: AtomClass = protor.getClass("CYS", "SG");  // .polar
```

### OONS Classifier

OONS radii (former FreeSASA default). Larger aliphatic carbon radii.

```zig
const oons = @import("classifier_oons.zig");

// Get radius (residue-specific → ANY → element estimation)
const radius: ?f64 = oons.getRadius("ALA", "CA");  // 2.00 (NACCESS is 1.87)

// In OONS, C_CAR and S are polar
const class: AtomClass = oons.getClass("ALA", "C");  // .polar
```

### Element-Based Estimation

```zig
const classifier = @import("classifier.zig");

// Estimate radius from atomic number (solves CA/Ca ambiguity)
const radius = classifier.guessRadiusFromAtomicNumber(6);   // 1.70 (Carbon)
const radius = classifier.guessRadiusFromAtomicNumber(20);  // 2.31 (Calcium)

// Estimate radius from element symbol
const radius = classifier.guessRadius("C");   // 1.70
const radius = classifier.guessRadius("FE");  // 1.26

// Extract element from atom name and estimate radius
const radius = classifier.guessRadiusFromAtomName(" CA "); // 1.70 (C)
const radius = classifier.guessRadiusFromAtomName("FE  "); // 1.26 (FE)

// Extract element symbol from atom name
const elem = classifier.extractElement(" CA "); // "C" (leading space → 1-char element)
const elem = classifier.extractElement("FE  "); // "FE" (no leading space → check 2-char)
```

### Explicit Element Identification via Atomic Number

Including an `element` field (atomic number array) in the input JSON resolves atom name ambiguity:

```json
{
  "x": [1.0, 2.0],
  "y": [3.0, 4.0],
  "z": [5.0, 6.0],
  "r": [1.7, 2.31],
  "atom_name": ["CA", "CA"],
  "element": [6, 20]
}
```

In the above example:
- 1st CA: Atomic number 6 = Carbon (Cα)
- 2nd CA: Atomic number 20 = Calcium (metal ion)

### Data Types

```zig
// Polarity class
pub const AtomClass = enum {
    polar,   // Hydrophilic
    apolar,  // Hydrophobic
    unknown, // Unknown
};

// Atom properties
pub const AtomProperties = struct {
    radius: f64,
    class: AtomClass,
};
```

## Classifier Comparison

| Property | NACCESS | ProtOr | OONS |
|----------|---------|--------|------|
| Aliphatic C | 1.87 Å | 1.88 Å | 2.00 Å |
| Aromatic C | 1.76 Å | 1.76 Å | 1.75 Å |
| Nitrogen | 1.65 Å | 1.64 Å | 1.55 Å |
| Oxygen | 1.40 Å | 1.42-1.46 Å | 1.40 Å |
| Sulfur | 1.85 Å (apolar) | 1.77 Å (polar) | 2.00 Å (polar) |
| ANY fallback | Yes | No | Yes |
| Classification | Atom type | Hybridization | Atom type |
| Reference | João Rodrigues | Tsai et al. 1999 | Ooi et al. |

## NACCESS Atom Types

Atom type definitions based on FreeSASA's naccess.config:

| Type | Radius (Å) | Class | Description | Examples |
|------|------------|-------|-------------|----------|
| C_ALI | 1.87 | apolar | Aliphatic carbon | CA, CB, CG (aliphatic) |
| C_CAR | 1.76 | apolar | Carbonyl/aromatic carbon | C (backbone), CG (aromatic) |
| C_NUC | 1.80 | apolar | Nucleic acid carbon | C1', C2', base carbons |
| N_AMD | 1.65 | polar | Amide nitrogen | N (backbone), NE, NH1/2 |
| N_AMN | 1.50 | polar | Amino nitrogen | NZ (LYS) |
| N_NUC | 1.60 | polar | Nucleic acid nitrogen | N1, N3, N7, N9 |
| O | 1.40 | polar | Oxygen | O, OG, OH |
| S | 1.85 | apolar | Sulfur | SG (CYS), SD (MET) |
| SE | 1.80 | apolar | Selenium | SE (SEC, MSE) |
| P | 1.90 | apolar | Phosphorus | P (nucleic acid backbone) |

## ProtOr Atom Types

Hybridization-based atom types based on Tsai et al. 1999:

| Type | Radius (Å) | Class | Description |
|------|------------|-------|-------------|
| C3H0 | 1.61 | apolar | sp2 carbon, no hydrogens |
| C3H1 | 1.76 | apolar | sp2 carbon, 1 hydrogen |
| C4H1 | 1.88 | apolar | sp3 carbon, 1 hydrogen |
| C4H2 | 1.88 | apolar | sp3 carbon, 2 hydrogens |
| C4H3 | 1.88 | apolar | sp3 carbon, 3 hydrogens (methyl) |
| N3H0 | 1.64 | polar | sp2 nitrogen, no hydrogens |
| N3H1 | 1.64 | polar | sp2 nitrogen, 1 hydrogen |
| N3H2 | 1.64 | polar | sp2 nitrogen, 2 hydrogens (amide) |
| N4H3 | 1.64 | polar | sp3 nitrogen, 3 hydrogens (amino) |
| O1H0 | 1.42 | polar | Carbonyl oxygen |
| O2H1 | 1.46 | polar | Hydroxyl oxygen |
| S2H0 | 1.77 | polar | Thioether sulfur |
| S2H1 | 1.77 | polar | Thiol sulfur |
| SE2H0 | 1.90 | polar | Selenoether (MSE) |
| SE2H1 | 1.90 | polar | Selenol (SEC) |
| P4H0 | 1.80 | polar | Phosphorus |

## OONS Atom Types

Atom types based on Ooi, Oobatake, Nemethy, Scheraga:

| Type | Radius (Å) | Class | Description |
|------|------------|-------|-------------|
| C_ALI | 2.00 | apolar | Aliphatic carbon |
| C_ARO | 1.75 | apolar | Aromatic carbon |
| C_CAR | 1.55 | polar | Carbonyl carbon |
| N | 1.55 | polar | Nitrogen |
| O | 1.40 | polar | Oxygen |
| S | 2.00 | polar | Sulfur |
| P | 1.80 | polar | Phosphorus |
| SE | 1.90 | polar | Selenium |
| U_POL | 1.50 | polar | Unknown polar (for ASX, GLX) |

## Supported Residues

### Amino Acids (20 standard + 2 non-standard)

| Residue | Description | Notes |
|---------|-------------|-------|
| ALA | Alanine | |
| ARG | Arginine | NH1/NH2 are N_AMD |
| ASN | Asparagine | |
| ASP | Aspartic acid | |
| CYS | Cysteine | SG is S |
| GLN | Glutamine | |
| GLU | Glutamic acid | |
| GLY | Glycine | |
| HIS | Histidine | Ring is C_CAR/N_AMD |
| ILE | Isoleucine | |
| LEU | Leucine | |
| LYS | Lysine | NZ is N_AMN |
| MET | Methionine | SD is S |
| PHE | Phenylalanine | Ring is C_CAR |
| PRO | Proline | |
| SER | Serine | OG is O |
| THR | Threonine | OG1 is O |
| TRP | Tryptophan | NE1 is N_AMD |
| TYR | Tyrosine | OH is O |
| VAL | Valine | |
| SEC | Selenocysteine | SE is SE |
| MSE | Selenomethionine | SE is SE |

### Nucleic Acids

| Residue | Description |
|---------|-------------|
| A, DA | Adenine (RNA/DNA) |
| C, DC | Cytosine (RNA/DNA) |
| G, DG | Guanine (RNA/DNA) |
| I, DI | Inosine (RNA/DNA) |
| T, DT | Thymine (RNA/DNA) |
| U, DU | Uracil (RNA/DNA) |

## PDB Atom Name Convention

PDB format atom names are 4 characters with positional significance:

```
Column 13-16: Atom name
- Leading space + 1-char element: " CA ", " N  ", " O  "
- 2-char element (no leading space): "FE  ", "ZN  ", "CA  " (Calcium)
```

The `extractElement` function follows this convention:
- Leading space → 1-char element (" CA " → "C")
- No leading space → Check for 2-char element ("CA  " → "CA", "CD1 " → "C")

## Element Radius Table

Van der Waals radii based on Mantina et al. 2009:

| Element | Radius (Å) | Element | Radius (Å) |
|---------|------------|---------|------------|
| H | 1.10 | FE | 1.26 |
| C | 1.70 | ZN | 1.39 |
| N | 1.55 | CU | 1.40 |
| O | 1.52 | MG | 1.73 |
| P | 1.80 | CA | 2.31 |
| S | 1.80 | NA | 2.27 |
| SE | 1.90 | K | 2.75 |

## Implementation Details

### O(1) Lookup

Uses `std.StaticStringMap` to build a compile-time hashmap:

```zig
const residue_atoms = std.StaticStringMap(AtomType).initComptime(.{
    .{ "ALA :CB  ", atom_types.C_ALI },
    .{ "ARG :CG  ", atom_types.C_ALI },
    // ...
});
```

Key format: `"RES :ATOM"` (9 characters, space-padded)

### Memory Usage

- All data is embedded at compile time
- No runtime allocation required
- Approximately 150 entries (amino acids + nucleic acids)

## Configuration File Format

Supports FreeSASA-compatible configuration file format. Allows defining custom atom radii.

### Basic Format

```
# Comment
name: MyClassifier

types:
C_ALI 1.87 apolar    # Type name, radius, class (polar/apolar)
C_CAR 1.76 apolar
O     1.40 polar

atoms:
ANY C   C_CAR        # Residue, atom name, type reference
ANY O   O
ANY CA  C_ALI
ALA CB  C_ALI        # Residue-specific definition
```

### Sections

| Section | Description |
|---------|-------------|
| `name:` | Classifier name (optional) |
| `types:` | Atom type definitions (name, radius, class) |
| `atoms:` | (residue, atom name) → type mapping |

### Library Usage

```zig
const parser = @import("classifier_parser.zig");

// Load from file
var classifier = try parser.parseConfigFile(allocator, "my_radii.config");
defer classifier.deinit();

// Use
const radius = classifier.getRadius("ALA", "CA");
```

## Completed Features

- [x] ProtOr classifier
- [x] OONS classifier
- [x] Configuration file parser (FreeSASA-compatible)
- [x] CLI integration (`--classifier`/`--config` options)

## References

- [FreeSASA naccess.config](https://github.com/mittinatten/freesasa/blob/master/share/naccess.config)
- [FreeSASA protor.config](https://github.com/mittinatten/freesasa/blob/master/share/protor.config)
- [FreeSASA oons.config](https://github.com/mittinatten/freesasa/blob/master/share/oons.config)
- Tsai, J., Taylor, R., Chothia, C., & Gerstein, M. (1999). The packing density in proteins: standard radii and volumes. J. Mol. Biol. 290:253-266.
- Mantina et al. (2009) Consistent van der Waals Radii for the Whole Main Group

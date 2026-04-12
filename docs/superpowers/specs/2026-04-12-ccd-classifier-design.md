# CCD Classifier Design Spec

## Summary

Introduce a new classifier (`ClassifierType.ccd`) for freesasa-zig that derives van der Waals radii from CCD (Chemical Component Dictionary) bond topology. By analyzing `_chem_comp_atom` and `_chem_comp_bond` data, the classifier determines each atom's hybridization and implicit hydrogen count, then maps these to ProtOr-compatible radii. This enables accurate radius assignment for any chemical component, not just standard amino acids.

## Motivation

Current classifiers (NACCESS, ProtOr, OONS) use hardcoded `(residue, atom_name) -> radius` lookup tables limited to standard amino acids and nucleic acids. Non-standard residues (modified AAs, ligands, cofactors) fall back to element-based VdW radii, losing hybridization-aware accuracy.

PDB's official mmCIF files now include inline `_chem_comp_bond`/`_chem_comp_atom` categories, making bond topology available without external CCD files.

## Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Radius granularity | Level 3: element + hybridization + H count | ProtOr/ChimeraX-compatible; CCD provides the data to derive this |
| CCD data supply | Hybrid: hardcoded + inline CCD + external CCD | Fast path for common residues; full coverage with CCD |
| Lookup priority | Hardcoded > inline CCD > external CCD > element fallback | Optimizes for the common case (stdAA-only structures) |
| New struct vs existing | New `CcdClassifier` struct (separate from `Classifier`) | Avoids forcing CCD semantics into flat HashMap; future flexibility |
| Integration method | `calc.zig` branches on `ClassifierType` | Simple; avoids premature abstraction of a shared interface |
| AtomClass assignment | Element-based: N,O=polar; C,S,P,Se=apolar | Sufficient for SASA; IDATM-level classification is overkill |
| Ionic radii | Out of scope | Requires coordination number from `_struct_conn`; defer |

## Architecture

### Data Flow

```
mmCIF file
  |
  +--> _atom_site (coordinates, residue names, atom names)
  |
  +--> _chem_comp_atom / _chem_comp_bond (inline CCD, if present)
         |
         v
CcdClassifier.getRadius(residue, atom_name)
  |
  1. Hardcoded table lookup (StaticStringMap, O(1))
  |    hit? -> return radius
  |
  2. Runtime component lookup (HashMap)
  |    - Populated from inline CCD or external CCD
  |    - Derive: bond orders -> hybridization -> implicit H count -> ProtOr radius
  |    hit? -> return radius
  |
  3. Element fallback (existing classifier.guessRadiusFromAtomName)
       -> return VdW radius
```

### Data Structures

#### CCD Component Representation

```zig
const CompAtom = struct {
    atom_id: FixedString4,
    type_symbol: FixedString4,
    aromatic: bool,
    leaving: bool,
};

const CompBond = struct {
    atom_idx_1: u16,
    atom_idx_2: u16,
    order: BondOrder,
    aromatic: bool,
};

const BondOrder = enum { single, double, triple, aromatic, delocalized, unknown };

const Component = struct {
    comp_id: FixedString5,
    atoms: []const CompAtom,
    bonds: []const CompBond,
};
```

#### Hybridization Derivation

```zig
const Hybridization = enum { sp3, sp2, sp, unknown };

const AtomTypeInfo = struct {
    hybridization: Hybridization,
    heavy_bond_count: u8,
    implicit_h_count: u8,
};
```

#### CcdClassifier

```zig
const CcdClassifier = struct {
    hardcoded: StaticStringMap(AtomProperties),
    runtime_components: HashMap(FixedString5, RuntimeComponent),
    allocator: Allocator,

    pub fn getRadius(self: *const CcdClassifier, residue: []const u8, atom: []const u8) ?f64;
    pub fn getClass(self: *const CcdClassifier, residue: []const u8, atom: []const u8) AtomClass;
    pub fn getProperties(self: *const CcdClassifier, residue: []const u8, atom: []const u8) ?AtomProperties;
    pub fn addComponent(self: *CcdClassifier, component: Component) !void;
    pub fn deinit(self: *CcdClassifier) void;
};
```

### Hybridization Derivation Pipeline

**Step 1: Bond order -> Hybridization** (reference: zreduce `ccd_derive.zig`)

```
For each atom in component:
  - Collect all bonds involving this atom
  - If any triple bond -> sp
  - If any double or aromatic bond -> sp2
  - Otherwise -> sp3
```

**Step 2: Implicit hydrogen count**

```
implicit_h = typical_valence(element) - heavy_bond_count

typical_valence: C=4, N=3, O=2, S=2, P=5
```

When CCD includes explicit hydrogen atoms in `_chem_comp_atom`, count H bonds directly instead.

**Step 3: ProtOr radius table lookup**

| Element | Hybridization | H count | Radius |
|---|---|---|---|
| C | sp2 | 0 | 1.61 |
| C | sp2 | 1 | 1.76 |
| C | sp3 | 0 | 1.61 (same as sp2/H0; quaternary C not in original ProtOr, inferred) |
| C | sp3 | 1-3 | 1.88 |
| N | any | any | 1.64 |
| O | sp2 | 0 | 1.42 |
| O | sp3 | 0 | 1.46 |
| O | sp3 | 1 | 1.46 |
| S | any | any | 1.77 |
| P | any | any | 1.80 |
| Se | any | any | 1.90 |

**Step 4: AtomClass**

- N, O -> `.polar`
- C, S, P, Se, others -> `.apolar`

### Hardcoded Components

The following comp_ids are pre-computed at compile time:

**Standard amino acids (20):** ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL

**Nucleic acids:** A, C, G, T, U, DA, DC, DG, DT, I, DI

**Modified residues:** MSE, SEC, PYL, HYP, MLY, SEP, TPO, ASX, GLX, PSU

**Caps/other:** ACE, NH2, HOH

Format: `StaticStringMap` keyed by `AtomKey(comp_id, atom_id)` -> `AtomProperties(radius, class)`.

### CCD Parsing

#### Inline CCD (mmCIF)

Extend `mmcif_parser.zig` to detect and parse `_chem_comp_atom`/`_chem_comp_bond` loops using the existing `cif_tokenizer.zig`. Only parse components whose comp_id is not in the hardcoded table.

#### External CCD

New `ccd_parser.zig` implementing a streaming parser (reference: zmmol `component_dict.zig`). Reads `_chem_comp_atom`/`_chem_comp_bond` loops from a multi-block CIF file without building a full document tree. Supports gzip-compressed input.

Future: zreduce-style binary CCD format for faster loading.

### Integration Points

#### classifier.zig

Add `.ccd` to `ClassifierType` enum.

#### calc.zig

Branch on classifier type:

```zig
if (classifier_type == .ccd) {
    var ccd_clf = CcdClassifier.init(allocator);
    defer ccd_clf.deinit();

    // Add inline CCD components (from mmCIF parse result)
    for (inline_components) |comp| {
        if (!ccd_clf.isHardcoded(comp.comp_id))
            try ccd_clf.addComponent(comp);
    }

    // Add external CCD components (if --ccd option)
    if (ccd_path) |path| {
        try ccd_clf.loadExternalCcd(path);
    }

    // Assign radii
    for (atoms, 0..) |atom, i| {
        radii[i] = ccd_clf.getRadius(atom.residue, atom.name)
            orelse classifier.guessRadiusFromAtomName(atom.name)
            orelse radii[i]; // keep element-based VdW
    }
} else {
    // existing classifier path
}
```

#### c_api.zig

Add `ZSASA_CLASSIFIER_CCD = 3`.

#### root.zig

Export `ccd_classifier`, `ccd_parser`, `hybridization` modules.

## File Plan

### New Files

| File | Purpose |
|---|---|
| `src/classifier_ccd.zig` | CcdClassifier struct, hardcoded table, public API |
| `src/ccd_parser.zig` | Streaming CCD `_chem_comp_atom`/`_chem_comp_bond` parser |
| `src/hybridization.zig` | Bond order -> hybridization -> H count -> radius derivation |

### Modified Files

| File | Change |
|---|---|
| `src/classifier.zig` | Add `.ccd` to `ClassifierType` |
| `src/mmcif_parser.zig` | Parse inline `_chem_comp_atom`/`_chem_comp_bond` |
| `src/calc.zig` | Add `.ccd` branch |
| `src/c_api.zig` | Add `ZSASA_CLASSIFIER_CCD` constant |
| `src/root.zig` | Export new modules |

## Test Strategy

1. **Hybridization unit tests** - Verify derivation from known CCD data (ALA, HEM, ATP) produces correct hybridization and H counts
2. **ProtOr compatibility tests** - For all stdAA atoms, CcdClassifier and ProtOr return identical radii (or document intentional differences)
3. **Inline CCD parse tests** - Parse `_chem_comp_bond` from real mmCIF files and verify correct Component construction
4. **E2E tests** - Compare `zsasa calc -c ccd` output with existing classifiers on test structures (1ubq.cif, 1crn.pdb)

## Out of Scope

- Ionic radii based on coordination number (requires `_struct_conn` parsing)
- Shared interface between `Classifier` and `CcdClassifier`
- External CCD binary format support
- Python bindings for `CcdClassifier`
- PDB format support for CCD (PDB files lack inline CCD data)

## References

- Tsai, J., Taylor, R., Chothia, C., & Gerstein, M. (1999). The packing density in proteins: standard radii and volumes. J Mol Biol, 290(1), 253-266. (ProtOr radii)
- ChimeraX radii: https://www.rbvi.ucsf.edu/chimerax/docs/user/radii.html
- zreduce `ccd_derive.zig`: hybridization derivation from CCD bond orders
- zmmol `component_dict.zig`: streaming CCD parser
- ChimeraX `stdresidues.cif`: template components for standard residues

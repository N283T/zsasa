---
sidebar_position: 3
---

# Output & Analysis

## Output Formats

### JSON (default)

Pretty-printed JSON with 2-space indentation:

```json
{
  "total_area": 18923.28,
  "atom_areas": [
    32.47,
    0.25,
    15.82
  ]
}
```

Use JSON for machine-readable output and long-term pipelines. Legacy text
formats such as `freesasa` and `rsa` are best-effort compatibility outputs for
tools or reports that expect FreeSASA/NACCESS-style text.

### Compact JSON

Single-line JSON without whitespace:

```json
{"total_area":18923.28,"atom_areas":[32.47,0.25,15.82]}
```

### CSV

Basic CSV with atom index and area:

```csv
atom_index,area
0,32.470000
1,0.250000
2,15.820000
total,18923.280000
```

When the input of `calc` has residue information (PDB, mmCIF, BinaryCIF and SDF/MOL input), `calc` writes the rich CSV instead, with one row per atom:

```csv
chain,residue,resnum,insertion_code,atom_name,x,y,z,radius,area
A,MET,1,,N,27.340,24.430,2.614,1.650,49.097438
A,MET,1,,CA,26.266,25.413,2.842,1.870,16.124513
H,SER,10,A,CB,28.523,15.820,8.182,1.870,45.686121
,,,,,,,,,699.097218
```

| Column | Content |
|--------|---------|
| `chain` | Chain ID, in full (mmCIF chain IDs can be longer than four characters) |
| `residue` | Residue name |
| `resnum` | Residue number |
| `insertion_code` | Insertion code; empty for a residue without one. Residues `10`, `10A` and `10B` differ only in this column |
| `atom_name` | Atom name |
| `x`, `y`, `z` | Coordinates in Å |
| `radius` | Atom radius in Å |
| `area` | SASA of the atom in Å² |

The last row holds the total area and leaves every other column empty. Text fields follow RFC 4180: a chain ID, residue name, insertion code or atom name that contains a comma, a double quote or a line break is enclosed in double quotes, with every double quote in it doubled (the PDB chain ID `,` is written as `","`). All other fields are written unquoted. Use a CSV parser instead of splitting lines at commas, and read the columns by name: `insertion_code` is new after zsasa 0.9.1, and it moved `atom_name` and the columns after it one position to the right.

`batch --format=csv` always writes the basic `atom_index,area` CSV per input file, also for structure input.

### FreeSASA-Compatible Text (`calc --format=freesasa`)

The single-structure `calc` command can write a FreeSASA-style text summary:

```text
## zsasa FreeSASA-compatible output ##

PARAMETERS
algorithm    : Shrake & Rupley
classifier   : naccess
probe-radius : 1.40
Test-points  : 100
input        : structure.pdb

RESULTS (A^2)
Total   :   18923.28
```

This format is intended for interoperability with workflows that expect a
FreeSASA-like human-readable report.

### RSA Text (`calc --format=rsa`)

The single-structure `calc` command can also write a residue table in the
`.rsa` format of NACCESS, which FreeSASA writes too:

```text
REM  zsasa FreeSASA/NACCESS-compatible RSA
REM  Absolute and relative SASAs for structure.pdb
REM  Atomic radii and polar/non-polar classes: NACCESS
REM  Reference values for relative SASA: Tien et al. 2013
REM  Algorithm: Shrake & Rupley
REM  Probe-radius: 1.40
REM  Test-points: 100
REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
RES MET A   1   241.62 107.9 149.97   N/A  91.65   N/A 169.86   N/A  71.76   N/A
RES GLN A   2   232.62 103.4 156.88   N/A  75.74   N/A  99.54   N/A 133.07   N/A
RES SER H  10A  224.86 145.1  83.12   N/A 141.74   N/A 100.89   N/A 123.97   N/A
END  Absolute sums over single chains surface
CHAIN  1 A      474.2        306.9        167.4        269.4        204.8
CHAIN  2 H      224.9         83.1        141.7        100.9        124.0
END  Absolute sums over all chains
TOTAL           699.1        390.0        309.1        370.3        328.8
```

`rsa` requires residue metadata, so use PDB/mmCIF input or another input format
that provides chain, residue name, residue number, and insertion code fields.
There is one `RES` row per residue as defined under
[Residue Identity](#residue-identity), and one `CHAIN` row per chain ID.

The all-atom relative value (`REL`, in percent) is the absolute value divided
by the maximum SASA of the residue type from Tien et al. (2013), the table
under [RSA Calculation](#rsa-calculation---rsa), whatever the classifier. It is
reported for the 20 standard amino acids. zsasa has no reference values for
the side-chain, main-chain, non-polar and polar columns, so their relative
values, and the all-atom relative value of any other residue, are printed as
`N/A`, FreeSASA's notation for a missing reference value.

The `Non-polar` and `All polar` columns split the area by atom, using the
polarity class that the active classifier gives each atom, the same classifier
that sets the radii (see
[Classifiers and CCD](../guide/classifiers.mdx#radii-reference-table)). The
classifiers disagree on some atoms: NACCESS classes sulfur as non-polar where
CCD, ProtOr and OONS class it as polar, and OONS classes carbonyl carbon as
polar. An atom that the classifier does not class (hydrogens, ligands outside
its tables) is classed by its element: N, O, P and S are polar, everything else
is non-polar. The `--polar` option prints the same split for the whole
structure.

#### RSA Column Layout

Rows follow the fixed columns of NACCESS, so readers that take fields by
position can read them. Columns are counted from 1:

| Row | Columns | Content |
|-----|---------|---------|
| `RES` | 1-3 | `RES` |
| | 5-7 | Residue name, right-justified |
| | 9 | Chain ID |
| | 10-13 | Residue number, right-justified |
| | 14 | Insertion code |
| | 16-80 | Five pairs of an absolute value (7 columns, 2 decimals) and a relative value (6 columns, 1 decimal): all atoms, side chain, main chain, non-polar, polar |
| `CHAIN` | 1-5 | `CHAIN` |
| | 6-8 | Number of the chain, right-justified |
| | 10 | Chain ID |
| | 12-21, 25-34, 38-47, 51-60, 64-73 | Absolute sums (10 columns, 1 decimal) in the order of the `RES` row |
| `TOTAL` | 1-5 | `TOTAL`, then the sums in the columns of the `CHAIN` row |

Chain ID, residue number and insertion code follow each other without a blank,
as in a PDB file (`A1000B`), so read `RES` rows by column, not by splitting at
blanks.

A row keeps these columns when the residue name has at most three characters,
the chain ID at most one, the residue number at most four (`-999` to `9999`)
and the insertion code at most one, absolute values are at most `999.99` and
relative values at most `999.9`. zsasa leaves the first column of every value
field blank, which keeps neighboring values apart. A label or value that is
too wide (a five-character residue name, a chain ID such as `AA`, a ligand with
more than 1000 Å²) is still written in full, with a blank between it and its
neighbors, so the rest of that row moves to the right. zsasa then prints a
warning that columns may be misaligned. Use `--format=json` for robust
machine-readable output.

:::note
Biopython's `Bio.PDB.NACCESS.process_rsa_data` reads these columns but converts
every relative value to a number, so it stops with `ValueError` at the first
`N/A`. This also happens with the RSA files of FreeSASA. Replace `   N/A` with
` -99.9` (the same width) before passing the lines to it.
:::

The `freesasa` and `rsa` formats are available for single `calc` runs only.
Batch output remains `json`, `compact`, `csv`, or `jsonl`.

### Trajectory Output (CSV)

The `traj` subcommand outputs CSV with per-frame total SASA:

```csv
frame,step,time,total_sasa
0,1,1.000,1840.88
1,2,2.000,1944.47
2,3,3.000,1848.46
...
```

## Analysis Features

### Per-Residue Aggregation (`--per-residue`)

Sums atom SASA per residue:

```
Per-residue SASA:
Chain  Res    Num       SASA  Atoms
----- ---- ------ ---------- ------
    A  MET      1     198.52     19
    A  LYS      2     142.31     22
    A  ALA      3      45.67      5
```

For mmCIF and BinaryCIF input, residue numbers are `label_seq_id` values by default and `auth_seq_id` values with `--auth-chain`; see [mmCIF Format](input.md#mmcif-format).

#### Residue Identity

The per-residue table, the [RSA text format](#rsa-text-calc---formatrsa) and the [JSONL residue map](../guide/batch.md#residue-maps-in-jsonl) use one definition of a residue, so they report the same residues with the same atoms and areas. Two atoms belong to the same residue when they are adjacent in the input and have the same

- chain ID (the full ID, also when it is longer than four characters),
- residue number,
- insertion code, and
- residue name.

A residue is therefore a run of consecutive atoms, and residues are listed in input order. Residues `10`, `10A` and `10B` are three residues, and so are `GLY 10` and `LYS 10` of one chain.

- **Non-contiguous residues.** If the atoms of a residue are not contiguous in the input, each run is reported as its own entry with the same labels. The entries are not merged, because the residue map describes a residue as an atom range (`residue_atom_start`, `residue_atom_count`).
- **Multi-model files.** `calc` and `batch` read all models superimposed (`calc --model=N` selects one), so every residue appears once per model. Each copy is its own entry, with the area that this copy has inside the superimposed structure; the copies are not summed. Use `calc --model=N` to get one entry per residue. (When each model consists of a single residue, consecutive models continue the same run and are reported as one entry.)

### RSA Calculation (`--rsa`)

Calculates Relative Solvent Accessibility (RSA = SASA / MaxSASA).

MaxSASA reference values from Tien et al. (2013):

| Residue | MaxSASA (Å²) | Residue | MaxSASA (Å²) |
|---------|-------------|---------|-------------|
| ALA | 129.0 | LEU | 201.0 |
| ARG | 274.0 | LYS | 236.0 |
| ASN | 195.0 | MET | 224.0 |
| ASP | 193.0 | PHE | 240.0 |
| CYS | 167.0 | PRO | 159.0 |
| GLN | 225.0 | SER | 155.0 |
| GLU | 223.0 | THR | 172.0 |
| GLY | 104.0 | TRP | 285.0 |
| HIS | 224.0 | TYR | 263.0 |
| ILE | 197.0 | VAL | 174.0 |

Output with RSA values (can exceed 1.0 for exposed terminal residues):

```
Per-residue SASA with RSA:
Chain  Res    Num       SASA    RSA  Atoms
----- ---- ------ ---------- ------ ------
    A  MET      1     198.52   0.89     19
    A  LYS      2     142.31   0.60     22
    A  ALA      3      45.67   0.35      5
```

### Polar/Nonpolar Summary (`--polar`)

Prints two summaries of the SASA, one by residue type and one by atom class. Automatically enables `--per-residue`.

```
Polar/Nonpolar SASA:
  Polar:       3478.56 Å² ( 72.1%) - 42 residues
  Nonpolar:    1344.73 Å² ( 27.9%) - 34 residues

Polar/Nonpolar SASA by atom class (classifier: NACCESS):
  Polar:       2353.88 Å² ( 48.8%) - 223 atoms
  Nonpolar:    2469.41 Å² ( 51.2%) - 379 atoms
```

**By residue type.** The first block adds up whole residues:

- **Polar**: ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR, TYR
- **Nonpolar**: ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, TRP, VAL
- **Unknown**: Non-standard residues (ligands, modified residues, etc.)

An `Unknown` line (`Unknown:  <area> Å² - <n> residues (excluded from %)`) is added only when the structure has non-standard residues with SASA; the percentages leave them out.

**By atom class.** The second block adds up atoms by the polarity class that the active classifier gives each of them, so it depends on `--classifier`. It is the split of the `Non-polar` and `All polar` columns of the [RSA text format](#rsa-text-calc---formatrsa): its two areas are the last two values of the `TOTAL` row. Atoms that the classifier does not class are classed by element (N, O, P and S polar, everything else nonpolar); when there are any, a line `(<n> atoms without a class from the classifier are classed by element)` follows. The two blocks answer different questions and do not agree: the carbon atoms of a polar residue count as polar in the first block and as nonpolar in the second.

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
H,GLY,10,,N,1.000,3.000,5.000,1.640,32.470000
H,GLY,10,,CA,2.000,4.000,6.000,1.880,0.250000
H,SER,10,A,N,3.000,5.000,7.000,1.640,15.820000
,,,,,,,,,18923.280000
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

The last row holds the total area and leaves every other column empty. Text fields follow RFC 4180: a chain ID, residue name, insertion code or atom name that contains a comma, a double quote or a line break is enclosed in double quotes, with every double quote in it doubled (the PDB chain ID `,` is written as `","`). All other fields are written unquoted, so use a CSV parser instead of splitting lines at commas. Read the columns by name: the `insertion_code` column was added after `resnum` in the release after 0.9.1, which moved `atom_name` and the columns after it one position to the right.

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

The single-structure `calc` command can also write a FreeSASA/NACCESS-style
RSA table:

```text
REM  zsasa FreeSASA/NACCESS-compatible RSA
REM  Absolute and relative SASAs for structure.pdb
REM  Atomic radii and reference values for relative SASA: naccess
REM  Algorithm: Shrake & Rupley
REM  Probe-radius: 1.40
REM  Test-points: 100
REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
RES ALA   A 1      30.00  23.3  20.00   N/A  10.00   N/A  20.00   N/A  10.00   N/A
END  Absolute sums over single chains surface
CHAIN  1   A       30.0         20.0         10.0         20.0         10.0
END  Absolute sums over all chains
TOTAL              30.0         20.0         10.0         20.0         10.0
```

`rsa` requires residue metadata, so use PDB/mmCIF input or another input format
that provides chain, residue name, residue number, and insertion code fields.
Relative all-atom RSA values are reported for standard amino acids; unavailable
relative values are printed as `N/A`, matching FreeSASA's convention.

The RSA text table follows legacy NACCESS-style fixed-width columns where
possible. If residue labels, residue numbers, chain IDs, or SASA/RSA values are
too wide for those columns, zsasa still writes the full values and prints a
warning that columns may be misaligned. Use `--format=json` for robust
machine-readable output.

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

Classifies residues and shows SASA breakdown. Automatically enables `--per-residue`.

- **Polar**: ARG, ASN, ASP, GLN, GLU, HIS, LYS, SER, THR, TYR
- **Nonpolar**: ALA, CYS, PHE, GLY, ILE, LEU, MET, PRO, TRP, VAL
- **Unknown**: Non-standard residues (ligands, modified residues, etc.)

```
Polar/Nonpolar SASA:
  Polar:       2345.67 Å² ( 45.2%) - 42 residues
  Nonpolar:    2845.23 Å² ( 54.8%) - 58 residues
```

An `Unknown` line (`Unknown:  <area> Å² - <n> residues (excluded from %)`) is added only when the structure has non-standard residues with SASA; the percentages leave them out.

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

When input has structural info (mmCIF/PDB), rich CSV is generated:

```csv
chain,residue,resnum,atom_name,x,y,z,radius,area
A,ALA,1,N,1.000,3.000,5.000,1.500,32.470000
A,ALA,1,CA,2.000,4.000,6.000,1.700,0.250000
,,,,,,,,18923.280000
```

### Trajectory Output (CSV)

The `traj` subcommand outputs CSV with per-frame total SASA:

```csv
frame,step,time,total_sasa
0,1,1.000,1866.44
1,2,2.000,1977.96
2,3,3.000,1884.93
...
```

## Analysis Features

### Per-Residue Aggregation (`--per-residue`)

Groups atom SASA by residue (chain + residue number + insertion code):

```
Per-residue SASA:
Chain  Res    Num       SASA  Atoms
----- ---- ------ ---------- ------
    A  MET      1     198.52     19
    A  LYS      2     142.31     22
    A  ALA      3      45.67      5
```

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
  Unknown:        0.00 Å² (  0.0%) -  0 residues
```

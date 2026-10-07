# Single-File Stress Benchmarks

Single-file benchmarks time one `calc` invocation per structure, from process start to written output. They cover eight structures from 10,919 to 4,506,416 atoms, in PDB and mmCIF form, at 100 sphere points.

## Runtime and memory

<div data-chart="single/atoms"></div>

At 10 threads `zsasa` f64 is the fastest tool on every structure in both formats, and it has the lowest peak memory. On the smallest structures the bitmask modes are slower than exact f64, because building the lookup table is a fixed cost that a single small structure does not repay.

<div data-chart="single/bars"></div>

### PDB input

<div data-table="single/t10-pdb"></div>

### mmCIF input

<div data-table="single/t10-mmcif"></div>

The "vs" columns are runtime ratios against `zsasa` f64: the comparator's runtime divided by the `zsasa` runtime. A dash means the comparator did not complete that input.

## Thread scaling

<div data-chart="single/threads"></div>

## Dataset

Inputs were normalized to protein-only files for fair comparator runs: hydrogens, alternate conformations, ligands, waters, and non-L-peptide chains were removed; atom and residue identifiers were wrapped into PDB field limits when needed. The mmCIF inputs are representations of the same eight workloads.

| Structure | Role | Atoms | Chains | Note |
| --- | --- | ---: | ---: | --- |
| AF-P49792-F10 | single-medium | 10,919 | 1 | AFDB single-chain case |
| AF-Q6ZS30-F1 | single-large | 21,611 | 1 | AFDB large single-chain case |
| AF-0000000066638622 | multi-medium | 14,618 | 2 | AFDB-derived two-chain case |
| AF-0000000065781219 | multi-large | 24,140 | 2 | AFDB-derived two-chain case |
| 3jc8 | 100k-atom | 107,500 | 3 | PDB assembly |
| 5vyc | RustSASA stress | 249,168 | 4 | Parser-heavy RustSASA case |
| 8rbs | FreeSASA stress | 164,605 | 5 | PDB coordinate-overflow parser case |
| 9fqr | maximum-size | 4,506,416 | 57 | Largest assembly in this subset |

Lahuta is excluded from this suite because the benchmarked SASA command targets AlphaFold-style chain-A inputs and cannot process this mixed multi-chain subset.

## Caveats

- `zsasa` and PDBTools.jl were measured in the 0.9.0 rerun. The FreeSASA and RustSASA values for PDB input are unchanged comparator measurements retained from the 0.6.0 suite.
- The 8rbs FreeSASA result for PDB input exposes fixed-width coordinate overflow: coordinates wider than the PDB 8.3 field are misread by the pinned FreeSASA parser, which produces a pathological setup time. The mmCIF form of 8rbs does not show it.
- The 5vyc RustSASA result for PDB input is dominated by parser time. The mmCIF form does not show it.
- FreeSASA did not complete the two AFDB-derived two-chain structures in mmCIF form, because it aborted without assigning atomic radii.
- PDBTools.jl timings include the full Julia wrapper invocation, including Julia startup.

## Evidence source

The charts and tables on this page are generated from `results/tables/single_file_t10_summary.csv` and `single_file_thread_scaling.csv` in [`N283T/zsasa-benchmarks`](https://github.com/N283T/zsasa-benchmarks).

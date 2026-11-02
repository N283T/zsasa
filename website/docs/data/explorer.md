---
sidebar_label: "SASA explorer"
layout: wide
description: "Search interface burial and polar/apolar SASA for 623,242 AFDB human dimer models and 550,122 SwissProt AFDB models, computed with zsasa."
---

# SASA explorer

Search the SASA, interface burial and polar/apolar split that `zsasa` computed for AFDB human heterodimers, homo-dimers and SwissProt models. Narrow the results with the sliders, click a column header to sort, and select a dimer row to see its full breakdown.

<nav class="jump" aria-label="On this page">
<a href="#search">Search</a>
<a href="#download">Download the tables</a>
<a href="#the-data">The data</a>
<a href="#what-the-columns-mean">Columns</a>
<a href="#how-the-values-were-computed">How the values were computed</a>
</nav>

<div id="search" data-explorer></div>

## Download {#download}

**Download CSV** above the table saves the rows that match the current search and filters (up to 100,000 rows). The full tables are Parquet files:

- [afdb_heterodimers.parquet](pathname:///zsasa/data/explorer/afdb_heterodimers.parquet): one row per heterodimer model (44 MB)
- [afdb_homodimers.parquet](pathname:///zsasa/data/explorer/afdb_homodimers.parquet): one row per homo-dimer model (9 MB)
- [swissprot.parquet](pathname:///zsasa/data/explorer/swissprot.parquet): one row per SwissProt model (17 MB)

The same files, with the runs and scripts behind them, are archived on Zenodo: [doi:10.5281/zenodo.23898582](https://doi.org/10.5281/zenodo.23898582).

### Licence and citation

The values are derived from AlphaFold Database models, which are available under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). When you reuse them, cite the data archive, `zsasa` and the sources of the models:

- Heterodimer and homo-dimer models and their scores (AFDB/NVIDIA complex release): Han, Y. *et al.* AlphaFold Database Expands to Proteome-Scale Quaternary Structures. *bioRxiv* (2026). [doi:10.64898/2026.03.27.714458](https://doi.org/10.64898/2026.03.27.714458)
- Data: Nagae, T. *Benchmark dataset and analysis artifacts for zsasa v0.12.1* (revision 3). Zenodo 2026. [doi:10.5281/zenodo.23898582](https://doi.org/10.5281/zenodo.23898582)
- `zsasa`: see [Citation](../citation.md).
- AlphaFold: Jumper, J. *et al.* Highly accurate protein structure prediction with AlphaFold. *Nature* **596**, 583–589 (2021). [doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2)
- AlphaFold Database: Varadi, M. *et al.* AlphaFold Protein Structure Database in 2024: providing structure coverage for over 214 million protein sequences. *Nucleic Acids Res.* **52**, D368–D375 (2024). [doi:10.1093/nar/gkad1011](https://doi.org/10.1093/nar/gkad1011); Bertoni, D. *et al.* AlphaFold Protein Structure Database 2025: a redesigned interface and updated structural coverage. *Nucleic Acids Res.* **54**, D358–D362 (2026). [doi:10.1093/nar/gkaf1226](https://doi.org/10.1093/nar/gkaf1226)

## The data

- **Heterodimers**: the 531,605 human–human heterodimer models of the AFDB/NVIDIA complex release.
- **Homo-dimers**: the 91,637 human homo-dimer models of the same release.
- **SwissProt**: the 550,122 SwissProt models of AFDB v6, with the SASA of each model and its polar/apolar split.

The model files come from:

| Set | Source |
| --- | --- |
| Heterodimers | AFDB/NVIDIA FTP release: [`heterodimers/`](https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/heterodimers/) (tar shards), scores in [`heterodimer_metadata.csv`](https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/) |
| Homo-dimers | AFDB/NVIDIA FTP release: the human (taxon 9606) models of [`homodimers/`](https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/homodimers/) (tar chunks of every species), scores in [`homodimer_metadata.csv`](https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/) |
| SwissProt | AFDB v6 bulk download: [`swissprot_pdb_v6.tar`](https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/) (27 GB); see the [AFDB download page](https://alphafold.ebi.ac.uk/download) |

The [README of the AFDB/NVIDIA release](https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/README.txt) describes its files and score columns. Most dimer models are only in this FTP release, not on the AFDB website.

The two dimer sets are the models of the large-scale workflow in the `zsasa` manuscript. Each dimer row has the SASA of chain A alone, chain B alone and the complex, the interface burial ΔSASA, the buried surface area (BSA), the polar/apolar split of each of these areas, and the AFDB/NVIDIA confidence scores.

The search runs in your browser on [DuckDB-WASM](https://duckdb.org/docs/api/wasm/overview). The first search of a set downloads its table (9 to 44 MB) into the browser; later searches run locally.

## What the columns mean

For a dimer with chains A and B:

| Column | Value |
| --- | --- |
| SASA A, SASA B | SASA of each chain on its own, cut out of the complex (Å²) |
| SASA AB | SASA of the complex (Å²) |
| ΔSASA | SASA A + SASA B − SASA AB: the area buried in the interface, counted on both sides (Å²) |
| BSA | ΔSASA / 2, the buried surface area per chain (Å²) |
| Polar, apolar | each area split by the polarity class of its atoms; ΔSASA polar + ΔSASA apolar = ΔSASA |
| Apolar % | apolar share of the interface (ΔSASA apolar / ΔSASA) |
| ipSAEmin | the smaller of the two directional ipSAE scores |
| ipTM, pDockQ | interface confidence scores from the AFDB/NVIDIA release |
| pLDDT | global pLDDT of the model, the mean over all residues |
| HQ | AFDB high-quality complex: ipSAEmin ≥ 0.6, pLDDT ≥ 70 and at most 10 backbone clashes |

The SwissProt table has the SASA of the whole model and its polar/apolar split, with the model's residue count and mean pLDDT.

:::warning[Predicted structures]
These are predicted models. A large ΔSASA is not evidence that two proteins interact: in the manuscript, burial and the confidence scores were only moderately correlated, and most heterodimers with a large interface are not high-quality models. Read ΔSASA together with ipSAEmin and the other scores.
:::

## How the values were computed

Every value comes from `zsasa` 0.13.0, the release that added the [polar/apolar split](../guide/batch.md#polar-split), with the settings of the manuscript's dimer workflow:

- Shrake–Rupley, bitmask mode, f32, 128 sphere points, 1.4 Å probe
- CCD classifier for radii and polarity classes ([Classifiers](../guide/classifiers.md)); atoms without a class are classed by element (N, O, P and S polar)
- dimers: a [BSA analysis](../guide/workflows.md#bsa-analysis) with chain A and chain B as the partners; SwissProt: a plain batch run

The bitmask mode at 128 points gives total SASA about 0.7% below the exact mode on static structures. The ΔSASA totals equal those of the `zsasa` 0.12.1 runs reported in the manuscript.

Chain accessions, genes, descriptions and pLDDT were read from the model files. ipTM, ipSAE, pDockQ and backbone clashes come from the score tables of the AFDB/NVIDIA release. The workflow files, scripts and tables are archived with the benchmark data; see [Citation](../citation.md).

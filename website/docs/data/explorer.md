---
sidebar_label: "SASA explorer"
description: "Search interface burial and polar/apolar SASA for 623,242 AFDB human dimer models and 550,122 SwissProt AFDB models, computed with zsasa."
---

# SASA explorer

Search the solvent accessible surface area that `zsasa` computed for two sets of AlphaFold Database (AFDB) models:

- **AFDB human dimers**: the 91,637 human homo-dimer models and 531,605 human–human heterodimer models of the AFDB/NVIDIA complex release. These are the models of the large-scale workflow in the `zsasa` manuscript. Each row has the SASA of chain A alone, chain B alone and the complex, the interface burial ΔSASA, the buried surface area (BSA), the polar/apolar split of each of these areas, and the AFDB/NVIDIA confidence scores.
- **SwissProt**: the 550,122 SwissProt models of AFDB v6, with the SASA of each model and its polar/apolar split.

The search runs in your browser on [DuckDB-WASM](https://duckdb.org/docs/api/wasm/overview), which reads only the parts of the data files that a query needs. Select a dimer row to see its full polar/apolar breakdown.

<div data-explorer></div>

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

## Data and licence

Download the full tables as Parquet:

- [afdb_dimers.parquet](pathname:///zsasa/data/explorer/afdb_dimers.parquet): one row per dimer model
- [swissprot.parquet](pathname:///zsasa/data/explorer/swissprot.parquet): one row per SwissProt model

The values are derived from AlphaFold Database models, which are available under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/). When you reuse them, cite the AlphaFold Database, the AFDB/NVIDIA complex release for the dimers, and `zsasa`.

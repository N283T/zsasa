# Paper Benchmark Website Refresh Design

## Goal

Refresh the `zsasa` website benchmark documentation so it matches the manuscript-oriented benchmark story from `zsasa-paper` and the pinned `v0.6.0` evidence in `zsasa-benchmarks`, while preserving only the older SwissProt benchmark as clearly labeled legacy context.

## Scope

Update benchmark-related website documentation and assets in the software repository:

- `website/docs/benchmarks/index.md`
- `website/docs/benchmarks/validation.md`
- `website/docs/benchmarks/batch.md`
- `website/docs/benchmarks/single-file.md`
- `website/docs/benchmarks/md.md`
- `website/docs/comparison.md`
- `website/src/pages/index.tsx`
- `website/static/assets/benchmarks/**`

Do not change Zig, Python, CLI, ABI, or benchmark runner behavior. Do not delete tracked non-website benchmark source data. Existing generated website benchmark images can be replaced when they are no longer referenced by the refreshed docs.

## Evidence Sources

Use these repositories as read-only sources:

- `/Users/nagaet/ghq/github.com/N283T/zsasa-paper`
  - Main benchmark narrative: `manuscript/manuscript.md`
  - Supplementary methodology and tables: `manuscript/supplement.md`
  - Paper figures: `manuscript/figures/benchmark/` and, if needed, `manuscript/figures/supplement/`
- `/Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks`
  - Summary CSVs: `results/tables/*.csv`
  - Website-ready exploratory figures: `results/figures/**/{png,svg}`
  - Tool and dataset metadata: `results/tables/tools.csv`, `results/tables/datasets.csv`

Numerical claims in the website should trace back to the CSV tables or manuscript text. Do not invent or manually adjust benchmark values.

## Information Architecture

Use a paper-first structure with five benchmark suites:

1. Static validation against FreeSASA.
2. Proteome-scale batch throughput.
3. Large and parser-heavy single-structure timing.
4. MD validation against MDTraj.
5. MD trajectory throughput.

Recommended page responsibilities:

- `benchmarks/index.md`: concise overview of the five suites, environment, pinned versions, and links to detail pages.
- `benchmarks/validation.md`: static FreeSASA agreement and MD MDTraj agreement, including exact versus bitmask interpretation.
- `benchmarks/batch.md`: E. coli and Human AFDB 10-thread batch results, thread-scaling summary for E. coli, and a legacy SwissProt section.
- `benchmarks/single-file.md`: the 8-structure curated subset, parser/large-structure stress story, runtime/RSS, and parse-versus-SASA timing.
- `benchmarks/md.md`: trajectory throughput and memory behavior for 5wvo_C, 6sup_A, and 5vz0_A.
- `comparison.md`: align performance summary and feature notes with current paper-era claims.
- `website/src/pages/index.tsx`: update homepage benchmark callout to the current v0.6.0 paper numbers.

## Figure and Table Design

Tables must not contain images. Every figure should be a standalone Markdown image block with nearby explanatory text.

Asset policy:

- Copy selected paper/benchmark images into a clean website asset directory such as `website/static/assets/benchmarks/paper/`.
- Prefer PNG for browser display. Keep SVG companions only when useful and when Docusaurus handles them cleanly.
- Replace old image references with independent figure blocks.
- Remove references to obsolete website images. Delete obsolete tracked website image files only after verifying they are unreferenced.

Suggested figure sources:

- Validation: `fig2_validation.png` from `zsasa-paper`, plus selected validation figures from `zsasa-benchmarks/results/figures/validation/png/` if a detail page needs more than the paper figure.
- Batch: `fig3_batch_ecoli.png`, `fig4_batch_human.png`, and optional supporting plots from `results/figures/batch_*`.
- Single-file: `fig5_single_runtime_rss.png`, `fig6_single_parse_sasa.png`, and optional supporting plots from `results/figures/single_file/png/`.
- MD: `fig7_md_throughput_rss.png` and optional supporting plots from `results/figures/md/png/`.

## Legacy SwissProt Handling

Keep SwissProt only as legacy reference context in `benchmarks/batch.md`:

- Label it clearly as a pre-pinned benchmark from the earlier website benchmark set.
- State that it was collected before the current `zsasa-benchmarks` pinned v0.6.0 harness and should not be mixed with the paper-era headline claims.
- Preserve the useful high-level observation: SwissProt-scale directory processing showed low memory and high throughput, but the current paper claims are based on E. coli and Human AFDB pinned reruns.
- Avoid presenting SwissProt as a current main result.

## Wording and Claims

Use careful language:

- State that exact f64/f32 modes reproduce FreeSASA under matched Shrake--Rupley settings.
- Describe bitmask as an approximate throughput-oriented mode with a quantified error envelope.
- Distinguish FreeSASA upstream from the `freesasa_batch` wrapper used for batch timing.
- Distinguish trajectory agreement from absolute ground truth.
- Avoid broad “up to” claims unless tied to a named workload, comparator, mode, and point count.

## Verification

After implementation:

1. Check for unreferenced or missing benchmark assets.
2. Run the Docusaurus build from `website/`. If dependencies are missing, run `npm ci` first.
3. Start the local website with `npm run start -- --host 127.0.0.1` or serve the production build, and show the user the local URL.
4. Report any checks that could not be completed.

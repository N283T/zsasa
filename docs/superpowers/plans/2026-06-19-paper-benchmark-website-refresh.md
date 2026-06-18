# Paper Benchmark Website Refresh Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refresh the website benchmark documentation and figures to match the `zsasa` manuscript and pinned `v0.6.0` benchmark evidence.

**Architecture:** Replace the older benchmark-story pages with paper-first Markdown pages backed by copied paper/benchmark figures and summary CSV values. Keep SwissProt only as a clearly labeled legacy pre-pinned result in the batch page. Verify the Docusaurus docs build and then serve locally for user review.

**Tech Stack:** Docusaurus Markdown docs, React/TypeScript homepage, static PNG/SVG assets, shell/Python helper scripts for asset and link checks.

---

### Task 1: Prepare Evidence and Assets

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/`
- Read-only source: `/Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/`
- Read-only source: `/Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/figures/`
- Read-only source: `/Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/tables/`

- [ ] **Step 1: Copy selected paper figures**

```bash
mkdir -p /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper
cp /Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/fig2_validation.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/fig3_batch_ecoli.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/fig4_batch_human.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/fig5_single_runtime_rss.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/fig6_single_parse_sasa.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-paper/manuscript/figures/benchmark/fig7_md_throughput_rss.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/
```

- [ ] **Step 2: Copy selected detailed benchmark figures**

```bash
mkdir -p /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/details
cp /Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/figures/validation/png/static_sr_mean_relative_error.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/details/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/figures/validation/png/md_r2.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/details/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/figures/batch_ecoli/png/ecoli_throughput_vs_threads.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/details/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/figures/single_file/png/single_runtime_vs_threads_grid.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/details/
cp /Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/figures/md/png/md_frames_per_sec_bar_grid.png /Users/nagaet/ghq/github.com/N283T/zsasa/website/static/assets/benchmarks/paper/details/
```

- [ ] **Step 3: Inspect source CSV values used by docs**

```bash
python3 - <<'PY'
import csv
from pathlib import Path
root = Path('/Users/nagaet/ghq/github.com/N283T/zsasa-benchmarks/results/tables')
for name in ['batch_t10_summary.csv', 'single_file_t10_summary.csv', 'md_summary.csv', 'validation_pairwise_summary.csv']:
    rows = list(csv.DictReader((root / name).open()))
    print(name, len(rows), rows[0])
PY
```

Expected: prints row counts and representative rows without errors.

### Task 2: Rewrite Benchmark Docs

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/benchmarks/index.md`
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/benchmarks/validation.md`
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/benchmarks/batch.md`
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/benchmarks/single-file.md`
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/benchmarks/md.md`

- [ ] **Step 1: Replace overview page with five-suite paper summary**

Write `index.md` with suite table, environment, evidence-source note, and links.

- [ ] **Step 2: Replace validation page**

Write `validation.md` with FreeSASA static validation, bitmask interpretation, MDTraj trajectory validation, figure blocks, and no images inside tables.

- [ ] **Step 3: Replace batch page**

Write `batch.md` with E. coli, Human, thread scaling, current paper-era summary, and a clearly labeled legacy SwissProt section.

- [ ] **Step 4: Replace single-file page**

Write `single-file.md` with the 8-structure subset, runtime/RSS, parse/SASA breakdown, parser caveats, and running notes.

- [ ] **Step 5: Replace MD page**

Write `md.md` with trajectory throughput, memory behavior, comparator interpretation, validation pointer, and running notes.

### Task 3: Align Related Website Copy

**Files:**
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs/comparison.md`
- Modify: `/Users/nagaet/ghq/github.com/N283T/zsasa/website/src/pages/index.tsx`

- [ ] **Step 1: Update comparison performance summary**

Replace stale 2,013-structure and 128-structure claims with current paper-era batch, single-file, and trajectory representative results.

- [ ] **Step 2: Update homepage benchmark callout**

Use current named workloads: E. coli 4,370 structures in 1.48 s with bitmask f32; Human 23,586 in 13.81 s; trajectory memory/speed claim tied to named comparators.

### Task 4: Verify Links, Assets, and Build

**Files:**
- Modify only if checks reveal broken references.

- [ ] **Step 1: Check image references**

```bash
python3 - <<'PY'
from pathlib import Path
import re
missing=[]
for p in Path('/Users/nagaet/ghq/github.com/N283T/zsasa/website/docs').rglob('*.md'):
    for img in re.findall(r'!\[[^\]]*\]\(([^)]+)\)', p.read_text()):
        if img.startswith('/zsasa/'):
            fs=Path('/Users/nagaet/ghq/github.com/N283T/zsasa/website/static') / img.removeprefix('/zsasa/')
        elif img.startswith('/assets/'):
            fs=Path('/Users/nagaet/ghq/github.com/N283T/zsasa/website/static') / img.removeprefix('/')
        elif img.startswith('pathname:///zsasa/'):
            fs=Path('/Users/nagaet/ghq/github.com/N283T/zsasa/website/static') / img.removeprefix('pathname:///zsasa/')
        else:
            continue
        if not fs.exists():
            missing.append((str(p), img, str(fs)))
print('\n'.join(map(str, missing)))
raise SystemExit(1 if missing else 0)
PY
```

Expected: no output and exit 0.

- [ ] **Step 2: Install website dependencies if needed**

```bash
cd /Users/nagaet/ghq/github.com/N283T/zsasa/website
if [ ! -d node_modules ]; then npm ci; fi
```

Expected: dependencies installed or already present.

- [ ] **Step 3: Build website**

```bash
cd /Users/nagaet/ghq/github.com/N283T/zsasa/website
npm run build
```

Expected: Docusaurus build succeeds.

- [ ] **Step 4: Serve locally for review**

```bash
cd /Users/nagaet/ghq/github.com/N283T/zsasa/website
npm run start -- --host 127.0.0.1
```

Expected: local URL printed for user review.

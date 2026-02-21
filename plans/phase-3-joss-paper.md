# Phase 3: JOSS Paper Preparation

## Goal

JOSS 投稿用の paper.md, paper.bib, figures を準備する。

## Prerequisites

- Phase 2 complete (v0.1.0 released, DOI obtained)
- ~6 months public history accumulated (~August 2026)

## Tasks

### 3.1 Paper Structure

- [ ] Create `paper/` directory
- [ ] `paper/paper.md` (750-1,750 words; RustSASA was ~1,100)
- [ ] `paper/paper.bib` (~20-30 references)
- [ ] `paper/figures/` directory

### 3.2 Paper Sections (paper.md)

1. **Summary** (~150 words) -- SASA importance, zsasa overview, key results
2. **Statement of Need** (~200 words) -- AlphaFold-scale datasets, existing tool limitations
3. **State of the Field** (~150 words) -- FreeSASA, RustSASA, MDTraj, Biopython
4. **Software Design** (~200 words) -- SoA layout, tiered SIMD, work-stealing thread pool, comptime generics, spatial hashing
5. **Research Impact** (~100 words) -- E. coli validation (R²=1.0), 2.3x speedup, 3.4x MD speedup
6. **AI Usage Disclosure** (~50-100 words) -- Tools used, scope, human design decisions
7. **Acknowledgements**

### 3.3 Key Differentiators (vs RustSASA -- "build vs contribute" justification)

| Feature | zsasa | RustSASA |
|---------|-------|----------|
| Algorithms | SR + LR | SR only |
| Precision | f64 default | f32 |
| Large structures | 2.3-3x faster | baseline |
| MD trajectories | 3.4x faster | baseline |
| Input formats | mmCIF, PDB, JSON, XTC | PDB |
| Thread control | explicit count | rayon global pool |
| Classifiers | NACCESS, ProtOr, OONS | ProtOr |
| Analysis | per-residue, RSA, polar/nonpolar | total SASA |

### 3.4 Figures (2-3 Composites)

- [ ] Figure 1: Validation
  - Panel A: SR validation scatter (`benchmarks/results/validation/ecoli/sr/validation_sr.png`)
  - Panel B: LR validation scatter (`benchmarks/results/validation/ecoli/lr/validation_lr.png`)
- [ ] Figure 2: Performance (3-4 panels)
  - Panel A: Large structure speedup (`benchmarks/results/plots/large/speedup_bar.png`)
  - Panel B: Thread scaling / efficiency (`benchmarks/results/plots/efficiency/summary.png`)
  - Panel C: Batch proteome (`benchmarks/results/plots/batch/ecoli_time_10t.png`)
  - Panel D: MD trajectory (`benchmarks/results/md/*/plots/bar.png`)
- [ ] 300+ DPI, panel labels (A, B, C, D)

### 3.5 Bibliography (paper.bib) -- Core References

**Algorithms:** Shrake & Rupley 1973, Lee & Richards 1971
**Tools:** FreeSASA (Mitternacht 2016), RustSASA (Campbell 2026), Biopython (Cock 2009), MDTraj (McGibbon 2015)
**MD frameworks:** MDAnalysis (Gowers 2016, Michaud-Agrawal 2011)
**Databases:** AlphaFold (Jumper 2021), AlphaFold DB (Varadi 2021)
**Classifiers:** ProtOr (Tsai 1999), OONS (Ooi 1987), NACCESS (Hubbard & Thornton 1993)
**Integrations:** Gemmi (Wojdyr 2022), Biotite (Kunzmann 2018)
**Other:** Zig language, Hyperfine (Peter 2023)

### 3.6 Paper CI

- [ ] `.github/workflows/paper.yml` -- `openjournals/paperdraft@master`
- [ ] `.github/workflows/ci.yml` paths-ignore に `paper/**` 追加

### 3.7 AI Usage Disclosure

JOSS 2026 new requirement. Must include:
- Tools/models used and versions
- Scope of assistance
- Human author reviewed/validated all outputs and made core design decisions
- Design decisions (SoA, tiered SIMD, thread pool, LR fast trig, comptime) are human-directed
- Validation (402 tests, E. coli proteome, 100k benchmarks) is human-verified

## Risks

| Risk | Impact | Mitigation |
|------|--------|------------|
| Short history (~29 days at creation) | HIGH | Submit after 6+ months, build community |
| AI-assisted dev questioned | HIGH | Transparent disclosure, design + validation evidence |
| Single author | MEDIUM | Community engagement, acknowledge contributors |
| FreeSASA author as reviewer | MEDIUM | Thorough honest validation |

## Reference: RustSASA Paper

- 61 lines, ~1,100 words, 25 references, 2 figures
- Sections: Summary, Statement of Need, Results, Methods, Acknowledgements
- DOI: 10.21105/joss.09537
- Local copy: `benchmarks/external/rustsasa-bench/paper/paper.md`

---
- [ ] **DONE** - Phase 3 complete

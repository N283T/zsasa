# Phase 1: MD Integration Documentation Update

## Goal
Update documentation to highlight MDTraj and MDAnalysis integrations as key features.

## Tasks

### 1. Main README.md Update
- [x] Add MD integrations to feature highlights
- [x] Add benchmark comparison (3.4x faster than mdsasa-bolt)
- [x] Add quick example for MDAnalysis usage

### 2. docs/python.md Update
- [x] Add MDTraj section with examples
- [x] Add MDAnalysis section with examples
- [x] Document `SASAAnalysis` class API
- [x] Document `compute_sasa()` function
- [x] Document precision parameter (`f32`/`f64`)
- [x] Add performance comparison table

### 3. New: docs/md-integration.md (Optional)
- [ ] Detailed MD trajectory analysis guide (skipped - docs/python.md is sufficient)
- [ ] Frame selection patterns
- [ ] Threading best practices
- [ ] Comparison with alternatives (mdakit-sasa, mdsasa-bolt)

## Key Points to Highlight

### Performance
- 3.4x faster than mdsasa-bolt (RustSASA) on real MD data
- Efficient frame-level parallelization
- Controllable thread count (unlike rayon's all-or-nothing)

### API Design
- Compatible with MDAnalysis `AnalysisBase` pattern
- MDTraj drop-in replacement for `mdtraj.shrake_rupley()`
- Consistent unit handling (Å² for MDAnalysis, nm² for MDTraj)

### Precision
- f64 precision by default (higher accuracy)
- f32 option for RustSASA/mdsasa-bolt comparison

## Benchmark Data to Include

| Implementation | Time (20k atoms, 1k frames) |
|----------------|------------------------------|
| freesasa-zig   | 8.8 s                        |
| mdsasa-bolt    | 30.3 s                       |
| Speedup        | **3.4x**                     |

## Current State

### Main README.md
- Has benchmark highlights for FreeSASA/RustSASA comparison
- **Missing**: MDTraj/MDAnalysis integrations
- **Missing**: MD trajectory benchmark (3.4x vs mdsasa-bolt)

### docs/python.md
- Has core API documentation
- Has gemmi/BioPython/Biotite integrations
- **Missing**: MDTraj section
- **Missing**: MDAnalysis section

### docs/README.md
- Index of all docs
- **Missing**: MD integration links

---
- [x] **DONE** - Phase complete

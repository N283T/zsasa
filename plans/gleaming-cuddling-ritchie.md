# Plan: Create `benchmarks/scripts/validation.py`

## Context

SASA validation (accuracy comparison between tools) is currently coupled with timing benchmarks in `bench.py` (single-file mode). Problems:

1. **Efficiency**: `bench.py` runs files one-by-one with timing overhead
2. **Interference**: Validation overhead affects timing measurements
3. **`bench_batch.py` can't help**: Uses hyperfine which discards all SASA output

New `validation.py` runs tools on a PDB directory, collects per-file SASA values, and compares across tools -- completely independent of timing benchmarks.

## Design

### CLI Interface (2 commands)

```bash
# Run tools and compare
./benchmarks/scripts/validation.py run \
    -i benchmarks/UP000000625_83333_ECOLI_v6/pdb \
    -n ecoli \
    --tool zig --tool freesasa --tool rustsasa \
    --algorithm sr --precision f64 --threads 1

# Re-analyze existing CSV
./benchmarks/scripts/validation.py compare \
    -d benchmarks/results/validation/ecoli
```

### Output

```
benchmarks/results/validation/<name>/
├── config.json          # System info, parameters
├── results.csv          # Per-file SASA values
└── validation_sr.png    # Scatter plot (always generated)
```

### CSV Format

```csv
structure,n_atoms,zsasa_f64,freesasa,rustsasa
1abc,1234,5678.90,5679.12,5677.45
```

Column name adapts to precision (`zsasa_f32`/`zsasa_f64`). RustSASA column omitted for LR.

### Data Flow

```
PDB directory
  │
  ├─ zsasa <input_dir> <tmp> --precision=f64 --algorithm=sr --threads=1
  │    → parse output JSONs: {"total_area": X}
  │
  ├─ freesasa <file.pdb> --shrake-rupley --resolution=100  (per file)
  │    → parse stdout: "Total   :  XXXXX"
  │
  └─ rust-sasa <input_dir> <tmp> --format json -t 1
  │    → parse output JSONs: {"total_area": X}
  │
  ▼
Merge by filename stem → CSV + scatter plot + statistics
```

## Tool Execution Strategy

### zsasa (batch mode)
- Command: `zsasa <input_dir> <output_dir> --algorithm=<a> --precision=<p> --threads=<n>`
- Output: one JSON per file in output_dir → `{"total_area": X, "atom_areas": [...]}`
- n_atoms: `len(atom_areas)` from JSON

### FreeSASA (per-file via CLI, no sasa_batch)
- Command per file: `freesasa <file.pdb> [--shrake-rupley --resolution=100 | --lee-richards --resolution=20]`
- Output: stdout default LOG format contains `Total   :  XXXXX`
- Parse: `re.search(r"[\d.]+", line)` on line starting with "Total" (same as `bench.py:302-308`)
- Supports both SR and LR algorithms (unlike `sasa_batch` which was SR-only)
- Binary: `benchmarks/external/freesasa-bench/src/freesasa`

### RustSASA (batch mode)
- Command: `rust-sasa <input_dir> <output_dir> --format json -t <n>`
- Output: same as zsasa (`{"total_area": X}`)
- SR only -- skip for LR with warning

## Implementation Steps

### Step 1: Script skeleton + CLI

Create `benchmarks/scripts/validation.py` with PEP 723 header.

**Reuse from `bench_batch.py`:**
- `Tool` enum, `get_root_dir()`, `get_binary_paths()`, `get_system_info()` (lines 57-111)
- `tempfile.TemporaryDirectory` pattern
- Binary existence check pattern

**Dependencies:** `polars`, `matplotlib`, `numpy`, `typer`, `rich`

**Binary paths:**
- `zsasa`: `root/zig-out/bin/zsasa`
- `freesasa`: `root/benchmarks/external/freesasa-bench/src/freesasa`
- `rustsasa`: `root/benchmarks/external/rustsasa-bench/target/release/rust-sasa`

### Step 2: Tool runners

**`run_zsasa(input_dir, tmp_dir, algorithm, precision, threads) → dict[str, tuple[float, int]]`**
- Run batch: `zsasa <input_dir> <tmp_dir> --algorithm=<a> --precision=<p> --threads=<n>`
- Scan `tmp_dir/*.json`, parse `total_area` and `len(atom_areas)`
- Returns `{stem: (total_sasa, n_atoms)}`

**`run_freesasa(input_dir, algorithm) → dict[str, float]`**
- Scan `input_dir/*.pdb`
- For each file: `freesasa <file> [--shrake-rupley --resolution=100 | --lee-richards --resolution=20]`
- Parse `Total` from stdout (regex from `bench.py:302-308`)
- Returns `{stem: total_sasa}`

**`run_rustsasa(input_dir, tmp_dir, threads) → dict[str, float]`**
- Run batch: `rust-sasa <input_dir> <tmp_dir> --format json -t <n>`
- Scan `tmp_dir/*.json`, parse `total_area`
- Returns `{stem: total_sasa}`

### Step 3: `run` command

1. Validate arguments (LR + rustsasa → warning+skip, binary existence)
2. Run each selected tool (with progress display via Rich)
3. Merge results into Polars DataFrame by stem
4. Write CSV + config.json
5. Compute and print statistics (R², mean/max relative error)
6. Generate scatter plot

### Step 4: `compare` command

1. Load existing `results.csv` with Polars
2. Compute statistics vs reference tool (default: freesasa)
3. Generate scatter plot

### Step 5: Statistics + scatter plot

**Reuse pattern from `analyze_plots.py:291-359`:**
- R² (Pearson correlation squared)
- Mean/max relative error %
- Scatter plot with y=x line and stats annotation box

## Key Files

| File | Role |
|------|------|
| `benchmarks/scripts/validation.py` | **NEW** - main script |
| `benchmarks/scripts/bench_batch.py` | Pattern: Tool enum, binaries, system info, tempdir |
| `benchmarks/scripts/bench.py:273-317` | FreeSASA output parsing (Total line regex) |
| `benchmarks/scripts/analyze_plots.py:291-359` | Validation stats/plot logic to adapt |

## Constraints

- Input: PDB directory
- FreeSASA: per-file execution (supports SR + LR)
- RustSASA: SR only, skip for LR with warning
- PEP 723 script metadata (`#!/usr/bin/env -S uv run --script`)

## Verification

1. Build zsasa: `zig build -Doptimize=ReleaseFast`
2. Run on small PDB set: `./benchmarks/scripts/validation.py run -i <pdb_dir> -n test`
3. Confirm CSV and scatter plot in `benchmarks/results/validation/test/`
4. Re-analyze: `./benchmarks/scripts/validation.py compare -d benchmarks/results/validation/test`
5. Verify R² ≈ 1.0 for f64 vs FreeSASA SR

---
- [x] **DONE** - Phase complete

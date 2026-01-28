# Plan: RustSASA Batch JSON Processing

## Goal

Extend RustSASA bench fork to support batch JSON processing for fair comparison.

## Current State

RustSASA bench (`benchmarks/external/rustsasa-bench`) has:
- `-J` flag for single JSON file input
- Directory processing for PDB files (with I/O timing mixed)
- `SASA_TIME_US` output for single-file benchmarking

## Requirements

1. Batch JSON directory processing
2. Separate SASA time from I/O time in output
3. Same output format as Zig for benchmark script compatibility
4. Use rayon for file-level parallelism (existing pattern)

## Design

### CLI Interface

```bash
# Single JSON (existing)
rust-sasa -J input.json -n 100 -t 8

# Batch JSON directory (new)
rust-sasa --json-dir ./input_dir/ --output-dir ./output_dir/ -n 100 -t 8
```

### Implementation Approach

Extend `main.rs` with new `--json-dir` option:

```rust
#[derive(Parser, Debug)]
struct Args {
    // ... existing args ...

    /// JSON input directory for batch benchmarking
    #[arg(long)]
    json_dir: Option<String>,

    /// Output directory (required with --json-dir)
    #[arg(long)]
    output_dir: Option<String>,
}
```

### Batch Processing Function

```rust
fn process_json_directory(
    input_dir: &str,
    output_dir: Option<&str>,
    n_points: usize,
    probe_radius: f32,
    threads: isize,
) -> Result<BatchResult, CLIError> {
    use rayon::prelude::*;

    // Scan for .json and .json.gz files
    let files: Vec<_> = std::fs::read_dir(input_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.file_name().to_string_lossy().to_string();
            name.ends_with(".json") || name.ends_with(".json.gz")
        })
        .collect();

    let results = Mutex::new(Vec::new());
    let total_sasa_time = AtomicU64::new(0);

    // Process in parallel
    let total_start = Instant::now();

    files.par_iter().for_each(|entry| {
        let path = entry.path();
        let path_str = path.to_str().unwrap();

        // Time SASA calculation only
        let (sasa_time_us, total_sasa) = process_json_file(
            path_str, n_points, probe_radius, 1  // single-threaded per file
        );

        total_sasa_time.fetch_add(sasa_time_us, Ordering::Relaxed);

        results.lock().unwrap().push(FileResult {
            filename: path_str.to_string(),
            sasa_time_us,
            total_sasa,
        });

        // Write output if output_dir specified
        if let Some(out_dir) = output_dir {
            // Write JSON output
        }
    });

    let total_time = total_start.elapsed();

    // Output in benchmark format
    eprintln!("BATCH_SASA_TIME_MS:{:.2}", total_sasa_time.load(Ordering::Relaxed) as f64 / 1000.0);
    eprintln!("BATCH_TOTAL_TIME_MS:{:.2}", total_time.as_millis());
    eprintln!("BATCH_FILES:{}", files.len());

    Ok(BatchResult { ... })
}
```

### Output Format

Match Zig output for benchmark script compatibility:

```
BATCH_SASA_TIME_MS:1234.56
BATCH_TOTAL_TIME_MS:1567.89
BATCH_FILES:100
```

## Implementation Steps

### Phase 1: Add CLI Options ✅

1. [x] Add `--json-dir` argument
2. [x] Add `--output-dir` argument
3. [x] Validate: json-dir requires output-dir or no output

### Phase 2: Batch Processing ✅

1. [x] Implement `process_json_directory()`
2. [x] Scan directory for JSON files
3. [x] Parallel processing with rayon
4. [x] Accumulate SASA time atomically

### Phase 3: Output & Testing ✅

1. [x] Output `BATCH_*` format to stderr
2. [x] Test with default dataset (9 files, 1356ms)
3. [x] Compare with Zig batch results (~0.3-0.6% difference, expected due to f32 vs f64)
4. [x] Validate SASA values match (within tolerance)

## Files to Modify

| File | Action |
|------|--------|
| `benchmarks/external/rustsasa-bench/src/main.rs` | Modify |

## Success Criteria

1. `rust-sasa --json-dir ./dir/` processes all JSON files
2. Output `BATCH_SASA_TIME_MS` matches Zig format
3. SASA results match single-file mode
4. Scales with thread count

---
- [x] **DONE** - Phase complete

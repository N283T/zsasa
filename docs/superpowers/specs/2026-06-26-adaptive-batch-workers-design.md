# Adaptive Batch Workers Design

## Summary

Add an experimental single-machine batch scheduler mode for very large `zsasa batch` runs that can become I/O-bound on ordinary PCs. The mode should automatically choose a file-worker count below the user-supplied `--threads` maximum when calibration shows that additional file-level concurrency no longer improves throughput. This is intended for large cold-cache datasets such as SwissProt-scale AFDB PDB directories where storage and page-cache behavior can dominate wall-clock time.

This design deliberately does not revive the old `--parallelism=atom|pipeline` feature. Batch mode remains file-level parallelism; the new mode only decides how many file workers to run.

## Motivation

The legacy SwissProt benchmark on the website notes that the 550,122-file SwissProt PDB v6 dataset becomes I/O-bound on a 32 GB M4 system because the 119 GB dataset exceeds available RAM. In that case, `zsasa` keeps low RSS but wall-clock time converges with another bitmask implementation because storage behavior dominates.

The project goal is not cluster-scale distribution. `zsasa` should remain attractive on a single laptop or workstation by avoiding excessive I/O pressure when the machine cannot keep the whole dataset hot in page cache.

## Prior Attempts and Non-Goals

Earlier history added `--parallelism=file|atom|pipeline` in commit `f621250`, with pipeline mode implemented as one I/O prefetch thread feeding atom-level SASA. That mode was later removed in commit `1c4c287`, and the changelog records the removal of `--parallelism=pipeline` and `--parallelism=atom`.

This feature should avoid those paths:

- Do not restore `--parallelism`.
- Do not process one file at a time with atom-level SASA as the primary batch strategy.
- Do not add cluster or multi-machine sharding as a first-class workflow.
- Do not change output schemas or default batch behavior.

## User Interface

Add an opt-in CLI flag:

```bash
zsasa batch input/ output.jsonl \
  --format=jsonl \
  --threads=10 \
  --adaptive-workers
```

Semantics:

- `--threads` remains the maximum file-worker count.
- Without `--adaptive-workers`, existing behavior is unchanged.
- With `--adaptive-workers`, `zsasa` runs a short calibration phase and then processes the remaining files with the selected worker count.
- The selected worker count and calibration summary are printed to stderr unless `--quiet` is set.

Future options such as `--adaptive-workers-sample=N` may be added later, but the first implementation should keep the interface minimal.

## Architecture

### Worker-count calibration

Implement a small calibration phase before the main run:

1. Scan and sort the input directory exactly as current batch mode does.
2. Build the same work item list used by normal parallel batch execution.
3. Choose candidate worker counts from powers of two up to the requested maximum, plus the exact maximum. For example, `--threads=10` tests `1, 2, 4, 8, 10`.
4. For each candidate, process a small disjoint calibration window from the beginning of the work list.
5. Measure wall-clock throughput in files/s for that candidate.
6. Select the smallest worker count whose throughput is close to the best observed throughput.
7. Process the remaining work items with the selected worker count.

The “smallest near-best” rule is important: if 4 workers and 10 workers are within a tolerance, choose 4 to reduce I/O pressure and improve stability on memory-constrained machines.

Initial constants can be conservative and internal:

- Per-candidate minimum: 64 work items.
- Per-candidate maximum: 256 work items.
- Near-best threshold: at least 95% of best observed throughput.
- If the dataset is too small to calibrate meaningfully, fall back to the existing requested worker count.

These constants should be easy to adjust in code after benchmark feedback.

### Reusing current batch code

The implementation should reuse the existing file-level processing path as much as possible:

- Keep `processOneFile`, `processOneSdfMolecule`, JSONL streaming, residue maps, classifiers, bitmask LUT reuse, and error handling semantics unchanged.
- Refactor `runBatchParallel` internally so a caller can run a contiguous work-item range with a fixed worker count.
- Calibration work items should produce normal outputs and normal `FileResult` entries. They are not dry runs and should not be recomputed.
- The final `BatchResult` should include both calibration and post-calibration results in deterministic work-item order.

This avoids doubling I/O for calibration and keeps output behavior predictable.

### Output and progress

For JSONL output, calibration results are streamed as they complete, just like current parallel results. JSONL order is already completion-order rather than sorted order; adaptive mode should not promise stronger ordering.

For per-file JSON/CSV/compact output, calibration writes files normally. The main phase skips already processed work items by using a start offset after the calibration windows.

Progress should cover the full number of work items. During calibration, the progress label may remain generic (`Processing items`) to avoid overcomplicating the first implementation.

### Error handling

A file that fails during calibration counts exactly like a file that fails during the normal parallel run. Throughput selection should use completed calibration items regardless of success status because failed parses still consume I/O and scheduling time.

If calibration cannot run because worker spawning or setup fails, return the underlying error rather than silently falling back. If the dataset is too small for calibration, fall back intentionally and report that decision when not quiet.

## Testing Plan

### Unit tests

- Candidate generation: `max=1`, `max=2`, `max=10`, and `max` not a power of two.
- Worker selection: choose the smallest candidate within 95% of best throughput.
- Small dataset fallback: dataset too small returns the requested worker count.

### Integration tests

Create a temporary directory with several small PDB files copied from existing fixtures. Run adaptive batch with `--threads=4 --adaptive-workers` and compare aggregate success count and total SASA values against normal `--threads=4` batch output.

For JSONL, verify adaptive mode writes parseable JSONL and the same number of successful rows as normal mode. Do not require line ordering.

### Manual performance checks

Use symlinked samples rather than copying large datasets:

- SwissProt sample: 1k, 5k, and optionally 20k files from `/Users/nagaet/pdb/afdb/swissprot_pdb_v6`.
- Human sample or full directory: `/Users/nagaet/pdb/afdb/UP000005640_9606_HUMAN_v6/pdb`.

Compare:

```bash
zsasa batch sample/ out-baseline.jsonl --format=jsonl --threads=10 --use-bitmask --precision=f32 --n-points=128 --quiet --timing
zsasa batch sample/ out-adaptive.jsonl --format=jsonl --threads=10 --adaptive-workers --use-bitmask --precision=f32 --n-points=128 --quiet --timing
```

Record total time, files/s, selected worker count, and peak RSS if available. Full SwissProt should be reserved for final stress testing because it is 119 GB and 550,122 files.

## Documentation

Update the batch guide and CLI reference when the implementation lands:

- Mark `--adaptive-workers` experimental.
- Explain that it is for single-machine, I/O-bound large directory runs.
- Clarify that it does not distribute work across machines.
- Mention that default batch behavior is unchanged unless the flag is set.

## Open Implementation Notes

- The current Human AFDB directory stores files under `pdb/` and `cif/`; existing batch scans are non-recursive, so manual checks should point at `.../pdb` or `.../cif`.
- Calibration windows should be disjoint to avoid recomputation.
- SDF expansion can produce more work items than files; adaptive mode should operate on work items, matching existing batch internals.

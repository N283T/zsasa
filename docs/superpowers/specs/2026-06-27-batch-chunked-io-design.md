# Batch Chunked I/O Design

## Goal

Evaluate and implement two large-batch strategies for SwissProt-scale directory processing:

1. **Simple chunk baseline**: process contiguous ranges of work items in chunks while keeping the existing per-file JSONL write path.
2. **Chunked JSONL writer**: process contiguous ranges in chunks and batch JSONL writes so each worker flushes once per chunk instead of once per file.

The feature is experimental and intended for performance measurement before deciding whether any CLI surface should become stable.

## Background

SwissProt full-run measurements showed that `zsasa batch` uses only about 2.4-2.5 CPU cores on a 10-core machine while processing 550,122 files, so the bottleneck is not SASA arithmetic. The current parallel JSONL path serializes one result at a time behind a mutex and flushes once per file. At SwissProt scale that means more than 550k lock/write/flush cycles.

Simple chunking might help if per-item scheduling or locality is a meaningful cost. Chunked JSONL writing should help if the writer flush/mutex path is the limiting cost. Both strategies should be benchmarkable separately.

## Non-goals

- No distributed or cluster execution.
- No output schema changes.
- No workflow TOML integration for the experimental flags.
- No attempt to preserve JSONL line order in chunked writer mode beyond valid one-line-per-successful-file output.
- No adaptive worker selection in this branch.

## CLI Surface

Add experimental options to `zsasa batch`:

- `--chunk-size=N`: enable chunk-range scheduling for parallel batch runs. `N=0` is invalid; absent means current behavior.
- `--chunked-jsonl`: when JSONL streaming is active, accumulate JSONL lines per worker chunk and flush once per chunk. Requires `--chunk-size=N` and `--format=jsonl`.

The options are intentionally simple because the exact public API can be revised after performance data.

## Architecture

The existing one-file-at-a-time worker remains the default path. New chunked execution uses the same `WorkItem` and `processOneFile` / `processOneSdfMolecule` logic but claims contiguous ranges from an atomic `next_index` by `chunk_size` instead of one item at a time.

For the simple chunk baseline, each file result is still sent through the existing `JsonlStreamWriter.writeResult()` call. This isolates the effect of chunk-range scheduling.

For chunked JSONL writing, each worker serializes successful file results into a thread-local `std.ArrayListUnmanaged(u8)` while processing its claimed range. After finishing the range, it acquires the JSONL stream mutex once, writes the accumulated bytes, flushes once, and resets the buffer.

## Error Handling and Ownership

- `FileResult` ownership remains unchanged: filenames and error messages belong to the result allocator; atom areas and residue maps are cleared after JSONL serialization when backed by the worker arena.
- JSONL serialization/write failures remain warning-only, matching current behavior.
- Chunked JSONL write failures set the same `write_failed` flag used by existing stream writer logic.
- If `--chunked-jsonl` is used without `--format=jsonl` or without `--chunk-size`, argument validation fails.

## Testing

Add tests for:

- CLI parsing of `--chunk-size=N` and `--chunked-jsonl`.
- Invalid CLI combinations.
- Chunk-size helper behavior.
- Simple chunk JSONL output remains parseable and has expected row count.
- Chunked JSONL output remains parseable and has expected row count.

## Benchmark Plan

After implementation, run SwissProt full benchmark with at least:

- Current baseline: `--threads=10`.
- Simple chunk baseline: `--threads=10 --chunk-size=256`.
- Chunked writer: `--threads=10 --chunk-size=256 --chunked-jsonl`.

Use the same legacy-compatible settings: `--format=jsonl --precision=f32 --n-points=100 --use-bitmask --timing`.

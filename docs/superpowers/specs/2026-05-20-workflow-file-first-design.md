# Workflow File-First Batch Execution Design

Date: 2026-05-20
Status: Approved design

## Problem

Workflow batch jobs currently execute job-first: each `[[jobs]]` entry calls the normal batch runner over the full input directory. For common dimer workflows such as chain A, chain B, and chain AB, the same files are scanned, parsed, classified, and prepared repeatedly. This is wasteful because only the chain subset changes between jobs.

## Goals

- Keep the existing workflow TOML schema unchanged.
- Preserve `[[jobs]].chains` semantics and all existing output layouts.
- Avoid repeated directory scans, structure parsing, classifier application, and SR bitmask/LUT preparation across jobs in one workflow run.
- Keep SASA calculations scientifically equivalent to the existing behavior: each job computes SASA on its selected atom set, not by subtracting or approximating from another job.
- Keep the non-workflow batch path as unchanged as practical.

## Non-goals

- Add a new `[[complexes]]` schema.
- Auto-generate pairwise or all-subset chain combinations.
- Share occlusion or neighbor-list state between different chain selections.
- Change JSON, JSONL, CSV, residue-map, or per-file output schemas.

## Proposed Architecture

Add a workflow-specific file-first runner for `zsasa batch --workflow`.

Current execution order:

```text
for each job:
  runBatch(input_dir, job_config)
```

New workflow execution order:

```text
scan input directory once
prepare shared calculation resources once
open one output stream per JSONL job when needed
for each input file:
  parse once
  apply classifier once
  build reusable chain-selection metadata
  for each workflow job:
    derive a selected AtomInput view/copy for the job chains
    run SASA for that selected atom set
    write to the same output destination as today
```

The runner should still produce per-job summaries compatible with the existing workflow completion counts.

## Component Boundaries

### Parsed input preparation

Create a helper that reads and classifies a structure once using a base `BatchConfig` without any job chain filter. For PDB, mmCIF, and BinaryCIF inputs, this preserves the full parsed structure for all workflow jobs. Classifier resources such as external CCD, SDF-derived CCD, and custom classifiers remain loaded once per workflow run.

### Chain selection

Create a helper that builds a job-specific `AtomInput` from a parsed/classified source `AtomInput` and a `[]const []const u8` chain filter.

The selected input owns its arrays and can be deinitialized normally. This is more allocation-heavy than a zero-copy indexed view, but it avoids changing the SASA algorithm interfaces and keeps the first implementation low-risk. If `job.chains` is null, the helper can either duplicate the full input or use a clearly documented fast path that does not mutate the source.

### SASA calculation

Reuse the existing `calculateSasaDispatch` path. Each job still computes SASA independently on its selected atoms. This preserves the current result for chain A, chain B, and chain AB jobs because each job's occluder set remains exactly the selected atom set.

### Output

Keep existing output layout:

- JSONL workflows write one file per job under the workflow output directory, such as `results/chain_a.jsonl`.
- Non-JSONL workflows write per-file outputs under one subdirectory per job, such as `results/chain_a/1abc.json`.
- Stdout behavior for single-job/no-output workflows should remain compatible where practical.

For JSONL, open and maintain one writer per job during the workflow run so the file-first loop can append each file's result to the correct job output.

### SDF handling

SDF does not naturally participate in chain-filtered dimer workflows. The first implementation should keep SDF behavior correct but does not need to optimize SDF molecule parsing across workflow jobs. A safe fallback is acceptable for SDF inputs if it preserves existing outputs.

## Error Handling

- A parse/classifier error for one input file should produce failed results for the affected file while allowing the workflow to continue, matching batch behavior.
- A job that selects no atoms should follow the existing behavior for an empty parsed input if one exists; otherwise return a clear per-file job error.
- Output write failures should mark the affected job/file as failed and include the error name in the result message.
- Invalid workflow configuration remains rejected before file processing begins.

## Testing Plan

- Add unit tests for the AtomInput chain-subset helper, including single chain, multiple chains, null chains, and no matching atoms.
- Add workflow-level tests comparing file-first outputs against the existing job-first behavior for chain A, chain B, and chain AB on a small fixture.
- Cover JSONL residue-map output for a chain-filtered workflow job.
- Run:
  - `zig fmt --check src/`
  - `zig build test`

## Compatibility and Performance Notes

This change should be behavior-preserving. The main expected speedup comes from doing parse/classifier work once per file instead of once per job. LUT preparation is also shared once per workflow run. SASA itself still runs once per job because each chain selection has a different occluder set.

The design intentionally avoids automatic complex generation. Users continue to list only the jobs they want, preventing combinatorial explosion for trimers and larger assemblies.

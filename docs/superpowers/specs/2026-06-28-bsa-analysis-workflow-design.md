# BSA Analysis Workflow Design

## Summary

Add a workflow-level BSA analysis mode for two-chain-group interfaces. The mode computes partner A SASA, partner B SASA, complex SASA, total ΔSASA, and BSA in one reproducible workflow analysis output. BSA and ΔSASA are intentionally separate concepts: ΔSASA is the unhalved loss of accessible area, while BSA is half of total ΔSASA for a two-partner interface.

## Goals

- Add an `[analysis]` section to workflow TOML for BSA-focused batch workflows.
- Keep the existing `[[jobs]]` workflow behavior and output schemas unchanged.
- Use a parser shape supported by the existing minimal TOML parser.
- Produce one JSONL row per input structure for analysis batch output.
- Make the first implementation correctness-first by reusing existing SASA calculation paths for partner A, partner B, and the AB complex.
- Leave room for a later direct interface calculator that skips full A/B/AB output materialization.

## Non-goals

- Do not add a general expression engine for derived workflow outputs.
- Do not implement multi-interface or N-partner automatic enumeration in the first change.
- Do not optimize BSA/ΔSASA by direct cross-interface point-difference calculation yet.
- Do not change normal `zsasa batch --workflow` job output behavior.
- Do not support BSA analysis for `calc --workflow` in this first change.

## Workflow schema

Use a single optional `[analysis]` table in sectioned workflow files:

```toml
[analysis]
type = "bsa"
name = "interface_ab"
partner_a = ["A"]
partner_b = ["B"]
level = "residue"
```

Fields:

- `type`: required for analysis mode. First supported value: `"bsa"`.
- `name`: optional safe output name. Defaults to `"bsa"`.
- `partner_a`: required string array of chain IDs for the first isolated partner.
- `partner_b`: required string array of chain IDs for the second isolated partner.
- `level`: optional output detail level. First implementation supports `"total"` and `"residue"`; default is `"total"`.

The schema uses `partner_a` and `partner_b` instead of nested arrays such as `partners = [["A"], ["B"]]` because the current TOML parser supports arrays of strings but not nested arrays.

Validation:

- `[analysis]` is only accepted in non-legacy sectioned workflows.
- `type = "bsa"` requires non-empty `partner_a` and `partner_b`.
- `name` must pass the same path-safety rules as job names.
- `level` must be `"total"` or `"residue"`.
- Batch BSA analysis requires `[input].dir` or positional input directory.
- JSONL analysis output requires `[output].dir` or positional output path/dir. The output file is `<output_dir>/<analysis.name>.jsonl`.
- If `[analysis]` is present, existing `[[jobs]]` entries are not required for the analysis output.

## Calculation semantics

For each input structure:

1. Parse and classify the full input once where practical.
2. Select atoms for partner A using `partner_a` chains.
3. Select atoms for partner B using `partner_b` chains.
4. Select atoms for the complex using `partner_a + partner_b` chains.
5. Compute SASA for all three selections with the workflow calculation settings.
6. Compute totals:
   - `delta_sasa_total = sasa_partner_a + sasa_partner_b - sasa_complex`
   - `bsa = delta_sasa_total / 2`
7. For `level = "residue"`, compute per-atom ΔSASA by mapping complex atom areas back to partner-selected atom order:
   - partner A atoms: `delta = sasa_partner_a_atom - sasa_complex_atom`
   - partner B atoms: `delta = sasa_partner_b_atom - sasa_complex_atom`
   Then aggregate consecutive atoms by residue using the existing residue-map conventions.

The first implementation may run only in the file-first workflow path and may fallback or reject unsupported combinations with explicit errors. Normal `[[jobs]]` behavior remains available for users who need existing per-job outputs.

## JSONL output schema

BSA analysis JSONL is a separate analysis schema, not an extension of normal SASA JSONL. It intentionally avoids `total_area` and `atom_areas` because those names are ambiguous for BSA.

One row per successfully processed input structure:

```json
{
  "filename": "1abc.cif",
  "analysis": "bsa",
  "name": "interface_ab",
  "partner_a": ["A"],
  "partner_b": ["B"],
  "sasa_partner_a": 1234.5,
  "sasa_partner_b": 1000.0,
  "sasa_complex": 1800.0,
  "delta_sasa_total": 434.5,
  "bsa": 217.25,
  "delta_sasa_level": "residue",
  "residue_chain": ["A", "B"],
  "residue_name": ["ARG", "TYR"],
  "residue_number": [12, 42],
  "residue_insertion_code": ["", ""],
  "residue_delta_sasa": [20.1, 14.5]
}
```

For `level = "total"`, residue arrays are omitted.

Rows for failed files are not written, matching current batch JSONL behavior. Failures contribute to the human-readable workflow summary.

## Architecture

### Workflow manifest

Extend `src/workflow_manifest.zig` with:

- `Analysis` struct stored on `Workflow`.
- Parsing for `[analysis]` via `parseAnalysis`.
- Validation helpers for BSA analysis fields and safe analysis names.
- Tests for valid BSA analysis, missing partners, invalid level, and legacy rejection.

### Analysis result writer

Extend `src/json_writer.zig` with a BSA JSONL serializer that accepts plain slices and optional residue ΔSASA map data. This keeps analysis output formatting out of workflow execution logic.

### Workflow execution

Extend `src/batch.zig` with an analysis path that runs before normal `[[jobs]]` output when `[analysis]` is present. The first implementation should:

- Reuse workflow config/classifier resource loading patterns.
- Scan input files using existing batch discovery.
- Parse/classify each file once.
- Use `copySelectedAtomInput` for partner and complex selections.
- Use `calculatePreparedInputResult` for SASA results.
- Serialize BSA rows to `<output_dir>/<analysis.name>.jsonl`.
- Keep regular job execution unchanged when no `[analysis]` is present.

If both `[analysis]` and `[[jobs]]` are present, the first implementation may run analysis output only and document that jobs are ignored in analysis mode, or run analysis before jobs if implementation remains simple. The initial implementation should prefer a small and clear analysis-only path.

## Testing

Add focused Zig tests:

- Workflow parser accepts `[analysis] type = "bsa"` and stores name, partners, and level.
- Workflow parser rejects BSA analysis without both partners.
- JSONL serializer outputs `bsa`, `delta_sasa_total`, and residue ΔSASA fields with expected names.
- Workflow integration test creates a tiny two-chain fixture, runs `batch --workflow`, and verifies `<output_dir>/<analysis.name>.jsonl` contains `"analysis":"bsa"`, `"bsa"`, `"delta_sasa_total"`, and residue ΔSASA fields for `level = "residue"`.

Run focused checks:

```bash
zig fmt --check src/
zig build test
zig build -Doptimize=ReleaseFast
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa --version
mkdir -p /tmp/zsasa-check
./zig-out/bin/zsasa calc examples/1ubq.pdb /tmp/zsasa-check/output.json
```

Docs check is not required unless website docs are changed. If website workflow docs are updated, run `cd website && npm run build` after dependencies are available.

## Documentation

Update `website/docs/guide/workflows.md` with a BSA analysis workflow example and the distinction between BSA and ΔSASA:

- ΔSASA is unhalved.
- BSA is `delta_sasa_total / 2` for the two-partner analysis.
- BSA analysis JSONL uses analysis-specific fields instead of normal SASA `total_area` and `atom_areas`.

## Future optimization

A later direct calculator can compute ΔSASA only for cross-interface candidate atoms and avoid materializing full partner and complex SASA outputs. This design deliberately keeps the public schema independent from the initial internal implementation so that optimization can be introduced without changing workflow files or JSONL consumers.

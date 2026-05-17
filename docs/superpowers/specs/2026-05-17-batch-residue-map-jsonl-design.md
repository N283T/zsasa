# Batch JSONL Residue Map Design

## Context

Issue #363 requests batch-compatible residue identifiers for large JSONL workflows. The motivating workflow is AFDB homo-dimer/interface analysis, where users need to join SASA output with residue-level data such as pLDDT, buried SASA, and PAE-derived confidence.

Current `zsasa batch --format=jsonl` output contains `filename`, `total_area`, and a bare `atom_areas` array. The array is numerically useful but lacks enough identifiers to map atom SASA values back to residues without reparsing the input structure.

## Decision

Add an opt-in compact residue map for batch JSONL output rather than emitting verbose per-atom JSON objects or a fully opinionated RSA report.

The new option is:

```bash
zsasa batch --format=jsonl --residue-map ...
```

Manifest equivalent:

```toml
residue_map = true
```

The JSONL shape extends the existing row with columnar residue arrays:

```json
{
  "filename": "AF-...-model_v1.cif.zst",
  "total_area": 16055.36,
  "atom_areas": [46.45, 10.81, 1.13],
  "residue_chain": ["A"],
  "residue_name": ["MET"],
  "residue_number": [1],
  "residue_insertion_code": [""],
  "residue_atom_start": [0],
  "residue_atom_count": [3],
  "residue_sasa": [58.39]
}
```

All residue arrays have identical length. `residue_atom_start` and `residue_atom_count` refer to positions in `atom_areas` for the same JSONL row.

## Scope

In scope:

- Add CLI flag `--residue-map` for `batch`.
- Add manifest global `residue_map = true`.
- Emit compact columnar residue arrays only for JSONL output when enabled.
- Preserve existing JSONL output exactly when the option is not enabled.
- Use the chain identifiers already selected by parser settings, including `--auth-chain` / manifest `auth_chain` for mmCIF chain matching and output.
- Include residue SASA as the sum of corresponding atom SASA values.

Out of scope for this change:

- RSA output in batch JSONL.
- Polar/nonpolar summaries in batch JSONL.
- Verbose per-atom metadata objects.
- Changing the default JSONL schema.
- Computing buried SASA or interface metrics directly.

## Data Model

Add a batch-owned residue map representation with columnar fields:

- `residue_chain: []FixedString4`
- `residue_name: []FixedString5`
- `residue_number: []i32`
- `residue_insertion_code: []FixedString4`
- `residue_atom_start: []usize`
- `residue_atom_count: []usize`
- `residue_sasa: []f64`

The map groups consecutive atoms with the same residue key:

```text
(chain_id, residue_num, insertion_code, residue_name)
```

This preserves compact ranges for normal PDB/mmCIF atom ordering. If the same residue appears in non-contiguous blocks, the output contains multiple residue rows with the same identifiers. That behavior keeps `atom_start/count` correct and avoids storing an atom-to-residue index array.

## Data Flow

1. `batch` parses `--residue-map` or manifest `residue_map`.
2. When `output_format == jsonl`, batch already stores `atom_areas` for streaming. With `residue_map` enabled, it also keeps enough residue map data for the current file before the per-file arena resets.
3. After SASA calculation, batch builds the residue map from `AtomInput` metadata and calculated atom areas.
4. The JSONL writer emits either the existing compact row or the extended residue-map row.

If input metadata is insufficient, the file should fail with a clear error only when `--residue-map` is requested. Without `--residue-map`, existing behavior remains unchanged.

## Error Handling

- `--residue-map` with non-JSONL output is rejected with a clear message: residue map is only supported for `--format=jsonl`.
- Missing chain/residue/insertion metadata while `--residue-map` is enabled marks the file as failed with an explanatory per-file error.
- Existing error handling for parse, classifier, SASA, and JSONL write failures remains unchanged.

## Testing

Add focused tests for:

- `batch` argument parsing accepts `--residue-map`.
- manifest parsing accepts `residue_map = true`.
- JSONL writer serializes columnar residue arrays.
- residue map construction groups consecutive atoms and computes `atom_start`, `atom_count`, and `residue_sasa`.
- non-contiguous repeated residue identifiers produce separate contiguous ranges.
- `--residue-map` without `--format=jsonl` is rejected.
- a small batch JSONL smoke test includes residue map fields when requested.

# Input Format Evolution

## Current Format (v1)

```json
{
  "atoms": [
    {"x": 1.0, "y": 2.0, "z": 3.0, "r": 1.8}
  ]
}
```

**Limitation:** 半径 `r` は事前に計算済みが前提。

## Extended Format (v2) - For Radius Classifier

Phase 9 (Radius Classifier) で半径自動判定を行うには、残基・原子情報が必要:

```json
{
  "atoms": [
    {
      "x": 1.0, "y": 2.0, "z": 3.0,
      "r": 1.8,
      "residue_name": "ALA",
      "atom_name": "CA",
      "chain_id": "A",
      "residue_number": 1,
      "element": "C"
    }
  ]
}
```

### Required Fields

| Field | Type | Purpose | Required For |
|-------|------|---------|--------------|
| `x`, `y`, `z` | f64 | 座標 | All |
| `r` | f64 | 半径（既知の場合） | v1 compatibility |

### Optional Fields (v2)

| Field | Type | Purpose | Required For |
|-------|------|---------|--------------|
| `residue_name` | string | 残基名 (ALA, GLY, ...) | Radius classifier |
| `atom_name` | string | 原子名 (CA, N, O, ...) | Radius classifier |
| `chain_id` | string | Chain ID | Output annotation |
| `residue_number` | int | 残基番号 | Output annotation |
| `element` | string | 元素記号 (C, N, O, ...) | Fallback radius |

### Radius Determination Priority

1. **Explicit `r`** - 入力に `r` があればそれを使用
2. **Classifier lookup** - `residue_name` + `atom_name` から半径を取得
3. **Element fallback** - `element` から汎用半径を取得
4. **Error** - 半径を決定できない場合はエラー

## Backward Compatibility

v2フォーマットはv1と後方互換:
- `r` が存在すればそれを使用（v1動作）
- `r` が無い場合のみclassifierを使用

```zig
// Pseudo-code
fn getRadius(atom: AtomJson, classifier: ?Classifier) !f64 {
    // Priority 1: explicit radius
    if (atom.r) |r| return r;

    // Priority 2: classifier lookup
    if (classifier) |c| {
        if (atom.residue_name) |res| {
            if (atom.atom_name) |atm| {
                if (c.getRadius(res, atm)) |r| return r;
            }
        }
    }

    // Priority 3: element fallback
    if (atom.element) |elem| {
        return getElementRadius(elem);
    }

    return error.CannotDetermineRadius;
}
```

## Phase Dependencies

```
Phase 9 (Radius Classifier)
    └── Requires: Extended input format (v2)
    └── Enables: Direct mmCIF/PDB → SASA calculation

Phase 10 (mmCIF Input)
    └── Requires: Extended input format (v2)
    └── Requires: Phase 9 (for radius assignment)
    └── Enables: Direct mmCIF file input

Phase 11 (Lee-Richards)
    └── No input format changes needed
    └── Uses existing AtomInput (x, y, z, r)
```

## Implementation Order (Suggested)

1. **Extend types.zig** - Add optional fields to AtomInput
2. **Update json_parser.zig** - Parse extended format
3. **Implement classifier.zig** - Radius lookup logic
4. **Update main.zig** - CLI options for classifier
5. **Test with extended JSON** - Verify backward compatibility

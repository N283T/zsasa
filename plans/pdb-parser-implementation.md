# Plan: PDB Parser Implementation

## Goal

既存の mmCIF/JSON パーサーと整合性のある PDB パーサーを実装する。

## Architecture

### 既存構造との整合性

```
Input Files → Parser → AtomInput → SASA Calculator
   .json   →  json_parser.zig  ─┐
   .cif    →  mmcif_parser.zig ─┼→ AtomInput (types.zig)
   .pdb    →  pdb_parser.zig  ─┘      ↓
                                  SASA calculation
```

### AtomInput 構造体 (出力)

```zig
pub const AtomInput = struct {
    x: []const f64,           // 必須: X座標
    y: []const f64,           // 必須: Y座標
    z: []const f64,           // 必須: Z座標
    r: []const f64,           // 必須: 原子半径 (VdW)

    residue: ?[]const []const u8,        // 残基名 (ALA, GLY...)
    atom_name: ?[]const []const u8,      // 原子名 (CA, CB...)
    element: ?[]const u8,                 // 元素番号
    chain_id: ?[]const []const u8,       // チェーンID
    residue_num: ?[]const i32,           // 残基番号
    insertion_code: ?[]const []const u8, // 挿入コード

    allocator: std.mem.Allocator,
};
```

## Implementation

### Phase 1: Core Parser

**File: `src/pdb_parser.zig`**

```zig
pub const PdbParser = struct {
    allocator: std.mem.Allocator,

    // Configuration (mmCIF と同じ)
    atom_only: bool = false,           // HETATM を除外
    first_alt_loc_only: bool = true,   // 最初の alt loc のみ
    model_num: ?u32 = null,            // 特定モデル番号
    chain_filter: ?[]const []const u8 = null,

    pub fn init(allocator: std.mem.Allocator) PdbParser
    pub fn parse(self: *PdbParser, source: []const u8) !AtomInput
    pub fn parseFile(self: *PdbParser, path: []const u8) !AtomInput
};
```

### Phase 2: PDB Record Parsing

**PDB フォーマット (固定幅)**

| Field | Columns | Offset | Length | Description |
|-------|---------|--------|--------|-------------|
| Record | 1-6 | 0 | 6 | ATOM/HETATM |
| Serial | 7-11 | 6 | 5 | 原子番号 |
| Name | 13-16 | 12 | 4 | 原子名 |
| AltLoc | 17 | 16 | 1 | 代替位置 |
| ResName | 18-20 | 17 | 3 | 残基名 |
| ChainID | 22 | 21 | 1 | チェーンID |
| ResSeq | 23-26 | 22 | 4 | 残基番号 |
| InsCode | 27 | 26 | 1 | 挿入コード |
| X | 31-38 | 30 | 8 | X座標 |
| Y | 39-46 | 38 | 8 | Y座標 |
| Z | 47-54 | 46 | 8 | Z座標 |
| Occupancy | 55-60 | 54 | 6 | 占有率 |
| TempFactor | 61-66 | 60 | 6 | B因子 |
| Element | 77-78 | 76 | 2 | 元素記号 |

**Parsing Functions**

```zig
fn parseAtomRecord(line: []const u8) !?AtomRecord {
    // 1. Record type check (ATOM/HETATM)
    // 2. Fixed-width field extraction
    // 3. Coordinate parsing (X, Y, Z)
    // 4. Element symbol → VdW radius
}

fn parseElement(line: []const u8) ?Element {
    // 1. Try columns 77-78 (element symbol)
    // 2. Fallback: infer from atom name (columns 13-16)
}
```

### Phase 3: Multi-Model Support

```zig
// MODEL/ENDMDL handling
fn isModelRecord(line: []const u8) bool
fn isEndModelRecord(line: []const u8) bool
fn parseModelNumber(line: []const u8) ?u32
```

### Phase 4: CLI Integration

**File: `src/main.zig`**

```zig
fn detectInputFormat(path: []const u8) InputFormat {
    if (endsWith(path, ".pdb") or endsWith(path, ".PDB") or
        endsWith(path, ".ent") or endsWith(path, ".ENT"))
        return .pdb;
    // ... existing logic
}

fn readInputFile(...) !AtomInput {
    return switch (format) {
        .pdb => {
            var parser = pdb_parser.PdbParser.init(allocator);
            parser.model_num = args.model_num;
            parser.chain_filter = parseChainFilter(...);
            return parser.parseFile(path);
        },
        // ... existing cases
    };
}
```

## Edge Cases

1. **短い行**: 80文字未満の行も許容 (元素記号フィールドがない場合あり)
2. **空白トリム**: 原子名・残基名の前後空白
3. **水素検出**: 元素記号 H/D または原子名から推定
4. **元素推定**: 原子名から元素を推定 (CA→C, FE→Fe)
5. **Alt Location**: 最初のもののみ採用 (FreeSASA 方式)
6. **Hybrid-36**: 将来対応 (99,999超の原子番号)

## Test Plan

1. **Unit tests**: 各フィールド抽出関数
2. **Integration tests**: 実際の PDB ファイル
3. **Comparison**: FreeSASA C との結果比較

## Files to Create/Modify

| File | Action |
|------|--------|
| `src/pdb_parser.zig` | **新規作成** |
| `src/main.zig` | 修正 (format detection, readInputFile) |
| `src/types.zig` | 修正 (InputFormat enum に .pdb 追加) |
| `build.zig` | 修正 (pdb_parser.zig をビルドに追加) |

## Verification

```bash
# Build
zig build -Doptimize=ReleaseFast

# Test with PDB file
./zig-out/bin/freesasa_zig test.pdb output.json

# Compare with FreeSASA C
freesasa test.pdb --format=json > freesasa_output.json
diff output.json freesasa_output.json
```

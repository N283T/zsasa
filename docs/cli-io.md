# CLI・入出力詳解

## CLI引数パース

### ファイル: `src/main.zig`

### Args構造体

```zig
const Args = struct {
    // 入出力
    input_file: ?[]const u8 = null,
    output_file: ?[]const u8 = null,

    // アルゴリズムパラメータ
    algorithm: Algorithm = .sr,     // sr or lr
    n_threads: usize = 0,           // 0 = 自動検出
    probe_radius: f64 = 1.4,        // デフォルト: 1.4 Å
    n_points: u32 = 100,            // デフォルト: 100 (SR)
    n_slices: u32 = 20,             // デフォルト: 20 (LR)

    // 分類器
    classifier: ?ClassifierType = null,  // naccess, protor, oons
    config_file: ?[]const u8 = null,     // カスタム設定

    // チェーン/モデル選択
    chain: ?[]const u8 = null,           // ラベルチェーンID
    model: u32 = 1,                      // モデル番号
    auth_chain: ?[]const u8 = null,      // authチェーンID

    // 解析オプション
    per_residue: bool = false,           // 残基単位出力
    rsa: bool = false,                   // RSA計算
    polar: bool = false,                 // 極性/非極性
    interface: bool = false,             // 界面SASA

    // 出力オプション
    output_format: OutputFormat = .json,
    quiet: bool = false,

    // モードフラグ
    validate_only: bool = false,
    show_help: bool = false,
    show_version: bool = false,
};
```

### パース処理

```zig
fn parseArgs(args: []const [:0]const u8) !Args {
    var result = Args{};
    var i: usize = 1;  // args[0]はプログラム名

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        // ロングオプション（--key=value形式）
        if (std.mem.startsWith(u8, arg, "--threads=")) {
            result.n_threads = try parseUsize(arg["--threads=".len..]);
        } else if (std.mem.startsWith(u8, arg, "--probe-radius=")) {
            result.probe_radius = try parseFloat(arg["--probe-radius=".len..]);
            if (result.probe_radius <= 0 or result.probe_radius > 10) {
                return error.InvalidProbeRadius;
            }
        } else if (std.mem.startsWith(u8, arg, "--n-points=")) {
            result.n_points = try parseU32(arg["--n-points=".len..]);
            if (result.n_points < 1 or result.n_points > 10000) {
                return error.InvalidNPoints;
            }
        } else if (std.mem.startsWith(u8, arg, "--format=")) {
            result.output_format = parseOutputFormat(arg["--format=".len..]);
        }
        // フラグオプション
        else if (std.mem.eql(u8, arg, "--validate")) {
            result.validate_only = true;
        } else if (std.mem.eql(u8, arg, "--quiet") or std.mem.eql(u8, arg, "-q")) {
            result.quiet = true;
        } else if (std.mem.eql(u8, arg, "--help") or std.mem.eql(u8, arg, "-h")) {
            result.show_help = true;
        } else if (std.mem.eql(u8, arg, "--version") or std.mem.eql(u8, arg, "-V")) {
            result.show_version = true;
        }
        // 位置引数
        else if (!std.mem.startsWith(u8, arg, "-")) {
            if (result.input_file == null) {
                result.input_file = arg;
            } else if (result.output_file == null) {
                result.output_file = arg;
            }
        } else {
            return error.UnknownOption;
        }
    }

    return result;
}
```

### ヘルプ出力

```zig
fn printHelp() void {
    const help_text =
        \\Usage: freesasa_zig [OPTIONS] <input> [output.json]
        \\
        \\Calculate Solvent Accessible Surface Area (SASA).
        \\Input formats: JSON, PDB (.pdb), mmCIF (.cif, .cif.gz)
        \\
        \\Algorithm:
        \\  --algorithm=ALGO   sr (Shrake-Rupley) or lr (Lee-Richards) (default: sr)
        \\  --threads=N        Number of threads (default: auto-detect)
        \\  --probe-radius=R   Probe radius in Angstroms (default: 1.4)
        \\  --n-points=N       Test points per atom, SR only (default: 100)
        \\  --n-slices=N       Slices per atom, LR only (default: 20)
        \\
        \\Classifier:
        \\  --classifier=TYPE  naccess, protor, or oons
        \\  --config=FILE      Custom classifier config file
        \\
        \\Chain/Model:
        \\  --chain=ID         Filter by label chain ID
        \\  --model=N          Select model number (default: 1)
        \\  --auth-chain=ID    Filter by auth chain ID
        \\
        \\Analysis:
        \\  --per-residue      Output per-residue SASA
        \\  --rsa              Calculate RSA (implies --per-residue)
        \\  --polar            Show polar/nonpolar summary
        \\  --interface        Calculate interface SASA
        \\
        \\Output:
        \\  --format=FORMAT    json, compact, or csv (default: json)
        \\  --validate         Validate input only
        \\  -q, --quiet        Suppress progress output
        \\  -h, --help         Show this help message
        \\  -V, --version      Show version
        \\
    ;
    std.debug.print("{s}", .{help_text});
}
```

---

## 入力フォーマット

### ファイル: `src/json_parser.zig`

### JSON構造

```json
{
  "x": [1.0, 2.0, 3.0, ...],
  "y": [1.0, 2.0, 3.0, ...],
  "z": [1.0, 2.0, 3.0, ...],
  "r": [1.7, 1.55, 1.52, ...]
}
```

**フィールド説明:**
- `x`, `y`, `z`: 原子座標（Å単位）
- `r`: van der Waals半径（Å単位）

### パース処理

```zig
pub fn parseAtomInput(allocator: Allocator, json_string: []const u8) !AtomInput {
    const parsed = try std.json.parseFromSlice(
        struct {
            x: []f64,
            y: []f64,
            z: []f64,
            r: []f64,
        },
        allocator,
        json_string,
        .{},
    );

    return AtomInput{
        .x = parsed.value.x,
        .y = parsed.value.y,
        .z = parsed.value.z,
        .r = parsed.value.r,
    };
}
```

### バリデーション

```zig
const MAX_RADIUS_ANGSTROMS: f64 = 100.0;

pub const ValidationError = struct {
    message: []const u8,
    index: ?usize = null,
    value: ?f64 = null,
};

pub const ValidationResult = struct {
    valid: bool,
    errors: []ValidationError,
    allocator: Allocator,

    pub fn deinit(self: *ValidationResult) void {
        self.allocator.free(self.errors);
    }
};

pub fn validateInput(allocator: Allocator, input: AtomInput) !ValidationResult {
    var errors = std.ArrayList(ValidationError).init(allocator);
    errdefer errors.deinit();

    const n = input.x.len;

    // 1. 空入力チェック
    if (n == 0) {
        try errors.append(.{ .message = "Input is empty" });
        return .{ .valid = false, .errors = try errors.toOwnedSlice(), .allocator = allocator };
    }

    // 2. 配列長チェック
    if (input.y.len != n or input.z.len != n or input.r.len != n) {
        try errors.append(.{ .message = "Array lengths do not match" });
    }

    // 3. 座標の有限性チェック
    for (input.x, 0..) |x, i| {
        if (!std.math.isFinite(x)) {
            try errors.append(.{
                .message = "x coordinate is not finite",
                .index = i,
                .value = x,
            });
        }
    }
    // y, zも同様...

    // 4. 半径チェック
    for (input.r, 0..) |r, i| {
        if (r <= 0) {
            try errors.append(.{
                .message = "radius must be positive",
                .index = i,
                .value = r,
            });
        } else if (r > MAX_RADIUS_ANGSTROMS) {
            try errors.append(.{
                .message = "radius exceeds maximum (100 Å)",
                .index = i,
                .value = r,
            });
        }
    }

    return .{
        .valid = errors.items.len == 0,
        .errors = try errors.toOwnedSlice(),
        .allocator = allocator,
    };
}
```

### エラーメッセージ出力

```zig
pub fn printValidationErrors(errors: []const ValidationError) void {
    std.debug.print("Input validation failed with {d} error(s):\n", .{errors.len});

    for (errors, 0..) |err, i| {
        // 10件までに制限
        if (i >= 10) {
            std.debug.print("  ... and {d} more errors\n", .{errors.len - 10});
            break;
        }

        if (err.index) |idx| {
            if (err.value) |val| {
                std.debug.print("  - Atom {d}: {s} (value: {d})\n", .{ idx, err.message, val });
            } else {
                std.debug.print("  - Atom {d}: {s}\n", .{ idx, err.message });
            }
        } else {
            std.debug.print("  - {s}\n", .{err.message});
        }
    }
}
```

**出力例:**
```
Input validation failed with 3 error(s):
  - Atom 42: radius must be positive (value: -1.5)
  - Atom 100: x coordinate is not finite (value: nan)
  - Atom 101: x coordinate is not finite (value: inf)
```

---

## PDB/mmCIF入力

### ファイル: `src/pdb_parser.zig`, `src/mmcif_parser.zig`

構造ファイルの拡張子で形式を自動判別:
- `.pdb` → PDBパーサー
- `.cif`, `.cif.gz` → mmCIFパーサー

```zig
fn detectInputFormat(path: []const u8) InputFormat {
    if (std.mem.endsWith(u8, path, ".pdb")) return .pdb;
    if (std.mem.endsWith(u8, path, ".cif")) return .mmcif;
    if (std.mem.endsWith(u8, path, ".cif.gz")) return .mmcif_gz;
    return .json;
}
```

### チェーン/モデル選択

```bash
# ラベルチェーンID（mmCIF label_asym_id）
--chain=A

# モデル番号（NMR構造など）
--model=1

# authチェーンID（PDB互換のchain ID）
--auth-chain=A
```

---

## 解析機能

### ファイル: `src/analysis.zig`

### 残基単位集計（--per-residue）

```zig
pub const ResidueSasa = struct {
    chain_id: []const u8,
    residue_name: []const u8,
    residue_num: i32,
    insertion_code: []const u8,
    sasa: f64,
    rsa: ?f64 = null,  // --rsa使用時
};
```

**出力例:**
```
Per-residue SASA:
Chain  Res   Num   SASA(Å²)
    A  MET     1     198.52
    A  LYS     2     142.31
    A  ALA     3      45.67
```

### RSA計算（--rsa）

RSA = SASA / MaxSASA

MaxSASA値はTien et al. (2013)の参照値を使用:
- ALA: 129.0 Å²
- ARG: 274.0 Å²
- etc.

**出力例:**
```
Per-residue SASA with RSA:
Chain  Res   Num   SASA(Å²)    RSA
    A  MET     1     198.52  0.958
    A  LYS     2     142.31  0.701
```

### 極性/非極性分類（--polar）

残基を極性/非極性に分類:
- **極性**: D, E, H, K, N, Q, R, S, T, Y
- **非極性**: A, C, F, G, I, L, M, P, V, W

**出力例:**
```
Polar/Nonpolar Summary:
  Polar:     2,345.67 Å² (45.2%)
  Nonpolar:  2,845.23 Å² (54.8%)
  Total:     5,190.90 Å²
```

### 界面SASA（--interface）

複合体界面 = Σ(個別チェーンSASA) - 複合体SASA

**出力例:**
```
Interface SASA Analysis:
Chain  Isolated(Å²)  Complex(Å²)  Buried(Å²)
    A      4,234.56      3,890.12      344.44
    B      4,567.89      4,123.45      444.44

Total buried at interface: 788.88 Å²
```

---

## 出力フォーマット

### ファイル: `src/json_writer.zig`

### OutputFormat enum

```zig
pub const OutputFormat = enum {
    json,     // Pretty-printed JSON（デフォルト）
    compact,  // Single-line JSON
    csv,      // CSV形式
};
```

### JSON出力（Pretty）

```zig
pub fn sasaResultToJsonPretty(allocator: Allocator, result: SasaResult) ![]u8 {
    const output = JsonOutput{
        .total_area = result.total_area,
        .atom_areas = result.atom_areas,
    };

    return std.json.stringify(allocator, output, .{
        .whitespace = .indent_2,
    });
}
```

**出力例:**
```json
{
  "total_area": 18923.28,
  "atom_areas": [
    32.47,
    0.25,
    15.82,
    ...
  ]
}
```

### JSON出力（Compact）

```zig
pub fn sasaResultToJson(allocator: Allocator, result: SasaResult) ![]u8 {
    const output = JsonOutput{
        .total_area = result.total_area,
        .atom_areas = result.atom_areas,
    };

    return std.json.stringify(allocator, output, .{});
}
```

**出力例:**
```json
{"total_area":18923.28,"atom_areas":[32.47,0.25,15.82,...]}
```

### CSV出力

```zig
pub fn sasaResultToCsv(allocator: Allocator, result: SasaResult) ![]u8 {
    var list = std.ArrayListUnmanaged(u8){};
    errdefer list.deinit(allocator);

    const writer = list.writer(allocator);

    // ヘッダー
    try writer.writeAll("atom_index,area\n");

    // 各原子のSASA
    for (result.atom_areas, 0..) |area, i| {
        try writer.print("{d},{d:.6}\n", .{ i, area });
    }

    // 合計
    try writer.print("total,{d:.6}\n", .{result.total_area});

    return list.toOwnedSlice(allocator);
}
```

**出力例:**
```csv
atom_index,area
0,32.470000
1,0.250000
2,15.820000
...
total,18923.280000
```

### ファイル書き込み

```zig
pub fn writeSasaResultWithFormat(
    allocator: Allocator,
    result: SasaResult,
    path: []const u8,
    format: OutputFormat,
) !void {
    const output_str = switch (format) {
        .json => try sasaResultToJsonPretty(allocator, result),
        .compact => try sasaResultToJson(allocator, result),
        .csv => try sasaResultToCsv(allocator, result),
    };
    defer allocator.free(output_str);

    const file = try std.fs.cwd().createFile(path, .{});
    defer file.close();

    try file.writeAll(output_str);
}
```

---

## メイン処理フロー

```zig
pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    // 1. 引数パース
    const args = try parseArgs(std.os.argv);

    // 2. ヘルプ/バージョン表示
    if (args.show_help) {
        printHelp();
        return;
    }
    if (args.show_version) {
        printVersion();
        return;
    }

    // 3. 入力ファイル読み込み
    const input_file = args.input_file orelse {
        std.debug.print("Error: No input file specified\n", .{});
        return error.NoInputFile;
    };
    const json_content = try std.fs.cwd().readFileAlloc(allocator, input_file, 100 * 1024 * 1024);
    defer allocator.free(json_content);

    // 4. JSONパース
    var atom_input = try json_parser.parseAtomInput(allocator, json_content);
    defer atom_input.deinit();

    // 5. バリデーション
    var validation = try json_parser.validateInput(allocator, atom_input);
    defer validation.deinit();

    if (!validation.valid) {
        json_parser.printValidationErrors(validation.errors);
        return error.ValidationFailed;
    }

    // 6. バリデーションのみモード
    if (args.validate_only) {
        if (!args.quiet) {
            std.debug.print("Input validation passed ({d} atoms)\n", .{atom_input.len()});
        }
        return;
    }

    // 7. SASA計算
    if (!args.quiet) {
        std.debug.print("Calculating SASA for {d} atoms...\n", .{atom_input.len()});
    }

    var result = try shrake_rupley.calculateSasa(allocator, atom_input, .{
        .n_points = args.n_points,
        .probe_radius = args.probe_radius,
        .n_threads = args.n_threads,
    });
    defer result.deinit();

    // 8. 結果出力
    if (args.output_file) |output_path| {
        try json_writer.writeSasaResultWithFormat(
            allocator,
            result,
            output_path,
            args.output_format,
        );
        if (!args.quiet) {
            std.debug.print("Results written to: {s}\n", .{output_path});
        }
    } else {
        // 標準出力
        const output = try json_writer.sasaResultToJsonPretty(allocator, result);
        defer allocator.free(output);
        std.io.getStdOut().writer().writeAll(output) catch {};
    }
}
```

---

## エラーハンドリング一覧

| エラー | 発生箇所 | 原因 |
|--------|----------|------|
| `InvalidJson` | json_parser | JSON構文エラー |
| `MissingField` | json_parser | 必須フィールドがない |
| `ArrayLengthMismatch` | json_parser | x,y,z,rの長さが不一致 |
| `EmptyInput` | json_parser | 原子数が0 |
| `InvalidRadius` | json_parser | 半径が0以下または100超 |
| `InvalidCoordinate` | json_parser | 座標がNaN/Inf |
| `InvalidProbeRadius` | main | プローブ半径が範囲外 |
| `InvalidNPoints` | main | ポイント数が範囲外 |
| `NoInputFile` | main | 入力ファイル未指定 |
| `FileNotFound` | main | 入力ファイルが存在しない |
| `OutOfMemory` | 各所 | メモリ不足 |

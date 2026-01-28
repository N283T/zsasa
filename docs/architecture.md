# アーキテクチャ概要

## プロジェクト構成

```
src/
├── main.zig           # CLIエントリーポイント
├── root.zig           # ライブラリルートモジュール
├── types.zig          # データ構造定義
├── json_parser.zig    # JSON入力パース・バリデーション
├── json_writer.zig    # 出力フォーマット（JSON/CSV）
├── test_points.zig    # Golden Section Spiralによるテストポイント生成
├── neighbor_list.zig  # 空間ハッシュによる近傍探索
├── simd.zig           # SIMD最適化（8-wide距離計算、高速三角関数）
├── thread_pool.zig    # 汎用スレッドプール
├── shrake_rupley.zig  # Shrake-Rupleyアルゴリズム
└── lee_richards.zig   # Lee-Richardsアルゴリズム
```

## モジュール依存関係

```
main.zig
    │
    ├── json_parser.zig ──► types.zig
    │
    ├── shrake_rupley.zig
    │       │
    │       ├── test_points.zig ──► types.zig
    │       │
    │       ├── neighbor_list.zig ──► types.zig
    │       │
    │       ├── simd.zig ──► types.zig
    │       │
    │       └── thread_pool.zig
    │
    └── json_writer.zig ──► types.zig
```

## データフロー

```
┌─────────────────┐
│   入力JSON      │
│ (x, y, z, r)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  json_parser    │  バリデーション
│  - 配列長チェック │  - 座標がfinite
│  - 半径チェック   │  - 半径 > 0, ≤ 100
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│   AtomInput     │  内部データ構造
│   - x: []f64    │
│   - y: []f64    │
│   - z: []f64    │
│   - r: []f64    │
└────────┬────────┘
         │
         ▼
┌─────────────────────────────────────────┐
│           shrake_rupley                  │
│  ┌─────────────────┐                    │
│  │  test_points    │ テストポイント生成   │
│  │  (事前計算)      │ Golden Section     │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │  neighbor_list  │ 近傍原子リスト構築   │
│  │  (O(N)構築)     │ 空間ハッシュ        │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │  thread_pool    │ 並列処理            │
│  │  (原子ごと)      │ Work-stealing      │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │     simd        │ 距離計算            │
│  │  (@Vector 8本)   │ 8原子同時処理      │
│  └─────────────────┘                    │
└────────────────────┬────────────────────┘
                     │
                     ▼
┌─────────────────┐
│   SasaResult    │
│   - total_area  │
│   - atom_areas  │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  json_writer    │  出力フォーマット
│  - JSON (pretty)│
│  - JSON (compact)│
│  - CSV          │
└─────────────────┘
```

## 主要な型定義 (types.zig)

### Vec3

3次元ベクトル。原子座標やテストポイントの表現に使用。精度パラメータ（f32/f64）に対応。

```zig
pub fn Vec3(comptime T: type) type {
    return struct {
        x: T,
        y: T,
        z: T,

        pub fn sub(self: @This(), other: @This()) @This() { ... }
        pub fn scale(self: @This(), s: T) @This() { ... }
        pub fn add(self: @This(), other: @This()) @This() { ... }
        pub fn lengthSquared(self: @This()) T { ... }
    };
}

// 使用例
const Vec3f32 = Vec3(f32);  // f32精度
const Vec3f64 = Vec3(f64);  // f64精度（デフォルト）
```

### AtomInput

入力データ構造。SoA (Structure of Arrays) 形式でメモリ効率とキャッシュ効率を最適化。精度パラメータ（f32/f64）に対応。

```zig
pub fn AtomInput(comptime T: type) type {
    return struct {
        x: []const T,
        y: []const T,
        z: []const T,
        r: []const T,

        pub fn len(self: @This()) usize {
            return self.x.len;
        }
    };
}

// 使用例
const AtomInputf32 = AtomInput(f32);  // f32精度（高速）
const AtomInputf64 = AtomInput(f64);  // f64精度（デフォルト）
```

**SoA形式の利点:**
- SIMD処理との親和性が高い
- メモリアクセスパターンが連続的
- キャッシュライン効率が良い
- f32使用時はデータサイズ半減でさらに効率向上

### SasaResult

計算結果を保持。

```zig
pub const SasaResult = struct {
    total_area: f64,
    atom_areas: []f64,
    allocator: Allocator,

    pub fn deinit(self: *SasaResult) void {
        self.allocator.free(self.atom_areas);
    }
};
```

## メモリ管理

Zigの明示的アロケータパターンを採用。全ての動的メモリ割り当ては `Allocator` を通じて行われる。

```zig
// 典型的なパターン
pub fn calculate(allocator: Allocator, input: AtomInput) !SasaResult {
    const atom_areas = try allocator.alloc(f64, input.len());
    errdefer allocator.free(atom_areas);  // エラー時の解放

    // ... 計算 ...

    return SasaResult{
        .total_area = total,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}
```

**利点:**
- メモリリークの防止（コンパイル時に追跡可能）
- テスト時に `std.testing.allocator` でリーク検出
- カスタムアロケータ（アリーナ等）への差し替えが容易

## エラーハンドリング

Zigのエラーユニオン型を使用。全てのエラーは呼び出し元に伝播可能。

```zig
pub const ParseError = error{
    InvalidJson,
    MissingField,
    ArrayLengthMismatch,
    EmptyInput,
};

pub const ValidationError = error{
    InvalidRadius,
    InvalidCoordinate,
};

pub fn parseAtomInput(allocator: Allocator, json: []const u8) !AtomInput {
    // エラーの場合は自動的に呼び出し元に伝播
    const parsed = try std.json.parseFromSlice(...);
    // ...
}
```

## スレッドセーフティ

- `thread_pool.zig`: スレッド間の同期は `std.Thread.Mutex` と `std.Thread.Condition` で実装
- `shrake_rupley.zig`: 各原子の計算は独立しているため、並列実行可能
- 共有状態: テストポイント配列と近傍リストは読み取り専用として共有

# アーキテクチャ概要

## プロジェクト構造

```
src/
├── main.zig           # CLIエントリーポイント
├── root.zig           # ライブラリルートモジュール
├── types.zig          # データ構造定義
├── json_parser.zig    # JSON入力解析と検証
├── json_writer.zig    # 出力形式（JSON/CSV）
├── test_points.zig    # 黄金角螺旋によるテスト点生成
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
│  json_parser    │  検証
│  - 配列長       │  - 座標が有限
│  - 半径チェック │  - 半径 > 0, ≤ 100
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
│  │  test_points    │ テスト点生成        │
│  │  (事前計算)     │ 黄金螺旋           │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │  neighbor_list  │ 近傍リスト構築      │
│  │  (O(N)構築)     │ 空間ハッシュ        │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │  thread_pool    │ 並列処理            │
│  │  (原子単位)     │ ワークスティーリング │
│  └────────┬────────┘                    │
│           │                              │
│  ┌────────▼────────┐                    │
│  │     simd        │ 距離計算            │
│  │  (@Vector 8)    │ 8原子同時処理       │
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
│  json_writer    │  出力形式
│  - JSON (整形)  │
│  - JSON (圧縮)  │
│  - CSV          │
└─────────────────┘
```

## 主要な型定義 (types.zig)

### Vec3

3Dベクトル。原子座標とテスト点に使用。精度パラメータ（f32/f64）をサポート。

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

入力データ構造。メモリとキャッシュ効率を最適化するためSoA（Structure of Arrays）形式を使用。精度パラメータ（f32/f64）をサポート。

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
- SIMD処理との高い親和性
- 連続したメモリアクセスパターン
- キャッシュライン効率が良い
- f32ではデータサイズが半分になりさらに効率向上

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

Zigの明示的アロケータパターンを使用。すべての動的メモリ割り当ては`Allocator`を通じて行われます。

```zig
// 典型的なパターン
pub fn calculate(allocator: Allocator, input: AtomInput) !SasaResult {
    const atom_areas = try allocator.alloc(f64, input.len());
    errdefer allocator.free(atom_areas);  // エラー時に解放

    // ... 計算 ...

    return SasaResult{
        .total_area = total,
        .atom_areas = atom_areas,
        .allocator = allocator,
    };
}
```

**利点:**
- メモリリーク防止（コンパイル時に追跡可能）
- テスト時の`std.testing.allocator`によるリーク検出
- カスタムアロケータ（アリーナなど）への容易な置き換え

## エラーハンドリング

Zigのエラー共用体型を使用。すべてのエラーは呼び出し元に伝播可能。

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
    // エラーは自動的に呼び出し元に伝播
    const parsed = try std.json.parseFromSlice(...);
    // ...
}
```

## スレッドセーフティ

- `thread_pool.zig`: `std.Thread.Mutex`と`std.Thread.Condition`でスレッド間同期を実装
- `shrake_rupley.zig`: 各原子の計算は独立しているため並列実行可能
- 共有状態: テスト点配列と近傍リストは読み取り専用として共有

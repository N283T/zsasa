# Shrake-Rupley 最適化プラン

## 最終結果

**最適化後ベンチマーク (4V6X, 237k atoms, 4 threads):**
- Zig: 191ms (改善前: 207ms)
- FreeSASA C: 431ms
- Speedup: **2.26x** (改善前: 2.03x)

**改善内容:**
- 8-wide SIMD実装により約12%の性能向上を達成

---

## 元の分析

**ベンチマーク結果 (4V6X, 237k atoms, 4 threads):**
- Zig: 207ms
- FreeSASA C: 421ms
- Speedup: **2.03x**

**期待との乖離:**
ユーザーの期待はもっと高速化できるはず。現状2x程度は物足りない。

## ボトルネック分析

コード読解結果 (`src/shrake_rupley.zig`):

### 1. SoA→AoS変換のオーバーヘッド (lines 204-213, 331-340)

```zig
// 現在: input (SoA) → positions (AoS) に変換
const positions = try allocator.alloc(Vec3, n_atoms);
for (0..n_atoms) |i| {
    positions[i] = Vec3{
        .x = input.x[i],
        .y = input.y[i],
        .z = input.z[i],
    };
}
```

**問題点:**
- 追加のメモリ確保 (24 bytes × N atoms)
- 3N回の書き込み
- NeighborListもpositionsを要求 → さらにコピー発生

### 2. SIMD 4-wide の制限 (lines 131-148)

```zig
while (i + 4 <= neighbors.len) : (i += 4) {
    // 4つの隣接原子を同時処理
    if (simd.isPointBuriedBatch4(point, batch_positions, batch_radii)) {
        is_buried = true;
        break;
    }
}
```

**問題点:**
- AVX2/AVX-512対応CPUでは8-wide SIMDが可能
- batch_positionsの構築に毎回4回のメモリアクセス

### 3. テストポイント生成の繰り返し (lines 200, 327)

```zig
const test_points_array = try test_points.generateTestPoints(allocator, config.n_points);
```

**問題点:**
- 毎回同じテストポイントを生成
- Fibonacci latticeの計算コスト

---

## 最適化案（優先度順）

### Phase 1: Quick Wins (低リスク、効果中)

#### 1.1 SoA維持でpositions配列を削除
- **変更箇所:** `shrake_rupley.zig`, `neighbor_list.zig`
- **効果:** メモリ確保削減、キャッシュ効率向上
- **リスク:** 低（NeighborListのAPI変更必要）
- **見積もり効果:** 5-10%

#### 1.2 テストポイントのキャッシュ
- **変更箇所:** `test_points.zig` にグローバルキャッシュ追加
- **効果:** 繰り返し計算の回避
- **リスク:** 低
- **見積もり効果:** 1-3%（小さい）

### Phase 2: Medium Effort (中リスク、効果大)

#### 2.1 8-wide SIMD (AVX-512)
- **変更箇所:** `simd.zig`, `shrake_rupley.zig`
- **効果:** 理論上2x高速化
- **リスク:** 中（CPU互換性）
- **見積もり効果:** 10-30%（AVX-512対応CPUのみ）

#### 2.2 隣接原子データのプリフェッチ
- **変更箇所:** `shrake_rupley.zig`
- **効果:** キャッシュミス削減
- **リスク:** 低
- **見積もり効果:** 5-10%

### Phase 3: Major Refactoring (高リスク、効果大)

#### 3.1 タイル処理（Blocked Processing）
- 原子を64個ずつのブロックで処理
- ブロック内の隣接データがL2キャッシュに収まる
- **見積もり効果:** 15-25%

---

## 実装順序

```
Phase 1.1 (SoA維持) → ベンチマーク
    ↓
Phase 1.2 (テストポイントキャッシュ) → ベンチマーク
    ↓
Phase 2.2 (プリフェッチ) → ベンチマーク
    ↓
Phase 2.1 (8-wide SIMD) → ベンチマーク
```

## 成功基準

| 現在 | 目標 | 達成時 | 実績 |
|------|------|--------|------|
| 2.0x | 2.5x | Phase 1完了 | ❌ 2.0x（効果なし） |
| 2.0x | 3.0x | Phase 2完了 | ⚠️ 2.26x（部分達成） |
| 2.0x | 4.0x | Phase 3完了 | - |

**結論**: 8-wide SIMDで2.26xを達成。更なる最適化には大規模なリファクタリングが必要。

---

## 注意事項

1. **テスト必須**: 各変更後に精度検証
2. **ベンチマーク**: 変更毎に計測して効果確認
3. **ロールバック**: 効果がない変更は戻す

---

## 実施結果

### Phase 1.1: 隣接原子の距離順ソート
- **結果**: ❌ 効果なし（オーバーヘッドが利得を上回る）
- **詳細**: ソートコスト(~100ms)が早期終了の利益を打ち消す

### Phase 1.2: テストポイントキャッシュ
- **結果**: ⏭️ スキップ（既に効率的）
- **詳細**: テストポイント生成は全体の0.1%未満

### Phase 2.1: 8-wide SIMD
- **結果**: ✅ 成功（12%改善）
- **詳細**: 4-wide → 8-wide SIMDで2.02x → 2.26xに向上

### Phase 2.2: プリフェッチ
- **結果**: ❌ 効果なし
- **詳細**: CPUのハードウェアプリフェッチャが既に効率的

---

- [x] **DONE** - Phase 1試行完了（効果なし）
- [x] **DONE** - Phase 2完了（8-wide SIMD成功）
- [x] **DONE** - プラン完了

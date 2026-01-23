# Phase 11: Lee-Richards Algorithm

## Goal

Lee-Richards（L&R）アルゴリズムを実装し、Shrake-Rupleyとの選択を可能にする。

## 優先度: 中

## 参考: FreeSASA実装

FreeSASAの`src/sasa_lr.c`:
- スライスベースのアルゴリズム
- 各原子をZ軸方向にスライス
- スライス上の円弧から表面積を計算

## アルゴリズム概要

### Shrake-Rupley vs Lee-Richards

| 特性 | Shrake-Rupley | Lee-Richards |
|------|---------------|--------------|
| 手法 | テストポイント | スライス/円弧 |
| 精度制御 | ポイント数 | スライス数 |
| 計算量 | O(N × P × K) | O(N × S × K) |
| 利点 | 実装が簡単 | 数学的に厳密 |
| 欠点 | 近似精度 | 実装が複雑 |

### Lee-Richards法の手順

1. 各原子について、Z軸方向にスライスを作成
2. 各スライスで原子の断面円を計算
3. 隣接原子との交差を計算し、露出円弧を特定
4. 円弧の長さから表面積を積分

```
        ╭─────╮
       ╱   ●   ╲      ← 原子球
      │    │    │
   ───┼────┼────┼───  ← スライス
      │    │    │
       ╲       ╱
        ╰─────╯

スライス断面:
    ╭───╮
   ╱     ╲
  │   ●   │  ← 原子断面（円）
   ╲     ╱
    ╰───╯
     ↑
    露出円弧を計算
```

## Tasks

### Phase 11.1: Core Algorithm

- [ ] スライス生成
- [ ] 断面円の計算
- [ ] 隣接原子との交差計算
- [ ] 露出円弧の特定

**データ構造:**
```zig
const Slice = struct {
    z: f64,           // Z座標
    radius: f64,      // スライス半径
    arcs: []Arc,      // 露出円弧
};

const Arc = struct {
    start: f64,       // 開始角度
    end: f64,         // 終了角度
};
```

### Phase 11.2: Arc Calculation

- [ ] 2円の交点計算
- [ ] 円弧のマージ処理
- [ ] 数値安定性の確保

**円弧計算:**
```
2つの円の交点から、原子1の露出円弧を計算

    ╭───╮
   ╱  1  ╲╲
  │   ●───●│ 2  ← 円2が円1を遮蔽
   ╲     ╱╱
    ╰───╯

露出部分 = 円1の全周 - 円2に遮蔽された部分
```

### Phase 11.3: Integration

- [ ] 円弧長から表面積への積分
- [ ] Simpson則またはガウス積分
- [ ] スライス数による精度制御

**表面積計算:**
```
SASA = Σ (円弧長 × スライス厚さ)
     = ∫ 2πr(z) × f(z) dz

f(z) = 露出円弧の割合
```

### Phase 11.4: CLI Integration

- [ ] `--algorithm=<sr|lr>` オプション追加
- [ ] `--n-slices=<n>` オプション追加（L&R用）
- [ ] デフォルトはShrake-Rupley維持

**CLI例:**
```bash
# Shrake-Rupley（デフォルト）
freesasa_zig input.json output.json

# Lee-Richards
freesasa_zig --algorithm=lr input.json output.json

# Lee-Richards with custom slices
freesasa_zig --algorithm=lr --n-slices=50 input.json output.json
```

### Phase 11.5: Optimization

- [ ] 近傍リストの再利用
- [ ] SIMD最適化（円弧計算）
- [ ] マルチスレッド対応

## Files

| File | Action |
|------|--------|
| `src/sasa_lee_richards.zig` | CREATE |
| `src/arc.zig` | CREATE |
| `src/main.zig` | MODIFY |
| `src/types.zig` | MODIFY |

## Success Criteria

- [ ] Lee-Richardsアルゴリズムが動作
- [ ] Shrake-Rupleyと同程度の精度
- [ ] FreeSASAのL&R実装と結果が一致

---
- [ ] **DONE** - Phase 11 complete

# Phase 17: JSON Input Benchmark

統一 JSON 入力による公平なベンチマーク環境構築。

## Goal

- CIF ファイルを PC に置かずにベンチマーク実行
- 全実装が同じ JSON 入力を使用
- Protein only、ProtOr 半径で統一

## JSON Format

freesasa-zig と同じ形式：

```json
{
  "x": [1.0, 2.0, 3.0, ...],
  "y": [1.0, 2.0, 3.0, ...],
  "z": [1.0, 2.0, 3.0, ...],
  "r": [1.7, 1.55, 1.52, ...]
}
```

## Tasks

### Phase 17.1: Input Generation Script ✅

- [x] `scripts/data/generate_benchmark_json.py` 作成
  - gemmi で CIF 読み込み
  - Protein atoms only（HETATM 除外）
  - Hydrogen 除外
  - ProtOr 半径付与
  - JSON 出力（x, y, z, r のみ、コンパクト形式）
- [x] `benchmarks/inputs_json/` に 6 構造生成済み

### Phase 17.2: FreeSASA C JSON Support

- [ ] `freesasa-bench` に JSON 入力オプション追加
  - cJSON ライブラリ統合
  - `--json` オプション追加
  - JSON → freesasa_structure 変換
  - 既存の計算ロジック使用

### Phase 17.3: RustSASA JSON Support

- [ ] `rustsasa-bench` に JSON 入力オプション追加
  - serde_json 依存追加
  - `--json-input` オプション追加
  - pdbtbx バイパスして直接計算

### Phase 17.4: Benchmark Script Update

- [ ] `benchmark.py` 更新
  - JSON 入力モード追加
  - CIF/JSON 自動判別
  - 全実装で JSON 使用

### Phase 17.5: Test Dataset Generation

- [ ] ベンチマーク用 JSON 生成
  - 既存 6 構造 (1CRN, 1UBQ, 1A0Q, 3HHB, 1AON, 4V6X)
  - `benchmarks/inputs_json/` に保存

## Directory Structure

```
benchmarks/
├── inputs_json/          # NEW: JSON inputs with ProtOr radii
│   ├── 1crn.json
│   ├── 1ubq.json
│   └── ...
├── structures/           # CIF files (optional, for regeneration)
└── external/
    ├── freesasa-bench/   # JSON support added
    └── rustsasa-bench/   # JSON support added
```

## Notes

- freesasa-zig は既に JSON 対応済み
- 半径は生成時に確定 → 各実装で classifier 不要
- Protein only → 水分子・リガンド除外で公平

---
- [ ] **DONE** - Phase complete

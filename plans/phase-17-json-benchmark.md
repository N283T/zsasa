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

- [x] `scripts/data/cif_to_protor_json.py` 作成（旧 generate_benchmark_json.py）
  - gemmi で CIF 読み込み
  - Protein atoms only（HETATM 除外）
  - Hydrogen 除外
  - ProtOr 半径付与
  - JSON 出力（x, y, z, r のみ、コンパクト形式、gzip 圧縮）
- [x] `benchmarks/inputs_json/` に 6 構造生成済み

### Phase 17.2: FreeSASA C JSON Support ✅

- [x] `freesasa-bench` に JSON 入力オプション追加
  - 独自 JSON パーサー実装（zlib で gzip 対応）
  - `--json-input` オプション追加
  - JSON → freesasa_structure 変換
  - 既存の計算ロジック使用

### Phase 17.3: RustSASA JSON Support ✅

- [x] `rustsasa-bench` に JSON 入力オプション追加
  - serde_json + flate2 依存追加
  - `-J/--json-input` オプション追加
  - pdbtbx バイパスして直接計算

### Phase 17.4: Benchmark Script ✅

- [x] 新規 `scripts/benchmark.py` 作成
  - JSON 入力のみ（公平な計測）
  - Zig/FreeSASA C/RustSASA 全対応
  - SASA 計算時間のみ計測
- [x] 旧 benchmark.py → `scripts/quick_compare.py` に改名（CIF 入力用）

### Phase 17.5: Test Dataset Generation ✅

- [x] ベンチマーク用 JSON 生成済み
  - 既存 6 構造 (1CRN, 1UBQ, 1A0Q, 3HHB, 1AON, 4V6X)
  - `benchmarks/inputs_json/` に .json.gz 形式で保存

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

## Results

初回ベンチマーク結果（M4 Mac, 10コア）:

| 構造 | 原子数 | Zig SR | RustSASA | FS-C | Zig vs FS-C | Zig vs Rust |
|------|--------|--------|----------|------|-------------|-------------|
| 1CRN | 327 | 0.43ms | 0.65ms | 0.85ms | 1.98x | 1.52x |
| 4V6X | 237,685 | 159ms | 501ms | 759ms | 4.77x | 3.16x |

---
- [x] **DONE** - Phase complete

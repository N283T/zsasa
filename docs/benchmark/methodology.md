# Benchmark Methodology / ベンチマーク手法

ベンチマークの測定方法と実行手順。結果は [results.md](results.md) を参照。

## SASA-Only Timing

公平な比較のため、**SASA計算時間のみ**を測定。ファイルI/Oは除外。

```
Total time = File I/O + SASA calculation + Output
                        ^^^^^^^^^^^^^^^^
                        これだけ測定
```

各実装での測定方法：

| Implementation | Method |
|----------------|--------|
| freesasa-zig | `--timing` オプションで内部計測（stderr出力） |
| FreeSASA C | パッチ版バイナリが SASA 計算時間を stderr に出力 |
| RustSASA | パッチ版バイナリが `SASA_TIME_US` を stderr に出力 |

## Parameters / パラメータ

| Parameter | Value | Notes |
|-----------|-------|-------|
| Algorithm | Shrake-Rupley | 全実装で対応（LRはZig/FreeSASAのみ） |
| n_points | 100 | テストポイント数 |
| probe_radius | 1.4 Å | 水分子半径 |
| Runs | 3 | 平均値を使用 |

## Stratified Sampling / 層化サンプリング

PDB 全構造（約 238K）から層化サンプリング：

| Bin | Range | Strategy |
|-----|-------|----------|
| 0-500 | 0 〜 500 | 比例配分 |
| 500-2k | 500 〜 2,000 | 比例配分 |
| 2k-10k | 2,000 〜 10,000 | 比例配分 |
| 10k-50k | 10,000 〜 50,000 | 比例配分 |
| 50k-200k | 50,000 〜 200,000 | **全件採用** |
| 200k+ | 200,000 〜 | **全件採用** |

希少な大規模構造（50k 以上）は全件採用し、残りを比例配分。

## Running Benchmarks / ベンチマーク実行

### Setup / セットアップ

```bash
# Zig バイナリをビルド
zig build -Doptimize=ReleaseFast

# 外部ツールのセットアップ (比較用)
cd benchmarks/external
git clone https://github.com/N283T/freesasa-bench.git
cd freesasa-bench && ./configure --enable-threads && make && cd ..
git clone https://github.com/N283T/rustsasa-bench.git
cd rustsasa-bench && cargo build --release --features cli && cd ..
```

### Index & Sample / インデックス＆サンプル生成

```bash
# インデックス作成（初回のみ）
./benchmarks/scripts/build_index.py benchmarks/inputs

# 分布確認
./benchmarks/scripts/sample.py benchmarks/inputs/index.json --analyze

# サンプル生成
./benchmarks/scripts/sample.py benchmarks/inputs/index.json \
    --target 100000 --seed 42 \
    -o benchmarks/samples/stratified_100k.json
```

### Single-File Mode / 個別ファイルモード

各ファイルを順番に処理し、アルゴリズム内並列化をテスト：

```bash
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/run.py --tool freesasa --algorithm lr --threads 1-10

# サンプルファイル指定
./benchmarks/scripts/run.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_100k.json \
    --threads 1,2,4,8,10
```

### Batch Mode / バッチモード

複数ファイルを並列処理し、スループットをテスト：

```bash
./benchmarks/scripts/run_batch.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_1k.json \
    --threads 1,4,8,10 --runs 3
```

| Tool | Batch Implementation |
|------|---------------------|
| Zig | Native directory input |
| Rust | `--json-dir` + rayon |
| FreeSASA | Shell script + background jobs |

### Analysis / 分析

```bash
./benchmarks/scripts/analyze.py summary    # サマリーテーブル
./benchmarks/scripts/analyze.py validate   # SASA検証
./benchmarks/scripts/analyze.py plot       # グラフ生成
./benchmarks/scripts/analyze.py all        # 全部

# CSVエクスポート
./benchmarks/scripts/analyze.py export_csv
```

## Scripts / スクリプト一覧

| Script | Purpose |
|--------|---------|
| `build_index.py` | 全ファイルから原子数インデックスを作成 |
| `sample.py` | インデックスから層化サンプリング |
| `run.py` | 個別ファイルモードでベンチマーク実行 |
| `run_batch.py` | バッチモードでベンチマーク実行 |
| `analyze.py` | 結果分析・グラフ生成 |
| `generate_json.py` | CIF→JSON変換 |

## Reproducibility / 再現性

同一シードで同一サンプル：

```bash
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o a.json
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o b.json
diff a.json b.json  # 差分なし
```

ベンチマーク実行時に表示されるシステム情報：

```
=== System Info ===
Platform: Darwin (arm64)
CPU: Apple M4
Cores: 10
Memory: 24.0 GB

=== Benchmark Config ===
Threads: 4
N-points: 100
Probe radius: 1.4 Å
Runs: 5
```

## Notes / 注意事項

1. **Initial run は遅い**: ファイルキャッシュ、warmup の影響
2. **Thread 数は CPU に依存**: 物理コア数に合わせるのが最適
3. **外部ツールは要パッチ**: SASA-only timing には改造版バイナリが必要

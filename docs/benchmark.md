# Benchmark / ベンチマーク

freesasa-zig のベンチマーク手法と結果。

## Methodology / 測定手法

### SASA-Only Timing

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

### Benchmark Dataset

6構造、小規模から超大規模まで：

| PDB ID | Atoms | Category | Description |
|--------|------:|----------|-------------|
| 1CRN | 327 | Tiny | Crambin |
| 1UBQ | 602 | Small | Ubiquitin |
| 1A0Q | 3,183 | Medium | Lipid transfer protein |
| 3HHB | 4,384 | Medium | Hemoglobin tetramer |
| 1AON | 58,674 | Large | GroEL-GroES complex |
| 4V6X | 237,685 | XLarge | Ribosome |

### Parameters

デフォルトパラメータ：

| Parameter | Value | Notes |
|-----------|-------|-------|
| Algorithm | Shrake-Rupley | 全実装で対応 |
| n_points | 100 | テストポイント数 |
| probe_radius | 1.4 Å | 水分子半径 |
| Threads | 4 | マルチスレッド比較時 |

### Runs

- Default: 3 runs (average)
- Configurable: `--runs N`

## Results / 結果

### vs FreeSASA C (4 threads)

Shrake-Rupley:

| Structure | Atoms | Zig (ms) | FS-C (ms) | Speedup |
|-----------|------:|--------:|---------:|--------:|
| 1CRN | 327 | 0.44 | 0.53 | 1.20x |
| 1UBQ | 602 | 0.63 | 0.88 | 1.40x |
| 1A0Q | 3,183 | 2.64 | 4.42 | 1.67x |
| 3HHB | 4,384 | 3.54 | 5.95 | 1.68x |
| 1AON | 58,674 | 45.61 | 98.28 | 2.15x |
| 4V6X | 237,685 | 189.06 | 424.53 | 2.25x |

Lee-Richards:

| Structure | Atoms | Zig (ms) | FS-C (ms) | Speedup |
|-----------|------:|--------:|---------:|--------:|
| 1CRN | 327 | 1.45 | 1.64 | 1.13x |
| 1UBQ | 602 | 2.19 | 2.99 | 1.37x |
| 1A0Q | 3,183 | 9.90 | 16.09 | 1.62x |
| 3HHB | 4,384 | 14.15 | 23.44 | 1.66x |
| 1AON | 58,674 | 182.83 | 317.41 | 1.74x |
| 4V6X | 237,685 | 741.86 | 1293.48 | 1.74x |

**Summary**: SR は **1.2x-2.3x 高速**、LR は **1.1x-1.7x 高速**。構造サイズが大きいほど差が開く。

### vs RustSASA

RustSASA は Shrake-Rupley のみ対応。

#### RustSASA の「5x 高速」について

RustSASA は「FreeSASA より 5x 高速」と謳っているが、これは**バッチ処理**での話：

| Benchmark | RustSASA | FreeSASA | Speedup |
|-----------|----------|----------|---------|
| Single protein | 4.0 ms | 4.0 ms | 1.0x |
| E. coli proteome (4,391 files) | 5.2 s | 28.0 s | 5.4x |

「5x」は RustSASA の CLI が GNU Parallel 相当の並列処理をしているため。
**アルゴリズム自体の速度は FreeSASA と同等**。

#### SASA-only comparison (single-threaded)

| Structure | Atoms | Zig (ms) | RustSASA (ms) | FreeSASA C (ms) | vs Rust | vs FS-C |
|-----------|------:|---------:|--------------:|----------------:|--------:|--------:|
| 1CRN | 327 | 0.40 | 0.69 | 0.79 | 1.7x | 2.0x |
| 1UBQ | 602 | 0.53 | 1.23 | 1.44 | 2.3x | 2.8x |
| 4V6X | 237,685 | 158.58 | 500.11 | 747.95 | 3.2x | 4.7x |

**Summary**: freesasa-zig は**アルゴリズムレベルで高速**。
- vs RustSASA: **1.7x-3.2x 高速**
- vs FreeSASA C: **2.0x-4.7x 高速**

### Batch Processing Benchmark / バッチ処理ベンチマーク

4,304 構造（層化サンプリング）を 10 スレッドで並列処理：

| Tool | Total Time | Files/sec | Speedup |
|------|------------|-----------|---------|
| **Zig** | 162.5 s | 26.5 | **8.3x** |
| **Rust** | 162.7 s | 26.4 | **8.3x** |
| FreeSASA | 1,355 s | 3.2 | 1.0x |

バッチモードでは各ファイルのSASA計算を1スレッドで実行し、ファイル処理を並列化。
Zig と Rust は同等の性能で、FreeSASA より約 **8x 高速**。

## Why is Zig Faster? / なぜ Zig が速いか

詳細は [optimizations.md](optimizations.md) を参照。

1. **Neighbor List**: O(N²) → O(N) の近傍探索
2. **8-wide SIMD**: `@Vector(8, f64)` による並列計算
3. **Fast Trigonometry**: LR 向け多項式近似（~37% 高速化）
4. **Work-stealing Thread Pool**: 効率的な並列処理

## Large-Scale Benchmark / 大規模ベンチマーク

238K 構造の全 PDB データセットから層化サンプリングでベンチマークを実行する方法。

### Stratified Sampling / 層化サンプリング

原子数に基づく対数スケール 6 ビンでサンプリング：

| Bin | Range | Strategy |
|-----|-------|----------|
| 0-500 | 0 〜 500 | 比例配分 |
| 500-2k | 500 〜 2,000 | 比例配分 |
| 2k-10k | 2,000 〜 10,000 | 比例配分 |
| 10k-50k | 10,000 〜 50,000 | 比例配分 |
| 50k-200k | 50,000 〜 200,000 | **全件採用** |
| 200k+ | 200,000 〜 | **全件採用** |

希少な大規模構造（50k 以上）は全件採用し、残りのターゲット数を他のビンから比例配分。

### Workflow / ワークフロー

```bash
# 1. インデックス作成（初回のみ、238K 件で約 5 分）
./benchmarks/scripts/build_index.py benchmarks/inputs
# → benchmarks/inputs/index.json

# 2. 分布確認
./benchmarks/scripts/sample.py benchmarks/inputs/index.json --analyze

# 3. サンプル生成（75,000 件、シード 42）
./benchmarks/scripts/sample.py benchmarks/inputs/index.json \
    --target 75000 --seed 42 \
    --output benchmarks/samples/stratified_75k.json

# 4. ベンチマーク実行
./benchmarks/scripts/run.py \
    --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_75k.json \
    --threads 1-10
```

### Scripts / スクリプト

| Script | Purpose |
|--------|---------|
| `build_index.py` | 全ファイルから原子数インデックスを作成 |
| `sample.py` | インデックスから層化サンプリング |
| `run.py --sample-file` | サンプルリストでフィルタして実行 |

### Reproducibility / 再現性

同一シードで同一結果が得られる：

```bash
# 同じシードなら同じサンプル
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o a.json
./benchmarks/scripts/sample.py index.json --target 75000 --seed 42 -o b.json
diff a.json b.json  # 差分なし
```

## Running Benchmarks / ベンチマーク実行

### Single-File Mode / 個別ファイルモード

各ファイルを順番に処理し、アルゴリズム内並列化をテスト：

```bash
# Run benchmark (single tool, single algorithm)
./benchmarks/scripts/run.py --tool zig --algorithm sr --threads 1-10
./benchmarks/scripts/run.py --tool freesasa --algorithm lr --threads 1-10

# With custom input directory
./benchmarks/scripts/run.py --tool zig --algorithm sr --input-dir benchmarks/inputs

# With stratified sample
./benchmarks/scripts/run.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_75k.json
```

### Batch Mode / バッチモード

複数ファイルを並列処理し、スループットをテスト：

```bash
# Batch benchmark (file-level parallelism)
./benchmarks/scripts/run_batch.py --tool zig --algorithm sr \
    --input-dir benchmarks/inputs \
    --sample-file benchmarks/samples/stratified_1k.json \
    --threads 1,4,8,10 --runs 3

# Each tool uses its native batch mode:
# - Zig: directory input with --timing
# - Rust: --json-dir with rayon
# - FreeSASA: freesasa_batch.sh with background jobs
```

バッチモードでは各ファイルのSASA計算は1スレッドで、ファイル処理を並列化：

| Tool | Batch Implementation |
|------|---------------------|
| Zig | Native directory input |
| Rust | `--json-dir` + rayon |
| FreeSASA | Shell script + background jobs |

### Analysis / 分析

```bash
# Analyze results
./benchmarks/scripts/analyze.py summary    # Summary tables
./benchmarks/scripts/analyze.py validate   # SASA validation
./benchmarks/scripts/analyze.py plot       # Generate graphs
./benchmarks/scripts/analyze.py all        # All of the above

# Generate JSON inputs from CIF files
./benchmarks/scripts/generate_json.py /path/to/cif benchmarks/dataset
./benchmarks/scripts/generate_json.py /path/to/cif out.tar.gz --archive  # Archive mode
```

### Requirements

- freesasa-zig binary (`zig build -Doptimize=ReleaseFast`)
- Python 3.11+ with dependencies (auto-installed via PEP 723)

### External Dependencies (Optional)

比較ベンチマーク用に、タイミングパッチ済みの fork を使用：

```bash
cd benchmarks/external

# FreeSASA C
git clone https://github.com/N283T/freesasa-bench.git
cd freesasa-bench && ./configure --enable-threads && make && cd ..

# RustSASA
git clone https://github.com/N283T/rustsasa-bench.git
cd rustsasa-bench && cargo build --release --features cli && cd ..
```

詳細は [benchmarks/external/README.md](../benchmarks/external/README.md) を参照。

## Reproducibility / 再現性

ベンチマーク実行時に表示される情報：

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

これらの情報を記録して結果の再現性を確保。

## Notes / 注意事項

1. **Initial run は遅い**: ファイルキャッシュ、JIT warmup などの影響
2. **Thread 数は CPU に依存**: 物理コア数に合わせるのが最適
3. **RustSASA は要パッチ**: SASA-only timing には改造版バイナリが必要

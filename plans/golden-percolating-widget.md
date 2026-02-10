# Phase 1: Refactor analyze.py

## Goal

`benchmarks/scripts/analyze.py` (2200行) を3ファイルに分割し、バグを修正する。

## 問題点

1. **バグ**: `speedup()` の `typer.Option()` デフォルト値が `all()` から呼び出し時に解決されない
2. **巨大ファイル**: 2200行の単一ファイル
3. **責務混在**: データ読み込み、プロット、CLI、テーブル生成が混在
4. **巨大関数**: `_collect_template_data()` が600行超

## 分割構成

```
benchmarks/scripts/
├── analyze.py        # CLI エントリポイント (コマンド定義)
├── analyze_plots.py  # プロット関数
└── analyze_data.py   # データ読み込み・ヘルパー・定数
```

### analyze_data.py (〜200行)
- 定数: `COLORS`, `LINESTYLES`, `BINS`, `RESULTS_DIR`, `PLOTS_DIR`
- `load_config()` - YAML設定読み込み
- `load_data()` - ベンチマークCSV読み込み
- `add_size_bin()` - サイズビン追加
- `compute_speedup_by_bin()` - ビン別スピードアップ計算
- `setup_style()` - matplotlib設定

### analyze_plots.py (〜800行)
- `_plot_scatter()` - 散布図描画
- `_plot_threads()` - スレッド図描画
- `_plot_speedup_single()` - スピードアップ図描画
- `plot_scatter()` - 散布図生成 (コマンド実装)
- `plot_threads()` - スレッド図生成
- `plot_grid()` - グリッド図生成
- `plot_validation()` - 検証図生成
- `plot_samples()` - サンプル図生成
- `plot_large()` - 大規模構造図生成
- `plot_efficiency()` - 効率図生成
- `plot_speedup()` - スピードアップ図生成

### analyze.py (〜600行)
- Typerアプリ定義
- コマンド登録 (プロット関数をラップ)
- `summary()` - サマリーテーブル
- `export_csv()` - CSV出力
- `render_docs()` - ドキュメント生成
- `all()` - 全実行
- `_generate_markdown_table()` - Markdownテーブル生成
- `_collect_template_data()` - テンプレートデータ収集

## バグ修正

```python
# Before (Line 2195)
speedup()  # typer.Option() が解決されない

# After
speedup(min_atoms=50000, top_n=5)  # デフォルト値を明示
```

## Tasks

1. `analyze_data.py` を作成
2. `analyze_plots.py` を作成
3. `analyze.py` を書き換え (インポート追加、コマンドラッパー)
4. バグ修正 (`all()` での `speedup()` 呼び出し)
5. 動作確認

## 検証

```bash
# 個別コマンド
./benchmarks/scripts/analyze.py summary
./benchmarks/scripts/analyze.py scatter
./benchmarks/scripts/analyze.py threads

# 全実行 (バグ修正確認)
./benchmarks/scripts/analyze.py all
```

## Critical Files

- `/Users/nagaet/freesasa-zig/benchmarks/scripts/analyze.py` (書き換え)
- `/Users/nagaet/freesasa-zig/benchmarks/scripts/analyze_data.py` (新規)
- `/Users/nagaet/freesasa-zig/benchmarks/scripts/analyze_plots.py` (新規)

---
- [x] **DONE** - Phase complete

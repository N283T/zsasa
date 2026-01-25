# Phase 17.4: Benchmark Script Redesign

## 現状の問題

現在の `benchmark.py` は:
- 結果をコンソールに表示するだけ
- 記録が残らない
- SASA値の検証がない
- グラフ生成機能がない

## 要件

### 1. 結果の記録

**実行条件 (JSON)**
```json
{
  "timestamp": "2025-01-25T12:00:00",
  "system": {
    "cpu": "Apple M4",
    "cores": 10,
    "memory": "32 GB",
    "os": "Darwin 24.6.0"
  },
  "parameters": {
    "n_runs": 3,
    "n_threads": 10,
    "n_points": 100,
    "n_slices": 20,
    "probe_radius": 1.4
  }
}
```

**結果 (CSV)**
```csv
pdb_id,n_atoms,algorithm,tool,run,sasa_time_ms,total_sasa
1crn,327,sr,zig,1,0.42,1234.56
1crn,327,sr,zig,2,0.41,1234.56
1crn,327,sr,freesasa-c,1,0.85,1234.78
...
```

### 2. SASA値の検証

- 各ツールの total_sasa を比較
- 許容誤差: 0.1% (同一入力なので高精度を期待)
- 不一致があれば警告表示
- 検証結果もCSVに記録

### 3. グラフ生成

**benchmark.py 内 (簡易)**
- 実行後に簡易グラフをターミナルに表示（rich使用）
- または PNG 1枚生成

**別スクリプト (詳細)**
- `scripts/plot_benchmark.py`
- CSV読み込み → 複数グラフ生成
  - 構造サイズ vs 実行時間 (log-log)
  - ツール間比較棒グラフ
  - スケーリング曲線
  - SASA値の一致確認プロット

### 4. 出力先

```
benchmarks/
├── results/
│   ├── 2025-01-25_120000_config.json
│   ├── 2025-01-25_120000_results.csv
│   └── 2025-01-25_120000_summary.png  (optional)
└── inputs_json/
    └── *.json.gz
```

## タスク

### Phase 17.4.1: benchmark.py リファクタリング ✅
- [x] 結果を CSV に記録
- [x] 実行条件を JSON に記録
- [x] 各ツールから total_sasa を取得
- [x] SASA 値の検証ロジック追加
- [x] 簡易サマリー表示（既存の rich テーブル改善）
- [x] 分離実行（--tool オプション）
- [x] スレッド範囲指定（--threads 1-10）

### Phase 17.4.2: グラフ生成スクリプト
- [ ] `scripts/plot_benchmark.py` 作成
- [ ] CSV 読み込み
- [ ] matplotlib で複数グラフ生成
  - 実行時間比較
  - スケーリング曲線
  - SASA 値比較

### Phase 17.4.3: テストと検証
- [ ] 全構造でベンチマーク実行
- [ ] SASA 値の一致確認
- [ ] グラフ生成確認

## 設計メモ

### SASA値取得方法

| ツール | 取得方法 |
|--------|----------|
| Zig | 出力JSON の `total_area` |
| FreeSASA C | stdout の `Total: XXX` |
| RustSASA | stdout の `Total SASA: XXX` |

### 許容誤差の根拠

- 同一座標・同一半径入力
- 同一アルゴリズム (SR 100点, LR 20スライス)
- 浮動小数点誤差 + 実装差異を考慮して 0.1%

---
- [ ] **DONE** - Phase complete

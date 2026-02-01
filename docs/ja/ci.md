# CI/CD設定

## 概要

GitHub Actionsによる継続的インテグレーション。

## トリガー

```yaml
on:
  push:
    branches: [main]
    paths-ignore:
      - '*.md'
      - 'docs/**'
      - 'plans/**'
      - 'benchmarks/**'
      - 'LICENSE'
  pull_request:
    branches: [main]
    paths-ignore: [同上]
  workflow_dispatch:  # 手動トリガー
```

### 除外パス（CIスキップ）

以下のファイルの変更ではCIは実行されません:

| パス | 理由 |
|------|--------|
| `*.md` | ドキュメントのみ |
| `docs/**` | ドキュメントのみ |
| `plans/**` | プランファイルのみ |
| `benchmarks/**` | ベンチマークデータとスクリプト |
| `LICENSE` | ライセンスファイル |

### 手動トリガー

`workflow_dispatch`によりGitHub UIからCIを手動実行可能。

Actions → CI → Run workflow

## ジョブ構成

```
┌──────────────┐  ┌─────────────────────────────────┐  ┌─────────────────────────────────┐
│ fmt          │  │ build (並列 x2)                 │  │ python (並列 x2)                │
│ (5min)       │  │ Linux / macOS                   │  │ Linux / macOS                   │
│              │  │ (15min)                         │  │ (15min)                         │
└──────┬───────┘  └────────────────┬────────────────┘  └────────────────┬────────────────┘
       │                           │                                    │
       └───────────┬───────────────┘                                    │
                   ▼                                                    │
         ┌─────────────────┐                                            │
         │ validate        │ ◄──────────────────────────────────────────┘
         │ (10min)         │   (pythonはfmtのみに依存)
         └─────────────────┘
```

## ジョブ詳細

### 1. フォーマットチェック (`fmt`)

**環境:** ubuntu-latest
**タイムアウト:** 5分
**実行時間:** 約20秒

```bash
zig fmt --check src/
```

`src/`ディレクトリ配下のすべての`.zig`ファイルを再帰的にチェック。

フォーマット違反があれば失敗。修正方法:

```bash
zig fmt src/
```

### 2. ビルド (`build`)

**環境:** Linux, macOS（並列実行）
**タイムアウト:** 15分
**実行時間:**
- Linux: 約45秒
- macOS: 約1分10秒

| ステップ | コマンド | 目的 |
|------|---------|---------|
| ビルド | `zig build` | デバッグビルド |
| テスト | `zig build test` | 全テスト実行 |
| リリース | `zig build -Doptimize=ReleaseFast` | 最適化ビルド |
| バイナリ確認 | - | バイナリ存在確認 |
| CLIテスト | `--help`, `--version` | CLI機能確認 |

**注:** WindowsはCIマトリクスから除外。WindowsユーザーはWSLを使用してください。

### 3. 検証 (`validate`)

**環境:** ubuntu-latest
**タイムアウト:** 10分
**依存:** fmt, build（両方成功後に実行）
**実行時間:** 約45秒

| ステップ | 説明 |
|------|-------------|
| サンプル実行 | `benchmarks/dataset/1a0q.json.gz`で計算 |
| JSON検証 | 出力フィールドを検証（Python） |
| CSVテスト | CSV出力をテスト |
| 検証モード | `--validate`フラグをテスト |

**JSON検証:**
- `total_area`フィールドが存在
- `atom_areas`フィールドが存在
- 原子数 == 3183（1A0Qの原子数）

### 4. Python (`python`)

**環境:** Linux, macOS（並列実行）
**タイムアウト:** 15分
**依存:** fmt
**実行時間:** 約1-2分

| ステップ | 説明 |
|------|-------------|
| 共有ライブラリビルド | `zig build -Doptimize=ReleaseFast` |
| ライブラリ確認 | 共有ライブラリ存在確認 |
| パッケージインストール | `uv pip install -e ".[dev]"` |
| テスト実行 | `pytest tests/ -v` |
| リント | `ruff check` / `ruff format --check` |
| 統合テスト | Pythonから基本機能を確認 |

**共有ライブラリ:**
- Linux: `libzsasa.so`
- macOS: `libzsasa.dylib`

**Pythonテストカバレッジ:**
- NumPy配列入出力
- SR/LRアルゴリズム
- エラーハンドリング
- パラメータ設定

## ローカルでのCI再現

```bash
# フォーマットチェック
zig fmt --check src/

# ビルドとテスト
zig build
zig build test
zig build -Doptimize=ReleaseFast

# CLIテスト
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa --version

# 検証テスト
./zig-out/bin/zsasa benchmarks/dataset/1a0q.json.gz /tmp/output.json
./zig-out/bin/zsasa --format=csv benchmarks/dataset/1a0q.json.gz /tmp/output.csv
./zig-out/bin/zsasa --validate benchmarks/dataset/1a0q.json.gz

# Pythonテスト
cd python
uv pip install -e ".[dev]"
pytest tests/ -v
ruff check zsasa/
ruff format --check zsasa/
```

## CI失敗時の対処

### フォーマットチェック失敗

```bash
# ローカルでフォーマット実行
zig fmt src/
git add -u && git commit -m "style: Format code"
git push
```

**注:** プッシュ済みのコミットには`--amend` + `git push -f`を使用しないでください。履歴が書き換わります。共有ブランチでは避けてください。

### ビルド失敗

1. エラーメッセージを確認
2. `zig build`をローカルで実行して再現
3. 修正してプッシュ

### テスト失敗

1. `zig build test`をローカルで実行
2. 失敗しているテスト名を特定
3. 該当テストを修正

### 検証失敗

1. `benchmarks/dataset/1a0q.json.gz`での実行を確認
2. 出力JSONの形式を確認
3. 原子数が3183であることを確認

### Pythonテスト失敗

1. 共有ライブラリがビルドされていることを確認
   ```bash
   ls -la zig-out/lib/
   ```
2. Pythonパッケージをインストール
   ```bash
   cd python && uv pip install -e ".[dev]"
   ```
3. テストを実行
   ```bash
   uv run --with pytest pytest tests/ -v
   ```

### Pythonリント失敗

```bash
cd python
ruff check zsasa/ --fix
ruff format zsasa/
git add -u && git commit -m "style: Fix Python lint issues"
git push
```

### タイムアウト

デフォルトのタイムアウト:
- fmt: 5分
- build: 15分
- validate: 10分
- python: 15分

通常これらを超えることはありません。超えた場合は無限ループを調査してください。

## Zigバージョン

CIで使用されるZigバージョン: **0.15.2**

`goto-bus-stop/setup-zig@v2`アクションで管理。

```yaml
- uses: goto-bus-stop/setup-zig@v2
  with:
    version: 0.15.2
```

## セキュリティ

- サードパーティアクション:
  - `actions/checkout@v4`
  - `goto-bus-stop/setup-zig@v2`
  - `actions/setup-python@v5`
  - `astral-sh/setup-uv@v4`
- ハードコードされたシークレット: なし
- 危険な権限: なし（デフォルト権限を使用）

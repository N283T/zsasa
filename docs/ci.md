# CI/CD 構成

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
      - 'examples/**'
      - 'scripts/**'
      - 'LICENSE'
  pull_request:
    branches: [main]
    paths-ignore: [同上]
  workflow_dispatch:  # 手動トリガー
```

### 除外パス（CIをスキップ）

以下のファイル変更ではCIは実行されない:

| パス | 理由 |
|------|------|
| `*.md` | ドキュメントのみ |
| `docs/**` | ドキュメントのみ |
| `plans/**` | 計画ファイルのみ |
| `examples/**` | サンプルデータ（テストフィクスチャ） |
| `scripts/**` | Pythonスクリプト（Zig無関係） |
| `LICENSE` | ライセンスファイル |

**注意:** `examples/input_1a0q.json`はvalidateジョブで使用されるテストフィクスチャです。このファイルを変更してもCIは実行されないため、変更時は手動でテストを実行してください。

### 手動トリガー

`workflow_dispatch`により、GitHub UIから手動でCIを実行可能。

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
         │ (10min)         │   (python は fmt のみに依存)
         └─────────────────┘
```

## ジョブ詳細

### 1. Format Check (`fmt`)

**実行環境:** ubuntu-latest
**タイムアウト:** 5分
**実行時間:** 約20秒

```bash
zig fmt --check src/
```

`src/`ディレクトリ以下の全`.zig`ファイルを再帰的にチェック。

フォーマット違反があると失敗。修正方法:

```bash
zig fmt src/
```

### 2. Build (`build`)

**実行環境:** Linux, macOS（並列実行）
**タイムアウト:** 15分
**実行時間:**
- Linux: 約45秒
- macOS: 約1分10秒

| ステップ | コマンド | 目的 |
|----------|----------|------|
| Build | `zig build` | Debug build |
| Test | `zig build test` | 全テスト実行 |
| Release | `zig build -Doptimize=ReleaseFast` | 最適化ビルド |
| Binary check | - | バイナリ存在確認 |
| CLI test | `--help`, `--version` | CLI動作確認 |

**注意:** WindowsはCIマトリックスから除外。WindowsユーザーはWSLを使用してください。

### 3. Validate (`validate`)

**実行環境:** ubuntu-latest
**タイムアウト:** 10分
**依存:** fmt, build（両方成功後に実行）
**実行時間:** 約45秒

| ステップ | 内容 |
|----------|------|
| Example run | `examples/input_1a0q.json`で計算 |
| JSON validation | 出力フィールド検証（Python） |
| CSV test | CSV出力テスト |
| Validate mode | `--validate`フラグテスト |

**JSON検証内容:**
- `total_area` フィールド存在
- `atom_areas` フィールド存在
- 原子数 == 3183（1A0Qの原子数）

### 4. Python (`python`)

**実行環境:** Linux, macOS（並列実行）
**タイムアウト:** 15分
**依存:** fmt
**実行時間:** 約1-2分

| ステップ | 内容 |
|----------|------|
| Build shared library | `zig build -Doptimize=ReleaseFast` |
| Library check | 共有ライブラリ存在確認 |
| Install package | `uv pip install -e ".[dev]"` |
| Run tests | `pytest tests/ -v` |
| Lint | `ruff check` / `ruff format --check` |
| Integration test | Pythonからの基本動作確認 |

**共有ライブラリ:**
- Linux: `libfreesasa_zig.so`
- macOS: `libfreesasa_zig.dylib`

**Pythonテスト内容:**
- NumPy配列入出力
- SR/LRアルゴリズム
- エラーハンドリング
- パラメータ設定

## ローカルでのCI再現

```bash
# フォーマットチェック
zig fmt --check src/

# ビルド・テスト
zig build
zig build test
zig build -Doptimize=ReleaseFast

# CLIテスト
./zig-out/bin/freesasa_zig --help
./zig-out/bin/freesasa_zig --version

# 検証テスト
./zig-out/bin/freesasa_zig examples/input_1a0q.json /tmp/output.json
./zig-out/bin/freesasa_zig --format=csv examples/input_1a0q.json /tmp/output.csv
./zig-out/bin/freesasa_zig --validate examples/input_1a0q.json

# Pythonテスト
cd python
uv pip install -e ".[dev]"
uv run --with pytest pytest tests/ -v
uv run --with ruff ruff check freesasa_zig/
uv run --with ruff ruff format --check freesasa_zig/
```

## CI失敗時の対処

### Format check failed

```bash
# ローカルでフォーマット実行
zig fmt src/
git add -u && git commit --amend
git push -f
```

### Build failed

1. エラーメッセージを確認
2. ローカルで `zig build` を実行して再現
3. 修正後にプッシュ

### Test failed

1. `zig build test` をローカルで実行
2. 失敗したテスト名を確認
3. 該当テストを修正

### Validate failed

1. `examples/input_1a0q.json` で実行を確認
2. 出力JSONの形式を確認
3. 原子数が3183であることを確認

### Python tests failed

1. 共有ライブラリがビルドされているか確認
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

### Python lint failed

```bash
cd python
uv run --with ruff ruff check freesasa_zig/ --fix
uv run --with ruff ruff format freesasa_zig/
git add -u && git commit --amend
git push -f
```

### Timeout

デフォルトタイムアウト:
- fmt: 5分
- build: 15分
- validate: 10分
- python: 15分

通常これらを超えることはない。超えた場合は無限ループの可能性を調査。

## Zigバージョン

CI で使用しているZigバージョン: **0.15.2**

`goto-bus-stop/setup-zig@v2` アクションで管理。

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

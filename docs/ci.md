# CI/CD 構成

## 概要

GitHub Actionsによる継続的インテグレーション。

## トリガー

```yaml
on:
  push:
    branches: [main]
  pull_request:
    branches: [main]
```

### 除外パス（CIをスキップ）

以下のファイル変更ではCIは実行されない:

| パス | 理由 |
|------|------|
| `*.md` | ドキュメントのみ |
| `docs/**` | ドキュメントのみ |
| `plans/**` | 計画ファイルのみ |
| `examples/**` | サンプルデータのみ |
| `scripts/**` | Pythonスクリプト（Zig無関係） |
| `LICENSE` | ライセンスファイル |

## ジョブ構成

```
┌──────────────┐  ┌─────────────────────────────────┐
│ fmt          │  │ build (並列 x3)                 │
│ (Linux)      │  │ Linux / macOS / Windows         │
└──────┬───────┘  └────────────────┬────────────────┘
       │                           │
       └───────────┬───────────────┘
                   ▼
         ┌─────────────────┐
         │ validate        │
         │ (Linux)         │
         └─────────────────┘
```

## ジョブ詳細

### 1. Format Check (`fmt`)

**実行環境:** ubuntu-latest
**実行時間:** 約20秒

```bash
zig fmt --check src/*.zig
```

フォーマット違反があると失敗。修正方法:

```bash
zig fmt src/*.zig
```

### 2. Build (`build`)

**実行環境:** Linux, macOS, Windows（並列実行）
**実行時間:**
- Linux: 約45秒
- macOS: 約1分10秒
- Windows: 約2分

| ステップ | コマンド | 目的 |
|----------|----------|------|
| Build | `zig build` | Debug build |
| Test | `zig build test` | 全テスト実行 |
| Release | `zig build -Doptimize=ReleaseFast` | 最適化ビルド |
| Binary check | - | バイナリ存在確認 |
| CLI test | `--help`, `--version` | CLI動作確認 |

### 3. Validate (`validate`)

**実行環境:** ubuntu-latest
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

## ローカルでのCI再現

```bash
# フォーマットチェック
zig fmt --check src/*.zig

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
```

## CI失敗時の対処

### Format check failed

```bash
# ローカルでフォーマット実行
zig fmt src/*.zig
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

## Zigバージョン

CI で使用しているZigバージョン: **0.15.2**

`goto-bus-stop/setup-zig@v2` アクションで管理。

```yaml
- uses: goto-bus-stop/setup-zig@v2
  with:
    version: 0.15.2
```

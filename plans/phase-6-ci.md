# Phase 6: CI/CD Setup

## Goal

GitHub Actionsによる継続的インテグレーションを設定し、品質を自動的に担保する。

## 優先度: 高

## 参考: FreeSASA CI

FreeSASAのCI構成（`.github/workflows/`）:
- Matrix build: gcc, clang
- オプション有無のテスト
- `make check`によるテスト実行

## Tasks

### Phase 6.1: Basic CI

- [x] `.github/workflows/ci.yml` 作成
- [x] Zig build テスト
- [x] Zig test 実行
- [x] 複数プラットフォーム対応（ubuntu, macos, windows）
- [x] `zig fmt --check` フォーマットチェック
- [x] Validate job（サンプル実行・出力検証）

**CI設計:**
```yaml
name: CI
on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  build:
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: goto-bus-stop/setup-zig@v2
        with:
          version: 0.15.2
      - run: zig build
      - run: zig build test
```

### Phase 6.2: Release Optimization

- [x] ReleaseFast build テスト
- [ ] ReleaseSafe build テスト（スキップ - ReleaseFastで十分）
- [ ] バイナリサイズチェック（スキップ - 必要に応じて追加）

### Phase 6.3: Benchmark CI (Optional)

- [ ] ベンチマーク自動実行
- [ ] 性能回帰検出
- [ ] 結果のアーティファクト保存

## Files

| File | Action |
|------|--------|
| `.github/workflows/ci.yml` | CREATE |
| `.github/workflows/release.yml` | CREATE (optional) |

## Success Criteria

- [x] PRごとにCIが自動実行される
- [x] build/test失敗時にマージをブロック
- [x] macOS/Linux/Windowsで動作確認

---
- [x] **DONE** - Phase 6 complete

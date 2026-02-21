# Phase 1: Release Preparation

## Goal

v0.1.0 リリース前に必要な community governance files とリポジトリ整備を完了する。

## Tasks

### 1.1 Community Governance Files (All Parallel)

- [ ] `CODE_OF_CONDUCT.md` -- Contributor Covenant v2.1
- [ ] `.github/ISSUE_TEMPLATE/bug_report.yml` -- YAML form (Zig version, OS, repro steps)
- [ ] `.github/ISSUE_TEMPLATE/feature_request.yml` -- YAML form (use case, proposal)
- [ ] `.github/pull_request_template.md` -- Description, type, testing, checklist
- [ ] `CITATION.cff` -- CFF 1.2.0 format (DOI は Phase 2 Zenodo 後に更新)

### 1.2 Repository Polish (All Parallel)

- [ ] GitHub topics 設定: `sasa`, `bioinformatics`, `zig`, `python`, `molecular-dynamics`, `structural-biology`, `simd`
- [ ] GitHub Discussions 有効化
- [ ] Homepage URL 設定
- [ ] `python/pyproject.toml` 更新: author name/email, Documentation/Issues URLs 追加
- [ ] Stray files 確認・削除 (`output.json`, `echo`, `time` が repo root にある可能性)

### 1.3 Community Building (Ongoing)

- [ ] 5-10 個の Issue 作成 (planned improvements, future features)
- [ ] Milestone 作成 (v0.2.0, v1.0.0)
- [ ] "good first issue" ラベル追加

## Key Files

| File | Action |
|------|--------|
| `CODE_OF_CONDUCT.md` | Create |
| `.github/ISSUE_TEMPLATE/bug_report.yml` | Create |
| `.github/ISSUE_TEMPLATE/feature_request.yml` | Create |
| `.github/pull_request_template.md` | Create |
| `CITATION.cff` | Create |
| `python/pyproject.toml` | Edit (author, URLs) |
| `.gitignore` | Edit (if stray files need ignoring) |

## Verification

- [ ] All new files render correctly on GitHub
- [ ] `zig build test` still passes
- [ ] CI workflow unaffected (new files are docs/config only)
- [ ] CITATION.cff validates (https://citation-file-format.github.io/cff-initializer-beta/)

## Notes

- CONTRIBUTING.md already exists and is comprehensive
- SECURITY.md already exists
- LICENSE (MIT) already exists
- README.md + README.ja.md already comprehensive

---
- [ ] **DONE** - Phase 1 complete

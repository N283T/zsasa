# Phase 2: v0.1.0 Release

## Goal

v0.1.0 の正式リリースを作成し、Zenodo DOI を取得する。

## Prerequisites

- Phase 1 complete (community files committed)

## Tasks (Sequential)

- [ ] 2.1: CHANGELOG.md の `[Unreleased]` セクション確認・整理
- [ ] 2.2: `git tag -a v0.1.0 -m "Release v0.1.0"` + push
- [ ] 2.3: `gh release create v0.1.0` with CHANGELOG excerpt
- [ ] 2.4: Zenodo integration (GitHub repo と連携、DOI 取得)
- [ ] 2.5: `CITATION.cff` に Zenodo DOI を追記
- [ ] 2.6: README.md にバッジ追加 (Zenodo DOI, CI status, license)
- [ ] 2.7: PyPI 公開 (`pip install zsasa`)

## Key Files

| File | Action |
|------|--------|
| `CHANGELOG.md` | Review/Edit |
| `CITATION.cff` | Edit (add DOI) |
| `README.md` | Edit (add badges) |

## Version Sync (Already Done)

All files already at `0.1.0`:
- `build.zig` line 9: `const version = "0.1.0";`
- `build.zig.zon` line 12: `.version = "0.1.0"`
- `python/pyproject.toml` line 7: `version = "0.1.0"`

## Verification

- [ ] `gh release view v0.1.0` shows correct release
- [ ] Zenodo DOI resolves correctly
- [ ] CITATION.cff DOI matches Zenodo
- [ ] README badges render on GitHub
- [ ] `pip install zsasa` works (if PyPI published)

---
- [ ] **DONE** - Phase 2 complete

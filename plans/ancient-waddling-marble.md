# Phase 1 Execution Plan: Release Preparation

## Context

zsasa (Zig SASA calculator) needs community governance files and repo polish before the v0.1.0 public release. The repo already has CONTRIBUTING.md, SECURITY.md, LICENSE (MIT), and comprehensive READMEs. This phase adds the remaining community files required for open-source best practices and JOSS submission.

**Branch:** `feature/release-prep` (already created)
**Author:** Tsubasa Nagae (no ORCID)

## Current State

- **Existing:** CONTRIBUTING.md, SECURITY.md, LICENSE, README.md, README.ja.md, .github/workflows/
- **Missing:** CODE_OF_CONDUCT.md, issue templates, PR template, CITATION.cff
- **Stray files:** `echo/` and `time/` (empty dirs), `output.json` (already gitignored)
- **pyproject.toml:** Needs author name/email update + Documentation/Issues URLs

## Batch 1: Community Files (parallel creation)

| # | File | Content |
|---|------|---------|
| 1 | `CODE_OF_CONDUCT.md` | Contributor Covenant v2.1, contact: GitHub Issues |
| 2 | `.github/ISSUE_TEMPLATE/bug_report.yml` | YAML form: Zig version, OS, repro steps, expected/actual |
| 3 | `.github/ISSUE_TEMPLATE/feature_request.yml` | YAML form: use case, proposal, alternatives |
| 4 | `.github/pull_request_template.md` | Description, type, testing checklist |
| 5 | `CITATION.cff` | CFF 1.2.0, author: Tsubasa Nagae, DOI placeholder for Phase 2 |

## Batch 2: Repo Polish

| # | Action | Detail |
|---|--------|--------|
| 6 | `python/pyproject.toml` | author: `Tsubasa Nagae`, add Documentation + Issues URLs |
| 7 | Remove stray dirs | `rm -rf echo/ time/`, add to .gitignore |
| 8 | GitHub repo settings | `gh api` calls: topics, discussions, homepage URL |

## Batch 3: Community Building

| # | Action | Detail |
|---|--------|--------|
| 9 | Issues + milestones + labels | 5-10 issues, v0.2.0/v1.0.0 milestones, "good first issue" label |

## Batch 4: Verification & Commit

| # | Action | Detail |
|---|--------|--------|
| 10 | `zig build test` | Ensure no regressions |
| 10 | Commit + push + PR | Single commit with all Phase 1 changes |

## Key Files to Modify/Create

- `CODE_OF_CONDUCT.md` (create)
- `.github/ISSUE_TEMPLATE/bug_report.yml` (create)
- `.github/ISSUE_TEMPLATE/feature_request.yml` (create)
- `.github/pull_request_template.md` (create)
- `CITATION.cff` (create)
- `python/pyproject.toml` (edit: lines 11, 41-43)
- `.gitignore` (edit: add echo/, time/)

## Verification

1. `zig build test` passes
2. `zig fmt --check src/` passes
3. All new files are valid (YAML syntax, CFF format)
4. CI unaffected (new files are docs/config only, CI skips `*.md`, `*.yml` templates)

---
- [ ] **DONE** - Phase 1 complete

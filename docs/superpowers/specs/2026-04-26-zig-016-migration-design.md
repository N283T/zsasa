# Zig 0.16 Migration & Native Flate Restoration Design

## Overview

Migrate zsasa from Zig 0.15.2 to Zig 0.16.0 and revert the C-zlib gzip
workaround introduced in #320 back to native `std.compress.flate`.

The work is split across two PRs released independently:

- **PR1**: toolchain bump + Zig 0.16 API migration (issue #342, all
  Zig-0.16-related items)
- **PR2**: gzip C-zlib → native flate revert (issue #342, last three items)

Each PR ships as its own patch release (PR1 → v0.2.9, PR2 → v0.2.10).

## Motivation

- Zig 0.16.0 (released 2026-04-15) introduces "I/O as an Interface" — `std.fs`,
  `std.net`, `std.process`, `std.Thread` sync primitives, and entropy/time APIs
  move behind a passed `Io` handle. Staying on 0.15 accumulates drift each
  release.
- Upstream `flate.Decompress` panic ([ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035))
  is confirmed fixed in 0.16, so the C-zlib workaround can be removed.
- Keep the toolchain current.

Source issues: #342 (this work), #343 (opt-in follow-up ideas, not in scope).

## Non-goals

- SASA algorithm changes unrelated to the migration.
- Performance work beyond what falls out of the 0.16 switch.
- Issue #343 follow-up ideas (`Io.Group`, `io.async`/`io.concurrent`,
  `ArenaAllocator` thread-safe rewrite, atomic file output API,
  `--error-style`/`--multiline-errors`).
- Threading model redesign. `std.Thread.*` calls migrate **mechanically** to
  the 0.16 equivalents; `thread_pool.zig`'s structure is preserved.
- Bench regression checks (separate effort, after both PRs land).
- C ABI changes for `c_api.zig`. Public signatures must not change so the
  Python bindings keep working.

## Scope

### PR1: toolchain bump + 0.16 API migration

Branch: `chore/zig-0.16-migration`

Steps (each one a single commit unless noted):

1. **toolchain bump**
   - `build.zig.zon`: `minimum_zig_version` `0.15.2` → `0.16.0`
   - `flake.nix`: `zig-overlay.packages.${system}."0.15.2"` → `"0.16.0"`,
     recompute `zigDeps` `outputHash`
   - `.github/workflows/{ci,docs,publish}.yml`: setup-zig version
     `0.15.2` → `0.16.0`
   - `README.md` version claim updates

2. **`std.fs.*` → `std.Io` migration** (~48 call sites)
   - Create one `Io.Threaded` instance in `main.zig` and thread `Io` through
     the call graph ("Juicy Main" style).
   - In `c_api.zig`, create `Io.Threaded` inside each entry function and keep
     it internal — public C signatures must not change.
   - Affected files: `mmap_reader.zig`, `compile_dict.zig`, `dcd.zig`,
     `calc.zig`, `batch.zig`, `c_api.zig`, and any other site surfaced by
     `grep -rn 'std\.fs\.'`.

3. **Reader/Writer API migration** (~3 call sites)
   - Three `std.fs.File.Writer.initStreaming` sites in `batch.zig` rewritten
     to the 0.16 API.

4. **`std.Thread.*` mechanical migration** (~20 call sites)
   - `std.Thread.Mutex` → 0.16 equivalent
   - `std.Thread.spawn` / `std.Thread.getCpuCount` / `std.Thread.sleep` → 0.16
     equivalents
   - `thread_pool.zig` structure unchanged. No `Io.Group` / `io.async`
     redesign.
   - Affected files: `thread_pool.zig`, `batch.zig`, `c_api.zig`,
     `shrake_rupley_bitmask.zig`.

5. **`mem.indexOf*` → `mem.find*` rename** (~18 call sites)
   - Mechanical rename.

6. **Unmanaged container migration** (~19 call sites)
   - `StringHashMap` / `HashMap` / `AutoHashMap` → Unmanaged variants with
     `.empty` initialization and op-time allocator.
   - Affected files include `ccd_parser.zig`, `classifier.zig`,
     `json_parser.zig`, `calc.zig`, `batch.zig`, `traj.zig`,
     `toml_classifier_parser.zig`, `classifier_parser.zig`,
     `classifier_ccd.zig`.

7. **`@cImport` → `build.zig`** (1 call site)
   - In 0.16, `@cImport` inside source files is removed. The single occurrence
     in `src/gzip.zig` migrates to a `b.addTranslateC`-based module wired in
     `build.zig`.
   - This is a placeholder change — PR2 deletes the C path entirely.

8. **`@floor`/`@ceil`/`@round`/`@trunc` audit** (~10 call sites)
   - 0.16 changes float→int conversion rules; `@intFromFloat` may need to be
     made explicit in more places.
   - Affected files: `shrake_rupley_bitmask.zig`, `neighbor_list.zig`,
     `bitmask_lut.zig`, `lee_richards.zig`.

9. **Tests, smoke runs, Python bindings fixes** (commits as needed)
   - `zig build test` until green.
   - 9-way smoke matrix (3 classifiers × 3 formats), `batch`, `traj`.
   - `pytest` in `python/` if `c_api.zig` internal changes broke any binding.

### PR2: native flate restoration

Branch: `chore/revert-gzip-to-native-flate` (cut from main after PR1 merges).

Steps:

1. **Rewrite `src/gzip.zig` to native `std.compress.flate`** (1 commit)
   - Public API unchanged: `readGzip`, `readGzipLimited`, `DEFAULT_MAX_SIZE`,
     `GzipError`. Callers (`mmcif_parser.zig`, `pdb_parser.zig`,
     `json_parser.zig`) are not touched.
   - Preserve the 4 GB decompression-bomb upper bound.

2. **Remove zlib from `build.zig`** (1 commit)
   - Drop `b.dependency("zlib", ...)` and both `zlib_artifact` /
     `zlib_pic_artifact`.
   - Drop `mod.linkLibrary(zlib_artifact)`,
     `exe_module.linkLibrary(zlib_artifact)`,
     `lib_module.linkLibrary(zlib_pic_artifact)`.
   - Drop `lib.linkLibC()` if no other C dependency remains (verify before
     deleting).
   - Drop the `b.addTranslateC` module added in PR1 step 7.

3. **Remove zlib from `build.zig.zon`** (1 commit)
   - Delete `.dependencies.zlib` entry.

4. **Verification smoke** (no extra commit unless fixes needed)

## Architecture

### `Io` propagation pattern

A single `Io.Threaded` instance is constructed at `main.zig` startup and
passed to all entry points that touch `std.Io.{File, Dir}`. This is the
"Juicy Main" pattern recommended by the 0.16 release notes.

For `c_api.zig`, which is the FFI boundary used by the Python bindings, each
exported function constructs its own `Io.Threaded` internally so the C ABI
stays unchanged. The lifetime is scoped to the FFI call.

### Thread pool

`thread_pool.zig` keeps its current chunked-queue worker design.
`std.Thread.spawn` / `getCpuCount` / `sleep` calls migrate one-for-one to the
0.16 equivalents (which the 0.16 release notes route through `std.Io`).

### Build system

PR1 step 7 introduces `b.addTranslateC` for `zlib.h`, replacing the
source-level `@cImport`. PR2 step 2 removes that path entirely along with the
zlib library link.

## Verification

Run for both PR1 and PR2:

```
1. zig build                                  # builds
2. zig build test                             # all tests pass
3. zig build -Doptimize=ReleaseFast           # release builds
4. zsasa calc 2oxd.cif.gz                     # 2713 atoms, 15546.61 Å²
5. for clf in ccd naccess oons:
     for fmt in pdb cif sdf:
       zsasa calc <sample.fmt> --classifier $clf
6. zsasa batch <test_dir>/                    # smoke
7. zsasa traj <test.dcd>                      # smoke
8. cd python && pytest                        # bindings pass
```

Step 4 is the headline confirmation case for PR2 (validates that native
flate handles a real-world gzip-compressed mmCIF without panicking).

## Risks

| Risk | Mitigation |
|---|---|
| Guessing 0.16 API names/signatures wrongly | Each step compiles before moving on. Use compiler errors as the source of truth, not memory of the release notes. |
| `Io.Threaded` propagation causes diff blowup at FFI boundary | `c_api.zig` constructs `Io.Threaded` internally per call, keeping the C ABI stable. |
| `flake.nix` `outputHash` recompute breaks the user's local environment | Use the standard `lib.fakeHash` → build → read correct hash → replace flow. Do not touch `~/dotfiles`. |
| Python bindings (`python/`) silently break | `c_api.zig` public signatures are not allowed to change. `pytest` is part of the verification gate. |
| PR2 native flate has edge-case behavior differences | Verify against issue #342's confirmation case (2oxd.cif.gz) plus existing `gzip.zig` tests and every gzipped fixture in `test_data/`. |
| `@floor`/`@ceil`/`@round`/`@trunc` integer-conversion changes silently shift results | Test suite (`lee_richards.zig`, `bitmask_lut.zig`, `shrake_rupley_bitmask.zig`) catches numeric drift. The 9-way smoke matrix also checks absolute SASA values. |
| Long stretches of `zig build test` red during step-by-step refactor | All commits get applied locally before pushing. The PR is only opened once the whole chain is green. |

## Release strategy

- PR1 merge → `/create-release` → v0.2.9 (chore: Zig 0.16 migration)
- PR2 merge → `/create-release` → v0.2.10 (fix: revert gzip to native flate)

CHANGELOG entries written when each release ships.

## References

- Migration issue: #342
- Follow-up ideas issue (out of scope): #343
- Zig 0.16 release notes: https://ziglang.org/download/0.16.0/release-notes.html
- Upstream flate bug confirmed fixed in 0.16: [ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035)
- Original zlib workaround PRs: #319, #320

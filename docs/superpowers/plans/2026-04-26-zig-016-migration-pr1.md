# Zig 0.16 Migration (PR1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Migrate zsasa from Zig 0.15.2 to 0.16.0. Toolchain bump + mechanical API migration (`std.fs` → `std.Io`, `std.Thread` → `std.Io`, `mem.indexOf*` → `mem.find*`, Managed → Unmanaged containers, `@cImport` → `build.zig`, float→int audit). The C-zlib gzip workaround stays in place; PR2 (separate plan) will revert it later.

**Architecture:** "Juicy Main" pattern — construct one `Io.Threaded` in `main.zig` and thread `Io` through the call graph. `c_api.zig` constructs `Io.Threaded` per FFI call to keep public C signatures stable. `thread_pool.zig` keeps its current chunked-queue design; threading primitives migrate one-for-one (no `Io.Group` redesign).

**Tech Stack:** Zig 0.16.0, mlugg/setup-zig CI action, mitchellh/zig-overlay for Nix, allyourcodebase/zlib (still present in PR1, removed in PR2).

**Spec:** `docs/superpowers/specs/2026-04-26-zig-016-migration-design.md`

**Reference:**
- Zig 0.16 release notes: https://ziglang.org/download/0.16.0/release-notes.html
- Issue #342 (migration), #343 (out-of-scope follow-ups), #320 (C-zlib workaround being preserved in PR1)

**Pre-flight:** Branch `chore/zig-0.16-migration` cut from `main`. Zig 0.16.0 binary available locally (verify with `zig version` after switching the toolchain).

---

## File Inventory

| File | Touched in tasks | Reason |
|---|---|---|
| `build.zig.zon` | 1 | `minimum_zig_version` bump |
| `flake.nix` | 1 | zig-overlay version + `zigDeps` outputHash |
| `.github/workflows/ci.yml` | 1 | setup-zig version |
| `.github/workflows/docs.yml` | 1 | setup-zig version |
| `.github/workflows/publish.yml` | 1 | setup-zig version |
| `README.md` | 1 | Version claim |
| `src/main.zig` | 2 | `Io.Threaded` init, propagation to subcommands |
| `src/root.zig` | 2 | Re-export `Io` if needed by library users |
| `src/calc.zig` | 2, 3, 7, 8 | Subcommand entry, `std.fs` 2 sites, find rename, Unmanaged map |
| `src/batch.zig` | 2, 3, 4, 5, 6, 7, 8 | Many `std.fs`, Reader/Writer, Mutex, spawn, find rename, Unmanaged map |
| `src/traj.zig` | 2, 3, 7, 8 | Subcommand entry, `std.fs`, find rename, Unmanaged map |
| `src/compile_dict.zig` | 2, 3 | Subcommand entry, `std.fs` 2 sites |
| `src/c_api.zig` | 3, 6 | Internal `Io.Threaded`, threading |
| `src/mmap_reader.zig` | 3 | `std.fs.cwd().openFile` |
| `src/dcd.zig` | 3 | `std.fs.File`, `std.fs.cwd().openFile` |
| `src/ccd_parser.zig` | 8 | `StringHashMap` → Unmanaged |
| `src/classifier.zig` | 8 | `HashMap` → Unmanaged |
| `src/classifier_parser.zig` | 8 | `StringHashMap` → Unmanaged |
| `src/classifier_ccd.zig` | 8 | `RuntimeMap` typedef |
| `src/toml_classifier_parser.zig` | 8 | `StringHashMap` → Unmanaged |
| `src/json_parser.zig` | 8 | `AutoHashMap` → Unmanaged |
| `src/thread_pool.zig` | 6 | `std.Thread.spawn`, `getCpuCount` |
| `src/shrake_rupley_bitmask.zig` | 6, 10 | `getCpuCount`, `@round` |
| `src/neighbor_list.zig` | 10 | `@ceil` → int |
| `src/bitmask_lut.zig` | 10 | `@round` |
| `src/lee_richards.zig` | 10 | `@floor`, `@mod` chain |
| `src/gzip.zig` | 9 | Remove `@cImport`, `@import("zlib_c")` |
| `build.zig` | 9 | Add `b.addTranslateC` for `zlib.h` |
| `python/tests/*` | 13 | Bindings smoke (only if breakage surfaces) |

---

### Task 1: Toolchain bump

Update version pins everywhere. Build will not yet pass (subsequent tasks fix the API breakage), but the toolchain is now 0.16.

**Files:**
- Modify: `build.zig.zon`
- Modify: `flake.nix`
- Modify: `.github/workflows/ci.yml`
- Modify: `.github/workflows/docs.yml`
- Modify: `.github/workflows/publish.yml`
- Modify: `README.md`

- [ ] **Step 1: Bump `minimum_zig_version` in `build.zig.zon`**

In `build.zig.zon`, change the line:
```zig
.minimum_zig_version = "0.15.2",
```
to:
```zig
.minimum_zig_version = "0.16.0",
```

- [ ] **Step 2: Bump zig-overlay version in `flake.nix`**

In `flake.nix`, change:
```nix
zig = zig-overlay.packages.${system}."0.15.2";
```
to:
```nix
zig = zig-overlay.packages.${system}."0.16.0";
```

- [ ] **Step 3: Stub `outputHash` for refresh**

Still in `flake.nix`, replace the existing `outputHash` line:
```nix
outputHash = "sha256-wxgdxQiNj7hOTagippwrDkeiZ5RZsagvJI67T28jf04=";
```
with:
```nix
outputHash = "sha256-AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=";
```

- [ ] **Step 4: Compute the real `zigDeps` hash**

Run:
```bash
nix build .#zsasa 2>&1 | grep -A2 "got:" | head -3
```
Expected output includes a line like `got: sha256-XXXX...XXXX=`. Copy the full `sha256-...=` value into `flake.nix` replacing the AAAA stub.

- [ ] **Step 5: Bump setup-zig version in CI workflows**

In each of `.github/workflows/ci.yml`, `.github/workflows/docs.yml`, `.github/workflows/publish.yml`, find every occurrence of:
```yaml
        with:
          version: 0.15.2
```
and change `0.15.2` to `0.16.0`.

Verify all sites updated:
```bash
grep -rn "0.15.2" .github/workflows/
```
Expected: no output.

- [ ] **Step 6: Update README version claim**

In `README.md`, find and replace any `Zig 0.15.2+` (or similar text claiming the minimum Zig version) with `Zig 0.16.0+`.

Verify:
```bash
grep -n "0.15" README.md
```
Expected: no output.

- [ ] **Step 7: Confirm Zig 0.16 binary is active**

Run:
```bash
zig version
```
Expected output: `0.16.0` (you may need `nix develop` or a shell update — coordinate with the user before installing system Zig).

- [ ] **Step 8: Commit**

```bash
git add build.zig.zon flake.nix .github/workflows/ README.md
git commit -m "chore: bump toolchain to Zig 0.16.0"
```

---

### Task 2: Add `Io.Threaded` to `main.zig` and propagate

Construct one `Io.Threaded` instance and pass `Io` to every subcommand handler. The handlers' signatures change to accept `io: std.Io`.

**Files:**
- Modify: `src/main.zig`
- Modify: `src/calc.zig` (function signature change)
- Modify: `src/batch.zig` (function signature change)
- Modify: `src/traj.zig` (function signature change)
- Modify: `src/compile_dict.zig` (function signature change)

- [ ] **Step 1: Construct `Io.Threaded` in `main()`**

In `src/main.zig`, near the top of `main()` (right after the GPA setup), insert:
```zig
var threaded: std.Io.Threaded = .init_single_threaded;
const io = threaded.io();
```

If parallel batch/traj need multi-threaded `Io`, switch to `std.Io.Threaded.init(.{ .allocator = allocator })` instead of `init_single_threaded` — verify by checking 0.16 release notes for the exact constructor.

- [ ] **Step 2: Pass `io` to each subcommand handler**

Update each subcommand call in `src/main.zig`:
```zig
calc.run(allocator, calc_args) catch |err| { ... };
```
to:
```zig
calc.run(allocator, io, calc_args) catch |err| { ... };
```

Apply the same pattern for `batch.run`, `traj.run`, `compile_dict.run` (find each call site).

- [ ] **Step 3: Update `calc.run` signature**

In `src/calc.zig`, find:
```zig
pub fn run(allocator: std.mem.Allocator, args: CalcArgs) !void {
```
Change to:
```zig
pub fn run(allocator: std.mem.Allocator, io: std.Io, args: CalcArgs) !void {
```

- [ ] **Step 4: Update `batch.run` signature**

In `src/batch.zig`, apply the same `io: std.Io` parameter addition to `pub fn run`.

- [ ] **Step 5: Update `traj.run` signature**

In `src/traj.zig`, apply the same.

- [ ] **Step 6: Update `compile_dict.run` signature**

In `src/compile_dict.zig`, apply the same.

- [ ] **Step 7: Compile (expect failures, but for the right reason)**

Run:
```bash
zig build 2>&1 | head -60
```
Expected: errors are now about missing `io` arguments at the call sites in `calc.zig`/`batch.zig`/etc. — not about `main.zig` itself. If `main.zig` errors persist, fix them first.

- [ ] **Step 8: Commit**

```bash
git add src/main.zig src/calc.zig src/batch.zig src/traj.zig src/compile_dict.zig
git commit -m "refactor: introduce Io.Threaded in main and thread io through subcommands"
```

---

### Task 3: Migrate `std.fs.cwd()` and `std.fs.{File,Dir}` operations

This is the largest single rewrite. ~48 call sites across many files. Apply the migration pattern uniformly.

**Migration patterns:**

| Old (0.15) | New (0.16) |
|---|---|
| `std.fs.cwd().openFile(path, .{})` | `std.Io.Dir.cwd().openFile(io, path, .{})` |
| `std.fs.cwd().createFile(path, .{})` | `std.Io.Dir.cwd().createFile(io, path, .{})` |
| `std.fs.cwd().openDir(path, .{...})` | `std.Io.Dir.cwd().openDir(io, path, .{...})` |
| `std.fs.cwd().makePath(path)` | `std.Io.Dir.cwd().createDirPath(io, path)` |
| `std.fs.File` (type) | `std.Io.File` |
| `std.fs.File.stdout()` | `std.Io.File.stdout()` (no `io` needed for static accessor; verify) |
| `std.fs.File.stderr()` | `std.Io.File.stderr()` (verify) |
| `std.fs.path.join(allocator, &.{...})` | `std.Io.Dir.path.join(allocator, &.{...})` (verify) |
| `file.close()` | `file.close(io)` |

**Files:**
- Modify: `src/mmap_reader.zig` (1 site)
- Modify: `src/compile_dict.zig` (2 sites)
- Modify: `src/dcd.zig` (2 sites — type + open)
- Modify: `src/calc.zig` (2 sites — pass `io` from `run` into helpers)
- Modify: `src/batch.zig` (~12 sites)
- Modify: `src/c_api.zig` (~6 sites — but construct `Io.Threaded` internally)
- Modify: `src/traj.zig` (verify with grep)

- [ ] **Step 1: List every site to migrate**

Run:
```bash
grep -rn "std\.fs\." src/ > /tmp/fs-sites.txt
wc -l /tmp/fs-sites.txt
```
Expected: ~48 lines.

- [ ] **Step 2: Migrate `mmap_reader.zig`**

In `src/mmap_reader.zig` find:
```zig
const file = try std.fs.cwd().openFile(path, .{});
```
Add `io: std.Io` parameter to the enclosing function (e.g., `pub fn read(allocator, path)` becomes `pub fn read(allocator, io, path)`), then change to:
```zig
const file = try std.Io.Dir.cwd().openFile(io, path, .{});
```
Update all callers of this function to pass `io`.

- [ ] **Step 3: Migrate `compile_dict.zig`**

In `src/compile_dict.zig` change both `std.fs.cwd().openFile(...)` and `std.fs.cwd().createFile(...)` lines to the `std.Io.Dir.cwd()...` form, passing the `io` parameter from `run()`.

- [ ] **Step 4: Migrate `dcd.zig`**

In `src/dcd.zig`:
- The struct field `file: std.fs.File` becomes `file: std.Io.File`.
- The `openFile` call gets `io` added.
- The struct's `init` function signature gains `io: std.Io`.
- Add `io: std.Io` field to the struct so `read`/`close` can use it without re-passing.

- [ ] **Step 5: Migrate `calc.zig`**

The two `std.fs.cwd().openFile(path, .{})` sites in `calc.zig` (line ~657 and ~1066) each gain `io` access via the `run` parameter. Pass `io` through helper functions as needed.

- [ ] **Step 6: Migrate `batch.zig`**

`batch.zig` has ~12 `std.fs.*` sites. Apply the same patterns. For struct fields like `file: std.fs.File` (line ~826), update to `std.Io.File`. For `std.fs.File.stdout()` / `stderr()`, keep the same call but verify the 0.16 path with the compiler.

For `std.fs.path.join(allocator, &.{...})`, change to `std.Io.Dir.path.join(allocator, &.{...})` (verify exact namespace via compiler).

- [ ] **Step 7: Migrate `c_api.zig` (FFI boundary, internal Io)**

In `src/c_api.zig`, every exported function that touches files must construct its own `Io.Threaded` (so the public C signature is unchanged):
```zig
export fn zsasa_calc_xxx(...) c_int {
    var threaded: std.Io.Threaded = .init_single_threaded;
    const io = threaded.io();
    // ... existing logic, now passing io to internal calls ...
}
```
The C ABI is preserved.

- [ ] **Step 8: Migrate `traj.zig` if it has any `std.fs.*` sites**

```bash
grep -n "std\.fs\." src/traj.zig
```
Apply the same patterns. Use `io` from the new `run` signature.

- [ ] **Step 9: Compile**

```bash
zig build 2>&1 | tee /tmp/build.log | head -80
```
Iterate on remaining errors. Common categories: callers of changed-signature helpers need updating; `file.close()` needs `io` argument.

- [ ] **Step 10: Commit**

```bash
git add src/
git commit -m "refactor: migrate std.fs to std.Io for Zig 0.16"
```

---

### Task 4: Reader/Writer API migration

Three `std.fs.File.Writer.initStreaming` sites in `batch.zig`. The 0.16 Reader/Writer API replaces `initStreaming` (verify exact name with the compiler).

**Files:**
- Modify: `src/batch.zig` (lines ~39, ~850, ~948)

- [ ] **Step 1: Find each site**

```bash
grep -n "Writer\.initStreaming\|Reader\.init" src/batch.zig
```
Expected: 3 lines.

- [ ] **Step 2: Migrate the stderr writer (line ~39)**

Locate:
```zig
const stderr_file = std.fs.File.stderr();
var w = std.fs.File.Writer.initStreaming(stderr_file, &buf);
```
Replace with the 0.16 idiom — most likely:
```zig
const stderr_file = std.Io.File.stderr();
var w = stderr_file.writer(&buf);
```
The exact name (`writer` vs `writerStreaming`) is whatever the compiler accepts. Try `writer` first; if that signature doesn't take a buffer, try `writerStreaming`.

- [ ] **Step 3: Migrate the JsonlWriter (line ~850)**

In the `JsonlWriter` struct's writer construction, apply the same pattern as Step 2, but using the `self.file` rather than `stderr()`.

- [ ] **Step 4: Migrate the jsonl file writer setup (line ~948)**

Same pattern, using the local `jf` file handle.

- [ ] **Step 5: Compile and iterate**

```bash
zig build 2>&1 | head -40
```
Fix any signature mismatches surfaced.

- [ ] **Step 6: Commit**

```bash
git add src/batch.zig
git commit -m "refactor: migrate File.Writer.initStreaming to 0.16 API"
```

---

### Task 5: `std.Thread.Mutex` → `std.Io.Mutex`

Single site to migrate.

**Files:**
- Modify: `src/batch.zig` (line ~825)

- [ ] **Step 1: Replace the type**

In `src/batch.zig` find:
```zig
mutex: std.Thread.Mutex = .{},
```
Replace with:
```zig
mutex: std.Io.Mutex = .{},
```

- [ ] **Step 2: Update lock/unlock sites if API differs**

Search for `mutex.lock()` and `mutex.unlock()` in `src/batch.zig`. If 0.16's `Io.Mutex` requires an `io` parameter for these calls, add it.

- [ ] **Step 3: Compile**

```bash
zig build 2>&1 | head -20
```

- [ ] **Step 4: Commit**

```bash
git add src/batch.zig
git commit -m "refactor: migrate Thread.Mutex to Io.Mutex"
```

---

### Task 6: `std.Thread.spawn` / `getCpuCount` / `sleep` mechanical migration

~20 call sites across `thread_pool.zig`, `batch.zig`, `c_api.zig`, `shrake_rupley_bitmask.zig`. Mechanical replacement. **No `Io.Group` redesign** — preserve the existing chunked-queue worker structure.

**Migration patterns (verify exact 0.16 names with the compiler):**

| Old (0.15) | New (0.16) candidate |
|---|---|
| `std.Thread.spawn(.{}, fn, args)` | `try io.async(fn, args)` (returns `Future`) — but this changes shape; if you need a join handle slice, use `std.Io.Thread.spawn(io, .{}, fn, args)` if available |
| `std.Thread.getCpuCount()` | `std.Io.cpuCount()` or `io.cpuCount()` (verify) |
| `std.Thread.sleep(ns)` | `try io.sleep(.fromNanoseconds(ns), .awake)` (verify) |
| `std.Thread` (type, e.g., `[]std.Thread`) | `std.Io.Thread` (verify) |

If `std.Thread.spawn` no longer exists in any form, the chunked-queue model needs to be implemented over `Io.Group`. **Stop and check with the user before doing that** — the spec explicitly defers `Io.Group` redesign to issue #343.

**Files:**
- Modify: `src/thread_pool.zig`
- Modify: `src/batch.zig`
- Modify: `src/c_api.zig`
- Modify: `src/shrake_rupley_bitmask.zig`

- [ ] **Step 1: List every site**

```bash
grep -rn "std\.Thread" src/
```

- [ ] **Step 2: Migrate `thread_pool.zig`**

`thread_pool.zig` has the core spawn pattern:
```zig
threads: []std.Thread,
...
const cpu_count = try std.Thread.getCpuCount();
const threads = try allocator.alloc(std.Thread, actual_threads);
thread.* = try std.Thread.spawn(.{}, workerLoop, .{ self, i });
```
Replace `std.Thread` → `std.Io.Thread` (or whatever the compiler accepts). Add `io: std.Io` to `ThreadPool.init` and use `io.cpuCount()` / `io.spawn(...)`. Update callers.

- [ ] **Step 3: Migrate `batch.zig` Thread sites**

`batch.zig` has spawn / getCpuCount / sleep sites. Apply the same pattern using `io` from the `run` signature (already added in Task 2).

For `std.Thread.sleep(50 * std.time.ns_per_ms)`, replace with `try io.sleep(.fromMillis(50), .awake)` (verify exact constructor).

- [ ] **Step 4: Migrate `c_api.zig` Thread sites**

`c_api.zig` has multiple spawn sites for batch worker functions. Apply the migration using the `Io.Threaded` already constructed inside each FFI entry. The threading allocator (`c_allocator`) stays.

- [ ] **Step 5: Migrate `shrake_rupley_bitmask.zig`**

Two `std.Thread.getCpuCount()` sites. If they're in helpers that don't currently take `io`, add the parameter and update callers (likely `calc.zig`).

- [ ] **Step 6: Compile**

```bash
zig build 2>&1 | head -40
```
If `std.Io.Thread.spawn` does not exist and the only async primitive is `Io.Group`, **stop and consult with user** — that crosses into issue #343 territory.

- [ ] **Step 7: Commit**

```bash
git add src/
git commit -m "refactor: migrate std.Thread.* primitives to std.Io.* equivalents"
```

---

### Task 7: Rename `mem.indexOf*` → `mem.find*`

Mechanical bulk rename. ~18 sites.

**Migration table (verify each new name with the compiler):**

| Old | New |
|---|---|
| `std.mem.indexOf` | `std.mem.find` |
| `std.mem.indexOfScalar` | `std.mem.findScalar` |
| `std.mem.indexOfAny` | `std.mem.findAny` |
| `std.mem.lastIndexOf` | `std.mem.findLast` |
| `std.mem.lastIndexOfScalar` | `std.mem.findLastScalar` |

**Files:**
- All files with matches from `grep -rn "indexOf\|lastIndexOf" src/`

- [ ] **Step 1: List every site**

```bash
grep -rn "std\.mem\.indexOf\|std\.mem\.lastIndexOf" src/
```

- [ ] **Step 2: Apply renames per file**

Use search-and-replace per affected file. Example for `src/calc.zig`:
```bash
# inspect, then apply
grep -n "std\.mem\.indexOf" src/calc.zig
```
Manually edit each file (do not blindly sed — names differ in arity).

For each match, check whether 0.16 calls it `find` (most likely), `findFirst`, or something else, by trying one and reading the compiler error.

- [ ] **Step 3: Compile**

```bash
zig build 2>&1 | head -40
```
If a name is wrong, the compiler will suggest a candidate. Apply.

- [ ] **Step 4: Commit**

```bash
git add src/
git commit -m "refactor: rename mem.indexOf* to mem.find* for Zig 0.16"
```

---

### Task 8: Migrate Managed → Unmanaged containers

Convert remaining managed `StringHashMap` / `HashMap` / `AutoHashMap` to Unmanaged variants with `.empty` initialization. ~19 sites.

**Migration pattern:**

Old (0.15 managed):
```zig
var map = std.StringHashMap(V).init(allocator);
defer map.deinit();
try map.put("k", v);
const got = map.get("k");
```

New (0.16 unmanaged with op-time allocator):
```zig
var map: std.StringHashMapUnmanaged(V) = .empty;
defer map.deinit(allocator);
try map.put(allocator, "k", v);
const got = map.get("k");
```

Note: 0.16 may have removed the managed variants entirely, in which case the rename is forced (compiler errors will surface every site).

**Files:**
- Modify: `src/calc.zig` (StringHashMap line ~794)
- Modify: `src/batch.zig` (StringHashMap line ~382)
- Modify: `src/ccd_parser.zig` (3 StringHashMap sites: struct field, builder, atom_name_map)
- Modify: `src/classifier_ccd.zig` (RuntimeMap typedef line ~786)
- Modify: `src/classifier_parser.zig` (StringHashMap line ~118)
- Modify: `src/classifier.zig` (HashMap line ~170)
- Modify: `src/toml_classifier_parser.zig` (StringHashMap line ~72)
- Modify: `src/json_parser.zig` (AutoHashMap line ~157)
- Modify: `src/traj.zig` (StringHashMap line ~1119)

- [ ] **Step 1: List every site**

```bash
grep -rn "std\.StringHashMap\|std\.HashMap\|std\.AutoHashMap" src/
```
Expected: ~19 lines (some are typedef + init pairs in same file).

- [ ] **Step 2: Migrate `ccd_parser.zig` (the most complex)**

`ccd_parser.zig` has 3 sites that share a builder pattern:
```zig
components: std.StringHashMap(StoredComponent),
...
.components = std.StringHashMap(StoredComponent).init(allocator),
```
Convert to:
```zig
components: std.StringHashMapUnmanaged(StoredComponent) = .empty,
...
// remove the .init(allocator) — field default `.empty` is sufficient
```
Update every `.put`, `.deinit`, `.get` call site to pass `allocator`. The function bodies that mutate the map need an allocator parameter (already present in most cases).

- [ ] **Step 3: Migrate `classifier.zig`**

The `HashMap(AtomKey, AtomProperties, AtomKeyContext, 80)` becomes `HashMapUnmanaged(AtomKey, AtomProperties, AtomKeyContext, 80)`. Apply same `.empty` + op-time allocator pattern.

- [ ] **Step 4: Migrate the remaining sites (`calc.zig`, `batch.zig`, `traj.zig`, `classifier_parser.zig`, `classifier_ccd.zig`, `toml_classifier_parser.zig`, `json_parser.zig`)**

Apply the same pattern uniformly. Each site: change type, replace `.init(allocator)` with `= .empty`, update `.put`/`.deinit` call sites.

- [ ] **Step 5: Compile**

```bash
zig build 2>&1 | head -40
```
Common errors: forgot to pass allocator to a `.put` or `.deinit`. Fix per error.

- [ ] **Step 6: Run tests for the affected modules**

```bash
zig build test 2>&1 | head -40
```
Map ordering / iteration changes are unlikely but possible — investigate if any test fails.

- [ ] **Step 7: Commit**

```bash
git add src/
git commit -m "refactor: migrate managed HashMaps to Unmanaged variants"
```

---

### Task 9: Move `@cImport` from `gzip.zig` to `build.zig`

In 0.16, `@cImport` inside source files is removed; translate-c is exposed as a build step. This is a **placeholder** change — PR2 will delete the entire C path.

**Files:**
- Modify: `build.zig` (add `b.addTranslateC` for `zlib.h`)
- Modify: `src/gzip.zig` (replace `@cImport(@cInclude("zlib.h"))` with `@import("zlib_c")`)

- [ ] **Step 1: Add the translate-c step in `build.zig`**

In `build.zig`, after `const zlib_artifact = zlib_dep.artifact("z");`, add:
```zig
const zlib_c = b.addTranslateC(.{
    .root_source_file = zlib_dep.path("zlib.h"),
    .target = target,
    .optimize = optimize,
});
const zlib_c_mod = zlib_c.createModule();
```

If `zlib_dep.path("zlib.h")` doesn't resolve (the upstream package may not expose the header at root), use a small wrapper header in the repo:
```bash
mkdir -p src/c
echo '#include <zlib.h>' > src/c/zlib_wrapper.h
```
Then:
```zig
const zlib_c = b.addTranslateC(.{
    .root_source_file = b.path("src/c/zlib_wrapper.h"),
    .target = target,
    .optimize = optimize,
});
zlib_c.addIncludePath(zlib_dep.path("."));
const zlib_c_mod = zlib_c.createModule();
```

- [ ] **Step 2: Wire `zlib_c_mod` into modules that need it**

In `build.zig`, add to `mod` (the library module), `exe_module`, and `lib_module`:
```zig
mod.addImport("zlib_c", zlib_c_mod);
exe_module.addImport("zlib_c", zlib_c_mod);
lib_module.addImport("zlib_c", zlib_c_mod);
```

- [ ] **Step 3: Update `gzip.zig` to import the module**

In `src/gzip.zig` change:
```zig
const c = @cImport(@cInclude("zlib.h"));
```
to:
```zig
const c = @import("zlib_c");
```

- [ ] **Step 4: Compile**

```bash
zig build 2>&1 | head -40
```

- [ ] **Step 5: Run gzip tests**

```bash
zig build test 2>&1 | grep -A3 "gzip"
```
Expected: gzip-related tests pass (PR2 will rewrite this module entirely).

- [ ] **Step 6: Commit**

```bash
git add build.zig src/gzip.zig src/c/ 2>/dev/null
git commit -m "refactor: move @cImport for zlib.h to build.zig translate-c"
```

---

### Task 10: `@floor`/`@ceil`/`@round`/`@trunc` integer-conversion audit

In 0.16, these intrinsics gain the ability to convert directly to integer types, replacing the `@intFromFloat(@ceil(x))` pattern.

**Files:**
- Modify: `src/shrake_rupley_bitmask.zig` (lines 58, 109)
- Modify: `src/neighbor_list.zig` (lines 81-83)
- Modify: `src/bitmask_lut.zig` (line 178)
- Modify: `src/lee_richards.zig` (lines 1361, 1362, 1546, 1547)

- [ ] **Step 1: List every site**

```bash
grep -rn "@floor\|@ceil\|@round\|@trunc" src/
```
Expected: ~10 lines.

- [ ] **Step 2: Audit `neighbor_list.zig` (the most likely candidate to simplify)**

In `src/neighbor_list.zig` lines 81-83, the current code is:
```zig
const nx = @max(1, @as(usize, @intFromFloat(@ceil((x_max - x_min) / cell_size))));
const ny = @max(1, @as(usize, @intFromFloat(@ceil((y_max - y_min) / cell_size))));
const nz = @max(1, @as(usize, @intFromFloat(@ceil((z_max - z_min) / cell_size))));
```
In 0.16, this can be simplified — try:
```zig
const nx: usize = @max(1, @ceil((x_max - x_min) / cell_size));
const ny: usize = @max(1, @ceil((y_max - y_min) / cell_size));
const nz: usize = @max(1, @ceil((z_max - z_min) / cell_size));
```
If the compiler accepts this (with the result type inferred from the usize annotation), it's correct. If not, fall back to the explicit form.

- [ ] **Step 3: Audit `shrake_rupley_bitmask.zig`, `bitmask_lut.zig`**

These use `@round` for binning. The result is then cast to integer downstream. Simplify if the 0.16 inference allows; otherwise leave unchanged (still compiles in 0.16).

- [ ] **Step 4: Audit `lee_richards.zig`**

The `@floor`/`@mod` chains compute float coordinates from integer indices — they stay float, no int conversion. These should just need to keep compiling. Verify with:
```bash
zig build 2>&1 | grep "lee_richards"
```

- [ ] **Step 5: Run targeted tests for numerical correctness**

```bash
zig build test 2>&1 | tail -40
```
Specifically watch tests in `lee_richards`, `shrake_rupley_bitmask`, `bitmask_lut`, `neighbor_list`. If any fail with numerical drift, revert that file's audit changes and use the explicit form.

- [ ] **Step 6: Commit**

```bash
git add src/
git commit -m "refactor: audit float-to-int conversions for Zig 0.16"
```

---

### Task 11: Full `zig build test` green

After all migrations land, the full test suite must pass.

- [ ] **Step 1: Build**

```bash
zig build 2>&1 | tail -20
```
Expected: clean build, no errors.

- [ ] **Step 2: Run all tests**

```bash
zig build test 2>&1 | tail -40
```
Expected: all tests pass.

- [ ] **Step 3: Build release**

```bash
zig build -Doptimize=ReleaseFast 2>&1 | tail -10
```
Expected: clean release build.

- [ ] **Step 4: If any test fails**

Investigate the failure. Common causes:
- HashMap iteration order changed → sort before comparing
- File API behavior change (e.g., `read` returning differently when EOF) → adjust expectations
- Threading: if migration accidentally serialized parallel code, tests should still pass but slowly

Do **not** mark tests `.skip()`. Fix or revert.

- [ ] **Step 5: Commit any fixes**

```bash
git add src/
git commit -m "fix: address test failures from Zig 0.16 migration"
```

---

### Task 12: Smoke test 9-way matrix + batch + traj

Manual smoke runs to verify real-world behavior beyond unit tests.

**Note on `2oxd.cif.gz`:** issue #342 specifies a `2oxd.cif.gz` confirmation case (2713 atoms, 15546.61 Å²). It's not in the repo. Either (a) download from RCSB:
```bash
curl -L -o /tmp/2oxd.cif.gz "https://files.rcsb.org/download/2oxd.cif.gz"
```
or (b) accept that PR2's verification is the canonical 2oxd check, and use repo fixtures here. **Use (a) for PR1's PR description evidence.**

- [ ] **Step 1: Build the binary**

```bash
zig build -Doptimize=ReleaseFast
ls -la zig-out/bin/zsasa
```

- [ ] **Step 2: 9-way classifier × format matrix on a small fixture**

```bash
ZSASA=./zig-out/bin/zsasa
for clf in ccd naccess oons; do
  for input in examples/1ubq.pdb examples/1ubq.cif test_data/ethanol_v2000.sdf; do
    echo "=== $clf × $input ==="
    $ZSASA calc "$input" --classifier "$clf" 2>&1 | tail -3
  done
done
```
Expected: all 9 runs print a "Total SASA" line without errors.

- [ ] **Step 3: Smoke `batch`**

```bash
./zig-out/bin/zsasa batch examples/ 2>&1 | tail -10
```
Expected: completes without error, prints per-file SASA.

- [ ] **Step 4: Smoke `traj`**

```bash
./zig-out/bin/zsasa traj test_data/1l2y.pdb test_data/1l2y.dcd 2>&1 | tail -10
```
Expected: completes, prints frame SASA values.

- [ ] **Step 5: Headline `2oxd.cif.gz` check**

```bash
curl -sL -o /tmp/2oxd.cif.gz "https://files.rcsb.org/download/2oxd.cif.gz"
./zig-out/bin/zsasa calc /tmp/2oxd.cif.gz 2>&1 | tail -5
```
Expected: 2713 atoms, total SASA 15546.61 Å² (allowing small rounding tolerance).
**Note:** PR1 still uses C-zlib for gzip, so a non-zero result confirms migration didn't break gzip plumbing. PR2 will redo this same check against native flate.

- [ ] **Step 6: Capture output for the PR description**

Save the output of Steps 2-5 to a file for inclusion in the PR body:
```bash
{
  echo "## Smoke matrix"
  ./zig-out/bin/zsasa calc examples/1ubq.pdb 2>&1 | tail -3
  # ... etc
} > /tmp/smoke-evidence.txt
```

---

### Task 13: Python bindings smoke

The Python package wraps `c_api.zig`. Since `c_api.zig` public signatures are unchanged, bindings should work — but verify.

**Files:**
- Modify (only if breakage): `python/zsasa/`, `python/tests/`

- [ ] **Step 1: Build the shared library**

```bash
zig build -Doptimize=ReleaseFast
ls -la zig-out/lib/libzsasa*
```

- [ ] **Step 2: Install Python package locally**

```bash
cd python
uv sync
```

- [ ] **Step 3: Run the test suite**

```bash
cd python
uv run pytest -v 2>&1 | tail -40
```
Expected: all tests pass.

- [ ] **Step 4: If failures surface**

If `c_api.zig` internals broke ABI subtly (e.g., a struct layout change because `std.Io.File` has different size than `std.fs.File`), investigate and fix in `c_api.zig`. **Do not change the Python wrapper unless the Zig side genuinely needs a new function.**

- [ ] **Step 5: Commit any binding fixes**

```bash
git add python/ src/c_api.zig
git commit -m "fix: restore python bindings compatibility on Zig 0.16"
```

---

### Task 14: Push and open PR1

- [ ] **Step 1: Push branch**

```bash
git push -u origin chore/zig-0.16-migration
```

- [ ] **Step 2: Open PR**

Use the smoke evidence captured in Task 12 Step 6 in the PR body.
```bash
gh pr create --title "chore: migrate to Zig 0.16" --body "$(cat <<'EOF'
## Summary
- Migrate toolchain to Zig 0.16.0 (build.zig.zon, flake.nix, CI workflows)
- Migrate `std.fs.*` to `std.Io.*` ("Juicy Main" pattern with `Io.Threaded`)
- Migrate `std.Thread.*` primitives to `std.Io.*` equivalents (mechanical, no `Io.Group` redesign)
- Rename `mem.indexOf*` → `mem.find*`
- Migrate managed HashMaps to Unmanaged variants
- Move `@cImport` for zlib.h to `build.zig` translate-c (placeholder; PR2 will remove the C path entirely)
- Audit `@floor`/`@ceil`/`@round`/`@trunc` integer conversions

C ABI preserved (`c_api.zig` constructs `Io.Threaded` internally per FFI call). Python bindings continue to work.

PR2 (separate) will revert the C-zlib gzip workaround to native `std.compress.flate`.

Refs: #342

## Smoke evidence
<paste from /tmp/smoke-evidence.txt>

## Test plan
- [x] `zig build` clean
- [x] `zig build test` all pass
- [x] `zig build -Doptimize=ReleaseFast` clean
- [x] `zsasa calc 2oxd.cif.gz` → 2713 atoms, 15546.61 Å² (still via C-zlib in PR1)
- [x] 9-way classifier × format smoke matrix
- [x] `zsasa batch` smoke
- [x] `zsasa traj` smoke
- [x] `python pytest` all pass
EOF
)"
```

- [ ] **Step 3: Wait for CI + review, address feedback per `git-workflow-rules.md`**

- [ ] **Step 4: Ask user before merging**

Per `git-workflow-rules.md`: ask explicit confirmation before `gh pr merge`.

- [ ] **Step 5: Post-merge cleanup**

```bash
git checkout main && git pull
git push origin --delete chore/zig-0.16-migration
git branch -d chore/zig-0.16-migration
```

- [ ] **Step 6: Release v0.2.9 via `/create-release` skill**

After PR1 merges, invoke `/create-release` to ship v0.2.9 (chore: Zig 0.16 migration).

- [ ] **Step 7: Hand off to PR2**

Write `docs/superpowers/plans/<date>-zig-016-migration-pr2.md` for the gzip native-flate revert (separate planning session).

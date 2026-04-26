# Zig 0.16 Migration (PR2) — Native Flate Restoration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Revert the C-zlib gzip workaround introduced in #320 back to native `std.compress.flate`, now that the upstream panic ([ziglang/zig#25035](https://github.com/ziglang/zig/issues/25035)) is fixed in Zig 0.16. Drop the entire `zlib` dependency chain (`build.zig.zon` entry, `b.addTranslateC` step, `linkLibrary` calls, wrapper header).

**Architecture:** `src/gzip.zig` is rewritten on top of `std.compress.flate.Decompress` with `Container.gzip`. The public API (`readGzip`, `readGzipLimited`, `DEFAULT_MAX_SIZE`, `GzipError`) stays byte-for-byte stable so callers (`mmcif_parser.zig`, `pdb_parser.zig`, `json_parser.zig`, `sdf_parser.zig`, `calc.zig`, `batch.zig`, `compile_dict.zig`) need zero changes. The function signature does **not** take `io: std.Io`, so the file is opened through a function-local `std.Io.Threaded.init_single_threaded` (the same pattern `c_api.zig` already uses for FFI entries). The 4 GB streaming decompression-bomb cap is preserved.

**Tech Stack:** Zig 0.16.0, `std.compress.flate` (native), `std.Io.Threaded` for file I/O.

**Spec:** `docs/superpowers/specs/2026-04-26-zig-016-migration-design.md` (PR2 section)

**Reference:**
- Zig 0.16 release notes: https://ziglang.org/download/0.16.0/release-notes.html
- Issue #342 (migration), #320 (the C-zlib workaround being reverted), confirmed-fixed upstream: ziglang/zig#25035
- PR1 plan: `docs/superpowers/plans/2026-04-26-zig-016-migration-pr1.md`

**Pre-flight:**
- Branch `chore/revert-gzip-to-native-flate` cut from latest `main` (already at v0.2.9, PR1 merged).
- `zig version` → `0.16.0`.
- `2oxd.cif.gz` available locally for the headline check (download command in Task 6 if not).
- Subagents executing tasks **must** use the `reference-zig` skill (zigdoc) to look up exact Zig 0.16 standard library API names — especially `std.compress.flate.Decompress`, `Container.gzip`, the new `std.Io.Reader`/`File.Reader` shapes. The compiler error message is the second source of truth when zigdoc returns `UnsupportedZigVersion`.

---

## File Inventory

| File | Touched in tasks | Reason |
|---|---|---|
| `src/gzip.zig` | 1 | Full rewrite: C-zlib → native flate. Public API and tests unchanged. |
| `src/c/zlib_wrapper.h` | 2 | Deleted — no longer needed once translate-c step is gone. |
| `build.zig` | 2 | Drop `zlib_dep`/`zlib_artifact`/`zlib_pic_artifact`, all `linkLibrary(zlib_*)` calls, the `b.addTranslateC` step, the `addImport("zlib_c", ...)` calls, and `lib_module.link_libc` if no other C dependency surfaces. |
| `build.zig.zon` | 3 | Drop `.dependencies.zlib` entry. |

No callers of `gzip.zig` are touched. No CHANGELOG edit in this PR — the release skill (`/create-release` v0.2.10) handles that after merge.

---

### Task 1: Rewrite `src/gzip.zig` with native `std.compress.flate`

Replace the C-zlib implementation with a streaming native-flate decompressor. The public surface (`readGzip`, `readGzipLimited`, `DEFAULT_MAX_SIZE`, `GzipError`) does not change. The 4 GB upper bound is enforced **before** the next chunk is read, matching the C-zlib version.

**Files:**
- Modify: `src/gzip.zig` (full rewrite of the implementation, tests preserved)

- [ ] **Step 1: Verify the Zig 0.16 flate API names with zigdoc**

Use the `reference-zig` skill to check the current names for:
- `std.compress.flate.Decompress` (struct)
- The `Container` enum variant `.gzip` (or whatever the current name is — could be `Container.gzip`, `flate.Container.gzip`, or `flate.Decompress.Container.gzip`)
- `Decompress.init` parameters (typically: input reader, container, dictionary slice)
- `Decompress.reader` field type (the public `std.Io.Reader` exposed by the decompressor)
- The `std.Io.Reader.read(dest)` method (or `readSliceShort`, whatever zero-copy short read returns `usize`)

If zigdoc returns `UnsupportedZigVersion`, fall back to `grep -rn 'pub fn init\|pub const Container' $(zig env | grep -i lib_dir)` against the locally installed std (`zig env --json | jq -r .std_dir`), or trust the compiler errors during Step 3.

- [ ] **Step 2: Read the current gzip.zig contents**

Run:
```bash
cat src/gzip.zig
```
Confirm the current implementation uses `@import("zlib_c")` and the `c.gzopen` / `c.gzread` / `c.gzclose` C API. Note the three existing tests at the bottom (`store block`, `nonexistent file`, `FileTooLarge`) — these stay verbatim.

- [ ] **Step 3: Rewrite `src/gzip.zig` with native flate**

Replace the file contents with the implementation below. The function-local `Io.Threaded` keeps the public signature unchanged (no `io: std.Io` parameter). API names marked **(verify)** must be confirmed via zigdoc or the compiler.

```zig
//! Gzip decompression using native std.compress.flate.
//!
//! Restored after #320's C-zlib workaround. The upstream panic
//! (https://github.com/ziglang/zig/issues/25035) is fixed in Zig 0.16.

const std = @import("std");

pub const GzipError = error{ GzipOpenFailed, GzipReadFailed, FileTooLarge, OutOfMemory };

const CHUNK_SIZE = 64 * 1024; // 64 KB read chunks

/// Default max decompressed size: 4 GB.
pub const DEFAULT_MAX_SIZE: usize = 4 * 1024 * 1024 * 1024;

/// Decompress a gzip file. Caller owns the returned slice.
pub fn readGzip(allocator: std.mem.Allocator, path: []const u8) GzipError![]u8 {
    return readGzipLimited(allocator, path, DEFAULT_MAX_SIZE);
}

/// Decompress a gzip file with a custom size limit. Caller owns the returned slice.
/// The size limit is checked before each chunk read so a malicious file cannot
/// allocate more than `max_size` bytes before the cap is detected.
pub fn readGzipLimited(allocator: std.mem.Allocator, path: []const u8, max_size: usize) GzipError![]u8 {
    var threaded: std.Io.Threaded = .init_single_threaded;
    const io = threaded.io();

    const file = std.Io.Dir.cwd().openFile(io, path, .{}) catch return error.GzipOpenFailed;
    defer file.close(io);

    var file_buf: [CHUNK_SIZE]u8 = undefined;
    var file_reader = file.reader(io, &file_buf); // (verify: file.reader signature)

    // Native gzip decompressor. (verify: Decompress.init parameter shape and Container.gzip path.)
    var decompress: std.compress.flate.Decompress = .init(file_reader.interface(), .gzip, &.{});
    var reader: *std.Io.Reader = &decompress.reader;

    var buf: std.ArrayListUnmanaged(u8) = .empty;
    errdefer buf.deinit(allocator);

    while (true) {
        if (buf.items.len >= max_size) return error.FileTooLarge;
        const room = max_size - buf.items.len;
        const want = @min(CHUNK_SIZE, room);
        try buf.ensureUnusedCapacity(allocator, want);
        const dest = buf.unusedCapacitySlice()[0..want];
        const n = reader.read(dest) catch |err| switch (err) {
            error.EndOfStream => break,
            else => return error.GzipReadFailed,
        };
        if (n == 0) break;
        buf.items.len += n;
    }

    return buf.toOwnedSlice(allocator);
}

// -- Tests --

test "readGzip decompresses gzip store block" {
    const allocator = std.testing.allocator;

    // Minimal gzip containing "Hello world\n" (store block, no compression)
    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    const content = try readGzip(allocator, tmp_path);
    defer allocator.free(content);

    try std.testing.expectEqualStrings("Hello world\n", content);
}

test "readGzip returns GzipOpenFailed for nonexistent file" {
    const allocator = std.testing.allocator;
    const result = readGzip(allocator, "/nonexistent/path/file.gz");
    try std.testing.expectError(error.GzipOpenFailed, result);
}

test "readGzipLimited returns FileTooLarge when limit exceeded" {
    const allocator = std.testing.allocator;

    // Minimal gzip containing "Hello world\n" (12 bytes decompressed)
    const gz_data = [_]u8{
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03,
        0x01, 0x0c, 0x00, 0xf3, 0xff, 0x48, 0x65, 0x6c, 0x6c, 0x6f,
        0x20, 0x77, 0x6f, 0x72, 0x6c, 0x64, 0x0a, 0xd5, 0xe0, 0x39,
        0xb7, 0x0c, 0x00, 0x00, 0x00,
    };

    var tmp_dir = std.testing.tmpDir(.{});
    defer tmp_dir.cleanup();
    try tmp_dir.dir.writeFile(std.testing.io, .{ .sub_path = "test.gz", .data = &gz_data });

    const tmp_path = try tmp_dir.dir.realPathFileAlloc(std.testing.io, "test.gz", allocator);
    defer allocator.free(tmp_path);

    // Limit to 5 bytes — "Hello world\n" is 12 bytes, should fail
    const result = readGzipLimited(allocator, tmp_path, 5);
    try std.testing.expectError(error.FileTooLarge, result);
}
```

**Verification fences (apply during compile):**
- `file.reader(io, &buf)` — if the actual signature is `file.reader(&buf)` (io carried by the file handle) or `file.readerStreaming(io, &buf)`, the compiler will say so. Adapt the call.
- `file_reader.interface()` — if the field is `.reader` instead of a method, change accordingly.
- `std.compress.flate.Decompress.init(reader, .gzip, &.{})` — the third arg may be optional in 0.16, or the container path may be `std.compress.flate.Container.gzip` accessed differently. If the compiler can't find `.gzip` directly, try `std.compress.flate.Container.gzip` first, then `std.compress.gzip.Container.gzip`. Use whichever the compiler accepts.
- `decompress.reader.read(dest)` — if `read` returns `error.EndOfStream` for clean EOF, the `switch (err)` branch above handles it. If 0.16 distinguishes it via `n == 0`, the switch is harmless dead code (still compiles). Either path is correct.

- [ ] **Step 4: Compile**

```bash
zig build 2>&1 | head -40
```
Expected: clean build, no errors. If there are errors, fix per the verification fences in Step 3 and recompile.

- [ ] **Step 5: Run gzip module tests**

```bash
zig build test 2>&1 | grep -A2 -i "gzip"
```
Expected: all 3 gzip tests pass. If any fails:
- `store block` failure → flate decompressor likely doesn't accept gzip wrapper. Wrong `Container` value. Fix.
- `nonexistent file` failure → `openFile` error not getting mapped to `GzipOpenFailed`. Verify the catch.
- `FileTooLarge` failure → cap check is in the wrong place (e.g., checked after read). Move it before `read`.

- [ ] **Step 6: Run full test suite**

```bash
zig build test 2>&1 | tail -20
```
Expected: same pass/skip count as PR1's final state (530 pass, 1 skip per the v0.2.9 record). No regressions.

- [ ] **Step 7: Build release**

```bash
zig build -Doptimize=ReleaseFast 2>&1 | tail -5
```
Expected: clean release build.

- [ ] **Step 8: Commit**

```bash
git add src/gzip.zig
git commit -m "$(cat <<'EOF'
fix: revert gzip to native std.compress.flate

The upstream panic (ziglang/zig#25035) that prompted #320's C-zlib workaround
is fixed in Zig 0.16. Restore the native flate decompressor while preserving
the 4 GB decompression-bomb cap and the public API (readGzip / readGzipLimited
/ DEFAULT_MAX_SIZE / GzipError).

Refs: #342
EOF
)"
```

---

### Task 2: Remove zlib from `build.zig` and delete the C wrapper

Drop every line referring to `zlib_dep`, `zlib_artifact`, `zlib_pic_artifact`, the translate-c step, and the `zlib_c` module. Verify whether `lib_module.link_libc = true` can also be removed (it can if `c_api.zig` has no remaining libc dependency).

**Files:**
- Modify: `build.zig`
- Delete: `src/c/zlib_wrapper.h`
- Delete: `src/c/` directory if it becomes empty after removing the wrapper

- [ ] **Step 1: Read current `build.zig` to confirm exact lines**

```bash
sed -n '8,100p' build.zig
```
Expected layout (from the post-PR1 state):
- Lines 9-13: `const zlib_dep = b.dependency("zlib", .{ ... });`
- Line 14: `const zlib_artifact = zlib_dep.artifact("z");`
- Lines 16-22: `const zlib_pic_dep = b.dependency(...)` and `const zlib_pic_artifact = ...;`
- Line 29: `mod.linkLibrary(zlib_artifact);`
- Line 53: `exe_module.linkLibrary(zlib_artifact);`
- Lines 66-73: `const zlib_c = b.addTranslateC(...);` block + `createModule()`
- Lines 76-77: `mod.addImport("zlib_c", zlib_c_mod);` and `exe_module.addImport(...)`
- Line 89: `lib_module.linkLibrary(zlib_pic_artifact);`
- Line 90: `lib_module.addImport("zlib_c", zlib_c_mod);`

If the line numbers have shifted, find each pattern with grep first:
```bash
grep -n 'zlib\|zlib_c\|linkLibrary\|addTranslateC' build.zig
```

- [ ] **Step 2: Delete the zlib dependency declarations**

In `build.zig`, remove the entire `zlib_dep`/`zlib_pic_dep` block (the lines from the comment `// zlib dependency...` through the end of the `zlib_pic_artifact` declaration). After removal, the build.zig should jump straight from the optimize/target setup into `const mod = b.addModule("zsasa", ...)`.

- [ ] **Step 3: Delete `mod.linkLibrary(zlib_artifact)`**

Remove the single line:
```zig
mod.linkLibrary(zlib_artifact);
```

- [ ] **Step 4: Delete `exe_module.linkLibrary(zlib_artifact)`**

Remove the single line:
```zig
exe_module.linkLibrary(zlib_artifact);
```

- [ ] **Step 5: Delete the translate-c block and the `zlib_c` module wiring**

Remove the entire block starting from the comment about `addTranslateC` through `mod.addImport("zlib_c", zlib_c_mod);` and `exe_module.addImport("zlib_c", zlib_c_mod);`. Specifically, delete:

```zig
// Translate zlib.h to a Zig module via addTranslateC ...
const zlib_c = b.addTranslateC(.{
    .root_source_file = b.path("src/c/zlib_wrapper.h"),
    .target = target,
    .optimize = optimize,
    .link_libc = true,
});
zlib_c.addIncludePath(zlib_dep.path("."));
const zlib_c_mod = zlib_c.createModule();

// Wire zlib_c into the modules that need it
mod.addImport("zlib_c", zlib_c_mod);
exe_module.addImport("zlib_c", zlib_c_mod);
```

- [ ] **Step 6: Delete `lib_module.linkLibrary(zlib_pic_artifact)` and `lib_module.addImport("zlib_c", ...)`**

Remove both lines:
```zig
lib_module.linkLibrary(zlib_pic_artifact);
lib_module.addImport("zlib_c", zlib_c_mod);
```

- [ ] **Step 7: Test whether `lib_module.link_libc = true` is still required**

The shared library exports symbols via `extern fn`, but that does **not** require libc. Try removing `link_libc = true` from the `lib_module` createModule call:

```zig
const lib_module = b.createModule(.{
    .root_source_file = b.path("src/c_api.zig"),
    .target = target,
    .optimize = optimize,
    // .link_libc = true,  ← attempt to delete
    .imports = &.{
        .{ .name = "zxdrfile", .module = zxdrfile_mod },
    },
});
```

Then:
```bash
zig build 2>&1 | tail -20
```

Decision rule:
- **If the build succeeds**: keep `link_libc` removed.
- **If the build fails with linker errors mentioning `__libc_*`, `malloc`, `__stack_chk_*`, `memcpy`, etc.**: restore `.link_libc = true,` and add a one-line comment explaining why (e.g., `// libc still needed for <symbol> in c_api.zig`).

- [ ] **Step 8: Delete the C wrapper header**

```bash
git rm src/c/zlib_wrapper.h
rmdir src/c 2>/dev/null || true   # only succeeds if dir is empty
```

- [ ] **Step 9: Verify clean build**

```bash
zig build 2>&1 | tail -20
```
Expected: clean build.

- [ ] **Step 10: Verify clean release build**

```bash
zig build -Doptimize=ReleaseFast 2>&1 | tail -10
```
Expected: clean release build, `zig-out/bin/zsasa` and `zig-out/lib/libzsasa.{dylib,so,dll}` produced.

- [ ] **Step 11: Run the full test suite**

```bash
zig build test 2>&1 | tail -10
```
Expected: all pass, no regressions vs Task 1's baseline.

- [ ] **Step 12: Commit**

```bash
git add build.zig
git add -A src/c/  # captures the deletion
git commit -m "$(cat <<'EOF'
chore: remove zlib dependency from build.zig

With gzip.zig back on native std.compress.flate, the zlib dependency, the
translate-c step, and the wrapper header are all unused. Drop them.

Refs: #342
EOF
)"
```

---

### Task 3: Remove zlib from `build.zig.zon`

The dependency is no longer fetched. Drop the entry so `zig fetch` and `zig build --fetch` skip it.

**Files:**
- Modify: `build.zig.zon`

- [ ] **Step 1: Confirm the current entry**

```bash
grep -A3 '\.zlib' build.zig.zon
```
Expected output:
```zig
        .zlib = .{
            .url = "git+https://github.com/allyourcodebase/zlib.git?ref=v1.3.1.1-2#c5115f4b69ef660f72a835c6638f80508ef284c7",
            .hash = "zlib-1.3.1-1-ZZQ7ldENAAA7qJjUXP6E6xnRuV-jDL9dyoJFc_eb3zQ6",
        },
```

- [ ] **Step 2: Delete the `.zlib` entry**

In `build.zig.zon` remove the four lines comprising the `.zlib = .{ ... },` entry. The `.dependencies` block should now contain only the `.zxdrfile` entry. The trailing comma on the last remaining entry must be present (Zig allows trailing comma).

After the edit, the dependencies block should look like:
```zig
    .dependencies = .{
        .zxdrfile = .{
            .url = "git+https://github.com/N283T/zxdrfile.git?ref=v0.4.0#47d8c42cff9b098febd1990cc50438c26987f85c",
            .hash = "zxdrfile-0.4.0-JIaQ9mIfAgAHB9bW-_p9g6_gfd0tyG3-MmA0LYJKB8rb",
        },
    },
```

- [ ] **Step 3: Verify clean fetch + build**

```bash
zig build --fetch 2>&1 | tail -5
zig build 2>&1 | tail -10
```
Expected: both succeed, no zlib entry pulled.

- [ ] **Step 4: Run tests**

```bash
zig build test 2>&1 | tail -10
```
Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add build.zig.zon
git commit -m "$(cat <<'EOF'
chore: drop zlib dependency from build.zig.zon

No longer referenced after Task 2's build.zig cleanup.

Refs: #342
EOF
)"
```

---

### Task 4: Headline `2oxd.cif.gz` verification

This is the canonical confirmation case from issue #342. Native flate must successfully decompress a real-world gzip-compressed mmCIF and the SASA result must match.

- [ ] **Step 1: Build the release binary**

```bash
zig build -Doptimize=ReleaseFast
ls -la zig-out/bin/zsasa
```

- [ ] **Step 2: Fetch `2oxd.cif.gz`**

```bash
mkdir -p /tmp/zsasa-pr2
curl -sL -o /tmp/zsasa-pr2/2oxd.cif.gz "https://files.rcsb.org/download/2oxd.cif.gz"
ls -la /tmp/zsasa-pr2/2oxd.cif.gz
```
Expected: file size > 100 KB (compressed mmCIF).

- [ ] **Step 3: Run the headline check**

```bash
./zig-out/bin/zsasa calc /tmp/zsasa-pr2/2oxd.cif.gz 2>&1 | tee /tmp/zsasa-pr2/2oxd-output.txt
```

- [ ] **Step 4: Verify the numbers**

Inspect `/tmp/zsasa-pr2/2oxd-output.txt`. The output **must** contain:
- 2713 atoms
- Total SASA: 15546.61 Å² (allowing rounding tolerance ±0.01 Å²)

If the atom count is wrong, the decompressed CIF is truncated → `readGzipLimited`'s cap or EOF handling is wrong, fix Task 1.
If the SASA is wrong, the cif is corrupted mid-stream → flate decompression silently dropped trailing bytes, fix Task 1's read loop.

- [ ] **Step 5: Save the evidence for the PR description**

```bash
{
  echo "## 2oxd.cif.gz headline check"
  echo ""
  echo '```'
  ./zig-out/bin/zsasa calc /tmp/zsasa-pr2/2oxd.cif.gz 2>&1 | tail -5
  echo '```'
} > /tmp/zsasa-pr2/headline-evidence.md
cat /tmp/zsasa-pr2/headline-evidence.md
```

(No commit — this is verification output for the PR body.)

---

### Task 5: 9-way smoke matrix + batch + traj

Same matrix that PR1 ran. Catches regressions outside the gzip path.

- [ ] **Step 1: Run the 9-way classifier × format matrix**

```bash
ZSASA=./zig-out/bin/zsasa
{
  for clf in ccd naccess oons; do
    for input in examples/1ubq.pdb examples/1ubq.cif test_data/ethanol_v2000.sdf; do
      echo "=== $clf × $input ==="
      "$ZSASA" calc "$input" --classifier "$clf" 2>&1 | tail -3
    done
  done
} | tee /tmp/zsasa-pr2/smoke-matrix.txt
```
Expected: all 9 invocations print a "Total SASA" line with no errors.

If any fixture path is missing, `ls examples/` and `ls test_data/` and substitute the closest-named file. Do not invent inputs.

- [ ] **Step 2: Smoke `batch`**

```bash
./zig-out/bin/zsasa batch examples/ 2>&1 | tail -10 | tee /tmp/zsasa-pr2/batch-smoke.txt
```
Expected: per-file SASA lines, no errors.

- [ ] **Step 3: Smoke `traj`**

```bash
./zig-out/bin/zsasa traj test_data/1l2y.pdb test_data/1l2y.dcd 2>&1 | tail -10 | tee /tmp/zsasa-pr2/traj-smoke.txt
```
Expected: per-frame SASA values, no errors.

If `test_data/1l2y.{pdb,dcd}` is missing, locate the actual trajectory fixtures with `ls test_data/*.dcd` and substitute. Do not invent inputs.

- [ ] **Step 4: Save evidence for the PR description**

```bash
{
  echo "## Smoke matrix"
  echo '```'
  cat /tmp/zsasa-pr2/smoke-matrix.txt
  echo '```'
  echo ""
  echo "## Batch smoke"
  echo '```'
  cat /tmp/zsasa-pr2/batch-smoke.txt
  echo '```'
  echo ""
  echo "## Traj smoke"
  echo '```'
  cat /tmp/zsasa-pr2/traj-smoke.txt
  echo '```'
} >> /tmp/zsasa-pr2/headline-evidence.md
```

(No commit.)

---

### Task 6: Python bindings smoke

`c_api.zig` was not edited, but the shared library no longer links zlib. Verify Python bindings still work.

- [ ] **Step 1: Build the shared library**

```bash
zig build -Doptimize=ReleaseFast
ls -la zig-out/lib/libzsasa*
```

- [ ] **Step 2: Sync Python deps and run tests**

```bash
cd python
uv sync
uv run pytest -v 2>&1 | tee /tmp/zsasa-pr2/pytest.txt | tail -40
```
Expected: all tests pass.

If any test fails because the dynamic loader can't find a symbol that came from zlib (unlikely, since `c_api.zig` doesn't call zlib), the wrapper header was still being pulled into the python build path. Re-check Task 2's deletions.

- [ ] **Step 3: Save evidence**

```bash
{
  echo "## Python bindings"
  echo '```'
  tail -20 /tmp/zsasa-pr2/pytest.txt
  echo '```'
} >> /tmp/zsasa-pr2/headline-evidence.md
```

(No commit.)

---

### Task 7: Push branch and open PR

- [ ] **Step 1: Confirm the commit history**

```bash
git log --oneline main..HEAD
```
Expected: exactly 3 commits, in order:
1. `fix: revert gzip to native std.compress.flate`
2. `chore: remove zlib dependency from build.zig`
3. `chore: drop zlib dependency from build.zig.zon`

If there are extra commits (e.g., a fixup that should have been squashed into Task 1), use `git rebase -i main` to squash before pushing.

- [ ] **Step 2: Push branch**

```bash
git push -u origin chore/revert-gzip-to-native-flate
```

- [ ] **Step 3: Open PR with the smoke evidence**

```bash
gh pr create --title "fix: revert gzip to native std.compress.flate" --body "$(cat <<EOF
## Summary
- Rewrite \`src/gzip.zig\` on top of \`std.compress.flate.Decompress\` with \`Container.gzip\`. Public API (\`readGzip\` / \`readGzipLimited\` / \`DEFAULT_MAX_SIZE\` / \`GzipError\`) and the 4 GB decompression-bomb cap unchanged.
- Drop the \`zlib\` dependency, the \`b.addTranslateC\` step, the \`linkLibrary\` calls, and \`src/c/zlib_wrapper.h\` from \`build.zig\` and \`build.zig.zon\`.
- Reverts the C-zlib workaround introduced in #320; the upstream panic ziglang/zig#25035 is fixed in Zig 0.16.

Callers (\`mmcif_parser.zig\`, \`pdb_parser.zig\`, \`json_parser.zig\`, \`sdf_parser.zig\`, \`calc.zig\`, \`batch.zig\`, \`compile_dict.zig\`) are not touched.

Refs: #342

$(cat /tmp/zsasa-pr2/headline-evidence.md)

## Test plan
- [x] \`zig build\` clean
- [x] \`zig build test\` all pass
- [x] \`zig build -Doptimize=ReleaseFast\` clean
- [x] \`zsasa calc 2oxd.cif.gz\` → 2713 atoms, 15546.61 Å² (native flate)
- [x] 9-way classifier × format smoke matrix
- [x] \`zsasa batch\` smoke
- [x] \`zsasa traj\` smoke
- [x] \`python pytest\` all pass
EOF
)"
```

- [ ] **Step 4: Wait for CI; address feedback per `git-workflow-rules.md`**

Run a final review pass via `pr-review-toolkit:review-pr` once CI is green. Fix any critical/important issues without asking. Push the fixes.

- [ ] **Step 5: Ask user before merging**

Per `git-workflow-rules.md`, post a single message to the user with CI status, review status, and an explicit merge ask. **Do not merge until the user replies "OK"/"merge"/"yes".**

- [ ] **Step 6: Post-merge cleanup**

```bash
git checkout main && git pull
git push origin --delete chore/revert-gzip-to-native-flate
git branch -d chore/revert-gzip-to-native-flate
git remote prune origin
```

- [ ] **Step 7: Release v0.2.10 via `/create-release`**

After the PR merges, invoke `/create-release v0.2.10` to ship the gzip native-flate restoration release. The release skill handles the changelog entry, version bumps (`build.zig`, `build.zig.zon`, `python/pyproject.toml`, `flake.nix`), release PR, and tag push.

- [ ] **Step 8: Close issue #342 (migration tracker)**

After v0.2.10 ships, both PR1 and PR2 are landed. Close #342 with a summary comment:

```bash
PR2_NUMBER=$(gh pr list --state merged --search "fix: revert gzip to native std.compress.flate" --json number --jq '.[0].number')
gh issue close 342 --comment "Migration complete.

- PR1 (Zig 0.16 migration): #345 → released as v0.2.9
- PR2 (gzip → native flate): #${PR2_NUMBER} → released as v0.2.10"
```

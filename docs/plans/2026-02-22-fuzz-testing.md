# Fuzz Testing Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add fuzz tests to cif_tokenizer, pdb_parser, and mmcif_parser using Zig 0.15.2's built-in `std.testing.fuzz()` to catch crashes and panics from malformed input (Issue #161).

**Architecture:** Each parser file gets a fuzz test appended to its existing test section. The fuzz test feeds arbitrary byte sequences from `std.testing.fuzz()` into the parser's `parse()` or tokenizer's `next()` function. Errors are expected and caught; only crashes/panics are failures. Seed corpus uses valid data from existing tests to give the fuzzer a useful starting point.

**Tech Stack:** Zig 0.15.2, `std.testing.fuzz()` built-in fuzzer

---

### Task 1: Add fuzz test to cif_tokenizer.zig

**Files:**
- Modify: `src/cif_tokenizer.zig` (append after line 462, end of test section)

**Context:** The CIF tokenizer (`Tokenizer` struct) takes a `[]const u8` via `Tokenizer.init(source)` and produces tokens via `tokenizer.next()`. It returns `.eof` when done. It can also return `.err` for malformed input (e.g., unterminated quoted strings). The fuzz test should consume all tokens to EOF without crashing.

**Step 1: Add the fuzz test**

Append this test at the end of `src/cif_tokenizer.zig` (after the `"quoted string with embedded quote char"` test):

```zig
test "fuzz tokenizer" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, input: []const u8) !void {
            var tok = Tokenizer.init(input);
            // Consume all tokens until EOF
            while (true) {
                const token = tok.next();
                if (token == .eof) break;
            }
        }
    }.testOne, .{
        .corpus = &.{
            "data_1ABC",
            "loop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n2 N",
            "'hello world'",
            "\"test value\"",
            ";first line\nsecond line\n;",
            "# comment\ndata_TEST",
            "DATA_test LOOP_",
            "'it''s ok'",
        },
    });
}
```

**Step 2: Run tests to verify it passes**

Run: `zig build test 2>&1 | tail -5`
Expected: All tests pass (fuzz test runs corpus entries as unit tests when not in `--fuzz` mode)

**Step 3: Commit**

```bash
git add src/cif_tokenizer.zig
git commit -m "test: Add fuzz test for CIF tokenizer (#161)"
```

---

### Task 2: Add fuzz test to pdb_parser.zig

**Files:**
- Modify: `src/pdb_parser.zig` (append after line 541, end of test section)

**Context:** `PdbParser` is initialized with `PdbParser.init(allocator)` and parses via `parser.parse(source)`. It returns `AtomInput` on success (must call `.deinit()`) or errors like `NoAtomsFound`, `InvalidCoordinate`, `LineTooShort`. Fuzz test must use `std.testing.allocator` to detect leaks.

**Step 1: Add the fuzz test**

Append this test at the end of `src/pdb_parser.zig` (after the `"PdbParser deuterium filter"` test):

```zig
test "fuzz pdb parser" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, input: []const u8) !void {
            var parser = PdbParser.init(std.testing.allocator);
            const result = parser.parse(input) catch return;
            result.deinit();
        }
    }.testOne, .{
        .corpus = &.{
            "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N\n",
            "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N\nATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  9.13           C\n",
            "HETATM  100  O   HOH A 101       5.000   5.000   5.000  1.00 20.00           O\n",
            "MODEL        1\nATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00 11.68           N\nENDMDL\n",
        },
    });
}
```

**Step 2: Run tests to verify it passes**

Run: `zig build test 2>&1 | tail -5`
Expected: All tests pass

**Step 3: Commit**

```bash
git add src/pdb_parser.zig
git commit -m "test: Add fuzz test for PDB parser (#161)"
```

---

### Task 3: Add fuzz test to mmcif_parser.zig

**Files:**
- Modify: `src/mmcif_parser.zig` (append after line 831, end of test section)

**Context:** `MmcifParser` is initialized with `MmcifParser.init(allocator)` and parses via `parser.parse(source)`. It delegates tokenization to `cif_tokenizer.Tokenizer`. It returns `AtomInput` on success (must call `.deinit()`) or errors like `NoAtomSiteLoop`, `MissingCoordinateField`, `InvalidCoordinate`. Fuzz test must use `std.testing.allocator` to detect leaks.

**Step 1: Add the fuzz test**

Append this test at the end of `src/mmcif_parser.zig` (after the `"startsWithIgnoreCase"` test):

```zig
test "fuzz mmcif parser" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, input: []const u8) !void {
            var parser = MmcifParser.init(std.testing.allocator);
            const result = parser.parse(input) catch return;
            result.deinit();
        }
    }.testOne, .{
        .corpus = &.{
            "data_TEST\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C 10.0 20.0 30.0\n#\n",
            "data_TEST\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C CA ALA 10.0 20.0 30.0\n2 N N ALA 11.0 21.0 31.0\n#\n",
            "data_TEST\nloop_\n_cell.length_a\n_cell.length_b\n10.0 20.0\n#\n",
        },
    });
}
```

**Step 2: Run tests to verify it passes**

Run: `zig build test 2>&1 | tail -5`
Expected: All tests pass

**Step 3: Commit**

```bash
git add src/mmcif_parser.zig
git commit -m "test: Add fuzz test for mmCIF parser (#161)"
```

---

### Task 4: Run fuzzer and verify no crashes

**Files:**
- No modifications (verification only)

**Context:** The `--fuzz` flag runs the fuzzer in continuous mode, generating mutations of the corpus. Let it run for ~60 seconds per target to look for crashes. This is a manual verification step.

**Step 1: Run the fuzzer briefly**

Run: `zig build test --fuzz 2>&1 &` (let it run ~30-60 seconds, then Ctrl+C)

Check output for any crash reports. The fuzzer opens an HTTP server showing coverage.

**Step 2: If crashes found, fix the parser**

If the fuzzer finds a crash:
1. Note the failing input (shown in the crash report)
2. Add a regression test with that input
3. Fix the parser to handle it gracefully (return error, not crash)
4. Re-run the fuzzer to verify the fix

**Step 3: Update CHANGELOG**

Add to the `[Unreleased]` section of `CHANGELOG.md` under `### Added`:

```markdown
- Fuzz tests for CIF tokenizer, PDB parser, and mmCIF parser using Zig's built-in `std.testing.fuzz()` (#161)
```

**Step 4: Commit CHANGELOG**

```bash
git add CHANGELOG.md
git commit -m "docs: Update CHANGELOG with fuzz testing (#161)"
```

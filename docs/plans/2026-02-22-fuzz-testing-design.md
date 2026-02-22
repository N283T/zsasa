# Fuzz Testing for PDB/mmCIF Parsers - Design

## Goal

Add fuzz tests to cif_tokenizer, pdb_parser, and mmcif_parser using Zig 0.15.2's built-in `std.testing.fuzz()` to catch crashes and panics from malformed input.

## Approach

Use Zig's built-in coverage-guided fuzzer. No external dependencies needed.

- Run with `zig build test --fuzz`
- Each parser gets a fuzz test that feeds arbitrary bytes and asserts no crash
- Errors (InvalidCoordinate, NoAtomsFound, etc.) are expected and acceptable
- Seed corpus from existing test strings for better initial coverage

## Targets

1. **cif_tokenizer.zig** - Tokenize arbitrary bytes to EOF without crash
2. **pdb_parser.zig** - Parse arbitrary bytes; errors OK, panics NG
3. **mmcif_parser.zig** - Parse arbitrary bytes; errors OK, panics NG

## Fuzz Test Pattern

```zig
test "fuzz <component>" {
    try std.testing.fuzz({}, struct {
        fn testOne(_: void, input: []const u8) !void {
            // Feed input to parser
            // catch errors (expected for malformed input)
            // only crashes/panics are failures
        }
    }.testOne, .{
        .corpus = &.{ /* valid seed data */ },
    });
}
```

## Non-Goals

- Property-based testing (correct output verification)
- Performance fuzzing
- File I/O fuzzing (parseFile paths)

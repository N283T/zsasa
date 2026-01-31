# Contributing to freesasa-zig

Thank you for your interest in contributing to freesasa-zig!

## Development Setup

### Prerequisites

- **Zig 0.15.2+** - Required for building the project
- **Python 3.11+** - Required for benchmarks and Python bindings
- **FreeSASA C** (optional) - For benchmark comparisons

### Getting the Source

```bash
git clone https://github.com/N283T/freesasa-zig.git
cd freesasa-zig
```

### Building

```bash
# Debug build
zig build

# Release build (recommended for benchmarks)
zig build -Doptimize=ReleaseFast

# Run tests
zig build test

# Verify installation
./zig-out/bin/freesasa_zig --version
```

### Python Bindings (Optional)

To use the Python bindings, first build the shared library:

```bash
zig build -Doptimize=ReleaseFast
```

Then install the Python package in development mode:

```bash
cd python
pip install -e ".[dev]"
pytest tests/ -v
```

For benchmark scripts:

```bash
cd benchmarks/scripts
./run.py --help
```

## Code Style

### Zig

- Run `zig fmt` before committing
- Follow Zig naming conventions (snake_case for functions/variables, PascalCase for types)
- Keep functions focused and under 50 lines when possible
- Use `comptime` for zero-cost abstractions

### Python

- Format with `ruff format .`
- Lint with `ruff check --fix .`
- Use type hints for function signatures

## Pull Request Process

1. **Fork** the repository and create a feature branch
2. **Write tests** for new functionality
3. **Run tests** locally: `zig build test`
4. **Format code**: `zig fmt src/*.zig`
5. **Create PR** with a clear description of changes
6. **Wait for review** and address any feedback

### Commit Messages

Use conventional commit format:

```
<type>: <description>

<optional body>
```

Types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `ci`

Examples:
- `feat: Add Lee-Richards algorithm support`
- `fix: Handle empty input gracefully`
- `perf: Optimize neighbor list construction`

## Reporting Issues

When reporting bugs, please include:

1. **Zig version** (`zig version`)
2. **Operating system** and architecture
3. **Minimal reproduction** - smallest input that triggers the issue
4. **Expected vs actual behavior**
5. **Error messages** or stack traces

For feature requests:
- Describe the use case
- Provide examples of expected API/behavior
- Consider performance implications for large structures

## Architecture Overview

```
src/
├── main.zig           # CLI entry point
├── types.zig          # Core data structures (AtomInput, SasaResult)
├── shrake_rupley.zig  # Shrake-Rupley algorithm (SIMD, multi-threaded)
├── lee_richards.zig   # Lee-Richards algorithm (SIMD, multi-threaded)
├── classifier*.zig    # Atom classifiers (NACCESS, ProtOr, OONS)
├── json_parser.zig    # JSON input parsing
├── pdb_parser.zig     # PDB format parsing
├── mmcif_parser.zig   # mmCIF format parsing
├── c_api.zig          # C ABI for Python bindings
└── thread_pool.zig    # Generic work-stealing thread pool
```

## Performance Considerations

- Always benchmark with `-Doptimize=ReleaseFast`
- Test with structures of varying sizes (small: ~300, medium: ~3000, large: ~50000+ atoms)
- SIMD optimizations should have scalar fallbacks
- Thread pool should scale reasonably from 1 to 10+ threads

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

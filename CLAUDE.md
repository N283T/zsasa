# Project Rules

## Build

**ALWAYS use ReleaseFast for Zig builds:**

```bash
zig build -Doptimize=ReleaseFast
```

Never use `zig build` without `-Doptimize=ReleaseFast` - debug builds are 3-10x slower and will give misleading benchmark results.

## Scripts

Scripts use typer + rich with PEP 723 metadata.

```bash
# Setup (first time)
cd scripts && uv sync

# Run scripts
./scripts/benchmark.py --help
./scripts/validate.py --help

# Lint/Type check
ruff check scripts/*.py
ty check --python scripts/.venv scripts/*.py
```

## Benchmark

Before benchmarking, ensure:
1. Zig binary is built with ReleaseFast
2. FreeSASA C is built with thread support

```bash
zig build -Doptimize=ReleaseFast
./scripts/benchmark.py
```

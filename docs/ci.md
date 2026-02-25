# CI/CD Configuration

## Overview

Continuous integration via GitHub Actions.

## Triggers

```yaml
on:
  push:
    branches: [main]
    paths-ignore:
      - '*.md'
      - 'docs/**'
      - 'plans/**'
      - 'benchmarks/**'
      - 'LICENSE'
  pull_request:
    branches: [main]
    paths-ignore: [same as above]
  workflow_dispatch:  # Manual trigger
```

### Excluded Paths (CI Skip)

CI does not run for changes to these files:

| Path | Reason |
|------|--------|
| `*.md` | Documentation only |
| `docs/**` | Documentation only |
| `plans/**` | Plan files only |
| `benchmarks/**` | Benchmark data and scripts |
| `LICENSE` | License file |

### Manual Trigger

`workflow_dispatch` allows manual CI execution from the GitHub UI.

Actions → CI → Run workflow

## Job Structure

```
┌──────────────┐  ┌─────────────────────────────────┐  ┌─────────────────────────────────┐
│ fmt          │  │ build (parallel x3)             │  │ python (parallel x3)            │
│ (5min)       │  │ Linux / macOS / Windows         │  │ Linux / macOS / Windows         │
│              │  │ (15min)                         │  │ (15min)                         │
└──────┬───────┘  └────────────────┬────────────────┘  └────────────────┬────────────────┘
       │                           │                                    │
       └───────────┬───────────────┘                                    │
                   ▼                                                    │
         ┌─────────────────┐                                            │
         │ validate        │ ◄──────────────────────────────────────────┘
         │ (10min)         │   (python depends only on fmt)
         └─────────────────┘
```

## Job Details

### 1. Format Check (`fmt`)

**Environment:** ubuntu-latest
**Timeout:** 5 minutes
**Runtime:** ~20 seconds

```bash
zig fmt --check src/
```

Recursively checks all `.zig` files under the `src/` directory.

Fails if any formatting violations exist. To fix:

```bash
zig fmt src/
```

### 2. Build (`build`)

**Environment:** Linux, macOS, Windows (parallel execution)
**Timeout:** 15 minutes
**Runtime:**
- Linux: ~45 seconds
- macOS: ~1 minute 10 seconds
- Windows: ~3-4 minutes

| Step | Command | Purpose |
|------|---------|---------|
| Build | `zig build` | Debug build |
| Test | `zig build test` | Run all tests |
| Release | `zig build -Doptimize=ReleaseFast` | Optimized build |
| Binary check | - | Verify binary exists |
| CLI test | `--help`, `--version` | CLI functionality check |

### 3. Validate (`validate`)

**Environment:** ubuntu-latest
**Timeout:** 10 minutes
**Dependencies:** fmt, build (runs after both succeed)
**Runtime:** ~45 seconds

| Step | Description |
|------|-------------|
| Example run | Compute with `examples/1ubq.pdb` |
| JSON validation | Validate output fields (Python) |
| CSV test | Test CSV output |
| Validate mode | Test `--validate` flag |

**JSON Validation:**
- `total_area` field exists
- `atom_areas` field exists
- Atom count > 500 (1UBQ has 602 atoms)

### 4. Python (`python`)

**Environment:** Linux, macOS, Windows (parallel execution)
**Timeout:** 15 minutes
**Dependencies:** fmt
**Runtime:** ~1-2 minutes

| Step | Description |
|------|-------------|
| Build shared library | `zig build -Doptimize=ReleaseFast` |
| Library check | Verify shared library exists |
| Install package | `uv pip install -e ".[dev]"` |
| Run tests | `pytest tests/ -v` |
| Lint | `ruff check` / `ruff format --check` |
| Integration test | Basic functionality from Python |

**Shared Libraries:**
- Linux: `libzsasa.so`
- macOS: `libzsasa.dylib`
- Windows: `zsasa.dll`

**Python Test Coverage:**
- NumPy array input/output
- SR/LR algorithms
- Error handling
- Parameter configuration

## Local CI Reproduction

```bash
# Format check
zig fmt --check src/

# Build and test
zig build
zig build test
zig build -Doptimize=ReleaseFast

# CLI test
./zig-out/bin/zsasa --help
./zig-out/bin/zsasa --version

# Validation test
./zig-out/bin/zsasa calc examples/1ubq.pdb /tmp/output.json
./zig-out/bin/zsasa calc --format=csv examples/1ubq.pdb /tmp/output.csv
./zig-out/bin/zsasa calc --validate examples/1ubq.pdb

# Python test
cd python
uv pip install -e ".[dev]"
pytest tests/ -v
ruff check zsasa/
ruff format --check zsasa/
```

## CI Failure Resolution

### Format check failed

```bash
# Run format locally
zig fmt src/
git add -u && git commit -m "style: Format code"
git push
```

**Note:** Do not use `--amend` + `git push -f` on pushed commits as this rewrites history. Avoid on shared branches.

### Build failed

1. Check error messages
2. Run `zig build` locally to reproduce
3. Fix and push

### Test failed

1. Run `zig build test` locally
2. Identify failing test name
3. Fix the relevant test

### Validate failed

1. Verify execution with `examples/1ubq.pdb`
2. Check output JSON format
3. Confirm atom count > 500

### Python tests failed

1. Verify shared library is built
   ```bash
   ls -la zig-out/lib/
   ```
2. Install Python package
   ```bash
   cd python && uv pip install -e ".[dev]"
   ```
3. Run tests
   ```bash
   uv run --with pytest pytest tests/ -v
   ```

### Python lint failed

```bash
cd python
ruff check zsasa/ --fix
ruff format zsasa/
git add -u && git commit -m "style: Fix Python lint issues"
git push
```

### Timeout

Default timeouts:
- fmt: 5 minutes
- build: 15 minutes
- validate: 10 minutes
- python: 15 minutes

These should not normally be exceeded. If exceeded, investigate for infinite loops.

## Zig Version

Zig version used in CI: **0.15.2**

Managed by the `mlugg/setup-zig@v2` action.

```yaml
- uses: mlugg/setup-zig@v2
  with:
    version: 0.15.2
```

## Security

- Third-party actions:
  - `actions/checkout@v4`
  - `mlugg/setup-zig@v2`
  - `actions/setup-python@v5`
  - `astral-sh/setup-uv@v4`
- Hardcoded secrets: None
- Dangerous permissions: None (uses default permissions)

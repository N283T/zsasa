---
sidebar_position: 1
---

# zsasa

<p align="center">
  <img src="/zsasa/img/logo.svg" alt="zsasa logo" width="420" />
</p>

High-performance Solvent Accessible Surface Area (SASA) calculator in Zig.
The current benchmark suite covers FreeSASA agreement, proteome-scale batch throughput, large single structures, and low-memory MD trajectory analysis.

## Quick Start

### Python

```bash
pip install zsasa
# or
uv add zsasa
```

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
```

### CLI

```bash
# One-line install
curl -fsSL https://raw.githubusercontent.com/N283T/zsasa/main/install.sh | sh

# Or with Nix
nix run github:N283T/zsasa -- calc structure.cif output.json
```

```bash
zsasa calc structure.cif output.json
```

## Documentation

| Section | Description |
|---------|-------------|
| [Getting Started](getting-started.md) | Install zsasa and run first CLI/Python calculations |
| [Choosing CLI or Python](guide/choosing-tool.md) | Pick the right interface for your data and workflow |
| [Batch Processing](guide/batch.md) | Process directories, JSONL output, residue maps, and chain filters |
| [Workflow Files](guide/workflows.md) | Reproducible TOML workflows for calc and batch jobs |
| [Classifiers and CCD](guide/classifiers.mdx) | Atom radius assignment, CCD dictionaries, and custom classifiers |
| [Trajectory Analysis](guide/trajectory.md) | MD trajectory SASA with CLI and Python integrations |
| [CLI Reference](cli/commands.md) | Command syntax and option reference |
| [Python API](python-api/) | Python bindings, core API, and integrations |
| [Benchmarks](benchmarks/) | Performance and accuracy benchmarks |
| [Comparison](comparison.md) | How zsasa compares to FreeSASA, RustSASA, and Lahuta |

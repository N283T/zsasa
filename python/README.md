# zsasa Python Bindings

Python bindings for [zsasa](https://github.com/N283T/zsasa) — a high-performance SASA calculator in Zig.

**[Full Documentation](https://n283t.github.io/zsasa/docs/python-api)**

## Installation

```bash
pip install zsasa
# or
uv add zsasa
```

The x86_64 wheels need a CPU with AVX2 and FMA (x86-64-v3, 2013 and later), and the macOS wheels need macOS 11 or later. On an older CPU, install from source with `pip install --no-binary zsasa zsasa` (needs Zig 0.16).

### Optional Dependencies

```bash
pip install zsasa[gemmi]      # Gemmi integration
pip install zsasa[biopython]  # BioPython integration
pip install zsasa[biotite]    # Biotite integration
pip install zsasa[all]        # All integrations
```

### Native Library

A wheel from PyPI bundles the `libzsasa` shared library and the `zsasa` binary. The library is looked for in the file `ZSASA_LIB` names, the package directory (and `zig-out` of a source checkout), the environment of the running interpreter (`<sys.prefix>/lib`, on Windows `<sys.prefix>\Library\bin`), and `/usr/local/lib` and `/usr/lib`, in that order; the current directory is never searched. Packagers that ship the library and the binary separately build a pure-Python wheel with `ZSASA_NO_BUNDLE=1`. See [How the Native Library Is Found](https://n283t.github.io/zsasa/docs/python-api#library-lookup).

## Quick Start

```python
import numpy as np
from zsasa import calculate_sasa

coords = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
radii = np.array([1.5, 1.5])
result = calculate_sasa(coords, radii)
print(f"Total SASA: {result.total_area:.2f} Å²")
```

```python
# With structure file (gemmi)
from zsasa.integrations.gemmi import calculate_sasa_from_structure
result = calculate_sasa_from_structure("protein.cif")
print(f"Total: {result.total_area:.1f} Å²")
```

## Features

- **Two algorithms**: Shrake-Rupley and Lee-Richards, with bitmask LUT optimization
- **Selectable precision**: `calculate_sasa_batch` takes `precision="f64"` (default) or `"f32"`
- **Multi-threading**: Automatic parallelization
- **Atom classification**: CCD, ProtOr, NACCESS, and OONS classifiers
- **Analysis**: Per-residue aggregation, RSA, polar/nonpolar classification
- **Batch processing**: `process_directory()` for proteome-scale datasets
- **MD trajectory**: [MDTraj](https://github.com/mdtraj/mdtraj) and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis) integrations; use [pyztraj](https://github.com/N283T/ztraj) for direct trajectory-file I/O
- **Integrations**: Gemmi, BioPython, Biotite

See the [full API reference](https://n283t.github.io/zsasa/docs/python-api) for details.

## Development

```bash
cd python
uv run --with pytest pytest tests/ -v    # Tests
ruff format . && ruff check --fix .      # Lint
```

## License

MIT

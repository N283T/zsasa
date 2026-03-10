# zsasa Python Bindings

Python bindings for [zsasa](https://github.com/N283T/zsasa) — a high-performance SASA calculator in Zig.

**[Full Documentation](https://n283t.github.io/zsasa/docs/python-api)**

## Installation

```bash
pip install zsasa
# or
uv add zsasa
```

### Optional Dependencies

```bash
pip install zsasa[gemmi]      # Gemmi integration
pip install zsasa[biopython]  # BioPython integration
pip install zsasa[biotite]    # Biotite integration
pip install zsasa[all]        # All integrations
```

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
- **Selectable precision**: f64 (default) or f32
- **Multi-threading**: Automatic parallelization
- **Atom classification**: NACCESS, ProtOr, OONS classifiers
- **Analysis**: Per-residue aggregation, RSA, polar/nonpolar classification
- **Batch processing**: `process_directory()` for proteome-scale datasets
- **MD trajectory**: Native XTC/DCD readers, [MDTraj](https://github.com/mdtraj/mdtraj) and [MDAnalysis](https://github.com/MDAnalysis/mdanalysis) integration
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

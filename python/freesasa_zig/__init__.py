"""freesasa-zig: Fast SASA calculation using Zig.

This package provides Python bindings for the freesasa-zig library,
a high-performance implementation of Solvent Accessible Surface Area (SASA)
calculation algorithms.

Example:
    >>> import numpy as np
    >>> from freesasa_zig import calculate_sasa
    >>>
    >>> # Single atom
    >>> coords = np.array([[0.0, 0.0, 0.0]])
    >>> radii = np.array([1.5])
    >>> result = calculate_sasa(coords, radii)
    >>> print(f"Total SASA: {result.total_area:.2f} Å²")
"""

from freesasa_zig.core import (
    SasaResult,
    calculate_sasa,
    get_version,
)

__all__ = [
    "calculate_sasa",
    "get_version",
    "SasaResult",
]

__version__ = get_version()

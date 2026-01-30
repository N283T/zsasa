"""Integration modules for external structure parsing libraries.

These modules provide convenience functions for using freesasa-zig
with popular structure parsing libraries like BioPython and gemmi.

Example:
    >>> # Using BioPython integration (requires: pip install freesasa-zig[biopython])
    >>> from freesasa_zig.integrations.biopython import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.pdb")

    >>> # Using gemmi integration (requires: pip install freesasa-zig[gemmi])
    >>> from freesasa_zig.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
"""

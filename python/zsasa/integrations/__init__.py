"""Integration modules for external structure parsing libraries.

These modules provide convenience functions for using zsasa
with popular structure parsing libraries like BioPython, Biotite, and gemmi.

Example:
    >>> # Using BioPython integration (requires: pip install zsasa[biopython])
    >>> from zsasa.integrations.biopython import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.pdb")

    >>> # Using Biotite integration (requires: pip install zsasa[biotite])
    >>> # Also works with AtomWorks (built on Biotite)
    >>> from zsasa.integrations.biotite import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.pdb")

    >>> # Using gemmi integration (requires: pip install zsasa[gemmi])
    >>> from zsasa.integrations.gemmi import calculate_sasa_from_structure
    >>> result = calculate_sasa_from_structure("protein.cif")
"""

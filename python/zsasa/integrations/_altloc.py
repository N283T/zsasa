"""Alternate-location (altloc) selection for the structure integrations.

Structure libraries hand over every alternate conformer of a residue. Counting
all of them puts several atoms on one site, so integrations that read atoms
one by one keep a single conformer per site with :func:`keep_auto_altloc`.

The rules are those of ``--altloc=auto`` in the zsasa command line
(``src/altloc.zig``), so a structure gives the same atoms either way.
"""

from __future__ import annotations

from collections.abc import Hashable, Sequence
from typing import NamedTuple

__all__ = ["SiteAtom", "keep_auto_altloc"]


class SiteAtom(NamedTuple):
    """What altloc selection needs to know about one atom.

    Attributes:
        position: Identifies the residue position: chain, residue number and
            insertion code (and the model, when atoms of several models are
            passed together).
        residue_name: Residue name. With microheterogeneity the alternates of
            one position have different residue names.
        atom_name: Atom name.
        altloc: Altloc ID, or an empty string for an atom without one.
        occupancy: Occupancy.
    """

    position: Hashable
    residue_name: str
    atom_name: str
    altloc: str
    occupancy: float


class _Survivor:
    """The residue that survives at one residue position."""

    __slots__ = ("has_alt_a", "occupancy", "residue_name")

    def __init__(self, residue_name: str, occupancy: float, has_alt_a: bool) -> None:
        self.residue_name = residue_name
        self.occupancy = occupancy
        self.has_alt_a = has_alt_a


class _Site:
    """What the atoms of one atom site amount to."""

    __slots__ = ("best", "best_occupancy", "decider")

    def __init__(self) -> None:
        # "blank" or "A": whichever came first in the file
        self.decider: str | None = None
        # Index of the alternate with the highest occupancy, other than A
        self.best: int | None = None
        self.best_occupancy = 0.0


def keep_auto_altloc(atoms: Sequence[SiteAtom]) -> list[bool]:
    """Decide which atoms survive altloc selection.

    Atoms without an altloc ID are always kept. For the others:

    1. One residue per residue position. When the alternates of a position
       carry more than one residue name (microheterogeneity), the residue of
       the first alternate ``A`` survives, or without an ``A`` the residue
       with the highest occupancy. The occupancy of a residue is that of its
       first alternate.
    2. One alternate per atom site of the surviving residue. Whichever comes
       first, an atom without an altloc ID or an alternate ``A``, decides: the
       first drops every alternate of the site, the second keeps the
       alternates ``A``. A site with neither keeps the alternate with the
       highest occupancy.

    Equal occupancies are a tie, and the residue or alternate that comes first
    wins it.

    Args:
        atoms: The atoms in file order.

    Returns:
        One flag per atom, True for the atoms to keep.
    """
    keep = [True] * len(atoms)
    if not any(atom.altloc for atom in atoms):
        return keep

    # Step 1: the residue that survives at each residue position
    survivors: dict[Hashable, _Survivor] = {}
    residues_with_alternates: set[tuple[Hashable, str]] = set()
    for atom in atoms:
        if not atom.altloc:
            continue
        residue_key = (atom.position, atom.residue_name)
        is_first_of_residue = residue_key not in residues_with_alternates
        residues_with_alternates.add(residue_key)

        survivor = survivors.get(atom.position)
        if survivor is not None and survivor.has_alt_a:
            continue
        if atom.altloc == "A":
            survivors[atom.position] = _Survivor(atom.residue_name, atom.occupancy, True)
        elif is_first_of_residue and (survivor is None or atom.occupancy > survivor.occupancy):
            survivors[atom.position] = _Survivor(atom.residue_name, atom.occupancy, False)

    # Step 2: what each atom site of a surviving residue holds. Atoms without
    # an altloc take part when their residue has alternates.
    sites: dict[tuple[Hashable, str, str], _Site] = {}
    site_of: list[_Site | None] = [None] * len(atoms)
    for i, atom in enumerate(atoms):
        if atom.altloc:
            if survivors[atom.position].residue_name != atom.residue_name:
                keep[i] = False
                continue
        elif (atom.position, atom.residue_name) not in residues_with_alternates:
            continue

        site = sites.setdefault((atom.position, atom.residue_name, atom.atom_name), _Site())
        site_of[i] = site
        if not atom.altloc:
            if site.decider is None:
                site.decider = "blank"
        elif atom.altloc == "A":
            if site.decider is None:
                site.decider = "A"
        elif site.best is None or atom.occupancy > site.best_occupancy:
            site.best = i
            site.best_occupancy = atom.occupancy

    for i, atom in enumerate(atoms):
        site = site_of[i]
        if not atom.altloc or site is None:
            continue
        if site.decider == "blank":
            keep[i] = False
        elif site.decider == "A":
            keep[i] = atom.altloc == "A"
        else:
            keep[i] = site.best == i

    return keep

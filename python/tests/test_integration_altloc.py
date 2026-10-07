"""Tests for the altloc selection shared by the structure integrations.

The expected results are those of ``--altloc=auto`` in the zsasa command line
(the tests of ``src/altloc.zig`` use the same cases).
"""

from __future__ import annotations

from zsasa.integrations._altloc import SiteAtom, keep_auto_altloc


def _atom(
    seq: int,
    residue_name: str,
    atom_name: str,
    altloc: str,
    occupancy: float,
    chain: str = "A",
) -> SiteAtom:
    return SiteAtom((chain, seq, " "), residue_name, atom_name, altloc, occupancy)


def test_atoms_without_altloc_are_all_kept() -> None:
    atoms = [_atom(1, "ALA", "N", "", 1.0), _atom(1, "ALA", "CA", "", 1.0)]
    assert keep_auto_altloc(atoms) == [True, True]
    assert keep_auto_altloc([]) == []


def test_blank_then_a_then_highest_occupancy() -> None:
    atoms = [
        # A wins whatever its occupancy and position
        _atom(1, "ALA", "CA", "B", 0.9),
        _atom(1, "ALA", "CA", "A", 0.1),
        # Without A the highest occupancy wins
        _atom(1, "ALA", "CB", "B", 0.3),
        _atom(1, "ALA", "CB", "C", 0.7),
        # An atom without an altloc drops the alternates of its site
        _atom(1, "ALA", "N", "", 1.0),
        _atom(1, "ALA", "N", "B", 0.5),
        # A tie goes to the first alternate
        _atom(1, "ALA", "O", "C", 0.5),
        _atom(1, "ALA", "O", "B", 0.5),
    ]
    assert keep_auto_altloc(atoms) == [False, True, False, True, True, False, True, False]


def test_one_alternate_per_site_when_occupancies_tie() -> None:
    atoms = [
        _atom(1, "ALA", "CA", "B", 0.5),
        _atom(1, "ALA", "CA", "C", 0.5),
        _atom(1, "ALA", "CB", "B", 0.0),
        _atom(1, "ALA", "CB", "C", 0.0),
        _atom(1, "ALA", "CB", "D", 0.0),
    ]
    assert keep_auto_altloc(atoms) == [True, False, True, False, False]


def test_sites_are_per_residue_position() -> None:
    atoms = [
        _atom(1, "ALA", "CA", "A", 0.6),
        _atom(1, "ALA", "CA", "B", 0.4),
        # A later residue that only has B keeps it
        _atom(2, "GLY", "CA", "B", 0.5),
        # The same number in another chain is another position
        _atom(1, "ALA", "CA", "B", 0.5, chain="B"),
    ]
    assert keep_auto_altloc(atoms) == [True, False, True, True]


def test_microheterogeneity_keeps_one_residue() -> None:
    # Position 2 is PRO as altloc A and SER as altloc B, with a shared N that
    # has no altloc. Position 3 is ILE as B and VAL as C, with no alternate A.
    atoms = [
        _atom(2, "PRO", "N", "", 1.0),
        _atom(2, "PRO", "CA", "A", 0.4),
        _atom(2, "PRO", "CD", "A", 0.4),
        _atom(2, "SER", "CA", "B", 0.6),
        _atom(2, "SER", "OG", "B", 0.6),
        _atom(3, "ILE", "CA", "B", 0.3),
        _atom(3, "ILE", "CD1", "B", 0.3),
        _atom(3, "VAL", "CA", "C", 0.7),
        _atom(3, "VAL", "CG1", "C", 0.7),
    ]
    # A, then the highest occupancy
    assert keep_auto_altloc(atoms) == [True, True, True, False, False, False, False, True, True]


def test_microheterogeneity_compares_first_alternates_and_keeps_first_on_tie() -> None:
    tie = [
        _atom(1, "LEU", "CA", "B", 0.5),
        _atom(1, "LEU", "CB", "B", 0.5),
        _atom(1, "ILE", "CA", "C", 0.5),
        _atom(1, "ILE", "CB", "C", 0.5),
    ]
    assert keep_auto_altloc(tie) == [True, True, False, False]

    # Only the first alternate of a residue counts: 0.4 against 0.5
    first_alternate = [
        _atom(1, "LEU", "CA", "B", 0.4),
        _atom(1, "LEU", "CB", "B", 0.9),
        _atom(1, "ILE", "CA", "C", 0.5),
        _atom(1, "ILE", "CB", "C", 0.1),
    ]
    assert keep_auto_altloc(first_alternate) == [False, False, True, True]


def test_microheterogeneity_with_interleaved_atoms() -> None:
    # The alternate A that decides for PRO comes last
    atoms = [
        _atom(2, "SER", "CA", "B", 0.6),
        _atom(2, "PRO", "CA", "C", 0.2),
        _atom(2, "SER", "CB", "B", 0.6),
        _atom(2, "PRO", "CB", "C", 0.2),
        _atom(2, "SER", "OG", "B", 0.6),
        _atom(2, "PRO", "CG", "A", 0.2),
    ]
    assert keep_auto_altloc(atoms) == [False, True, False, True, False, True]


def test_residues_without_alternates_are_left_alone() -> None:
    # Two residues without altlocs share chain and number, next to a residue
    # with alternates at the same position
    atoms = [
        _atom(1, "ALA", "CA", "", 1.0),
        _atom(1, "GLY", "CA", "", 1.0),
        _atom(1, "SER", "CA", "A", 0.5),
        _atom(1, "PRO", "CA", "B", 0.5),
    ]
    assert keep_auto_altloc(atoms) == [True, True, True, False]

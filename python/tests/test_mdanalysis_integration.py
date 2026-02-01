"""Tests for MDAnalysis integration."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    import MDAnalysis as mda

MDAnalysis = pytest.importorskip("MDAnalysis")

from zsasa.mdanalysis import (  # noqa: E402
    SASAAnalysis,
    _get_element,
    _get_radii_from_atomgroup,
    _get_radius,
    compute_sasa,
)


def make_universe(
    n_atoms: int = 10,
    n_residues: int = 2,
    n_frames: int = 3,
) -> mda.Universe:
    """Create a simple test universe."""
    from MDAnalysis.coordinates.memory import MemoryReader
    from MDAnalysis.core.topologyattrs import (
        Atomnames,
        Atomtypes,
        Elements,
        Resids,
        Resnames,
        Resnums,
    )

    # Create empty universe
    u = MDAnalysis.Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_residues,
        n_segments=1,
        atom_resindex=np.repeat(np.arange(n_residues), n_atoms // n_residues),
        trajectory=True,
    )

    # Add topology attributes
    u.add_TopologyAttr(Atomnames(["CA"] * n_atoms))
    u.add_TopologyAttr(Atomtypes(["C"] * n_atoms))
    u.add_TopologyAttr(Elements(["C"] * n_atoms))
    u.add_TopologyAttr(Resnames(["ALA"] * n_residues))
    u.add_TopologyAttr(Resids(list(range(1, n_residues + 1))))
    u.add_TopologyAttr(Resnums(list(range(1, n_residues + 1))))

    # Create coordinates - atoms spread out to avoid overlap
    coords = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    for frame in range(n_frames):
        for i in range(n_atoms):
            # Spread atoms along x-axis with slight variation per frame
            coords[frame, i] = [i * 5.0 + frame * 0.1, 0.0, 0.0]

    # Attach trajectory
    u.trajectory = MemoryReader(coords)

    return u


class TestElementInference:
    """Test element inference from atoms."""

    def test_get_element_from_element_attr(self) -> None:
        """Test element extraction from element attribute."""
        u = make_universe(n_atoms=1, n_residues=1, n_frames=1)
        atom = u.atoms[0]
        element = _get_element(atom)
        assert element == "C"

    def test_get_radius_carbon(self) -> None:
        """Test radius for carbon atom."""
        u = make_universe(n_atoms=1, n_residues=1, n_frames=1)
        atom = u.atoms[0]
        radius = _get_radius(atom)
        assert radius == pytest.approx(1.70, rel=0.01)

    def test_get_radii_from_atomgroup(self) -> None:
        """Test radii extraction from atomgroup."""
        u = make_universe(n_atoms=5, n_residues=1, n_frames=1)
        radii = _get_radii_from_atomgroup(u.atoms)
        assert len(radii) == 5
        assert all(r == pytest.approx(1.70, rel=0.01) for r in radii)


class TestSASAAnalysis:
    """Test SASAAnalysis class."""

    @pytest.fixture
    def universe(self) -> mda.Universe:
        """Create test universe."""
        return make_universe(n_atoms=10, n_residues=2, n_frames=3)

    def test_init_with_universe(self, universe: mda.Universe) -> None:
        """Test initialization with Universe."""
        analysis = SASAAnalysis(universe)
        assert analysis.atomgroup.n_atoms == 10

    def test_init_with_atomgroup(self, universe: mda.Universe) -> None:
        """Test initialization with AtomGroup."""
        analysis = SASAAnalysis(universe.atoms)
        assert analysis.atomgroup.n_atoms == 10

    def test_init_with_selection(self, universe: mda.Universe) -> None:
        """Test initialization with selection string."""
        analysis = SASAAnalysis(universe, select="resid 1")
        assert analysis.atomgroup.n_atoms == 5

    def test_run_basic(self, universe: mda.Universe) -> None:
        """Test basic run."""
        analysis = SASAAnalysis(universe)
        analysis.run()

        assert analysis.n_frames == 3
        assert analysis.results.atom_area.shape == (3, 10)
        assert analysis.results.residue_area.shape == (3, 2)
        assert analysis.results.total_area.shape == (3,)
        assert analysis.results.mean_total_area > 0

    def test_run_with_frame_range(self, universe: mda.Universe) -> None:
        """Test run with start/stop/step."""
        analysis = SASAAnalysis(universe)
        analysis.run(start=0, stop=2, step=1)

        assert analysis.n_frames == 2

    def test_run_with_step(self, universe: mda.Universe) -> None:
        """Test run with step > 1."""
        analysis = SASAAnalysis(universe)
        analysis.run(start=0, stop=3, step=2)

        assert analysis.n_frames == 2  # frames 0 and 2

    def test_results_dict_access(self, universe: mda.Universe) -> None:
        """Test dict-like access to results."""
        analysis = SASAAnalysis(universe)
        analysis.run()

        # Dict-like access for mdsasa-bolt compatibility
        assert analysis.results["total_area"] is not None
        assert analysis.results["residue_area"] is not None
        assert analysis.results["atom_area"] is not None

    def test_all_sasa_positive(self, universe: mda.Universe) -> None:
        """Test that all SASA values are positive."""
        analysis = SASAAnalysis(universe)
        analysis.run()

        assert np.all(analysis.results.atom_area >= 0)
        assert np.all(analysis.results.residue_area >= 0)
        assert np.all(analysis.results.total_area >= 0)

    def test_residue_aggregation(self, universe: mda.Universe) -> None:
        """Test that residue SASA matches atom SASA sum."""
        analysis = SASAAnalysis(universe)
        analysis.run()

        # Each residue has 5 atoms
        for frame in range(analysis.n_frames):
            res0_atoms = analysis.results.atom_area[frame, :5].sum()
            res1_atoms = analysis.results.atom_area[frame, 5:].sum()

            assert analysis.results.residue_area[frame, 0] == pytest.approx(
                res0_atoms, rel=0.01
            )
            assert analysis.results.residue_area[frame, 1] == pytest.approx(
                res1_atoms, rel=0.01
            )

    def test_total_matches_sum(self, universe: mda.Universe) -> None:
        """Test that total SASA matches atom SASA sum."""
        analysis = SASAAnalysis(universe)
        analysis.run()

        for frame in range(analysis.n_frames):
            expected = analysis.results.atom_area[frame].sum()
            assert analysis.results.total_area[frame] == pytest.approx(
                expected, rel=0.01
            )

    def test_method_chaining(self, universe: mda.Universe) -> None:
        """Test that run() returns self for chaining."""
        analysis = SASAAnalysis(universe)
        result = analysis.run()
        assert result is analysis


class TestSASAAnalysisAlgorithms:
    """Test different algorithms."""

    @pytest.fixture
    def universe(self) -> mda.Universe:
        """Create test universe."""
        return make_universe(n_atoms=10, n_residues=2, n_frames=2)

    def test_sr_algorithm(self, universe: mda.Universe) -> None:
        """Test Shrake-Rupley algorithm."""
        analysis = SASAAnalysis(universe)
        analysis.run(algorithm="sr", n_points=500)

        assert analysis.results.atom_area is not None
        assert np.all(analysis.results.atom_area > 0)

    def test_lr_algorithm(self, universe: mda.Universe) -> None:
        """Test Lee-Richards algorithm."""
        analysis = SASAAnalysis(universe)
        analysis.run(algorithm="lr", n_slices=20)

        assert analysis.results.atom_area is not None
        assert np.all(analysis.results.atom_area > 0)

    def test_sr_lr_similar(self, universe: mda.Universe) -> None:
        """Test that SR and LR give similar results."""
        analysis_sr = SASAAnalysis(universe)
        analysis_sr.run(algorithm="sr", n_points=500)

        analysis_lr = SASAAnalysis(universe)
        analysis_lr.run(algorithm="lr", n_slices=50)

        np.testing.assert_allclose(
            analysis_sr.results.total_area,
            analysis_lr.results.total_area,
            rtol=0.05,
        )

    def test_invalid_algorithm(self, universe: mda.Universe) -> None:
        """Test that invalid algorithm raises error."""
        analysis = SASAAnalysis(universe)
        with pytest.raises(ValueError, match="algorithm"):
            analysis.run(algorithm="invalid")


class TestSASAAnalysisThreading:
    """Test threading behavior."""

    @pytest.fixture
    def universe(self) -> mda.Universe:
        """Create test universe."""
        return make_universe(n_atoms=20, n_residues=4, n_frames=5)

    def test_threading_consistency(self, universe: mda.Universe) -> None:
        """Test that threading doesn't affect results."""
        analysis_1 = SASAAnalysis(universe)
        analysis_1.run(n_threads=1)

        analysis_4 = SASAAnalysis(universe)
        analysis_4.run(n_threads=4)

        np.testing.assert_array_almost_equal(
            analysis_1.results.atom_area,
            analysis_4.results.atom_area,
            decimal=4,
        )


class TestComputeSasaFunction:
    """Test convenience function."""

    @pytest.fixture
    def universe(self) -> mda.Universe:
        """Create test universe."""
        return make_universe(n_atoms=10, n_residues=2, n_frames=3)

    def test_atom_mode(self, universe: mda.Universe) -> None:
        """Test atom mode output."""
        sasa = compute_sasa(universe, mode="atom")
        assert sasa.shape == (3, 10)

    def test_residue_mode(self, universe: mda.Universe) -> None:
        """Test residue mode output."""
        sasa = compute_sasa(universe, mode="residue")
        assert sasa.shape == (3, 2)

    def test_total_mode(self, universe: mda.Universe) -> None:
        """Test total mode output."""
        sasa = compute_sasa(universe, mode="total")
        assert sasa.shape == (3,)

    def test_with_selection(self, universe: mda.Universe) -> None:
        """Test with atom selection."""
        sasa = compute_sasa(universe, select="resid 1", mode="atom")
        assert sasa.shape == (3, 5)  # Only first residue

    def test_with_frame_range(self, universe: mda.Universe) -> None:
        """Test with frame range."""
        sasa = compute_sasa(universe, start=1, stop=3, mode="total")
        assert sasa.shape == (2,)


class TestUnitsAndValues:
    """Test that output units and values are correct."""

    def test_isolated_atom_sasa(self) -> None:
        """Test SASA of isolated atom."""
        from MDAnalysis.coordinates.memory import MemoryReader
        from MDAnalysis.core.topologyattrs import (
            Atomnames,
            Atomtypes,
            Elements,
            Resids,
            Resnames,
            Resnums,
        )

        # Create universe with single atom
        u = MDAnalysis.Universe.empty(
            n_atoms=1,
            n_residues=1,
            n_segments=1,
            atom_resindex=[0],
            trajectory=True,
        )
        u.add_TopologyAttr(Atomnames(["CA"]))
        u.add_TopologyAttr(Atomtypes(["C"]))
        u.add_TopologyAttr(Elements(["C"]))
        u.add_TopologyAttr(Resnames(["ALA"]))
        u.add_TopologyAttr(Resids([1]))
        u.add_TopologyAttr(Resnums([1]))

        coords = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
        u.trajectory = MemoryReader(coords)

        analysis = SASAAnalysis(u)
        analysis.run(probe_radius=1.4)

        # Carbon radius = 1.70 Å, probe = 1.4 Å
        # Effective radius = 3.1 Å
        # Area = 4 * pi * 3.1^2 = 120.76 Å²
        expected = 4 * np.pi * (1.70 + 1.4) ** 2
        assert analysis.results.total_area[0] == pytest.approx(expected, rel=0.05)

    def test_output_in_angstrom_squared(self) -> None:
        """Test that output is in Å² (not nm²)."""
        u = make_universe(n_atoms=1, n_residues=1, n_frames=1)
        analysis = SASAAnalysis(u)
        analysis.run()

        # SASA should be ~100-150 Å² for single carbon
        # If it were nm², it would be ~1-1.5
        assert analysis.results.total_area[0] > 50
        assert analysis.results.total_area[0] < 200

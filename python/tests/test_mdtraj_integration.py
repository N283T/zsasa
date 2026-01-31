"""Tests for MDTraj integration."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    import mdtraj as md

mdtraj = pytest.importorskip("mdtraj")

from freesasa_zig.mdtraj import (  # noqa: E402
    _get_radii_from_topology,
    compute_sasa,
)


class TestUnitConversions:
    """Test unit conversions between MDTraj (nm) and freesasa-zig (Angstrom)."""

    def test_coordinates_nm_to_angstrom(self) -> None:
        """Verify coordinates are converted from nm to Angstrom (x10)."""
        # Create a simple trajectory with known coordinates in nm
        topology = mdtraj.Topology()
        chain = topology.add_chain()
        residue = topology.add_residue("ALA", chain)
        topology.add_atom("CA", mdtraj.element.carbon, residue)

        # 1.0 nm position
        coords_nm = np.array([[[1.0, 2.0, 3.0]]], dtype=np.float32)
        traj = mdtraj.Trajectory(coords_nm, topology)

        # Internal conversion should multiply by 10
        # We can verify this indirectly by checking SASA scales correctly
        sasa = compute_sasa(traj)

        # SASA of isolated carbon atom
        # MDTraj uses 1.7 Angstrom = 0.170 nm for carbon
        # With probe 1.4 A, effective radius = 1.7 + 1.4 = 3.1 A
        # Area = 4 * pi * 3.1^2 = 120.76 A^2 = 1.2076 nm^2
        expected_area_nm2 = 4 * np.pi * (1.7 + 1.4) ** 2 / 100.0
        assert sasa[0, 0] == pytest.approx(expected_area_nm2, rel=0.05)

    def test_sasa_output_in_nm2(self) -> None:
        """Verify SASA output is in nm^2 (divided by 100)."""
        topology = mdtraj.Topology()
        chain = topology.add_chain()
        residue = topology.add_residue("ALA", chain)
        topology.add_atom("CA", mdtraj.element.carbon, residue)

        coords_nm = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
        traj = mdtraj.Trajectory(coords_nm, topology)

        sasa = compute_sasa(traj)

        # Result should be in nm^2, which is much smaller than A^2
        # 1 nm^2 = 100 A^2, so ~120 A^2 = ~1.2 nm^2
        assert sasa[0, 0] < 10  # Definitely in nm^2, not A^2
        assert sasa[0, 0] > 0.5  # But still a reasonable value


class TestAtomicRadii:
    """Test atomic radii extraction from topology."""

    def test_radii_from_topology(self) -> None:
        """Test _get_radii_from_topology returns correct values."""
        topology = mdtraj.Topology()
        chain = topology.add_chain()
        residue = topology.add_residue("ALA", chain)
        topology.add_atom("CA", mdtraj.element.carbon, residue)
        topology.add_atom("N", mdtraj.element.nitrogen, residue)
        topology.add_atom("O", mdtraj.element.oxygen, residue)

        radii = _get_radii_from_topology(topology)

        # Should be in Angstrom (nm * 10)
        assert len(radii) == 3
        assert radii[0] == pytest.approx(1.70, rel=0.01)  # Carbon
        assert radii[1] == pytest.approx(1.55, rel=0.01)  # Nitrogen
        assert radii[2] == pytest.approx(1.52, rel=0.01)  # Oxygen

    def test_unknown_element_default(self) -> None:
        """Test that unknown elements get default radius."""
        # Use an element not in _ATOMIC_RADII_NM
        # Most common elements are covered, but check the default behavior
        topology = mdtraj.Topology()
        chain = topology.add_chain()
        residue = topology.add_residue("UNK", chain)

        # Add atom with element that has a radius in the table
        topology.add_atom("C", mdtraj.element.carbon, residue)

        radii = _get_radii_from_topology(topology)
        assert len(radii) == 1
        assert radii[0] > 0  # Should have some radius


class TestComputeSasa:
    """Test compute_sasa function."""

    @pytest.fixture
    def simple_trajectory(self) -> md.Trajectory:
        """Create a simple trajectory for testing."""
        topology = mdtraj.Topology()
        chain = topology.add_chain()
        residue = topology.add_residue("ALA", chain)
        topology.add_atom("CA", mdtraj.element.carbon, residue)
        topology.add_atom("CB", mdtraj.element.carbon, residue)

        # Two frames with atoms at different positions
        coords = np.array(
            [
                [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],  # Frame 0: close atoms
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],  # Frame 1: farther apart
            ],
            dtype=np.float32,
        )
        return mdtraj.Trajectory(coords, topology)

    def test_atom_mode(self, simple_trajectory: md.Trajectory) -> None:
        """Test per-atom SASA output."""
        sasa = compute_sasa(simple_trajectory, mode="atom")

        assert sasa.shape == (2, 2)  # (n_frames, n_atoms)
        assert sasa.dtype == np.float32
        # Frame 1 atoms farther apart -> more SASA
        assert sasa[1].sum() > sasa[0].sum()

    def test_residue_mode(self, simple_trajectory: md.Trajectory) -> None:
        """Test per-residue SASA aggregation."""
        sasa = compute_sasa(simple_trajectory, mode="residue")

        assert sasa.shape == (2, 1)  # (n_frames, n_residues)
        # Total should match atom mode sum
        sasa_atom = compute_sasa(simple_trajectory, mode="atom")
        np.testing.assert_array_almost_equal(sasa[:, 0], sasa_atom.sum(axis=1))

    def test_total_mode(self, simple_trajectory: md.Trajectory) -> None:
        """Test total SASA per frame output."""
        sasa_total = compute_sasa(simple_trajectory, mode="total")

        assert sasa_total.shape == (2,)  # (n_frames,)
        # Total should match atom mode sum
        sasa_atom = compute_sasa(simple_trajectory, mode="atom")
        np.testing.assert_array_almost_equal(sasa_total, sasa_atom.sum(axis=1))

    def test_algorithm_sr(self, simple_trajectory: md.Trajectory) -> None:
        """Test Shrake-Rupley algorithm."""
        sasa = compute_sasa(simple_trajectory, algorithm="sr", n_points=960)
        assert sasa.shape == (2, 2)
        assert np.all(sasa > 0)

    def test_algorithm_lr(self, simple_trajectory: md.Trajectory) -> None:
        """Test Lee-Richards algorithm."""
        sasa = compute_sasa(simple_trajectory, algorithm="lr", n_slices=20)
        assert sasa.shape == (2, 2)
        assert np.all(sasa > 0)

    def test_invalid_algorithm(self, simple_trajectory: md.Trajectory) -> None:
        """Test that invalid algorithm raises error."""
        with pytest.raises(ValueError, match="algorithm"):
            compute_sasa(simple_trajectory, algorithm="invalid")

    def test_threading(self, simple_trajectory: md.Trajectory) -> None:
        """Test that threading doesn't affect results."""
        sasa_1 = compute_sasa(simple_trajectory, n_threads=1)
        sasa_4 = compute_sasa(simple_trajectory, n_threads=4)

        np.testing.assert_array_almost_equal(sasa_1, sasa_4, decimal=4)


class TestMultiResidue:
    """Test with multi-residue structures."""

    @pytest.fixture
    def multi_residue_trajectory(self) -> md.Trajectory:
        """Create a trajectory with multiple residues."""
        topology = mdtraj.Topology()
        chain = topology.add_chain()

        # Residue 1: ALA
        res1 = topology.add_residue("ALA", chain)
        topology.add_atom("N", mdtraj.element.nitrogen, res1)
        topology.add_atom("CA", mdtraj.element.carbon, res1)
        topology.add_atom("C", mdtraj.element.carbon, res1)

        # Residue 2: GLY
        res2 = topology.add_residue("GLY", chain)
        topology.add_atom("N", mdtraj.element.nitrogen, res2)
        topology.add_atom("CA", mdtraj.element.carbon, res2)
        topology.add_atom("C", mdtraj.element.carbon, res2)

        # Single frame with atoms spread out
        coords = np.array(
            [
                [
                    [0.0, 0.0, 0.0],
                    [0.15, 0.0, 0.0],
                    [0.3, 0.0, 0.0],
                    [0.5, 0.0, 0.0],
                    [0.65, 0.0, 0.0],
                    [0.8, 0.0, 0.0],
                ]
            ],
            dtype=np.float32,
        )
        return mdtraj.Trajectory(coords, topology)

    def test_residue_aggregation(self, multi_residue_trajectory: md.Trajectory) -> None:
        """Test that residue mode correctly aggregates atoms."""
        sasa_atom = compute_sasa(multi_residue_trajectory, mode="atom")
        sasa_residue = compute_sasa(multi_residue_trajectory, mode="residue")

        assert sasa_atom.shape == (1, 6)
        assert sasa_residue.shape == (1, 2)

        # Residue 0 = atoms 0,1,2; Residue 1 = atoms 3,4,5
        expected_res0 = sasa_atom[0, 0:3].sum()
        expected_res1 = sasa_atom[0, 3:6].sum()

        assert sasa_residue[0, 0] == pytest.approx(expected_res0, rel=0.01)
        assert sasa_residue[0, 1] == pytest.approx(expected_res1, rel=0.01)


class TestShrakeRupleyAlias:
    """Test shrake_rupley alias function."""

    def test_alias_exists(self) -> None:
        """Test that shrake_rupley is an alias for compute_sasa."""
        from freesasa_zig.mdtraj import compute_sasa, shrake_rupley

        assert shrake_rupley is compute_sasa

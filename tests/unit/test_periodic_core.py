"""Unit tests for the unified periodic-geometry layer (core/periodic.py).

These pin the *primitives*: an exact toroidal distance transform, torus-closing
labelling, periodic SASA, and whole-molecule reconstruction. Every check is
against an independent brute-force reference computed on a torus by hand.
"""

from __future__ import annotations

import numpy as np
import pytest

from lat.core.errors import LatError
from lat.core.pbc import Box
from lat.core.periodic import (
    bond_graph_consistent,
    component_percolating,
    grid_axes,
    label_periodic,
    periodic_distance_transform,
    periodic_nearest,
    periodic_neighbors,
    residue_com,
    sasa_shrake_rupley,
    selection_is_whole,
    voxel_occupancy,
    whole_positions,
)
from lat.core.topology import Topology


def _topo(n, mols=None, bonds=None, types=None, n_atoms=None):
    return Topology(
        n_atoms=n_atoms if n_atoms is not None else n,
        atom_id=np.arange(1, (n_atoms if n_atoms is not None else n) + 1),
        types=np.zeros(n, dtype=int) if types is None else np.asarray(types, int),
        box=Box.orthorhombic([10.0, 10.0, 10.0]),
        mols=None if mols is None else np.asarray(mols),
        bonds=None if bonds is None else np.asarray(bonds, int).reshape(-1, 2),
    )


# ---------------------------------------------------------------------------
# distance transform
# ---------------------------------------------------------------------------
def _brute_toroidal_edt(mask, spacing):
    """O(N^2) toroidal EDT by explicit minimisation over all voxel pairs."""
    idx = np.array(np.nonzero(mask)).T
    dims = np.array(mask.shape)
    out = np.full(mask.shape, np.inf)
    for p in np.ndindex(*mask.shape):
        d = np.abs(idx - np.array(p))
        d = np.minimum(d, dims - d)          # circular per axis
        dist = np.sqrt(((d * spacing) ** 2).sum(axis=1))
        out[p] = dist.min()
    return out


def test_periodic_edt_matches_brute_force():
    rng = np.random.default_rng(7)
    mask = rng.random((9, 8, 7)) < 0.15
    spacing = (1.3, 0.7, 1.1)
    got = periodic_distance_transform(mask, spacing)
    want = _brute_toroidal_edt(mask, spacing)
    assert np.allclose(got, want, atol=1e-9)


def test_periodic_edt_is_periodic_not_euclidean():
    """A single occupied corner voxel must be closer to the *far* corner on the
    torus than the plain Euclidean distance suggests."""
    mask = np.zeros((20, 20, 20), dtype=bool)
    mask[0, 0, 0] = True
    d = periodic_distance_transform(mask, (1.0, 1.0, 1.0))
    assert d[0, 0, 0] == 0.0
    assert np.isclose(d[19, 19, 19], np.sqrt(3.0))          # wraps one voxel
    # non-periodic EDT would give |19,19,19| = 32.9
    assert d[19, 19, 19] < 5.0
    assert np.isclose(d[10, 0, 0], 10.0)                    # equidistant both ways


def test_periodic_edt_rejects_empty_mask():
    with pytest.raises(LatError, match="empty"):
        periodic_distance_transform(np.zeros((4, 4, 4), dtype=bool), (1.0, 1.0, 1.0))


# ---------------------------------------------------------------------------
# labelling
# ---------------------------------------------------------------------------
def test_label_periodic_merges_across_faces():
    """Two blobs sitting on opposite faces are one component on the torus."""
    m = np.zeros((12, 12, 12), dtype=bool)
    m[0, :, :] = True
    m[-1, :, :] = True
    labels, n = label_periodic(m, connectivity=6)
    assert n == 1


def test_label_periodic_keeps_disconnected_blobs_apart():
    m = np.zeros((12, 12, 12), dtype=bool)
    m[2, 2, 2] = True
    m[8, 8, 8] = True
    labels, n = label_periodic(m, connectivity=6)
    assert n == 2
    assert labels[2, 2, 2] != labels[8, 8, 8]


def test_label_periodic_connectivity_6_vs_26():
    """Face-touching cubes are one component under 26 but two under 6 -- but
    here the diagonal bridge is only joined through a corner."""
    m = np.zeros((8, 8, 8), dtype=bool)
    m[1, 1, 1] = True
    m[2, 2, 2] = True
    assert label_periodic(m, connectivity=6)[1] == 2
    assert label_periodic(m, connectivity=26)[1] == 1


def test_component_percolating():
    m = np.zeros((10, 10, 10), dtype=bool)
    m[:, 5, 5] = True                     # spans x only
    labels, n = label_periodic(m, 6)
    assert n == 1
    assert component_percolating(labels, 1)
    m2 = np.zeros((10, 10, 10), dtype=bool)
    m2[4:6, 4:6, 4:6] = True
    labels2, _ = label_periodic(m2, 6)
    assert not component_percolating(labels2, 1)


# ---------------------------------------------------------------------------
# voxel occupancy
# ---------------------------------------------------------------------------
def test_voxel_occupancy_is_periodic_and_complete():
    """Voxel occupancy must equal a brute-force test of every grid centre, and
    an atom sitting on a face must mark voxels on the opposite side too."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    axes, spacing, _ = grid_axes(box, 0.25)
    pos = np.array([[0.0, 5.0, 5.0], [7.3, 2.1, 9.9]])
    radii = np.array([1.0, 0.6])
    occ = voxel_occupancy(pos, box, axes, radii)
    assert occ[0].any() and occ[-1].any()        # wraps in x
    # brute force: mark every grid centre within the radius (minimum image)
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    want = np.zeros(occ.shape, dtype=bool)
    for p, r in zip(pos, radii):
        want |= np.linalg.norm(box.minimum_image(grid - p), axis=-1) <= r
    assert np.array_equal(occ, want)


def test_voxel_occupancy_volume_is_within_sampling_bias():
    """A single sphere's voxel volume is the analytic volume plus the usual
    half-voxel centre-sampling bias (declared, not a defect)."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    axes, spacing, _ = grid_axes(box, 0.25)
    occ = voxel_occupancy(np.array([[0.0, 5.0, 5.0]]), box, axes, np.array([1.0]))
    want = 4.0 / 3.0 * np.pi
    got = occ.sum() * np.prod(spacing)
    assert abs(got - want) / want < 0.06


def test_voxel_grid_rejects_triclinic():
    box = Box.from_lammps_data([0, 0, 0], [10, 10, 10], [1.0, 0.0, 0.0])
    with pytest.raises(LatError, match="triclinic"):
        grid_axes(box, 0.5)


# ---------------------------------------------------------------------------
# molecule level
# ---------------------------------------------------------------------------
def test_whole_positions_reconstructs_a_chain_crossing_the_face():
    """A 3-atom chain with 4 Å bonds laid across the x face of a 10 Å box: each
    bond is restored to its 4 Å minimum-image length, so the chain is 8 Å long
    even though its wrapped coordinates span only 6 Å."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[8.5, 5.0, 5.0], [2.5, 5.0, 5.0], [6.5, 5.0, 5.0]])
    topo = _topo(3, mols=[1, 1, 1], bonds=[[0, 1], [1, 2]])
    got = whole_positions(pos, box, topo)
    assert np.allclose(got[:, 0], [8.5, 12.5, 16.5])
    assert np.isclose(np.ptp(got[:, 0]), 8.0)
    # a chain longer than half the box must still work (the first-atom MIC trick
    # would fold the far end back)
    pos2 = np.array([[0.0, 0, 0], [6.0, 0, 0], [2.0, 0, 0]])
    topo2 = _topo(3, mols=[1, 1, 1], bonds=[[0, 1], [1, 2]])
    got2 = whole_positions(pos2, box, topo2)
    assert np.isclose(np.ptp(got2[:, 0]), 8.0)          # 0 -> -4 -> -8


def test_whole_positions_subset_only_touches_selected_molecule():
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[9.5, 0, 0], [0.5, 0, 0], [1.5, 0, 0], [2.5, 0, 0]])
    topo = _topo(4, mols=[1, 1, 2, 2], bonds=[[0, 1], [2, 3]])
    out = whole_positions(pos, box, topo, indices=[0, 1])
    assert out.shape == (2, 3)
    assert np.isclose(abs(out[0, 0] - out[1, 0]), 1.0)
    full = whole_positions(pos, box, topo)
    assert np.allclose(full[2:], pos[2:])       # untouched molecule unchanged


def test_selection_is_whole():
    topo = _topo(4, mols=[1, 1, 2, 2], bonds=[[0, 1], [2, 3]])
    assert selection_is_whole(topo, [0, 1]) is True
    assert selection_is_whole(topo, [0, 1, 2, 3]) is True
    assert selection_is_whole(topo, [0]) is False
    assert selection_is_whole(topo, [2, 3, 0]) is False
    assert selection_is_whole(_topo(2), [0, 1]) is None


def test_residue_com_is_exact_for_a_whole_molecule_spanning_the_box():
    """Atoms spread over 8 Å of a 10 Å box have a well-defined centre only when
    the molecule is reconstructed; the naive first-atom MIC centre gives 0.67
    instead of 4.0."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[0.0, 0, 0], [4.0, 0, 0], [8.0, 0, 0],
                    [1.0, 5, 0], [2.0, 5, 0], [3.0, 5, 0]])
    topo = _topo(6, mols=[1, 1, 1, 2, 2, 2], bonds=[[0, 1], [1, 2], [3, 4], [4, 5]])
    com, warn = residue_com(pos, box, topo, [0, 1, 2])
    assert warn is None
    assert com[0] == pytest.approx(4.0)


def test_residue_com_warns_when_the_selection_is_not_whole_molecules():
    """The same geometry split across three partially selected molecules has no
    well-defined centre: WARN with the size of the ambiguity."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[0.0, 0, 0], [4.0, 0, 0], [8.0, 0, 0],
                    [1.0, 5, 0], [2.0, 5, 0], [3.0, 5, 0]])
    topo = _topo(6, mols=[1, 1, 2, 2, 3, 3], bonds=[[0, 3], [1, 4], [2, 5]])
    _, warn = residue_com(pos, box, topo, [0, 1, 2], what="ligand")
    assert warn is not None and "whole molecules" in warn and "ligand" in warn
    # single atoms are always exact, even inside a partially selected molecule
    _, warn1 = residue_com(pos, box, topo, [0])
    assert warn1 is None


# ---------------------------------------------------------------------------
# pair level / periodicity of the surface layer
# ---------------------------------------------------------------------------
def test_single_atom_sasa_is_analytic():
    box = Box.orthorhombic([20.0, 20.0, 20.0])
    pos = np.array([[10.0, 10.0, 10.0]])
    total, per = sasa_shrake_rupley(pos, np.array([1.7]), 1.4, box, n_points=4096)
    want = 4 * np.pi * (1.7 + 1.4) ** 2
    assert abs(total - want) / want < 1e-3
    assert per.shape == (1,)


def test_sasa_is_invariant_under_placement_across_the_boundary():
    """The same dimer must give the same SASA wherever it sits in the cell."""
    box = Box.orthorhombic([14.0, 14.0, 14.0])
    radii = np.array([1.7, 1.7])
    mid = np.array([[7.0, 7.0, 7.0], [8.4, 7.0, 7.0]])
    across = np.array([[13.7, 7.0, 7.0], [1.1, 7.0, 7.0]])
    a_mid, _ = sasa_shrake_rupley(mid, radii, 1.4, box, n_points=1024)
    a_acr, _ = sasa_shrake_rupley(across, radii, 1.4, box, n_points=1024)
    assert abs(a_mid - a_acr) / a_mid < 1e-6


def test_sasa_shielding_across_the_face():
    """A neighbour placed just across the face must shield, not be ignored."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    radii = np.array([1.7, 1.7])
    near = np.array([[5.0, 5.0, 5.0], [5.0 + 3.0, 5.0, 5.0]])
    far = np.array([[5.0, 5.0, 5.0], [5.0 - 3.0, 5.0, 5.0]])   # wraps to the same gap
    a_near, _ = sasa_shrake_rupley(near, radii, 1.4, box, n_points=2048)
    a_far, _ = sasa_shrake_rupley(far, radii, 1.4, box, n_points=2048)
    assert abs(a_near - a_far) / a_near < 1e-6
    lone, _ = sasa_shrake_rupley(near[:1], radii[:1], 1.4, box, n_points=2048)
    assert a_near < 2 * lone                                 # mutual shielding


def test_periodic_neighbors_ball_and_query():
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[0.2, 5.0, 5.0], [9.8, 5.0, 5.0], [5.0, 5.0, 5.0]])
    tree = periodic_neighbors(pos, box, 1.0)
    hits = tree.ball([0.2, 5.0, 5.0], 1.0)
    assert set(hits.tolist()) == {0, 1}                      # 0.4 A apart via the face
    d, i = tree.query([0.2, 5.0, 5.0], k=1)
    assert i[0] == 0 and d[0] == pytest.approx(0.0)


def test_periodic_neighbors_rejects_radius_above_cutoff():
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    tree = periodic_neighbors(np.zeros((1, 3)), box, 1.0)
    with pytest.raises(LatError, match="exceeds the cutoff"):
        tree.ball([0.0, 0.0, 0.0], 2.0)


def test_periodic_nearest_uses_minimum_image():
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    centers = np.array([[9.5, 5, 5], [5.0, 5, 5]])
    idx, d = periodic_nearest(np.array([[0.2, 5, 5]]), centers, box)
    assert idx[0] == 0 and d[0] == pytest.approx(0.7)


def test_periodic_nearest_rejects_triclinic():
    box = Box.from_lammps_data([0, 0, 0], [10, 10, 10], [2.0, 0, 0])
    with pytest.raises(LatError, match="triclinic"):
        periodic_nearest(np.zeros((1, 3)), np.zeros((1, 3)), box)


def test_bond_graph_consistency_detects_a_cell_encircling_ring():
    """A ring that goes around the cell cannot be made whole: the closure bond
    needs a lattice vector, which is exactly what gmx rejects with
    'inconsistent shifts over periodic boundaries'."""
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    # chain: 1 -> 5 -> 9 (consistent); ring closure 9 -> 1 needs the wrap
    pos = np.array([[1.0, 1.0, 1.0], [5.0, 1.0, 1.0], [9.0, 1.0, 1.0],
                    [1.0, 5.0, 1.0]])
    topo = _topo(4, mols=[1, 1, 1, 1], bonds=[[0, 1], [1, 2], [2, 0], [0, 3]])
    assert bond_graph_consistent(pos, box, topo, [0, 1, 2, 3]) is False
    # the same chain without the closure is fine
    topo2 = _topo(4, mols=[1, 1, 1, 1], bonds=[[0, 1], [1, 2], [0, 3]])
    assert bond_graph_consistent(pos, box, topo2, [0, 1, 2, 3]) is True
    # no bond information at all -> unknown
    assert bond_graph_consistent(pos, box, _topo(4), None) is None


def test_residue_com_warns_and_falls_back_for_a_wrapping_network():
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[1.0, 1.0, 1.0], [5.0, 1.0, 1.0], [9.0, 1.0, 1.0],
                    [1.0, 5.0, 1.0]])
    topo = _topo(4, mols=[1, 1, 1, 1], bonds=[[0, 1], [1, 2], [2, 0], [0, 3]])
    com, warn = residue_com(pos, box, topo, [0, 1, 2, 3])
    assert warn is not None and "wraps the periodic cell" in warn
    assert np.allclose(com, box.wrap(pos.mean(axis=0)))    # wrapped mean

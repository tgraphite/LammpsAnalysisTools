"""Single entry point for periodic (PBC) geometry.

Policy
------
Modules **must not** build their own minimum-image arithmetic, KD-tree,
distance transform or connected-component labelling. Every periodic decision is
made here, on top of three thin primitives that this module wraps:

* :mod:`lat.core.pbc`      -- the cell (orthorhombic + triclinic), ``Box``
* :mod:`lat.core.neighbor` -- periodic atom-atom pair search
* :mod:`lat.core.unwrap`   -- bond-graph unwrapping (``gmx rmpbc`` analogue)

``tests/unit/test_no_ad_hoc_pbc.py`` enforces the policy.

Four layers
-----------
1. molecule level  : ``whole_positions``, ``residue_com``, ``residue_coms``
2. pair level      : ``periodic_neighbors``
3. grid level      : ``grid_axes``, ``voxel_occupancy``,
                     ``periodic_distance_transform``, ``label_periodic``
4. surface level   : ``sasa_shrake_rupley``

Why this exists
---------------
Before 0.2.0 there were four parallel "periodic geometry" implementations in
the tree, only two of them principled, and three modules silently assumed
wrapped coordinates were physical. See ``docs/07_decision_log.md`` §6.2 and
``docs/MANUAL_ALGORITHMS.md`` §0.3.

Grid layer convention
---------------------
A Cartesian voxel grid only tiles a cell periodically when the cell is
orthorhombic. The grid functions therefore raise :class:`LatError` on a
triclinic cell instead of returning subtly wrong numbers. The molecule, pair
and surface layers work for both cell types.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from .errors import LatError
from .pbc import Box
from .unwrap import build_adjacency, make_whole, molecule_groups

__all__ = [
    "topology_cache",
    "whole_positions",
    "split_molecules",
    "selection_is_whole",
    "bond_graph_consistent",
    "com_reference_spread",
    "residue_com",
    "residue_coms",
    "PeriodicNeighbors",
    "periodic_neighbors",
    "periodic_nearest",
    "periodic_distance",
    "periodic_vector",
    "grid_axes",
    "voxel_occupancy",
    "periodic_distance_transform",
    "label_periodic",
    "component_percolating",
    "percolating_components",
    "sasa_shrake_rupley",
]

# COM ambiguity below this (Å) is ignored when warning about a selection that
# is not a union of whole molecules.
COM_SPREAD_TOL = 0.1


# ---------------------------------------------------------------------------
# per-topology cache (bond adjacency + molecule grouping)
# ---------------------------------------------------------------------------
class _TopologyCache:
    __slots__ = ("topo", "adj", "adj_done", "mol_ids", "mol_groups")

    def __init__(self, topo):
        self.topo = topo
        self.adj = None
        self.adj_done = False
        self.mol_ids = None
        self.mol_groups = None


_CACHE: dict[int, _TopologyCache] = {}
_CACHE_MAX = 16


def topology_cache(topo) -> _TopologyCache:
    """Cached bond adjacency / molecule grouping for one topology."""
    key = id(topo)
    entry = _CACHE.get(key)
    if entry is None or entry.topo is not topo:
        if len(_CACHE) >= _CACHE_MAX:
            _CACHE.clear()
        entry = _TopologyCache(topo)
        _CACHE[key] = entry
    return entry


def _adjacency(topo):
    entry = topology_cache(topo)
    if not entry.adj_done:
        entry.adj = build_adjacency(topo.n_atoms, topo.bonds)
        entry.adj_done = True
    return entry.adj


def _molecules(topo):
    """(mol ids per atom, list of index arrays) or (None, None)."""
    entry = topology_cache(topo)
    if entry.mol_groups is None:
        if topo.mols is None:
            return None, None
        ids = np.asarray(topo.mols)
        entry.mol_ids = ids
        entry.mol_groups = molecule_groups(ids)
    return entry.mol_ids, entry.mol_groups


# ---------------------------------------------------------------------------
# 1. molecule level
# ---------------------------------------------------------------------------
def whole_positions(pos, box: Box, topo, indices=None) -> np.ndarray:
    """Coordinates with every *selected* molecule made whole across the cell.

    Uses the bond graph when available (correct for chains longer than half a
    box) and molecule IDs as a fallback. With ``indices`` given, only the
    connected components containing a selected atom are traversed, and only
    those rows are returned.
    """
    pos = np.asarray(pos, float)
    idx = None if indices is None else np.asarray(indices, dtype=int)
    adj = _adjacency(topo) if topo is not None else None
    if adj is not None:
        out = make_whole(pos, box, adj=adj, seeds=idx)
        return out if idx is None else out[idx]
    if topo is not None:
        mol_ids, _ = _molecules(topo)
        if mol_ids is not None:
            sub = pos if idx is None else pos[idx]
            groups = molecule_groups(mol_ids if idx is None else mol_ids[idx])
            return make_whole(sub, box, groups=groups)
    out = pos.copy()
    return out if idx is None else out[idx]


def split_molecules(topo, indices):
    """Split `indices` into whole molecules.

    Returns a list of index arrays when molecule identity is known, else None.
    Molecules only partially covered by `indices` are still returned as the
    partial group (callers detect that with `selection_is_whole`).
    """
    mol_ids, _ = _molecules(topo) if topo is not None else (None, None)
    if mol_ids is None:
        return None
    idx = np.asarray(indices, dtype=int)
    if len(idx) == 0:
        return []
    order = np.argsort(mol_ids[idx], kind="stable")
    sorted_idx = idx[order]
    sorted_mol = mol_ids[order]
    uniq, starts = np.unique(sorted_mol, return_index=True)
    ends = list(starts[1:]) + [len(sorted_mol)]
    return [sorted_idx[s:e] for s, e in zip(starts, ends)]


def selection_is_whole(topo, indices) -> bool | None:
    """True/False when molecule identity is known, None when it is not.

    A selection is "whole" when it contains either all atoms or no atoms of
    every molecule it touches, so a single reference atom is enough to
    reconstruct each molecule's physical centre exactly.
    """
    mol_ids, _ = _molecules(topo) if topo is not None else (None, None)
    if mol_ids is None:
        return None
    idx = np.asarray(indices, dtype=int)
    touched, counts = np.unique(mol_ids[idx], return_counts=True)
    _, total = np.unique(mol_ids, return_counts=True)
    # total is aligned to np.unique(mol_ids) which is the same sorted order
    return bool(np.array_equal(counts, total[np.searchsorted(np.unique(mol_ids), touched)]))


def bond_graph_consistent(pos, box: Box, topo, indices=None, tol: float = 1e-3):
    """Can the bond graph covering `indices` be unwrapped consistently?

    A spanning-tree reconstruction (like ``gmx rmpbc``) is only valid when no
    ring of the graph encircles the cell. In a 2D-periodic covalent network --
    an electrode slab, a gel, a cross-linked film -- walking around one ring
    returns to the starting atom displaced by a lattice vector, so *no* choice
    of images is right everywhere. gmx rejects this case outright
    ("inconsistent shifts over periodic boundaries"); we detect it the same way,
    by comparing every reconstructed bond length with its minimum-image length.

    Returns True/False, or None when there is no bond graph to check.
    """
    if topo is None:
        return None
    bonds = topo.bonds
    adj = _adjacency(topo)
    if adj is None or bonds is None or len(bonds) == 0:
        return None
    bonds = np.asarray(bonds, dtype=int)
    mol_ids, _ = _molecules(topo)
    if mol_ids is not None and indices is not None:
        touched = np.isin(mol_ids, np.unique(mol_ids[np.asarray(indices, dtype=int)]))
    else:
        touched = np.ones(topo.n_atoms, dtype=bool)
    seeds = None if indices is None else np.asarray(indices, dtype=int)
    keep = touched[bonds[:, 0]] & touched[bonds[:, 1]]
    b = bonds[keep]
    if len(b) == 0:
        return True
    p = np.asarray(pos, float)
    out = make_whole(p, box, adj=adj, seeds=seeds)
    d_rec = np.linalg.norm(out[b[:, 0]] - out[b[:, 1]], axis=1)
    d_mic = np.linalg.norm(box.minimum_image(p[b[:, 0]] - p[b[:, 1]]), axis=1)
    return bool(np.all(np.abs(d_rec - d_mic) <= tol))


def _com_of(p, box, masses=None):
    """Weighted centre of already-selected coordinates, wrapped into the cell."""
    p = np.asarray(p, float)
    if masses is None:
        c = p.mean(axis=0)
    else:
        w = np.asarray(masses, float)
        c = (p * w[:, None]).sum(axis=0) / w.sum()
    return box.wrap(c)


def _com_unwrapped(pos, box, atoms, masses=None):
    """Centre of the given atoms, unwrapped relative to their first atom.

    Only exact when the atoms are connected and span less than half a box; the
    caller checks that with `com_reference_spread`.
    """
    p = np.asarray(pos, float)[atoms]
    ref = p[0]
    rel = ref + box.minimum_image(p - ref)
    if masses is None:
        return rel.mean(axis=0)
    w = np.asarray(masses, float)[atoms]
    return (rel * w[:, None]).sum(axis=0) / w.sum()


def com_reference_spread(pos, box, indices, n_refs: int = 3) -> float:
    """How much the centre depends on the reference atom (Å).

    Non-zero only when the selection is *not* a union of whole molecules; the
    value bounds the error of the naive "first atom + minimum image" centre.
    """
    idx = np.asarray(indices, dtype=int)
    if len(idx) < 2:
        return 0.0
    picks = np.unique(np.linspace(0, len(idx) - 1, min(n_refs, len(idx))).astype(int))
    # rolling the atom list changes which atom _com_unwrapped uses as reference
    coms = [_com_unwrapped(pos, box, np.roll(idx, -int(k))) for k in picks]
    if len(coms) < 2:
        return 0.0
    base = coms[0]
    return max(float(np.linalg.norm(box.minimum_image(c - base))) for c in coms[1:])


def residue_com(pos, box, topo, indices, masses=None, what="selection"):
    """Centre of a residue / selection.

    Returns ``(com, warning)``. The centre is exact when the selection is a
    union of whole molecules whose bond graph can be unwrapped. Otherwise the
    naive wrapped / first-atom centre is used and a WARN explains why (gmx
    refuses both cases -- see docs/11_examples_recheck.md).
    """
    idx = np.asarray(indices, dtype=int)
    if len(idx) == 0:
        raise LatError("residue_com: empty selection")
    mol_ids, _ = _molecules(topo) if topo is not None else (None, None)
    if mol_ids is None:
        return _com_unwrapped(pos, box, idx, masses), None
    if selection_is_whole(topo, idx):
        if bond_graph_consistent(pos, box, topo, idx) is False:
            # A 2D-periodic covalent network (electrode slab, gel) cannot be
            # made whole: the spanning tree is internally inconsistent. The
            # wrapped mean (what 0.1.0 always used) is unbiased for a selection
            # whose atoms tile the cell, so it is the best available answer --
            # but say so loudly instead of reporting a displaced centre.
            return _com_of(np.asarray(pos, float)[idx], box,
                           masses=None if masses is None else np.asarray(masses, float)[idx]), (
                f"{what} contains a bond graph that wraps the periodic cell (a ring "
                f"encircles the box, e.g. a 2D network or cross-linked film), so its "
                f"molecules cannot be made whole; using the wrapped centre, whose "
                f"position along the wrapped axes is meaningless (gmx refuses this "
                f"case with 'inconsistent shifts over periodic boundaries')"
            )
        # exact: reconstruct each molecule along its bond graph, then average
        return _com_of(whole_positions(pos, box, topo, indices=idx), box,
                       masses=None if masses is None else np.asarray(masses, float)[idx]), None
    spread = com_reference_spread(pos, box, idx)
    warning = None
    if spread > COM_SPREAD_TOL:
        warning = (
            f"{what} is not a union of whole molecules and its centre depends on "
            f"the periodic-image convention by up to {spread:.3g} A; results are "
            f"approximate -- select whole molecules to make it exact"
        )
    return _com_unwrapped(pos, box, idx, masses), warning


def residue_coms(pos, box, topo, groups, masses=None, what="selection"):
    """Centres of several residue groups. Returns ``(coms, warnings)``."""
    coms = np.empty((len(groups), 3))
    warnings: list[str] = []
    seen: set[str] = set()
    for k, g in enumerate(groups):
        c, w = residue_com(pos, box, topo, g, masses=masses, what=what)
        coms[k] = c
        if w and w not in seen:
            seen.add(w)
            warnings.append(w)
    return coms, warnings


# ---------------------------------------------------------------------------
# 2. pair level
# ---------------------------------------------------------------------------
class PeriodicNeighbors:
    """Periodic neighbour queries against one fixed coordinate set.

    Orthorhombic cells use ``cKDTree(boxsize=...)``; triclinic cells use an
    image tiling sized by the largest query radius. All returned distances are
    minimum-image distances.
    """

    def __init__(self, pos, box: Box, cutoff: float):
        pos = box.wrap(np.asarray(pos, float))
        self.box = box
        self.pos = pos
        self.n = len(pos)
        self.cutoff = float(cutoff)
        if box.is_orthorhombic:
            self._tree = cKDTree(pos - box.origin, boxsize=box.lengths)
            self._map = None
        else:
            shifts = box.image_shifts(self.cutoff)
            ghost = (pos[None, :, :] + (shifts[:, None, :] @ box.H.T)).reshape(-1, 3)
            self._tree = cKDTree(ghost)
            self._map = np.tile(np.arange(self.n), len(shifts))

    def _query_point(self, x):
        x = np.asarray(x, float).reshape(3)
        if self._map is None:
            return self.box.wrap(x) - self.box.origin
        return x

    def ball(self, x, radius: float) -> np.ndarray:
        """Indices within `radius` of `x` (minimum image), self included."""
        if radius > self.cutoff + 1e-12:
            raise LatError(
                f"periodic_neighbors: query radius {radius:g} A exceeds the "
                f"cutoff {self.cutoff:g} A the tree was built with"
            )
        hits = self._tree.query_ball_point(self._query_point(x), float(radius))
        if self._map is None:
            return np.asarray(hits, dtype=int)
        return np.unique(self._map[np.asarray(hits, dtype=int)])

    def query(self, x, k: int = 1):
        """Distances and indices of the k nearest periodic images of atoms."""
        if self._map is None:
            d, i = self._tree.query(self._query_point(x), k=int(k))
            return np.atleast_1d(d), np.atleast_1d(i).astype(int)
        d, g = self._tree.query(np.asarray(x, float).reshape(3), k=int(k))
        d = np.atleast_1d(d)
        g = np.atleast_1d(g).astype(int)
        idx = self._map[g]
        keep = np.ones(len(idx), dtype=bool)
        for m in range(1, len(idx)):
            if np.any(idx[:m] == idx[m]):
                keep[m] = False
        return d[keep], idx[keep]


def periodic_neighbors(pos, box: Box, cutoff: float) -> PeriodicNeighbors:
    return PeriodicNeighbors(pos, box, cutoff)


def periodic_nearest(points, centers, box: Box):
    """Nearest `centers` entry (minimum image) for each row of `points`.

    Returns ``(index, distance)``. Orthorhombic only: a Cartesian point cloud
    has no triclinic periodic meaning.
    """
    if not box.is_orthorhombic:
        raise LatError("periodic_nearest: a triclinic cell has no periodic "
                       "Cartesian point cloud (grid analyses need an orthorhombic box)")
    c = box.wrap(np.asarray(centers, float)) - box.origin
    p = box.wrap(np.asarray(points, float)) - box.origin
    tree = cKDTree(c, boxsize=box.lengths)
    d, i = tree.query(p)
    return np.atleast_1d(i).astype(int), np.atleast_1d(d)


def periodic_vector(a, b, box: Box) -> np.ndarray:
    """Minimum-image displacement ``a - b`` (row vectors or single points)."""
    return box.minimum_image(np.asarray(a, float) - np.asarray(b, float))


def periodic_distance(a, b, box: Box) -> float:
    """Minimum-image distance between one point and one point (Å)."""
    return float(np.linalg.norm(periodic_vector(a, b, box)))


# ---------------------------------------------------------------------------
# 3. grid level (orthorhombic cells only)
# ---------------------------------------------------------------------------
def _require_ortho(box: Box, who: str) -> None:
    if not box.is_orthorhombic:
        raise LatError(
            f"{who}: a Cartesian voxel grid only tiles an orthorhombic cell "
            f"periodically; this box is triclinic (tilt {box.tilt}). Use a "
            f"molecule/pair-based analysis instead, or make the cell orthorhombic"
        )


def grid_axes(box: Box, dr: float, max_points: float = 5e7):
    """Regular grid cell centres covering the box.

    Returns ``(axes, spacing, coarsened)``. `dr` is increased just enough to
    keep the voxel count under `max_points`.
    """
    _require_ortho(box, "grid_axes")
    dr = float(dr)
    if dr <= 0:
        raise LatError(f"grid spacing must be positive, got {dr}")
    L = box.lengths
    n = np.maximum(np.round(L / dr), 1)
    if float(np.prod(n)) > max_points:
        n = np.maximum(np.floor(n * (max_points / float(np.prod(n))) ** (1.0 / 3.0)), 1)
        coarsened = True
    else:
        coarsened = False
    spacing = L / n
    axes = [box.lo[k] + (np.arange(int(n[k])) + 0.5) * spacing[k] for k in range(3)]
    return axes, spacing, coarsened


def voxel_occupancy(pos, box: Box, axes, radii) -> np.ndarray:
    """Mark voxels whose centre is within `radii` of an atom (periodic).

    Vectorised per atom over its local voxel block.
    """
    _require_ortho(box, "voxel_occupancy")
    dims = np.array([len(a) for a in axes], dtype=int)
    spacing = np.array([a[1] - a[0] for a in axes]) if all(len(a) > 1 for a in axes) \
        else box.lengths / dims
    occ = np.zeros(tuple(int(d) for d in dims), dtype=bool)
    L = box.lengths
    pos = box.wrap(np.asarray(pos, float))
    los = np.array([a[0] for a in axes])
    radii = np.asarray(radii, float)

    def _axis(k, x, r):
        i0 = int(np.floor((x - r - los[k]) / spacing[k]))
        i1 = int(np.ceil((x + r - los[k]) / spacing[k]))
        ii = np.arange(i0, i1 + 1)
        iw = np.mod(ii, dims[k])
        ax = axes[k][iw]
        d = x - ax
        d -= L[k] * np.round(d / L[k])          # minimum image
        d = np.abs(d)
        if len(np.unique(iw)) != len(iw):
            uw, inv = np.unique(iw, return_inverse=True)
            dm = np.full(len(uw), np.inf)
            np.minimum.at(dm, inv.ravel(), d)
            return uw, dm
        return iw, d

    for atom in range(len(pos)):
        r = radii[atom]
        if not np.isfinite(r) or r <= 0:
            continue
        ux, dx = _axis(0, pos[atom, 0], r)
        uy, dy = _axis(1, pos[atom, 1], r)
        uz, dz = _axis(2, pos[atom, 2], r)
        mask = (dx[:, None, None] ** 2 + dy[None, :, None] ** 2
                + dz[None, None, :] ** 2) <= r * r
        occ[np.ix_(ux, uy, uz)] |= mask
    return occ


_SQ_CACHE: dict[tuple[int, float], np.ndarray] = {}


def _circ_sqdist(n: int, h: float) -> np.ndarray:
    """(n, n) squared *circular* distances between grid lines, step `h`."""
    key = (int(n), round(float(h), 9))
    D = _SQ_CACHE.get(key)
    if D is None:
        i = np.arange(int(n))
        d = np.abs(i[:, None] - i[None, :])
        d = np.minimum(d, int(n) - d)
        D = (d * float(h)) ** 2
        if len(_SQ_CACHE) > 32:
            _SQ_CACHE.clear()
        _SQ_CACHE[key] = D
    return D


def _minplus_axis(f, axis: int, D: np.ndarray, budget: int = 64 << 20):
    f = np.moveaxis(f, axis, -1)
    shape = f.shape
    lines = np.ascontiguousarray(f).reshape(-1, shape[-1])
    n = shape[-1]
    out = np.empty_like(lines)
    chunk = max(1, int(budget // max(n * n * 8, 1)))
    for s in range(0, len(lines), chunk):
        blk = lines[s:s + chunk]
        out[s:s + chunk] = (blk[:, None, :] + D[None, :, :]).min(axis=-1)
    return np.moveaxis(out.reshape(shape), -1, axis)


def periodic_distance_transform(mask, spacing) -> np.ndarray:
    """Exact toroidal Euclidean distance (Å) from every voxel to the nearest
    True voxel of `mask`.

    `mask` tiles the cell and `spacing` is the per-axis voxel pitch. The
    minimum-image EDT is separable: each axis gets a circular 1-D
    (min,+) squared-distance transform.
    """
    mask = np.ascontiguousarray(mask, dtype=bool)
    if mask.ndim != 3:
        raise LatError(f"periodic_distance_transform: expected a 3-D mask, got {mask.ndim}-D")
    if not mask.any():
        raise LatError("periodic_distance_transform: the mask is empty")
    spacing = np.asarray(spacing, float).ravel()
    if spacing.size != 3 or np.any(spacing <= 0):
        raise LatError(f"periodic_distance_transform: bad spacing {spacing!r}")
    f = np.where(mask, 0.0, np.inf)
    for axis in range(3):
        f = _minplus_axis(f, axis, _circ_sqdist(f.shape[axis], spacing[axis]))
    return np.sqrt(f)


def label_periodic(mask, connectivity: int = 6):
    """Connected components of `mask` on the torus.

    Returns ``(labels, n_components)``; ``labels`` has the mask shape with 0
    for background and 1..n for components. Components that leave through one
    face and re-enter through the opposite one are a single component.

    Connectivity 6 = faces, 26 = faces+edges+corners.

    Implementation: label a wrap-padded copy of the mask, then merge label
    pairs that denote the *same physical voxel* through the pad (a padded slab
    and the face it copies can end up in different components, because the pad
    only makes one of the two wrap adjacencies traversable).
    """
    mask = np.ascontiguousarray(mask, dtype=bool)
    if mask.ndim != 3:
        raise LatError(f"label_periodic: expected a 3-D mask, got {mask.ndim}-D")
    if int(connectivity) not in (6, 26):
        raise LatError(f"label_periodic: connectivity must be 6 or 26, got {connectivity}")
    if not mask.any():
        return np.zeros(mask.shape, dtype=np.int32), 0
    structure = ndimage.generate_binary_structure(3, 1 if int(connectivity) == 6 else 3)
    padded = np.pad(mask, 1, mode="wrap")
    lab, n_lab = ndimage.label(padded, structure=structure)

    parent = list(range(n_lab + 1))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(int(a)), find(int(b))
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    for axis in range(3):
        n = mask.shape[axis]
        # padded index 0 copies the last plane, index n+1 copies the first
        for i, j in ((0, n), (n + 1, 1)):
            lo = np.take(lab, i, axis=axis).ravel()
            hi = np.take(lab, j, axis=axis).ravel()
            pairs = np.unique(np.stack([lo, hi], axis=1), axis=0)
            pairs = pairs[(pairs[:, 0] > 0) | (pairs[:, 1] > 0)]
            for a, b in pairs:
                union(a, b)

    core = lab[1:-1, 1:-1, 1:-1]
    roots = np.array([find(r) for r in range(n_lab + 1)])
    uniq, inv = np.unique(roots[core], return_inverse=True)
    inv = inv.reshape(core.shape)
    if uniq.size and uniq[0] == 0:
        return inv.astype(np.int32), int(uniq.size - 1)
    return (inv + 1).astype(np.int32), int(uniq.size)


def percolating_components(labels) -> set[int]:
    """Component ids that span the cell along at least one axis.

    Percolation is a property of the *cell*: the component must reach both
    opposite faces, not merely be connected through the periodic image.
    """
    out: set[int] = set()
    for lo, hi in ((labels[0], labels[-1]),
                   (labels[:, 0], labels[:, -1]),
                   (labels[:, :, 0], labels[:, :, -1])):
        a = set(int(v) for v in np.unique(lo))
        b = set(int(v) for v in np.unique(hi))
        out |= (a & b) - {0}
    return out


def component_percolating(labels, cid: int) -> bool:
    m = labels == int(cid)
    return bool(
        (m[0].any() and m[-1].any())
        or (m[:, 0].any() and m[:, -1].any())
        or (m[:, :, 0].any() and m[:, :, -1].any())
    )


# ---------------------------------------------------------------------------
# 4. surface level
# ---------------------------------------------------------------------------
_DIR_CACHE: dict[int, np.ndarray] = {}


def _sphere_dirs(n_points: int) -> np.ndarray:
    """Fibonacci-sphere unit directions (deterministic, no RNG)."""
    n = int(n_points)
    if n < 1:
        raise LatError(f"sasa: n_points must be >= 1, got {n}")
    D = _DIR_CACHE.get(n)
    if D is None:
        k = np.arange(n) + 0.5
        phi = np.arccos(1 - 2 * k / n)
        theta = np.pi * (1 + 5 ** 0.5) * k
        D = np.column_stack(
            [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
        )
        if len(_DIR_CACHE) > 16:
            _DIR_CACHE.clear()
        _DIR_CACHE[n] = D
    return D


def sasa_shrake_rupley(pos, radii, probe_radius, box: Box, n_points: int = 92):
    """Periodic Shrake-Rupley solvent accessible surface area.

    `radii` are the vdW radii of exactly the atoms in `pos`, and the returned
    total/per-atom areas are for those atoms only. Neighbour search and every
    sphere-point test are minimum-image, so atoms near a face are treated the
    same as atoms in the middle of the cell.
    """
    pos = np.asarray(pos, float)
    radii = np.asarray(radii, float)
    if len(pos) != len(radii):
        raise LatError("sasa_shrake_rupley: pos and radii must have the same length")
    if probe_radius < 0:
        raise LatError(f"sasa_shrake_rupley: probe radius must be >= 0, got {probe_radius}")
    if len(pos) == 0:
        return 0.0, np.zeros(0)
    dirs = _sphere_dirs(n_points)
    rr = radii + float(probe_radius)
    rmax = float(rr.max())
    tree = periodic_neighbors(pos, box, 2.0 * rmax)
    pos_w = box.wrap(pos)
    total = 0.0
    per_atom = np.zeros(len(pos))
    for i in range(len(pos)):
        Ri = rr[i]
        neigh = tree.ball(pos_w[i], Ri + rmax)
        pts = pos_w[i] + dirs * Ri
        exposed = np.ones(len(dirs), dtype=bool)
        for j in neigh:
            if int(j) == i:
                continue
            d2 = np.sum(box.minimum_image(pts - pos_w[int(j)]) ** 2, axis=1)
            exposed &= d2 > rr[int(j)] ** 2
        area = 4 * np.pi * Ri ** 2 * float(exposed.mean())
        per_atom[i] = area
        total += area
    return total, per_atom

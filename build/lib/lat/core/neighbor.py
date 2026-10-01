"""Neighbor search: periodic pairs within a cutoff (orthorhombic / triclinic)."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .pbc import Box


def _empty():
    return (
        np.empty(0, dtype=int),
        np.empty(0, dtype=int),
        np.empty(0, dtype=float),
    )


def pairs_within(pos_a, pos_b, box: Box, cutoff: float, same: bool = False):
    """Return (indices_a, indices_b, distances) for pairs within `cutoff` (Å).

    `same=True` treats A and B as the same set: self-pairs are excluded and each
    unordered pair is returned once.
    """
    # cKDTree(boxsize=...) requires coordinates in [0, L); wrap() returns them
    # in [origin, origin+L). Shift by origin for the periodic tree (differences
    # are unaffected, so MIC distances stay correct).
    pos_a = box.wrap(np.asarray(pos_a, float)) - box.origin
    pos_b = pos_a if same else box.wrap(np.asarray(pos_b, float)) - box.origin

    if box.is_orthorhombic:
        L = box.lengths
        ta = cKDTree(pos_a, boxsize=L)
        if same:
            pr = ta.query_pairs(float(cutoff), output_type="ndarray")
            if pr.size == 0:
                return _empty()
            i, j = pr[:, 0], pr[:, 1]
            d = np.linalg.norm(box.minimum_image(pos_a[i] - pos_a[j]), axis=1)
            return i, j, d
        lists = ta.query_ball_point(pos_b, float(cutoff))
        i, j = [], []
        for jb, lst in enumerate(lists):
            for ia in lst:
                i.append(ia)
                j.append(jb)
        if not i:
            return _empty()
        i = np.asarray(i, int)
        j = np.asarray(j, int)
        d = np.linalg.norm(box.minimum_image(pos_a[i] - pos_b[j]), axis=1)
        return i, j, d

    # triclinic: replicate B over neighboring images
    shifts = box.image_shifts(float(cutoff))
    ghost = (pos_b[None, :, :] + (shifts[:, None, :] @ box.H.T)).reshape(-1, 3)
    map_j = np.tile(np.arange(len(pos_b)), len(shifts))
    tree = cKDTree(ghost)
    lists = tree.query_ball_point(pos_a, float(cutoff))

    if same:
        best: dict[tuple[int, int], float] = {}
        for ia, lst in enumerate(lists):
            for g in lst:
                j = int(map_j[g])
                if j == ia:
                    continue
                d = float(np.linalg.norm(pos_a[ia] - ghost[g]))
                key = (ia, j) if ia < j else (j, ia)
                if key not in best or d < best[key]:
                    best[key] = d
        if not best:
            return _empty()
        i = np.asarray([k[0] for k in best], int)
        j = np.asarray([k[1] for k in best], int)
        d = np.asarray(list(best.values()), float)
        return i, j, d

    i, j, d = [], [], []
    for ia, lst in enumerate(lists):
        for g in lst:
            i.append(ia)
            j.append(int(map_j[g]))
            d.append(float(np.linalg.norm(pos_a[ia] - ghost[g])))
    if not i:
        return _empty()
    return np.asarray(i, int), np.asarray(j, int), np.asarray(d, float)


def neighbor_counts(pos_a, pos_b, box: Box, cutoff: float) -> np.ndarray:
    """Per-A count of B atoms within cutoff (self excluded when same set)."""
    i, j, _ = pairs_within(pos_a, pos_b, box, cutoff, same=False)
    counts = np.zeros(len(pos_a), dtype=int)
    if len(i):
        np.add.at(counts, i, 1)
    return counts


def distance_matrix(pos_a, pos_b, box: Box) -> np.ndarray:
    """Full (Na, Nb) MIC distance matrix (small systems only)."""
    a = np.asarray(pos_a, float)
    b = np.asarray(pos_b, float)
    d = a[:, None, :] - b[None, :, :]
    return np.linalg.norm(box.minimum_image(d), axis=2)

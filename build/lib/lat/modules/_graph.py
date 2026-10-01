"""Residue-level bipartite graph engine shared by agg / cluster.

Neighbors -> graph -> connected components -> endpoint-label classification.
One edge per residue pair (a polydentate anion counts once).
"""

from __future__ import annotations

from collections import Counter

import numpy as np

from ..core.errors import LatError
from ..core.neighbor import pairs_within


def residues(topo, indices):
    """Group atom indices into residues (molecules). Returns (mols, atoms_by_mol)."""
    if topo.mols is None:
        raise LatError("this analysis needs molecule IDs (mol or residue information)")
    mols = topo.mols[indices]
    order = np.argsort(mols, kind="stable")
    sorted_idx = indices[order]
    sorted_mols = mols[order]
    uniq, starts = np.unique(sorted_mols, return_index=True)
    atoms_by_mol = [sorted_idx[s:e] for s, e in zip(starts, list(starts[1:]) + [len(sorted_idx)])]
    return uniq, atoms_by_mol


def _pack(atoms_by_mol):
    all_atoms = np.concatenate(atoms_by_mol) if atoms_by_mol else np.empty(0, int)
    res_of_atom = np.concatenate([np.full(len(a), k) for k, a in enumerate(atoms_by_mol)]) \
        if atoms_by_mol else np.empty(0, int)
    return all_atoms, res_of_atom


def edges_min(pos, mols_a, atoms_a, mols_b, atoms_b, box, cutoff):
    a_all, a_res = _pack(atoms_a)
    b_all, b_res = _pack(atoms_b)
    if len(a_all) == 0 or len(b_all) == 0:
        return set()
    ia, ib, _ = pairs_within(pos[a_all], pos[b_all], box, cutoff, same=False)
    return set(zip(a_res[ia].tolist(), b_res[ib].tolist()))


def _coms(pos, atoms_by_mol, box):
    """Per-residue COM, unwrapped within the residue (PBC-safe)."""
    out = np.empty((len(atoms_by_mol), 3))
    for k, a in enumerate(atoms_by_mol):
        ref = pos[a[0]]
        d = box.minimum_image(pos[a] - ref)
        out[k] = ref + d.mean(axis=0)
    return out


def edges_center(pos, atoms_a, atoms_b, box, cutoff):
    if not atoms_a or not atoms_b:
        return set()
    com_a = _coms(pos, atoms_a, box)
    com_b = _coms(pos, atoms_b, box)
    ia, ib, _ = pairs_within(com_a, com_b, box, cutoff, same=False)
    return set(zip(ia.tolist(), ib.tolist()))


def edges(pos, atoms_a, atoms_b, box, cutoff, mode):
    if mode == "min":
        return edges_min(pos, None, atoms_a, None, atoms_b, box, cutoff)
    if mode == "center":
        return edges_center(pos, atoms_a, atoms_b, box, cutoff)
    raise LatError(f"-mode must be min|center, got {mode!r}")


class Components:
    def __init__(self, n_a, n_b, edge_set):
        self.n_a = n_a
        self.n_b = n_b
        self.edge_set = edge_set
        parent = list(range(n_a + n_b))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for a, b in edge_set:
            ra, rb = find(a), find(n_a + b)
            if ra != rb:
                parent[ra] = rb

        comp: dict[int, dict] = {}
        for a in range(n_a):
            comp.setdefault(find(a), {"a": 0, "b": 0})["a"] += 1
        for b in range(n_b):
            comp.setdefault(find(n_a + b), {"a": 0, "b": 0})["b"] += 1
        self.comp = comp
        self.root_of_a = [find(a) for a in range(n_a)]
        self.root_of_b = [find(n_a + b) for b in range(n_b)]
        self.deg_a = Counter(a for a, _ in edge_set)
        self.deg_b = Counter(b for _, b in edge_set)
        members: dict[int, dict] = {}
        for a in range(n_a):
            members.setdefault(self.root_of_a[a], {"a": [], "b": []})["a"].append(a)
        for b in range(n_b):
            members.setdefault(self.root_of_b[b], {"a": [], "b": []})["b"].append(b)
        self.members = members

    def type_counts(self) -> Counter:
        c = Counter()
        for info in self.comp.values():
            c[(info["a"], info["b"])] += 1
        return c

    def label(self, side: str) -> str:
        key = "a" if side == "a" else "b"
        other = "b" if side == "a" else "a"
        roots = self.root_of_a if side == "a" else self.root_of_b
        counts = Counter()
        for r in roots:
            info = self.comp[r]
            if info[key] == 1 and info[other] == 0:
                counts["SSIP"] += 1
            elif info[key] == 1 and info[other] == 1:
                counts["CIP"] += 1
            else:
                counts["AGG"] += 1
        return counts

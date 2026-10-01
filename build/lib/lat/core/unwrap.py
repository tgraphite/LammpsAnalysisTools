"""Make molecules whole across the periodic boundary (gmx rmpbc analogue).

Wrapped coordinates split any molecule that spans a box face. Analyses such as
Rg and RMSD must reconstruct whole molecules first, otherwise they are biased
low. We follow bonded connectivity (spanning forest) rather than a naive
COM-to-first-atom minimum image, which fails for chains longer than half a box.
"""

from __future__ import annotations

import numpy as np


def build_adjacency(n_atoms: int, bonds):
    """Adjacency list from a (n_bonds, 2) integer bond array (or None)."""
    if bonds is None or len(bonds) == 0:
        return None
    adj = [[] for _ in range(n_atoms)]
    for i, j in np.asarray(bonds, dtype=int):
        adj[int(i)].append(int(j))
        adj[int(j)].append(int(i))
    return adj


def molecule_groups(mols):
    order = np.argsort(mols, kind="stable")
    sm = mols[order]
    uniq, starts = np.unique(sm, return_index=True)
    return [order[s:e] for s, e in zip(starts, list(starts[1:]) + [len(sm)])]


def make_whole(pos, box, adj=None, groups=None):
    """Return a copy of `pos` with each connected molecule made whole.

    `adj` (from build_adjacency) is preferred; `groups` (from molecule_groups)
    is the fallback when bond information is unavailable.
    """
    out = np.asarray(pos, float).copy()
    if adj is not None:
        n = len(out)
        seen = np.zeros(n, dtype=bool)
        for s in range(n):
            if seen[s]:
                continue
            seen[s] = True
            stack = [s]
            while stack:
                a = stack.pop()
                for b in adj[a]:
                    if not seen[b]:
                        out[b] = out[a] + box.minimum_image(pos[b] - pos[a])
                        seen[b] = True
                        stack.append(b)
        return out
    if groups is not None:
        for g in groups:
            ref = pos[g[0]]
            out[g] = ref + box.minimum_image(pos[g] - ref)
    return out

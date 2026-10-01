"""lat rdf — radial distribution function and coordination number."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.grid import rdf_bins
from ..core.neighbor import pairs_within
from ..core.table import Field, Table
from .base import Module, Param
from . import register


def _sanitize(text: str) -> str:
    out = text.replace(" ", "_")
    for ch in "(),":
        out = out.replace(ch, "")
    return out


def _label(topo, selection) -> str:
    types = np.unique(topo.types[selection.indices])
    els = {topo.type_element.get(int(t)) for t in types}
    if len(els) == 1 and None not in els:
        return els.pop()
    return _sanitize(selection.literal)


def _first_shell_min(r, g):
    n = len(g)
    peak = None
    for i in range(1, n - 1):
        if g[i] >= g[i - 1] and g[i] >= g[i + 1] and g[i] > 1.0:
            peak = i
            break
    if peak is None:
        return None
    for j in range(peak + 1, n - 1):
        if g[j] <= g[j - 1] and g[j] <= g[j + 1] and g[j] < 0.5:
            return j
    return None


@register
class RDF(Module):
    name = "rdf"
    summary = "radial distribution function g(r) and coordination number"
    help_doc = "algorithms/rdf.md"
    inputs = ["-s", "-sel1", "-sel2"]
    outputs = ["<base>.rdf.out.csv", "<base>.rdf.cn.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), required=True, help="first selection"),
        Param("sel2", flags=("-sel2",), required=True, help="second selection"),
        Param("nbins", flags=("-nbins",), config="modules.rdf.nbins", type="int"),
        Param("cutoff", flags=("-cutoff",), config="modules.rdf.cutoff", type="float"),
        Param("frame_range", flags=("-frame_range",), help="start:stop[:stride]"),
        Param("stride", flags=("-stride",), type="int"),
        Param("intramol", flags=("-intramol",), config="modules.rdf.exclude_intramolecular",
              type="flag", help="exclude pairs in the same molecule (needs mol/Bonds)"),
    ]

    def run(self, ctx, *, sel1, sel2, nbins, cutoff, frame_range, stride, intramol):
        topo = ctx.topo
        s1 = ctx.selection(sel1)
        s2 = ctx.selection(sel2)
        if s1.n_atoms == 0:
            raise LatError(f"selection {sel1!r} matched 0 atoms")
        if s2.n_atoms == 0:
            raise LatError(f"selection {sel2!r} matched 0 atoms")

        nbins = int(nbins)
        cutoff = float(cutoff)
        edges, r = rdf_bins(nbins, cutoff)

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        nframes = len(frames)
        same = (s1.n_atoms == s2.n_atoms
                and np.array_equal(np.sort(s1.indices), np.sort(s2.indices)))
        factor = 2.0 if same else 1.0
        n_a, n_b = s1.n_atoms, s2.n_atoms

        warnings: list[str] = []
        mols = topo.mols if intramol else None
        if intramol and mols is None:
            warnings.append("intramolecular exclusion requested but no molecule IDs available; ignored")
            mols = None

        hist_sum = np.zeros(nbins)
        hist_v_sum = np.zeros(nbins)
        idx1, idx2 = s1.indices, s2.indices
        min_edge = None
        for fr in ctx.frames(frames):
            pos = fr.positions
            box = fr.box
            min_edge = float(box.lengths.min()) if min_edge is None else min(min_edge, float(box.lengths.min()))
            i, j, d = pairs_within(pos[idx1], None if same else pos[idx2], box, cutoff, same=same)
            if not same and len(i):
                # overlapping (but not identical) selections: drop self pairs
                keep = idx1[i] != idx2[j]
                i, j, d = i[keep], j[keep], d[keep]
            if mols is not None and len(i):
                mi = mols[idx1][i] if same else mols[idx1][i]
                mj = mols[idx1][j] if same else mols[idx2][j]
                keep = mi != mj
                i, j, d = i[keep], j[keep], d[keep]
            counts, _ = np.histogram(d, bins=edges)
            hist_sum += counts
            hist_v_sum += counts * box.volume

        if min_edge is not None and cutoff > 0.5 * min_edge:
            warnings.append(f"cutoff {cutoff:g} A exceeds half the minimum box edge ({0.5*min_edge:g} A)")
        if nframes < 2:
            warnings.append("fewer than 2 frames; g(r) statistics are poor")

        # exact spherical shell volume (4/3)pi(r_hi^3 - r_lo^3)
        shell = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
        with np.errstate(divide="ignore", invalid="ignore"):
            g = factor * hist_v_sum / (nframes * n_a * n_b * shell)
        g = np.nan_to_num(g)
        cn = factor * np.cumsum(hist_sum) / (nframes * n_a)

        l1, l2 = _label(topo, s1), _label(topo, s2)
        selection = {"sel1": s1.literal, "sel2": s2.literal}
        fr_str = ctx.frame_range_str(frame_range, frames)

        scalars = {}
        if len(g):
            ipk = int(np.argmax(g))
            scalars["r_first_peak/A"] = float(r[ipk])
            scalars["g_first_peak"] = float(g[ipk])
        imin = _first_shell_min(r, g)
        if imin is not None:
            scalars["r_first_min/A"] = float(r[imin])
            scalars["cn_first_shell"] = float(cn[imin])
        else:
            warnings.append(
                "no first-shell minimum found (g never returned below 0.5 after the "
                "first peak); r_first_min/cn_first_shell omitted"
            )

        g_table = Table(
            title=f"RDF {s1.literal} - {s2.literal}",
            columns=[Field("r", "Angstrom"), Field(f"g_{l1}_{l2}", "1")],
            data=np.column_stack([r, g]),
            tag="rdf",
            scalars=scalars,
            warnings=warnings,
            meta={"module": "rdf"},
            selection=selection,
            frame_range=fr_str,
        )
        cn_table = Table(
            title=f"Coordination number {s1.literal} - {s2.literal}",
            columns=[Field("r", "Angstrom"), Field(f"cn_{l1}_{l2}", "1")],
            data=np.column_stack([r, cn]),
            tag="rdf.cn",
            scalars=scalars,
            warnings=warnings,
            meta={"module": "rdf"},
            selection=selection,
            frame_range=fr_str,
        )
        return [g_table, cn_table]

"""lat cluster — coordination-number and cluster-size distributions."""

from __future__ import annotations

from collections import Counter

import numpy as np

from ..core.errors import LatError
from ..core.table import Field, Table
from ._graph import Components, edges, residues
from .base import Module, Param
from . import register


@register
class Cluster(Module):
    name = "cluster"
    summary = "coordination-number distribution (and cluster sizes)"
    help_doc = "algorithms/agg.md"
    inputs = ["-s", "-sel1", "-sel2", "-cutoff"]
    outputs = ["<base>.cluster.cn.out.csv", "<base>.cluster.size.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), required=True),
        Param("sel2", flags=("-sel2",), required=True),
        Param("cutoff", flags=("-cutoff",), config="modules.cluster.cutoff", type="float",
              required=True),
        Param("mode", flags=("-mode",), config="modules.cluster.mode", choices=("min", "center")),
        Param("focus", flags=("-focus",), config="modules.cluster.focus", type="flag",
              help="skip graph, only per-center coordination numbers"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel1, sel2, cutoff, mode, focus, frame_range, stride):
        topo = ctx.topo
        if cutoff is None:
            raise LatError("cluster: -cutoff is required (no default)")
        s1 = ctx.selection(sel1)
        s2 = ctx.selection(sel2)
        am, aa = residues(topo, s1.indices)
        bm, ba = residues(topo, s2.indices)

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        cnt_dist = Counter()
        size_dist = Counter()
        n_a_total = 0
        warnings: list[str] = []
        for fr in ctx.frames(frames):
            pos = fr.positions
            e = edges(pos, aa, ba, fr.box, float(cutoff), mode)
            comp = Components(len(am), len(bm), e)
            for a in range(len(am)):
                cnt_dist[comp.deg_a.get(a, 0)] += 1
            n_a_total += len(am)
            if not focus:
                for mem in comp.members.values():
                    size_dist[len(mem["a"]) + len(mem["b"])] += 1

        cn_rows = [[k, cnt_dist[k], cnt_dist[k] / n_a_total] for k in sorted(cnt_dist)]
        tables = [Table(
            title=f"Coordination number distribution (A={sel1}, B={sel2})",
            columns=[Field("cn", "1"), Field("count", "1"), Field("p", "1")],
            data=np.array(cn_rows, float), tag="cluster.cn", warnings=warnings,
            scalars={"mean_cn": float(sum(k * v for k, v in cnt_dist.items()) / n_a_total)},
            meta={"module": "cluster"},
            selection={"sel1": sel1, "sel2": sel2, "cutoff": f"{cutoff} A", "mode": mode},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]
        if not focus:
            total = sum(size_dist.values())
            size_rows = [[k, size_dist[k], size_dist[k] / total] for k in sorted(size_dist)]
            tables.append(Table(
                title="Cluster size distribution",
                columns=[Field("size", "1"), Field("count", "1"), Field("p", "1")],
                data=np.array(size_rows, float), tag="cluster.size", warnings=warnings,
                meta={"module": "cluster"},
                selection={"sel1": sel1, "sel2": sel2, "cutoff": f"{cutoff} A", "mode": mode},
                frame_range=ctx.frame_range_str(frame_range, frames),
            ))
        return tables

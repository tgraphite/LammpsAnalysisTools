"""lat agg — SSIP / CIP / AGG speciation via a residue-level bipartite graph."""

from __future__ import annotations

from collections import Counter

import numpy as np

from ..core.errors import LatError
from ..core.table import Field, Table
from ._graph import Components, edges, residues
from .base import Module, Param
from . import register


def _sorted_types(keys):
    return sorted(keys, key=lambda k: (k[0] + k[1], k[0], k[1]))


@register
class Agg(Module):
    name = "agg"
    summary = "ion aggregation speciation (SSIP/CIP/AGG) from a residue graph"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel1", "-sel2", "-cutoff"]
    outputs = ["<base>.agg.out.csv", "<base>.agg.labels.out.csv", "<base>.agg.speciation.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), required=True, help="center group (e.g. Li, type 2)"),
        Param("sel2", flags=("-sel2",), required=True, help="counter group (e.g. FSI)",),
        Param("cutoff", flags=("-cutoff",), config="modules.agg.cutoff", type="float",
              required=True, help="required; no default"),
        Param("mode", flags=("-mode",), config="modules.agg.mode", choices=("min", "center"),
              help="edge criterion: min atom-atom or residue COM"),
        Param("solv", flags=("-solv",), config="modules.agg.solv", help="optional solvent group"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel1, sel2, cutoff, mode, solv, frame_range, stride):
        topo = ctx.topo
        if cutoff is None:
            raise LatError("agg: -cutoff is required (no default)")
        s1 = ctx.selection(sel1)
        s2 = ctx.selection(sel2)
        if s1.n_atoms == 0 or s2.n_atoms == 0:
            raise LatError("agg: empty selection")
        am, aa = residues(topo, s1.indices)
        bm, ba = residues(topo, s2.indices)

        solv_res = None
        if solv:
            ss = ctx.selection(solv)
            solv_res = residues(topo, ss.indices)

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        if topo.charges is not None:
            qa = float(np.sum(topo.charges[s1.indices]))
            qb = float(np.sum(topo.charges[s2.indices]))
            scale = max(abs(qa), abs(qb), 1.0)
            if abs(qa + qb) > 0.1 * scale:
                warnings.append(f"charge imbalance: sel1 sum={qa:+.3f}, sel2 sum={qb:+.3f}")

        typekeys: set = set()
        speckeys: set = set()
        type_rows = []
        label_rows = []
        spec_rows = []
        for fr in ctx.frames(frames):
            pos = fr.positions
            e = edges(pos, aa, ba, fr.box, float(cutoff), mode)
            comp = Components(len(am), len(bm), e)
            tc = comp.type_counts()
            typekeys |= set(tc)
            type_rows.append((fr.index, tc))
            if len(comp.comp) == 1 and (len(am) + len(bm)) > 2:
                warnings.append(
                    "all residues are in a single connected component; cutoff may be too large"
                )
            la = comp.label("a")
            lb = comp.label("b")
            na, nb = len(am), len(bm)
            label_rows.append((fr.index,
                               la["SSIP"] / na, la["CIP"] / na, la["AGG"] / na,
                               lb["SSIP"] / nb, lb["CIP"] / nb, lb["AGG"] / nb))
            if solv_res is not None:
                sm, sa = solv_res
                sk = Counter()
                for root, mem in comp.members.items():
                    info = comp.comp[root]
                    atoms = np.concatenate(
                        [aa[k] for k in mem["a"]] + [ba[k] for k in mem["b"]]
                    ) if (mem["a"] or mem["b"]) else np.empty(0, int)
                    s = _count_solvent(pos, atoms, sa, sm, fr.box, float(cutoff))
                    sk[(info["a"], info["b"], s)] += 1
                speckeys |= set(sk)
                spec_rows.append((fr.index, sk))

        warnings = list(dict.fromkeys(warnings))
        tables = []
        cols = [Field("frame", "frame")] + [Field(f"A{i}B{j}", "1") for (i, j) in _sorted_types(typekeys)]
        data = np.array([[fr] + [tc.get((i, j), 0) for (i, j) in _sorted_types(typekeys)]
                         for fr, tc in type_rows], float)
        definitions = {"cutoff": f"{cutoff:g} A", "mode": mode}
        tables.append(Table(
            title="Cluster species count",
            columns=cols, data=data, tag="agg",
            xlabel="frame",
            ylabel="count",
            definitions=definitions, warnings=warnings,
            meta={"module": "agg"},
            selection={"sel1": sel1, "sel2": sel2},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))

        lab_cols = [Field("frame", "frame"), Field("A_SSIP", "1"), Field("A_CIP", "1"),
                    Field("A_AGG", "1"), Field("B_SSIP", "1"), Field("B_CIP", "1"),
                    Field("B_AGG", "1")]
        tables.append(Table(
            title="Label fractions",
            columns=lab_cols, data=np.array(label_rows, float), tag="agg.labels",
            xlabel="frame",
            ylabel="fraction",
            definitions=definitions, warnings=warnings, meta={"module": "agg"},
            selection={"sel1": sel1, "sel2": sel2},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))

        if solv_res is not None:
            keys = sorted(speckeys, key=lambda k: (sum(k), k))
            srows = []
            for fr, sk in spec_rows:
                for k in keys:
                    srows.append((fr, k[0], k[1], k[2], sk.get(k, 0)))
            tables.append(Table(
                title="Speciation A_iB_jS_k",
                columns=[Field("frame", "frame"), Field("i", "1"), Field("j", "1"),
                         Field("k", "1"), Field("count", "1")],
                data=np.array(srows, float), tag="agg.speciation",
                xlabel="frame",
                ylabel="count",
                definitions=definitions, warnings=warnings,
                meta={"module": "agg"},
                selection={"sel1": sel1, "sel2": sel2, "solv": solv},
                frame_range=ctx.frame_range_str(frame_range, frames),
            ))
        return tables


def _count_solvent(pos, comp_atoms, solv_atoms, solv_mols, box, cutoff):
    from ..core.neighbor import pairs_within
    all_solv = np.concatenate(solv_atoms) if solv_atoms else np.empty(0, int)
    if len(all_solv) == 0 or len(comp_atoms) == 0:
        return 0
    ia, ib, _ = pairs_within(pos[comp_atoms], pos[all_solv], box, cutoff, same=False)
    if len(ib) == 0:
        return 0
    res_of = np.concatenate([np.full(len(a), m) for m, a in zip(solv_mols, solv_atoms)])
    return len(set(res_of[ib].tolist()))

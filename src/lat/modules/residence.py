"""lat residence — ligand residence times around a reference group."""

from __future__ import annotations

from collections import Counter

import numpy as np

from ..core.correlate import linear_fit
from ..core.errors import LatError
from ..core.neighbor import pairs_within
from ..core.periodic import periodic_distance, residue_com, residue_coms
from ..core.table import Field, Table
from ._graph import residues
from .base import Module, Param
from . import register


@register
class Residence(Module):
    name = "residence"
    summary = "ligand residence time (intermittent/continuous ACF + intervals)"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-ref", "-lig", "-cutoff"]
    outputs = ["<base>.residence.out.csv", "<base>.residence.acf.out.csv",
               "<base>.residence.intervals.out.csv"]
    params = [
        Param("ref", flags=("-ref",), required=True, help="reference group"),
        Param("lig", flags=("-lig",), required=True, help="ligand group"),
        Param("cutoff", flags=("-cutoff",), config="modules.residence.cutoff", type="float",
              required=True),
        Param("mode", flags=("-mode",), config="modules.residence.mode",
              choices=("min", "center")),
        Param("definition", flags=("-definition",), config="modules.residence.definition",
              choices=("cont", "inter", "both")),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, ref, lig, cutoff, mode, definition, frame_range, stride):
        topo = ctx.topo
        if cutoff is None:
            raise LatError("residence: -cutoff is required (no default)")
        rsel = ctx.selection(ref)
        if rsel.n_atoms == 0:
            raise LatError("residence: empty reference group")
        lm, la = residues(topo, lig and ctx.selection(lig).indices)
        n_lig = len(lm)
        if n_lig == 0:
            raise LatError("residence: no ligand residues")

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        ref_idx = rsel.indices
        lig_all = np.concatenate(la) if la else np.empty(0, int)
        res_of = np.concatenate([np.full(len(a), k) for k, a in enumerate(la)]) \
            if la else np.empty(0, int)

        present = np.zeros((n_lig, len(frames)), dtype=bool)
        occ_rows = []
        for c, fr in enumerate(ctx.frames(frames)):
            pos = fr.positions
            if mode == "center":
                com_l, ws = residue_coms(pos, fr.box, topo, la, what="ligand residue")
                com_r, wr = residue_com(pos, fr.box, topo, ref_idx, what="reference group")
                if wr and wr not in warnings:
                    warnings.append(wr)
                for w in ws:
                    if w not in warnings:
                        warnings.append(w)
                d = np.array([periodic_distance(c_l, com_r, fr.box) for c_l in com_l])
                present[d < float(cutoff), c] = True
            elif mode == "min":
                ia, ib, _ = pairs_within(pos[lig_all], pos[ref_idx], fr.box, float(cutoff))
                if len(ia):
                    present[np.unique(res_of[ia]), c] = True
            else:
                raise LatError(f"residence: -mode must be min|center, got {mode!r}")
            occ_rows.append((fr.index, int(present[:, c].sum())))

        if present.sum() == 0:
            warnings.append("no ligand bound under the given cutoff")

        C_inter = _acf_inter(present)
        C_cont = _acf_cont(present)
        tau_int = float(C_inter[1:].sum())
        y = np.clip(C_inter[1:], 1e-12, None)
        ok = y > 0
        if ok.sum() >= 3:
            slope, _, r2 = linear_fit(np.arange(1, len(C_inter))[ok], np.log(y[ok]))
            tau_fit = float(-1.0 / slope) if slope < 0 else float("nan")
        else:
            tau_fit, r2 = float("nan"), float("nan")
            warnings.append("too few points for lifetime fit; integral only")
        if len(frames) > 1 and tau_int > 0.5 * (len(frames) - 1):
            warnings.append("ACF not decayed within half the trajectory; lifetime is truncated")

        scalars = {"tau_integral_frames": tau_int, "tau_fit_frames": tau_fit, "r2": r2}
        dt = ctx.frame_interval
        if dt:
            scalars["tau_integral"] = f"{tau_int} x {dt} ps"
            scalars["tau_fit"] = f"{tau_fit} x {dt} ps"

        definitions = {"cutoff": f"{cutoff:g} A", "mode": mode}
        tables = [Table(
            title="Bound ligand count",
            columns=[Field("frame", "frame"), Field("n_bound", "1")],
            data=np.array(occ_rows, float), tag="residence",
            xlabel="frame",
            ylabel="n(bound)",
            definitions=definitions, warnings=warnings,
            scalars=scalars, meta={"module": "residence"},
            selection={"ref": ref, "lig": lig},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]
        acf = np.column_stack([np.arange(len(frames)), C_inter, C_cont])
        tables.append(Table(
            title="Residence ACF",
            columns=[Field("lag", "frame"), Field("C_inter", "1"), Field("C_cont", "1")],
            data=acf, tag="residence.acf",
            xlabel="lag (frame)",
            ylabel="C(t)",
            definitions=definitions, warnings=warnings, scalars=scalars,
            meta={"module": "residence"},
            selection={"ref": ref, "lig": lig},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))
        lengths = Counter(_run_lengths(present))
        if lengths:
            rows = [[k, v] for k, v in sorted(lengths.items())]
        else:
            rows = []
        tables.append(Table(
            title="Residence interval distribution",
            columns=[Field("interval_frames", "frame"), Field("count", "1")],
            data=np.array(rows, float) if rows else np.empty((0, 2)), tag="residence.intervals",
            xlabel="interval (frame)",
            ylabel="count",
            definitions=definitions, warnings=warnings, meta={"module": "residence"},
            selection={"ref": ref, "lig": lig},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))
        return tables


def _acf_inter(present):
    T = present.shape[1]
    total = present.sum()
    C = np.zeros(T)
    if total == 0:
        return C
    for lag in range(T):
        num = (present[:, : T - lag] & present[:, lag:]).sum()
        C[lag] = num * T / ((T - lag) * total)
    return C


def _acf_cont(present):
    T = present.shape[1]
    total = present.sum()
    C = np.zeros(T)
    if total == 0:
        return C
    C[0] = 1.0
    for lag in range(1, T):
        run = present[:, : T - lag].copy()
        for k in range(1, lag + 1):
            run &= present[:, k: T - lag + k]
        C[lag] = run.sum() * T / ((T - lag) * total)
    return C


def _run_lengths(present):
    out = []
    for j in range(present.shape[0]):
        series = present[j]
        run = 0
        for v in series:
            if v:
                run += 1
            elif run:
                out.append(run)
                run = 0
        if run:
            out.append(run)
    return out

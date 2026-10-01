"""lat rms — Kabsch-aligned RMSD per frame and RMSF per atom."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.table import Field, Table
from ..core.unwrap import build_adjacency, make_whole, molecule_groups
from .base import Module, Param
from . import register


def _kabsch(mobile, ref):
    """Row-vector rotation R minimizing |mobile @ R - ref| (inputs centered)."""
    h = mobile.T @ ref
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, d]) @ vt


@register
class RMS(Module):
    name = "rms"
    summary = "RMSD vs a reference frame and per-atom RMSF"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel"]
    outputs = ["<base>.rmsd.out.csv", "<base>.rmsf.out.csv"]
    params = [
        Param("sel", flags=("-sel",), required=True),
        Param("fit", flags=("-fit",), help="atoms used for alignment (default: -sel)"),
        Param("ref", flags=("-ref",), config="modules.rms.ref", help="first | average | frame index"),
        Param("rmsf", flags=("-rmsf",), config="modules.rms.rmsf", type="flag",
              help="also write per-atom RMSF"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel, fit, ref, rmsf, frame_range, stride):
        topo = ctx.topo
        selection = ctx.selection(sel)
        if selection.n_atoms == 0:
            raise LatError(f"selection {sel!r} matched 0 atoms")
        fit_sel = ctx.selection(fit) if fit else selection
        if fit_sel.n_atoms == 0:
            raise LatError(f"fit selection {fit!r} matched 0 atoms")

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        if len(frames) < 2:
            raise LatError("rms needs at least 2 frames")

        warnings: list[str] = []
        adj = build_adjacency(topo.n_atoms, topo.bonds)
        groups = None
        if adj is None:
            mols = topo.mols if topo.mols is not None else None
            if mols is None:
                warnings.append(
                    "no molecule information (mol/Bonds); RMSD may jump across the "
                    "periodic boundary for molecules that span the box"
                )
            else:
                groups = molecule_groups(mols)

        def _whole(fr):
            return make_whole(fr.positions, fr.box, adj=adj, groups=groups)

        # reference
        if ref in (None, "first", "0"):
            ref_index = frames[0]
        elif str(ref).isdigit():
            ref_index = int(ref)
        elif ref == "average":
            ref_index = None
        else:
            raise LatError(f"-ref must be first|average|<frame index>, got {ref!r}")

        si, fi = selection.indices, fit_sel.indices
        if ref_index is None:
            acc = []
            for fr in ctx.frames(frames):
                acc.append(_whole(fr))
            mean_all = np.mean(acc, axis=0)
            ref_calc = mean_all[si]
            ref_fit = mean_all[fi]
        else:
            if ref_index not in frames:
                raise LatError(f"reference frame {ref_index} not in the selected range")
            for fr in ctx.frames([ref_index]):
                pos = _whole(fr)
            ref_calc = pos[si]
            ref_fit = pos[fi]

        fit_center = ref_fit.mean(axis=0)
        ref_fit_c = ref_fit - fit_center

        rmsd_rows = []
        aligned = []
        for fr in ctx.frames(frames):
            pos = _whole(fr)
            mobile_fit_c = pos[fi] - pos[fi].mean(axis=0)
            rot = _kabsch(mobile_fit_c, ref_fit_c)
            calc_c = pos[si] - pos[fi].mean(axis=0)
            moved = calc_c @ rot + fit_center
            d = moved - ref_calc
            rmsd = float(np.sqrt((d * d).sum(axis=1).mean()))
            rmsd_rows.append((int(fr.index), rmsd))
            if rmsf:
                aligned.append(moved)

        rmsd_data = np.array(rmsd_rows, float)
        tables = [Table(
            title="RMSD",
            columns=[Field("frame", "frame"), Field("rmsd", "Angstrom")],
            data=rmsd_data,
            tag="rmsd",
            xlabel="frame",
            ylabel="RMSD (Angstrom)",
            definitions={"ref": str(ref), "fit": fit or selection.literal},
            warnings=warnings,
            meta={"module": "rms"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

        if rmsf:
            arr = np.stack(aligned)  # (T, Na, 3)
            mean = arr.mean(axis=0)
            fluct = np.sqrt(((arr - mean) ** 2).sum(axis=2).mean(axis=0))
            ids = topo.atom_id[si]
            tables.append(Table(
                title="RMSF",
                columns=[Field("atom", "id"), Field("rmsf", "Angstrom")],
                data=np.column_stack([ids, fluct]),
                tag="rmsf",
                xlabel="atom id",
                ylabel="RMSF (Angstrom)",
                warnings=warnings,
                meta={"module": "rms"},
                selection={"sel": selection.literal},
                frame_range=ctx.frame_range_str(frame_range, frames),
            ))
        else:
            warnings.append("RMSF not requested (-rmsf)")
        return tables

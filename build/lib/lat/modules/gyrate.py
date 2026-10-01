"""lat gyrate — radius of gyration of a selection."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.table import Field, Table
from ..core.unwrap import build_adjacency, make_whole, molecule_groups
from .base import Module, Param
from . import register


@register
class Gyrate(Module):
    name = "gyrate"
    summary = "radius of gyration Rg per frame"
    help_doc = "algorithms/density_msd.md"
    inputs = ["-s", "-sel"]
    outputs = ["<base>.gyrate.out.csv"]
    params = [
        Param("sel", flags=("-sel",), required=True),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel, frame_range, stride):
        topo = ctx.topo
        selection = ctx.selection(sel)
        if selection.n_atoms == 0:
            raise LatError(f"selection {sel!r} matched 0 atoms")
        idx = selection.indices

        warnings: list[str] = []
        try:
            masses = topo.atom_mass()[idx]
        except LatError:
            masses = np.ones(len(idx))
            warnings.append("no masses available; using uniform weights for Rg")

        adj = build_adjacency(topo.n_atoms, topo.bonds)
        groups = None
        if adj is None:
            mols = topo.mols[idx] if topo.mols is not None else None
            if mols is None:
                warnings.append(
                    "no molecule information (mol/Bonds); Rg may be underestimated "
                    "for molecules that span the periodic boundary"
                )
            else:
                groups = molecule_groups(mols)

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        rows = []
        for fr in ctx.frames(frames):
            if adj is not None:
                pos = make_whole(fr.positions, fr.box, adj=adj)[idx]
            else:
                pos = make_whole(fr.positions[idx], fr.box, groups=groups)
            total = masses.sum()
            com = (masses[:, None] * pos).sum(axis=0) / total
            d2 = ((pos - com) ** 2).sum(axis=1)
            rg = float(np.sqrt(np.sum(masses * d2) / total))
            rows.append((fr.index, rg))
        data = np.array(rows, float)
        return [Table(
            title=f"Radius of gyration {selection.literal}",
            columns=[Field("frame", "frame"), Field("rg", "Angstrom")],
            data=data,
            tag="gyrate",
            scalars={"rg_mean": float(data[:, 1].mean()),
                     "rg_min": float(data[:, 1].min()),
                     "rg_max": float(data[:, 1].max())},
            warnings=warnings,
            meta={"module": "gyrate"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

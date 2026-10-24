"""lat contact — inter-group contact counts and distance measures."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.neighbor import pairs_within
from ..core.periodic import periodic_distance, residue_com
from ..core.table import Field, Table
from .base import Module, Param
from . import register


@register
class Contact(Module):
    name = "contact"
    summary = "contact counts and minimum/COM distances between two groups"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel1", "-sel2", "-cutoff"]
    outputs = ["<base>.contact.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), required=True),
        Param("sel2", flags=("-sel2",), required=True),
        Param("cutoff", flags=("-cutoff",), config="modules.contact.cutoff", type="float",
              required=True),
        Param("mode", flags=("-mode",), config="modules.contact.mode", choices=("min", "center")),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel1, sel2, cutoff, mode, frame_range, stride):
        if cutoff is None:
            raise LatError("contact: -cutoff is required (no default)")
        s1 = ctx.selection(sel1)
        s2 = ctx.selection(sel2)
        if s1.n_atoms == 0 or s2.n_atoms == 0:
            raise LatError("contact: empty selection")

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        rows = []
        for fr in ctx.frames(frames):
            pos = fr.positions
            i, j, d = pairs_within(pos[s1.indices], pos[s2.indices], fr.box, float(cutoff))
            n_contact = len(d)
            min_d = float(d.min()) if n_contact else float("nan")
            com1, w1 = residue_com(pos, fr.box, ctx.topo, s1.indices,
                                   what=f"sel1={s1.literal}")
            com2, w2 = residue_com(pos, fr.box, ctx.topo, s2.indices,
                                   what=f"sel2={s2.literal}")
            for w in (w1, w2):
                if w and w not in warnings:
                    warnings.append(w)
            com_d = periodic_distance(com1, com2, fr.box)
            rows.append((fr.index, n_contact, min_d, com_d))
        data = np.array(rows, float)
        if np.all(np.isnan(data[:, 2])):
            warnings.append("no contacts within cutoff; minimum distance is undefined")
        return [Table(
            title="Contacts",
            columns=[Field("frame", "frame"), Field("n_contacts", "1"),
                     Field("min_dist", "Angstrom"), Field("com_dist", "Angstrom")],
            data=data,
            tag="contact",
            xlabel="frame",
            ylabel="count / distance",
            definitions={"cutoff": f"{cutoff:g} A"},
            scalars={"mean_contacts": float(np.nanmean(data[:, 1]))},
            warnings=warnings,
            meta={"module": "contact"},
            selection={"sel1": s1.literal, "sel2": s2.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

"""lat contact — inter-group contact counts and distance measures."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.neighbor import pairs_within
from ..core.periodic import periodic_distance, residue_com, residue_coms
from ..core.table import Field, Table
from ._graph import residues
from .base import Module, Param
from . import register


@register
class Contact(Module):
    name = "contact"
    summary = "contact counts and minimum/COM distances between two groups"
    help_doc = "manual/contact.md"
    inputs = ["-s", "-sel1", "-sel2", "-cutoff"]
    outputs = ["<base>.contact.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), required=True, help="first group"),
        Param("sel2", flags=("-sel2",), required=True, help="second group"),
        Param("cutoff", flags=("-cutoff",), config="modules.contact.cutoff", type="float",
              required=True, help="contact cutoff in Angstrom (required; no default)"),
        Param("mode", flags=("-mode",), config="modules.contact.mode", choices=("min", "center"),
              help="count atom pairs (min) or residue pairs by COM (center)"),
        Param("frame_range", flags=("-frame_range",), help="frame window start:stop[:stride]"),
        Param("stride", flags=("-stride",), type="int", help="sample every N frames (default: 1)"),
    ]

    def run(self, ctx, *, sel1, sel2, cutoff, mode, frame_range, stride):
        if cutoff is None:
            raise LatError("contact: -cutoff is required (no default)")
        if mode not in ("min", "center"):
            raise LatError(f"contact: -mode must be min|center, got {mode!r}")
        s1 = ctx.selection(sel1)
        s2 = ctx.selection(sel2)
        if s1.n_atoms == 0 or s2.n_atoms == 0:
            raise LatError("contact: empty selection")

        same = (s1.n_atoms == s2.n_atoms
                and np.array_equal(np.sort(s1.indices), np.sort(s2.indices)))
        # For -mode center the contact graph is built on residue centres; the
        # residue partition is topology-only, so it is computed once.
        groups1 = groups2 = None
        if mode == "center":
            _, groups1 = residues(ctx.topo, s1.indices)
            groups2 = groups1 if same else residues(ctx.topo, s2.indices)[1]

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        rows = []
        for fr in ctx.frames(frames):
            pos = fr.positions
            if mode == "center":
                c1, w1 = residue_coms(pos, fr.box, ctx.topo, groups1,
                                      what=f"sel1={s1.literal}")
                c2 = c1 if same else residue_coms(
                    pos, fr.box, ctx.topo, groups2, what=f"sel2={s2.literal}")[0]
                i, j, d = pairs_within(c1, c2, fr.box, float(cutoff), same=same)
                n_contact = len(d)
                min_d = float(d.min()) if n_contact else float("nan")
                com1, cw1 = residue_com(pos, fr.box, ctx.topo, s1.indices,
                                        what=f"sel1={s1.literal}")
                com2, cw2 = residue_com(pos, fr.box, ctx.topo, s2.indices,
                                        what=f"sel2={s2.literal}")
                for w in (w1, cw1, cw2):
                    for item in (w if isinstance(w, list) else [w]):
                        if item and item not in warnings:
                            warnings.append(item)
            else:
                i, j, d = pairs_within(pos[s1.indices], pos[s2.indices], fr.box,
                                       float(cutoff), same=same)
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
            definitions={"cutoff": f"{cutoff:g} A",
                         "mode": "atom pairs (min)" if mode == "min" else "residue COM pairs (center)"},
            scalars={"mean_contacts": float(np.nanmean(data[:, 1]))},
            warnings=warnings,
            meta={"module": "contact"},
            selection={"sel1": s1.literal, "sel2": s2.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

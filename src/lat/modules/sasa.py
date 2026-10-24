"""lat sasa — solvent accessible surface area (Shrake-Rupley, periodic).

The surface is always built from *all* atoms, so a selected atom is shielded by
its neighbours even when they are outside the selection; `-sel` only chooses
whose contribution is reported separately.

All neighbour search and sphere-point tests are minimum-image (see
lat.core.periodic.sasa_shrake_rupley), so atoms near a cell face are treated
exactly like atoms in the middle of the cell.
"""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.periodic import sasa_shrake_rupley
from ..core.table import Field, Table
from ..core.topology import vdw_radii
from .base import Module, Param
from . import register


@register
class Sasa(Module):
    name = "sasa"
    summary = "solvent accessible surface area (Shrake-Rupley, PBC-aware)"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel", "-radius", "-sasa_points"]
    outputs = ["<base>.sasa.out.csv"]
    params = [
        Param("radius", flags=("-radius",), config="modules.sasa.radius", type="float"),
        Param("sasa_points", flags=("-sasa_points",), config="modules.sasa.sasa_points",
              type="int"),
        Param("sel", flags=("-sel",), help="atoms whose area is reported (surface uses all)"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, radius, sasa_points, sel, frame_range, stride):
        if radius is None or sasa_points is None:
            raise LatError("sasa: -radius and -sasa_points need defaults (check config)")
        selection = ctx.selection(sel or "all")
        idx = selection.indices
        if len(idx) == 0:
            raise LatError(f"sasa: selection {(sel or 'all')!r} matched 0 atoms")

        radii = vdw_radii(ctx.topo, ctx.cfg, ctx.element_map())
        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        rows = []
        for fr in ctx.frames(frames):
            total, per_atom = sasa_shrake_rupley(fr.positions, radii, float(radius), fr.box,
                                                 int(sasa_points))
            rows.append((fr.index, total, float(per_atom[idx].sum())))
        data = np.array(rows, float)
        definitions = {
            "radius": f"{radius} A",
            "sasa_points": int(sasa_points),
            "surface": "all atoms (shielding is always computed over the whole system)",
        }
        return [Table(
            title="SASA",
            columns=[Field("frame", "frame"), Field("sasa_total", "Angstrom^2"),
                     Field("sasa_sel", "Angstrom^2")],
            data=data,
            tag="sasa",
            xlabel="frame",
            ylabel="SASA (Angstrom^2)",
            definitions=definitions,
            scalars={"sasa_mean": float(data[:, 1].mean()),
                     "sasa_sel_mean": float(data[:, 2].mean())},
            warnings=warnings,
            meta={"module": "sasa"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

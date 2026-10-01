"""lat pistack — pi-pi stacking via ring centers and normals.

Ported from /mnt/d/dev/83_pipi_stack_code (geometry.py + analyzer.py).
Ring definition is explicit: each residue in the selection contributes one ring
(use -sel1/-sel2 to pick the aromatic atoms).
"""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.neighbor import pairs_within
from ..core.table import Field, Table
from ._graph import residues
from .base import Module, Param
from . import register


def _ring_geometry(pos, atoms_by_ring):
    centers = np.empty((len(atoms_by_ring), 3))
    normals = np.empty((len(atoms_by_ring), 3))
    for k, atoms in enumerate(atoms_by_ring):
        c = pos[atoms]
        centers[k] = c.mean(axis=0)
        centered = c - centers[k]
        if len(c) >= 3:
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            n = vt[-1]
            normals[k] = n / (np.linalg.norm(n) + 1e-12)
        else:
            normals[k] = np.array([0.0, 0.0, 1.0])
    return centers, normals


@register
class PiStack(Module):
    name = "pistack"
    summary = "pi-pi stacking count and per-pair geometry"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel1", "-sel2"]
    outputs = ["<base>.pistack.out.csv", "<base>.pistack.pairs.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), required=True, help="ring atoms of group A"),
        Param("sel2", flags=("-sel2",), help="ring atoms of group B (default: self-stacking)"),
        Param("distance", flags=("-distance",), config="modules.pistack.distance", type="float"),
        Param("angle", flags=("-angle",), config="modules.pistack.angle", type="float"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel1, sel2, distance, angle, frame_range, stride):
        topo = ctx.topo
        s1 = ctx.selection(sel1)
        _, rings_a = residues(topo, s1.indices)
        rings_a = [a for a in rings_a if len(a) >= 3]
        if not rings_a:
            raise LatError("pistack: sel1 must select >=3 atoms per residue (one ring per residue)")
        self_mode = sel2 is None
        if self_mode:
            rings_b = rings_a
        else:
            s2 = ctx.selection(sel2)
            _, rings_b = residues(topo, s2.indices)
            rings_b = [a for a in rings_b if len(a) >= 3]
            if not rings_b:
                raise LatError("pistack: sel2 must select >=3 atoms per residue")

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        if len(rings_a) < 2 and self_mode:
            warnings.append("fewer than 2 rings in sel1; no self-stacking possible")

        counts = []
        pair_rows = []
        for fr in ctx.frames(frames):
            pos = fr.positions
            ca, na = _ring_geometry(pos, rings_a)
            cb, nb = _ring_geometry(pos, rings_b)
            if self_mode:
                ia, ib, d = pairs_within(ca, None, fr.box, float(distance), same=True)
            else:
                ia, ib, d = pairs_within(ca, cb, fr.box, float(distance), same=False)
            n = 0
            for k in range(len(ia)):
                a, b = int(ia[k]), int(ib[k])
                dot = abs(float(np.dot(na[a], nb[b])))
                ang = float(np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))))
                if ang <= float(angle):
                    n += 1
                    pair_rows.append((fr.index, a, b, float(d[k]), ang))
            counts.append((fr.index, n))

        definitions = {"distance": f"{distance} A", "angle": f"{angle} deg"}
        tables = [Table(
            title="Pi-pi stacking",
            columns=[Field("frame", "frame"), Field("n_pistack", "1")],
            data=np.array(counts, float), tag="pistack",
            xlabel="frame",
            ylabel="n(pistack)",
            definitions=definitions,
            scalars={"mean_n_pistack": float(np.mean([c for _, c in counts]))},
            warnings=warnings, meta={"module": "pistack"},
            selection={"sel1": sel1, "sel2": sel2 or "self"},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]
        pair_data = np.array(pair_rows, float) if pair_rows else np.empty((0, 5))
        tables.append(Table(
            title="Pi-pi stacking pairs",
            columns=[Field("frame", "frame"), Field("ring_i", "1"), Field("ring_j", "1"),
                     Field("dist", "Angstrom"), Field("angle", "deg")],
            data=pair_data, tag="pistack.pairs",
            xlabel="frame",
            ylabel="pair count / geometry",
            definitions=definitions, warnings=warnings, meta={"module": "pistack"},
            selection={"sel1": sel1, "sel2": sel2 or "self"},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))
        if not pair_rows:
            warnings.append("no pi-pi stacking pairs found under the given criteria")
        return tables

"""lat select — inspect/export atom selections (screen output, optional ndx)."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.selection import write_ndx
from .base import Module, Param
from . import register


@register
class Select(Module):
    name = "select"
    summary = "inspect a selection and/or export it to a GROMACS .ndx"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel", "-export-ndx", "-name"]
    outputs = [".ndx (optional)"]
    screen_only = True
    params = [
        Param("sel", flags=("-sel",), help="selection expression"),
        Param("export_ndx", flags=("-export-ndx",), help="write selection to this .ndx"),
        Param("name", flags=("-name",), help="group name in the exported ndx"),
    ]

    def run(self, ctx, *, sel, export_ndx, name):
        topo = ctx.topo
        print(f"topology: {topo.source}")
        print(
            f"  atoms={topo.n_atoms}  types={len(topo.unique_types)}  "
            f"bonds={0 if topo.bonds is None else len(topo.bonds)}  frames={ctx.n_frames}"
        )
        print(f"  atom_style={topo.meta.get('style', 'n/a')}  box={np.round(topo.box.lengths, 4).tolist()}")

        if sel is None:
            for t in topo.unique_types:
                m = topo.mask_of_types([int(t)])
                el = topo.type_element.get(int(t), "?")
                mass = topo.type_mass.get(int(t), float("nan"))
                print(f"  type {int(t):>4}  element={el:<3}  count={int(m.sum()):>8}  mass={mass:.4f}")
            return []

        selection = ctx.selection(sel)
        idx = selection.indices
        if selection.n_atoms == 0:
            raise LatError(f"selection {sel!r} matched 0 atoms")
        types, counts = np.unique(topo.types[idx], return_counts=True)
        print(f"selection: {sel}")
        print(f"  matched atoms: {selection.n_atoms}")
        print("  type distribution:")
        for t, c in zip(types, counts):
            el = topo.type_element.get(int(t), "?")
            print(f"    type {int(t):>4}  element={str(el or '?'):<3}  count={int(c):>8}")
        if topo.charges is not None:
            print(f"  charge sum: {float(np.sum(topo.charges[idx])):.6g}")
        if topo.mols is not None:
            print(f"  molecules:  {len(np.unique(topo.mols[idx]))}")
        if topo.resnames is not None:
            rn = sorted(set(map(str, topo.resnames[idx])))
            print(f"  resnames:   {', '.join(rn)}")

        if export_ndx:
            group = name or "lat"
            write_ndx(export_ndx, {group: (idx + 1).tolist()})
            print(f"  wrote {export_ndx} [{group}] with {len(idx)} atoms")
        return []

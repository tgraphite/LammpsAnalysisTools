"""lat density — 1D / 2D density profile of a selection."""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.grid import smooth as smooth_fn
from ..core.table import Field, Table
from .base import Module, Param
from . import register

_AXES = {"x": 0, "y": 1, "z": 2}
_MASS_FACTOR = 1e27 / 6.02214076e23  # (g/mol)/A^3 -> kg/m^3


@register
class Density(Module):
    name = "density"
    summary = "density profile along an axis (1D) or in a plane (2D)"
    help_doc = "algorithms/density_msd.md"
    inputs = ["-s", "-sel", "-axis"]
    outputs = ["<base>.density.out.csv", "<base>.density.2d.out.csv"]
    params = [
        Param("sel", flags=("-sel",), help="selection (default: all)"),
        Param("axis", flags=("-axis",), config="modules.density.axis", choices=("x", "y", "z")),
        Param("nbins", flags=("-nbins",), config="modules.density.nbins", type="int"),
        Param("dim", flags=("-dim",), config="modules.density.dim", type="int", choices=(1, 2)),
        Param("quantity", flags=("-quantity",), config="modules.density.quantity",
              choices=("mass", "number")),
        Param("smooth", flags=("-smooth",), config="modules.density.smooth",
              choices=("none", "savgol", "gaussian")),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel, axis, nbins, dim, quantity, smooth, frame_range, stride):
        topo = ctx.topo
        if axis not in _AXES:
            raise LatError(f"-axis must be x/y/z, got {axis!r}")
        if int(dim) not in (1, 2):
            raise LatError(f"-dim must be 1 or 2, got {dim}")
        if quantity not in ("mass", "number"):
            raise LatError(f"-quantity must be mass|number, got {quantity!r}")

        selection = ctx.selection(sel or "all")
        if selection.n_atoms == 0:
            raise LatError(f"selection {sel!r} matched 0 atoms")
        idx = selection.indices
        weights = None
        if quantity == "mass":
            weights = topo.atom_mass()[idx]

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        nframes = len(frames)
        ax = _AXES[axis]

        # geometry from the first selected frame
        first = next(iter(ctx.frames(frames)))
        box = first.box
        lo, hi = box.lo[ax], box.hi[ax]
        nbins = int(nbins)

        warnings: list[str] = []
        if hi - lo <= 0:
            raise LatError(f"degenerate axis {axis} (length {hi - lo})")

        if int(dim) == 1:
            value = self._profile_1d(ctx, frames, idx, weights, box, ax, lo, hi, nbins, quantity)
            sm = self._grid_smooth(ctx, smooth)
            value = smooth_fn(value, smooth, **sm) if smooth and smooth != "none" else value
            centers = 0.5 * (np.linspace(lo, hi, nbins + 1)[:-1] + np.linspace(lo, hi, nbins + 1)[1:])
            unit = "kg/m^3" if quantity == "mass" else "1/Angstrom^3"
            label = "rho" if quantity == "mass" else "n"
            return [Table(
                title=f"Density profile ({quantity}) {selection.literal} along {axis}",
                columns=[Field(axis, "Angstrom"), Field(f"{label}_{axis}", unit)],
                data=np.column_stack([centers, value]),
                tag="density",
                warnings=warnings,
                meta={"module": "density"},
                selection={"sel": selection.literal, "quantity": quantity},
                frame_range=ctx.frame_range_str(frame_range, frames),
            )]

        # 2D long table
        other = [i for i in (0, 1, 2) if i != ax]
        a, b = other
        loa, hia = box.lo[a], box.hi[a]
        lob, hib = box.lo[b], box.hi[b]
        edges_a = np.linspace(loa, hia, nbins + 1)
        edges_b = np.linspace(lob, hib, nbins + 1)
        da = (hia - loa) / nbins
        db = (hib - lob) / nbins
        acc = np.zeros((nbins, nbins))
        for fr in ctx.frames(frames):
            pos = fr.positions[idx]
            h, _, _ = np.histogram2d(pos[:, a], pos[:, b], bins=[edges_a, edges_b], weights=weights)
            cross_f = fr.box.volume / fr.box.lengths[ax]
            acc += h / (da * db * cross_f)
        acc /= nframes
        if quantity == "mass":
            acc = acc * _MASS_FACTOR
            unit = "kg/m^3"
        else:
            unit = "1/Angstrom^3"
        ca = 0.5 * (edges_a[:-1] + edges_a[1:])
        cb = 0.5 * (edges_b[:-1] + edges_b[1:])
        aa, bb = np.meshgrid(ca, cb, indexing="ij")
        data = np.column_stack([aa.ravel(), bb.ravel(), acc.ravel()])
        names = ("xyz"[a], "xyz"[b])
        return [Table(
            title=f"2D density ({quantity}) {selection.literal} in {names[0]}{names[1]}",
            columns=[Field(names[0], "Angstrom"), Field(names[1], "Angstrom"),
                     Field("rho", unit)],
            data=data,
            tag="density",
            mode="2d",
            kind="grid",
            warnings=warnings,
            meta={"module": "density"},
            selection={"sel": selection.literal, "quantity": quantity},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

    @staticmethod
    def _profile_1d(ctx, frames, idx, weights, box, ax, lo, hi, nbins, quantity):
        edges = np.linspace(lo, hi, nbins + 1)
        dr = (hi - lo) / nbins
        acc = np.zeros(nbins)
        for fr in ctx.frames(frames):
            coords = fr.positions[idx, ax]
            h, _ = np.histogram(coords, bins=edges, weights=weights)
            cross_f = fr.box.volume / fr.box.lengths[ax]
            acc += h / (cross_f * dr)
        acc /= len(frames)
        if quantity == "mass":
            acc = acc * _MASS_FACTOR
        return acc

    @staticmethod
    def _grid_smooth(ctx, method):
        g = ctx.cfg.get("grid", default={}) or {}
        return {
            "window": int(g.get("savgol_window", 5)),
            "order": int(g.get("savgol_order", 2)),
            "sigma": float(g.get("gaussian_sigma", 1.0)),
        }

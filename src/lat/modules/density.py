"""lat density — 1D / 2D density profile of a selection.

Binning follows the gmx density convention: positions are binned in *fractional*
(box-relative) coordinates each frame, so the bins track a changing (NPT) box;
the reported axis is scaled back to the mean box length.
"""

from __future__ import annotations

import numpy as np

from ..core.errors import LatError
from ..core.grid import box_change_warning, fractional_bins
from ..core.table import Field, Table
from .base import Module, Param
from . import register

_AXES = {"x": 0, "y": 1, "z": 2}
_MASS_FACTOR = 1e27 / 6.02214076e23  # (g/mol)/A^3 -> kg/m^3


@register
class Density(Module):
    name = "density"
    summary = "density profile along an axis (1D) or in a plane (2D)"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel", "-axis"]
    outputs = ["<base>.density.out.csv", "<base>.density.2d.out.csv"]
    params = [
        Param("sel", flags=("-sel",), help="selection (default: all)"),
        Param("axis", flags=("-axis",), config="modules.density.axis", choices=("x", "y", "z")),
        Param("nbins", flags=("-nbins",), config="modules.density.nbins", type="int"),
        Param("dim", flags=("-dim",), config="modules.density.dim", type="int", choices=(1, 2)),
        Param("quantity", flags=("-quantity",), config="modules.density.quantity",
              choices=("mass", "number")),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel, axis, nbins, dim, quantity, frame_range, stride):
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
        weights = topo.atom_mass()[idx] if quantity == "mass" else None

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        nbins = int(nbins)
        ax = _AXES[axis]
        unit = "kg/m^3" if quantity == "mass" else "1/Angstrom^3"
        fr_str = ctx.frame_range_str(frame_range, frames)
        definitions = {"quantity": quantity, "nbins": nbins, "dim": int(dim)}

        if int(dim) == 1:
            value, centers, lengths = self._profile_1d(ctx, frames, idx, weights, ax, nbins, quantity)
            warnings = []
            w = box_change_warning(lengths, [ax])
            if w:
                warnings.append(w)
            label = "rho" if quantity == "mass" else "n"
            ylabel = ("density (kg/m^3)" if quantity == "mass"
                      else "number density (1/Angstrom^3)")
            definitions["axis"] = axis
            return [Table(
                title="Density profile",
                columns=[Field(axis, "Angstrom"), Field(f"{label}_{axis}", unit)],
                data=np.column_stack([centers, value]),
                tag="density",
                ylabel=ylabel,
                definitions=definitions,
                warnings=warnings,
                meta={"module": "density"},
                selection={"sel": selection.literal},
                frame_range=fr_str,
            )]

        other = [i for i in (0, 1, 2) if i != ax]
        a, b = other
        acc, ca, cb, lengths = self._profile_2d(ctx, frames, idx, weights, a, b, nbins, quantity)
        warnings = []
        w = box_change_warning(lengths, [a, b])
        if w:
            warnings.append(w)
        names = ("xyz"[a], "xyz"[b])
        aa, bb = np.meshgrid(ca, cb, indexing="ij")
        data = np.column_stack([aa.ravel(), bb.ravel(), acc.ravel()])
        return [Table(
            title="2D density",
            columns=[Field(names[0], "Angstrom"), Field(names[1], "Angstrom"),
                     Field("rho", unit)],
            data=data,
            tag="density",
            mode="2d",
            kind="grid",
            definitions=definitions,
            warnings=warnings,
            meta={"module": "density"},
            selection={"sel": selection.literal},
            frame_range=fr_str,
        )]

    @staticmethod
    def _profile_1d(ctx, frames, idx, weights, ax, nbins, quantity):
        edges, centers = fractional_bins(nbins)
        acc = np.zeros(nbins)
        lengths = np.zeros((len(frames), 3))
        origin_sum = np.zeros(3)
        n = 0
        for fr in ctx.frames(frames):
            L = fr.box.lengths
            lo = fr.box.origin
            s = (fr.positions[idx, ax] - lo[ax]) / L[ax]
            s -= np.floor(s)                       # wrap fractional coord into [0, 1)
            h, _ = np.histogram(s, bins=edges, weights=weights)
            acc += h * nbins / fr.box.volume       # fractional slab volume = V/nbins
            lengths[n] = L
            origin_sum += lo
            n += 1
        acc /= max(n, 1)
        mean_L = lengths[:n].mean(axis=0)
        mean_lo = origin_sum / max(n, 1)
        x = mean_lo[ax] + centers * mean_L[ax]
        if quantity == "mass":
            acc = acc * _MASS_FACTOR
        return acc, x, lengths[:n]

    @staticmethod
    def _profile_2d(ctx, frames, idx, weights, a, b, nbins, quantity):
        edges, centers = fractional_bins(nbins)
        acc = np.zeros((nbins, nbins))
        lengths = np.zeros((len(frames), 3))
        origin_sum = np.zeros(3)
        n = 0
        for fr in ctx.frames(frames):
            L = fr.box.lengths
            lo = fr.box.origin
            sa = (fr.positions[idx, a] - lo[a]) / L[a]
            sb = (fr.positions[idx, b] - lo[b]) / L[b]
            sa -= np.floor(sa)
            sb -= np.floor(sb)
            h, _, _ = np.histogram2d(sa, sb, bins=[edges, edges], weights=weights)
            acc += h * (nbins * nbins) / fr.box.volume
            lengths[n] = L
            origin_sum += lo
            n += 1
        acc /= max(n, 1)
        mean_L = lengths[:n].mean(axis=0)
        mean_lo = origin_sum / max(n, 1)
        ca = mean_lo[a] + centers * mean_L[a]
        cb = mean_lo[b] + centers * mean_L[b]
        if quantity == "mass":
            acc = acc * _MASS_FACTOR
        return acc, ca, cb, lengths[:n]

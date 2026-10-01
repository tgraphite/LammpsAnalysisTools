"""lat probe — free volume, SASA, PSD, porosity and connectivity.

Shared core: probe sphere + vdW spheres + spatial occupancy.
PSD/porosity use a probe-grid approximation (declared); SASA is Shrake-Rupley.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from ..core.errors import LatError
from ..core.table import Field, Table
from .base import Module, Param
from . import register


def _vdw_radii(topo, cfg, element_map=None):
    table = cfg.get("vdw", default={}) or {}
    emap = element_map if element_map is not None else topo.type_element
    radii = np.empty(topo.n_atoms)
    missing = []
    for i, t in enumerate(topo.types):
        el = emap.get(int(t))
        if el is None or el not in table:
            missing.append(el or f"type {int(t)}")
            radii[i] = np.nan
        else:
            radii[i] = float(table[el])
    if missing:
        raise LatError(
            f"probe: no vdW radius for {sorted(set(missing))}; provide element "
            f"identity (Masses section, an `element` dump column, or -types file) "
            f"or add the element to config vdw"
        )
    return radii


def _grid_axes(box, dr):
    axes = []
    for k in range(3):
        lo, hi = box.lo[k], box.hi[k]
        n = max(int(round((hi - lo) / dr)), 1)
        c = lo + (np.arange(n) + 0.5) * (hi - lo) / n
        axes.append(c)
    return axes


def _occupancy(pos, box, axes, radii, dr):
    """Mark grid voxels within (vdW + probe) of any atom, periodic in every axis.

    Vectorized per atom over its local voxel block: O(sum (2r/dr+1)^3) with
    numpy, instead of a Python triple loop over every voxel.
    """
    dims = np.array([len(a) for a in axes], dtype=int)
    occ = np.zeros(tuple(int(d) for d in dims), dtype=bool)
    L = box.lengths
    dr = float(dr)
    pos = box.wrap(np.asarray(pos, float))
    los = np.array([a[0] for a in axes])

    def _axis(k, x):
        i0 = int(np.floor((x - radii_atom - los[k]) / dr))
        i1 = int(np.ceil((x + radii_atom - los[k]) / dr))
        ii = np.arange(i0, i1 + 1)
        iw = np.mod(ii, dims[k])
        ax = axes[k][iw]
        d = x - ax
        d -= L[k] * np.round(d / L[k])  # minimum image
        d = np.abs(d)
        if len(np.unique(iw)) != len(iw):
            uw, inv = np.unique(iw, return_inverse=True)
            dm = np.full(len(uw), np.inf)
            np.minimum.at(dm, inv, d)
            return uw, dm
        return iw, d

    for atom in range(len(pos)):
        radii_atom = radii[atom]
        if not np.isfinite(radii_atom) or radii_atom <= 0:
            continue
        r = radii_atom
        ux, dx = _axis(0, pos[atom, 0])
        uy, dy = _axis(1, pos[atom, 1])
        uz, dz = _axis(2, pos[atom, 2])
        mask = (dx[:, None, None] ** 2 + dy[None, :, None] ** 2
                + dz[None, None, :] ** 2) <= r * r
        occ[np.ix_(ux, uy, uz)] |= mask
    return occ


def _sasa(pos, radii, probe, n_points=92):
    if not hasattr(_sasa, "_dirs"):
        k = np.arange(n_points) + 0.5
        phi = np.arccos(1 - 2 * k / n_points)
        theta = np.pi * (1 + 5 ** 0.5) * k
        _sasa._dirs = np.column_stack([
            np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)
        ])
    dirs = _sasa._dirs
    rr = radii + probe
    tree = cKDTree(pos)
    total = 0.0
    per_atom = np.zeros(len(pos))
    for i in range(len(pos)):
        Ri = rr[i]
        neigh = tree.query_ball_point(pos[i], Ri + rr.max())
        neigh = [j for j in neigh if j != i]
        pts = pos[i] + dirs * Ri
        exposed = np.ones(len(dirs), dtype=bool)
        for j in neigh:
            d2 = np.sum((pts - pos[j]) ** 2, axis=1)
            exposed &= d2 > rr[j] ** 2
        area = 4 * np.pi * Ri ** 2 * exposed.mean()
        per_atom[i] = area
        total += area
    return total, per_atom


@register
class Probe(Module):
    name = "probe"
    summary = "free volume / porosity / SASA / PSD / connectivity"
    help_doc = "algorithms/probe.md"
    inputs = ["-s", "-radius", "-mode", "-grid_dr"]
    outputs = ["<base>.probe.<mode>.out.csv"]
    params = [
        Param("radius", flags=("-radius",), config="modules.probe.radius", type="float"),
        Param("mode", flags=("-mode",), config="modules.probe.mode",
              choices=("ffv", "psd", "porosity", "connectivity", "sasa", "all")),
        Param("grid_dr", flags=("-grid_dr",), config="modules.probe.grid_dr", type="float"),
        Param("connectivity", flags=("-connectivity",), config="modules.probe.connectivity",
              type="int"),
        Param("sel", flags=("-sel",), help="for SASA: subset (default: all)"),
        Param("sasa_points", flags=("-sasa_points",), config="modules.probe.sasa_points",
              type="int"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, radius, mode, grid_dr, connectivity, sel, sasa_points,
            frame_range, stride):
        if mode not in ("ffv", "psd", "porosity", "connectivity", "sasa", "all"):
            raise LatError(f"probe: unknown -mode {mode!r}")
        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        tables = []

        want_grid = mode in ("ffv", "porosity", "psd", "connectivity", "all")
        want_sasa = mode in ("sasa", "all")

        if want_grid:
            radii = _vdw_radii(ctx.topo, ctx.cfg, ctx.element_map()) + float(radius)
            ffv_rows = []
            psd_acc = None
            conn_rows = []
            first = next(iter(ctx.frames(frames)))
            dr = float(grid_dr)
            axes = _grid_axes(first.box, dr)
            ntotal = len(axes[0]) * len(axes[1]) * len(axes[2])
            cap = 5e7
            if ntotal > cap:
                # coarsen automatically rather than failing on large boxes
                dr = dr * (ntotal / cap) ** (1.0 / 3.0)
                axes = _grid_axes(first.box, dr)
                ntotal = len(axes[0]) * len(axes[1]) * len(axes[2])
                warnings.append(
                    f"probe grid would need {ntotal} points; coarsened grid_dr "
                    f"{float(grid_dr):.3g} -> {dr:.3g} A to stay under {cap:.0e}"
                )
            for fr in ctx.frames(frames):
                occ = _occupancy(fr.positions, fr.box, axes, radii, dr)
                free = ~occ
                ffv = float(free.mean())
                ffv_rows.append((fr.index, ffv))
                if mode in ("psd", "all") and psd_acc is None:
                    dt = ndimage.distance_transform_edt(free, sampling=dr)
                    r_pore = dt - float(radius)
                    local = (r_pore == ndimage.maximum_filter(r_pore, size=3)) & (r_pore > 0)
                    vals = r_pore[local]
                    voxel = dr ** 3
                    total_free = free.sum() * voxel
                    hist, edges = np.histogram(vals, bins=40,
                                               weights=np.full(vals.size, voxel))
                    rmid = 0.5 * (edges[:-1] + edges[1:])
                    psd = hist / total_free / (edges[1] - edges[0])
                    psd_acc = np.column_stack([rmid, psd])
                if mode in ("connectivity", "all"):
                    structure = ndimage.generate_binary_structure(
                        3, 1 if int(connectivity) == 6 else 3)
                    lab, ncomp = ndimage.label(free, structure=structure)
                    for cid in range(1, ncomp + 1):
                        vol = float((lab == cid).sum()) * dr ** 3
                        percolating = _percolating(lab == cid)
                        conn_rows.append((fr.index, cid, vol, int(percolating)))

            ffv_data = np.array(ffv_rows, float)
            ffv_table = Table(
                title=f"Free volume fraction (probe {radius} A)",
                columns=[Field("frame", "frame"), Field("ffv", "1")],
                data=ffv_data, tag="probe.ffv",
                scalars={"ffv_mean": float(ffv_data[:, 1].mean()),
                         "porosity_mean": float(ffv_data[:, 1].mean())},
                warnings=warnings, meta={"module": "probe"},
                selection={"radius": f"{radius} A", "grid_dr": f"{dr} A"},
                frame_range=ctx.frame_range_str(frame_range, frames),
            )
            tables.append(ffv_table)
            if psd_acc is not None:
                tables.append(Table(
                    title="Pore size distribution (probe-grid approximation)",
                    columns=[Field("r", "Angstrom"), Field("dv_dr", "1/Angstrom")],
                    data=psd_acc, tag="probe.psd", warnings=warnings,
                    meta={"module": "probe"},
                    selection={"radius": f"{radius} A", "grid_dr": f"{dr} A"},
                    frame_range=ctx.frame_range_str(frame_range, frames),
                ))
            if conn_rows:
                tables.append(Table(
                    title="Free-space connectivity",
                    columns=[Field("frame", "frame"), Field("component", "1"),
                             Field("volume", "Angstrom^3"), Field("percolating", "1")],
                    data=np.array(conn_rows, float), tag="probe.connectivity",
                    warnings=warnings, meta={"module": "probe"},
                    selection={"mode": "connectivity", "grid_dr": f"{dr} A",
                               "connectivity": int(connectivity)},
                    frame_range=ctx.frame_range_str(frame_range, frames),
                ))

        if want_sasa:
            selection = ctx.selection(sel or "all")
            idx = selection.indices
            radii = _vdw_radii(ctx.topo, ctx.cfg, ctx.element_map())[idx]
            rows = []
            for fr in ctx.frames(frames):
                total, _ = _sasa(fr.positions[idx], radii, float(radius), int(sasa_points))
                rows.append((fr.index, total))
            data = np.array(rows, float)
            tables.append(Table(
                title=f"SASA {selection.literal} (probe {radius} A)",
                columns=[Field("frame", "frame"), Field("sasa", "Angstrom^2")],
                data=data, tag="probe.sasa",
                scalars={"sasa_mean": float(data[:, 1].mean())},
                warnings=warnings, meta={"module": "probe"},
                selection={"sel": selection.literal, "radius": f"{radius} A"},
                frame_range=ctx.frame_range_str(frame_range, frames),
            ))

        if not tables:
            raise LatError(f"probe: mode {mode!r} produced no output")
        return tables


def _percolating(mask):
    b = np.zeros(6, bool)
    b[0] = mask[0].any(); b[1] = mask[-1].any()
    b[2] = mask[:, 0].any(); b[3] = mask[:, -1].any()
    b[4] = mask[:, :, 0].any(); b[5] = mask[:, :, -1].any()
    return (b[0] and b[1]) or (b[2] and b[3]) or (b[4] and b[5])

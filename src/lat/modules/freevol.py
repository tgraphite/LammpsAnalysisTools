"""lat freevol — probe-accessible free volume from a periodic voxel grid.

Modes (default `ffv`):

* ``ffv``           free volume fraction per frame. A voxel is *occupied* when
                    its centre is within ``r_vdW + r_probe`` of any atom.
* ``psd``           pore size distribution. The distance transform is taken to
                    the *vdW* surface (Gelb-Gubbins convention, so the pore
                    radius is the radius of the largest sphere that fits), the
                    probe-inflated grid decides which voxels are accessible,
                    and every accessible voxel is assigned to its nearest pore
                    centre (periodic Voronoi) so that ``integral dV/dr dr = 1``.
* ``connectivity``  connected components of the free space, closed across the
                    periodic faces, summarised per frame.

Grid analyses need a fixed orthorhombic cell: a Cartesian voxel grid only tiles
an orthorhombic cell periodically, so a triclinic cell or a cell that drifts
over the selected frames is an error rather than a silently wrong number.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from ..core.errors import LatError
from ..core.periodic import (
    grid_axes,
    label_periodic,
    percolating_components,
    periodic_distance_transform,
    periodic_nearest,
    voxel_occupancy,
)
from ..core.table import Field, Table
from ..core.topology import vdw_radii
from .base import Module, Param
from . import register

MODES = ("ffv", "psd", "connectivity")

# A bool occupancy grid is cheap; the distance transform needs float64 arrays,
# so PSD/connectivity coarsen earlier than a plain FFV run.
MAX_VOXELS = {m: (2.5e7 if m in ("psd", "connectivity") else 5e7) for m in MODES}

BOX_WARN_TOL = 1e-3     # relative cell drift that earns a WARN
BOX_ERR_TOL = 2e-2      # ... and that is a hard error
VORONOI_CHUNK = 1_000_000


def _flat_to_cart(idx, dims, lo, spacing):
    """Cartesian coordinates of flat voxel indices on the grid."""
    nz = dims[2]
    iz = idx % nz
    iy = (idx // nz) % dims[1]
    ix = idx // (nz * dims[1])
    return lo + (np.column_stack([ix, iy, iz]) + 0.5) * spacing


def _pore_centers(r_field, free, spacing):
    """One representative voxel per local-maximum plateau of the distance field.

    Returns the flat indices of the pore centres. Local maxima are taken over
    26-connected plateaus (a cubic-grid EDT is full of ridges, so a plain
    equality test would count every voxel of a plateau).
    """
    size = 3
    local = ndimage.maximum_filter(r_field, size=size, mode="nearest")
    plateau = free & (r_field >= local - 1e-12) & (r_field > 0)
    if not plateau.any():
        return np.empty(0, dtype=int)
    labels, ncomp = label_periodic(plateau, connectivity=26)
    lab_p = labels.ravel()
    rv = r_field.ravel()
    plat_idx = np.flatnonzero(plateau.ravel())
    lab_sel = lab_p[plat_idx]
    mx = np.full(ncomp + 1, -np.inf)
    np.maximum.at(mx, lab_sel, rv[plat_idx])
    is_max = rv[plat_idx] >= mx[lab_sel] - 1e-12
    cand = plat_idx[is_max]
    _, first = np.unique(lab_p[cand], return_index=True)
    return cand[first]


@register
class FreeVol(Module):
    name = "freevol"
    summary = "probe-accessible free volume: FFV, PSD, free-space connectivity"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-mode", "-radius", "-grid_dr", "-connectivity"]
    outputs = ["<base>.freevol.out.csv", "<base>.freevol.psd.out.csv",
               "<base>.freevol.connectivity.out.csv"]
    params = [
        Param("radius", flags=("-radius",), config="modules.freevol.radius", type="float"),
        Param("mode", flags=("-mode",), config="modules.freevol.mode", choices=MODES),
        Param("grid_dr", flags=("-grid_dr",), config="modules.freevol.grid_dr", type="float"),
        Param("connectivity", flags=("-connectivity",), config="modules.freevol.connectivity",
              type="int", help="6 = faces, 26 = faces+edges+corners"),
        Param("psd_rmax", flags=("-psd_rmax",), config="modules.freevol.psd_rmax", type="float",
              help="upper PSD radius (A); default: data-driven"),
        Param("components", flags=("-components",), config="modules.freevol.components",
              type="flag", help="also write the per-component connectivity table"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, radius, mode, grid_dr, connectivity, psd_rmax, components,
            frame_range, stride):
        if radius is None or grid_dr is None:
            raise LatError("freevol: -radius and -grid_dr need defaults (check config)")
        if mode not in MODES:
            raise LatError(f"freevol: unknown -mode {mode!r}; choose from {', '.join(MODES)}")
        if int(connectivity) not in (6, 26):
            raise LatError(f"freevol: -connectivity must be 6 or 26, got {connectivity}")

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        r_vdw = vdw_radii(ctx.topo, ctx.cfg, ctx.element_map())
        r_probe = r_vdw + float(radius)

        ref_box = None
        axes = spacing = None
        dims = None
        voxel_volume = 0.0
        bin_edges = None
        psd_hist = None
        psd_centers_total = 0

        ffv_rows = []
        conn_rows = []
        comp_rows = []
        n_seen = 0
        n_warned_box = False

        for fr in ctx.frames(frames):
            if ref_box is None:
                ref_box = fr.box
                axes, spacing, coarsened = grid_axes(
                    ref_box, float(grid_dr), max_points=MAX_VOXELS[mode])
                dims = np.array([len(a) for a in axes], dtype=int)
                voxel_volume = float(np.prod(spacing))
                if coarsened:
                    warnings.append(
                        f"voxel grid coarsened from grid_dr {float(grid_dr):.3g} A to "
                        f"{float(np.mean(spacing)):.3g} A to stay under "
                        f"{MAX_VOXELS[mode]:.0e} voxels ({int(dims.prod())} used)"
                    )
            else:
                rel = float(np.abs(fr.box.H - ref_box.H).max() / ref_box.H.max())
                if rel > BOX_ERR_TOL:
                    raise LatError(
                        "freevol: the cell changes by "
                        f"{rel:.2%} over the selected frames (frame {fr.index}: lengths "
                        f"{np.round(fr.box.lengths, 4)} vs {np.round(ref_box.lengths, 4)}); "
                        f"a voxel grid is tied to one cell -- restrict -frame_range to a "
                        f"constant-volume segment"
                    )
                if rel > BOX_WARN_TOL and not n_warned_box:
                    warnings.append(
                        f"cell drifts by up to {rel:.2%} over the selected frames; the "
                        f"voxel grid is fixed to the first selected frame, so FFV/PSD are "
                        f"approximate (use a constant-volume segment for exact numbers)"
                    )
                    n_warned_box = True

            occ = voxel_occupancy(fr.positions, ref_box, axes, r_probe)
            free = ~occ
            n_seen += 1

            if mode == "ffv":
                ffv_rows.append((fr.index, float(free.mean())))
                continue

            if not free.any():
                warnings.append(f"frame {fr.index}: no probe-accessible voxel")
                if mode == "connectivity":
                    conn_rows.append((fr.index, 0, 0, 0.0, 0.0, 0.0))
                continue

            if mode == "psd":
                occ_vdw = voxel_occupancy(fr.positions, ref_box, axes, r_vdw)
                dt = periodic_distance_transform(occ_vdw, spacing)
                r_field = np.where(free, dt, 0.0)
                rmax = float(r_field.max())
                width = float(np.mean(spacing))
                if bin_edges is None:
                    top = float(psd_rmax) if psd_rmax else max(
                        np.ceil(rmax / width) * width, 2 * width)
                    bin_edges = np.arange(0.0, top + width, width)
                    psd_hist = np.zeros(len(bin_edges) - 1)
                rep = _pore_centers(r_field, free, spacing)
                if len(rep) == 0:
                    warnings.append(f"frame {fr.index}: no pore centre found")
                    continue
                psd_centers_total += len(rep)
                dims_lo = ref_box.lo
                centres = _flat_to_cart(rep, dims, dims_lo, spacing)
                r_centres = r_field.ravel()[rep]
                free_idx = np.flatnonzero(free.ravel())
                counts = np.zeros(len(rep), dtype=float)
                for s in range(0, len(free_idx), VORONOI_CHUNK):
                    chunk = free_idx[s:s + VORONOI_CHUNK]
                    nearest, _ = periodic_nearest(
                        _flat_to_cart(chunk, dims, dims_lo, spacing), centres, ref_box)
                    counts += np.bincount(nearest, minlength=len(rep))
                weights = counts * voxel_volume
                h, _ = np.histogram(r_centres, bins=bin_edges, weights=weights)
                if h.sum() > 0:
                    psd_hist += h
                    overflow = float(weights.sum() - h.sum())
                    if overflow / weights.sum() > 1e-3:
                        warnings.append(
                            f"frame {fr.index}: {overflow / weights.sum():.1%} of the "
                            f"free volume sits above -psd_rmax; raise it for a complete PSD"
                        )
                continue

            labels, ncomp = label_periodic(free, int(connectivity))
            percolating = percolating_components(labels)
            counts = np.bincount(labels.ravel(), minlength=ncomp + 1).astype(float)
            counts[0] = 0.0
            vols = counts * voxel_volume
            v_largest = float(vols.max()) if ncomp else 0.0
            v_total = float(vols.sum())
            v_perc = float(sum(vols[c] for c in percolating))
            conn_rows.append((fr.index, ncomp, len(percolating), v_largest,
                              v_largest / v_total if v_total > 0 else 0.0, v_perc))
            if components:
                order = np.argsort(vols)[::-1]
                for cid in order:
                    if cid == 0 or vols[cid] <= 0:
                        continue
                    comp_rows.append((fr.index, float(vols[cid]),
                                      1.0 if int(cid) in percolating else 0.0))

        if n_seen == 0:
            raise LatError("freevol: no frames selected")

        tables = []
        definitions = {
            "radius": f"{radius} A",
            "grid_dr": (f"{float(np.mean(spacing)):.4g} A (requested {float(grid_dr):g} A)"
                        if ref_box is not None else f"{float(grid_dr):g} A"),
            "voxels": str(int(dims.prod())) if dims is not None else "n/a",
            "cell": "first selected frame (fixed)",
        }
        fr_str = ctx.frame_range_str(frame_range, frames)

        if mode == "ffv":
            data = np.array(ffv_rows, float)
            tables.append(Table(
                title="Free volume fraction",
                columns=[Field("frame", "frame"), Field("ffv", "1")],
                data=data, tag="freevol",
                xlabel="frame",
                ylabel="free volume fraction",
                definitions=definitions,
                scalars={"ffv_mean": float(data[:, 1].mean()),
                         "ffv_std": float(data[:, 1].std())},
                warnings=warnings, meta={"module": "freevol"},
                selection={"atoms": "all"},
                frame_range=fr_str,
            ))

        if mode == "psd":
            if psd_hist is None or psd_hist.sum() <= 0:
                raise LatError("freevol -mode psd: no pore centre was found in any frame")
            width = float(bin_edges[1] - bin_edges[0])
            rmid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            dv = psd_hist / (psd_hist.sum() * width)
            tables.append(Table(
                title="Pore size distribution (probe-grid, periodic)",
                columns=[Field("r", "Angstrom"), Field("dv_dr", "1/Angstrom")],
                data=np.column_stack([rmid, dv]), tag="freevol.psd",
                xlabel="pore radius (Angstrom)",
                ylabel="dV/dr (1/Angstrom)",
                definitions={**definitions,
                             "definition": "r = distance to the nearest vdW surface at a "
                                           "local maximum; Voronoi-normalised over the "
                                           "probe-accessible volume"},
                scalars={"integral": float(dv.sum() * width),
                         "r_peak": float(rmid[int(np.argmax(dv))]),
                         "n_pore_centers": int(psd_centers_total),
                         "frames": int(n_seen)},
                warnings=warnings, meta={"module": "freevol"},
                selection={"atoms": "all"},
                frame_range=fr_str,
            ))

        if mode == "connectivity":
            data = np.array(conn_rows, float)
            tables.append(Table(
                title="Free-space connectivity",
                columns=[Field("frame", "frame"), Field("n_components", "1"),
                         Field("n_percolating", "1"), Field("v_largest", "Angstrom^3"),
                         Field("frac_largest", "1"),
                         Field("v_percolating", "Angstrom^3")],
                data=data, tag="freevol.connectivity",
                xlabel="frame",
                ylabel="component volume",
                definitions={**definitions, "connectivity": int(connectivity),
                             "percolating": "component spans the cell along >=1 axis"},
                scalars={"n_components_mean": float(data[:, 1].mean()),
                         "frac_largest_mean": float(data[:, 4].mean())},
                warnings=warnings, meta={"module": "freevol"},
                selection={"atoms": "all"},
                frame_range=fr_str,
            ))
            if components and comp_rows:
                tables.append(Table(
                    title="Free-space components",
                    columns=[Field("frame", "frame"), Field("volume", "Angstrom^3"),
                             Field("percolating", "1")],
                    data=np.array(comp_rows, float), tag="freevol.connectivity.components",
                    xlabel="frame",
                    ylabel="component volume",
                    definitions={**definitions, "connectivity": int(connectivity)},
                    warnings=warnings, meta={"module": "freevol"},
                    selection={"atoms": "all"},
                    frame_range=fr_str,
                ))

        if not tables:
            raise LatError(f"freevol: mode {mode!r} produced no output")
        return tables

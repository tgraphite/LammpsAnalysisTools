"""lat hbond — hydrogen bonds aligned with gmx hbond conventions.

Geometry: D...A (default) or H...A distance; D-H...A angle > 180 - angle.
Requires explicit bonds (LAMMPS Bonds section) to identify D-H; no inference.
"""

from __future__ import annotations

import numpy as np

from ..core.correlate import linear_fit
from ..core.errors import LatError
from ..core.grid import box_change_warning, fractional_bins
from ..core.neighbor import pairs_within
from ..core.table import Field, Table
from .base import Module, Param
from . import register

_ACCEPTORS = {"N", "O", "F", "S", "Cl", "Br", "I"}


def _adjacency(n_atoms, bonds):
    if bonds is None or len(bonds) == 0:
        raise LatError("hbond requires bond information (LAMMPS Bonds section); not found")
    adj = [[] for _ in range(n_atoms)]
    for a, b in bonds:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    return adj


@register
class HBond(Module):
    name = "hbond"
    summary = "hydrogen bonds: count, ACF, lifetime and geometry distributions"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel1", "-sel2"]
    outputs = ["<base>.hbond.out.csv", "<base>.hbond.acf.out.csv", "<base>.hbond.dist_*.out.csv"]
    params = [
        Param("sel1", flags=("-sel1",), help="donor group (default: all donors)"),
        Param("sel2", flags=("-sel2",), help="acceptor group (default: all acceptors)"),
        Param("distance", flags=("-distance",), config="modules.hbond.distance", type="float"),
        Param("distance_mode", flags=("-distance_mode",), config="modules.hbond.distance_mode",
              choices=("DA", "HA"), help="DA (gmx default) or HA"),
        Param("angle", flags=("-angle",), config="modules.hbond.angle", type="float"),
        Param("nbins", flags=("-nbins",), config="modules.hbond.nbins", type="int"),
        Param("lifetime", flags=("-lifetime",), config="modules.hbond.compute_lifetime",
              type="flag", help="compute ACF and lifetime"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel1, sel2, distance, distance_mode, angle, nbins,
            lifetime, frame_range, stride):
        topo = ctx.topo
        elem = np.array([topo.type_element.get(int(t)) or "X" for t in topo.types])
        adj = _adjacency(topo.n_atoms, topo.bonds)

        if distance_mode not in ("DA", "HA"):
            raise LatError(f"-distance_mode must be DA|HA, got {distance_mode!r}")

        # group masks
        if sel1 is None and sel2 is None:
            g1 = np.ones(topo.n_atoms, dtype=bool)
            g2 = g1.copy()
        elif sel1 is not None and sel2 is None:
            g1 = _mask(ctx, sel1, topo)
            g2 = g1.copy()
        else:
            g1 = _mask(ctx, sel1, topo) if sel1 is not None else np.ones(topo.n_atoms, bool)
            g2 = _mask(ctx, sel2, topo)
        if not (np.array_equal(g1, g2) or not np.any(g1 & g2)):
            raise LatError("hbond: sel1 and sel2 must be identical or completely non-overlapping")

        donors, hydrogens = _donors(elem, adj, g1)
        acceptors = np.flatnonzero(g2 & np.isin(elem, list(_ACCEPTORS)))
        if not len(donors):
            raise LatError("hbond: no donor (heavy atom bonded to H) found in sel1")
        if not len(acceptors):
            raise LatError("hbond: no acceptor atom found in sel2")

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        if len(frames) < 2:
            raise LatError("hbond needs at least 2 frames")
        warnings: list[str] = []
        same = np.array_equal(g1, g2)

        triple_frames: list[set] = []
        counts = []
        centers = []
        bonds_set = _bonded_pairs(topo.bonds)
        for fr in ctx.frames(frames):
            pos = fr.positions
            triples, mid = _frame_hbonds(donors, acceptors, pos, fr.box,
                                         float(distance), distance_mode,
                                         float(angle), same=same,
                                         bonds_set=bonds_set)
            triple_frames.append(triples)
            counts.append(len(triples))
            centers.append(mid)

        if max(counts) == 0:
            warnings.append("no hydrogen bonds found under the given criteria")

        tables = [Table(
            title="Hydrogen bond count",
            columns=[Field("frame", "frame"), Field("n_hbond", "1")],
            data=np.column_stack([np.arange(len(frames)), np.asarray(counts, float)]),
            tag="hbond",
            xlabel="frame",
            ylabel="n(hbond)",
            definitions=_def_dict(distance, distance_mode, angle),
            warnings=warnings,
            meta={"module": "hbond"},
            selection=_sel_dict(sel1, sel2, distance, distance_mode, angle),
            frame_range=ctx.frame_range_str(frame_range, frames),
            scalars={"mean_n_hbond": float(np.mean(counts))},
        )]

        if lifetime:
            acf_table, tscal, acf_warn = _acf(triple_frames)
            warnings.extend(acf_warn)
            acf_table.warnings = list(warnings)
            tables.append(acf_table)
            tables[0].scalars.update(tscal)

        tables.extend(_geom_tables(centers, frames, int(nbins), warnings, sel1, sel2,
                                   distance, distance_mode, angle, ctx, frame_range))
        return tables


def _mask(ctx, literal, topo):
    sel = ctx.selection(literal)
    m = np.zeros(topo.n_atoms, dtype=bool)
    m[sel.indices] = True
    return m


def _donors(elem, adj, group):
    donors, hydrogens = {}, {}
    for i in range(len(elem)):
        if not group[i] or elem[i] == "H":
            continue
        hs = [j for j in adj[i] if elem[j] == "H" and group[j]]
        if hs:
            donors[i] = hs
    return donors, hydrogens


def _bonded_pairs(bonds):
    s = set()
    if bonds is None:
        return s
    for a, b in bonds:
        s.add((int(a), int(b)))
        s.add((int(b), int(a)))
    return s


def _frame_hbonds(donors, acceptors, pos, box, distance, mode, angle, same, bonds_set):
    """gmx convention: D...A (or H...A) <= distance and angle H-D...A <= angle.

    The angle is evaluated at the DONOR D between D->H and D->A, as in gmx
    hbond (default 30 deg). All vectors use the minimum image.
    """
    triples = set()
    mids = []
    cos_max = np.cos(np.deg2rad(angle))
    acc_pos = pos[acceptors]
    if not len(acceptors):
        return triples, mids

    def _emit(d, h, a):
        vdh = box.minimum_image(pos[h] - pos[d])
        ndh = float(np.linalg.norm(vdh))
        vda = box.minimum_image(pos[a] - pos[d])
        nda = float(np.linalg.norm(vda))
        if ndh == 0 or nda == 0:
            return
        if float(np.dot(vdh, vda)) / (ndh * nda) < cos_max:
            return
        triples.add((d, h, a))
        # midpoint of D and its periodic image nearest to A (then wrapped back
        # into the cell so a boundary-straddling bond is not dropped by the
        # histogram). The image of D nearest A is A + MI(D - A).
        mid = pos[a] + 0.5 * box.minimum_image(pos[d] - pos[a])
        mids.append(box.wrap(mid))

    if mode == "DA":
        # candidate D...A pairs within the cutoff, via periodic neighbor search
        d_ids = np.fromiter(donors.keys(), dtype=int)
        if not len(d_ids):
            return triples, mids
        ia, ja, _ = pairs_within(pos[d_ids], acc_pos, box, distance, same=False)
        for k in range(len(ia)):
            d = int(d_ids[ia[k]])
            a = int(acceptors[ja[k]])
            if a == d or (d, a) in bonds_set:
                continue
            for h in donors[d]:
                if h == a:
                    continue
                _emit(d, h, a)
    else:
        # candidate H...A pairs within the cutoff
        h2d = {h: d for d, hs in donors.items() for h in hs}
        h_ids = np.fromiter(h2d.keys(), dtype=int)
        if not len(h_ids):
            return triples, mids
        ia, ja, _ = pairs_within(pos[h_ids], acc_pos, box, distance, same=False)
        for k in range(len(ia)):
            h = int(h_ids[ia[k]])
            a = int(acceptors[ja[k]])
            d = int(h2d[h])
            if a == d or a == h or (d, a) in bonds_set:
                continue
            _emit(d, h, a)
    return triples, mids


def _acf(triple_frames):
    T = len(triple_frames)
    universe = sorted(set().union(*triple_frames)) if triple_frames else []
    warns = []
    if not universe:
        return (Table(title="Hbond ACF", columns=[Field("lag", "frame"),
                     Field("C_inter", "1"), Field("C_cont", "1")],
                     data=np.zeros((0, 3)), tag="hbond.acf",
                     warnings=["undefined"], meta={"module": "hbond"}), {}, [])

    idx = {t: i for i, t in enumerate(universe)}
    present = np.zeros((len(universe), T), dtype=bool)
    for k, s in enumerate(triple_frames):
        for t in s:
            present[idx[t], k] = True
    total = present.sum()
    # unbiased over the T-lag valid origins: C(tau) = <h(0)h(tau)> / <h>
    C_inter = np.zeros(T)
    for lag in range(T):
        num = (present[:, : T - lag] & present[:, lag:]).sum()
        C_inter[lag] = num * T / ((T - lag) * total) if T > lag else 0.0
    C_cont = np.zeros(T)
    for lag in range(T):
        if lag == 0:
            C_cont[lag] = 1.0
            continue
        run = present[:, : T - lag].copy()
        for k in range(1, lag + 1):
            run &= present[:, k: T - lag + k]
        C_cont[lag] = run.sum() * T / ((T - lag) * total)

    # lifetime: integral and exponential fit
    tau_integral = float(C_inter[1:].sum()) if T > 1 else 0.0
    fit_start = 1
    y = np.clip(C_inter[fit_start:], 1e-12, None)
    ok = y > 0
    if ok.sum() >= 3:
        slope, intercept, r2 = linear_fit(np.arange(fit_start, T)[ok], np.log(y[ok]))
        tau_fit = float(-1.0 / slope) if slope < 0 else float("nan")
        r2_fit = r2
        if not (0 < tau_fit <= T):
            tau_fit = float("nan")
            warns.append("ACF does not decay within the trajectory; lifetime fit unreliable")
    else:
        tau_fit, r2_fit = float("nan"), float("nan")
        warns.append("too few points for a lifetime fit; only integral reported")

    data = np.column_stack([np.arange(T), C_inter, C_cont])
    table = Table(
        title="Hydrogen bond ACF",
        columns=[Field("lag", "frame"), Field("C_inter", "1"), Field("C_cont", "1")],
        data=data,
        tag="hbond.acf",
        xlabel="lag (frame)",
        ylabel="C(t)",
        scalars={"tau_integral": tau_integral, "tau_fit": tau_fit, "r2": r2_fit,
                 "acf_unit": "frame"},
        warnings=warns,
        meta={"module": "hbond"},
    )
    return table, {"tau_integral": tau_integral, "tau_fit": tau_fit}, warns


def _geom_tables(centers, frames, nbins, warnings, sel1, sel2, distance, mode, angle,
                 ctx, frame_range):
    # gmx-density-style binning: fractional (box-relative) bins per frame, axis
    # scaled back to the mean box length. The donor-acceptor midpoints are in
    # absolute coordinates but each frame's box may differ under NPT.
    edges, centers_s = fractional_bins(nbins)
    names = [(0, 2, "xz"), (1, 2, "yz"), (0, 1, "xy")]
    acc = {name: np.zeros((nbins, nbins)) for _, _, name in names}
    arrays = [np.asarray(c, float).reshape(-1, 3) for c in centers]
    lengths = np.zeros((len(frames), 3))
    origin_sum = np.zeros(3)
    n = 0
    for fr, arr in zip(ctx.frames(frames), arrays):
        L = fr.box.lengths
        lo = fr.box.origin
        for (a, b, name) in names:
            if len(arr) == 0:
                continue
            sa = (arr[:, a] - lo[a]) / L[a]
            sb = (arr[:, b] - lo[b]) / L[b]
            sa -= np.floor(sa)
            sb -= np.floor(sb)
            h, _, _ = np.histogram2d(sa, sb, bins=[edges, edges])
            acc[name] += h
        lengths[n] = L
        origin_sum += lo
        n += 1
    n = max(n, 1)
    mean_L = lengths.mean(axis=0)
    mean_lo = origin_sum / n
    w = box_change_warning(lengths, [0, 1, 2])
    if w:
        warnings.append(w)

    tables = []
    for (a, b, name) in names:
        ca = mean_lo[a] + centers_s * mean_L[a]
        cb = mean_lo[b] + centers_s * mean_L[b]
        aa, bb = np.meshgrid(ca, cb, indexing="ij")
        tables.append(Table(
            title="Hbond geometry",
            columns=[Field("xyz"[a], "Angstrom"), Field("xyz"[b], "Angstrom"),
                     Field("count", "1/frame")],
            data=np.column_stack([aa.ravel(), bb.ravel(), (acc[name] / n).ravel()]),
            tag=f"hbond.dist_{name}",
            kind="grid",
            definitions=_def_dict(distance, mode, angle),
            warnings=warnings,
            meta={"module": "hbond"},
            selection=_sel_dict(sel1, sel2, distance, mode, angle),
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))
    return tables


def _sel_dict(sel1, sel2, distance, mode, angle):
    d = {}
    if sel1:
        d["sel1"] = sel1
    if sel2:
        d["sel2"] = sel2
    return d


def _def_dict(distance, mode, angle):
    return {"distance": f"{distance} A", "distance_mode": mode, "angle": f"{angle} deg"}

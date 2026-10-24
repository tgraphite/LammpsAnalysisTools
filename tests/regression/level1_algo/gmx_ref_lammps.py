#!/usr/bin/env python3
"""Level-1 regression for **LAMMPS input**: lat vs native gmx ground truth.

Complements `compare_gmx.py` (GROMACS-native input) and covers the half of the
IO path that had no native reference at all: a LAMMPS `data` + `lammpstrj`.

It converts the system with `../convert/lammps_to_gmx.py`, then runs, on the
same frames and the same converted coordinates:

    module   gmx command                        lat module
    ------   --------------------------------   ---------
    gyrate   gmx gyrate        (System, MolN)   gyrate
    rms      gmx rms    -nomw -fit rot+trans    rms
    rmsf     gmx rmsf   -fit yes                rms  -rmsf
    sasa     gmx sasa   -probe -ndots           sasa
    freevol  (no gmx analogue; integral + connectivity reported for the record)
    rdf      gmx rdf    -bin -rmax              rdf  -intramol
    hbond    gmx hbond  -num                    hbond

Notes that cost real debugging time and are therefore encoded here:
* gmx applies `rmpbc` (make molecules whole) itself when given a .tpr, so the
  comparison is make-whole against make-whole.  Feeding gmx a .gro as `-s`
  silently disables that.
* `gmx rms` defaults to **mass-weighted** superposition (`-mw yes`); lat's
  Kabsch is unweighted, so `-nomw` is mandatory for a like-for-like number
  (otherwise gmx is ~2% lower for a C/H/O polymer).
* SASA radii: gmx uses r = 0.5*sigma*2^(1/6), which the converter matches to
  lat's vdW table.
* gmx rdf is selected by index-group name (a `name C*` wildcard selection is
  rejected by gmx 2022's parser in this context).

Usage
-----
  gmx_ref_lammps.py --data system.in.data --traj run.lammpstrj \
      --out gmx_ref --rdf-pairs C:O,C:C,O:O
  gmx_ref_lammps.py ... --frames 0:101:1 --quick
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "tests" / "regression" / "convert"))

import lammps_to_gmx as conv  # noqa: E402
from lat.config import load as load_config  # noqa: E402
from lat.core.context import Context  # noqa: E402
from lat.core.io import open_system  # noqa: E402
from lat.modules import get  # noqa: E402


def read_xvg(path):
    d = np.loadtxt(path, comments=["@", "#"])
    return d if d.ndim == 2 else d.reshape(1, -1)


def gmx(cmd, stdin=None):
    return conv.gmx(cmd, stdin=stdin)


class Case:
    """One lat-vs-gmx comparison on aligned 1-D series."""

    def __init__(self, name, unit):
        self.name, self.unit = name, unit
        self.n = self.mean_lat = self.mean_gmx = None
        self.mean_abs = self.max_abs = self.rel = None
        self.note = ""

    def series(self, lat, gmx):
        self.n = min(len(lat), len(gmx))
        a, b = np.asarray(lat[:self.n], float), np.asarray(gmx[:self.n], float)
        d = a - b
        self.mean_lat, self.mean_gmx = float(a.mean()), float(b.mean())
        self.mean_abs, self.max_abs = float(np.abs(d).mean()), float(np.abs(d).max())
        den = abs(self.mean_gmx) if self.mean_gmx else 1.0
        self.rel = 100.0 * self.mean_abs / den
        return self


def compare_rdf(x_lat, g_lat, cn_lat, g_gmx, cn_gmx, r_min):
    # lat reports r in Angstrom, gmx in nm: convert the gmx grid to Angstrom
    # before interpolating (interpolating A against an nm grid silently clamps
    # everything past 1 nm to a constant).
    gi = np.interp(x_lat, g_gmx[:, 0] * 10.0, g_gmx[:, 1])
    ci = np.interp(x_lat, cn_gmx[:, 0] * 10.0, cn_gmx[:, 1])
    m = x_lat > r_min
    d = (gi - g_lat)[m]
    c = Case("rdf", "g(r)")
    c.n = int(m.sum())
    c.mean_lat, c.mean_gmx = float(g_lat[m].mean()), float(gi[m].mean())
    c.mean_abs = float(np.abs(d).mean())
    c.max_abs = float(np.abs(d).max())
    c.rel = 100.0 * abs(cn_lat[-1] - ci[-1]) / abs(ci[-1]) if ci[-1] else np.nan
    c.note = f"cn lat={cn_lat[-1]:.3f} gmx={ci[-1]:.3f} (cn_rel column)"
    return c


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--out", default="gmx_ref", help="working/output directory")
    ap.add_argument("--frames", default=None, help="start:stop (exclusive) of the dump")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--rdf-pairs", default="C:O,C:C,O:O")
    ap.add_argument("--rdf-bin", type=float, default=0.002, help="nm (0.02 A, as in compare_gmx.py)")
    ap.add_argument("--rdf-rmax", type=float, default=1.0, help="nm")
    ap.add_argument("--rdf-min", type=float, default=0.12, help="nm, ignore below")
    ap.add_argument("--json", help="also write the metrics here")
    a = ap.parse_args(argv)

    if not conv.shutil.which(conv.GMX):
        print(f"SKIP: gmx ({conv.GMX}) not found", file=sys.stderr)
        return 2

    work = Path(a.out)
    win = None
    if a.frames:
        p = a.frames.split(":")
        win = (int(p[0]), int(p[1]))
    conv.convert(a.data, a.traj, work, stride=a.stride, frames=win)
    tpr, ref, xtc = work / "sys.tpr", work / "ref.tpr", work / "whole.xtc"

    cfg = load_config()
    topo, source = open_system(a.data, a.traj, cfg=cfg)
    ctx = Context(topo, source, cfg=cfg)
    nf = ctx.n_frames
    fr = a.frames if a.frames else f"0:{nf}:{a.stride}"
    probe = float(cfg.get("modules.sasa.radius"))
    ndots = int(cfg.get("modules.sasa.sasa_points"))
    hb = dict(distance=float(cfg.get("modules.hbond.distance")),
              angle=float(cfg.get("modules.hbond.angle")),
              mode=cfg.get("modules.hbond.distance_mode"))

    results, order = {}, []

    def add(case):
        results[case.name] = case.__dict__
        order.append(case.name)
        print(f"  {case.name:22s} n={case.n:<5d} lat={case.mean_lat:12.4f} "
              f"gmx={case.mean_gmx:12.4f} mean|d|={case.mean_abs:10.4f} "
              f"max|d|={case.max_abs:10.4f} ({case.rel:6.2f}%) {case.note}")

    # ---- gyrate: whole system + per molecule -----------------------------
    print("[gyrate]")
    o = work / "rg.xvg"
    gmx(["gyrate", "-f", str(xtc), "-s", str(tpr), "-n", str(work / "sys.ndx"),
         "-o", str(o)], stdin="System\n")
    lat = get("gyrate")().run(ctx, sel="all", frame_range=fr, stride=None)
    add(Case("gyrate_all", "A").series(lat[0].data[:, 1], read_xvg(o)[:, 1] * 10.0))
    n_mol = len(np.unique(topo.mols))
    lat_mol = get("gyrate")().run(ctx, sel=f"mol 1", frame_range=fr, stride=None)
    for m in range(1, n_mol + 1):
        o = work / f"rg_m{m}.xvg"
        gmx(["gyrate", "-f", str(xtc), "-s", str(tpr), "-n", str(work / "sys.ndx"),
             "-o", str(o)], stdin=f"Mol{m}\n")
        if m > 1:
            lat_mol = get("gyrate")().run(ctx, sel=f"mol {m}", frame_range=fr, stride=None)
        add(Case(f"gyrate_mol{m}", "A").series(lat_mol[0].data[:, 1],
                                               read_xvg(o)[:, 1] * 10.0))

    # ---- rms / rmsf ------------------------------------------------------
    print("[rms / rmsf]")
    o = work / "rmsd.xvg"
    gmx(["rms", "-s", str(ref), "-f", str(xtc), "-n", str(work / "sys.ndx"),
         "-fit", "rot+trans", "-nomw", "-o", str(o)], stdin="System\nSystem\n")
    rms_t = get("rms")().run(ctx, sel="all", fit=None, ref="first", rmsf=True,
                             frame_range=fr, stride=None)
    add(Case("rmsd", "A").series(rms_t[0].data[:, 1], read_xvg(o)[:, 1] * 10.0))

    o = work / "rmsf.xvg"
    gmx(["rmsf", "-s", str(ref), "-f", str(xtc), "-n", str(work / "sys.ndx"),
         "-fit", "yes", "-o", str(o)], stdin="System\n")
    add(Case("rmsf", "A").series(rms_t[1].data[:, 1], read_xvg(o)[:, 1] * 10.0))

    # ---- sasa ------------------------------------------------------------
    print("[sasa]")
    o = work / "sasa.xvg"
    gmx(["sasa", "-s", str(tpr), "-f", str(xtc), "-surface", "System",
         "-probe", str(probe / 10.0), "-ndots", str(ndots), "-o", str(o)])
    lat = get("sasa")().run(ctx, radius=probe, sasa_points=ndots, sel=None,
                            frame_range=fr, stride=None)
    add(Case("sasa", "A^2").series(lat[0].data[:, 1], read_xvg(o)[:, 1] * 100.0))

    # ---- freevol: no native reference, so report internally-checkable facts --
    print("[freevol]")
    dr = float(cfg.get("modules.freevol.grid_dr"))
    psd = get("freevol")().run(ctx, radius=probe, mode="psd", grid_dr=dr,
                               connectivity=26, psd_rmax=None, components=False,
                               frame_range="0:1:1", stride=None)[0]
    conn = get("freevol")().run(ctx, radius=probe, mode="connectivity", grid_dr=dr,
                                connectivity=26, psd_rmax=None, components=False,
                                frame_range=fr, stride=None)[0]
    sanity = {
        "psd_integral": float(psd.scalars["integral"]),
        "psd_r_peak": float(psd.scalars["r_peak"]),
        "psd_n_pore_centers": int(psd.scalars["n_pore_centers"]),
        "psd_grid_dr": psd.definitions["grid_dr"],
        "connectivity_n_components_mean": float(conn.scalars["n_components_mean"]),
        "connectivity_frac_largest_mean": float(conn.scalars["frac_largest_mean"]),
        "connectivity": conn.definitions["connectivity"],
    }
    results["_freevol_sanity"] = sanity
    for k, v in sanity.items():
        print(f"  {k:34s} {v}")

    # ---- rdf (intramolecular included, like gmx) -------------------------
    print("[rdf]")
    nbins = int(round(a.rdf_rmax / a.rdf_bin))
    for pair in a.rdf_pairs.split(","):
        e1, e2 = pair.split(":")
        g, cn = work / "g.xvg", work / "cn.xvg"
        gmx(["rdf", "-f", str(xtc), "-s", str(tpr), "-n", str(work / "sys.ndx"),
             "-o", str(g), "-cn", str(cn), "-bin", str(a.rdf_bin),
             "-rmax", str(a.rdf_rmax), "-ref", f"El_{e1}", "-sel", f"El_{e2}"])
        t = get("rdf")().run(ctx, sel1=f"element {e1}", sel2=f"element {e2}",
                             nbins=nbins, cutoff=a.rdf_rmax * 10.0,
                             frame_range=fr, stride=None, intramol=True)
        c = compare_rdf(t[0].data[:, 0], t[0].data[:, 1], t[1].data[:, 1],
                        read_xvg(g), read_xvg(cn), a.rdf_min * 10.0)
        c.name = f"rdf_{e1}-{e2}"
        add(c)

    # ---- hbond -----------------------------------------------------------
    print("[hbond]")
    o = work / "hbnum.xvg"
    ok, log = conv.run_gmx(["hbond", "-f", str(xtc), "-s", str(tpr),
                            "-n", str(work / "sys.ndx"), "-num", str(o)],
                           stdin="System\nSystem\n")
    if not ok:
        c = Case("hbond", "count")
        c.note = "gmx hbond failed: " + log.strip().splitlines()[-1][:80]
        add(c)
    else:
        t = get("hbond")().run(ctx, sel1=None, sel2=None, distance=hb["distance"],
                               distance_mode=hb["mode"], angle=hb["angle"],
                               nbins=50, lifetime=False, frame_range=fr, stride=None)
        h = read_xvg(o)
        add(Case("hbond", "count").series(t[0].data[:, 1], h[:, 1]))
        results[order[-1]]["note"] = (
            f"lat donors: {t[0].scalars.get('mean_n_hbond', float('nan')):.3f} mean; "
            f"gmx within-0.35nm pairs mean {h[:, 2].mean():.1f}")

    print("\n| quantity | n | lat | gmx | mean\\|d\\| | max\\|d\\| | rel | note |")
    print("|---|---|---|---|---|---|---|---|")
    for k in order:
        r = results[k]
        if r["rel"] is None:
            print(f"| {k} | - | - | - | - | - | - | {r['note']} |")
            continue
        note = r.get("note", "")
        if k.startswith("rdf_"):
            note = "rel = coordination-number rel. diff. " + note
        print(f"| {k} | {r['n']} | {r['mean_lat']:.4f} | {r['mean_gmx']:.4f} | "
              f"{r['mean_abs']:.4f} | {r['max_abs']:.4f} | {r['rel']:.2f}% | {note} |")

    if a.json:
        Path(a.json).write_text(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())

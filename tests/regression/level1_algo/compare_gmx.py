#!/usr/bin/env python3
"""Level-1 regression: lat vs GROMACS on real electrolyte/polymer systems.

Loads the gmx topology once through lat's API, runs RDF for each (ref, sel)
pair, and compares against a freshly generated `gmx rdf` ground truth on a fine
grid. Reports max/RMS/max-relative error, peak position/height and shell CN.

Usage:
  compare_gmx.py --system sys1 --base <11-systems> [--pairs LI:FSI,LI:FEC]
                 [--rmax 3.0] [--bin 0.001] [--frames 90:101:1] [--work DIR]
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from lat.config import load as load_config  # noqa: E402
from lat.core.context import Context  # noqa: E402
from lat.core.io import open_system  # noqa: E402
from lat.modules import get  # noqa: E402

GMX = "/usr/local/gromacs/bin/gmx"


def read_xvg(path):
    data = np.loadtxt(path, comments=["@", "#"])
    return data if data.ndim == 2 else data.reshape(1, -1)


def run_gmx(base_xtc, base_tpr, ref, sels, bin_nm, rmax, b, e, outdir):
    out = outdir / "gmx_g.xvg"
    cn = outdir / "gmx_cn.xvg"
    env = {"GMX_MAXBACKUP": "-1", "PATH": "/usr/local/gromacs/bin:/usr/bin:/bin"}
    cmd = [GMX, "rdf", "-f", str(base_xtc), "-s", str(base_tpr),
           "-o", str(out), "-cn", str(cn), "-bin", str(bin_nm),
           "-sel", *sels, "-ref", ref, "-b", str(b), "-e", str(e), "-rmax", str(rmax)]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"gmx rdf failed:\n{proc.stderr[-2000:]}")
    return out, cn


def compare(x, gmx_g, lat_g, gmx_cn=None, lat_cn=None, label=""):
    gi = np.interp(x, gmx_g[:, 0], gmx_g[:, 1])
    m = (x > 0.12) & (x < x.max() - 0.02)
    d = (gi - lat_g)[m]
    res = {
        "max_abs": float(np.abs(d).max()),
        "rms": float(np.sqrt((d ** 2).mean())),
        "mean": float(d.mean()),
        "peak_gmx_r": float(gmx_g[np.argmax(gmx_g[:, 1]), 0]),
        "peak_lat_r": float(x[np.argmax(lat_g)]),
        "peak_gmx": float(gmx_g[:, 1].max()),
        "peak_lat": float(lat_g.max()),
    }
    if gmx_cn is not None and lat_cn is not None:
        ci = np.interp(x, gmx_cn[:, 0], gmx_cn[:, 1])
        dc = (ci - lat_cn)[m]
        res["cn_max_abs"] = float(np.abs(dc).max())
        res["cn_rel"] = float(np.abs(dc).max() / max(gmx_cn[:, 1].max(), 1e-12))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--system", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--pairs", default="LI:LI,LI:FSI,LI:FEC,LI:FEMC,LI:TTE,LI:MBA,LI:PMBA")
    ap.add_argument("--rmax", type=float, default=3.0)
    ap.add_argument("--bin", type=float, default=0.001)
    ap.add_argument("--frames", default="90:101:1")
    ap.add_argument("--work", default="/tmp/opencode/gt_run")
    args = ap.parse_args()

    base = Path(args.base) / f"{args.system}_npt"
    tpr, xtc = Path(str(base) + ".tpr"), Path(str(base) + ".xtc")
    work = Path(args.work) / args.system
    work.mkdir(parents=True, exist_ok=True)

    pairs = [tuple(p.split(":")) for p in args.pairs.split(",")]
    pairs = [(r, s) for r, s in pairs if s in _available_resnames(tpr)]
    ref = pairs[0][0]
    sels = sorted({s for _, s in pairs})

    start, stop, stride = (int(v) for v in args.frames.split(":"))
    b = 100.0 * start
    e = 100.0 * (stop - 1)
    gmx_g_path, gmx_cn_path = run_gmx(xtc, tpr, ref, sels, args.bin, args.rmax, b, e, work)
    gmx_g = read_xvg(gmx_g_path)
    gmx_cn = read_xvg(gmx_cn_path)
    col = {name: i + 1 for i, name in enumerate(sels)}
    cncol = {name: i + 1 for i, name in enumerate(sels)}

    cfg = load_config()
    topo, source = open_system(str(tpr), str(xtc), cfg=cfg)
    ctx = Context(topo, source, cfg=cfg)
    nbins = int(round(args.rmax / args.bin))  # rmax[nm]/bin[nm], lat uses Angstrom internally

    print(f"# {args.system}: atoms={topo.n_atoms} frames={ctx.n_frames} nbins={nbins} bin={args.bin} nm")
    print(f"{'pair':12s} {'max|dg|':>8s} {'rms':>8s} {'mean':>9s} "
          f"{'peak_gmx':>9s} {'peak_lat':>9s} {'cn_max|d|':>10s} {'cn_rel':>8s}")
    worst = 0.0
    for r, s in pairs:
        mod = get("rdf")()
        tables = mod.run(ctx, sel1=f"resname {r}", sel2=f"resname {s}",
                         nbins=nbins, cutoff=args.rmax * 10.0,
                         frame_range=args.frames, stride=None, intramol=True)
        g_tab, cn_tab = tables
        x = g_tab.data[:, 0] / 10.0
        stats = compare(x, gmx_g[:, [0, col[s]]], g_tab.data[:, 1],
                        gmx_cn=gmx_cn[:, [0, cncol[s]]], lat_cn=cn_tab.data[:, 1],
                        label=f"{r}-{s}")
        worst = max(worst, stats["max_abs"])
        print(f"{r+'-'+s:12s} {stats['max_abs']:8.4g} {stats['rms']:8.4g} "
              f"{stats['mean']:+9.2g} {stats['peak_gmx']:9.3f} {stats['peak_lat']:9.3f} "
              f"{stats.get('cn_max_abs', float('nan')):10.4g} "
              f"{stats.get('cn_rel', float('nan')):8.3g}")
    print(f"# worst max|dg| = {worst:.4g}")


def _available_resnames(tpr):
    import MDAnalysis as mda
    import warnings
    warnings.filterwarnings("ignore")
    u = mda.Universe(str(tpr))
    return set(map(str, u.atoms.resnames))


if __name__ == "__main__":
    main()

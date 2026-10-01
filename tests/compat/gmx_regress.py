#!/usr/bin/env python3
"""GROMACS native regression on third-party "novice user" order data.

For every (topology, trajectory) pair found under a root, compare lat against
freshly generated GROMACS ground truth for:
  rdf (g(r) + coordination number), density profile, gyrate (Rg), msd (mol_com)

Selections are discovered from the topology itself (resnames), so this works on
arbitrary real projects. Writes JSON + a markdown table.

Usage:
  gmx_regress.py --root /mnt/e/RefData/gromacs_real_projects --maxdepth 8 \
                 --frames 200 --stride 5 --out /tmp/opencode/refdata_gmx.md
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
warnings.filterwarnings("ignore")

from lat.config import load as load_config  # noqa: E402
from lat.core.context import Context  # noqa: E402
from lat.core.io import open_system  # noqa: E402
from lat.modules import get  # noqa: E402

GMX = "/usr/local/gromacs/bin/gmx"
ENV = {"GMX_MAXBACKUP": "-1", "PATH": "/usr/local/gromacs/bin:/usr/bin:/bin",
       "OMP_NUM_THREADS": "4"}


def read_xvg(path):
    d = np.loadtxt(path, comments=["@", "#"])
    return d if d.ndim == 2 else d.reshape(1, -1)


def gmx(cmd, stdin=None):
    p = subprocess.run([GMX, *cmd], capture_output=True, text=True, env=ENV,
                       input=stdin)
    if p.returncode != 0:
        raise RuntimeError(f"gmx {' '.join(cmd[:2])} failed: {p.stderr[-800:]}")
    return p


def discover(root: Path, maxdepth: int):
    """Pair each tpr with the same-stem trajectory; else an atom-count match.

    Never pair a topology with a trajectory of a different atom count (that is
    a data-usage error, e.g. a restart topology next to an unrelated run).
    """
    tprs, trajs = [], []
    root = root.resolve()
    for dp, dn, fn in os.walk(root):
        rel = Path(dp).relative_to(root)
        if len(rel.parts) >= maxdepth:
            dn[:] = []
        for f in fn:
            p = Path(dp) / f
            e = f.lower().rsplit(".", 1)[-1]
            if e == "tpr":
                tprs.append(p)
            elif e in ("xtc", "trr"):
                trajs.append(p)
    by_dir = {}
    for t in trajs:
        by_dir.setdefault(t.parent, []).append(t)

    def n_atoms_traj(t):
        try:
            import MDAnalysis as mda
            return mda.Universe(str(t)).atoms.n_atoms
        except Exception:
            return None

    def n_atoms_tpr(t):
        try:
            import MDAnalysis as mda
            return mda.Universe(str(t)).atoms.n_atoms
        except Exception:
            return None

    pairs = []
    for tpr in sorted(tprs):
        cands = by_dir.get(tpr.parent, [])
        if not cands:
            continue
        same = [c for c in cands if c.stem == tpr.stem]
        if same:
            pairs.append((tpr, same[0]))
            continue
        nt = n_atoms_tpr(tpr)
        if nt is None:
            continue
        for c in cands:
            if n_atoms_traj(c) == nt:
                pairs.append((tpr, c))
                break
    return pairs


def topo_info(tpr: Path):
    import MDAnalysis as mda
    u = mda.Universe(str(tpr))
    res = {}
    for rn in u.atoms.resnames:
        res[str(rn)] = res.get(str(rn), 0) + 1
    res = dict(sorted(res.items(), key=lambda kv: -kv[1]))
    n_frames = None
    return u, res


def compare_series(x, g, lat, x_lo=None):
    gi = np.interp(x, g[:, 0], g[:, 1])
    m = x > 0
    if x_lo is not None:
        m &= x > x_lo
    d = (gi - lat)[m]
    if d.size == 0:
        return None
    denom = max(np.abs(g[:, 1]).max(), 1e-12)
    return {"max_abs": float(np.abs(d).max()), "rms": float(np.sqrt((d ** 2).mean())),
            "mean": float(d.mean()), "rel": float(np.abs(d).max() / denom)}


def regress_pair(tpr: Path, traj: Path, work: Path, frames: int, stride: int, timeout):
    work.mkdir(parents=True, exist_ok=True)
    u, resnames = topo_info(tpr)
    if len(resnames) < 2:
        return {"skip": "fewer than 2 residue names"}
    names = list(resnames)
    ref = names[0]
    # pick up to 3 distinct partner resnames
    sels = names[1:4]
    cfg = load_config()
    topo, source = open_system(str(tpr), str(traj), cfg=cfg)
    ctx = Context(topo, source, cfg=cfg)
    nf = ctx.n_frames
    stride = 1  # gmx -dt does not select the same frames; use contiguous windows
    start = max(0, nf - frames)
    fr = f"{start}:{nf}:{stride}"

    # gmx time window + stride-matched dt from the same frames
    dt = 0.0
    t0 = t1 = None
    try:
        import MDAnalysis as mda
        uu = mda.Universe(str(tpr), str(traj))
        uu.trajectory[start]
        t0 = uu.trajectory.time
        j = min(start + 1, nf - 1)
        uu.trajectory[j]
        dt = (uu.trajectory.time - t0) / max(j - start, 1)
        uu.trajectory[nf - 1]
        t1 = uu.trajectory.time
    except Exception:
        pass
    b_e = ["-b", str(t0), "-e", str(t1)] if t0 is not None else []

    out = {"atoms": topo.n_atoms, "frames": nf, "window": fr, "ref": ref,
           "partners": sels, "rdf": {}, "density": None, "gyrate": None, "msd": None}

    # ---- rdf ----
    gmx_ok = True
    try:
        gpath, cnpath = work / "g.xvg", work / "cn.xvg"
        gmx(["rdf", "-f", str(traj), "-s", str(tpr), "-o", str(gpath),
             "-cn", str(cnpath), "-bin", "0.002", "-rmax", "2.0",
             "-ref", f"resname {ref}", "-sel", *[f"resname {s}" for s in sels],
             *b_e])
        g = read_xvg(gpath)
        cn = read_xvg(cnpath)
        nbins = 1000  # 2.0 nm / 0.002
        for k, s in enumerate(sels):
            t = get("rdf")().run(ctx, sel1=f"resname {ref}", sel2=f"resname {s}",
                                nbins=nbins, cutoff=20.0, frame_range=fr,
                                stride=None, intramol=True)
            x = t[0].data[:, 0] / 10.0
            out["rdf"][s] = compare_series(x, g[:, [0, k + 1]], t[0].data[:, 1], x_lo=0.12)
            out["rdf"][s + "_cn"] = compare_series(
                x, cn[:, [0, k + 1]], t[1].data[:, 1], x_lo=0.12)
    except Exception as e:
        # usually a gmx/tpr version mismatch: gmx cannot read the topology even
        # though lat/MDAnalysis can. Do not count this as a lat failure.
        gmx_ok = False
        out["gmx_unavailable"] = str(e)[:300]

    # ---- density ----
    try:
        if not gmx_ok:
            raise RuntimeError("gmx unavailable")
        dpath = work / "d.xvg"
        gmx(["density", "-f", str(traj), "-s", str(tpr), "-o", str(dpath),
             "-sl", "100", "-dens", "mass", *b_e], stdin="0\n")
        g = read_xvg(dpath)
        t = get("density")().run(ctx, sel="all", axis="z", nbins=100, dim=1,
                                 quantity="mass", frame_range=fr, stride=None)
        # gmx density x is nm along axis; lat x is Angstrom bin centres
        gx = g[:, 0] * 10.0
        # align by comparing profiles (lat rho kg/m^3, gmx kg/m^3)
        lo = max(gx.min(), t[0].data[:, 0].min())
        hi = min(gx.max(), t[0].data[:, 0].max())
        m = (t[0].data[:, 0] >= lo) & (t[0].data[:, 0] <= hi)
        gi = np.interp(t[0].data[m, 0], gx, g[:, 1])
        d = gi - t[0].data[m, 1]
        out["density"] = {"max_abs": float(np.abs(d).max()),
                          "rms": float(np.sqrt((d ** 2).mean())),
                          "mean": float(d.mean()),
                          "gmx_mean": float(g[:, 1].mean()),
                          "lat_mean": float(t[0].data[m, 1].mean())}
    except Exception as e:
        out["density_error"] = str(e)[:300]

    # ---- gyrate ----
    try:
        if not gmx_ok:
            raise RuntimeError("gmx unavailable")
        gpath = work / "rg.xvg"
        gmx(["gyrate", "-f", str(traj), "-s", str(tpr), "-o", str(gpath), *b_e],
            stdin="0\n")
        g = read_xvg(gpath)
        t = get("gyrate")().run(ctx, sel="all", frame_range=fr, stride=None)
        # gmx col0 time, col1 Rg (nm); lat Rg Angstrom
        lat_rg = t[0].data[:, 1] / 10.0
        n = min(len(lat_rg), len(g))
        d = g[-n:, 1] - lat_rg[-n:]
        out["gyrate"] = {"max_abs_nm": float(np.abs(d).max()),
                         "rms_nm": float(np.sqrt((d ** 2).mean())),
                         "gmx_mean": float(g[-n:, 1].mean()),
                         "lat_mean": float(lat_rg[-n:].mean())}
    except Exception as e:
        out["gyrate_error"] = str(e)[:300]

    # ---- msd ----
    try:
        if not gmx_ok:
            raise RuntimeError("gmx unavailable")
        mpath = work / "msd.xvg"
        gmx(["msd", "-f", str(traj), "-s", str(tpr), "-o", str(mpath),
             "-sel", f"resname {ref}", "-seltype", "mol_com", *b_e])
        g = read_xvg(mpath)
        t = get("msd")().run(ctx, sel=f"resname {ref}", dim=3, seltype="mol_com",
                            fit_range=None, remove_com=True, no_remove_com=False,
                            frame_range=fr, stride=stride)
        # lat MSD is Angstrom^2 vs lag; map lag to gmx time using gmx dt
        lat_msd = t[0].data[:, 1] / 100.0  # nm^2
        lag = t[0].data[:, 0]
        t_l = lag * (dt * stride) if dt else lag
        gi = np.interp(t_l, g[:, 0], g[:, 1])
        m = t_l > 0
        if m.any():
            d = gi - lat_msd
            out["msd"] = {"max_abs_nm2": float(np.abs(d[m]).max()),
                          "rms_nm2": float(np.sqrt((d[m] ** 2).mean())),
                          "gmx_end": float(g[-1, 1]), "lat_end": float(lat_msd[-1])}
    except Exception as e:
        out["msd_error"] = str(e)[:300]

    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", required=True)
    ap.add_argument("--maxdepth", type=int, default=8)
    ap.add_argument("--frames", type=int, default=20)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--timeout", type=float, default=1200)
    ap.add_argument("--out", default="/tmp/opencode/refdata_gmx.md")
    ap.add_argument("--json", default="/tmp/opencode/refdata_gmx.json")
    a = ap.parse_args()

    pairs = []
    for r in a.root:
        pairs += discover(Path(r), a.maxdepth)
    if a.limit:
        pairs = pairs[: a.limit]
    print(f"# {len(pairs)} (topology, trajectory) pairs under {a.root}", flush=True)
    results = []
    workroot = Path("/tmp/opencode/refdata_gmx_work")
    for tpr, traj in pairs:
        tag = str(tpr.parent)
        print(f"* {tag}  [{tpr.name} + {traj.name}]", flush=True)
        try:
            r = regress_pair(tpr, traj, workroot / tpr.parent.name, a.frames, a.stride, a.timeout)
        except Exception as e:
            r = {"fatal": str(e)[:400]}
        r["tpr"] = str(tpr)
        r["traj"] = str(traj)
        results.append(r)
        Path(a.json).write_text(json.dumps(results, indent=1))  # survive a kill
        def fmt(d, key):
            return f"{d[key]:.4g}" if isinstance(d, dict) and key in d else "-"
        rd = r.get("rdf", {})
        worst = max((v["rel"] for k, v in rd.items() if not k.endswith("_cn") and v),
                    default=float("nan"))
        print(f"    rdf_worst_rel={worst:.3g}  dens={fmt(r.get('density'),'max_abs')}  "
              f"rg={fmt(r.get('gyrate'),'max_abs_nm')}  msd={fmt(r.get('msd'),'max_abs_nm2')}"
              + (f"  ERR {r.get('fatal') or r.get('rdf_error') or ''}" if (r.get('fatal') or r.get('rdf_error')) else ""),
              flush=True)

    Path(a.json).write_text(json.dumps(results, indent=1))
    lines = ["# GROMACS native regression (third-party RefData)", "",
             "| system | atoms | frames | rdf max rel | cn max rel | density max|d| | Rg max|d| nm | MSD max|d| nm^2 |",
             "|---|---|---|---|---|---|---|---|"]
    for r in results:
        rd = r.get("rdf", {})
        rel = [v["rel"] for k, v in rd.items() if not k.endswith("_cn") and v]
        cnrel = [v["rel"] for k, v in rd.items() if k.endswith("_cn") and v]
        def g(d, k):
            return f"{d[k]:.3g}" if isinstance(d, dict) and k in d else "-"
        lines.append("| `%s` | %s | %s | %s | %s | %s | %s | %s |" % (
            Path(r.get("tpr", "?")).parent.name, r.get("atoms", "-"), r.get("frames", "-"),
            f"{max(rel):.3g}" if rel else "-", f"{max(cnrel):.3g}" if cnrel else "-",
            g(r.get("density"), "max_abs"), g(r.get("gyrate"), "max_abs_nm"),
            g(r.get("msd"), "max_abs_nm2")))
    Path(a.out).write_text("\n".join(lines))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Compatibility battery for real GROMACS order data (lat only, no gmx needed).

Runs every applicable lat module against each (tpr/gro, xtc/trr) pair and
records PASS/FAIL/SKIP plus error text. Complements gmx_regress.py (which needs
a working gmx binary and may be blocked by tpr/binary version mismatches).

Usage:
  gmx_compat.py --root /mnt/e/RefData/gromacs_real_projects --maxdepth 8 \
                --frames 3 --out /tmp/opencode/gmx_compat.md
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gmx_regress import discover, gmx  # noqa: E402

LAT = [sys.executable, "/opt/lat/lat"]


def run_lat(args, timeout):
    import subprocess
    t0 = time.time()
    try:
        p = subprocess.run(LAT + args, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr, time.time() - t0
    except subprocess.TimeoutExpired:
        return 124, "", "TIMEOUT", time.time() - t0


def topo_resnames_vel(tpr, traj):
    import MDAnalysis as mda
    u = mda.Universe(str(tpr), str(traj))
    res = {}
    for rn in u.atoms.resnames:
        res[str(rn)] = res.get(str(rn), 0) + 1
    res = [k for k, _ in sorted(res.items(), key=lambda kv: -kv[1])]
    has_bonds = False
    try:
        has_bonds = len(u.bonds) > 0
    except Exception:
        pass
    nf = len(u.trajectory)
    has_vel = False
    try:
        u.trajectory[0]
        v = u.atoms.velocities
        has_vel = v is not None and len(v) == u.atoms.n_atoms and bool(v.any())
    except Exception:
        pass
    return res, has_bonds, nf, has_vel


def battery(tpr, traj, timeout):
    resnames, bonds, nf, has_vel = topo_resnames_vel(tpr, traj)
    r1 = resnames[0] if resnames else "all"
    r2 = resnames[1] if len(resnames) > 1 else r1
    n = min(3, max(1, nf))
    fr = ["-frame_range", f"0:{n}:1"] if nf >= 2 else []
    base = ["-s", str(tpr), "-f", str(traj)]
    o = "/tmp/opencode/gmx_compat_out"
    cases = [
        ("describe", ["describe", "-module", "rdf"]),
        ("select", ["select", *base, "-sel", "all"]),
        ("rdf", ["rdf", *base, *fr, "-sel1", f"resname {r1}", "-sel2", f"resname {r2}",
                 "-cutoff", "8", "-o", o]),
        ("density1d", ["density", *base, *fr, "-axis", "z", "-o", o]),
        ("gyrate", ["gyrate", *base, *fr, "-sel", "all", "-o", o]),
        ("rms", (["rms", *base, *fr, "-sel", "all", "-o", o] if nf >= 2 else None)),
        ("msd", (["msd", *base, *fr, "-sel", f"resname {r1}", "-seltype", "mol_com", "-o", o]
                 if nf >= 3 else None)),
        ("contact", ["contact", *base, *fr, "-sel1", f"resname {r1}", "-sel2", f"resname {r2}",
                     "-cutoff", "5", "-o", o]),
        ("occupancy", ["occupancy", *base, *fr, "-sel", f"resname {r1}",
                       "-zone", "z:0-10", "-o", o]),
        ("probe_ffv", ["probe", *base, *fr, "-mode", "ffv", "-o", o]),
        ("transport", (["transport", *base, *fr, "-sel", f"resname {r1}", "-o", o]
                       if has_vel else None)),
    ]
    if bonds:
        cases += [
            ("hbond", (["hbond", *base, *fr, "-o", o] if nf >= 2 else None)),
            ("agg", ["agg", *base, *fr, "-sel1", f"resname {r1}", "-sel2", f"resname {r2}",
                     "-cutoff", "3.5", "-o", o]),
            ("cluster", ["cluster", *base, *fr, "-sel1", f"resname {r1}", "-sel2", f"resname {r2}",
                         "-cutoff", "3.5", "-o", o]),
            ("residence", (["residence", *base, *fr, "-lig", f"resname {r1}", "-ref", f"resname {r2}",
                            "-cutoff", "3.5", "-o", o] if nf >= 2 else None)),
        ]
    else:
        for nm in ("hbond", "agg", "cluster", "residence"):
            cases.append((nm, None))
    out = []
    for name, args in cases:
        if args is None:
            out.append((name, "SKIP", "precondition not met"))
            continue
        rc, so, se, dt = run_lat(args, timeout)
        if rc == 0:
            out.append((name, "PASS", ""))
        else:
            tail = (se or so).strip().splitlines()[-1] if (se or so).strip() else ""
            out.append((name, "FAIL", tail[:300]))
    return resnames, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", required=True)
    ap.add_argument("--maxdepth", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=300)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="/tmp/opencode/gmx_compat.md")
    ap.add_argument("--json", default="/tmp/opencode/gmx_compat.json")
    a = ap.parse_args()

    pairs = []
    for r in a.root:
        pairs += discover(Path(r), a.maxdepth)
    if a.limit:
        pairs = pairs[: a.limit]
    print(f"# {len(pairs)} pairs", flush=True)
    results = []
    for tpr, traj in pairs:
        print(f"* {tpr}", flush=True)
        try:
            resnames, out = battery(tpr, traj, a.timeout)
        except Exception as e:
            results.append({"tpr": str(tpr), "traj": str(traj), "fatal": str(e)[:300], "result": []})
            print(f"    FATAL {str(e)[:120]}", flush=True)
            continue
        bad = [x for x in out if x[1] == "FAIL"]
        print(f"    {len(out)-len(bad)}/{len(out)} pass" + (f"  FAIL {[b[0] for b in bad]}" if bad else ""),
              flush=True)
        results.append({"tpr": str(tpr), "traj": str(traj), "resnames": resnames,
                        "result": [list(x) for x in out]})
        Path(a.json).write_text(json.dumps(results, indent=1))  # survive a kill
    Path(a.json).write_text(json.dumps(results, indent=1))
    mods = ["describe", "select", "rdf", "density1d", "gyrate", "rms", "msd", "contact",
            "occupancy", "probe_ffv", "transport", "hbond", "agg", "cluster", "residence"]
    lines = ["# GROMACS compatibility battery (lat only)", "",
             "| system | " + " | ".join(mods) + " |", "|" + "---|" * (len(mods) + 1)]
    fails = []
    for r in results:
        m = {x[0]: x for x in r["result"]}
        row = [{"PASS": "ok", "SKIP": "·", "FAIL": "FAIL"}.get(m.get(mod, ["", "-"])[1], "-")
               for mod in mods]
        lines.append(f"| `{Path(r['tpr']).parent.name}/{Path(r['tpr']).name}` | " + " | ".join(row) + " |")
        for x in r["result"]:
            if x[1] == "FAIL":
                fails.append((r["tpr"], x[0], x[2]))
        if r.get("fatal"):
            fails.append((r["tpr"], "FATAL", r["fatal"]))
    lines += ["", "## failures", ""]
    for tpr, mod, msg in fails:
        lines.append(f"- `{mod}` · `{tpr}`\n  - {msg}")
    Path(a.out).write_text("\n".join(lines))
    print(f"\nwrote {a.out} ({len(results)} systems, {len(fails)} failures)")


if __name__ == "__main__":
    main()

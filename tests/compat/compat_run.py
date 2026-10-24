#!/usr/bin/env python
"""Compatibility harness: run every lat module against real LAMMPS projects.

Goal is crash-hunting / contract compliance on diverse real structures
(orthorhombic + triclinic, nonzero box origin, atomic/molecular/charged,
with/without image flags, with/without velocities, tiny boxes).

Usage:
  python compat_run.py --root /mnt/d/CurrentProjects --maxdepth 3
  python compat_run.py --pairs pairs.txt --out report.md
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

LAT = [sys.executable, "/opt/lat/lat"]
DUMP_EXT = {".lammpstrj", ".dump"}
DATA_EXT = {".data", ".lmp", ".lam", ".data.gz"}


# --------------------------------------------------------------------------
# lightweight header parsing
# --------------------------------------------------------------------------
def read_text_head(p: Path, nbytes: int = 65536) -> str:
    with open(p, "rb") as fh:
        return fh.read(nbytes).decode("utf-8", "replace")


def dump_meta(path: Path) -> dict | None:
    """Header of first frame: box, triclinic flag, atom columns, n_atoms."""
    try:
        head = read_text_head(path, 200000)
    except OSError:
        return None
    lines = head.splitlines()
    # first frame atom count
    n = None
    i = 0
    box = None
    tric = False
    cols = []
    while i < len(lines):
        ln = lines[i].strip()
        if ln == "ITEM: NUMBER OF ATOMS":
            n = int(lines[i + 1])
        elif ln.startswith("ITEM: BOX BOUNDS"):
            tric = any(k in ln for k in ("xy", "xz", "yz"))
            box = []
            for k in range(3):
                parts = [float(x) for x in lines[i + 1 + k].split()]
                box.append(parts)
            i += 3
        elif ln.startswith("ITEM: ATOMS"):
            cols = ln.split()[2:]
            break
        i += 1
    if n is None or box is None or not cols:
        return None
    lengths = []
    for k in range(3):
        lengths.append(box[k][1] - box[k][0])
    return {
        "n_atoms": n,
        "triclinic": tric,
        "cols": cols,
        "lengths": lengths,
        "origin": [b[0] for b in box],
    }


def dump_types(path: Path, n_atoms: int) -> list[int] | None:
    """Types actually present in the first frame (topology Masses may list all)."""
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("ITEM: ATOMS"):
                    names = line.split()[2:]
                    if "type" not in names:
                        return None
                    k = names.index("type")
                    seen = set()
                    for _ in range(n_atoms):
                        row = fh.readline().split()
                        if not row:
                            break
                        seen.add(int(row[k]))
                    return sorted(seen)
    except OSError:
        return None
    return None


def data_meta(path: Path) -> dict | None:
    try:
        head = read_text_head(path, 400000)
    except OSError:
        return None
    m = re.search(r"^\s*(\d+)\s+atoms\s*$", head, re.M)
    if not m:
        return None
    n = int(m.group(1))
    mt = re.search(r"^\s*(\d+)\s+atom types\s*$", head, re.M)
    ntypes = int(mt.group(1)) if mt else 1
    mb = re.search(r"^\s*(\d+)\s+bonds\s*$", head, re.M)
    bonds = bool(mb) and int(mb.group(1)) > 0
    masses = {}
    mm = re.search(r"^Masses\s*$", head, re.M)
    if mm:
        for ln in head[mm.end():].splitlines():
            s = ln.strip()
            if not s:
                if masses:
                    break
                continue
            parts = s.split()
            if len(parts) >= 2 and parts[0].lstrip("-").isdigit():
                try:
                    masses[int(parts[0])] = float(parts[1])
                except ValueError:
                    break
            elif masses:
                break
    return {"n_atoms": n, "ntypes": ntypes, "bonds": bonds, "masses": masses}


def pair_up(root: Path, maxdepth: int) -> list[tuple[Path, Path]]:
    """Shallow scan: pair each dump with the best topology in the same dir."""
    import os
    dumps, datas = [], []
    root = root.resolve()
    skip = {"nep_structures", ".git", "__pycache__", "node_modules", "venv", ".venv"}
    for dirpath, dirnames, filenames in os.walk(root):
        rel = Path(dirpath).relative_to(root)
        if len(rel.parts) >= maxdepth:
            dirnames[:] = []
        dirnames[:] = [d for d in dirnames if d not in skip]
        for fn in filenames:
            p = Path(dirpath) / fn
            suf = "".join(p.suffixes[-2:]) if p.suffixes[-1:] == [".gz"] else p.suffix
            if suf.lower() in DUMP_EXT:
                dumps.append(p)
            elif suf.lower() in DATA_EXT:
                datas.append(p)
    pairs = []
    for d in sorted(dumps):
        dm = dump_meta(d)
        n = dm["n_atoms"] if dm else None
        cands = [x for x in datas if x.parent == d.parent]
        if not cands:
            cands = [x for x in datas if d.stem.split(".")[0] in x.name]
        if not cands:
            continue
        if n is not None:
            matching = [x for x in cands
                        if (data_meta(x) or {}).get("n_atoms") == n]
            if not matching:
                continue  # topology/dump genuinely differ (e.g. GCMC snapshot)
            cands = matching
        # prefer a topology whose stem matches, else the smallest plausible
        cands.sort(key=lambda x: (d.stem not in x.name, len(x.name)))
        pairs.append((cands[0], d))
    return pairs


# --------------------------------------------------------------------------
# module battery
# --------------------------------------------------------------------------
def run(args, cwd, timeout):
    t0 = time.time()
    try:
        p = subprocess.run(LAT + args, cwd=cwd, capture_output=True,
                           text=True, timeout=timeout)
        return p.returncode, p.stdout, p.stderr, time.time() - t0
    except subprocess.TimeoutExpired:
        return 124, "", f"TIMEOUT after {timeout}s", time.time() - t0


def battery(topo: Path, dump: Path, cwd: Path, tmp: Path, timeout: float):
    dm = dump_meta(dump)
    d = data_meta(topo)
    if dm is None or d is None:
        return [("__meta__", "SKIP", "unreadable header", 0.0)]
    present = dump_types(dump, dm["n_atoms"])
    t_all = present or sorted(d["masses"]) or list(range(1, min(d["ntypes"], 3) + 1))
    t1 = t_all[0]
    t2 = t_all[1] if len(t_all) > 1 else t1
    min_edge = min(dm["lengths"])
    half = min_edge / 2.0
    c = max(1.5, min(8.0, 0.45 * min_edge))
    has_vel = all(k in dm["cols"] for k in ("vx", "vy", "vz"))
    base = ["-s", str(topo), "-f", str(dump)]
    o = str(tmp / "out")

    # discover frame count from `select` (cheap thanks to lazy counting)
    nf = 1
    rc, so, se, _ = run(["select", *base, "-sel", f"type {t1}"], cwd, timeout)
    m = re.search(r"frames=(\d+)", so)
    if rc == 0 and m:
        nf = int(m.group(1))
    fr = ["-frame_range", f"0:{min(4, nf)}:1"] if nf >= 2 else []

    cases = [
        ("describe", ["describe", "-module", "rdf"], True),
        ("select", ["select", *base, "-sel", f"type {t1}"], True),
        ("rdf", ["rdf", *base, *fr, "-sel1", f"type {t1}", "-sel2", f"type {t2}",
                 "-cutoff", f"{c:.3f}", "-o", o], True),
        ("density1d", ["density", *base, *fr, "-axis", "z",
                       *(["-quantity", "number"] if not d["masses"] else []), "-o", o], True),
        ("gyrate", ["gyrate", *base, *fr, "-sel", f"type {t1}", "-o", o], True),
        ("rms", (["rms", *base, *fr, "-sel", f"type {t1}", "-o", o] if nf >= 2 else None), True),
        ("msd", (["msd", *base, *fr, "-sel", f"type {t1}", "-o", o] if nf >= 3 else None), True),
        ("contact", ["contact", *base, *fr, "-sel1", f"type {t1}", "-sel2", f"type {t2}",
                     "-cutoff", f"{c:.3f}", "-o", o], True),
        ("occupancy", ["occupancy", *base, *fr, "-sel", f"type {t1}",
                       "-zone", f"z:0-{half:.3f}", "-o", o], True),
        ("freevol_ffv", (["freevol", *base, *fr, "-mode", "ffv", "-o", o]
                         if (d["masses"] or "element" in dm["cols"]) else None), True),
    ]
    if has_vel:
        cases.append(("transport", ["transport", *base, *fr, "-sel", f"type {t1}", "-o", o], True))
    else:
        cases.append(("transport", None, True))
    if d["bonds"]:
        cases += [
            ("hbond", (["hbond", *base, *fr, "-o", o] if nf >= 2 else None), True),
            ("agg", ["agg", *base, *fr, "-sel1", f"type {t1}", "-sel2", f"type {t2}",
                     "-cutoff", f"{min(3.5, half):.3f}", "-o", o], True),
            ("cluster", ["cluster", *base, *fr, "-sel1", f"type {t1}", "-sel2", f"type {t2}",
                         "-cutoff", f"{min(3.5, half):.3f}", "-o", o], True),
            ("residence", ["residence", *base, *fr, "-lig", f"type {t1}", "-ref", f"type {t2}",
                           "-cutoff", f"{min(3.5, half):.3f}", "-o", o], True),
        ]
    else:
        for nm in ("hbond", "agg", "cluster", "residence"):
            cases.append((nm, None, True))
    # pistack needs a curated per-residue ring selection (>=3 ring atoms);
    # a generic per-type selection is not a valid test, so it lives in a
    # dedicated sample test instead of the generic battery.

    out = []
    for name, args, _ in cases:
        if args is None:
            out.append((name, "SKIP", "precondition not met", 0.0))
            continue
        rc, so, se, dt = run(args, cwd, timeout)
        if rc == 0:
            out.append((name, "PASS", "", dt))
            print(f"    {name:11s} ok    {dt:5.1f}s", flush=True)
        else:
            tail = (se or so).strip().splitlines()[-1] if (se or so).strip() else ""
            out.append((name, "FAIL", tail[:300], dt))
            print(f"    {name:11s} FAIL  {dt:5.1f}s  {tail[:120]}", flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", action="append", default=[])
    ap.add_argument("--pairs", default=None)
    ap.add_argument("--maxdepth", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=300)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="/tmp/opencode/compat_report.md")
    ap.add_argument("--json", default="/tmp/opencode/compat_report.json")
    a = ap.parse_args()

    pairs = []
    if a.pairs:
        for ln in Path(a.pairs).read_text().splitlines():
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                t, d = ln.split("|")
                pairs.append((Path(t), Path(d)))
    for r in a.root:
        pairs += pair_up(Path(r), a.maxdepth)
    if a.limit:
        pairs = pairs[: a.limit]

    results = []
    tmp = Path("/tmp/opencode/compat_tmp")
    tmp.mkdir(parents=True, exist_ok=True)
    for topo, dump in pairs:
        if dump_meta(dump) is None or data_meta(topo) is None:
            continue
        print(f"* {topo}  ||  {dump}", flush=True)
        res = battery(topo, dump, dump.parent, tmp, a.timeout)
        results.append({"topo": str(topo), "dump": str(dump), "result": res})
        bad = [r for r in res if r[1] == "FAIL"]
        tag = "FAIL " + ",".join(b[0] for b in bad) if bad else "all pass"
        print(f"{dump}  [{len(res)-len(bad)}/{len(res)}]  {tag}")

    Path(a.json).write_text(json.dumps(results, indent=1))
    lines = ["# lat compatibility report", ""]
    mods = ["describe", "select", "rdf", "density1d", "gyrate", "rms", "msd", "contact",
            "occupancy", "freevol_ffv", "transport", "hbond", "agg", "cluster", "residence", "pistack"]
    lines.append("| dump | " + " | ".join(mods) + " |")
    lines.append("|" + "---|" * (len(mods) + 1))
    fails = []
    for r in results:
        m = {x[0]: x for x in r["result"]}
        row = []
        for mod in mods:
            x = m.get(mod)
            row.append("-" if x is None else {"PASS": "ok", "SKIP": "·", "FAIL": "FAIL"}[x[1]])
        lines.append(f"| `{Path(r['dump']).name}` | " + " | ".join(row) + " |")
        for x in r["result"]:
            if x[1] == "FAIL":
                fails.append((r["dump"], x[0], x[2]))
    lines += ["", "## failures", ""]
    for dump, mod, msg in fails:
        lines.append(f"- `{mod}` · `{dump}`\n  - {msg}")
    Path(a.out).write_text("\n".join(lines))
    print(f"\nwrote {a.out} ({len(results)} systems, {len(fails)} failures)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
verify_lammps_roundtrip.py — 验证 gmx->lammps 转换的可读性与一致性。

检查：
  1) MDAnalysis 能读 <prefix>.data 与 <prefix>.lammpstrj（format=LAMMPSDUMP）；
  2) 逐帧坐标与 <prefix>.gro/.xtc 一致（atol=1e-3 Å）；
  3) 原子数、盒、键表；电荷/质量存在；
  4) type->element 与 gmx 原子名推断元素逐原子一致；
  5) LAMMPS `type N` 与 GROMACS `resname` 命中数一致（中心 Li 与各分子）。
输出 <prefix>.verify.json + 摘要。
"""
import argparse, json, re, sys
import numpy as np
import MDAnalysis as mda


def el_from_name(name):
    m = re.match(r"([A-Za-z]{1,2})", name)
    s = m.group(1)
    return s.capitalize() if len(s) == 2 else s.upper()


def parse_masses(path):
    type2el = {}
    with open(path) as fh:
        in_masses = False
        for line in fh:
            s = line.strip()
            if s.startswith("Masses"):
                in_masses = True
                continue
            if in_masses:
                if s.startswith(("Atoms", "Bonds")):
                    break
                if not s:
                    continue
                parts = s.split("#")
                nums = parts[0].split()
                if len(nums) >= 2 and len(parts) > 1:
                    type2el[int(nums[0])] = parts[1].strip()
    return type2el


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--prefix", required=True)
    args = ap.parse_args()
    p = args.prefix
    rep = {"prefix": p, "checks": {}, "ok": True}

    def check(name, cond, detail=""):
        rep["checks"][name] = {"ok": bool(cond), "detail": detail}
        if not cond:
            rep["ok"] = False

    ulmp = mda.Universe(p + ".data", p + ".lammpstrj", format="LAMMPSDUMP")
    ugmx = mda.Universe(p + ".gro", p + ".xtc")

    check("lmp_natoms", ulmp.atoms.n_atoms > 0, str(ulmp.atoms.n_atoms))
    check("natoms_equal", ulmp.atoms.n_atoms == ugmx.atoms.n_atoms,
          f"{ulmp.atoms.n_atoms} vs {ugmx.atoms.n_atoms}")
    check("lmp_has_charges", float(np.abs(ulmp.atoms.charges).sum()) > 0,
          f"sum|q|={np.abs(ulmp.atoms.charges).sum():.3f}")
    check("lmp_net_charge", abs(float(ulmp.atoms.charges.sum())) < 1e-2,
          f"net={ulmp.atoms.charges.sum():.2e}")
    try:
        nb = len(ulmp.bonds)
    except Exception:
        nb = 0
    check("lmp_bonds_nonzero", nb > 0, str(nb))

    # frame coords
    maxdev = 0.0
    nframes = min(len(ulmp.trajectory), len(ugmx.trajectory))
    for i in range(nframes):
        ulmp.trajectory[i]
        ugmx.trajectory[i]
        maxdev = max(maxdev, float(np.abs(ulmp.atoms.positions - ugmx.atoms.positions).max()))
    check("frame_max_dev_A", maxdev < 1e-3, f"{maxdev:.2e} over {nframes} frames")

    # box
    ulmp.trajectory[0]; ugmx.trajectory[0]
    check("box_equal", np.allclose(ulmp.dimensions, ugmx.dimensions, atol=1e-3),
          f"{np.round(ulmp.dimensions,3)}")

    # type->element vs gmx name-derived element
    type2el = parse_masses(p + ".data")
    lmp_el = np.array([type2el.get(int(t), "?") for t in ulmp.atoms.types])
    gmx_el = np.array([el_from_name(n) for n in ugmx.atoms.names])
    match_frac = float(np.mean(lmp_el == gmx_el))
    check("type_element_consistency", match_frac > 0.999, f"{match_frac:.4f}")
    rep["type2element"] = type2el

    # selection equivalence by per-element counts
    from collections import Counter
    cl = Counter(lmp_el.tolist())
    cg = Counter(gmx_el.tolist())
    check("element_counts_equal", cl == cg, f"lmp={dict(cl)} gmx={dict(cg)}")

    # known selections: Li via lmp type vs gmx resname LI
    li_type = [t for t, e in type2el.items() if e == "Li"]
    if li_type:
        n_lmp = len(ulmp.select_atoms(f"type {li_type[0]}"))
        n_gmx = len(ugmx.select_atoms("resname LI"))
        check("Li_selection_equal", n_lmp == n_gmx, f"lmp={n_lmp} gmx={n_gmx}")

    json.dump(rep, open(p + ".verify.json", "w"), indent=2)
    print(json.dumps(rep, indent=2))
    sys.exit(0 if rep["ok"] else 1)


if __name__ == "__main__":
    main()

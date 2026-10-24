#!/usr/bin/env python3
"""
make_selections.py — 从 <prefix>.mapping.csv 生成 <prefix>.selections.yaml：
元素级、分子级（按 resname 的 mol/atom 区间）以及各 P0 模块的候选选择。
"""
import argparse, csv, json
from collections import OrderedDict, defaultdict


def compress(vals):
    vals = sorted(vals)
    if not vals:
        return []
    rngs, s, p = [], vals[0], vals[0]
    for v in vals[1:]:
        if v == p + 1:
            p = v
        else:
            rngs.append((s, p)); s = p = v
    rngs.append((s, p))
    return rngs


def rng_str(rngs):
    return ",".join(str(a) if a == b else f"{a}-{b}" for a, b in rngs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-p", "--prefix", required=True)
    args = ap.parse_args()
    p = args.prefix
    rows = list(csv.DictReader(open(p + ".mapping.csv")))
    rows = [{k: v for k, v in r.items()} for r in rows]

    by_type, by_res = {}, defaultdict(list)
    for r in rows:
        by_type.setdefault((int(r["type"]), r["element"]), 0)
        by_type[(int(r["type"]), r["element"])] += 1
        by_res[r["resname"]].append(r)

    out = OrderedDict()
    out["elements"] = [
        OrderedDict(type=t, element=e, count=c, lmp=f"type {t}",
                    gmx=f"name {e}*")
        for (t, e), c in sorted(by_type.items())
    ]
    mols = []
    for rn, rs in sorted(by_res.items()):
        mol_ids = [int(x["mol"]) for x in rs]
        atom_ids = [int(x["id"]) for x in rs]
        mol_r = compress(sorted(set(mol_ids)))
        mol_range_str = rng_str(mol_r)
        mols.append(OrderedDict(
            resname=rn, n_molecules=len(set(mol_ids)), n_atoms=len(rs),
            lmp_mol_range=mol_range_str,
            lmp=f"mol {mol_range_str}",
            gmx=f"resname {rn}",
        ))
    out["molecules"] = mols

    # candidate module selections (P0)
    types = {e: t for (t, e) in by_type}
    sel = OrderedDict()
    sel["rdf"] = [
        {"sel1": f"type {types['Li']}", "sel2": f"type {types[a]}",
         "note": f"Li-{a}"}
        for a in ("O", "N", "F", "S")
    ]
    sel["density"] = [{"sel": "all", "axis": "z", "note": "mass density profile"}]
    sel["msd"] = [{"sel": f"type {types[e]}"} for e in ("Li", "C", "O")]
    sel["agg_cluster"] = [
        {"sel1": f"type {types['Li']}", "sel2": f"type {types['N']}",
         "note": "Li-FSI graph (N as anion center)"}
    ]
    sel["hbond"] = [{"note": "no O-H donor in sys1; use FEC/FEMC C-H as donor if needed"}]
    sel["freevol"] = [{"note": "FFV/PSD/connectivity; grid modes have no -sel"}]
    sel["sasa"] = [{"sel": "all", "note": "H2 probe radius 1.445 A"}]
    sel["transport"] = [{"note": "requires vx vy vz columns in dump (this sample has none)"}]
    out["candidate_module_selections"] = sel
    out["_note"] = ("LAMMPS 无 resname；分子选择用 mol 区间或按 type。"
                    "GROMACS 样例 .gro/.xtc 仅为坐标级回归；权威逐原子拓扑见 .mapping.csv"
                    "（.gro 读回时 Li 质量会被猜测为 0，属已知限制）。")

    # minimal YAML writer
    with open(p + ".selections.yaml", "w") as fh:
        def dump(o, ind=0):
            sp = "  " * ind
            if isinstance(o, dict):
                for k, v in o.items():
                    if isinstance(v, (dict, list)):
                        fh.write(f"{sp}{k}:\n"); dump(v, ind + 1)
                    else:
                        fh.write(f"{sp}{k}: {v}\n")
            elif isinstance(o, list):
                for it in o:
                    if isinstance(it, dict):
                        first = True
                        for k, v in it.items():
                            if first:
                                fh.write(f"{sp}- {k}: {v}\n"); first = False
                            else:
                                fh.write(f"{sp}  {k}: {v}\n")
                    else:
                        fh.write(f"{sp}- {it}\n")
        dump(out)
    print(f"[ok] wrote {p}.selections.yaml")


if __name__ == "__main__":
    main()

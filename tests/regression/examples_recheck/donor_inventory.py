#!/usr/bin/env python3
"""Donor inventory of a LAMMPS data file: which heavy atoms carry H, by element.

This is the evidence behind the 0.2.0 hbond donor rule (a C-H is not a donor,
mirroring `gmx hbond`, which prints "Found 0 donors ... No Donors found" for a
donor-free group).

    python3 donor_inventory.py <data> [sel1 ...]
"""
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/opt/lat/src")
from lat.config import load
from lat.core.io import open_system
from lat.core.selection import parse

POLAR = {"N", "O", "F", "S", "Cl", "Br", "I"}


def main():
    path = Path(sys.argv[1])
    sel = " ".join(sys.argv[2:]) or "all"
    cfg = load()
    topo, _ = open_system(str(path), None, cfg=cfg)
    elem = topo.type_element
    el = [elem.get(int(t)) for t in topo.types]
    adj = [[] for _ in range(topo.n_atoms)]
    for a, b in topo.bonds:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    group = set(parse(topo, sel).indices.tolist())
    heavy_h = Counter()
    polar_donor = Counter()
    for i in range(topo.n_atoms):
        if el[i] == "H" or not adj[i]:
            continue
        hs = [j for j in adj[i] if el[j] == "H" and j in group]
        if hs and i in group:
            heavy_h[el[i]] += 1
            if el[i] in POLAR:
                polar_donor[el[i]] += 1
    print(f"file     : {path}")
    print(f"selection: {sel}   ({len(group)} atoms)")
    print(f"atoms    : {topo.n_atoms}, bonds {0 if topo.bonds is None else len(topo.bonds)}")
    print("heavy atom with bonded H (inside the selection), by element:")
    for k, v in sorted(heavy_h.items()):
        mark = "donor" if k in POLAR else "NOT a donor (gmx/lat)"
        print(f"  {k:>3s}: {v:6d}   {mark}")
    n_polar = sum(polar_donor.values())
    n_any = sum(heavy_h.values())
    print(f"total with H: {n_any}; polar donors: {n_polar}"
          + ("  -> hbond refuses this group" if n_polar == 0 else ""))


if __name__ == "__main__":
    main()

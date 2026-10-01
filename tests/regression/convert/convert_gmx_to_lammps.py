#!/usr/bin/env python3
"""
convert_gmx_to_lammps.py — 把 GROMACS 体系（tpr/gro + xtc）转成 lat 可读的
LAMMPS 格式：<prefix>.data（atom_style full）+ <prefix>.lammpstrj，并可选写出
缩减版 gmx 样例（.gro/.xtc）与映射/选择清单。

设计目标：为 lat 的两级输入回归提供自包含样例，替代 InterMol/无头 VMD。

用法示例：
  python convert_gmx_to_lammps.py \
      -s sys1_npt.tpr -f sys1_npt.xtc \
      -o samples/sys1_reduced \
      --per-resname 200 --frames 5

约定：
  * 元素由质量按 ±0.2% 匹配 IUPAC 标准原子量（ase.data），无解/多解 -> 报错。
  * 一个元素一个 LAMMPS type；type 按原子序数升序编号。
  * atom_style full: id mol type q x y z。
  * 时间=帧号；dump 的 TIMESTEP 写帧序号。
"""
import argparse, csv, json, os, sys
from collections import OrderedDict

import numpy as np
import MDAnalysis as mda

try:
    from ase.data import chemical_symbols, atomic_masses
except Exception:
    chemical_symbols, atomic_masses = None, None

TOL = 0.002  # ±0.2%


def build_element_mass():
    if chemical_symbols is None:
        raise RuntimeError("ase is required to build the element mass table")
    d = OrderedDict()
    for z, sym in enumerate(chemical_symbols):
        if z == 0:
            continue
        d[sym] = float(atomic_masses[z])
    return d


def mass_to_element(mass, table, tol=TOL):
    hits = [e for e, m in table.items() if abs(m - mass) <= tol * m]
    if len(hits) == 1:
        return hits[0]
    return None if not hits else "AMBIGUOUS:" + ",".join(hits)


def lammps_box_from_dimensions(dim, tol=0.1):
    """dim = [a,b,c,alpha,beta,gamma] (MDAnalysis), return box lines + flags."""
    a, b, c, al, be, ga = dim
    ortho = (abs(al - 90) < 1e-3 and abs(be - 90) < 1e-3 and abs(ga - 90) < 1e-3)
    if ortho:
        return (0.0, a), (0.0, b), (0.0, c), None, True
    al, be, ga = np.radians([al, be, ga])
    lx = a
    xy = b * np.cos(ga)
    xz = c * np.cos(be)
    ly = np.sqrt(b * b - xy * xy)
    yz = (b * c * np.cos(al) - xy * xz) / ly
    lz = np.sqrt(c * c - xz * xz - yz * yz)
    return (0.0, lx), (0.0, ly), (0.0, lz), (xy, xz, yz), False


def pick_residues(u, per_resname, seed=0):
    """Pick first N residues per resname (N=None -> all). Returns sorted resindices."""
    chosen = []
    resnames = u.residues.resnames
    uniq = sorted(set(resnames))
    for rn in uniq:
        idx = np.where(resnames == rn)[0]
        if per_resname is not None and per_resname < len(idx):
            idx = idx[:per_resname]
        chosen.append(idx)
    return np.sort(np.concatenate(chosen)) if chosen else np.array([], dtype=int)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-s", "--top", required=True, help="GROMACS topology (tpr/gro/pdb)")
    ap.add_argument("-f", "--traj", default=None, help="GROMACS trajectory (xtc/trr/dcd)")
    ap.add_argument("-o", "--prefix", required=True, help="output prefix")
    ap.add_argument("--per-resname", type=int, default=None, help="keep first N residues per resname (default all)")
    ap.add_argument("--frames", type=int, default=1, help="number of frames to write (default 1)")
    ap.add_argument("--start", type=int, default=0, help="start frame index")
    ap.add_argument("--stride", type=int, default=None, help="frame stride (default spread over trajectory)")
    ap.add_argument("--no-bonds", action="store_true", help="omit Bonds section")
    ap.add_argument("--no-gmx-sample", action="store_true", help="do not write reduced .gro/.xtc")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.prefix)) or ".", exist_ok=True)
    out = args.prefix

    u = mda.Universe(args.top, args.traj) if args.traj else mda.Universe(args.top)
    nframes = len(u.trajectory)
    u.trajectory[0]
    dim = u.dimensions
    if dim is None:
        raise RuntimeError("no box info available from topology/trajectory")

    # ---- select residues (whole molecules) ----
    resind = pick_residues(u, args.per_resname)
    mask = np.isin(u.atoms.resindices, resind)
    atoms = u.atoms[mask]
    old2new = -np.ones(u.atoms.n_atoms, dtype=int)
    old2new[atoms.indices] = np.arange(atoms.n_atoms)

    # ---- element / type ----
    mass_table = build_element_mass()
    elements = []
    for m in atoms.masses:
        e = mass_to_element(m, mass_table)
        if e is None or e.startswith("AMBIGUOUS"):
            raise RuntimeError(f"cannot uniquely map mass {m:.5f} to an element ({e})")
        elements.append(e)
    uniq_el = sorted(set(elements), key=lambda e: chemical_symbols.index(e))
    el2type = {e: i + 1 for i, e in enumerate(uniq_el)}
    type2el = {v: k for k, v in el2type.items()}
    atom_types = np.array([el2type[e] for e in elements])
    atom_els = np.array(elements)

    # ---- molecule ids: residue index+1 (renumbered) ----
    resid_uniq = sorted(set(atoms.residues.resindices))
    res2mol = {r: i + 1 for i, r in enumerate(resid_uniq)}
    atom_mol = np.array([res2mol[r] for r in atoms.resindices])

    # ---- bonds (intra-subset) ----
    bonds = []
    for b in u.bonds:
        i, j = b.indices
        if old2new[i] >= 0 and old2new[j] >= 0:
            bonds.append((old2new[i] + 1, old2new[j] + 1))
    bonds = sorted(set(tuple(sorted(x)) for x in bonds))

    # ---- write data file (first frame coords) ----
    xlo, ylo, zlo = None, None, None
    xr, yr, zr, tilt, ortho = lammps_box_from_dimensions(dim)
    pos = atoms.positions  # angstrom
    q = atoms.charges if hasattr(atoms, "charges") else np.zeros(atoms.n_atoms)

    with open(out + ".data", "w") as fh:
        fh.write("LAMMPS data file converted from GROMACS by lat convert (frame 0)\n\n")
        fh.write(f"{atoms.n_atoms} atoms\n")
        fh.write(f"{len(bonds)} bonds\n\n")
        fh.write(f"{len(uniq_el)} atom types\n")
        fh.write(f"{1 if bonds else 0} bond types\n\n")
        fh.write(f"{xr[0]:.6f} {xr[1]:.6f} xlo xhi\n")
        fh.write(f"{yr[0]:.6f} {yr[1]:.6f} ylo yhi\n")
        fh.write(f"{zr[0]:.6f} {zr[1]:.6f} zlo zhi\n")
        if tilt is not None:
            fh.write(f"{tilt[0]:.6f} {tilt[1]:.6f} {tilt[2]:.6f} xy xz yz\n")
        fh.write("\nMasses\n\n")
        for e in uniq_el:
            fh.write(f"{el2type[e]} {mass_table[e]:.6f}  # {e}\n")
        fh.write("\nAtoms # full\n\n")
        for k in range(atoms.n_atoms):
            fh.write(f"{k+1} {atom_mol[k]} {atom_types[k]} {q[k]:.6f} "
                     f"{pos[k,0]:.6f} {pos[k,1]:.6f} {pos[k,2]:.6f}\n")
        if bonds:
            fh.write("\nBonds\n\n")
            for bid, (i, j) in enumerate(bonds, 1):
                fh.write(f"{bid} 1 {i} {j}\n")

    # ---- write dump (lammpstrj) ----
    if args.stride is None:
        args.stride = max(1, nframes // max(1, args.frames))
    frame_ids = [min(nframes - 1, args.start + k * args.stride) for k in range(args.frames)]
    with open(out + ".lammpstrj", "w") as fh:
        for fi in frame_ids:
            u.trajectory[fi]
            box = u.dimensions
            xr, yr, zr, tilt, ortho = lammps_box_from_dimensions(box)
            fh.write("ITEM: TIMESTEP\n")
            fh.write(f"{fi}\n")
            fh.write("ITEM: NUMBER OF ATOMS\n")
            fh.write(f"{atoms.n_atoms}\n")
            fh.write("ITEM: BOX BOUNDS pp pp pp\n")
            fh.write(f"{xr[0]:.6f} {xr[1]:.6f}\n")
            fh.write(f"{yr[0]:.6f} {yr[1]:.6f}\n")
            fh.write(f"{zr[0]:.6f} {zr[1]:.6f}\n")
            if tilt is not None:
                fh.write(f"{tilt[0]:.6f} {tilt[1]:.6f} {tilt[2]:.6f}\n")
            fh.write("ITEM: ATOMS id mol type q x y z\n")
            p = atoms.positions
            for k in range(atoms.n_atoms):
                fh.write(f"{k+1} {atom_mol[k]} {atom_types[k]} {q[k]:.6f} "
                         f"{p[k,0]:.6f} {p[k,1]:.6f} {p[k,2]:.6f}\n")

    # ---- reduced gmx sample (.gro + .xtc) ----
    if not args.no_gmx_sample:
        u.trajectory[0]
        atoms.write(out + ".gro")
        with mda.Writer(out + ".xtc", atoms.n_atoms) as W:
            for fi in frame_ids:
                u.trajectory[fi]
                W.write(atoms)

    # ---- mapping csv ----
    with open(out + ".mapping.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "mol", "type", "element", "mass", "charge", "resname", "resid", "atomname"])
        for k in range(atoms.n_atoms):
            a = atoms[k]
            w.writerow([k + 1, atom_mol[k], atom_types[k], atom_els[k],
                        f"{mass_table[atom_els[k]]:.6f}", f"{q[k]:.6f}",
                        a.resname, a.resid, a.name])

    # ---- selections ----
    counts = {}
    for e in uniq_el:
        counts[e] = int(np.sum(atom_els == e))
    resname_counts = {}
    for rn in sorted(set(atoms.residues.resnames)):
        resname_counts[rn] = int(np.sum(atoms.residues.resnames == rn))
    sel = OrderedDict()
    sel["_meta"] = OrderedDict(
        source_top=args.top, source_traj=args.traj, frames=frame_ids,
        natoms=int(atoms.n_atoms), nmol=int(len(resid_uniq)),
        ntypes=len(uniq_el), type2element={str(k): v for k, v in type2el.items()},
        element_mass_source="ase.data / IUPAC standard atomic weights (tolerance 0.2%)",
        vdw_source="reax_tools ELEMENT_ATOMIC_RADII (Bondi)",
    )
    sel["selections"] = []
    for e in uniq_el:
        sel["selections"].append(OrderedDict(
            name=e, element=e, type=int(el2type[e]),
            mda=f"element {e}" if False else f"resname {e.upper()}" if e.upper() in resname_counts else f"name {e.upper()}",
            lmp=f"type {el2type[e]}",
            count=counts[e],
        ))
    with open(out + ".selections.yaml", "w") as fh:
        # minimal YAML (no pyyaml dependency)
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
        dump(sel)

    print(f"[ok] wrote {out}.data, {out}.lammpstrj")
    print(f"     natoms={atoms.n_atoms} nmol={len(resid_uniq)} ntypes={len(uniq_el)} frames={frame_ids}")
    print(f"     type->element: {type2el}")
    print(f"     element counts: {counts}")
    print(f"     resname counts: {resname_counts}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""LAMMPS data + dump  ->  GROMACS topology + trajectory (for gmx ground truth).

Why this exists: `tests/compat/gmx_regress.py` can only compare lat against gmx
on GROMACS-native inputs (tpr/xtc).  Nothing turned a LAMMPS system into
something `/usr/bin/gmx` could read, so LAMMPS-only systems (and therefore the
LAMMPS data-file parser, its atom-order handling and its bond list) had no
native-gmx ground truth at all.

Key conventions
---------------
* GROMACS atom order = **id order** (id 1..N).  The LAMMPS Bonds section
  references atom ids, so with ids exactly 1..N the bond list maps straight to
  GROMACS 1-based indices and no bond bookkeeping is needed.
* Each dump frame is sorted by `id` (LAMMPS reorders atoms spatially).
* One `[ moleculetype ]` per LAMMPS molecule, so `--pbc whole` and
  `-nmol`-style molecule splitting behave; `[ molecules ]` lists one of each.
* Charges are zero and only masses / bonds / coordinates are meaningful: this
  topology is for analysis (gyrate, rms, rmsf, sasa, rdf, hbond), not dynamics.
* vdW sigma is derived from lat's own vdw table via r = 0.5*sigma*2^(1/6),
  which is the radius `gmx sasa` uses, so SASA radii match lat's Bondi-like set.
* Trajectory is written as a multi-frame .gro in a temp file (a .gro is nm by
  definition, so there is no unit ambiguity) and immediately converted to .xtc
  by gmx trjconv.

Usage
-----
  lammps_to_gmx.py --data sys.data --traj sys.lammpstrj --out gmxdir
  lammps_to_gmx.py --data sys.data --out gmxdir --frames 0:101:1
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

GMX = os.environ.get("GMX", shutil.which("gmx") or "gmx")
_TWO_SIXTH = 2.0 ** (1.0 / 6.0)

# Independent mass -> element table (only what this converter needs to name
# atoms and to pick vdW radii).  Deliberately NOT lat's inference, so the
# reference stays independent of the code under test.
_ELEMENT_MASS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "P": 30.974, "S": 32.06, "Cl": 35.45, "Na": 22.99, "Li": 6.94,
}

_MDP = """\
integrator      = md
nsteps          = 0
dt              = 0.001
continuation    = yes
cutoff-scheme   = Verlet
coulombtype     = cutoff
rcoulomb        = {rcut}
vdwtype         = cutoff
rvdw            = {rcut}
nstlist         = 1
verlet-buffer-tolerance = -1
pbc             = xyz
constraints     = none
"""


class GmxError(RuntimeError):
    pass


def gmx(cmd, stdin=None, cwd=None):
    env = dict(os.environ, GMX_MAXBACKUP="-1")
    p = subprocess.run([GMX, *cmd], input=stdin, capture_output=True, text=True,
                       env=env, cwd=cwd)
    if p.returncode != 0:
        raise GmxError(f"gmx {' '.join(str(c) for c in cmd[:3])} failed:\n"
                       f"{p.stderr[-1500:]}")
    return p


def run_gmx(cmd, stdin=None, cwd=None):
    """Like gmx() but tolerates failure (returns (ok, stdout+stderr))."""
    try:
        p = gmx(cmd, stdin=stdin, cwd=cwd)
        return True, p.stdout + p.stderr
    except GmxError as exc:
        return False, str(exc)


# --------------------------------------------------------------------- parse
# column layout per LAMMPS atom_style, taken from the `Atoms # <style>` comment
_ATOM_LAYOUT = {
    "atomic":    [(0, "id"), (1, "type"), (2, "x"), (3, "y"), (4, "z")],
    "charge":    [(0, "id"), (1, "type"), (2, "q"), (3, "x"), (4, "y"), (5, "z")],
    "molecular": [(0, "id"), (1, "mol"), (2, "type"), (3, "x"), (4, "y"), (5, "z")],
    "bond":      [(0, "id"), (1, "mol"), (2, "type"), (3, "x"), (4, "y"), (5, "z")],
    "full":      [(0, "id"), (1, "mol"), (2, "type"), (3, "q"), (4, "x"), (5, "y"), (6, "z")],
}

# Sections this converter ignores but must not mistake for masses/atoms/bonds.
_OTHER_SECTIONS = ("pair", "pairij", "bond", "angle", "dihedral", "improper",
                   "atom", "velocities", "ellipsoids", "lines", "triangles",
                   "bodies", "angles", "dihedrals", "impropers")


def read_lammps_data(path):
    sec = None
    style = None
    atoms, bond_rows, masses = [], [], {}
    lo, hi, tilt = [None] * 3, [None] * 3, None
    for line in open(path):
        parts = line.split("#")
        s = parts[0].strip()
        comment = parts[1].strip() if len(parts) > 1 else ""
        if not s:
            continue
        head = s.split()[0]
        first_numeric = head[:1].isdigit() or (head[:1] in "+-." and head[1:2].isdigit())
        t = s.split()
        # box bounds can appear before or between sections; match on the pattern
        if len(t) >= 4 and t[2:4] == ["xlo", "xhi"]:
            lo[0], hi[0] = float(t[0]), float(t[1])
            continue
        if len(t) >= 4 and t[2:4] == ["ylo", "yhi"]:
            lo[1], hi[1] = float(t[0]), float(t[1])
            continue
        if len(t) >= 4 and t[2:4] == ["zlo", "zhi"]:
            lo[2], hi[2] = float(t[0]), float(t[1])
            continue
        if len(t) >= 6 and t[3:6] == ["xy", "xz", "yz"]:
            tilt = [float(t[0]), float(t[1]), float(t[2])]
            continue
        if not first_numeric:
            low = head.lower()
            if low in ("masses", "atoms", "bonds"):
                sec = low
                if low == "atoms" and comment:
                    style = comment.split()[0].lower()
            else:
                sec = "other"          # e.g. `Pair Coeffs`, `Velocities`, a title line
            continue
        if sec == "masses":
            masses[int(t[0])] = float(t[1])
        elif sec == "atoms":
            style_l = style or ("full" if len(t) >= 7 else "atomic")
            layout = _ATOM_LAYOUT.get(style_l)
            if layout is None:
                raise SystemExit(f"{path}: unsupported Atoms style {style!r}")
            vals = {name: t[i] for i, name in layout}
            atoms.append((int(vals["id"]), int(vals.get("mol", 1)), int(vals["type"]),
                          float(vals["x"]), float(vals["y"]), float(vals["z"])))
        elif sec == "bonds":
            if len(t) >= 4:
                bond_rows.append((int(t[2]), int(t[3])))
    if not atoms:
        raise SystemExit(f"{path}: no Atoms section")
    if any(v is None for v in lo + hi):
        raise SystemExit(f"{path}: incomplete box bounds")
    if tilt is not None and any(abs(v) > 1e-8 for v in tilt):
        raise SystemExit(f"{path}: triclinic boxes are not supported by this converter")
    return dict(atoms=atoms, bonds=np.array(bond_rows) - 1 if bond_rows else None,
                masses=masses, lo=np.array(lo, float), hi=np.array(hi, float))


def elements_of_types(masses, tol=0.002):
    out = {}
    for t, m in masses.items():
        hits = [e for e, mm in _ELEMENT_MASS.items() if abs(m - mm) / mm <= tol]
        out[t] = hits[0] if len(hits) == 1 else None
    return out


def vdw_radii(cfg_path=None):
    """lat's vdW table (Angstrom) if available, else {}."""
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))
        from lat.config import load as load_config
        cfg = load_config(cfg_path) if cfg_path else load_config()
        return {str(k): float(v) for k, v in (cfg.get("vdw", default={}) or {}).items()}
    except Exception:
        return {}


# --------------------------------------------------------------------- write
def _frames(path):
    """Yield (ids, positions, box) per frame from a LAMMPS dump, in file order.

    `box` is (lo, lengths, tilt). The cell comes from the trajectory itself: a
    data file written early in a run can carry a different box (NPT), and
    mixing the two silently rescales the system.
    """
    with open(path) as fh:
        n = None
        box = None
        while True:
            line = fh.readline()
            if not line:
                return
            if line.startswith("ITEM: TIMESTEP"):
                fh.readline()
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                n = int(fh.readline())
            elif line.startswith("ITEM: BOX BOUNDS"):
                rows = [[float(v) for v in fh.readline().split()] for _ in range(3)]
                lo = np.array([r[0] for r in rows])
                hi = np.array([r[1] for r in rows])
                tilt = (np.array([rows[0][2], rows[0][3], rows[0][4]])
                        if len(rows[0]) >= 5 else np.zeros(3))
                if np.any(np.abs(tilt) > 1e-8):
                    raise SystemExit(f"{path}: triclinic trajectory boxes are not "
                                     f"supported by this converter")
                box = (lo, hi - lo, tilt)
            elif line.startswith("ITEM: ATOMS"):
                cols = line.split()[2:]
                ci = cols.index("id")
                xyz = [cols.index(c) for c in ("x", "y", "z")]
                raw = [fh.readline().split() for _ in range(n)]
                yield (np.array([r[ci] for r in raw], dtype=int),
                       np.array([[float(r[k]) for k in xyz] for r in raw]),
                       box)


def convert(data, traj, outdir, cfg_path=None, stride=1, frames=None):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sysd = read_lammps_data(data)

    ids = np.array([a[0] for a in sysd["atoms"]])
    order = np.argsort(ids, kind="stable")
    if not np.array_equal(ids[order], np.arange(1, len(ids) + 1)):
        raise SystemExit(f"{data}: this converter expects atom ids to be exactly 1..N")
    mols = np.array([a[1] for a in sysd["atoms"]])[order]
    types = np.array([a[2] for a in sysd["atoms"]])[order]
    xyz = np.array([[a[3], a[4], a[5]] for a in sysd["atoms"]])[order]
    bonds = sysd["bonds"]
    lo, L = sysd["lo"], sysd["hi"] - sysd["lo"]
    n = len(ids)

    elem = elements_of_types(sysd["masses"])
    unknown = [t for t, e in elem.items() if e is None]
    if unknown:
        raise SystemExit(f"{data}: cannot map type(s) {unknown} to an element "
                         f"(masses {[sysd['masses'][t] for t in unknown]})")
    print("[lammps_to_gmx] type->element:",
          {int(t): e for t, e in sorted(elem.items())})

    rad = vdw_radii(cfg_path)
    missing = sorted({elem[int(t)] for t in types} - set(rad))
    if missing:
        raise SystemExit(f"no vdW radius for {missing}; pass --cfg or check lat config")

    # Atom names are the bare element symbol: the .gro name field is five
    # characters, and appending the atom id overflows it past 9999 atoms,
    # which corrupts the fixed-format line. Analysis selects by index group.
    names = [elem[int(types[i])] for i in range(n)]

    def atom_lines(pos_nm):
        return ["%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (
            int(mols[i]), "MOL", names[i], i + 1,
            pos_nm[i, 0], pos_nm[i, 1], pos_nm[i, 2]) for i in range(n)]

    def box_line(box_len):
        return "%10.5f%10.5f%10.5f\n" % tuple(np.asarray(box_len, float) / 10.0)

    # The single-frame sys.gro (which grompp turns into the tpr, and therefore
    # the reference cell for every analysis) must use the trajectory's own
    # cell: a data file written before an NPT run carries a different box, and
    # mixing the two silently rescales the system.
    sys_lo, sys_L, sys_xyz = lo, L, xyz
    if traj:
        for fid, p, box in _frames(traj):
            if box is not None:
                sys_lo, sys_L = box[0], box[1]
            sys_xyz = p[np.argsort(fid)]
            break
    with open(outdir / "sys.gro", "w") as fh:
        fh.write("LAMMPS -> GROMACS\n%d\n" % n)
        fh.writelines(atom_lines((sys_xyz - sys_lo) / 10.0))
        fh.write(box_line(sys_L))
    L = sys_L

    blocks = {int(m): np.where(mols == m)[0] for m in np.unique(mols)}
    for m, rows in blocks.items():
        if not np.array_equal(rows, np.arange(rows[0], rows[0] + len(rows))):
            raise SystemExit(f"molecule {m} is not contiguous in id order")
    with open(outdir / "sys.top", "w") as fh:
        fh.write("[ defaults ]\n1 1 no 1.0 1.0\n\n[ atomtypes ]\n")
        for t in sorted(sysd["masses"]):
            e = elem[t]
            sigma = 2.0 * rad[e] / 10.0 / _TWO_SIXTH
            fh.write(f"T{t}  {sysd['masses'][t]:.5f}  0.0000  A  {sigma:.5f}  1.0e-01\n")
        for m, rows in blocks.items():
            off = int(rows[0])
            fh.write(f"\n[ moleculetype ]\nMOL{m}  3\n\n[ atoms ]\n")
            for i in rows:
                t = int(types[i])
                fh.write("%5d  T%-3d %5d  MOL  %-4s %5d  %8.5f  %8.5f\n" % (
                    i - off + 1, t, 1, names[i], i - off + 1, 0.0, sysd["masses"][t]))
            fh.write("\n[ bonds ]\n")
            for a, b in (bonds if bonds is not None else []):
                if off <= a < off + len(rows):
                    if not (off <= b < off + len(rows)):
                        raise SystemExit(f"bond {a}-{b} crosses molecule blocks")
                    fh.write("%5d %5d 1  0.150000  300000.0\n" % (a - off + 1, b - off + 1))
        fh.write("\n[ system ]\nLAMMPS system\n\n[ molecules ]\n")
        for m in blocks:
            fh.write(f"MOL{m}  1\n")

    with open(outdir / "sys.ndx", "w") as fh:
        def grp(name, idx):
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(idx), 15):
                fh.write(" ".join(str(int(k) + 1) for k in idx[i:i + 15]) + "\n")
            fh.write("\n")
        grp("System", np.arange(n))
        for m, rows in blocks.items():
            grp(f"Mol{m}", rows)
        # element groups: gmx rdf selects by group name far more reliably than
        # by a wildcard selection string
        for e in sorted({elem[int(t)] for t in types}):
            grp(f"El_{e}", np.where(np.array([elem[int(t)] for t in types]) == e)[0])

    mdp = outdir / "grompp.mdp"
    mdp.write_text(_MDP.format(rcut=min(1.0, float(L.min()) / 10.0 / 2.0)))
    gmx(["grompp", "-f", str(mdp), "-c", str(outdir / "sys.gro"),
         "-p", str(outdir / "sys.top"), "-o", str(outdir / "sys.tpr"),
         "-maxwarn", "20"])

    if traj:
        with tempfile.TemporaryDirectory() as tmp:
            gro = Path(tmp) / "traj.gro"
            nfr = 0
            with open(gro, "w") as fh:
                for k, (fid, p, fbox) in enumerate(_frames(traj)):
                    if k % stride:
                        continue
                    if frames and not (frames[0] <= k < frames[1]):
                        continue
                    flo, flen = (fbox[0], fbox[1]) if fbox is not None else (lo, L)
                    pos = (p[np.argsort(fid)] - flo) / 10.0
                    fh.write("frame %d\n%d\n" % (k, n))
                    fh.writelines(atom_lines(pos))
                    fh.write(box_line(flen))
                    nfr += 1
            if nfr == 0:
                raise SystemExit(f"{traj}: no frames selected")
            gmx(["trjconv", "-s", str(outdir / "sys.tpr"), "-f", str(gro),
                 "-n", str(outdir / "sys.ndx"), "-o", str(outdir / "traj.xtc")],
                stdin="0\n")
            print(f"[lammps_to_gmx] traj.xtc: {nfr} frames")
        # `-pbc whole` fails on a data file whose molecule decomposition puts a
        # box-spanning object into one moleculetype (gmx: "inconsistent shifts
        # over periodic boundaries"). Analyses that need whole molecules
        # (gyrate/rms/rmsf) are impossible on such a system anyway; distance
        # and hbond analyses only need minimum-image coordinates, so fall back
        # to the wrapped trajectory instead of failing the whole conversion.
        try:
            gmx(["trjconv", "-s", str(outdir / "sys.tpr"), "-f", str(outdir / "traj.xtc"),
                 "-n", str(outdir / "sys.ndx"), "-pbc", "whole",
                 "-o", str(outdir / "whole.xtc")], stdin="0\n")
            xtc = outdir / "whole.xtc"
        except GmxError as exc:
            print("[lammps_to_gmx] WARNING: gmx trjconv -pbc whole failed "
                  f"({str(exc).strip().splitlines()[-1][:90]}); continuing with the "
                  "wrapped trajectory -- gyrate/rms/rmsf are NOT valid here")
            wh = outdir / "whole.xtc"
            if wh.exists():
                wh.unlink()
            xtc = outdir / "traj.xtc"
        gmx(["trjconv", "-s", str(outdir / "sys.tpr"), "-f", str(xtc),
             "-n", str(outdir / "sys.ndx"), "-dump", "0", "-o", str(outdir / "ref.gro")],
            stdin="0\n")
        gmx(["grompp", "-f", str(mdp), "-c", str(outdir / "ref.gro"),
             "-p", str(outdir / "sys.top"), "-o", str(outdir / "ref.tpr"),
             "-maxwarn", "20"])
    print(f"[lammps_to_gmx] wrote {outdir}/sys.gro sys.top sys.ndx sys.tpr"
          + (f" rec.tpr {xtc.name}" if traj else ""))
    return outdir


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="LAMMPS data file")
    ap.add_argument("--traj", help="LAMMPS dump (optional)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--cfg", help="LAT_CONFIG to take vdW radii from")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--frames", help="start:stop (exclusive) frame window")
    a = ap.parse_args(argv)
    win = None
    if a.frames:
        p = a.frames.split(":")
        win = (int(p[0]), int(p[1]))
    convert(a.data, a.traj, a.out, cfg_path=a.cfg, stride=a.stride, frames=win)
    return 0


if __name__ == "__main__":
    sys.exit(main())

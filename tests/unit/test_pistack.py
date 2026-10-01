"""Unit test for pistack: two parallel rings stack, a perpendicular one does not."""

import subprocess
import sys

import numpy as np

R = 1.4
HEX = np.array([[R * np.cos(t), R * np.sin(t), 0.0]
                for t in np.linspace(0, 2 * np.pi, 6, endpoint=False)])


def build(path, z_shift=3.5, perpendicular=False):
    ring_b = HEX.copy()
    if perpendicular:
        ring_b = np.column_stack([HEX[:, 0], np.zeros(6), HEX[:, 1]])
    ring_b = ring_b + np.array([0, 0, z_shift])
    coords = np.vstack([HEX, ring_b]) + 10.0
    mols = [1] * 6 + [2] * 6
    lines = ["LAMMPS data file", "", "12 atoms", "1 atom types", "",
             "0.0 20.0 xlo xhi", "0.0 20.0 ylo yhi", "0.0 20.0 zlo zhi", "",
             "Masses", "", "1 12.011", "", "Atoms # full", ""]
    for i in range(12):
        lines.append(f"{i+1} {mols[i]} 1 0.0 {coords[i,0]} {coords[i,1]} {coords[i,2]}")
    path.write_text("\n".join(lines) + "\n")
    return coords


def write_dump(path, coords):
    out = ["ITEM: TIMESTEP", "0", "ITEM: NUMBER OF ATOMS", "12",
           "ITEM: BOX BOUNDS pp pp pp", "0.0 20.0", "0.0 20.0", "0.0 20.0",
           "ITEM: ATOMS id mol type q x y z"]
    mols = [1] * 6 + [2] * 6
    for i in range(12):
        out.append(f"{i+1} {mols[i]} 1 0.0 {coords[i,0]} {coords[i,1]} {coords[i,2]}")
    path.write_text("\n".join(out) + "\n")


def run_lat(repo_root, args):
    return subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)


def test_parallel_vs_perpendicular(tmp_path, repo_root):
    for perpendicular, expected in ((False, 1), (True, 0)):
        d = tmp_path / f"p{int(perpendicular)}"
        d.mkdir()
        coords = build(d / "p.data", perpendicular=perpendicular)
        write_dump(d / "p.lammpstrj", coords)
        proc = run_lat(repo_root, ["pistack", "-s", str(d / "p.data"), "-f", str(d / "p.lammpstrj"),
                                   "-sel1", "type 1", "-distance", "5.5", "-angle", "30",
                                   "-o", str(d / "out")])
        assert proc.returncode == 0, proc.stderr
        table = np.atleast_2d(np.genfromtxt(d / "out.pistack.out.csv", delimiter=",", comments="#"))
        assert int(table[0, 1]) == expected

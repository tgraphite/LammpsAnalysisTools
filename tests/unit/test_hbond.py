"""Unit test for hbond on a hand-built O-H...O dimer (expect exactly 1 bond)."""

import subprocess
import sys

import numpy as np
import pytest

DATA = """LAMMPS data file

6 atoms
4 bonds

2 atom types
1 bond types

0.0 20.0 xlo xhi
0.0 20.0 ylo yhi
0.0 20.0 zlo zhi

Masses

1 15.999
2 1.008

Atoms # full

1 1 1 0.0 0.0 0.0 0.0
2 1 2 0.0 0.96 0.0 0.0
3 1 2 0.0 -0.24 0.93 0.0
4 2 1 0.0 2.9 0.0 0.0
5 2 2 0.0 3.86 0.0 0.0
6 2 2 0.0 2.66 0.93 0.0

Bonds

1 1 1 2
2 1 1 3
3 1 4 5
4 1 4 6
"""


def write_dump(path, n_frames=3):
    out = []
    for k in range(n_frames):
        out += ["ITEM: TIMESTEP", str(k), "ITEM: NUMBER OF ATOMS", "6",
                "ITEM: BOX BOUNDS pp pp pp", "0.0 20.0", "0.0 20.0", "0.0 20.0",
                "ITEM: ATOMS id type x y z",
                "1 1 0.0 0.0 0.0", "2 2 0.0 0.96 0.0", "3 2 0.0 -0.24 0.93",
                "4 1 0.0 2.9 0.0", "5 2 0.0 3.86 0.0", "6 2 0.0 2.66 0.93"]
    path.write_text("\n".join(out) + "\n")


def run_lat(repo_root, args):
    return subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)


@pytest.fixture
def setup(tmp_path):
    data = tmp_path / "h.data"
    dump = tmp_path / "h.lammpstrj"
    data.write_text(DATA)
    write_dump(dump)
    return tmp_path, data, dump


def test_one_hbond_per_frame(repo_root, setup):
    tmp, data, dump = setup
    proc = run_lat(repo_root, ["hbond", "-s", str(data), "-f", str(dump),
                               "-o", str(tmp / "h")])
    assert proc.returncode == 0, proc.stderr
    out = np.genfromtxt(tmp / "h.hbond.out.csv", delimiter=",", comments="#")
    assert list(out[:, 1]) == [1.0, 1.0, 1.0]
    acf = np.genfromtxt(tmp / "h.hbond.acf.out.csv", delimiter=",", comments="#")
    assert acf[0, 1] == pytest.approx(1.0)


def test_overlapping_groups_error(repo_root, setup):
    tmp, data, dump = setup
    proc = run_lat(repo_root, ["hbond", "-s", str(data), "-f", str(dump),
                               "-sel1", "type 1", "-sel2", "type 1 2",
                               "-o", str(tmp / "h2")])
    assert proc.returncode == 1
    assert "identical or completely non-overlapping" in proc.stderr


def test_hbond_midpoint_is_actual_midpoint():
    from lat.core.pbc import Box
    from lat.modules.hbond import _frame_hbonds

    box = Box.orthorhombic([10.0, 10.0, 10.0])
    pos = np.array([[1.0, 0.0, 0.0], [1.5, 0.0, 0.0], [2.0, 0.0, 0.0]])
    triples, mids = _frame_hbonds({0: [1]}, np.array([2]), pos, box, 3.5, "DA",
                                  30.0, same=False, bonds_set=set())
    assert triples == {(0, 1, 2)}
    assert np.allclose(mids[0], [1.5, 0.0, 0.0])


def test_hbond_midpoint_wraps_across_boundary():
    from lat.core.pbc import Box
    from lat.modules.hbond import _frame_hbonds

    box = Box.orthorhombic([10.0, 10.0, 10.0])
    # D at 0.5, H at 0.1 (points -x), A at 9.5: D's nearest image is at 10.5,
    # so the true midpoint (10.0) wraps to 0.0.
    pos = np.array([[0.5, 0.0, 0.0], [0.1, 0.0, 0.0], [9.5, 0.0, 0.0]])
    triples, mids = _frame_hbonds({0: [1]}, np.array([2]), pos, box, 3.5, "DA",
                                  30.0, same=False, bonds_set=set())
    assert triples == {(0, 1, 2)}
    assert np.allclose(mids[0], [0.0, 0.0, 0.0], atol=1e-9)

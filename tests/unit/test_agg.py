"""Unit test for agg / cluster on a hand-built 1D ion arrangement."""

import subprocess
import sys

import numpy as np
import pytest

# A1-A2 bridged by B1 (AGG); A3-B2 (CIP); A4 isolated (SSIP)
COORDS = {
    1: (0.0, 0.0, 0.0),
    2: (3.5, 0.0, 0.0),
    3: (1.75, 0.0, 0.0),
    4: (20.0, 0.0, 0.0),
    5: (23.0, 0.0, 0.0),
    6: (35.0, 0.0, 0.0),
}
TYPES = {1: 1, 2: 1, 3: 2, 4: 1, 5: 2, 6: 1}
CHARGES = {1: 1.0, 2: 1.0, 3: -1.0, 4: 1.0, 5: -1.0, 6: 1.0}


def write_data(path):
    lines = ["LAMMPS data file", "", "6 atoms", "2 atom types", "",
             "0.0 40.0 xlo xhi", "0.0 10.0 ylo yhi", "0.0 10.0 zlo zhi", "",
             "Masses", "", "1 6.94", "2 18.998", "", "Atoms # full", ""]
    for i in range(1, 7):
        x, y, z = COORDS[i]
        lines.append(f"{i} {i} {TYPES[i]} {CHARGES[i]} {x} {y} {z}")
    path.write_text("\n".join(lines) + "\n")


def write_dump(path, n_frames=2):
    out = []
    for k in range(n_frames):
        out += ["ITEM: TIMESTEP", str(k), "ITEM: NUMBER OF ATOMS", "6",
                "ITEM: BOX BOUNDS pp pp pp", "0.0 40.0", "0.0 10.0", "0.0 10.0",
                "ITEM: ATOMS id mol type q x y z"]
        for i in range(1, 7):
            x, y, z = COORDS[i]
            out.append(f"{i} {i} {TYPES[i]} {CHARGES[i]} {x} {y} {z}")
    path.write_text("\n".join(out) + "\n")


def run_lat(repo_root, args):
    return subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)


@pytest.fixture
def setup(tmp_path):
    data = tmp_path / "a.data"
    dump = tmp_path / "a.lammpstrj"
    write_data(data)
    write_dump(dump)
    return tmp_path, data, dump


def test_agg_speciation(repo_root, setup):
    tmp, data, dump = setup
    proc = run_lat(repo_root, ["agg", "-s", str(data), "-f", str(dump),
                               "-sel1", "type 1", "-sel2", "type 2",
                               "-cutoff", "5", "-o", str(tmp / "a")])
    assert proc.returncode == 0, proc.stderr
    out = np.genfromtxt(tmp / "a.agg.out.csv", delimiter=",", comments="#")
    header = [l for l in open(tmp / "a.agg.out.csv") if l.startswith("# columns")][0]
    cols = [c.strip() for c in header.split(":")[1].split(",")]
    row = dict(zip(cols, out[0]))
    assert row["A2B1"] == 1 and row["A1B1"] == 1 and row["A1B0"] == 1

    lab = np.genfromtxt(tmp / "a.agg.labels.out.csv", delimiter=",", comments="#")
    # frame, A_SSIP, A_CIP, A_AGG, B_SSIP, B_CIP, B_AGG
    assert lab[0, 1:4].tolist() == [0.25, 0.25, 0.5]
    assert lab[0, 4:7].tolist() == [0.0, 0.5, 0.5]


def test_cluster_distributions(repo_root, setup):
    tmp, data, dump = setup
    proc = run_lat(repo_root, ["cluster", "-s", str(data), "-f", str(dump),
                               "-sel1", "type 1", "-sel2", "type 2",
                               "-cutoff", "5", "-o", str(tmp / "c")])
    assert proc.returncode == 0, proc.stderr
    cn = np.genfromtxt(tmp / "c.cluster.cn.out.csv", delimiter=",", comments="#")
    p = {int(k): round(float(v), 6) for k, v in cn[:, [0, 2]]}
    assert p == {0: 0.25, 1: 0.75}
    size = np.genfromtxt(tmp / "c.cluster.size.out.csv", delimiter=",", comments="#")
    ps = {int(k): round(float(v), 6) for k, v in size[:, [0, 2]]}
    assert ps == {1: pytest.approx(1 / 3), 2: pytest.approx(1 / 3), 3: pytest.approx(1 / 3)}


def test_cutoff_required(repo_root, setup):
    tmp, data, dump = setup
    proc = run_lat(repo_root, ["agg", "-s", str(data), "-f", str(dump),
                               "-sel1", "type 1", "-sel2", "type 2", "-o", str(tmp / "x")])
    assert proc.returncode == 1
    assert "cutoff" in proc.stderr.lower()

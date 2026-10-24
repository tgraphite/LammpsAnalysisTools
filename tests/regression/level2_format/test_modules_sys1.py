"""End-to-end smoke: every module runs on the real reduced sample.

This is the integration layer that unit tests on synthetic systems miss
(e.g. frames with zero hydrogen bonds).
"""

import subprocess
import sys

import numpy as np
import pytest


def lat(repo_root, args):
    return subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)


CASES = [
    ("select", ["select", "-sel", "type 2"]),
    ("describe", ["describe", "-module", "rdf"]),
    ("rdf", ["rdf", "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "10"]),
    ("density1d", ["density", "-sel", "type 2", "-axis", "z"]),
    ("density2d", ["density", "-sel", "type 2", "-dim", "2", "-nbins", "20"]),
    ("gyrate", ["gyrate", "-sel", "mol 401-600"]),
    ("contact", ["contact", "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "3"]),
    ("rms", ["rms", "-sel", "type 3"]),
    ("msd", ["msd", "-sel", "type 2"]),
    # hbond is deliberately absent: sys1_reduced is a Li/FSI/FEC electrolyte
    # whose only H atoms sit on carbon, so it has no polar-H donor.  gmx hbond
    # refuses such a group with "No Donors found" and lat mirrors that on
    # purpose -- see test_hbond_requires_polar_donor.
    ("agg", ["agg", "-sel1", "type 2", "-sel2", "type 4 5 6 7", "-cutoff", "2.5"]),
    ("cluster", ["cluster", "-sel1", "type 2", "-sel2", "type 4 5 6 7",
                 "-cutoff", "2.5", "-focus"]),
    ("pistack", ["pistack", "-sel1", "type 3"]),
    ("freevol_ffv", ["freevol", "-mode", "ffv", "-grid_dr", "3.0"]),
    ("sasa", ["sasa", "-sel", "type 2"]),
    ("residence", ["residence", "-ref", "type 2", "-lig", "type 4 5 6 7", "-cutoff", "3"]),
    ("occupancy", ["occupancy", "-sel", "type 2", "-zone", "z:0-45,45-90,90-135"]),
]


@pytest.mark.parametrize("name,args", CASES, ids=[c[0] for c in CASES])
def test_module_runs_on_sample(name, args, repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    cmd = list(args)
    if name != "describe":
        cmd += ["-o", str(tmp_path / name), "-s", str(base) + ".data"]
        if name != "select":
            cmd += ["-f", str(base) + ".lammpstrj"]
    proc = lat(repo_root, cmd)
    assert proc.returncode == 0, f"{name} failed:\n{proc.stderr}"
    if name not in ("select", "describe"):
        assert list(tmp_path.glob(f"{name}.*.out.csv")), f"{name} produced no CSV"


def test_hbond_requires_polar_donor(repo_root, sample_dir, tmp_path):
    """A C-H is not a hydrogen-bond donor (gmx: "No Donors found")."""
    base = sample_dir / "sys1_reduced"
    proc = lat(repo_root, ["hbond", "-s", str(base) + ".data",
                           "-f", str(base) + ".lammpstrj", "-o", str(tmp_path / "h")])
    assert proc.returncode == 1
    assert "no donor" in proc.stderr.lower()


# A two-molecule O-H...O system whose acceptor is parked 8 A away: donors exist
# but no frame satisfies the criteria, so the geometry tables must still be
# written (this is the "frames with zero hydrogen bonds" case).
_OH_DATA = """LAMMPS data file

6 atoms
4 bonds

2 atom types

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
4 2 1 0.0 8.0 0.0 0.0
5 2 2 0.0 8.96 0.0 0.0
6 2 2 0.0 7.76 0.93 0.0

Bonds

1 1 1 2
2 1 1 3
3 1 4 5
4 1 4 6
"""


def test_hbond_geometry_with_empty_frames(repo_root, tmp_path):
    data = tmp_path / "oh.data"
    dump = tmp_path / "oh.lammpstrj"
    data.write_text(_OH_DATA)
    frames = []
    for k in range(3):
        frames += ["ITEM: TIMESTEP", str(k), "ITEM: NUMBER OF ATOMS", "6",
                   "ITEM: BOX BOUNDS pp pp pp", "0.0 20.0", "0.0 20.0", "0.0 20.0",
                   "ITEM: ATOMS id type x y z",
                   "1 1 0.0 0.0 0.0", "2 2 0.0 0.96 0.0", "3 2 0.0 -0.24 0.93",
                   "4 1 0.0 8.0 0.0", "5 2 0.0 8.96 0.0", "6 2 0.0 7.76 0.93"]
    dump.write_text("\n".join(frames) + "\n")

    proc = lat(repo_root, ["hbond", "-s", str(data), "-f", str(dump),
                           "-o", str(tmp_path / "h")])
    assert proc.returncode == 0, proc.stderr
    for tag in ("hbond", "hbond.acf", "hbond.dist_xz", "hbond.dist_yz", "hbond.dist_xy"):
        assert (tmp_path / f"h.{tag}.out.csv").is_file()
    counts = np.genfromtxt(tmp_path / "h.hbond.out.csv", delimiter=",", comments="#")
    assert list(counts[:, 1]) == [0.0, 0.0, 0.0]


def test_transport_guard_on_sample(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    proc = lat(repo_root, ["transport", "-s", str(base) + ".data",
                           "-f", str(base) + ".lammpstrj", "-sel", "type 2",
                           "-o", str(tmp_path / "t")])
    assert proc.returncode == 1
    assert "velocity" in proc.stderr.lower()

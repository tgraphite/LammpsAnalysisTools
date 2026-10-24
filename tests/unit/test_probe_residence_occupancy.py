"""Smoke/unit tests for freevol, residence, occupancy and transport."""

import subprocess
import sys

import numpy as np
import pytest


def run_lat(repo_root, args):
    return subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)


def write_full_data(path, atoms, box=10.0):
    """atoms: list of (mol, type, q, x, y, z)."""
    lines = ["LAMMPS data file", "", f"{len(atoms)} atoms", "2 atom types", "",
             f"0.0 {box} xlo xhi", f"0.0 {box} ylo yhi", f"0.0 {box} zlo zhi", "",
             "Masses", "", "1 12.011", "2 1.008", "", "Atoms # full", ""]
    for i, (m, t, q, x, y, z) in enumerate(atoms, 1):
        lines.append(f"{i} {m} {t} {q} {x} {y} {z}")
    path.write_text("\n".join(lines) + "\n")


def write_dump_full(path, frames, cols="x y z"):
    out = []
    for k, (atoms, boxes) in enumerate(frames):
        out += ["ITEM: TIMESTEP", str(k), "ITEM: NUMBER OF ATOMS", str(len(atoms)),
                "ITEM: BOX BOUNDS pp pp pp"]
        out += [f"0.0 {boxes[0]}", f"0.0 {boxes[1]}", f"0.0 {boxes[2]}"]
        out += [f"ITEM: ATOMS id mol type q {cols}"]
        for i, (m, t, q, x, y, z) in enumerate(atoms, 1):
            out.append(f"{i} {m} {t} {q} {x} {y} {z}")
    path.write_text("\n".join(out) + "\n")


def test_freevol_ffv(tmp_path, repo_root):
    atoms = [(1, 1, 0.0, 5.0, 5.0, 5.0)]
    data = tmp_path / "p.data"
    dump = tmp_path / "p.lammpstrj"
    write_full_data(data, atoms)
    write_dump_full(dump, [(atoms, (10, 10, 10))] * 2)
    proc = run_lat(repo_root, ["freevol", "-s", str(data), "-f", str(dump), "-mode", "ffv",
                               "-grid_dr", "0.5", "-radius", "0.3", "-o", str(tmp_path / "out")])
    assert proc.returncode == 0, proc.stderr
    header = open(tmp_path / "out.freevol.out.csv").read()
    import re
    m = re.search(r"ffv_mean=([\d.eE+-]+)", header)
    expected_occ = (4 / 3) * np.pi * (1.70 + 0.3) ** 3 / 1000
    assert float(m.group(1)) == pytest.approx(1 - expected_occ, abs=0.03)


def test_residence_always_bound(tmp_path, repo_root):
    atoms = [(1, 2, 0.0, 5.0, 5.0, 5.0), (2, 1, 0.0, 5.0, 7.0, 5.0)]
    data = tmp_path / "r.data"
    dump = tmp_path / "r.lammpstrj"
    write_full_data(data, atoms)
    write_dump_full(dump, [(atoms, (10, 10, 10))] * 4)
    proc = run_lat(repo_root, ["residence", "-s", str(data), "-f", str(dump),
                               "-ref", "id 1", "-lig", "id 2", "-cutoff", "3",
                               "-mode", "min", "-o", str(tmp_path / "r")])
    assert proc.returncode == 0, proc.stderr
    header = open(tmp_path / "r.residence.acf.out.csv").read()
    import re
    tau = float(re.search(r"tau_integral_frames=([\d.eE+-]+)", header).group(1))
    # permanently bound: unbiased ACF C(lag)=1 -> tau = sum_{lag>=1} = T-1 = 3
    assert tau == pytest.approx(3.0)


def test_occupancy_constant_zones(tmp_path, repo_root):
    atoms = [(1, 1, 0.0, 5.0, 5.0, 2.0), (2, 1, 0.0, 5.0, 5.0, 8.0)]
    data = tmp_path / "o.data"
    dump = tmp_path / "o.lammpstrj"
    write_full_data(data, atoms)
    write_dump_full(dump, [(atoms, (10, 10, 10))] * 3)
    proc = run_lat(repo_root, ["occupancy", "-s", str(data), "-f", str(dump),
                               "-sel", "type 1", "-zone", "z:0-5,5-10", "-o", str(tmp_path / "o")])
    assert proc.returncode == 0, proc.stderr
    out = np.genfromtxt(tmp_path / "o.occupancy.out.csv", delimiter=",", comments="#")
    assert list(out[:, 1]) == [1.0, 1.0, 1.0]
    assert list(out[:, 2]) == [1.0, 1.0, 1.0]


def test_transport_diffusion_with_velocity(tmp_path, repo_root):
    rng = np.random.default_rng(0)
    atoms = [(1, 1, 0.0, 5.0, 5.0, 5.0)]
    data = tmp_path / "tv.data"
    write_full_data(data, atoms)
    lines = []
    steps = np.cumsum(rng.normal(0, 0.2, size=(5, 3)) + rng.normal(0, 0.4, size=(5, 3)), axis=0)
    for k in range(5):
        x, y, z = 5 + steps[k]
        vx, vy, vz = rng.normal(0, 0.1, 3)
        lines += ["ITEM: TIMESTEP", str(k), "ITEM: NUMBER OF ATOMS", "1",
                  "ITEM: BOX BOUNDS pp pp pp", "0.0 10.0", "0.0 10.0", "0.0 10.0",
                  "ITEM: ATOMS id mol type q vx vy vz xu yu zu",
                  f"1 1 1 0.0 {vx} {vy} {vz} {x} {y} {z}"]
    (tmp_path / "tv.lammpstrj").write_text("\n".join(lines) + "\n")
    proc = run_lat(repo_root, ["transport", "-s", str(data), "-f", str(tmp_path / "tv.lammpstrj"),
                               "-sel", "all", "-o", str(tmp_path / "tv")])
    assert proc.returncode == 0, proc.stderr
    d = np.atleast_2d(np.genfromtxt(tmp_path / "tv.transport.d.out.csv", delimiter=",", comments="#"))
    assert d.shape[0] == 1 and np.isfinite(d[0, 2])


def test_transport_requires_velocity(tmp_path, repo_root):
    atoms = [(1, 1, 0.0, 5.0, 5.0, 5.0)]
    data = tmp_path / "t.data"
    dump = tmp_path / "t.lammpstrj"
    write_full_data(data, atoms)
    write_dump_full(dump, [(atoms, (10, 10, 10))] * 3)
    proc = run_lat(repo_root, ["transport", "-s", str(data), "-f", str(dump),
                               "-sel", "all", "-o", str(tmp_path / "t")])
    assert proc.returncode == 1
    assert "velocity" in proc.stderr.lower()

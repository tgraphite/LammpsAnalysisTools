"""Unit tests for density / msd / rms using synthetic LAMMPS inputs."""

import subprocess
import sys

import numpy as np
import pytest


def write_data(path, positions, box=20.0):
    n = len(positions)
    lines = [
        "LAMMPS data file", "", f"{n} atoms", "1 atom types", "",
        f"0.0 {box} xlo xhi", f"0.0 {box} ylo yhi", f"0.0 {box} zlo zhi", "",
        "Masses", "", "1 1.0", "", "Atoms # atomic", "",
    ]
    for i, p in enumerate(positions, 1):
        lines.append(f"{i} 1 {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    path.write_text("\n".join(lines) + "\n")


def write_dump(path, frames, coord="x"):
    names = "x y z" if coord == "x" else "xu yu zu"
    out = []
    for k, pos in enumerate(frames):
        out += [
            "ITEM: TIMESTEP", str(k),
            "ITEM: NUMBER OF ATOMS", str(len(pos)),
            "ITEM: BOX BOUNDS pp pp pp",
            "0.0 20.0", "0.0 20.0", "0.0 20.0",
            f"ITEM: ATOMS id type {names}",
        ]
        out += [f"{i} 1 {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}" for i, p in enumerate(pos, 1)]
    path.write_text("\n".join(out) + "\n")


def run_lat(repo_root, args):
    proc = subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return proc


def load(path):
    return np.genfromtxt(path, delimiter=",", comments="#")


@pytest.fixture
def synth(tmp_path):
    rng = np.random.default_rng(0)
    n, t = 30, 8
    steps = rng.normal(0, 0.3, size=(t, n, 3))
    steps[0] = 0.0
    frames = np.cumsum(steps, axis=0) + 10.0
    data = tmp_path / "s.data"
    dump = tmp_path / "s.lammpstrj"
    write_data(data, frames[0])
    write_dump(dump, frames, coord="xu")
    return tmp_path, data, dump, frames


def test_msd_matches_direct(repo_root, synth):
    tmp, data, dump, frames = synth
    run_lat(repo_root, ["msd", "-s", str(data), "-f", str(dump),
                        "-sel", "type 1", "-o", str(tmp / "m")])
    out = load(tmp / "m.msd.out.csv")
    # direct MSD with COM removal (matches config default)
    traj = frames - frames.mean(axis=1, keepdims=True)
    t = len(traj)
    ref = np.array([
        np.mean(np.sum((traj[lag:] - traj[:t - lag]) ** 2, axis=2))
        for lag in range(t)
    ])
    assert np.allclose(out[:, 1], ref, atol=1e-8)


def test_density_uniform(repo_root, tmp_path):
    tmp = tmp_path
    rng = np.random.default_rng(1)
    n = 4000
    pos = rng.uniform(0, 20, size=(n, 3))
    data = tmp_path / "u.data"
    dump = tmp_path / "u.lammpstrj"
    write_data(data, pos)
    write_dump(dump, [pos], coord="x")
    run_lat(repo_root, ["density", "-s", str(data), "-f", str(dump),
                        "-sel", "all", "-axis", "z", "-nbins", "20", "-quantity", "number",
                        "-o", str(tmp / "d")])
    out = load(tmp / "d.density.out.csv")
    expected = n / (20.0 ** 3)
    assert out[:, 1].mean() == pytest.approx(expected, rel=1e-9)
    assert np.allclose(out[:, 1], expected, rtol=0.2)


def test_rms_rigid_rotation(repo_root, tmp_path):
    tmp = tmp_path
    base = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], float) + 10.0
    frames = []
    for k in range(4):
        ang = 0.1 * k
        c, s = np.cos(ang), np.sin(ang)
        r = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        p = (base - 10.0) @ r.T + 10.0
        frames.append(p)
    data = tmp_path / "r.data"
    dump = tmp_path / "r.lammpstrj"
    write_data(data, frames[0])
    write_dump(dump, frames, coord="x")
    run_lat(repo_root, ["rms", "-s", str(data), "-f", str(dump),
                        "-sel", "all", "-ref", "first", "-o", str(tmp / "r")])
    out = load(tmp / "r.rmsd.out.csv")
    assert out[:, 1].max() < 1e-6


def test_density_fractional_binning_tracks_shifting_box(repo_root, tmp_path):
    """Atoms at fixed fractional coords must land in the same bin across frames
    even when the NPT box origin shifts (gmx-density-style fractional binning)."""
    data = tmp_path / "f.data"
    dump = tmp_path / "f.lammpstrj"
    # 4 atoms all at z-fraction 0.26, spread in x/y
    write_data(data, [[5.0, 5.0, 5.2], [15.0, 5.0, 5.2],
                      [5.0, 15.0, 5.2], [15.0, 15.0, 5.2]])
    # frame 1: box origin shifts by 2 and length grows to 21; atoms keep the
    # same fractional z (0.26) -> 2 + 0.26*21 = 7.46
    frames = [
        ([[5.0, 5.0, 5.2], [15.0, 5.0, 5.2], [5.0, 15.0, 5.2], [15.0, 15.0, 5.2]], 0.0, 20.0),
        ([[7.0, 5.0, 7.46], [17.0, 5.0, 7.46], [7.0, 15.0, 7.46], [17.0, 15.0, 7.46]], 2.0, 21.0),
    ]
    lines = []
    for k, (pos, lo, length) in enumerate(frames):
        hi = lo + length
        lines += ["ITEM: TIMESTEP", str(k), "ITEM: NUMBER OF ATOMS", "4",
                  "ITEM: BOX BOUNDS pp pp pp",
                  f"{lo} {hi}", f"{lo} {hi}", f"{lo} {hi}",
                  "ITEM: ATOMS id type x y z"]
        lines += [f"{i} 1 {p[0]} {p[1]} {p[2]}" for i, p in enumerate(pos, 1)]
    dump.write_text("\n".join(lines) + "\n")

    run_lat(repo_root, ["density", "-s", str(data), "-f", str(dump),
                        "-axis", "z", "-nbins", "20", "-o", str(tmp_path / "f")])
    prof = load(tmp_path / "f.density.out.csv")
    nonzero = prof[prof[:, 1] > 0]
    # both frames land in the same fractional bin -> exactly one populated bin
    assert len(nonzero) == 1
    # axis scaled to mean box (origin 1, length 20.5); bin centre fraction 0.275
    assert nonzero[0, 0] == pytest.approx(1.0 + 0.275 * 20.5, abs=0.6)
    assert "box length varies" in (tmp_path / "f.density.out.csv").read_text()


def test_density_rejects_smooth_flag(repo_root, tmp_path):
    """Smoothing was moved to plot time; density must no longer accept it."""
    data = tmp_path / "g.data"
    dump = tmp_path / "g.lammpstrj"
    write_data(data, [[5.0, 5.0, 5.0]])
    write_dump(dump, [[[5.0, 5.0, 5.0]]], coord="x")
    proc = subprocess.run(
        [sys.executable, str(repo_root / "lat"), "density", "-s", str(data),
         "-f", str(dump), "-axis", "z", "-smooth", "gaussian", "-o", str(tmp_path / "g")],
        capture_output=True, text=True)
    assert proc.returncode != 0
    assert "error" in proc.stderr.lower()

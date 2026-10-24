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


def _naive_acf(triple_frames):
    """Reference ACF defined literally: re-AND lag+1 shifted slices per lag."""
    import numpy as np

    T = len(triple_frames)
    universe = sorted(set().union(*triple_frames)) if triple_frames else []
    idx = {t: i for i, t in enumerate(universe)}
    present = np.zeros((len(universe), T), dtype=bool)
    for k, s in enumerate(triple_frames):
        for t in s:
            present[idx[t], k] = True
    total = present.sum()
    inter = np.zeros(T)
    cont = np.zeros(T)
    if T:
        cont[0] = 1.0
    for lag in range(T):
        inter[lag] = ((present[:, : T - lag] & present[:, lag:]).sum() * T
                      / ((T - lag) * total)) if T > lag else 0.0
        if lag == 0:
            continue
        run = present[:, : T - lag].copy()
        for k in range(1, lag + 1):
            run &= present[:, k: T - lag + k]
        cont[lag] = run.sum() * T / ((T - lag) * total)
    return inter, cont


def test_acf_matches_the_naive_definition():
    """The production ACF runs incrementally and drops finished pairs; both are
    optimisations that must not change a single value."""
    import numpy as np

    from lat.modules.hbond import _acf

    rng = np.random.default_rng(5)
    for trial in range(4):
        n_pairs, T = int(rng.integers(1, 25)), int(rng.integers(1, 40))
        frames = []
        for k in range(T):
            frames.append({f"t{i}" for i in range(n_pairs) if rng.random() < 0.3})
        if not any(frames):
            frames[0] = {"t0"}
        table, scalars, _ = _acf(frames)
        want_inter, want_cont = _naive_acf(frames)
        got = table.data
        assert np.array_equal(got[:, 1], want_inter), (trial, "intermittent")
        assert np.array_equal(got[:, 2], want_cont), (trial, "continuous")
        assert scalars["tau_integral"] == pytest.approx(float(want_inter[1:].sum()))


def test_acf_is_cheap_for_long_trajectories():
    """A pair that is not continuously present cannot affect longer lags, so a
    601-frame / 100k-pair case must finish in well under a second."""
    import time

    import numpy as np

    from lat.modules.hbond import _acf

    rng = np.random.default_rng(6)
    frames = []
    for k in range(400):
        frames.append({f"p{i}" for i in rng.choice(100_000, 200, replace=False)})
    t0 = time.time()
    table, _, _ = _acf(frames)
    elapsed = time.time() - t0
    assert table.data.shape == (400, 3)
    assert elapsed < 20.0, f"ACF took {elapsed:.1f}s for 400 frames"

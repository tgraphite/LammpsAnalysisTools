"""Level-2 regression: lat rdf on the reduced sample.

Validates (a) LAMMPS and internal-GROMACS input paths agree bit-for-bit and
(b) g(r) agrees with an independent implementation (MDAnalysis InterRDF).
"""

import subprocess
import sys

import numpy as np
import pytest


def run_lat(repo_root, args):
    proc = subprocess.run(
        [sys.executable, str(repo_root / "lat"), *args],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"lat failed:\n{proc.stderr}"
    return proc


def load(path):
    return np.genfromtxt(path, delimiter=",", comments="#")


def test_rdf_two_paths_and_mda(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    # MDAnalysis InterRDF includes intramolecular pairs, so request them explicitly.
    common = ["rdf", "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "10", "-intramol"]

    run_lat(repo_root, [*common, "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                        "-o", str(tmp_path / "lmp")])
    run_lat(repo_root, [*common, "-s", str(base) + ".gro", "-f", str(base) + ".xtc",
                        "-o", str(tmp_path / "gmx")])

    lmp = load(tmp_path / "lmp.rdf.out.csv")
    gmx = load(tmp_path / "gmx.rdf.out.csv")
    lmp_cn = load(tmp_path / "lmp.rdf.cn.out.csv")

    # (a) format independence
    assert lmp.shape == gmx.shape
    assert np.array_equal(lmp[:, 0], gmx[:, 0])
    assert np.allclose(lmp[:, 1], gmx[:, 1], atol=0.0, rtol=0.0)

    # (b) independent reference
    mda = pytest.importorskip("MDAnalysis")
    from MDAnalysis.analysis.rdf import InterRDF

    u = mda.Universe(str(base) + ".data")
    u.load_new(str(base) + ".lammpstrj", format="LAMMPSDUMP")
    rdf = InterRDF(u.select_atoms("type 2"), u.select_atoms("type 5"),
                   nbins=200, range=(0.0, 10.0))
    rdf.run()
    res = rdf.results
    assert np.allclose(np.asarray(res.bins), lmp[:, 0])
    assert np.abs(np.asarray(res.rdf) - lmp[:, 1]).max() < 0.05

    # CN must be exactly reconstructible from the reference pair counts
    n_frames = len(u.trajectory)
    cn_ref = np.cumsum(np.asarray(res.count)) / (n_frames * len(u.select_atoms("type 2")))
    assert np.allclose(cn_ref, lmp_cn[:, 1], atol=1e-12)


def test_rdf_intramolecular_default_excludes(repo_root, sample_dir, tmp_path):
    """Default excludes same-molecule pairs; -intramol includes them."""
    base = sample_dir / "sys1_reduced"
    common = ["rdf", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
              "-sel1", "type 1", "-sel2", "type 1", "-cutoff", "10"]

    run_lat(repo_root, [*common, "-o", str(tmp_path / "default")])
    run_lat(repo_root, [*common, "-intramol", "-o", str(tmp_path / "intra")])

    default = load(tmp_path / "default.rdf.out.csv")[:, 1]
    intra = load(tmp_path / "intra.rdf.out.csv")[:, 1]

    # intramolecular pairs exist for type 1, so the two must differ
    assert np.abs(default - intra).max() > 1e-6
    # excluding a subset of pairs can only lower g(r) pointwise
    assert np.all(default <= intra + 1e-12)

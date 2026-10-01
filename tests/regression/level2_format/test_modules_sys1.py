"""End-to-end smoke: every module runs on the real reduced sample.

This is the integration layer that unit tests on synthetic systems miss
(e.g. frames with zero hydrogen bonds).
"""

import subprocess
import sys

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
    ("hbond", ["hbond"]),
    ("agg", ["agg", "-sel1", "type 2", "-sel2", "type 4 5 6 7", "-cutoff", "2.5"]),
    ("cluster", ["cluster", "-sel1", "type 2", "-sel2", "type 4 5 6 7",
                 "-cutoff", "2.5", "-focus"]),
    ("pistack", ["pistack", "-sel1", "type 3"]),
    ("probe_ffv", ["probe", "-mode", "ffv", "-grid_dr", "3.0"]),
    ("probe_sasa", ["probe", "-mode", "sasa", "-sel", "type 2"]),
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


def test_hbond_geometry_with_empty_frames(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    proc = lat(repo_root, ["hbond", "-s", str(base) + ".data", "-f", str(base) + ".lammpstrj",
                           "-o", str(tmp_path / "h")])
    assert proc.returncode == 0, proc.stderr
    for tag in ("hbond", "hbond.acf", "hbond.dist_xz", "hbond.dist_yz", "hbond.dist_xy"):
        assert (tmp_path / f"h.{tag}.out.csv").is_file()


def test_transport_guard_on_sample(repo_root, sample_dir, tmp_path):
    base = sample_dir / "sys1_reduced"
    proc = lat(repo_root, ["transport", "-s", str(base) + ".data",
                           "-f", str(base) + ".lammpstrj", "-sel", "type 2",
                           "-o", str(tmp_path / "t")])
    assert proc.returncode == 1
    assert "velocity" in proc.stderr.lower()

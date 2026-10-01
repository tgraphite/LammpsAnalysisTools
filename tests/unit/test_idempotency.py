"""Idempotency: running any subcommand twice must give identical results.

The only permitted difference is the audit timestamp line (# generated:).
Triggering the second run also exercises output backup rotation (.bak.1).
"""

import re
import subprocess
import sys

import pytest

SAMPLE = "tests/regression/samples/sys1_reduced/sys1_reduced"

CASES = [
    ("rdf", ["rdf", "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "10"]),
    ("density", ["density", "-sel", "type 2", "-axis", "z"]),
    ("density2d", ["density", "-sel", "type 2", "-dim", "2", "-nbins", "20"]),
    ("gyrate", ["gyrate", "-sel", "mol 401-600"]),
    ("contact", ["contact", "-sel1", "type 2", "-sel2", "type 5", "-cutoff", "3"]),
    ("rms", ["rms", "-sel", "type 3"]),
    ("msd", ["msd", "-sel", "type 2"]),
    ("hbond", ["hbond"]),
    ("agg", ["agg", "-sel1", "type 2", "-sel2", "type 4 5 6 7", "-cutoff", "2.5"]),
    ("cluster", ["cluster", "-sel1", "type 2", "-sel2", "type 4 5 6 7", "-cutoff", "2.5", "-focus"]),
    ("pistack", ["pistack", "-sel1", "type 3"]),
    ("probe_ffv", ["probe", "-mode", "ffv", "-grid_dr", "4.0"]),
    ("probe_sasa", ["probe", "-mode", "sasa", "-sel", "type 2"]),
    ("residence", ["residence", "-ref", "type 2", "-lig", "type 4 5 6 7", "-cutoff", "3"]),
    ("occupancy", ["occupancy", "-sel", "type 2", "-zone", "z:0-45,45-90,90-135"]),
]

_IGNORE = re.compile(r"^#\s*(generated|command):")


def _stable(path):
    lines = []
    for line in open(path):
        if _IGNORE.match(line):
            continue
        lines.append(line)
    return "".join(lines)


def _run(repo_root, args):
    proc = subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return proc


@pytest.mark.parametrize("name,args", CASES, ids=[c[0] for c in CASES])
def test_idempotent(name, args, repo_root, tmp_path):
    base = list(args) + ["-o", str(tmp_path / name), "-s", SAMPLE + ".data"]
    if name != "select":
        base += ["-f", SAMPLE + ".lammpstrj"]
    _run(repo_root, base)
    files = sorted(tmp_path.glob(f"{name}.*.out.csv"))
    assert files, f"{name}: no output"
    first = {f.name: _stable(f) for f in files}

    _run(repo_root, base)  # same -o -> rotates previous outputs to .bak.1
    for f in files:
        assert _stable(f) == first[f.name], f"{name}: {f.name} changed on re-run"
    assert list(tmp_path.glob("*.bak.1")), f"{name}: expected .bak.1 rotation"


def test_screen_modules_idempotent(repo_root, sample_dir):
    data = str(sample_dir / "sys1_reduced.data")
    for args in (["select", "-sel", "type 2"], ["describe", "-module", "rdf"],
                 ["describe"]):
        p1 = _run(repo_root, args + ([] if "select" not in args else ["-s", data]))
        p2 = _run(repo_root, args + ([] if "select" not in args else ["-s", data]))
        assert p1.stdout == p2.stdout

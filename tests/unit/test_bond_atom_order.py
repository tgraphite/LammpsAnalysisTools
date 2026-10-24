"""The Bonds section references atom *ids*, so it must follow the id re-sort.

A LAMMPS data file's Atoms section is not required to be sorted by id (OVITO
writes it in its own order).  `lammps_data.parse` re-sorts every per-atom array
into id order; the bond list, built from `id2idx`, therefore has to be permuted
by the same map.  When it is not, the bond graph becomes a random relabelling
and `make_whole` walks non-bonded pairs, which silently corrupts gyrate / rms /
hbond.  Regression: an OVITO-written PET file had 880 of 6476 edges matching
real bonds and a 3.6x overestimated Rg.
"""

import subprocess
import sys

import numpy as np
import pytest


def write_scrambled_data(path, box=20.0):
    """4 atoms / 2 two-atom molecules, Atoms rows deliberately not id-sorted.

    Molecule 1 (ids 1,2) straddles the x face: |0.5 - 18.5| wraps to 2.0 A,
    so its true Rg is 1.0 A whereas the wrapped coordinates give 9.0 A.
    """
    path.write_text(
        "\n".join([
            "LAMMPS data file (unsorted Atoms section)", "",
            "4 atoms", "2 bonds", "1 atom types", "",
            f"0.0 {box} xlo xhi", f"0.0 {box} ylo yhi", f"0.0 {box} zlo zhi", "",
            "Masses", "", "1 12.011", "",
            "Atoms # molecular", "",
            "3 2 1 5.0 5.0 5.0",
            "1 1 1 0.5 5.0 5.0",
            "4 2 1 7.0 5.0 5.0",
            "2 1 1 18.5 5.0 5.0",
            "",
            "Bonds", "",
            "1 1 1 2",
            "2 1 3 4",
            "",
        ]) + "\n"
    )


def test_bonds_are_permuted_with_atom_order(tmp_path):
    from lat.core.io import lammps_data

    f = tmp_path / "scrambled.data"
    write_scrambled_data(f)
    topo = lammps_data.parse(str(f))

    assert topo.atom_id.tolist() == [1, 2, 3, 4]
    got = {tuple(sorted(b)) for b in topo.bonds}
    assert got == {(0, 1), (2, 3)}, got
    # each bond must stay inside one molecule
    for a, b in topo.bonds:
        assert topo.mols[a] == topo.mols[b]


def test_gyrate_makes_straddling_molecule_whole(tmp_path, repo_root):
    f = tmp_path / "scrambled.data"
    write_scrambled_data(f)
    proc = subprocess.run(
        [sys.executable, str(repo_root / "lat"), "gyrate",
         "-s", str(f), "-sel", "mol 1", "-o", str(tmp_path / "mol1")],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    rg = np.genfromtxt(tmp_path / "mol1.gyrate.out.csv", delimiter=",", comments="#")[1]
    # two equal masses 2.0 A apart -> Rg = 1.0 A; wrapped coords would give 9.0
    assert rg == pytest.approx(1.0, abs=1e-6)

"""Translation invariance of every PBC-sensitive module.

A rigid translation of the whole system by a non-lattice vector leaves every
physical distance unchanged; it only changes *where* molecules sit relative to
the cell faces. Any observable that handles PBC correctly must therefore be
identical in both placements. Before 0.2.0 this failed for

* ``contact -mode center`` and ``residence -mode center`` (averaged wrapped
  coordinates, so a molecule on a face got the wrong centre),
* ``pistack`` (ring centroid/normal from wrapped coordinates),
* ``agg``/``cluster -mode center`` (first-atom minimum image).

``rdf -intramol`` is included because it is *pairwise* minimum image and is
therefore already invariant; the test pins that so a future "optimisation"
cannot silently break it.
"""

from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

BOX = 20.0
# hexagons: 6-membered rings, one straddling the x face, one in open space
RING_R = 1.4


def _hexagon(cx, cy, cz, n=6):
    a = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([cx + RING_R * np.cos(a), cy + RING_R * np.sin(a),
                            np.full(n, cz)])


# The hexagon at x=0 straddles the face, so its wrapped centroid lands at
# x=10 while its physical centre is x=0. Cuts chosen between the two:
# the true centre separation is 4 A, the wrapped-coordinate one is 6 A.
def _system(shift):
    """Two rigid hexagons; `shift` translates the whole system along x."""
    mols, pos = [], []
    for mol_id, (cx, cy) in ((1, (0.0, 10.0)), (2, (4.0, 10.0))):
        pts = _hexagon(cx, cy, 10.0)
        mols += [mol_id] * len(pts)
        pos.append(pts)
    pos = np.vstack(pos)
    pos[:, 0] = np.mod(pos[:, 0] + shift, BOX)
    mols = np.array(mols)
    bonds = []
    for a in range(0, 12, 6):
        for k in range(6):
            bonds.append((a + k, a + (k + 1) % 6))
    return pos, mols, np.array(bonds)


def write_system(data_path, dump_path, pos, mols, bonds):
    n = len(pos)
    lines = ["LAMMPS data file", "", f"{n} atoms", f"{len(bonds)} bonds",
             "2 atom types", "",
             f"0.0 {BOX} xlo xhi", f"0.0 {BOX} ylo yhi", f"0.0 {BOX} zlo zhi", "",
             "Masses", "", "1 12.011", "2 15.999", "",
             "Atoms # molecular", ""]
    for i, (m, p) in enumerate(zip(mols, pos), 1):
        lines.append(f"{i} {m} 1 {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    lines += ["", "Bonds", ""]
    for k, (a, b) in enumerate(bonds, 1):
        lines.append(f"{k} 1 {a + 1} {b + 1}")
    data_path.write_text("\n".join(lines) + "\n")

    out = ["ITEM: TIMESTEP", "0", "ITEM: NUMBER OF ATOMS", str(n),
           "ITEM: BOX BOUNDS pp pp pp",
           f"0.0 {BOX}", f"0.0 {BOX}", f"0.0 {BOX}",
           "ITEM: ATOMS id mol type x y z"]
    for i, (m, p) in enumerate(zip(mols, pos), 1):
        out.append(f"{i} {m} 1 {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    dump_path.write_text("\n".join(out) + "\n")


def run_lat(repo_root, args):
    proc = subprocess.run([sys.executable, str(repo_root / "lat"), *args],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    return proc


def read_table(path):
    """(column names, data) from a LAT CSV; names come from the `# columns:` header."""
    names = None
    for line in path.read_text().splitlines():
        if line.startswith("# columns:"):
            names = [s.strip() for s in line.split(":", 1)[1].split(",")]
    data = np.atleast_2d(np.genfromtxt(path, delimiter=",", comments="#"))
    assert names is not None, f"no # columns: header in {path}"
    assert len(names) == data.shape[1], f"{path}: {names} vs {data.shape}"
    return names, data


@pytest.fixture
def placements(repo_root, tmp_path):
    """(out_dir_A, out_dir_B) for the same system shifted by 7.3 Å in x."""
    dirs = []
    for tag, shift in (("a", 0.0), ("b", 7.3)):
        d = tmp_path / tag
        d.mkdir()
        pos, mols, bonds = _system(shift)
        data, dump = d / "sys.data", d / "sys.lammpstrj"
        write_system(data, dump, pos, mols, bonds)
        dirs.append((d, data, dump))
    return dirs


@pytest.mark.parametrize(
    "module,extra,prefix,tag,cols",
    [
        ("contact", ["-sel1", "mol 1", "-sel2", "mol 2", "-cutoff", "7.0", "-mode", "center"],
         "c", "c.contact.out.csv", ("n_contacts", "min_dist", "com_dist")),
        ("contact", ["-sel1", "mol 1", "-sel2", "mol 2", "-cutoff", "7.0", "-mode", "min"],
         "cmin", "cmin.contact.out.csv", ("n_contacts", "min_dist", "com_dist")),
        ("residence", ["-ref", "mol 2", "-lig", "mol 1", "-cutoff", "5.0", "-mode", "center"],
         "r", "r.residence.out.csv", ("n_bound",)),
        ("residence", ["-ref", "mol 2", "-lig", "mol 1", "-cutoff", "5.0", "-mode", "min"],
         "rmin", "rmin.residence.out.csv", ("n_bound",)),
        ("pistack", ["-sel1", "mol 1", "-sel2", "mol 2", "-distance", "5.0", "-angle", "30"],
         "p", "p.pistack.out.csv", ("n_pistack",)),
        ("agg", ["-sel1", "mol 1", "-sel2", "mol 2", "-cutoff", "7.0", "-mode", "center"],
         "g", "g.agg.out.csv", None),
        ("cluster", ["-sel1", "mol 1", "-sel2", "mol 2", "-cutoff", "7.0", "-mode", "center"],
         "cl", "cl.cluster.cn.out.csv", None),
        ("rdf", ["-sel1", "mol 1", "-sel2", "mol 1", "-nbins", "20", "-cutoff", "8.0",
                 "-intramol"], "d", "d.rdf.out.csv", None),
    ],
)
def test_module_is_translation_invariant(repo_root, placements, module, extra, prefix, tag, cols):
    results = []
    for d, data, dump in placements:
        run_lat(repo_root, [module, "-s", str(data), "-f", str(dump),
                            "-o", str(d / prefix), *extra])
        assert (d / tag).exists(), f"{module} did not write {tag}"
        results.append(read_table(d / tag))
    (names_a, a), (names_b, b) = results
    assert names_a == names_b
    assert len(a) == len(b)
    for k, name in enumerate(names_a):
        if cols is not None and name not in cols:
            continue
        va, vb = a[:, k], b[:, k]
        bad = ~(np.isclose(va, vb, rtol=1e-9, atol=1e-9, equal_nan=True))
        assert not bad.any(), (
            f"{module} {tag} column {name!r} changed under a rigid translation: "
            f"{va[bad]} vs {vb[bad]} -- the module treats wrapped coordinates as physical"
        )

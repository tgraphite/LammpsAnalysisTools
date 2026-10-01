"""Unit tests for the core primitives (no external references needed)."""

import numpy as np
import pytest

from lat.config import load as load_config
from lat.core.pbc import Box
from lat.core.selection import parse as parse_selection
from lat.core.table import Field, Table
from lat.core.output import write_table
from lat.core.topology import element_mass_table, infer_type_elements


def test_pbc_orthorhombic_minimum_image():
    box = Box.orthorhombic([10.0, 10.0, 10.0])
    d = box.minimum_image(np.array([[9.5, 0.0, 0.0], [-9.5, 0.0, 0.0]]))
    assert np.allclose(d[:, 0], [-0.5, 0.5])
    assert box.volume == pytest.approx(1000.0)
    assert box.is_orthorhombic


def test_pbc_triclinic_minimum_image():
    box = Box.from_lammps_data([0, 0, 0], [10, 10, 10], [2.0, 1.0, 0.5])
    assert not box.is_orthorhombic
    d = np.array([7.3, -4.1, 2.2])
    # brute-force minimum over neighboring images
    best = None
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                img = d + np.array([i, j, k]) @ box.H.T
                dist = np.linalg.norm(img)
                best = dist if best is None else min(best, dist)
    assert np.linalg.norm(box.minimum_image(d[None, :])[0]) == pytest.approx(best)


def test_element_inference_tolerance():
    table = element_mass_table(load_config())
    te, problems = infer_type_elements({1: 6.94, 2: 12.011}, table, 0.002)
    assert te[1] == "Li" and te[2] == "C"
    # reax_tools' Li=7.0 must not match within 0.2%
    te2, problems2 = infer_type_elements({1: 7.0}, table, 0.002)
    assert te2[1] is None and 1 in problems2


def _toy_topo():
    from lat.core.topology import Topology

    n = 6
    return Topology(
        n_atoms=n,
        atom_id=np.arange(1, n + 1),
        types=np.array([1, 1, 1, 2, 2, 2]),
        box=Box.orthorhombic([10, 10, 10]),
        positions=np.zeros((n, 3)),
        mols=np.array([1, 1, 1, 2, 2, 2]),
        type_element={1: "C", 2: "O"},
    )


def test_selection_grammar():
    topo = _toy_topo()
    assert parse_selection(topo, "all").n_atoms == 6
    assert parse_selection(topo, "type 2").n_atoms == 3
    assert parse_selection(topo, "type 1 2").n_atoms == 6
    assert parse_selection(topo, "element O").n_atoms == 3
    assert parse_selection(topo, "mol 2").n_atoms == 3
    assert parse_selection(topo, "element O and mol 1").n_atoms == 0
    assert parse_selection(topo, "type 1 and not mol 1").n_atoms == 0
    assert parse_selection(topo, "id 1-3").n_atoms == 3
    assert parse_selection(topo, "index 0 1").n_atoms == 2


def test_output_backup_rotation(tmp_path):
    t = Table(
        title="t", columns=[Field("x", "A"), Field("y", "1")],
        data=np.array([[1.0, 2.0]]), tag="demo", meta={"module": "demo"},
    )
    p1 = write_table(t, basename="s", outdir=tmp_path, backup_on_conflict=True)
    p2 = write_table(t, basename="s", outdir=tmp_path, backup_on_conflict=True)
    assert p1 == p2
    assert (tmp_path / "s.demo.out.csv.bak.1").is_file()

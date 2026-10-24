"""Compatibility fixes found by the real-project campaign (tests/compat)."""

import numpy as np
import pytest

from lat.core.io import lammps_data, lammps_dump
from lat.core.pbc import Box


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_atom_style_accelerator_suffix(tmp_path):
    data = """LAMMPS data file

4 atoms
2 atom types

0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Masses

1 12.0
2 1.0

Atoms # full/kk

1 1 1 0.0 0.0 0.0 0.0
2 1 1 0.0 1.0 0.0 0.0
3 2 1 0.0 0.0 1.0 0.0
4 2 1 0.0 0.0 0.0 1.0
"""
    topo = lammps_data.parse(_write(tmp_path, "a.data", data))
    assert topo.n_atoms == 4
    assert np.allclose(topo.positions[1], [1.0, 0.0, 0.0])


def test_dump_string_element_column(tmp_path):
    dump = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z element
1 1 1.0 2.0 3.0 Ni
2 2 4.0 5.0 6.0 N
"""
    p = _write(tmp_path, "d.lammpstrj", dump)
    frames = list(lammps_dump.DumpReader(p, n_atoms=2))
    assert len(frames) == 1
    fr = frames[0]
    assert np.allclose(fr.data["x"], [1.0, 4.0])
    assert list(fr.data["element"]) == ["Ni", "N"]


def test_element_map_falls_back_to_dump_column(tmp_path):
    from lat.core.context import Context
    from lat.core.io import DumpFrameSource
    data = """LAMMPS data file

2 atoms
2 atom types

0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Atoms # atomic

1 1 0 0 0
2 2 1 1 1
"""
    dump = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z element
1 1 1.0 2.0 3.0 Cu
2 2 4.0 5.0 6.0 O
"""
    topo = lammps_data.parse(_write(tmp_path, "n.data", data))
    src = DumpFrameSource(_write(tmp_path, "n.lammpstrj", dump), topo.n_atoms)
    ctx = Context(topo, src)
    assert ctx.element_map() == {1: "Cu", 2: "O"}


def test_dump_atom_order_realigned_by_id(tmp_path):
    """A dump may list atoms out of id order; frames must be realigned."""
    from lat.core.context import Context
    from lat.core.io import open_system
    import lat.modules as M

    M._load_builtin()
    data = """LAMMPS data file

4 atoms
2 atom types

0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Atoms # atomic

1 1 0 0 0
2 1 1 0 0
3 2 5 0 0
4 2 5 1.2 0
"""
    dump = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
4
ITEM: BOX BOUNDS pp pp pp
0 10
0 10
0 10
ITEM: ATOMS id type x y z
4 2 5 1.2 0
3 2 5 0 0
2 1 1 0 0
1 1 0 0 0
"""
    topo, src = open_system(_write(tmp_path, "o.data", data),
                            _write(tmp_path, "o.lammpstrj", dump))
    ctx = Context(topo, src)
    table = M.get("rdf")().run(ctx, sel1="type 2", sel2="type 2", nbins=200,
                              cutoff=4.0, frame_range=None, stride=None,
                              intramol=False)[0]
    r, g = table.data[:, 0], table.data[:, 1]
    # the two type-2 atoms are 1.2 A apart; a mis-aligned order would peak at 1.0
    assert abs(r[int(np.argmax(g))] - 1.2) < 0.05


def test_vdw_table_covers_metals_and_metalloids():
    import yaml
    from pathlib import Path
    cfg = yaml.safe_load((Path(__file__).resolve().parents[2] / "config" / "lat.yaml").read_text())
    vdw = cfg["vdw"]
    for el in ("Te", "Ti", "Mo", "Fe", "U"):
        assert el in vdw and float(vdw[el]) > 0
    assert len(vdw) > 90


def test_freevol_grid_auto_coarsens_large_box(tmp_path):
    from lat.config import load as load_config
    from lat.core.context import Context
    from lat.core.io import open_system
    from lat.modules import get
    data = """LAMMPS data file

1 atoms
1 atom types

0 2000 xlo xhi
0 2000 ylo yhi
0 2000 zlo zhi

Masses

1 12.0

Atoms # atomic

1 1 1000 1000 1000
"""
    dump = """ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0 2000
0 2000
0 2000
ITEM: ATOMS id type x y z
1 1 1000 1000 1000
"""
    topo = _write(tmp_path, "big.data", data)
    traj = _write(tmp_path, "big.lammpstrj", dump)
    cfg = load_config()
    topo_o, src = open_system(str(topo), str(traj), cfg=cfg)
    ctx = Context(topo_o, src, cfg=cfg)
    # grid_dr = 0.5 A over a 2000 A box would be 6.4e10 voxels; must coarsen
    tables = get("freevol")().run(ctx, radius=1.445, mode="ffv", grid_dr=0.5,
                                  connectivity=26, psd_rmax=None, components=False,
                                  frame_range="0:1:1", stride=None)
    assert tables and tables[0].data.shape[0] == 1
    assert any("coarsened" in w for w in tables[0].warnings)
    # one atom smaller than a coarsened voxel: the grid cannot see it
    assert 0.99 <= tables[0].data[0, 1] <= 1.0

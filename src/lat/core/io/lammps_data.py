"""LAMMPS data-file parser (atom_style inference included)."""

from __future__ import annotations

import numpy as np

from ..errors import LatError, get_logger
from ..pbc import Box
from ..topology import Topology, element_mass_table, infer_type_elements

STYLE_COLUMNS = {
    "atomic": ["id", "type", "x", "y", "z"],
    "charge": ["id", "type", "q", "x", "y", "z"],
    "molecular": ["id", "mol", "type", "x", "y", "z"],
    "bond": ["id", "mol", "type", "x", "y", "z"],
    "angle": ["id", "mol", "type", "x", "y", "z"],
    "full": ["id", "mol", "type", "q", "x", "y", "z"],
}

_SECTION_WORDS = {
    "Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
    "Ellipsoids", "Lines", "Triangles", "Bodies", "Atom", "Pair", "Bond", "Angle",
    "Dihedral", "Improper",
}
_COUNT_WORDS = {"atoms", "bonds", "angles", "dihedrals", "impropers", "types", "ellipsoids", "lines", "triangles", "bodies"}


def _float_ok(tok: str) -> bool:
    return any(c in tok for c in ".eE")


def _infer_style(tokens: list[str], has_bonds: bool) -> str:
    n = len(tokens)
    if n == len(STYLE_COLUMNS["atomic"]):
        return "atomic"
    if n == len(STYLE_COLUMNS["full"]):
        return "full"
    if n == 6:
        if _float_ok(tokens[2]):
            return "charge"
        if has_bonds:
            return "molecular"
        raise LatError(
            "ambiguous 6-column Atoms section (charge vs molecular); pass -style explicitly"
        )
    raise LatError(f"unsupported Atoms section with {n} columns; pass -style explicitly")


def parse(path, *, style: str | None = None, cfg=None) -> Topology:
    with open(path) as fh:
        raw_lines = fh.readlines()

    natoms = None
    box_lo = [None, None, None]
    box_hi = [None, None, None]
    tilt = None
    style_hint = None
    sections: dict[str, list[str]] = {}
    current = None

    for raw in raw_lines:
        line = raw.rstrip("\n")
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        head = tokens[0]

        if not head[:1].isdigit() and head in _SECTION_WORDS:
            low = head.lower()
            if low == "masses":
                current = "masses"
            elif low == "atoms":
                current = "atoms"
                if "#" in tokens:
                    style_hint = tokens[tokens.index("#") + 1]
            elif low == "velocities":
                current = "velocities"
            elif low == "bonds":
                current = "bonds"
            else:
                current = None
            if current:
                sections.setdefault(current, [])
            continue

        if current is not None:
            data = line.split("#")[0].strip()
            if data:
                sections[current].append(data)
            continue

        # pre-section header
        if tokens[0].isdigit() and tokens[-1].lower() in _COUNT_WORDS:
            if tokens[-1] == "atoms":
                natoms = int(tokens[0])
            continue
        if len(tokens) >= 4 and tokens[2:4] == ["xlo", "xhi"]:
            box_lo[0], box_hi[0] = float(tokens[0]), float(tokens[1])
        elif len(tokens) >= 4 and tokens[2:4] == ["ylo", "yhi"]:
            box_lo[1], box_hi[1] = float(tokens[0]), float(tokens[1])
        elif len(tokens) >= 4 and tokens[2:4] == ["zlo", "zhi"]:
            box_lo[2], box_hi[2] = float(tokens[0]), float(tokens[1])
        elif len(tokens) >= 6 and tokens[3:6] == ["xy", "xz", "yz"]:
            tilt = [float(tokens[0]), float(tokens[1]), float(tokens[2])]

    if natoms is None:
        raise LatError(f"{path}: 'N atoms' header not found")
    if any(v is None for v in box_lo + box_hi):
        raise LatError(f"{path}: incomplete box bounds")

    atom_rows = [r.split() for r in sections.get("atoms", [])]
    if not atom_rows:
        raise LatError(f"{path}: Atoms section missing or empty")
    if len(atom_rows) != natoms:
        raise LatError(f"{path}: header says {natoms} atoms but Atoms section has {len(atom_rows)}")

    has_bonds = bool(sections.get("bonds"))
    st = style or style_hint or _infer_style(atom_rows[0], has_bonds)
    # strip accelerator suffixes: full/kk, charge/omp, atomic/gpu, ...
    st = st.split("/")[0].strip().lower()
    if st not in STYLE_COLUMNS:
        raise LatError(f"{path}: unknown atom_style {st!r}; known: {sorted(STYLE_COLUMNS)}")
    if style and style_hint and style.split("/")[0] != style_hint.split("/")[0]:
        get_logger().warning(f"-style {style} overrides Atoms header style {style_hint}")
    cols = STYLE_COLUMNS[st]

    table = np.array([[float(v) for v in row[: len(cols)]] for row in atom_rows], float)
    col = {name: table[:, k] for k, name in enumerate(cols)}

    ids = col["id"].astype(int)
    id2idx = {int(a): i for i, a in enumerate(ids)}
    types = col["type"].astype(int)
    positions = np.column_stack([col["x"], col["y"], col["z"]])
    charges = col.get("q")
    mols = col["mol"].astype(int) if "mol" in col else None

    type_mass: dict[int, float] = {}
    for row in sections.get("masses", []):
        t = row.split()
        if len(t) >= 2:
            type_mass[int(t[0])] = float(t[1])

    velocities = None
    vel_rows = sections.get("velocities", [])
    if vel_rows:
        vid = np.array([int(r.split()[0]) for r in vel_rows])
        vv = np.array([[float(x) for x in r.split()[1:4]] for r in vel_rows])
        velocities = np.full((natoms, 3), np.nan)
        for k, aid in enumerate(vid):
            if int(aid) in id2idx:
                velocities[id2idx[int(aid)]] = vv[k]

    bonds = None
    bond_rows = sections.get("bonds", [])
    if bond_rows:
        b = []
        for r in bond_rows:
            t = r.split()
            if len(t) >= 4:
                a1, a2 = int(t[2]), int(t[3])
                if a1 in id2idx and a2 in id2idx:
                    b.append((id2idx[a1], id2idx[a2]))
        bonds = np.array(b, dtype=int) if b else None

    box = Box.from_lammps_data(box_lo, box_hi, tilt)

    cfg = cfg
    mass_table = element_mass_table(cfg) if cfg is not None else {}
    type_element, problems = ({}, {})
    if type_mass and mass_table:
        tol = float(cfg.get("element_match_tolerance", default=0.002)) if cfg is not None else 0.002
        type_element, problems = infer_type_elements(type_mass, mass_table, tol)
        if problems:
            for t, msg in problems.items():
                get_logger().info(f"type {t}: {msg}")

    order = np.argsort(ids, kind="stable")
    def _re(x):
        return None if x is None else x[order]

    return Topology(
        n_atoms=natoms,
        atom_id=ids[order],
        types=types[order],
        box=box,
        positions=positions[order],
        velocities=_re(velocities),
        masses=None,
        charges=_re(charges),
        mols=_re(mols),
        bonds=bonds,
        type_mass=type_mass,
        type_element=type_element,
        element_problems=problems,
        source=str(path),
        meta={"style": st, "format": "lammps-data"},
    )

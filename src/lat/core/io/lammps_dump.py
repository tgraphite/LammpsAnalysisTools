"""Streaming LAMMPS dump reader (frame index = order, LAMMPS timestep ignored)."""

from __future__ import annotations

import numpy as np

from ..errors import LatError
from ..frame import Frame
from ..pbc import Box

_INT_COLS = {"id", "type", "mol", "ix", "iy", "iz"}
_STR_COLS = {"element", "name", "resname", "resid"}


def _columns(names, raw):
    """Build per-column arrays; tolerate non-numeric columns (e.g. element)."""
    try:
        arr = np.array(raw, dtype=float)
    except ValueError:
        strarr = np.array(raw, dtype=object)
        data = {}
        for k, name in enumerate(names):
            col = strarr[:, k]
            if name in _INT_COLS:
                data[name] = col.astype(int)
            elif name in _STR_COLS:
                data[name] = col
            else:
                try:
                    data[name] = col.astype(float)
                except ValueError:
                    data[name] = col
        return data
    data = {}
    for k, name in enumerate(names):
        col = arr[:, k]
        data[name] = col.astype(int) if name in _INT_COLS else col
    return data


class DumpReader:
    """Iterable of Frame; safe to iterate multiple times (reopens the file)."""

    def __init__(self, path, n_atoms: int | None = None):
        self.path = str(path)
        self.expected = n_atoms

    def __iter__(self):
        with open(self.path) as fh:
            index = -1
            n_atoms = None
            triclinic = False
            while True:
                line = fh.readline()
                if not line:
                    break
                if line.startswith("ITEM: TIMESTEP"):
                    fh.readline()
                elif line.startswith("ITEM: NUMBER OF ATOMS"):
                    n_atoms = int(fh.readline())
                    if self.expected is not None and n_atoms != self.expected:
                        raise LatError(
                            f"{self.path}: frame {index + 1} has {n_atoms} atoms but the "
                            f"topology has {self.expected} (atom count changed across "
                            f"frames, e.g. reactive insertion/deletion); restrict the "
                            f"analysis with -frame_range/-stride to a consistent window"
                        )
                elif line.startswith("ITEM: BOX BOUNDS"):
                    triclinic = "xy" in line.split()
                    bounds = []
                    for _ in range(3):
                        bounds.extend(float(x) for x in fh.readline().split())
                    box = Box.from_lammps_dump(bounds, triclinic)
                elif line.startswith("ITEM: ATOMS"):
                    names = line.split()[2:]
                    raw = [fh.readline().split() for _ in range(n_atoms)]
                    data = _columns(names, raw)
                    index += 1
                    yield Frame(index=index, box=box, n_atoms=n_atoms, data=data)
        if index < 0:
            raise LatError(f"{self.path}: no frames found")

"""IO backends. LAMMPS data/dump are first-class; GROMACS is internal (regression)."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..errors import LatError
from ..frame import Frame
from . import gmx, lammps_data, lammps_dump

__all__ = ["gmx", "lammps_data", "lammps_dump", "detect_kind", "open_system"]


def detect_kind(path) -> str:
    p = str(path).lower()
    if p.endswith(".lammpstrj") or p.endswith(".dump"):
        return "dump"
    s = Path(p).suffix
    if s in (".data", ".lmp"):
        return "data"
    if s in (".tpr", ".gro", ".xtc", ".trr", ".pdb"):
        return s[1:]
    raise LatError(f"cannot determine format of {path!r}; expected LAMMPS data/dump or GROMACS")


class SingleFrameSource:
    n_frames_known = True

    def __init__(self, topo):
        self.topo = topo
        self.n_frames = 1

    def frames(self, indices):
        pos = self.topo.positions
        if pos is None:
            raise LatError("topology carries no coordinates and no trajectory was given (-f)")
        data = {"x": pos[:, 0], "y": pos[:, 1], "z": pos[:, 2]}
        if self.topo.velocities is not None:
            v = self.topo.velocities
            data.update({"vx": v[:, 0], "vy": v[:, 1], "vz": v[:, 2]})
        yield Frame(index=0, box=self.topo.box, n_atoms=self.topo.n_atoms, data=data)


class DumpFrameSource:
    def __init__(self, path, n_atoms):
        self.path = str(path)
        self.n_atoms = n_atoms
        self._n_frames = None

    @property
    def n_frames_known(self):
        return self._n_frames is not None

    @property
    def n_frames(self):
        """Lazily count frames: avoids a full multi-GB scan for bounded reads."""
        if self._n_frames is None:
            self._n_frames = _count_frames(self.path)
        return self._n_frames

    def frames(self, indices):
        wanted = sorted(set(int(i) for i in indices))
        if not wanted:
            return
        last = wanted[-1]
        ws = set(wanted)
        # frames stream in order: stop at the last requested frame instead of
        # reading a multi-GB trajectory to EOF
        for frame in lammps_dump.DumpReader(self.path, self.n_atoms):
            if frame.index in ws:
                yield frame
            if frame.index >= last:
                break


class MdaFrameSource:
    n_frames_known = True

    def __init__(self, universe: "gmx.GmxUniverse"):
        self.u = universe
        self.n_frames = max(universe.n_frames, 1)

    def frames(self, indices):
        yield from self.u.frames(indices)


def _count_frames(path) -> int:
    """Count frames with chunked binary reads (much faster than per-line)."""
    n = 0
    prev = b""
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            data = prev + chunk
            n += data.count(b"ITEM: TIMESTEP")
            prev = data[-16:]
    if n == 0:
        raise LatError(f"{path}: no frames found")
    return n


def open_system(topo_path, traj_path, *, style: str | None = None, cfg=None):
    """Return (Topology, FrameSource)."""
    if topo_path is None:
        raise LatError("a topology file is required (-s/--topo)")

    k_topo = detect_kind(topo_path)
    k_traj = detect_kind(traj_path) if traj_path is not None else None

    if k_topo == "data":
        topo = lammps_data.parse(topo_path, style=style, cfg=cfg)
        if k_traj in (None,):
            return topo, SingleFrameSource(topo)
        if k_traj == "dump":
            src = DumpFrameSource(traj_path, topo.n_atoms)
            return topo, src
        if k_traj in ("xtc", "trr", "gro"):
            uni = gmx.GmxUniverse(topo_path, traj_path)
            topo.n_frames = max(uni.n_frames, 1)
            return topo, MdaFrameSource(uni)
        raise LatError(f"unsupported trajectory format {k_traj!r} with a LAMMPS data topology")

    if k_topo in ("tpr", "gro", "pdb"):
        uni = gmx.GmxUniverse(topo_path, traj_path)
        topo = uni.topology(cfg=cfg)
        topo.n_frames = max(uni.n_frames, 1)
        return topo, MdaFrameSource(uni)

    raise LatError(f"unsupported topology format {k_topo!r}")

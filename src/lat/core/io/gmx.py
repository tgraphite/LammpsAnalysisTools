"""Internal GROMACS reader via MDAnalysis (regression only; not user-facing)."""

from __future__ import annotations

import re

import numpy as np

from ..errors import LatError
from ..frame import Frame
from ..pbc import Box
from ..topology import Topology, element_mass_table, infer_type_elements


def _symbols() -> set[str]:
    from ase.data import chemical_symbols

    return set(chemical_symbols[1:])


def guess_element(name: str) -> str | None:
    s = re.sub(r"[^A-Za-z]", "", str(name))
    if not s:
        return None
    syms = _symbols()
    if len(s) >= 2:
        cand = s[:2].capitalize()
        if cand in syms:
            return cand
    cand = s[:1].upper()
    return cand if cand in syms else None


class GmxUniverse:
    def __init__(self, topo_path, traj_path=None):
        try:
            import MDAnalysis as mda
        except Exception as exc:  # pragma: no cover
            raise LatError(f"MDAnalysis required for GROMACS input: {exc}")
        try:
            self.u = (mda.Universe(str(topo_path), str(traj_path)) if traj_path
                      else mda.Universe(str(topo_path)))
            self.n_frames = len(self.u.trajectory)
        except Exception as exc:
            raise LatError(
                f"cannot read GROMACS input ({topo_path}"
                + (f", {traj_path}" if traj_path else "")
                + f"): {exc}. The file may be truncated/corrupt or from a newer "
                "GROMACS version than the installed reader supports."
            )
        try:  # .tpr boxes exist only after a frame is loaded
            self.u.trajectory[0]
        except Exception as exc:
            raise LatError(f"cannot read the first frame of {traj_path}): {exc}")

    def topology(self, cfg=None) -> Topology:
        u = self.u
        ag = u.atoms
        if ag.n_atoms == 0:
            raise LatError("GROMACS universe has no atoms")

        elements = []
        try:
            elements = [str(e).capitalize() for e in ag.elements]
        except Exception:
            elements = []
        if not elements or any(not e for e in elements):
            elements = [guess_element(n) for n in ag.names]

        syms = _symbols()
        valid = [e if (e and e in syms) else None for e in elements]
        if any(e is None for e in valid):
            # fall back to mass inference per atom (tpr carries masses)
            mass_table = element_mass_table(cfg) if cfg is not None else {}
            try:
                raw_masses = np.array(ag.masses, dtype=float)
            except Exception:
                raw_masses = None
            if raw_masses is None or not np.isfinite(raw_masses).any() or not raw_masses.any():
                raise LatError(
                    "cannot determine elements: GROMACS input lacks element names and masses"
                )
            types, type_element, problems = _types_from_masses(raw_masses, mass_table, cfg)
        else:
            from ase.data import chemical_symbols

            znum = {s: i for i, s in enumerate(chemical_symbols)}
            order = sorted(set(valid), key=lambda s: znum.get(s, 999))
            zmap = {el: i + 1 for i, el in enumerate(order)}
            types = np.array([zmap[e] for e in valid], dtype=int)
            type_element = {zmap[e]: e for e in order}
            problems = {}

        box = _box_from_dimensions(getattr(u, "dimensions", None))
        positions = np.array(ag.positions, dtype=float)
        try:
            masses = np.array(ag.masses, dtype=float)
        except Exception:
            masses = None
        try:
            charges = np.array(ag.charges, dtype=float)
            if charges.size and not np.isfinite(charges).any():
                charges = None
        except Exception:
            charges = None

        mols = np.array(ag.resindices, dtype=int) + 1
        resnames = np.array([str(r) for r in ag.resnames])
        names = np.array([str(n) for n in ag.names])
        resids = np.array(ag.resids, dtype=int)

        bonds = None
        try:
            if hasattr(u, "bonds") and len(u.bonds) > 0:
                bonds = np.array([b.indices for b in u.bonds], dtype=int)
        except Exception:
            bonds = None

        # per-type mass from the first atom of each type (fall back to elements)
        mass_table = element_mass_table(cfg) if cfg is not None else {}
        type_mass = {}
        for t in np.unique(types):
            if masses is not None:
                sel = masses[types == t]
                if len(sel) and np.isfinite(sel[0]) and sel[0] > 0:
                    type_mass[int(t)] = float(sel[0])
                    continue
            el = type_element.get(int(t))
            if el and el in mass_table:
                type_mass[int(t)] = float(mass_table[el])

        velocities = None
        try:
            v = np.array(ag.velocities, dtype=float)
            if v.size and np.isfinite(v).any():
                velocities = v
        except Exception:
            velocities = None

        return Topology(
            n_atoms=ag.n_atoms,
            atom_id=np.arange(1, ag.n_atoms + 1),
            types=types,
            box=box,
            positions=positions,
            velocities=velocities,
            masses=masses,
            charges=charges,
            mols=mols,
            resnames=resnames,
            names=names,
            resids=resids,
            bonds=bonds,
            type_mass=type_mass,
            type_element=type_element,
            element_problems=problems,
            n_frames=max(self.n_frames, 1),
            source=str(self.u.filename or "gromacs"),
            meta={"format": "gromacs"},
        )

    def frames(self, indices):
        """Yield frames; also provide continuous (no-jump) coordinates as xu/yu/zu.

        xtc/trr do not store image flags, so GROMACS tools remove PBC jumps
        internally (msd/transport rely on this). We reconstruct continuous
        coordinates frame-to-frame with the minimum-image displacement, which
        is valid as long as the per-step motion is below half the box.
        """
        u = self.u
        prev = None
        unw = None
        idxs = [int(i) for i in indices]
        # xtc random access re-scans the file from the start; iterate in order
        # with a single slice so reading N selected frames is O(N), not O(N^2).
        if len(idxs) >= 2 and all(idxs[k + 1] - idxs[k] == idxs[1] - idxs[0]
                                  for k in range(len(idxs) - 1)) and idxs[1] > idxs[0]:
            frames_iter = u.trajectory[idxs[0]:idxs[-1] + 1:idxs[1] - idxs[0]]
        else:
            frames_iter = (u.trajectory[i] for i in idxs)
        for ts in frames_iter:
            idx = ts.frame
            dims = getattr(u, "dimensions", None)
            box = _box_from_dimensions(dims)
            pos = np.array(u.atoms.positions)
            if prev is None:
                unw = pos.copy()
            else:
                unw = unw + box.minimum_image(pos - prev)
            prev = pos
            data = {
                "x": pos[:, 0], "y": pos[:, 1], "z": pos[:, 2],
                "xu": unw[:, 0], "yu": unw[:, 1], "zu": unw[:, 2],
            }
            try:
                v = np.array(u.atoms.velocities)
                if v.size and np.isfinite(v).any():
                    data["vx"], data["vy"], data["vz"] = v[:, 0], v[:, 1], v[:, 2]
            except Exception:
                pass
            yield Frame(index=int(idx), box=box, n_atoms=u.atoms.n_atoms, data=data)


def _box_from_dimensions(dims):
    if dims is None or not np.all(np.isfinite(dims)) or np.allclose(dims[:3], 0):
        raise LatError("GROMACS universe has no box dimensions; load a frame first")
    lx, ly, lz, alpha, beta, gamma = [float(x) for x in dims]
    return Box.from_lengths_angles([lx, ly, lz], [alpha, beta, gamma])


def _types_from_masses(masses, mass_table, cfg):
    tol = float(cfg.get("element_match_tolerance", default=0.002)) if cfg is not None else 0.002
    type_mass, type_ids, problems = {}, np.zeros(len(masses), dtype=int), {}
    unique = {}
    for m in masses:
        key = round(float(m), 3)
        unique.setdefault(key, len(unique) + 1)
    for key, tid in unique.items():
        type_mass[tid] = key
    type_element, problems = infer_type_elements(type_mass, mass_table, tol)
    for i, m in enumerate(masses):
        type_ids[i] = unique[round(float(m), 3)]
    return type_ids, type_element, problems

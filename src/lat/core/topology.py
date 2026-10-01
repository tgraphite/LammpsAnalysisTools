"""Topology: types / elements / masses / charges / molecules / bonds.

No molecule inference (M3): `mol` or explicit `Bonds` come from the input.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .errors import LatError
from .pbc import Box


def element_mass_table(cfg) -> dict:
    """Full IUPAC masses from ase.data, overridden by config entries."""
    table: dict[str, float] = {}
    try:
        from ase.data import atomic_masses, chemical_symbols

        for z in range(1, len(chemical_symbols)):
            table[chemical_symbols[z]] = float(atomic_masses[z])
    except Exception:  # ase optional at runtime
        pass
    cfg_mass = cfg.get("element_mass", default={}) or {}
    table.update({k: float(v) for k, v in cfg_mass.items()})
    return table


def infer_type_elements(type_mass: dict, mass_table: dict, tol: float):
    """Return (type_element, problems). Ambiguous/unknown types get None."""
    result: dict[int, str | None] = {}
    problems: dict[int, str] = {}
    for t, m in sorted(type_mass.items()):
        matches = [
            el for el, mm in mass_table.items() if mm > 0 and abs(m - mm) / mm <= tol
        ]
        if len(matches) == 1:
            result[t] = matches[0]
        elif not matches:
            result[t] = None
            problems[t] = f"mass {m:.4f} matches no element within +/-{tol*100:.1f}%"
        else:
            result[t] = None
            problems[t] = f"mass {m:.4f} matches multiple elements {matches}"
    return result, problems


@dataclass
class Topology:
    n_atoms: int
    atom_id: np.ndarray
    types: np.ndarray
    box: Box
    positions: np.ndarray | None = None
    velocities: np.ndarray | None = None
    masses: np.ndarray | None = None
    charges: np.ndarray | None = None
    mols: np.ndarray | None = None
    resnames: np.ndarray | None = None
    names: np.ndarray | None = None
    resids: np.ndarray | None = None
    bonds: np.ndarray | None = None
    type_mass: dict = field(default_factory=dict)
    type_element: dict = field(default_factory=dict)
    element_problems: dict = field(default_factory=dict)
    n_frames: int = 1
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        self.n_atoms = int(self.n_atoms)
        self.atom_id = np.asarray(self.atom_id)
        self.types = np.asarray(self.types, dtype=int)

    # ---- helpers ------------------------------------------------------
    @property
    def unique_types(self) -> np.ndarray:
        return np.unique(self.types)

    def mask_of_types(self, ts) -> np.ndarray:
        ts = np.atleast_1d(np.asarray(ts, dtype=int))
        return np.isin(self.types, ts)

    def elements_of_types(self, ts) -> list[str]:
        out = []
        for t in np.atleast_1d(ts):
            el = self.type_element.get(int(t))
            if el is None:
                raise LatError(
                    f"type {int(t)} has no element mapping"
                    + (f" ({self.element_problems.get(int(t), '')})" if self.element_problems else "")
                    + "; provide -types or a complete -style/mass table"
                )
            out.append(el)
        return out

    def atom_mass(self) -> np.ndarray:
        if self.masses is not None:
            return self.masses
        if not self.type_mass:
            raise LatError(
                "no per-type masses available (data file has no Masses section); "
                "mass weighting is impossible — use -quantity number for density, "
                "or supply masses"
            )
        return np.array([self.type_mass.get(int(t), 0.0) for t in self.types])

    def require_elements(self) -> None:
        bad = [t for t, e in self.type_element.items() if e is None]
        if bad:
            detail = "; ".join(f"type {t}: {self.element_problems.get(t, 'unknown')}" for t in bad)
            raise LatError(f"cannot map type->element for {bad}: {detail}")


def reorder_by_id(arrays: dict, ids: np.ndarray):
    order = np.argsort(ids, kind="stable")
    return order, {k: (v[order] if v is not None else None) for k, v in arrays.items()}

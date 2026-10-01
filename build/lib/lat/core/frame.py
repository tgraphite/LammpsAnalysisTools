"""A single trajectory frame with lazy coordinate access."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .errors import LatError
from .pbc import Box


@dataclass
class Frame:
    index: int
    box: Box
    n_atoms: int
    data: dict = field(default_factory=dict)

    def has(self, *names) -> bool:
        return all(n in self.data for n in names)

    def vector(self, name) -> np.ndarray:
        if name not in self.data:
            raise LatError(f"trajectory frame {self.index} lacks column {name!r}")
        return np.asarray(self.data[name], float)

    def _axes(self, pos_names, scaled_names, base_names):
        if self.has(*pos_names):
            return np.column_stack([self.vector(n) for n in pos_names])
        if self.has(*scaled_names):
            L = self.box.lengths
            return np.column_stack([self.vector(n) for n in scaled_names]) * L + self.box.origin
        if self.has(*base_names):
            x = np.column_stack([self.vector(n) for n in base_names])
            if self.has("ix", "iy", "iz"):
                img = np.column_stack([self.vector(n) for n in ("ix", "iy", "iz")])
                x = x + img * self.box.lengths
            return x
        raise LatError(
            f"trajectory frame {self.index} lacks coordinates "
            f"(need x/y/z, xs/ys/zs or xu/yu/zu)"
        )

    @property
    def positions(self) -> np.ndarray:
        """Wrapped (or as-given) cartesian positions, shape (N,3)."""
        if self.has("x", "y", "z") or self.has("xs", "ys", "zs"):
            return self._axes(("x", "y", "z"), ("xs", "ys", "zs"), ("x", "y", "z"))
        return self.unwrapped

    @property
    def unwrapped(self) -> np.ndarray:
        """Best-effort continuous coordinates for MSD/transport."""
        if self.has("xu", "yu", "zu"):
            return np.column_stack([self.vector(n) for n in ("xu", "yu", "zu")])
        if self.has("xsu", "ysu", "zsu"):
            L = self.box.lengths
            return np.column_stack([self.vector(n) for n in ("xsu", "ysu", "zsu")]) * L + self.box.origin
        if self.has("x", "y", "z", "ix", "iy", "iz"):
            x = np.column_stack([self.vector(n) for n in ("x", "y", "z")])
            img = np.column_stack([self.vector(n) for n in ("ix", "iy", "iz")])
            return x + img * self.box.lengths
        return self.positions

    @property
    def velocities(self) -> np.ndarray | None:
        if self.has("vx", "vy", "vz"):
            return np.column_stack([self.vector(n) for n in ("vx", "vy", "vz")])
        return None

    @property
    def types(self) -> np.ndarray | None:
        return np.asarray(self.data["type"], dtype=int) if "type" in self.data else None

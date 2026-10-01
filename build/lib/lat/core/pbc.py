"""Periodic boundary conditions: orthorhombic and triclinic minimum image.

Convention: H is a 3x3 matrix whose *columns* are the box vectors a, b, c.
Cartesian row vectors r = s @ H.T ; fractional s = r @ inv(H).T.
Triclinic MIC rounds fractional coordinates to the nearest image (standard
minimum-image convention).
"""

from __future__ import annotations

import numpy as np

from .errors import LatError

_ORTHO_TOL = 1e-8


class Box:
    def __init__(self, H: np.ndarray, origin=None):
        self.H = np.asarray(H, dtype=float).reshape(3, 3)
        if abs(np.linalg.det(self.H)) < 1e-12:
            raise LatError("degenerate box (zero volume)")
        self.Hinv = np.linalg.inv(self.H)
        self.origin = np.zeros(3) if origin is None else np.asarray(origin, float)

    # ---- constructors -------------------------------------------------
    @classmethod
    def orthorhombic(cls, lengths, origin=None):
        L = np.asarray(lengths, float)
        return cls(np.diag(L), origin=origin)

    @classmethod
    def from_lammps_data(cls, lo, hi, tilt=None):
        lo = np.asarray(lo, float)
        hi = np.asarray(hi, float)
        t = np.zeros(3) if tilt is None else np.asarray(tilt, float)
        H = np.array(
            [
                [hi[0] - lo[0], t[0], t[1]],
                [0.0, hi[1] - lo[1], t[2]],
                [0.0, 0.0, hi[2] - lo[2]],
            ]
        )
        return cls(H, origin=lo)

    @classmethod
    def from_lammps_dump(cls, bounds, triclinic: bool):
        """bounds: (xlo_b, xhi_b, ylo_b, yhi_b, zlo_b, zhi_b[, xy, xz, yz])."""
        bounds = [float(b) for b in bounds]
        if triclinic:
            xlo_b, xhi_b, ylo_b, yhi_b, zlo_b, zhi_b, xy, xz, yz = bounds
            xlo = xlo_b - min(0.0, xy, xz, xy + xz)
            xhi = xhi_b - max(0.0, xy, xz, xy + xz)
            ylo = ylo_b - min(0.0, yz)
            yhi = yhi_b - max(0.0, yz)
            return cls.from_lammps_data([xlo, ylo, zlo_b], [xhi, yhi, zhi_b], [xy, xz, yz])
        xlo, xhi, ylo, yhi, zlo, zhi = bounds
        return cls.from_lammps_data([xlo, ylo, zlo], [xhi, yhi, zhi])

    @classmethod
    def from_lengths_angles(cls, lengths, angles, origin=None):
        """GROMACS-style box: lengths (Å) and angles (deg)."""
        lx, ly, lz = [float(v) for v in lengths]
        alpha, beta, gamma = np.deg2rad([float(v) for v in angles])
        sg = np.sin(gamma)
        if abs(sg) < 1e-12:
            raise LatError("invalid box angle gamma=0")
        cx = lz * np.cos(beta)
        cy = lz * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sg
        cz = np.sqrt(max(lz * lz - cx * cx - cy * cy, 0.0))
        H = np.array(
            [
                [lx, ly * np.cos(gamma), cx],
                [0.0, ly * sg, cy],
                [0.0, 0.0, cz],
            ]
        )
        return cls(H, origin=origin)

    # ---- properties ---------------------------------------------------
    @property
    def lengths(self) -> np.ndarray:
        return np.linalg.norm(self.H, axis=0)

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.H)))

    @property
    def tilt(self) -> np.ndarray:
        return np.array([self.H[0, 1], self.H[0, 2], self.H[1, 2]])

    @property
    def is_orthorhombic(self) -> bool:
        return bool(np.all(np.abs(self.tilt) < _ORTHO_TOL))

    @property
    def lo(self) -> np.ndarray:
        return self.origin.copy()

    @property
    def hi(self) -> np.ndarray:
        return self.origin + self.lengths

    # ---- transforms ---------------------------------------------------
    def cart2frac(self, r: np.ndarray) -> np.ndarray:
        return (np.asarray(r, float) - self.origin) @ self.Hinv.T

    def frac2cart(self, s: np.ndarray) -> np.ndarray:
        return np.asarray(s, float) @ self.H.T + self.origin

    def minimum_image(self, delta: np.ndarray) -> np.ndarray:
        """Vectorized MIC for an array of row-vector displacements."""
        d = np.asarray(delta, float)
        s = d @ self.Hinv.T
        s -= np.round(s)
        return s @ self.H.T

    def wrap(self, r: np.ndarray) -> np.ndarray:
        s = self.cart2frac(r)
        s -= np.floor(s)
        # keep strictly inside the cell so KDTree(boxsize=...) never sees the
        # upper boundary (floating-point frac can land on 1.0)
        s = np.clip(s, 0.0, 1.0 - 1e-9)
        return self.frac2cart(s)

    def image_shifts(self, cutoff: float) -> np.ndarray:
        """Integer image shifts needed so no pair within `cutoff` is missed."""
        lengths = self.lengths
        # perpendicular height of each face: V / area(opposite face)
        H = self.H
        heights = np.empty(3)
        faces = [(1, 2), (0, 2), (0, 1)]
        for k, (i, j) in enumerate(faces):
            area = np.linalg.norm(np.cross(H[:, i], H[:, j]))
            heights[k] = self.volume / area if area > 0 else lengths[k]
        n = np.maximum(np.ceil(cutoff / heights).astype(int), 1)
        rng = [np.arange(-n[k], n[k] + 1) for k in range(3)]
        grids = np.meshgrid(*rng, indexing="ij")
        return np.stack([g.ravel() for g in grids], axis=1).astype(float)

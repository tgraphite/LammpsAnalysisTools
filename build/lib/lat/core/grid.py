"""Automatic 1D/3D grids: bin edges, accumulation, optional smoothing."""

from __future__ import annotations

import numpy as np

from .errors import LatError
from .pbc import Box


def axis_bins(length: float, nbins: int, max_bins: int) -> int:
    if nbins < 1:
        raise LatError(f"nbins must be >= 1, got {nbins}")
    if length <= 0:
        raise LatError(f"cannot grid a non-positive axis length ({length})")
    if nbins > max_bins:
        raise LatError(f"nbins={nbins} exceeds max_bins_axis={max_bins}")
    return int(nbins)


def rdf_bins(nbins: int, cutoff: float):
    edges = np.linspace(0.0, float(cutoff), int(nbins) + 1)
    r = 0.5 * (edges[:-1] + edges[1:])
    return edges, r


def profile_bins(box: Box, axis: str, nbins: int):
    idx = {"x": 0, "y": 1, "z": 2}.get(axis)
    if idx is None:
        raise LatError(f"axis must be x/y/z, got {axis!r}")
    lo = box.lo[idx]
    hi = box.hi[idx]
    edges = np.linspace(lo, hi, int(nbins) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers


def smooth(y: np.ndarray, method: str, *, window: int, order: int, sigma: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import savgol_filter

    y = np.asarray(y, float)
    if method == "none":
        return y.copy()
    if method == "savgol":
        w = int(window)
        if w % 2 == 0:
            w += 1
        if w > len(y):
            raise LatError(f"savgol window {window} larger than data length {len(y)}")
        return savgol_filter(y, w, int(order))
    if method == "gaussian":
        return gaussian_filter1d(y, float(sigma), mode="nearest")
    raise LatError(f"unknown smoothing method {method!r}")

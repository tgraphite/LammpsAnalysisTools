"""Automatic 1D/3D grids: bin edges, accumulation, optional smoothing."""

from __future__ import annotations

import numpy as np

from .errors import LatError


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


def fractional_bins(nbins: int):
    """Unit-interval edges and centers for box-relative (fractional) binning."""
    edges = np.linspace(0.0, 1.0, int(nbins) + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers


def box_change_warning(lengths, axes, tol: float = 0.02):
    """WARN text if any of `axes` deviates from its mean length by more than tol.

    `lengths` is (nframes, 3). Returns None when every requested axis stays
    within tolerance of its mean.
    """
    L = np.asarray(lengths, float)
    if L.size == 0:
        return None
    mean = L.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = np.abs(L - mean).max(axis=0) / mean
    bad = [(ax, float(dev[ax])) for ax in axes if mean[ax] > 0 and dev[ax] > tol]
    if not bad:
        return None
    parts = ", ".join(f"{'xyz'[ax]} ±{100 * d:.1f}%" for ax, d in bad)
    return (f"box length varies ({parts}) beyond ±{tol * 100:.0f}% of the mean; "
            f"bins are fractional and the axis is scaled to the mean box length")


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

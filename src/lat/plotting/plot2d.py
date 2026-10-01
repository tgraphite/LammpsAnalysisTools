"""Plot a lat 2D long-table (x, y, value) as a color map."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from ..core.errors import LatError  # noqa: E402
from .plot1d import read_lat_csv  # noqa: E402


def _smooth_grid(grid, method, window, order, sigma):
    if not method or method == "none":
        return grid
    finite = np.isfinite(grid)
    if not finite.any():
        return grid
    fill = grid.copy()
    fill[~finite] = grid[finite].mean()
    if method == "gaussian":
        from scipy.ndimage import gaussian_filter
        out = gaussian_filter(fill, float(sigma))
    elif method == "savgol":
        from scipy.signal import savgol_filter
        w = int(window)
        if w % 2 == 0:
            w += 1
        if w > min(grid.shape):
            raise LatError(f"savgol window {window} larger than grid size {grid.shape}")
        out = savgol_filter(fill, w, int(order), axis=0)
        out = savgol_filter(out, w, int(order), axis=1)
    else:
        raise LatError(f"unknown smoothing method {method!r}")
    out[~finite] = np.nan
    return out


def plot_grid(path, out=None, cmap="viridis", xlabel=None, ylabel=None, title=None,
              smooth="none", window=5, order=2, sigma=1.0):
    parsed = read_lat_csv(path)
    header = parsed["header"]
    cols = parsed["columns"]
    data = parsed["data"]
    if data.shape[1] != 3:
        raise LatError(f"{path}: 2D plot needs exactly 3 columns (x, y, value), got {data.shape[1]}")

    x, y, v = data[:, 0], data[:, 1], data[:, 2]
    xu, yu = np.unique(x), np.unique(y)
    if len(xu) * len(yu) == len(x):
        xi = {val: i for i, val in enumerate(xu)}
        yi = {val: i for i, val in enumerate(yu)}
        grid = np.full((len(xu), len(yu)), np.nan)
        for a, b, c in zip(x, y, v):
            grid[xi[a], yi[b]] = c
        grid = _smooth_grid(grid, smooth, window, order, sigma)
        fig, ax = plt.subplots(figsize=(6, 5))
        mesh = ax.pcolormesh(xu, yu, grid.T, cmap=cmap, shading="auto")
    else:
        fig, ax = plt.subplots(figsize=(6, 5))
        mesh = ax.tricontourf(x, y, v, levels=20, cmap=cmap)
    ax.set_xlabel(xlabel or header.get("xlabel", cols[0]))
    ax.set_ylabel(ylabel or header.get("ylabel", cols[1]))
    ax.set_title(title or header.get("title", ""))
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.locator = MaxNLocator(nbins=8)
    cbar.update_ticks()
    cbar.set_label(cols[2])
    fig.tight_layout()
    out = Path(out) if out else Path(path).with_suffix(".png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out

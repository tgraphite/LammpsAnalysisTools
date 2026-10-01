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


def plot_grid(path, out=None, cmap="viridis", title=None):
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
        fig, ax = plt.subplots(figsize=(6, 5))
        mesh = ax.pcolormesh(xu, yu, grid.T, cmap=cmap, shading="auto")
    else:
        fig, ax = plt.subplots(figsize=(6, 5))
        mesh = ax.tricontourf(x, y, v, levels=20, cmap=cmap)
    ax.set_xlabel(f"{cols[0]} ({header.get('units', '').split(',')[0].split('=')[-1].strip()})")
    ax.set_ylabel(cols[1])
    ax.set_title(title or header.get("title", Path(path).stem))
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.locator = MaxNLocator(nbins=8)
    cbar.update_ticks()
    cbar.set_label(cols[2])
    fig.tight_layout()
    out = Path(out) if out else Path(path).with_suffix(".png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out

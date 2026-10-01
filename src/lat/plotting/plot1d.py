"""Plot y columns of a lat 1D CSV against its x column.

Reads the plotting contract from the CSV header (`title`, `xlabel`, `ylabel`,
`columns`); the y names in `columns` are used as the legend labels, so a user
can rename e.g. `g_type_1_2` to `Water-Oil` and re-run to get a clean figure.

The x axis is always the first column. `-columns` restricts which y columns are
drawn; `-bar` draws grouped bars and `-stack` draws a normalized (occupancy)
stacked bar chart.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.errors import LatError  # noqa: E402
from ..core.grid import smooth as smooth_fn  # noqa: E402


def read_lat_csv(path):
    """Return dict(header=..., columns=[...], data=np.ndarray)."""
    path = Path(path)
    if not path.is_file():
        raise LatError(f"CSV not found: {path}")
    header: dict[str, str] = {}
    columns = None
    rows = []
    in_data = False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not in_data:
                if line.strip().startswith("# --- data ---"):
                    in_data = True
                    continue
                if line.startswith("#") and ":" in line:
                    key, val = line.lstrip("#").split(":", 1)
                    header[key.strip()] = val.strip()
                    if key.strip() == "columns":
                        columns = [c.strip() for c in val.split(",")]
                continue
            if line:
                rows.append([float(x) for x in line.split(",")])
    if columns is None:
        raise LatError(f"{path}: no 'columns' header found; not a lat CSV")
    if not rows:
        raise LatError(f"{path}: no data rows")
    return {"header": header, "columns": columns, "data": np.array(rows)}


def _select_y(columns, requested):
    """Resolve requested y-column names to [(index, name), ...].

    The x column is always `columns[0]` and is always drawn, so it is never a
    selectable y column. `requested=None` means "all y columns".
    """
    if not requested:
        return [(k, columns[k]) for k in range(1, len(columns))]
    x_name = columns[0]
    lookup = {columns[k]: k for k in range(1, len(columns))}
    picked: list[tuple[int, str]] = []
    seen: set[str] = set()
    for name in requested:
        if name == x_name:
            continue  # x is always drawn; ignore it as a y selection
        if name not in lookup:
            raise LatError(
                f"column {name!r} not found; available y columns: "
                f"{', '.join(columns[1:])}")
        if name in seen:
            continue
        seen.add(name)
        picked.append((lookup[name], name))
    if not picked:
        raise LatError(
            "-columns selected no y column (the first column is always the x axis)")
    return picked


def _normalize(ys):
    """Normalize each column to its group sum (occupancy fraction)."""
    total = ys.sum(axis=0)
    return np.divide(ys, total, out=np.zeros_like(ys, float), where=total > 0)


def _bar_geometry(x, n_series):
    """Return (bar_width, offsets) for grouped bars at positions `x`."""
    x = np.asarray(x, float)
    if x.size > 1:
        dx = float(np.min(np.diff(np.sort(x))))
        if not np.isfinite(dx) or dx <= 0:
            dx = 1.0
    else:
        dx = 1.0
    group = dx * 0.8
    width = group / max(n_series, 1)
    offsets = (np.arange(n_series) - (n_series - 1) / 2.0) * width
    return width, offsets


def plot_csv(path, out=None, xlabel=None, ylabel=None, title=None,
             smooth="none", window=5, order=2, sigma=1.0,
             columns=None, mode="line"):
    parsed = read_lat_csv(path)
    header = parsed["header"]
    cols = parsed["columns"]
    data = parsed["data"]
    if data.shape[1] < 2:
        raise LatError(f"{path}: needs at least x and one y column")

    picked = _select_y(cols, columns)
    if mode not in ("line", "bar", "stack"):
        raise LatError(f"unknown plot mode {mode!r}")
    if mode != "line" and smooth not in (None, "none"):
        raise LatError("bar/stack charts are not smoothed; drop -smooth or use a line plot")

    x = data[:, 0]
    fig, ax = plt.subplots(figsize=(6, 4))

    if mode == "line":
        for k, name in picked:
            y = data[:, k]
            if smooth and smooth != "none":
                y = smooth_fn(y, smooth, window=int(window), order=int(order),
                              sigma=float(sigma))
            ax.plot(x, y, label=name)
        y_default = header.get("ylabel", "value")

    elif mode == "bar":
        width, offsets = _bar_geometry(x, len(picked))
        for (k, name), off in zip(picked, offsets):
            ax.bar(x + off, data[:, k], width=width, label=name)
        y_default = header.get("ylabel", "value")

    else:  # stack: normalize each x group to an occupancy fraction, then stack
        ys = np.vstack([data[:, k] for k, _ in picked])
        frac = _normalize(ys)
        width, _ = _bar_geometry(x, 1)
        bottom = np.zeros(len(x))
        for (k, name), fr in zip(picked, frac):
            ax.bar(x, fr, width=width, bottom=bottom, label=name)
            bottom = bottom + fr
        ax.set_ylim(0.0, 1.0)
        y_default = "fraction"

    ax.set_xlabel(xlabel or header.get("xlabel", cols[0]))
    ax.set_ylabel(ylabel or y_default)
    ax.set_title(title or header.get("title", ""))
    # `columns` doubles as the legend; always show it, without a box.
    ax.legend(frameon=False)
    fig.tight_layout()
    out = Path(out) if out else Path(path).with_suffix(".png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out

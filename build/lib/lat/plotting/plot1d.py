"""Plot every y column of a lat 1D CSV against its x column."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ..core.errors import LatError  # noqa: E402


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


def plot_csv(path, out=None, xlabel=None, ylabel=None, title=None):
    parsed = read_lat_csv(path)
    header = parsed["header"]
    cols = parsed["columns"]
    data = parsed["data"]
    if data.shape[1] < 2:
        raise LatError(f"{path}: needs at least x and one y column")

    units = dict(
        (kv.split("=", 1)[0].strip(), kv.split("=", 1)[1].strip())
        for kv in header.get("units", "").split(",")
        if "=" in kv
    )
    fig, ax = plt.subplots(figsize=(6, 4))
    for k in range(1, data.shape[1]):
        name = cols[k]
        label = name if name not in units else f"{name} ({units[name]})"
        ax.plot(data[:, 0], data[:, k], label=label)
    xname = xlabel or cols[0]
    ax.set_xlabel(f"{xname} ({units.get(cols[0], '?')})")
    ax.set_ylabel(ylabel or (cols[1] if len(cols) == 2 else "value"))
    ax.set_title(title or header.get("title", Path(path).stem))
    if data.shape[1] > 2:
        ax.legend()
    fig.tight_layout()
    out = Path(out) if out else Path(path).with_suffix(".png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out

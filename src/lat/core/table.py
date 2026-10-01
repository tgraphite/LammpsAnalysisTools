"""Table contract: modules return Tables; output.py serializes them."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .errors import LatError


@dataclass(frozen=True)
class Field:
    name: str
    unit: str

    def label(self) -> str:
        return f"{self.name} ({self.unit})"


@dataclass
class Table:
    """One CSV file: a single table with a shared axis/axes.

    kind="curve":  columns = [x, y1, y2, ...]
    kind="grid":   columns = [x, y, value]   (long form 2D)
    """

    title: str
    columns: list[Field]
    data: np.ndarray
    tag: str
    kind: str = "curve"
    mode: str | None = None
    xlabel: str | None = None
    ylabel: str | None = None
    definitions: dict = field(default_factory=dict)
    scalars: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)
    selection: dict = field(default_factory=dict)
    frame_range: str = "all"

    def __post_init__(self) -> None:
        self.data = np.asarray(self.data)
        if self.data.ndim != 2 or self.data.shape[1] != len(self.columns):
            raise LatError(
                f"Table {self.tag!r}: data shape {self.data.shape} inconsistent with "
                f"{len(self.columns)} columns"
            )
        if self.kind not in ("curve", "grid"):
            raise LatError(f"Table {self.tag!r}: unknown kind {self.kind!r}")

    @property
    def x(self) -> Field:
        return self.columns[0]

    @property
    def y(self) -> list[Field]:
        return list(self.columns[1:])

    @property
    def x_label(self) -> str:
        """Axis label for x: explicit, else the x column's `name (unit)`."""
        return self.xlabel or self.columns[0].label()

    @property
    def y_label(self) -> str:
        """Axis label for y: explicit, else derived from the y columns."""
        if self.ylabel:
            return self.ylabel
        if self.kind == "grid":
            return self.columns[1].label() if len(self.columns) > 1 else "value"
        ys = self.columns[1:]
        units = {c.unit for c in ys}
        return f"value ({units.pop()})" if len(units) == 1 else "value"

    @property
    def n_rows(self) -> int:
        return int(self.data.shape[0])

    def add_warning(self, msg: str) -> None:
        self.warnings.append(msg)

    def filename(self, basename: str, ext: str = "csv", suffix: str = "out") -> str:
        parts = [basename, self.tag]
        if self.mode:
            parts.append(self.mode)
        parts.append(suffix)
        return ".".join(parts) + f".{ext}"

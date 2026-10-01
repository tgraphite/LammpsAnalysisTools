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

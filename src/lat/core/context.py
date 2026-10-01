"""Context: topology + streaming frame accessor + units + frame bookkeeping."""

from __future__ import annotations

import time

import numpy as np

from .errors import LatError, get_logger
from .selection import Selection, parse as parse_selection
from .topology import Topology


def resolve_frames(spec: str | None, n_frames: int | None, stride: int = 1) -> list[int]:
    """frame_range syntax start:stop:stride (stop exclusive); default all frames.

    `n_frames` may be None when the total is not yet known (lazy dump counting):
    an explicit start:stop is then honoured without an upper-bound check, so a
    bounded read never has to scan a multi-GB trajectory just to count frames.
    """
    need_total = spec in (None, "", "all") or (str(spec).count(":") == 0)
    if need_total and n_frames is None:
        raise LatError(
            "total frame count is unknown; give an explicit frame_range start:stop "
            "or let the reader count frames"
        )
    if spec in (None, "", "all"):
        start, stop, st = 0, n_frames, int(stride)
    else:
        parts = str(spec).split(":")
        if len(parts) == 1:
            start, stop, st = int(parts[0]), n_frames, int(stride)
        elif len(parts) == 2:
            start, stop, st = int(parts[0]), int(parts[1]), int(stride)
        elif len(parts) == 3:
            start, stop, st = int(parts[0]), int(parts[1]), int(parts[2])
        else:
            raise LatError(f"bad frame_range {spec!r}; expected start:stop[:stride]")
    if st < 1:
        raise LatError(f"frame stride must be >= 1, got {st}")
    if start < 0 or start >= stop:
        raise LatError(f"frame_range {spec!r} is empty (start must be >= 0 and < stop)")
    if n_frames is not None and stop > n_frames:
        raise LatError(
            f"frame_range {spec!r} out of bounds for {n_frames} frames (0..{n_frames - 1})"
        )
    return list(range(start, stop, st))


class Context:
    def __init__(self, topo: Topology, source, *, cfg=None, unit: str = "real",
                 command: str = "", debug: bool = False, ndx: dict | None = None):
        self.topo = topo
        self.source = source
        self.cfg = cfg
        self.unit = unit
        self.command = command
        self.debug = debug
        self.ndx = ndx or {}
        self._n_frames = None
        self.frame_interval = cfg.get("meta.frame_interval") if cfg is not None else None

    @property
    def n_frames(self) -> int:
        """Total frames; counted lazily so explicit frame_range stays cheap."""
        if self._n_frames is None:
            self._n_frames = int(self.source.n_frames)
        return self._n_frames

    @property
    def n_frames_known(self) -> bool:
        """True when the total is already known without scanning the file."""
        if self._n_frames is not None:
            return True
        return bool(getattr(self.source, "n_frames_known", False))

    def selection(self, literal: str | None) -> Selection:
        return parse_selection(self.topo, literal, self.ndx)

    def element_map(self) -> dict:
        """type->element, falling back to an `element` column in the trajectory.

        Some topology writers (e.g. OVITO) omit the Masses section, so element
        identity is unavailable from the data file. If the dump carries an
        `element` column, recover the mapping from its first frame.
        """
        m = {int(k): v for k, v in (self.topo.type_element or {}).items()}
        types = self.topo.types
        n_types = len(set(int(t) for t in types)) if types is not None else 0
        if n_types and len(m) >= n_types:
            return m
        try:
            for fr in self.source.frames([0]):
                if "type" in fr.data and "element" in fr.data:
                    for t, e in zip(fr.data["type"], fr.data["element"]):
                        m.setdefault(int(t), str(e))
                break
        except Exception:
            pass
        return m

    def frames(self, indices):
        """Stream frames with periodic progress on stderr (always shown)."""
        indices = list(indices)
        total = len(indices)
        if total == 0:
            return
        log = get_logger()
        src = getattr(self.source, "path", type(self.source).__name__)
        log.notice(f"reading {total} frame(s) from {src}")
        t0 = time.time()
        n = 0
        last = t0
        for fr in self.source.frames(indices):
            n += 1
            now = time.time()
            if n == 1 or n == total or now - last >= 2.0:
                log.notice(f"  processed frame {n}/{total} ({now - t0:.1f}s)")
                last = now
            yield fr
        log.notice(f"  done: {n}/{total} frames in {time.time() - t0:.1f}s")

    def resolve_frames(self, spec: str | None, stride: int = 1) -> list[int]:
        # force a count only when the spec needs the total, or when the source
        # already knows it for free (so out-of-range explicit specs still error)
        need_total = spec in (None, "", "all") or str(spec).count(":") == 0
        nf = self.n_frames if (need_total or self.n_frames_known) else None
        return resolve_frames(spec, nf, stride)

    def effective_dt(self, indices):
        """Time between two *selected* frames (ps), accounting for stride."""
        if not self.frame_interval or len(indices) < 2:
            return self.frame_interval
        step = abs(int(indices[1]) - int(indices[0]))
        return float(self.frame_interval) * max(step, 1)

    def frame_range_str(self, spec: str | None, indices) -> str:
        if spec in (None, "", "all"):
            if len(indices) == self.n_frames:
                return "0:%d:1" % self.n_frames
        return str(spec) if spec not in (None, "") else "all"

    def box_from_frame(self, frame):
        return frame.box

    def average_box(self, indices) -> np.ndarray:
        """Mean box lengths over the selected frames (for density normalization)."""
        acc = np.zeros(3)
        n = 0
        for fr in self.frames(indices):
            acc += fr.box.lengths
            n += 1
        if n == 0:
            raise LatError("no frames to average box over")
        return acc / n

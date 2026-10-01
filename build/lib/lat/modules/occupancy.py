"""lat occupancy — zone populations, crossings and first-passage."""

from __future__ import annotations

from collections import Counter

import numpy as np

from ..core.errors import LatError
from ..core.table import Field, Table
from .base import Module, Param
from . import register

_AXES = {"x": 0, "y": 1, "z": 2}


def parse_zone(spec: str):
    """Parse 'z:0-10,10-20 x:0-5' -> (axes, per_axis_ranges, n_zones)."""
    if not spec:
        raise LatError("occupancy: -zone is required")
    per_axis: dict[int, list[tuple[float, float]]] = {}
    for token in spec.split():
        if ":" not in token:
            raise LatError(f"occupancy: bad zone token {token!r}; expected axis:min-max[,min-max]")
        axis, ranges = token.split(":", 1)
        if axis not in _AXES:
            raise LatError(f"occupancy: unknown axis {axis!r}")
        k = _AXES[axis]
        blocks = []
        for r in ranges.split(","):
            if "-" not in r[1:]:
                raise LatError(f"occupancy: bad range {r!r}")
            a, b = r.split("-", 1)
            a, b = float(a), float(b)
            if b <= a:
                raise LatError(f"occupancy: empty range [{a},{b})")
            blocks.append((a, b))
        blocks.sort()
        for (a1, b1), (a2, b2) in zip(blocks, blocks[1:]):
            if a2 < b1:
                raise LatError(f"occupancy: overlapping ranges on axis {axis}")
        per_axis[k] = blocks
    counts = [len(v) for v in per_axis.values()]
    n_zones = 1
    for c in counts:
        n_zones *= c
    axes = sorted(per_axis)
    return axes, per_axis, n_zones


def _labels(pos, axes, per_axis, n_zones):
    n = len(pos)
    label = np.zeros(n, dtype=int)
    zone = np.ones(n, dtype=bool)
    zone_id = np.ones(n, dtype=int)
    zone_index = 0
    # iterate cartesian product of per-axis block choices
    choices = [per_axis[k] for k in axes]
    import itertools

    for combo in itertools.product(*choices):
        inside = np.ones(n, dtype=bool)
        for k, (a, b) in zip(axes, combo):
            c = pos[:, k]
            inside &= (c >= a) & (c < b)
        zone_index += 1
        label[inside] = zone_index
    return label


@register
class Occupancy(Module):
    name = "occupancy"
    summary = "zone occupancy, boundary crossings and first-passage times"
    help_doc = "algorithms/occupancy.md"
    inputs = ["-s", "-sel", "-zone"]
    outputs = ["<base>.occupancy.out.csv", "<base>.occupancy.crossings.out.csv",
               "<base>.occupancy.first_passage.out.csv", "<base>.occupancy.enrichment.out.csv"]
    params = [
        Param("sel", flags=("-sel",), required=True),
        Param("zone", flags=("-zone",), config="modules.occupancy.zone", required=True),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel, zone, frame_range, stride):
        axes, per_axis, n_zones = parse_zone(zone)
        selection = ctx.selection(sel)
        if selection.n_atoms == 0:
            raise LatError("occupancy: empty selection")
        idx = selection.indices
        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        warnings: list[str] = []
        if len(frames) < 2:
            warnings.append("fewer than 2 frames; crossings/first-passage are trivial")

        occ_rows = []
        cross_rows = []
        first_passage = np.full(len(idx), -1, dtype=int)
        prev = None
        for c, fr in enumerate(ctx.frames(frames)):
            lab = _labels(fr.positions[idx], axes, per_axis, n_zones)
            counts = [int((lab == z).sum()) for z in range(1, n_zones + 1)]
            occ_rows.append([fr.index] + counts)
            if prev is not None:
                changed = np.flatnonzero(lab != prev)
                cc = Counter((int(prev[i]), int(lab[i])) for i in changed)
                for (a, b), n in sorted(cc.items()):
                    cross_rows.append((fr.index, a, b, n))
            arrived = (lab > 0) & (first_passage < 0)
            first_passage[arrived] = fr.index
            prev = lab

        tables = [Table(
            title=f"Zone occupancy {zone}",
            columns=[Field("frame", "frame")] + [Field(f"zone{z}", "1") for z in range(1, n_zones + 1)],
            data=np.array(occ_rows, float), tag="occupancy", warnings=warnings,
            meta={"module": "occupancy"},
            selection={"sel": selection.literal, "zone": zone},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]
        tables.append(Table(
            title="Zone crossings",
            columns=[Field("frame", "frame"), Field("from_zone", "1"), Field("to_zone", "1"),
                     Field("count", "1")],
            data=np.array(cross_rows, float) if cross_rows else np.empty((0, 4)),
            tag="occupancy.crossings", warnings=warnings, meta={"module": "occupancy"},
            selection={"sel": selection.literal, "zone": zone},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))
        fp = Counter(int(v) for v in first_passage if v >= 0)
        fp_rows = [[k, fp[k]] for k in sorted(fp)]
        tables.append(Table(
            title="First passage distribution",
            columns=[Field("frame", "frame"), Field("count", "1")],
            data=np.array(fp_rows, float) if fp_rows else np.empty((0, 2)),
            tag="occupancy.first_passage", warnings=warnings, meta={"module": "occupancy"},
            selection={"sel": selection.literal, "zone": zone},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))

        box = next(iter(ctx.frames(frames))).box
        total_vol = box.volume
        n_frames = len(frames)
        enr_rows = []
        import itertools

        zi = 0
        for combo in itertools.product(*[per_axis[k] for k in axes]):
            zi += 1
            vol = 1.0
            for k, (a, b) in zip(axes, combo):
                lo = max(a, box.lo[k]); hi = min(b, box.hi[k])
                vol *= max(hi - lo, 0.0)
            occ_fraction = np.array(occ_rows, float)[:, zi].mean() / len(idx)
            vol_fraction = vol / total_vol
            enr_rows.append((zi, vol_fraction, occ_fraction,
                             occ_fraction / vol_fraction if vol_fraction > 0 else float("nan")))
        tables.append(Table(
            title=f"Zone enrichment {zone}",
            columns=[Field("zone", "1"), Field("volume_fraction", "1"),
                     Field("occupancy_fraction", "1"), Field("enrichment", "1")],
            data=np.array(enr_rows, float), tag="occupancy.enrichment", warnings=warnings,
            meta={"module": "occupancy"},
            selection={"sel": selection.literal, "zone": zone},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))
        return tables

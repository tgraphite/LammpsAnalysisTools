"""lat msd — mean squared displacement and Einstein diffusion coefficient."""

from __future__ import annotations

import numpy as np

from ..core.correlate import linear_fit
from ..core.errors import LatError
from ..core.table import Field, Table
from .base import Module, Param
from . import register


def _parse_fit_range(spec, n_lags):
    if spec in (None, "", "auto"):
        start = max(1, int(0.1 * n_lags))
        stop = max(start + 2, int(0.5 * n_lags))
        return start, min(stop, n_lags)
    parts = str(spec).split(":")
    if len(parts) != 2:
        raise LatError(f"bad -fit_range {spec!r}; expected start:stop")
    return int(parts[0]), int(parts[1])


@register
class MSD(Module):
    name = "msd"
    summary = "mean squared displacement and diffusion coefficient"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel"]
    outputs = ["<base>.msd.out.csv"]
    params = [
        Param("sel", flags=("-sel",), required=True),
        Param("dim", flags=("-dim",), config="modules.msd.dim", type="int", choices=(2, 3)),
        Param("seltype", flags=("-seltype",), config="modules.msd.seltype",
              choices=("atom", "mol_com"), help="atom positions or per-molecule COM"),
        Param("fit_range", flags=("-fit_range",), config="modules.msd.fit_range"),
        Param("remove_com", flags=("-remove_com",), config="modules.msd.remove_com",
              type="flag", help="subtract the selection COM motion"),
        Param("no_remove_com", flags=("-noremove_com",), type="flag"),
        Param("stride", flags=("-stride",), type="int"),
        Param("frame_range", flags=("-frame_range",)),
    ]

    def run(self, ctx, *, sel, dim, fit_range, remove_com, no_remove_com, stride,
            frame_range, seltype):
        selection = ctx.selection(sel)
        if selection.n_atoms == 0:
            raise LatError(f"selection {sel!r} matched 0 atoms")
        idx = selection.indices
        dim = int(dim)
        if dim not in (2, 3):
            raise LatError(f"-dim must be 2 or 3, got {dim}")
        seltype = seltype or "atom"

        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        if len(frames) < 3:
            raise LatError("msd needs at least 3 frames")

        warnings: list[str] = []
        first = next(iter(ctx.frames(frames)))
        unwrapped_ok = first.has("xu", "yu", "zu") or first.has("ix", "iy", "iz")

        if seltype == "mol_com":
            if ctx.topo.mols is None:
                raise LatError("msd -seltype mol_com needs molecule IDs")
            mol_ids = ctx.topo.mols[idx]
            uniq = np.unique(mol_ids)
            groups = [idx[mol_ids == m] for m in uniq]

            traj = np.stack([
                np.array([fr.unwrapped[a].mean(axis=0) for a in groups])
                for fr in ctx.frames(frames)
            ])  # (T, M, 3); unwrapped coords are continuous
        else:
            traj = np.stack([fr.unwrapped[idx] for fr in ctx.frames(frames)])  # (T, N, 3)

        do_com = bool(remove_com) and not no_remove_com
        if do_com:
            traj = traj - traj.mean(axis=1, keepdims=True)
        else:
            warnings.append("COM drift not removed; D may be biased if the system drifts")

        if not unwrapped_ok:
            warnings.append(
                "no image flags (ix/iy/iz) or unwrapped coords (xu/yu/zu); "
                "MSD computed on wrapped coordinates may saturate"
            )

        t, n, _ = traj.shape
        msd = np.zeros((t, 2 + dim))
        msd[:, 0] = np.arange(t)
        for lag in range(t):
            d = traj[lag:] - traj[: t - lag]
            sq = d * d
            msd[lag, 1] = sq[:, :, :dim].sum(axis=2).mean()
            for k in range(dim):
                msd[lag, 2 + k] = sq[:, :, k].mean()

        start, stop = _parse_fit_range(fit_range, t)
        x = msd[start:stop, 0]
        y = msd[start:stop, 1]
        slope, intercept, r2 = linear_fit(x, y)
        factor = 6.0 if dim == 3 else 4.0
        D = slope / factor

        if r2 < 0.8:
            warnings.append(f"low fit quality R^2={r2:.3f}; D is unreliable")

        scalars = {
            "D_A2_per_frame": float(D),
            "R2": float(r2),
            "fit_start": int(start),
            "fit_stop": int(stop),
            "n_origins": int(t),
        }
        dt = ctx.effective_dt(frames)
        if dt:
            # Angstrom^2/frame / (ps/frame) -> cm^2/s
            D_cm2_s = D * (1e-16 / (float(dt) * 1e-12))
            scalars["D_cm2_per_s"] = float(D_cm2_s)
            scalars["frame_interval_ps"] = float(dt)

        comp = "xyz"[:dim]
        cols = [Field("lag", "frame"), Field("msd", "Angstrom^2")] + [
            Field(f"msd_{c}", "Angstrom^2") for c in comp
        ]
        return [Table(
            title="Mean squared displacement",
            columns=cols,
            data=msd,
            tag="msd",
            xlabel="lag (frame)",
            ylabel="MSD (Angstrom^2)",
            definitions={"dim": dim},
            scalars=scalars,
            warnings=warnings,
            meta={"module": "msd"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

"""lat transport — diffusion, Nernst-Einstein and Green-Kubo conductivity."""

from __future__ import annotations

import numpy as np

from ..core.correlate import linear_fit
from ..core.errors import LatError
from ..core.table import Field, Table
from .base import Module, Param
from . import register

F_CONST = 96485.33212  # C/mol
R_CONST = 8.314462618  # J/mol/K
E_CONST = 1.602176634e-19
KB_CONST = 1.380649e-23
NA_CONST = 6.02214076e23


def _msd_slope(traj):
    t = len(traj)
    y = np.array([np.mean(np.sum((traj[lag:] - traj[: t - lag]) ** 2, axis=2))
                  for lag in range(t)])
    start = max(1, int(0.1 * t))
    stop = max(start + 2, int(0.5 * t))
    slope, _, r2 = linear_fit(np.arange(start, stop), y[start:stop])
    return slope, r2


@register
class Transport(Module):
    name = "transport"
    summary = "diffusion, Nernst-Einstein and Green-Kubo transport"
    help_doc = "MANUAL_ALGORITHMS.md"
    inputs = ["-s", "-sel"]
    outputs = ["<base>.transport.d.out.csv", "<base>.transport.ne.out.csv",
               "<base>.transport.gk.out.csv"]
    params = [
        Param("sel", flags=("-sel",), required=True),
        Param("temperature", flags=("-temperature",), config="modules.transport.temperature",
              type="float"),
        Param("velocity_unit", flags=("-velocity_unit",),
              config="modules.transport.velocity_unit"),
        Param("frame_range", flags=("-frame_range",)),
        Param("stride", flags=("-stride",), type="int"),
    ]

    def run(self, ctx, *, sel, temperature, velocity_unit, frame_range, stride):
        topo = ctx.topo
        selection = ctx.selection(sel)
        if selection.n_atoms == 0:
            raise LatError("transport: empty selection")
        idx = selection.indices
        frames = ctx.resolve_frames(frame_range, int(stride) if stride else 1)
        if len(frames) < 3:
            raise LatError("transport needs at least 3 frames")

        first = next(iter(ctx.frames(frames)))
        has_vel = first.has("vx", "vy", "vz")
        warnings: list[str] = []
        if not has_vel:
            raise LatError("transport: trajectory must contain velocity columns vx, vy, vz")

        traj = np.stack([fr.unwrapped[idx] for fr in ctx.frames(frames)])
        traj = traj - traj.mean(axis=1, keepdims=True)

        types = topo.types[idx]
        charges = topo.charges[idx] if topo.charges is not None else None

        dt = ctx.effective_dt(frames)
        box = first.box
        V_m3 = box.volume * 1e-30

        d_rows = []
        species = []
        elem_map = {}
        for t in np.unique(types):
            sub = traj[:, types == t, :]
            slope, r2 = _msd_slope(sub)
            D = slope / 6.0
            el = topo.type_element.get(int(t), str(int(t)))
            elem_map[int(t)] = el
            # D in Angstrom^2/frame -> cm^2/s (1e-4) and m^2/s (1e-8) per ps
            D_cm2_s = D / float(dt) * 1e-4 if dt else float("nan")
            d_rows.append((int(t), int((types == t).sum()), D, D_cm2_s, r2))
            species.append((int(t), el, int((types == t).sum()), D))
        if dt is None:
            warnings.append("frame_interval not set; D reported only in Angstrom^2/frame")

        definitions = {"type_element": "; ".join(f"{k}={v}" for k, v in elem_map.items())}
        if dt is not None:
            definitions["frame_interval"] = f"{dt} ps"
        tables = [Table(
            title="Diffusion coefficients",
            columns=[Field("type", "1"), Field("n", "1"),
                     Field("D_A2_per_frame", "Angstrom^2/frame"),
                     Field("D_cm2_s", "cm^2/s"), Field("r2", "1")],
            data=np.array(d_rows, float),
            tag="transport.d",
            xlabel="type",
            ylabel="D / r^2",
            definitions=definitions, warnings=warnings, meta={"module": "transport"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )]

        if charges is None:
            warnings.append("no charges available; NE/GK conductivity skipped")
            return tables
        if not temperature:
            warnings.append("no -temperature; NE/GK conductivity skipped")
            return tables
        if dt is None:
            warnings.append("no frame_interval; NE/GK conductivity skipped")
            return tables

        ne_rows = []
        denom = 0.0
        contrib = []
        for (t, el, n, D) in species:
            q = float(np.mean(charges[types == t]))
            c = n / (NA_CONST * V_m3)
            D_si = D / float(dt) * 1e-8  # Angstrom^2/frame -> m^2/s
            w = q * q * c * D_si
            contrib.append(w)
            denom += w
            ne_rows.append((el, q, c, D_si, w))
        sigma_ne = F_CONST ** 2 / (R_CONST * float(temperature)) * denom
        ne_definitions = dict(definitions)
        ne_definitions["temperature"] = f"{temperature} K"
        tables.append(Table(
            title="Nernst-Einstein conductivity",
            columns=[Field("type", "1"), Field("charge", "e"),
                     Field("conc", "mol/m^3"), Field("D", "m^2/s"), Field("t_i", "1")],
            data=np.array([[i, q, c, D, (w / denom if denom else 0.0)]
                           for i, (_, q, c, D, w) in enumerate(ne_rows)], float),
            tag="transport.ne",
            xlabel="type",
            ylabel="D / transport number",
            definitions=ne_definitions,
            scalars={"sigma_NE_S_m": float(sigma_ne), "temperature_K": float(temperature)},
            warnings=warnings + ["NE ignores ion-ion correlation; typically overestimates sigma"],
            meta={"module": "transport"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        ))

        # Green-Kubo
        v_unit = velocity_unit or "A/ps"
        if v_unit == "A/ps":
            scale = 1e-10 / 1e-12
        elif v_unit == "A/fs":
            scale = 1e-10 / 1e-15
        else:
            raise LatError(f"transport: velocity_unit must be A/ps or A/fs, got {v_unit!r}")
        vel = np.stack([fr.velocities[idx] for fr in ctx.frames(frames)]) * scale
        q = charges
        J = np.einsum("tac,a->tc", vel, q) * E_CONST  # C m/s
        T = len(J)
        acf = np.zeros(T)
        for lag in range(T):
            acf[lag] = np.sum(J[: T - lag] * J[lag:], axis=1).mean()
        dt_s = float(dt) * 1e-12
        running = np.cumsum(acf) * dt_s
        # trapezoidal integration removes the missing half weight of acf[0]
        running = running - 0.5 * acf[0] * dt_s
        sigma_gk = running / (3.0 * V_m3 * KB_CONST * float(temperature))
        gk_definitions = dict(ne_definitions)
        gk_definitions["velocity_unit"] = v_unit
        gk = Table(
            title="Green-Kubo running conductivity",
            columns=[Field("lag", "frame"), Field("sigma_running", "S/m")],
            data=np.column_stack([np.arange(T), sigma_gk]),
            tag="transport.gk",
            xlabel="lag (frame)",
            ylabel="sigma (S/m)",
            definitions=gk_definitions,
            scalars={"sigma_GK_S_m": float(sigma_gk[-1]),
                     "sigma_NE_S_m": float(sigma_ne),
                     "temperature_K": float(temperature)},
            warnings=warnings, meta={"module": "transport"},
            selection={"sel": selection.literal},
            frame_range=ctx.frame_range_str(frame_range, frames),
        )
        tables.append(gk)
        return tables

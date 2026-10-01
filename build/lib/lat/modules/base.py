"""Module contract: a pure function module(ctx, **params) -> list[Table]."""

from __future__ import annotations

from dataclasses import dataclass, field

from ..core.errors import LatError


@dataclass
class Param:
    name: str
    flags: tuple = ()
    config: str | None = None
    type: str = "str"  # str|int|float|flag
    required: bool = False
    nargs: str | int | None = None
    help: str = ""
    choices: tuple = ()


class Module:
    name: str = ""
    summary: str = ""
    help_doc: str = ""  # relative path under docs/, used by `lat describe`
    params: list[Param] = []
    inputs: list[str] = []
    outputs: list[str] = []
    requires_topology: bool = True
    screen_only: bool = False

    def run(self, ctx, **params):
        raise NotImplementedError

    # ---- helpers ------------------------------------------------------
    @classmethod
    def resolve_defaults(cls, cfg, provided: dict) -> dict:
        """Fill every declared param from CLI-provided values, else config.

        No code-level physical defaults: a missing value is an ERROR.
        """
        out = {}
        for p in cls.params:
            if provided.get(p.name) is not None:
                out[p.name] = provided[p.name]
                continue
            if p.config:
                value = cfg.get(p.config)
                if value is not None:
                    out[p.name] = value
                    continue
            if p.required:
                raise LatError(
                    f"{cls.name}: required parameter {p.name!r} "
                    f"({', '.join(p.flags)}) was not provided"
                )
            out[p.name] = None
        return out

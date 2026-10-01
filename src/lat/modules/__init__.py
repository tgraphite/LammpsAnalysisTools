"""Module registry."""

from __future__ import annotations

from ..core.errors import LatError
from .base import Module

_REGISTRY: dict[str, type[Module]] = {}


def register(cls: type[Module]) -> type[Module]:
    if not cls.name:
        raise LatError(f"module {cls.__name__} has no name")
    _REGISTRY[cls.name] = cls
    return cls


def get(name: str) -> type[Module]:
    if name not in _REGISTRY:
        raise LatError(f"unknown module {name!r}; known: {', '.join(sorted(_REGISTRY))}")
    return _REGISTRY[name]


def names() -> list[str]:
    return sorted(_REGISTRY)


def all_modules() -> dict[str, type[Module]]:
    return dict(_REGISTRY)


def _load_builtin() -> None:
    from . import (  # noqa: F401
        agg,
        cluster,
        contact,
        density,
        describe,
        gyrate,
        hbond,
        msd,
        occupancy,
        pistack,
        probe,
        rdf,
        residence,
        rms,
        select,
        transport,
    )


_load_builtin()

"""Load the single program-root config/lat.yaml.

Rules (M11 + decision 18):
  * exactly one config file, at <program-root>/config/lat.yaml;
  * the codebase must not carry fallback defaults for physical parameters;
    a missing key is an ERROR;
  * CLI arguments override the config.
"""

from __future__ import annotations

import os
from pathlib import Path

import yaml

from . import __version__
from .core.errors import LatError, get_logger

_MISSING = object()

# src/lat/config.py -> src/lat -> src -> <program root>  (source checkout)
_PROGRAM_ROOT = Path(__file__).resolve().parents[2]
# config bundled inside the installed package (for `pip install` on a new device)
_PACKAGE_DATA = Path(__file__).resolve().parent / "data" / "lat.yaml"
_SOURCE_CONFIG = _PROGRAM_ROOT / "config" / "lat.yaml"


def program_root() -> Path:
    return _PROGRAM_ROOT


def config_path() -> Path:
    """Resolution order: $LAT_CONFIG > source checkout > bundled package data.

    Set LAT_CONFIG=/path/to/lat.yaml to use a project-specific configuration.
    """
    env = os.environ.get("LAT_CONFIG")
    if env:
        return Path(env)
    if _SOURCE_CONFIG.is_file():
        return _SOURCE_CONFIG
    return _PACKAGE_DATA


class Config:
    def __init__(self, raw: dict, path: Path):
        self.raw = raw
        self.path = path

    def get(self, dotted: str, default=_MISSING):
        node = self.raw
        for key in dotted.split("."):
            if not isinstance(node, dict) or key not in node:
                if default is _MISSING:
                    return _MISSING
                return default
            node = node[key]
        return node

    def require(self, dotted: str):
        value = self.get(dotted)
        if value is _MISSING:
            raise LatError(
                f"config key {dotted!r} missing in {self.path}; no code-level fallback is allowed"
            )
        return value


def load(path=None) -> Config:
    p = Path(path) if path else config_path()
    if not p.is_file():
        raise LatError(f"config file not found: {p}")
    with open(p) as fh:
        raw = yaml.safe_load(fh)
    if not isinstance(raw, dict):
        raise LatError(f"config file {p} must contain a mapping at top level")
    declared = (raw.get("meta") or {}).get("version")
    if declared is not None and str(declared) != __version__:
        get_logger().warning(
            f"config {p} declares meta.version {declared} but this lat is "
            f"{__version__}; module parameters may not match this release"
        )
    return Config(raw, p)

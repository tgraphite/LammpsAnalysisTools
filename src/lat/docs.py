"""Helpers shared by the per-module help manuals and their static checks.

The parameter table in every manual is generated from the module class and the
active configuration.  ``tests/unit/test_manuals.py`` re-renders the same table
and fails when a manual drifts from the code.
"""

from __future__ import annotations

from pathlib import Path

from .config import Config, manual_path

_TYPE_LABEL = {
    "str": "字符串",
    "int": "整数",
    "float": "浮点数",
    "flag": "开关",
}


def _fmt_default(param, cfg: Config | None) -> str:
    if param.config and cfg is not None:
        value = cfg.get(param.config)
        if value is not None:
            if isinstance(value, bool):
                return "true" if value else "false"
            return str(value)
    if param.required:
        return "必填"
    return "—"


def param_rows(cls, cfg: Config | None = None) -> list[list[str]]:
    rows = []
    for p in cls.params:
        flags = ", ".join(f"`{f}`" for f in (p.flags or (p.name,)))
        rows.append([
            flags,
            _TYPE_LABEL.get(p.type, p.type),
            _fmt_default(p, cfg),
            p.help or "—",
        ])
    return rows


def render_param_table(cls, cfg: Config | None = None) -> str:
    lines = ["| 参数 | 类型 | 默认 | 说明 |", "|---|---|---|---|"]
    for row in param_rows(cls, cfg):
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def manual_for(name: str) -> Path | None:
    return manual_path(f"manual/{name}.md")

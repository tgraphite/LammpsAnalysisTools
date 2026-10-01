"""The only place that writes analysis files.

Implements the frozen CSV header contract, atomic writes, and conflict
backup (.bak.1/.bak.2...).
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from .errors import LatError, get_logger
from .table import Table


def _fmt(value) -> str:
    if value is None:
        return "none"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def _scalar_str(scalars: dict) -> str:
    if not scalars:
        return "none"
    return "; ".join(f"{k}={_fmt(v)}" for k, v in scalars.items())


def _kv_str(d: dict) -> str:
    if not d:
        return "none"
    return "; ".join(f"{k}={v}" for k, v in d.items())


def _header(table: Table, *, include_audit: bool, command: str, version: str,
            prefix: str = "#") -> list[str]:
    module = table.meta.get("module", table.tag.split(".")[0])
    cols = ", ".join(c.name for c in table.columns)
    warnings = "; ".join(table.warnings) if table.warnings else "none"
    lines = [
        f"{prefix} program: LAT version {version}",
        f"{prefix} module: {module}",
    ]
    if include_audit:
        lines.append(f"{prefix} command: {command}")
    lines += [
        prefix,
        f"{prefix} selections: {_kv_str(table.selection)}",
        f"{prefix} definitions: {_kv_str(table.definitions)}",
        f"{prefix} frame_range: {table.frame_range}",
        f"{prefix} scalars: {_scalar_str(table.scalars)}",
        f"{prefix} warning: {warnings}",
        prefix,
        f"{prefix} title: {table.title}",
        f"{prefix} xlabel: {table.x_label}",
        f"{prefix} ylabel: {table.y_label}",
        prefix,
        f"{prefix} columns: {cols}",
        f"{prefix} --- data ---",
    ]
    return lines


def backup_existing(path: Path) -> Path | None:
    """Rotate an existing target to .bak.1, shifting older backups up."""
    if not path.exists():
        return None
    n = 1
    while path.with_name(path.name + f".bak.{n}").exists():
        n += 1
    for i in range(n - 1, 0, -1):
        src = path.with_name(path.name + f".bak.{i}")
        dst = path.with_name(path.name + f".bak.{i + 1}")
        os.replace(src, dst)
    bak = path.with_name(path.name + ".bak.1")
    os.replace(path, bak)
    return bak


def write_table(
    table: Table,
    *,
    basename: str,
    outdir: Path | str = ".",
    ext: str = "csv",
    suffix: str = "out",
    prefix: str = "#",
    float_format: str = "%.6g",
    atomic: bool = True,
    backup_on_conflict: bool = True,
    include_audit: bool = True,
    command: str = "",
    version: str = "",
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    target = outdir / table.filename(basename, ext=ext, suffix=suffix)

    if backup_on_conflict:
        bak = backup_existing(target)
        if bak is not None:
            get_logger().info(f"existing {target.name} moved to {bak.name}")

    fmt = "%s" if float_format in ("%s", "str") else float_format
    body_lines = []
    for row in table.data:
        body_lines.append(",".join(fmt % v for v in row))
    header = _header(table, include_audit=include_audit, command=command,
                     version=version, prefix=prefix)
    text = "\n".join(header + body_lines) + "\n"

    if atomic:
        fd, tmp = tempfile.mkstemp(dir=str(outdir), prefix=".lat_tmp_", suffix="." + ext)
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write(text)
            os.replace(tmp, target)
        except BaseException:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise
    else:
        with open(target, "w") as fh:
            fh.write(text)
    return target


def write_all(tables, **kwargs) -> list[Path]:
    if not tables:
        raise LatError("module produced no tables to write")
    return [write_table(t, **kwargs) for t in tables]

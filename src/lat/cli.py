"""lat command line: lat <module> -s <data> -f <traj> [...]."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .config import load as load_config
from .core.context import Context
from .core.errors import LatError, Logger, get_logger, set_logger
from .core.io import open_system
from .core.output import write_all
from .core.selection import read_ndx
from .modules import all_modules, get

_ALLOWED_UNITS = {"real", "metal", "lj"}


def _type_of(type_name):
    return {"str": str, "int": int, "float": float}.get(type_name)


def build_parser(cls, cfg=None) -> argparse.ArgumentParser:
    epilog = (f"manual: docs/MANUAL_ALGORITHMS.md  |  `lat describe -module {cls.name}` "
              "for parameters; defaults come from the lat config")
    parser = argparse.ArgumentParser(
        prog=f"lat {cls.name}", description=cls.summary, epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter, add_help=True
    )
    if getattr(cls, "requires_topology", True):
        parser.add_argument("-s", "-topo", dest="topo", required=True, help="LAMMPS data file")
        parser.add_argument("-f", "-traj", dest="traj", help="LAMMPS dump / GROMACS traj")
        parser.add_argument("-o", dest="output", help="output basename")
        parser.add_argument("-unit", dest="unit", help="length unit (default: real/Angstrom)")
        parser.add_argument("-style", dest="style", help="force atom_style")
        parser.add_argument("-types", dest="types", help="explicit type->element file")
        parser.add_argument("-index", "-n", dest="index", help="GROMACS .ndx with named groups")
        parser.add_argument("-pbc", dest="pbc", action="store_true", default=None)
        parser.add_argument("-nopbc", dest="nopbc", action="store_true")
        parser.add_argument("-debug", dest="debug", action="store_true")
        parser.add_argument("-v", dest="verbose", action="store_true")
    for p in cls.params:
        help_text = p.help
        if p.config and cfg is not None:
            val = cfg.get(p.config)
            if val is not None:
                help_text = f"{help_text} [default: {val}]".strip()
        kwargs = {"help": help_text}
        if p.choices:
            kwargs["choices"] = list(p.choices)
        if p.type == "flag":
            parser.add_argument(*p.flags, dest=p.name, action="store_true", default=None, **kwargs)
        else:
            parser.add_argument(*p.flags, dest=p.name, type=_type_of(p.type), default=None, **kwargs)
    return parser


def _parse_types_file(path) -> dict[int, str]:
    mapping = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 2:
                continue
            a, b = parts[0], parts[1]
            if a.isdigit():
                mapping[int(a)] = b.capitalize()
            elif b.isdigit():
                mapping[int(b)] = a.capitalize()
    if not mapping:
        raise LatError(f"{path}: no type->element pairs parsed")
    return mapping


def _run_plot(kind: str, argv) -> int:
    parser = argparse.ArgumentParser(prog=f"lat {kind}")
    parser.add_argument("-f", dest="csv", required=True,
                        help="lat output CSV (header carries title/xlabel/ylabel/columns)")
    parser.add_argument("-o", dest="out", help="output image path")
    parser.add_argument("-title", dest="title", help="override the header title")
    parser.add_argument("-smooth", dest="smooth", choices=("none", "gaussian", "savgol"),
                        help="smooth the plotted curve(s); the CSV stays raw")
    if kind == "plot1d":
        parser.add_argument("-columns", dest="columns", nargs="+", metavar="NAME",
                            help="only draw these y columns (comma- or space-separated); "
                                 "the x axis is always the first column")
        parser.add_argument("-bar", dest="bar", action="store_true",
                            help="grouped bar chart, one bar per y column")
        parser.add_argument("-stack", dest="stack", action="store_true",
                            help="stacked bar chart; y columns are normalized to "
                                 "occupancy fractions")
    args = parser.parse_args(argv)
    from .plotting import plot1d, plot2d

    cfg = load_config()
    smooth = args.smooth or str(cfg.get("plot.smooth", default="none") or "none")
    sm = {
        "window": int(cfg.get("plot.savgol_window", default=5) or 5),
        "order": int(cfg.get("plot.savgol_order", default=2) or 2),
        "sigma": float(cfg.get("plot.gaussian_sigma", default=1.0) or 1.0),
    }
    if kind == "plot1d":
        if args.bar and args.stack:
            raise LatError("-bar and -stack are mutually exclusive")
        col_sel = None
        if args.columns:
            col_sel = [c for tok in args.columns
                       for c in (t.strip() for t in tok.split(",")) if c]
        mode = "bar" if args.bar else "stack" if args.stack else "line"
        out = plot1d(args.csv, out=args.out, title=args.title, smooth=smooth,
                     columns=col_sel, mode=mode, **sm)
    else:
        out = plot2d(args.csv, out=args.out, title=args.title, smooth=smooth, **sm)
    print(f"wrote {out}", file=sys.stderr)
    return 0


def _apply_config_paths(cfg) -> None:
    pass


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in ("-h", "--help", "help"):
        _top_help()
        return 0
    if argv[0] in ("-V", "--version"):
        print(f"lat {__version__}")
        return 0

    module_name = argv[0]
    if module_name in ("plot1d", "plot2d"):
        try:
            return _run_plot(module_name, argv[1:])
        except LatError as exc:
            print(f"[lat ERROR] {exc}", file=sys.stderr)
            return 1

    try:
        cfg = load_config()
        logger = Logger(level="warning")
        set_logger(logger)
        cls = get(module_name)
        parser = build_parser(cls, cfg)
        args = parser.parse_args(argv[1:])

        if getattr(args, "verbose", None):
            logger.level = "info"
        elif cfg.get("logger.level"):
            logger.level = str(cfg.get("logger.level"))

        command = "lat " + " ".join(argv)
        _execute(cls, args, cfg, command)
        return 0
    except LatError as exc:
        print(f"[lat ERROR] {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("[lat ERROR] interrupted", file=sys.stderr)
        return 130


def _execute(cls, args, cfg, command):
    logger = get_logger()
    provided = {}
    for p in cls.params:
        provided[p.name] = getattr(args, p.name, None)
    params = cls.resolve_defaults(cfg, provided)

    mod = cls()
    if not getattr(cls, "requires_topology", True):
        mod.run(None, **params)
        return []

    unit = args.unit or cfg.get("io.default_unit")
    if args.unit is None:
        logger.info(f"no -unit given; assuming {unit} (Angstrom)")
    if unit not in _ALLOWED_UNITS:
        raise LatError(f"unit {unit!r} not supported (allowed: {sorted(_ALLOWED_UNITS)})")

    ndx = read_ndx(args.index) if args.index else None
    import time as _time
    _t0 = _time.time()
    logger.notice(f"reading topology {args.topo}")
    topo, source = open_system(args.topo, args.traj, style=args.style, cfg=cfg)
    logger.notice(f"  parsed {topo.n_atoms} atoms in {_time.time() - _t0:.1f}s"
                  + (f"; trajectory {args.traj}" if args.traj else ""))

    if args.types:
        topo.type_element.update(_parse_types_file(args.types))
        topo.element_problems = {}

    ctx = Context(topo, source, cfg=cfg, unit=unit, command=command,
                  debug=bool(args.debug), ndx=ndx)

    tables = mod.run(ctx, **params)
    emitted = set()
    for t in tables:
        for w in t.warnings:
            if w not in emitted:
                logger.warning(w)
                emitted.add(w)
    if getattr(cls, "screen_only", False) or not tables:
        return tables

    basename = args.output or Path(args.topo).stem
    outdir = Path(basename).parent if Path(basename).parent != Path("") else Path(".")
    stem = Path(basename).name
    out_cfg = cfg.get("output", default={}) or {}
    paths = write_all(
        tables,
        basename=stem,
        outdir=outdir,
        ext=out_cfg.get("ext", "csv"),
        suffix=out_cfg.get("suffix", "out"),
        prefix=out_cfg.get("header_prefix", "#"),
        float_format=out_cfg.get("float_format", "%.6g"),
        atomic=bool(out_cfg.get("atomic_write", True)),
        backup_on_conflict=bool(cfg.get("io.backup_on_conflict", default=True)),
        include_audit=bool(out_cfg.get("include_audit", True)),
        command=command,
        version=__version__,
    )
    for path in paths:
        logger.notice(f"wrote {path}")
    if args.debug and out_cfg.get("debug_json", False):
        audit = {
            "command": command,
            "version": __version__,
            "config": str(cfg.path),
            "tables": [{"tag": t.tag, "rows": t.n_rows, "scalars": t.scalars,
                        "selection": t.selection, "warnings": t.warnings} for t in tables],
            "warnings": logger.warnings,
        }
        apath = Path(outdir) / f"{stem}.debug.json"
        with open(apath, "w") as fh:
            json.dump(audit, fh, indent=2)
        logger.info(f"wrote {apath}")
    return tables


def _top_help() -> None:
    print("lat — LAMMPS Analysis Tool")
    print(f"version {__version__}")
    print("\nusage: lat <module> -s <data> [-f <traj>] [options]\n")
    print("modules:")
    for name in sorted(all_modules()):
        print(f"  {name:<12} {get(name).summary}")
    print("\nplotting:")
    print("  plot1d       plot a 1D CSV")
    print("  plot2d       plot a 2D CSV")
    print("\n  describe -module <name>   show a module's parameters and contract")


if __name__ == "__main__":
    sys.exit(main())

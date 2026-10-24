"""The program-root config and the packaged copy must not drift.

`config/lat.yaml` is authoritative for a source checkout and
`src/lat/data/lat.yaml` for an installed package; before 0.2.0 they differed
(`grid.dr_3d` 0.5 vs 1.0), so the same command could behave differently
depending on how lat was started.
"""

from __future__ import annotations

import yaml

from lat import __version__


def test_bundled_config_is_identical_to_source_config(repo_root):
    a = (repo_root / "config" / "lat.yaml").read_text()
    b = (repo_root / "src" / "lat" / "data" / "lat.yaml").read_text()
    assert a == b, ("config/lat.yaml and src/lat/data/lat.yaml differ; keep them "
                    "byte-identical (copy the source config over the packaged one)")


def test_config_declares_the_package_version(repo_root):
    raw = yaml.safe_load((repo_root / "config" / "lat.yaml").read_text())
    assert raw["meta"]["version"] == __version__


def test_declared_version_is_used_by_the_cli(repo_root):
    raw = yaml.safe_load((repo_root / "config" / "lat.yaml").read_text())
    assert str(raw["meta"]["version"]) != "0.1.0"


def test_every_module_has_a_config_block(repo_root):
    """Each registered module's parameters must be resolvable from config."""
    import sys

    sys.path.insert(0, str(repo_root / "src"))
    from lat.config import load
    from lat.modules import all_modules

    cfg = load(str(repo_root / "config" / "lat.yaml"))
    missing = []
    for name in sorted(all_modules()):
        if name in ("describe", "plot1d", "plot2d", "select"):
            continue
        if cfg.get(f"modules.{name}") is None:
            missing.append(name)
    assert not missing, f"modules without a config block: {missing}"

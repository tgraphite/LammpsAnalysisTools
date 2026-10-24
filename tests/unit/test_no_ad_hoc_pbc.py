"""Architecture guard: periodic geometry lives in exactly one place.

Before 0.2.0 there were four parallel "periodic geometry" implementations in
``src/lat/modules`` and three modules silently treated wrapped coordinates as
physical (SASA lost every neighbour across a face; contact/residence ``-mode
center`` averaged wrapped coordinates; agg/cluster/pistack folded chains longer
than half a box). See ``docs/07_decision_log.md`` §6.2.

This test fails when a module goes back to hand-rolling PBC instead of calling
:mod:`lat.core.periodic`.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

MODULES_DIR = Path(__file__).resolve().parents[2] / "src" / "lat" / "modules"

# Tokens that must never appear in a module: each one is a way to compute a
# periodic quantity that this project has already got wrong once.
BLOCKED_TOKENS = {
    "cKDTree": "use lat.core.periodic.periodic_neighbors / pairs_within",
    "KDTree(": "use lat.core.periodic.periodic_neighbors",
    "distance_transform_edt": "use lat.core.periodic.periodic_distance_transform",
    "ndimage.label": "use lat.core.periodic.label_periodic",
    "generate_binary_structure": "use lat.core.periodic.label_periodic",
    "minimum_image": "use lat.core.periodic.whole_positions / periodic_distance",
    "np.pad": "wrap-aware padding belongs in lat.core.periodic.label_periodic",
}

# Modules that must obtain a specific periodic primitive, by name.
REQUIRED = {
    "_graph.py": ("residue_coms",),
    "contact.py": ("residue_com", "periodic_distance"),
    "pistack.py": ("whole_positions",),
    "residence.py": ("residue_coms",),
}

# A bare positional mean is the "wrapped centroid" bug. These modules average
# something that is not a wrapped position: unwrapped coordinates (msd), box
# lengths (density, hbond), coordinates already made whole by periodic
# (rms, pistack), or cluster-superposed frames (rms).
POSITIONAL_MEAN_ALLOWED = {"msd.py", "rms.py", "density.py", "hbond.py", "pistack.py"}
# rdf is deliberately absent from REQUIRED: `-intramol` counts *pairwise*
# minimum-image distances, which are already exact on wrapped coordinates.
# tests/regression/level2_format/test_pbc_invariance.py pins that.

_MEAN_AXIS0 = re.compile(r"\.mean\(\s*axis\s*=\s*0")


def _module_files():
    return sorted(p for p in MODULES_DIR.glob("*.py"))


def test_modules_directory_exists():
    assert MODULES_DIR.is_dir(), f"module directory not found: {MODULES_DIR}"


@pytest.mark.parametrize("path", _module_files(), ids=lambda p: p.name)
def test_no_hand_rolled_periodic_geometry(path):
    text = path.read_text()
    hits = [f"{tok} ({why})" for tok, why in BLOCKED_TOKENS.items() if tok in text]
    assert not hits, f"{path.name} hand-rolls periodic geometry: " + "; ".join(hits)


@pytest.mark.parametrize("name,symbols", sorted(REQUIRED.items()))
def test_pbc_sensitive_module_routes_through_periodic(name, symbols):
    path = MODULES_DIR / name
    assert path.exists(), f"{name} is missing"
    text = path.read_text()
    assert "core.periodic import" in text or "core import periodic" in text, (
        f"{name} must import from lat.core.periodic"
    )
    for sym in symbols:
        assert sym in text, f"{name} must use lat.core.periodic.{sym}"


@pytest.mark.parametrize("path", _module_files(), ids=lambda p: p.name)
def test_no_wrapped_positional_mean_outside_allowlist(path):
    if path.name in POSITIONAL_MEAN_ALLOWED:
        pytest.skip("averages unwrapped coordinates or box lengths, not positions")
    hits = _MEAN_AXIS0.findall(path.read_text())
    assert not hits, (
        f"{path.name} averages positions with a bare .mean(axis=0); a wrapped "
        f"centroid is wrong for anything crossing a face -- use "
        f"lat.core.periodic.residue_com"
    )


def test_probe_module_is_gone():
    """0.2.0 split `probe` into `freevol` + `sasa`; no compatibility shim."""
    assert not (MODULES_DIR / "probe.py").exists(), (
        "probe.py must stay deleted: its grid modes are lat freevol and its "
        "surface mode is lat sasa"
    )


def test_freevol_and_sasa_modules_exist():
    assert (MODULES_DIR / "freevol.py").exists()
    assert (MODULES_DIR / "sasa.py").exists()

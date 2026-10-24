#!/usr/bin/env python3
"""Re-run every LAT artifact under examples/ that this release changes.

Why this exists
---------------
lat 0.2.0 fixes three defects that changed numbers in already-published example
outputs, and splits one module (`probe` -> `freevol` + `sasa`). The owners of
those outputs must be able to see, per file, *what* changed and *why* before the
new numbers overwrite the old ones -- so nothing is overwritten here: every run
writes into this directory, mirroring the original relative path, and
MAPPING.tsv records the old path, the new path and the command.

Usage
-----
    python3 recheck.py --list                 # show the job table
    python3 recheck.py --affected             # re-run only the changed modules
    python3 recheck.py --verify               # re-run unchanged modules, expect identical
    python3 recheck.py --report               # compare mirror vs originals

`--report` then writes MAPPING.tsv and prints the per-file differences.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent                      # _lat_v020_recheck
EXAMPLES = HERE.parent                                      # examples/
LAT = os.environ.get("LAT", "/opt/lat/lat")
PY = sys.executable

# ---------------------------------------------------------------------------
# job table: (relative workdir, reason, status, [lat argv...])
#   status "affected" -> numbers are expected to change (defect fix / module split)
#   status "verify"   -> numbers must be bit-identical (code path untouched)
# `-o` is retargeted to this directory automatically.
# ---------------------------------------------------------------------------
D31 = "31_linear_polymer/analysis"
D22 = "22_surfactant/8x_system"
D24 = "24_electrolyte_electrode/sandwich_licl"
D24C = "24_electrolyte_electrode/sandwich_licl_cpm"
D21 = "21_ethanol_water"
D23 = "23_electrolyte"
D22A = "22_surfactant/1x_system"
D32 = "32_molecular_crystal"

JOBS = [
    # ---- 31: bond-graph unwrapping (Bonds reference atom ids, not file rows) --
    (D31, "affected", "bond-id order fix (gyrate read Rg of the wrong atoms)",
     ["gyrate", "-s", "system.in.data", "-f", "run.langevin.out.lammpstrj",
      "-sel", "all", "-o", "run"]),
    (D31, "affected", "bond-id order fix (rmsd/rmsf read the wrong atoms)",
     ["rms", "-s", "system.in.data", "-f", "run.langevin.out.lammpstrj",
      "-sel", "all", "-rmsf", "-o", "run"]),
    # ---- 31: probe split + periodic SASA + periodic grid --------------------
    (D31, "affected", "probe -> freevol; PBC-correct occupancy unchanged, new tag",
     ["freevol", "-s", "system.in.data", "-f", "run.langevin.out.lammpstrj",
      "-mode", "ffv", "-o", "run"]),
    (D31, "affected", "probe -> freevol; PSD redefined (vdW surface + Voronoi norm)",
     ["freevol", "-s", "system.in.data", "-f", "run.langevin.out.lammpstrj",
      "-mode", "psd", "-stride", "10", "-o", "run"]),
    (D31, "affected", "probe -> freevol; torus-closed connectivity, 26-connectivity",
     ["freevol", "-s", "system.in.data", "-f", "run.langevin.out.lammpstrj",
      "-mode", "connectivity", "-o", "run"]),
    (D31, "affected", "probe -> sasa; SASA is now minimum-image (was +14%)",
     ["sasa", "-s", "system.in.data", "-f", "run.langevin.out.lammpstrj",
      "-o", "run"]),
    # ---- 22: hbond donor definition (C-H was counted as a donor) -------------
    (D22, "affected", "hbond donor fix (only N/O/F/S/Cl/Br/I bonded to H donate)",
     ["hbond", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-sel1", "all", "-sel2", "all", "-o", "hbond_any", "-frame_range", "160:200"]),
    (D22, "affected", "hbond donor fix (only N/O/F/S/Cl/Br/I bonded to H donate)",
     ["hbond", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-sel1", "type 2 3 4 5 6", "-sel2", "type 9 10",
      "-o", "hbond_surfactant_water", "-frame_range", "160:200"]),
    # ---- 24: hbond donor + contact COM --------------------------------------
    (D24, "affected", "hbond donor fix",
     ["hbond", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-sel1", "all", "-sel2", "all", "-o", "any"]),
    (D24, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 10", "-sel2", "type 1-8", "-o", "li_go"]),
    (D24, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 10", "-sel2", "type 9", "-o", "li_gr"]),
    (D24, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 11", "-sel2", "type 1-8", "-o", "cl_go"]),
    (D24, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 11", "-sel2", "type 9", "-o", "cl_gr"]),
    (D24C, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.cpm.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 10", "-sel2", "type 1-8", "-o", "li_go"]),
    (D24C, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.cpm.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 10", "-sel2", "type 9", "-o", "li_gr"]),
    (D24C, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.cpm.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 11", "-sel2", "type 1-8", "-o", "cl_go"]),
    (D24C, "affected", "contact -mode min: com_dist now uses whole-molecule centres",
     ["contact", "-s", "prerun.out.data", "-f", "run.cpm.out.lammpstrj",
      "-cutoff", "5", "-sel1", "type 11", "-sel2", "type 9", "-o", "cl_gr"]),
    # ---- remaining artifacts: unchanged code paths, re-run for bit-identity --
    # ---- unchanged code paths, re-run to prove bit-identity ------------------
    (D21, "verify", "rdf untouched", ["rdf", "-s", "system.data", "-f", "equil.out.lammpstrj",
      "-sel1", "type 6", "-sel2", "type 6", "-o", "owow", "-frame_range", "100:200",
      "-cutoff", "15"]),
    (D21, "verify", "rdf untouched", ["rdf", "-s", "system.data", "-f", "equil.out.lammpstrj",
      "-sel1", "type 5", "-sel2", "type 6", "-o", "oeow", "-frame_range", "100:200",
      "-cutoff", "15"]),
    (D22A, "verify", "rdf untouched", ["rdf", "-s", "system.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-cutoff", "15", "-sel1", "type 7 8",
      "-sel2", "type 9 10", "-o", "oil_water"]),
    (D23, "verify", "agg -mode min untouched (config default is min)",
     ["agg", "-s", "run.out.data", "-f", "run.out.lammpstrj", "-sel1", "type 1",
      "-sel2", "type 2-5", "-solv", "type 6-18", "-o", "li_fsi", "-cutoff", "2.9",
      "-frame_range", "0:200:10"]),
    (D23, "verify", "cluster -mode min untouched",
     ["cluster", "-s", "run.out.data", "-f", "run.out.lammpstrj", "-sel1", "type 1",
      "-sel2", "type 6-18", "-cutoff", "2.9", "-frame_range", "100:200", "-o", "li_solv"]),
    (D23, "verify", "rdf untouched", ["rdf", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-sel1", "type 1", "-sel2", "type 6-18", "-frame_range", "100:200", "-o", "li_solv"]),
    (D24, "verify", "density untouched", ["density", "-s", "prerun.out.data",
      "-f", "run.efield.out.lammpstrj", "-sel", "type 10", "-frame_range", "500:600",
      "-o", "li_end", "-nbins", "200"]),
    (D24, "verify", "msd untouched", ["msd", "-s", "prerun.out.data",
      "-f", "run.efield.out.lammpstrj", "-o", "li", "-sel", "type 10"]),
    ("32_molecular_crystal/run_constraincof", "verify", "density untouched",
     ["density", "-s", "prerun.freezecof.out.data", "-f", "run.out.lammpstrj",
      "-sel", "type 10", "-o", "li", "-dim", "2", "-axis", "z"]),
]
JOBS += [
    # 21: the third rdf command
    (D21, "verify", "rdf untouched", ["rdf", "-s", "system.data", "-f", "equil.out.lammpstrj",
      "-sel1", "type 1 2 3 4 5", "-sel2", "type 6 7", "-o", "ew", "-frame_range", "100:200",
      "-cutoff", "15"]),
    # 22/1x: the two remaining rdf commands
    (D22A, "verify", "rdf untouched", ["rdf", "-s", "system.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-cutoff", "15", "-sel1", "type 5 6", "-sel2", "type 6",
      "-o", "surf_so_water_o"]),
    (D22A, "verify", "rdf untouched", ["rdf", "-s", "system.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-cutoff", "15", "-sel1", "type 2 3 4", "-sel2", "type 7 8",
      "-o", "surf_org_oil"]),
    # 22/8x: rdf + density (not touched by 0.2.0)
    (D22, "verify", "rdf untouched", ["rdf", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-cutoff", "15", "-sel1", "type 7 8", "-sel2", "type 9 10",
      "-o", "oil_water"]),
    (D22, "verify", "rdf untouched", ["rdf", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-cutoff", "15", "-sel1", "type 5 6", "-sel2", "type 6",
      "-o", "surf_so_water_o"]),
    (D22, "verify", "rdf untouched", ["rdf", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-cutoff", "15", "-sel1", "type 2 3 4", "-sel2", "type 7 8",
      "-o", "surf_org_oil"]),
    (D22, "verify", "stored value is 1/L too low (70x); density itself is untouched",
     ["density", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-axis", "z", "-dim", "2", "-sel", "type 7 8",
      "-o", "oil_density_2d"]),
    (D22, "verify", "stored value is 1/L too low (70x); density itself is untouched",
     ["density", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-frame_range", "160:200", "-axis", "z", "-dim", "2", "-sel", "type 2 3 4 5 6",
      "-o", "surfactant_density_2d"]),
    # 23: rdf/agg/cluster/msd variants
    (D23, "verify", "rdf untouched", ["rdf", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-sel1", "type 1", "-sel2", "type 2-5", "-frame_range", "100:200", "-o", "li_fsi"]),
    (D23, "verify", "agg -mode min untouched", ["agg", "-s", "run.out.data",
      "-f", "run.out.lammpstrj", "-sel1", "type 1", "-sel2", "type 2-5", "-o",
      "li_fsi_lowercutoff", "-cutoff", "2.4", "-frame_range", "100:200"]),
    (D23, "verify", "cluster untouched (reports size-1 components)", ["cluster", "-s",
      "run.out.data", "-f", "run.out.lammpstrj", "-sel1", "type 1", "-sel2", "type 2-5",
      "-cutoff", "2.95", "-frame_range", "100:200", "-o", "li_fsi"]),
    (D23, "verify", "msd untouched", ["msd", "-s", "run.out.data", "-f", "run.out.lammpstrj",
      "-sel", "type 1", "-o", "Li"]),
    # 24 licl / cpm: density variants
    (D24, "verify", "density untouched", ["density", "-s", "prerun.out.data",
      "-f", "run.efield.out.lammpstrj", "-sel", "type 11", "-frame_range", "500:600",
      "-o", "cl_end", "-nbins", "200"]),
    (D24C, "verify", "density untouched", ["density", "-s", "prerun.out.data",
      "-f", "run.cpm.out.lammpstrj", "-sel", "type 10", "-frame_range", "500:600",
      "-o", "li_end", "-nbins", "200"]),
    (D24C, "verify", "density untouched", ["density", "-s", "prerun.out.data",
      "-f", "run.cpm.out.lammpstrj", "-sel", "type 11", "-frame_range", "500:600",
      "-o", "cl_end", "-nbins", "200"]),
    (D24C, "verify", "msd untouched", ["msd", "-s", "prerun.out.data",
      "-f", "run.cpm.out.lammpstrj", "-o", "li", "-sel", "type 10"]),
    # 24 lifsi / vaccum
    ("24_electrolyte_electrode/sandwich_lifsi", "verify", "density untouched",
     ["density", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-sel", "type 11-14", "-frame_range", "300:350", "-o", "tmp_end"]),
    ("24_electrolyte_electrode/sandwich_lifsi", "verify", "density untouched",
     ["density", "-s", "prerun.out.data", "-f", "run.efield.out.lammpstrj",
      "-sel", "type 11-14", "-frame_range", "0:20", "-o", "tmp_start"]),
    ("24_electrolyte_electrode/vaccum_lifsi", "verify", "density untouched",
     ["density", "-s", "prerun.out.data", "-f", "prerun.out.lammpstrj",
      "-sel", "type 9-26", "-frame_range", "90:100", "-o", "prerun_final", "-nbins", "200"]),
    ("24_electrolyte_electrode/vaccum_lifsi", "verify", "density untouched",
     ["density", "-s", "prerun.out.data", "-f", "prerun.out.lammpstrj",
      "-sel", "type 9-26", "-frame_range", "0:10", "-o", "prerun_init", "-nbins", "200"]),
    ("24_electrolyte_electrode/vaccum_lifsi", "verify", "density untouched",
     ["density", "-s", "prerun.out.data", "-f", "run.out.lammpstrj",
      "-sel", "type 9-26", "-frame_range", "180:200", "-o", "run_final", "-nbins", "200"]),
    ("24_electrolyte_electrode/vaccum_lifsi", "verify", "density untouched",
     ["density", "-s", "prerun.out.data", "-f", "run.out.lammpstrj",
      "-sel", "type 9-26", "-frame_range", "0:20", "-o", "run_init", "-nbins", "200"]),
    # 32: the second density map
    ("32_molecular_crystal/run_constraincof", "verify", "density untouched",
     ["density", "-s", "prerun.freezecof.out.data", "-f", "run.out.lammpstrj",
      "-sel", "type 11", "-o", "cl", "-dim", "2", "-axis", "z"]),
]

# `probe` was renamed, so the old file name is not derivable from the new tag.
OLD_NAME_OVERRIDE = {
    (D31, "freevol", "run"): {"freevol": "run.probe.ffv.out.csv",
                              "freevol.psd": "run.probe.psd.out.csv",
                              "freevol.connectivity": "run.probe.connectivity.out.csv"},
    (D31, "sasa", "run"): {"sasa": "run.probe.sasa.out.csv"},
    (D22, "density", "oil_density_2d"): {"density.2d": "oil.density.2d.out.csv"},
}

# Old artifacts that were renamed by hand after the run, and the mirror file
# that now stands in for them (kept under the *original* name so the copy-back
# is one-to-one).
# values are paths relative to this directory
MIRROR_ALIAS = {
    (D22, "density", "oil_density_2d"):
        {"density.2d": f"{D22}/oil.density.2d.out.csv"},
    ("32_molecular_crystal/run_constraincof", "density", "li"):
        {"density.2d": "32_molecular_crystal/li.density.2d.out.csv"},
    ("32_molecular_crystal/run_constraincof", "density", "cl"):
        {"density.2d": "32_molecular_crystal/cl.density.2d.out.csv"},
}

# Artifacts whose command is chemically void (no polar-H donor in -sel1): lat
# now refuses, and gmx reports "No Donors found" for the same groups. The old
# files are the C-H...O artifact and must be deleted rather than regenerated.
INVALID = [
    (D22, f"hbond_surfactant_water.hbond{t}.out.csv",
     "no polar-H donor in sel1 (surfactant is C/H/O/S only, 1920 C-H): 0.1.0 counted "
     "C-H...O artifacts; lat 0.2.0 refuses and gmx prints 'No Donors found'")
    for t in ("", ".acf", ".dist_xy", ".dist_xz", ".dist_yz")
]

# The 32 density artifact was moved up out of the run directory after the fact,
# so its old copy lives one level above the directory the command runs in.
OLD_DIR_OVERRIDE = {
    ("32_molecular_crystal/run_constraincof", "li"): "32_molecular_crystal",
    ("32_molecular_crystal/run_constraincof", "cl"): "32_molecular_crystal",
}

# Round-1 exploration artifacts that are not deliverables at all.
SCAFFOLDING = [
    f"{D31}/latfix.gyrate.out.csv", f"{D31}/latfix.rmsd.out.csv", f"{D31}/latfix.rmsf.out.csv",
    f"{D31}/run.ffv.out.csv", f"{D31}/run.psd.out.csv", f"{D31}/run.sasa.out.csv",
    f"{D31}/run.ffv.out.png", f"{D31}/run.psd.out.png", f"{D31}/run.sasa.out.png",
    f"{D31}/latcheck.fix.gyrate.csv", f"{D31}/latcheck.fix.maxedgejump.csv",
    f"{D31}/latcheck.fix.mol_rg.csv", f"{D31}/latcheck.fix.rmsd.csv",
    f"{D31}/latcheck.fix.rmsf.csv", f"{D31}/latcheck.gmx.csv",
    f"{D31}/latcheck.lat.gyrate.csv", f"{D31}/latcheck.lat.maxedgejump.csv",
    f"{D31}/latcheck.lat.mol_rg.csv", f"{D31}/latcheck.lat.rmsd.csv",
    f"{D31}/latcheck.lat.rmsf.csv", f"{D31}/latcheck_bond_order.patch",
    f"{D31}/latcheck_gmx_compare.py", f"{D31}/latcheck_jumps.py",
    f"{D31}/latcheck_probe.py", f"{D31}/latcheck_verify.py",
    f"{D31}/gmx_convert3.py", f"{D31}/bonds_fixed.tmp.npy",
    f"{D31}/run.ffv.out.png", f"{D31}/run.psd.out.png", f"{D31}/run.sasa.out.png",
]


def mirror_dir(wd: str) -> Path:
    return HERE / wd


def retarget(argv: list[str], wd: str) -> list[str]:
    out = list(argv)
    i = out.index("-o")
    out[i + 1] = str(mirror_dir(wd) / out[i + 1])
    return out


def job_id(wd: str, argv: list[str]) -> str:
    """Unique per command: workdir, module, -o prefix and any -mode value."""
    parts = [wd.replace("/", "_"), argv[0], argv[argv.index("-o") + 1]]
    if "-mode" in argv:
        parts.append(argv[argv.index("-mode") + 1])
    return "__".join(parts)


def run_job(wd: str, status: str, why: str, argv: list[str], force: bool) -> dict:
    jid = job_id(wd, argv)
    log = HERE / "logs" / f"{jid}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    mirror_dir(wd).mkdir(parents=True, exist_ok=True)
    cmd = [PY, str(LAT), *retarget(argv, wd)]
    # a `.done` marker with rc=0 is the only proof of a completed run; an
    # interrupted job leaves a log behind but must not be treated as cached
    done = HERE / "logs" / f"{jid}.done"
    if force or not (done.exists() and done.read_text().strip() == "0"):
        with open(log, "w") as fh:
            proc = subprocess.run(cmd, cwd=EXAMPLES / wd, stdout=fh, stderr=subprocess.STDOUT)
        rc = proc.returncode
        (HERE / "logs" / f"{jid}.cmd").write_text(" ".join(cmd) + "\n")
        done.write_text(str(rc) + "\n")
    else:
        rc = 0
        print(f"  [skip: cached] {jid}")
    return {"id": jid, "wd": wd, "status": status, "why": why,
            "argv": " ".join(argv), "rc": rc, "log": str(log.relative_to(HERE))}


def compare(old: Path, new: Path) -> dict:
    """Value comparison of two LAT CSVs, by *position* (NaN-aware).

    Stored example outputs were written by earlier releases whose column
    *names* and, in two cases, row counts differ; comparing by position is what
    answers "did the numbers change under the current code".
    """
    import numpy as np

    def load(p):
        names = None
        for line in p.read_text().splitlines():
            if line.startswith("# columns:"):
                names = [s.strip() for s in line.split(":", 1)[1].split(",")]
        data = np.atleast_2d(np.genfromtxt(p, delimiter=",", comments="#", ndmin=2))
        if data.size == 0:
            data = np.empty((0, len(names or [])))
        return names or [], data

    na, da = load(old)
    nb, db = load(new)
    notes = []
    if na != nb:
        notes.append(f"column names changed: {'+'.join(na)} -> {'+'.join(nb)}")
    out = {"cols_old": ",".join(na), "cols_new": ",".join(nb),
           "rows_old": da.shape[0], "rows_new": db.shape[0], "max_rel": None,
           "note": "; ".join(notes)}
    if da.shape != db.shape:
        if da.shape[0] != db.shape[0]:
            out["note"] = "; ".join(notes + [f"shape changed {da.shape} -> {db.shape}"])
            return out
        notes.append(f"shape changed {da.shape} -> {db.shape}; compared the common columns")
        k = min(da.shape[1], db.shape[1])
        da, db = da[:, :k], db[:, :k]
        out["note"] = "; ".join(notes)
    if da.size == 0:
        out["max_rel"] = 0.0
        return out
    with np.errstate(invalid="ignore", divide="ignore"):
        rel = np.abs(db - da) / np.maximum(np.abs(da), 1e-12)
    rel = np.where(np.isnan(da) & np.isnan(db), 0.0, rel)
    rel = np.where(np.isnan(rel), np.inf, rel)
    out["max_rel"] = float(rel.max()) if rel.size else 0.0
    out["changed_cols"] = ",".join(name for k, name in enumerate(na)
                                   if rel[:, k].max() > 1e-9) if rel.size else ""
    return out


def output_names(wd: str, argv: list[str], prefix: str) -> list[str]:
    """New output file names, from the module's declared tags."""
    sys.path.insert(0, "/opt/lat/src")
    from lat.modules import get

    mod = get(argv[0])
    return [f"{prefix}.{t}.out.csv" for t in _tags(wd, argv, mod)]


def _tags(wd: str, argv: list[str], mod) -> list[str]:
    """Tags a module writes for this argv (only the ones this project regenerates)."""
    name = argv[0]
    if name == "rms":
        return ["rmsd", "rmsf"]
    if name == "freevol":
        mode = "ffv" if "-mode" not in argv else argv[argv.index("-mode") + 1]
        return {"ffv": ["freevol"], "psd": ["freevol.psd"],
                "connectivity": ["freevol.connectivity"]}[mode]
    if name == "sasa":
        return ["sasa"]
    return {
        "gyrate": ["gyrate"], "rdf": ["rdf", "rdf.cn"], "contact": ["contact"],
        "hbond": ["hbond", "hbond.acf", "hbond.dist_xy", "hbond.dist_xz", "hbond.dist_yz"],
        "density": ["density.2d" if "-dim" in argv and argv[argv.index("-dim") + 1] == "2"
                    else "density"],
        "msd": ["msd"], "agg": ["agg", "agg.labels", "agg.speciation"],
        "cluster": ["cluster.cn", "cluster.size"],
    }[name]


def main() -> int:
    ids = [job_id(wd, argv) for wd, _, _, argv in JOBS]
    assert len(set(ids)) == len(ids), (
        "job ids collide, two jobs would share one log and one would be skipped: "
        + ", ".join(sorted({i for i in ids if ids.count(i) > 1}))
    )
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--affected", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--force", action="store_true", help="re-run even if a log exists")
    ap.add_argument("--filter", default=None, help="only jobs whose workdir contains this")
    a = ap.parse_args()

    if a.list:
        for wd, status, why, argv in JOBS:
            print(f"{status:8s} {wd:55s} {' '.join(argv)}")
            print(f"{'':8s} {'':55s} -> {why}")
        return 0

    selected = [j for j in JOBS
                if ((a.affected and j[1] == "affected")
                    or (a.verify and j[1] == "verify"))
                and (a.filter is None or a.filter in j[0])]
    if not selected:
        selected = [j for j in JOBS if j[1] == "affected"] if not a.report else []

    rows = []
    for wd, status, why, argv in selected:
        print(f"[run] {status} {wd}: {' '.join(argv)}")
        info = run_job(wd, status, why, argv, a.force)
        print(f"      rc={info['rc']} log={info['log']}")
        rows.append(info)

    if a.report or rows:
        import csv

        lines = []
        for wd, status, why, argv in JOBS:
            prefix = argv[argv.index("-o") + 1]
            old_wd = EXAMPLES / OLD_DIR_OVERRIDE.get((wd, prefix), wd)
            old_name = OLD_NAME_OVERRIDE.get((wd, argv[0], prefix), {})
            alias = MIRROR_ALIAS.get((wd, argv[0], prefix), {})
            for tag in _tags(wd, argv, None):
                new = (HERE / alias[tag]) if tag in alias else \
                    (mirror_dir(wd) / f"{prefix}.{tag}.out.csv")
                old = old_wd / old_name.get(tag, f"{prefix}.{tag}.out.csv")
                if new.exists() and old.exists():
                    cmp = compare(old, new)
                    lines.append([wd, status, str(old.relative_to(EXAMPLES)),
                                  str(new.relative_to(HERE)),
                                  argv[0], why, cmp["rows_old"], cmp["rows_new"],
                                  f"{cmp['max_rel']:.3e}" if cmp["max_rel"] is not None else "n/a",
                                  cmp.get("changed_cols", ""),
                                  cmp.get("note", "") or "mirror-only"])
                elif new.exists():
                    lines.append([wd, status, "", str(new.relative_to(HERE)), argv[0], why,
                                  "", "", "", "", "new file (no old counterpart)"])
        with open(HERE / "MAPPING.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["workdir", "status", "old_file", "new_file", "module", "reason",
                        "rows_old", "rows_new", "max_rel_diff", "changed_cols", "note"])
            w.writerows(lines)
            for s in SCAFFOLDING:
                w.writerow([str(Path(s).parent), "scaffolding", s, "", "",
                            "round-1 exploration artifact: delete", "", "", "", "", "delete"])
            for wd, name, why in INVALID:
                w.writerow([wd, "invalid", f"{wd}/{name}", "", "", why, "", "", "", "", "delete"])
        print(f"\nwrote {HERE / 'MAPPING.tsv'}")
        for ln in lines:
            if ln[1] in ("verify",) and ln[8] not in ("0.000e+00", "n/a", ""):
                print(f"  !! {ln[0]}/{ln[2]}: max_rel_diff={ln[8]} (expected 0)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

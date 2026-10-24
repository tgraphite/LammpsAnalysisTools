#!/usr/bin/env python3
"""Before/after comparison figures for a LAT behaviour change.

Given a tree of *old* LAT CSV outputs (already in the project) and a tree of
*new* ones produced by the changed build (see `recheck.py`), draw one figure per
quantity that materially changed, so the change can be reviewed at a glance
instead of by diffing numbers.

    make_before_after.py --examples /path/examples \
        --mirror /path/examples/_lat_v020_recheck \
        --out docs/figures/0.2.0

Only quantities listed in CASES are drawn; the mapping between old and new file
names (renames) is part of the case table so the table stays the single source
of truth for "what changed".
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# CJK labels: pick whichever Noto CJK face is installed, else fall back to
# English-only text (matplotlib would otherwise draw tofu boxes).
for _f in ("Noto Sans CJK JP", "Noto Sans CJK SC", "WenQuanYi Zen Hei"):
    if any(_f == f.name for f in matplotlib.font_manager.fontManager.ttflist):
        plt.rcParams["font.sans-serif"] = [_f, "DejaVu Sans"]
        break
plt.rcParams["axes.unicode_minus"] = False

# (relative path, [old name, new name], kind, title, unit, xlabel)
CASES = [
    ("31_linear_polymer/analysis", ["run.gyrate.out.csv"], "line",
     "Rg（4 条 PET 链整体）", "Å", "frame",
     "键表按 atom id 重排：旧值把非键合原子当键，Rg 高估 3.60×"),
    ("31_linear_polymer/analysis", ["run.rmsd.out.csv"], "line",
     "RMSD（相对第 0 帧）", "Å", "frame",
     "同一键序缺陷：RMSD 高估 3.50×"),
    ("31_linear_polymer/analysis", ["run.rmsf.out.csv"], "line",
     "RMSF（逐原子，6200 个）", "Å", "atom id",
     "同一键序缺陷：RMSF 高估 4.53×，并出现与结构无关的锯齿"),
    ("31_linear_polymer/analysis",
     ["run.probe.ffv.out.csv", "run.freevol.out.csv"], "line",
     "自由体积分数 FFV", "1", "frame",
     "占用判据本来就按轴 MIC 回绕 → 逐位不变（对照组）"),
    ("31_linear_polymer/analysis",
     ["run.probe.sasa.out.csv", "run.sasa.out.csv"], "line",
     "SASA（Shrake–Rupley）", "Å²", "frame",
     "邻域/遮挡改用最小镜像：对 gmx sasa 从 +14.24% 降到 +0.59%"),
    ("31_linear_polymer/analysis",
     ["run.probe.psd.out.csv", "run.freevol.psd.out.csv"], "line",
     "孔径分布 PSD", "1/Å", "孔半径 (Å)",
     "重定义：vdW 距离变换 + 平台去重 + Voronoi 归一（∫dV/dr dr = 1）"),
    ("31_linear_polymer/analysis",
     ["run.probe.connectivity.out.csv", "run.freevol.connectivity.out.csv"], "connectivity",
     "自由空间连通性", "1", "frame",
     "环面闭合 + 26 邻接 + 逐帧汇总（旧表 30812 行里 112 个分量是单体素粉尘）"),
    ("22_surfactant/8x_system", ["hbond_any.hbond.out.csv"], "line",
     "氢键数（全体×全体，32800 原子）", "1", "frame",
     "存储文件来自更早实现；gmx hbond 给 13808.5，0.2.0 给 13908.8（+0.7%）"),
    ("22_surfactant/8x_system", ["hbond_any.hbond.acf.out.csv"], "line",
     "氢键自相关 ACF（C_inter / C_cont）", "1", "lag (frame)",
     "旧表只有 11 行（40 帧窗口），现为 40 行"),
    ("22_surfactant/8x_system", ["hbond_any.hbond.dist_xy.out.csv"], "map2d",
     "氢键中点 xy 分布", "count", "x (Å)",
     "计数整体上移（同一批帧、同一批给体）"),
    ("22_surfactant/8x_system", ["oil.density.2d.out.csv"], "map2d",
     "油相 2D 密度", "kg/m³", "x (Å)",
     "旧文件少一个盒长因子：0.635 vs 解析值 44.33 kg/m³（×70）"),
    ("22_surfactant/8x_system", ["surfactant_density_2d.density.2d.out.csv"], "map2d",
     "表面活性剂 2D 密度", "kg/m³", "x (Å)",
     "同上：2.772 vs 解析值 193.56 kg/m³（×70）"),
    ("23_electrolyte", ["li_solv.cluster.size.out.csv"], "bars",
     "Li–溶剂 团簇尺寸分布", "p", "团簇尺寸",
     "旧表漏掉 size=1 一档：p 之和 0.2337（新表 1.0000）"),
    ("23_electrolyte", ["li_solv.rdf.cn.out.csv"], "line",
     "Li–溶剂 配位数 n(r)", "1", "r (Å)",
     "旧表只有 100 行（头部却写 nbins=200），现为 200 行"),
    ("23_electrolyte", ["li_fsi.rdf.cn.out.csv"], "line",
     "Li–FSI 配位数 n(r)", "1", "r (Å)",
     "旧表只有 99 行，现为 200 行"),
    ("24_electrolyte_electrode/sandwich_licl", ["any.hbond.out.csv"], "line",
     "氢键数（LiCl/电极体系，14456 原子）", "1", "frame",
     "存储文件来自更早实现；0.2.0 与 gmx 的差（+13.8%）全部来自受体集合（lat 认 Cl⁻）"),
]


def load_csv(path: Path):
    """(column names, data) from a LAT CSV."""
    names = None
    for line in path.read_text().splitlines():
        if line.startswith("# columns:"):
            names = [s.strip() for s in line.split(":", 1)[1].split(",")]
    data = np.atleast_2d(np.genfromtxt(path, delimiter=",", comments="#", ndmin=2))
    return names or [], data


def grid_from_xyz(data):
    """Rebuild a 2D map from a flat (x, y, value) table."""
    xs = np.unique(data[:, 0])
    ys = np.unique(data[:, 1])
    m = np.full((len(xs), len(ys)), np.nan)
    ix = np.searchsorted(xs, data[:, 0])
    iy = np.searchsorted(ys, data[:, 1])
    m[ix, iy] = data[:, 2]
    return xs, ys, m


def draw(args, rel_dir, names, kind, title, unit, xlabel, note):
    old_p = args.examples / rel_dir / names[0]
    new_p = args.mirror / rel_dir / (names[1] if len(names) > 1 else names[0])
    if not old_p.is_file() or not new_p.is_file():
        print(f"  skip (missing): {rel_dir}/{names[0]}")
        return None
    _, a = load_csv(old_p)
    _, b = load_csv(new_p)
    fig = plt.figure(figsize=(11, 5.2))
    fig.suptitle(f"{title}\n{note}", fontsize=11, y=0.985)
    ax = fig.add_subplot(111)

    if kind == "line":
        ax.plot(a[:, 0], a[:, 1], lw=1.0, color="#c0392b", label=f"0.1.0  (mean {np.nanmean(a[:, 1]):.4g})")
        ax.plot(b[:, 0], b[:, 1], lw=1.0, color="#1f6feb", label=f"0.2.0  (mean {np.nanmean(b[:, 1]):.4g})")
        if a.shape[1] > 2 and b.shape[1] > 2:
            ax.plot(a[:, 0], a[:, 2], lw=0.8, ls=":", color="#c0392b", label="0.1.0 col3")
            ax.plot(b[:, 0], b[:, 2], lw=0.8, ls=":", color="#1f6feb", label="0.2.0 col3")
        ax.set_xlabel(xlabel)
    elif kind == "bars":
        sizes = sorted(set(a[:, 0]) | set(b[:, 0]))
        xa = np.arange(len(sizes))
        va = {int(s): v for s, v in zip(a[:, 0], a[:, 2])}
        vb = {int(s): v for s, v in zip(b[:, 0], b[:, 2])}
        ax.bar(xa - 0.2, [va.get(int(s), 0) for s in sizes], 0.4, color="#c0392b", label="0.1.0")
        ax.bar(xa + 0.2, [vb.get(int(s), 0) for s in sizes], 0.4, color="#1f6feb", label="0.2.0")
        ax.set_xticks(xa, [int(s) for s in sizes])
        ax.set_yscale("log")
        ax.set_xlabel(xlabel)
    elif kind == "map2d":
        xs_a, ys_a, ma = grid_from_xyz(a)
        xs_b, ys_b, mb = grid_from_xyz(b)
        fig.clf()
        vmax = float(np.nanmax(mb))
        for k, (m, lab, xs, ys) in enumerate(((ma, "0.1.0", xs_a, ys_a),
                                              (mb, "0.2.0", xs_b, ys_b))):
            ax = fig.add_subplot(1, 2, k + 1)
            im = ax.imshow(m.T, origin="lower", vmin=0, vmax=vmax, aspect="equal",
                           extent=[xs[0], xs[-1], ys[0], ys[-1]])
            ax.set_title(f"{lab}  (max {np.nanmax(m):.3g})")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("y (Å)")
            if k == 0:
                notes = "两图共用色标；若左图几乎空白即为旧归一化缺陷的直观表现"
        fig.suptitle(f"{title}\n{note}\n{notes}", fontsize=11, y=0.99)
        fig.colorbar(im, ax=fig.axes, shrink=0.85, label=unit)
        out = args.out / (rel_dir.replace("/", "_") + "__" + names[0].replace(".out.csv", "") + ".png")
        fig.savefig(out, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"  wrote {out.relative_to(args.out.parent)}")
        return out
    elif kind == "connectivity":
        # old: one row per component; new: one row per frame
        fr_a, cnt_a = np.unique(a[:, 0], return_counts=True)
        ax.plot(fr_a, cnt_a, lw=1.0, color="#c0392b",
                label=f"0.1.0 分量数/帧 (max {cnt_a.max()})")
        ax.plot(b[:, 0], b[:, 1], lw=1.0, color="#1f6feb",
                label=f"0.2.0 分量数/帧 (max {int(b[:, 1].max())})")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("连通分量数/帧")
        ax2 = ax.twinx()
        ax2.plot(b[:, 0], b[:, 4], lw=1.2, ls="--", color="#2c7a7b",
                 label="0.2.0 最大分量体积占比 frac_largest")
        ax2.set_ylabel("frac_largest")
        ax2.set_ylim(0.998, 1.0005)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")
    ax.grid(alpha=0.25)
    if kind != "connectivity":
        ax.legend(fontsize=8)
    ax.set_ylabel(f"{title} ({unit})" if kind == "line" else unit)
    out = args.out / (rel_dir.replace("/", "_") + "__" + names[0].replace(".out.csv", "") + ".png")
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out.relative_to(args.out.parent)}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--examples", required=True)
    ap.add_argument("--mirror", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    a.examples, a.mirror, a.out = Path(a.examples), Path(a.mirror), Path(a.out)
    a.out.mkdir(parents=True, exist_ok=True)
    made = []
    for rel_dir, names, kind, title, unit, xlabel, note in CASES:
        print(f"[{title}]")
        p = draw(a, rel_dir, names, kind, title, unit, xlabel, note)
        if p:
            made.append(p)
    print(f"\n{len(made)} figure(s) in {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

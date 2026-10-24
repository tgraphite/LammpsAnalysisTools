# 08 · 交接 / 恢复笔记（compact 后用）

> 目的：让新会话在零上下文下快速接续。先读本文件，再按需读 `07_decision_log.md`、
> `MANUAL_TESTER.md`、`MANUAL_ALGORITHMS.md`、`03_io_selection.md`、`05_regression.md`。

## 0. 一句话状态

**全部 P0 模块已实现，`pytest tests` 161 passed（全部）。** 2026-10 完成了一轮「LAMMPS→GROMACS 原生对拍」：
新增 `tests/regression/convert/lammps_to_gmx.py` + `tests/regression/level1_algo/gmx_ref_lammps.py`，
在 4 链 PET（6200 原子 / 501 帧）上逐帧比 gyrate / rms / rmsf / sasa / rdf / hbond。
**修掉两个真 bug**：① `lammps_data.parse` 的 `bonds` 没跟着 id 排序重排（Rg 高估 3.6×、
RMSD 3.5×、RMSF 4.5×）；② `hbond._donors` 把 C–H 当给体（计数高一个量级）。
修正后 gyrate/RMSD/RMSF 与 gmx 一致到 1e-4 Å 量级、hbond 逐帧基本相同、RDF 配位数差 0.26%。
**0.2.0 已随后落地**（用户拍板：不留旧版本、接受数值变化并重新生成产物）：新增
`src/lat/core/periodic.py` 作为周期几何唯一入口（架构守卫 `test_no_ad_hoc_pbc.py`），
把 `probe` 拆成 `freevol` + `sasa`（`probe` 删除，无兼容壳），修掉 SASA 的 PBC（+14.2% → −0.2%…−0.8%）、
`contact`/`residence -mode center` 的 wrapped 质心、`pistack` 的裸环心、PSD 的周期 EDT/bin/归一化、
connectivity 的环面与 26 邻接，并统一到「整体平移后物理量不变」这一判据
（`tests/regression/level2_format/test_pbc_invariance.py`）。完整清单见 `07_decision_log.md` §7 与
`12_0.2.0_release.md`；examples 旧产物的重跑映射见 `examples/_lat_v020_recheck/MAPPING.tsv`。
下一步 = P1 审计 / `current` 与 gmx 对拍 / examples 映射确认后覆盖回原目录。

## 1. 环境（本机）

- Python：`/opt/miniconda3/bin/python`（MDAnalysis 2.10.0, ase 3.28, numpy 1.26, scipy 1.17,
  matplotlib 3.10, pyyaml 6.0, sklearn 1.8）。
- GROMACS：`/usr/bin/gmx` 2022.5（Debian；另一台曾用 `/usr/local/gromacs/bin/gmx` 2023.2）。
  回归脚本现在按 `$GMX` → `PATH` → 旧路径解析，不再硬编码。仅内部回归用，不对外暴露。
- LAMMPS：`/opt/lammps/build/lmp`（2023.2；`stress/atom`、`voronoi/atom`、REAXFF、USER-NEP 等）。
- 其它：VMD `/opt/vmd193/vmd`、Multiwfn、GPUMD、reacnetgenerator、packmol。

## 2. 仓库地图（`/opt/lat`）

- `config/lat.yaml` — 唯一配置来源（默认值/单位/vdW/元素质量/各模块默认）。
- `docs/00–10` — 来龙去脉、元要求、架构、IO、输出契约、回归、参考、决策日志、本文件、测试战役。
- `docs/MANUAL_ALGORITHMS.md` + `docs/module_contracts.md` — 算法原理与模块输入/输出/警报。
- `src/lat/` — 实现（core + 15 个模块 + 2 个绘图命令）。
- `tests/regression/convert/` — 双向转换器：`convert_gmx_to_lammps.py`（gmx→LAMMPS）与
  `lammps_to_gmx.py`（LAMMPS→gmx，为原生 gmx 对拍而建）。
- `tests/regression/level1_algo/` — `compare_gmx.py`（GROMACS 原生输入）与
  `gmx_ref_lammps.py`（LAMMPS 输入；含 gyrate/rms/rmsf/sasa/rdf/hbond）。
- `tests/regression/samples/sys1_reduced/` — 自包含回归样例（已落盘，自检全通过）。
  注意：它是 Li/FSI/FEC 电解液，只有 C–H，**没有电负性给体**，因此 hbond 在样本级必然报错。

## 3. 冻结决策（要点，完整见 07）

- 时间=**帧号**，忽略 LAMMPS 时间；寿命/滞留标「k × 帧间距」。
- 默认长度单位 **Å**；未给 `-unit` 报 **INFO**（不中断）。
- **不推断分子**；分子用 `mol`/`Bonds`；无交互输入。
- `type→element`：质量 ±0.2% 匹配 **IUPAC（ase.data）**；无解/多解 ERROR；**不能用 reax_tools 的 Li=7.0**。
- vdW 半径：`reax_tools::ELEMENT_ATOMIC_RADII`（Bondi，进 `lat.yaml:vdw`）；probe 用纯 vdW。
- rdf `cutoff` 默认 **10 Å**；agg/cluster `cutoff` **必填**；hbond 默认 **D···A 3.5 Å**、H–D–A **>150°**。
- probe 默认 **H2 探针 1.445 Å**（H2 动力学直径 2.89 Å / 2）。
- transport 需 dump 含 `vx vy vz`，否则 ERROR；默认去质心漂移。
- occupancy 语法 `-zone x/y/z:min-max[,min-max][ ...]`；分子体积自动 zone = **P2**。
- 输出命名 `<basename>.<tag>.<mode>.out.csv`；CSV 头已冻结；同名冲突改名 `.bak.1/.bak.2…`。
- **不用退出码**；INFO/WARNING/ERROR；全部 stderr；`-debug` 落 JSON 审计（对外不暴露）。
- 唯一 `lat.yaml` 在程序根目录（开发人员用）；所有物理参数有 CLI 开关。
- **不支持多文件/glob**（一轨迹一模拟）。
- 优先级：除 `order` 与分子体积 zone = **P2** 外，其余全部 **P0**；P1 待 P0 验收后审计。

## 4. 已验证事实 / 实现陷阱

- LAMMPS `data`/`dump` **都不携带 units 与 atom_style**；`atom_style` 由 `Atoms` 列数推断（见 `03`）。
- **6 列歧义**（charge vs molecular/bond/angle）需启发式，仍不确定 ERROR 要 `-style`。
- MDAnalysis 读 `.lammpstrj` **必须显式 `format="LAMMPSDUMP"`**。
- MDAnalysis 2.10 用 `AtomGroup.indices`（`.inds` 已废弃）。
- MDAnalysis 写不了 LAMMPS → 自研 writer 已解决。
- `.gro` 不带质量/电荷；其 Li 质量会被猜为 0 → 质量相关量以 `*.mapping.csv`/data 为准。
- `mda.Universe` 读 tpr 无 `dimensions`，需先 `u.trajectory[0]`。
- 离子缔合：**真·残基级二部图**；禁止 42 项目的贪心最近邻（低估 AGG）。

## 5. 回归样例

路径：`tests/regression/samples/sys1_reduced/`（800 分子/6600 原子，帧 `[0,20,40,60,80]`）。
- 生成：`convert_gmx_to_lammps.py -s <tpr> -f <xtc> -o sys1_reduced --per-resname 200 --frames 5` 再 `make_selections.py`。
- 校验：`verify_lammps_roundtrip.py -p sys1_reduced` → `ok: true`（坐标偏差 9.5e-7 Å）。
- 源数据：`/mnt/f/CRO/40_multi_electrolytes/11-systems/sys1_npt.{tpr,xtc}`；sys4/sys6 含 PMBA 聚合物。
- 缺口：无固-液/金属/强氢键/ππ；这些用 `tests/data/` 自造。

## 6. P0 实现顺序

1. ✅ `core`（io/topology/selection/pbc/neighbor/grid/correlate/output/context/config/errors）
   + `select` + `describe` + `plot1d/2d` + `cli`（`./lat`）——`config/lat.yaml` 单源、备份轮转、
   INFO/WARNING/ERROR、原子写、类型质量推断（IUPAC）。
2. ✅ `rdf`（第一条数值回归：lammps vs gmx 完全相同；vs MDAnalysis ≤0.05）。
3. ✅ `density`（1D/2D，kg/m³ 或 1/Å³）→ ✅ `msd`（解包优先 `xu`>`ix`>wrapped+WARN，Einstein D）→
   ✅ `rms`（Kabsch RMSD + RMSF）。单测：均匀密度、合成随机游走 MSD、刚体旋转 RMSD≈0。
4. ✅ `hbond`（D···A 3.5 Å / >150°，ACF+寿命+3 投影分布，需 Bonds）→ ✅ `agg`（残基二部图、
   SSIP/CIP/AGG + 可选 `-solv` speciation）→ ✅ `cluster`（CN/尺寸分布，`-focus` 退化 CN）→
   ✅ `pistack`（移植 83_pipi_stack：SVD 环法向 + 中心距离 + 夹角）。
5. ✅ `gyrate` `contact` `freevol`（ffv/psd/connectivity）`sasa` `residence`（inter/cont ACF、
   区间分布）`transport`（每物种 D + NE + GK，强制速度列）`occupancy`（zone 语法、穿越、首次通过、富集）。

**全部 P0 模块已实现；0.2.0 时 `pytest tests` 161 passed。** 含：
- 真实样例模块 E2E 冒烟（17 模块 + transport 守卫）；
- **level1 gmx 逐值回归**（`compare_gmx.py`，sys1–6 含 PMBA：RDF/CN ≤0.2%、密度 <1%、MSD ~5%）；
- 幂等测试 16 项；PBC 边界崩溃 / hbond 空帧 / transport 单位 / ACF 无偏等缺陷已修。
见 `docs/09_test_campaign.md`。剩余：P1 审计；agg 的 CN 口径开关（P1）；hbond/gyrate/sasa/current 逐值对拍。

工程约束（继续遵守）：模块代码不得出现物理默认值，全部来自 `config/lat.yaml`；`output.py` 是唯一写盘点；
新模块须加 pytest（`tests/unit` 或 `tests/regression/level2_format`）。
新增 CLI 子命令只需在 `src/lat/modules/<name>.py` 用 `@register`；`modules/__init__._load_builtin` 里补 import。
`tests/regression/level1_algo` 待接全量 11-systems 的 gmx `.xvg` golden（183k 原子，暂缓）。

## 7. 关键外部路径

- 回归数据：`/mnt/f/CRO/40_multi_electrolytes/11-systems`（不要复制全量进仓库，太大）。
- ππ 已验证代码：`/mnt/d/dev/83_pipi_stack_code`、`/mnt/d/dev/pi-stack/gmx_pistack`。
- reax_tools：`/mnt/d/CurrentProjects/reax_tools`（vdW 表 `src/cpp/argparser.cpp:62`）。
- 反面参考：`/mnt/f/CRO/42_multi_electrolytes_3/90-code/coord_distributions.ipynb`。
- 笔记/临时：`/tmp/opencode/`（element_mass.json、element_radii.json 等）。

## 8. 明确不做

能量/自由能、PMF、MMPBSA、对接、反应网络、固体/金属力学、DXA/位错、SDF、谱、
phonon、相场、DPD/CG、AIMD/QM、PKA。

## 9. 待办/不确定

- CI 规模（全量 sys vs 抽帧子集）：开发时实测，只要前后行为意图一致。
- P1 审计：P0 验收后判断是否有功能重复/必要性不足的模块。
- `order` 的参考轴定义（P2 再定）。

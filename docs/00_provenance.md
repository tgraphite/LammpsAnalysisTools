# 00 · 来龙去脉与决策记录

## 1. 这个项目为什么存在

起因是一次「MD 计算项清单」的梳理：最初有 **135 个商品名**，但逐一核对底层实现后，
发现它们并不是 135 种能力，而是少数算法内核 + 若干输入契约 + 大量输出切片/叙事。

关键结论（分析过程见对话与 06_reference）：

- 真正的**算法内核约 13 个**（近邻、直方图/体素、时间相关、Einstein/GK、图连通分量、
  拟合/降维、几何模板、压力张量、探针/表面、泊松、FFT、逐帧键、界面几何）。
- 收敛到**最小 8–9 个可复用子程序**即可覆盖约 80% 高频需求。
- 现有工程（GROMACS 侧）分析以 `gmx` 子命令为主，LAMMPS 侧则长期缺少一个统一、易用、
  可回归、输出契约固定的后处理工具——零散脚本各行其是（例如配位分布用贪心最近邻，破坏边数守恒）。

`lat` 就是为填补这个空缺而设计：**LAMMPS 专用、gmx 习惯、收敛后的子程序集、文档即教材、数据逻辑分离。**

## 2. 范围决策（已确认）

| 项 | 决策 |
|---|---|
| 输入 | 完全兼容 LAMMPS `data` + `lammpstrj`；内部可读 GROMACS（不对外暴露） |
| 回归 | **两级输入回归**：① 算法 GT（gmx xvg 经 MDAnalysis 直接比对）；② 格式兼容（转 LAMMPS） |
| 时间语义 | 全部用**帧号**；忽略 LAMMPS 时间；寿命/滞留以帧为单位，标注「k × 帧间距」 |
| 拓扑 | **不推断分子**；用户运行前用静态 MDAnalysis 语义或 `ndx` 明确计算目标 |
| 交互 | **无交互式输入**（避免 gmx 式交互破坏幂等：同一 id 10 可能是 protein 也可能是 SOL） |
| 固体/金属力学 | **不开发**（不属于本程序职责） |
| 能量 | **不开发**（只能在 LAMMPS on-the-fly 计算） |
| defect / sdf / spectra | **不开发**（交给 OVITO / 专用谱工具） |
| ππ 堆积 | 直接**移植**已有 100% 验证的代码（`/mnt/d/dev/83_pipi_stack_code`、`gmx_pistack`） |
| 氢键回归 | 只要求**数值一致**：自造计算对象，能算出点数即可（回归不考虑物理） |
| freevol / sasa | 0.1.0 里 `freevolume` + `pore` + `sasa` 集成在单一 `lat probe`；0.2.0 按算法族拆回 `lat freevol`（体素网格）与 `lat sasa`（曲面），见 MANUAL_ALGORITHMS.md §12 |
| zone/event | 更名为 `lat occupancy`（区域占据与穿越；见 MANUAL_ALGORITHMS.md） |
| transport | 显式要求 `lammpstrj` 含速度列（`vx vy vz`），否则 FATAL |

## 3. 元要求决策（M1–M15）

见 `01_meta_requirements.md`，其中 M1–M15 已全部确认；补充确认：

- M2：时间全部帧号；hbond 自相关-寿命输出为「k × 帧间距」。
- M5：保留 debug 审计能力，但**成品不对用户暴露无关审计信息**。
- M10：`lat select` 额外可导出 `ndx` 供校验。
- M3：选择器只接受静态 MDAnalysis 语义 / ndx；不解析复杂手工语义（交由 `group2ndx`）。

## 4. 命名

`lat <>`：LAMMPS Analysis Tool。子命令即模块名，短横线参数，单横线。

## 5. 参考代码索引（详见 06_reference.md）

- GROMACS 2023.2（本机 `/usr/local/gromacs/bin/gmx`）：`rdf -cn`、`density`、`hbond`、`msd`、`rms/rmsf`、`freevolume`、`sasa`。
- MDAnalysis 2.10（本机 miniconda）：RDF/MSD/HBond/RMS/PCA/Density、`lib.distances`。
- 已验证 ππ 代码：`/mnt/d/dev/83_pipi_stack_code`、`/mnt/d/dev/pi-stack/gmx_pistack`。
- 回归数据：`/mnt/f/CRO/40_multi_electrolytes/11-systems`（含 gmx GT：rdf/cn、density、msd）。
- 反面参考：`/mnt/f/CRO/42_multi_electrolytes_3/90-code/coord_distributions.ipynb`（贪心最近邻，不采纳其算法）。

# 06 · 参考代码与权威来源

> 原则：优先 **GROMACS 语义** 与 **权威文献/高质量开源**；in-house 代码仅作「能跑通」参考，不默认正确。

## 1. GROMACS（算法语义基准，本机 `/usr/local/gromacs/bin/gmx` 2023.2）

| 命令 | 借鉴点 |
|---|---|
| `gmx rdf -cn` | g(r) 归一化、运行配位数、`-rmax`、`-selrpos/-seltype`（排除分子内需自行处理） |
| `gmx density` | 沿轴质量/数量密度、分组 |
| `gmx hbond` | D–H···A 几何判据、`-ac` 存在函数自相关、`-life` 寿命（间断/连续两种解读）、`-dist/-ang` |
| `gmx msd` | 去 PBC 解包、Einstein 拟合、`-type`（3D/xy） |
| `gmx rms/rmsf` | Kabsch 对齐与逐原子涨落 |
| `gmx freevolume` | **探针插入**：以 vdW 半径判占据、`-radius` 探针半径（`lat probe` 的算法基础） |
| `gmx sasa` | 探针可及表面积（Shrake–Rupley 类） |

注意：GROMACS 的分析输出为 `.xvg`；lat 输出为自有 CSV 契约（见 04）。

## 2. MDAnalysis 2.10（本机 miniconda，可作交叉校验与 IO）

- `MDAnalysis.analysis.rdf`、`msd`、`hydrogenbonds`、`rms`、`pca`、`density`、`lineardensity`、
  `distances`、`dihedrals`、`encore`、`polymer`、`contacts`。
- `MDAnalysis.lib.distances`：PBC/三斜距离。
- 读取器：`XTCReader/TRRReader/TPRReader`（gmx 内部用）；`LAMMPS.DATAReader/DumpReader`（对外主路径）。

## 3. 已验证的 ππ 代码（只移植，不改算法）

- `/mnt/d/dev/83_pipi_stack_code/`（`src/analyzer.py`、`geometry.py`、`cell_list.py`、`topology_parser.py`、`tests/`）
- `/mnt/d/dev/pi-stack/gmx_pistack` 与 `pi_pi_stack_report.html`（含 result_loose/medium.csv 参考输出）

## 4. 反面参考（不采纳其算法）

- `/mnt/f/CRO/42_multi_electrolytes_3/90-code/coord_distributions.ipynb`：
  配位分布用「贪心最近邻一对一 + 配位数上限封顶」，破坏边数守恒，会系统性低估 AGG。
  `lat agg` 必须建**真正的残基级二部图**（见 MANUAL_ALGORITHMS.md）。

## 5. 离子缔合（SSIP/CIP/AGG）文献口径

- 阳离子中心：Li⁺ 第一壳阴离子数 0/1/≥2（Hou et al. *Chem. Sci.* 2021；Efaw et al. *Nat. Mater.* 2023）。
- 阴离子中心：阴离子配位 0/1/≥2 个阳离子（Henderson & Seo, *JES* 2022）。
- 1:1 电解质：**平均配位数相等（边数守恒），但 SSIP/CIP/AGG 分数一般不对称**。
- 实现必须报中心、截断、残基级计数，并做边数守恒自检。

## 6. LAMMPS 参考

- 本机 `/opt/lammps`（2023.2，packages: REAXFF/DPD/VORONOI/QEQ/MC/MANYBODY/USER-NEP…）。
- 相关 compute/fix：`compute rdf/msd/vacf/cluster/atom/coord/atom/voronoi/atom/gyration/stress/atom`、
  `fix ave/time/ave/chunk`、`fix gcmc`。
- 本工具只消费 `data`/`lammpstrj`，不驱动引擎。

## 7. vdW 半径表（采用）

- 来源：`reax_tools`（本机 `/mnt/d/CurrentProjects/reax_tools`，公认反应网络基准）。
- 文件：`src/cpp/argparser.cpp:62` 的 `ELEMENT_ATOMIC_RADII`（元素符号 → Å，Bondi 型）。
- `reax_tools` 用 `d < 0.6·(r_i+r_j)`（`RVDW_FACTOR=1.2`）做键感知；lat 的 **probe 探针占据用纯 vdW 半径，不加因子**。
- lat 的 vdW 表以 nm 为单位整表复制进 `config/lat.yaml:vdw`；按元素符号索引，需先建立 `type→element` 映射。

## 8. 回归数据

- `/mnt/f/CRO/40_multi_electrolytes/11-systems`（gmx GT：rdf/cn、density、msd）。
- 自包含缩减样例：`tests/regression/samples/sys1_reduced/`。

## 9. 转换器（本项目自研）

- `tests/regression/convert/convert_gmx_to_lammps.py`：MDAnalysis 读 gmx（tpr/xtc）→ 写 LAMMPS
  `data`（atom_style full）+ `lammpstrj`，可选缩减分子数与帧数，并输出 `.mapping.csv` / `.selections.yaml`。
- `tests/regression/convert/make_selections.py`：生成分子/元素/模块候选选择。
- `tests/regression/convert/verify_lammps_roundtrip.py`：坐标/盒/键/元素/选择一致性自检。
- 元素质量匹配用 **IUPAC 标准原子量（ase.data）**；vdW 用 reax_tools 表（**单位 Å**，非 nm——`06` 旧文误写 nm）。
- MDAnalysis 读取 LAMMPS dump 需显式 `format="LAMMPSDUMP"`。

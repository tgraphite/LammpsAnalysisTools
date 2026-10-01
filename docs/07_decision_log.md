# 07 · 决策日志与优先级

> 本文件记录已确认决策（开发前冻结）与优先级。所有「待定」项已在此关闭。

## 1. 优先级定义（用户口径）

| 级别 | 含义 |
|---|---|
| **P0** | 一会就要开发（所有已确认项） |
| **P1** | 待验收 P0 后开发（必要性 / 是否重复功能尚需判断） |
| **P2** | 以后验证了关键点再开发（算法或定义不确定） |

**全部已确认功能 = P0。** 唯一调整：
- `gyrate`(rg)、`contact` = **P0**；
- `order` = **P2**（参考轴定义不确定）；
- 以分子体积自动生成 `occupancy` zone = **P2**（多孔材料控制区指认问题）；
- 其余「必要性/重复性」未决项 = P1（当前为空，待 P0 验收后审计）。

## 2. P0 模块清单

`rdf` `density` `hbond` `rms` `msd` `agg` `cluster` `pistack` `probe` `residence`
`transport` `occupancy` `gyrate` `contact` `select` `plot1d` `plot2d` `describe`
+ `core`（io/topology/selection/pbc/neighbor/grid/correlate/output/errors/context）。

## 3. 已确认决策

### 单位 / 时间
- 默认长度 **Å**（real 与 metal 都是 Å；lj 太小众不默认）。未显式给 `-unit` → 报 **INFO**，按默认 real/Å。
- 时间 **全部帧号**，忽略 LAMMPS 时间。寿命/滞留以帧输出，标「k × 帧间距」。

### 输入 / 格式（含源码依据）
- `atom_style` 数据列（LAMMPS 2023.2 `fields_data_atom`）：
  `atomic`=id type x y z(5)；`charge`=id type q x y z(6)；
  `molecular/bond/angle`=id molecule type x y z(6)；`full`=id molecule type q x y z(7)。
  **5→atomic、7→full 可确定；6 列歧义**（charge vs molecular）→ 启发式判断，仍不确定 ERROR 要求 `-style`。
- **不支持多文件/glob 拼接**：一轨迹一模拟，用户自行保证。
- `frame_range`：`start:stop:stride`，默认全帧。
- MSD/transport 解包：`xu yu zu` > `ix iy iz` > 最小镜像增量+WARN。
- **dump 原子顺序按 `id` 对齐拓扑**：dump 顺序是 LAMMPS 内部顺序，与 data 文件（按 `id` 升序）不一致；
  逐帧按 `id` 重排，否则 `type` 与坐标错配、g(r) 失真（见 13_ethanol_water 案例）。
  缺 `id` 列 → 假定一致并 WARN；`id` 集合不匹配 → ERROR。
- **hbond 几何中点**：`mid = A + 0.5·MIC(D−A)` 再 wrap；曾误写为 `0.5·(MIC(D−A)+A)`，
  导致二维分布被压到 1/4 盒（见 14_surfactant 案例）。
- PBC 支持正交/三斜；`-pbc` 隐藏默认开。
- **GROMACS 可内部直读，不对外暴露**（供回归）。

### 拓扑 / 元素
- **不推断分子**；分子内对需用户提供 `mol` 列或 GROMACS 残基编号。
- **rdf 默认排除同分子内对**（不计算分子内 RDF）；`-intramol` 打开（包含分子内）。
  有键但无分子 ID 时无法排除 → WARN，须显式 `-intramol`。
- `type→element`：默认按质量 **±0.2%** 匹配元素表；多解/无解 → **ERROR 中断**。
  行尾注释不可靠，不采用；可用 `-types` 文件显式覆盖。
- vdW 半径表采用 `reax_tools` 的 `ELEMENT_ATOMIC_RADII`（Bondi 型），整表进 `lat.yaml:vdw`；
  probe 用纯 vdW，键感知才用 `0.6*(ri+rj)`。

### 选择器
- 静态 MDAnalysis 语义 + `ndx`；无交互输入。
- 读入 `-index`/`-n` 时，`-sel1 -sel2` 等**允许是 ndx 中的组名**。
- **selection 字面值只存在于输出文件头**（因其可能为 `file.ndx -> water` 或 `type 1 2 3`）。
- `lat select` 打印命中统计，可导出 ndx。

### 算法口径
- **绝对轴分布（density 1D/2D、hbond 几何）**：照 **gmx density 口径**——逐帧按**分数坐标**分箱、
  轴按**平均盒长**标定；相关轴盒长偏离均值 **>±2% WARN**。非正交/畸变盒的截面近似由用户负责。
- **平滑只在绘图**：任何计算命令都输出**原始结果**（含 density；ACF 算不出寿命也照常输出）；
  `-smooth {none,gaussian,savgol}` 只在 `plot1d`/`plot2d` 生效。
- **rdf**：`cutoff` 默认 **10 Å**（不随盒子膨胀），n bins 200；`sel1≡sel2` 默认不算自对。
- **hbond**：对齐 gmx；**D···A 3.5 Å**（`-da` 默认）、H–D–A **>150°**；两组须完全相同或完全不重叠。
- **agg/cluster**：`cutoff` **必填**（不设默认）；真·残基级二部图；边数守恒 WARN。
- **probe**：默认 **H2 探针，半径 1.445 Å**（H2 动力学直径 2.89 Å 的一半）；PSD 用探针-网格近似并声明。
- **transport**：强制 `lammpstrj` 含速度列，否则 ERROR；默认去质心漂移；NE + GK 都给。
- **occupancy**：区域语法 `-zone x/y/z:min-max[,min-max][ ...]`，支持多维度约束；边界半开 `[a,b)`。
- **msd**：默认 3D、去质心漂移。

### 输出
- 命名：`<basename>.<tag>.<mode>.out.csv`，命名即意图。
  例：`mysystem.rdf.cn.out.csv`、`mysystem.rg.out.csv`、`mysystem.rdf_Li_FSI.out.csv`。
- **CSV 头冻结**：`program/module/command` + `selections/definitions/frame_range/scalars/warning`
  + **绘图文本契约** `title/xlabel/ylabel/columns`（`columns` 的 y 名即 legend）。
- **`plot1d` 输入用 `-f`**，不保留位置参数；多图/分组模式已放弃（不记入文档）。
- **`plot1d` 列选择/柱状**：`-columns` 按名只画部分 y 列（x 恒为第一列，永远绘制，不可选为 y）；
  `-bar` 分组柱、`-stack` 归一化堆积（占据比、y 轴 0–1）；`-bar`/`-stack` 互斥、仅限 `plot1d`、不平滑。
- **同名冲突**：旧文件改名 `.bak.1`、`.bak.2`…（不覆盖、不报错）。
- `-debug`：额外落完整 JSON 审计；**对外文档不暴露 debug 模式**。
- 不用退出码；分级 **INFO/WARNING/ERROR**（ERROR 中断），全部到 stderr。

### 配置 / 工程
- 有且仅有程序根目录下的**一个 `lat.yaml`**（开发人员用）；用户可配置的物理参数全部有命令行开关。
- 计算/绘图解耦（`plot1d`/`plot2d` 只读 CSV）。
- 回归：两级输入（算法 GT + 格式兼容）；**开发期重装 InterMol 或无头 VMD** 做 gmx→lammps 转换。

## 4. 遗留技术前置

- **gmx→lammps 转换链：已完成**，改为自研 MDAnalysis writer，无需 InterMol/VMD。
  脚本 `tests/regression/convert/`，样例 `tests/regression/samples/sys1_reduced/`（自检全通过）。
- 转换后原子顺序/type 编号的选择器映射：已随样例落盘（`*.mapping.csv`、`*.selections.yaml`）。
- CI 规模（全量 vs 抽帧子集）——开发时实测后定，只要前后行为意图一致。

## 5. 明确不做

能量/自由能、PMF、MMPBSA、对接、反应网络、固体/金属力学、DXA/位错、SDF、谱(VDOS/SAXS)、
phonon、相场、DPD/CG、AIMD/QM、PKA。
